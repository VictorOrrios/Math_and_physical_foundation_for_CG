/*
 * meshTexture-Sparse.cpp
 *
 * Written by Jose Miguel Espadero <josemiguel.espadero@urjc.es>
 *
 * This code is written as material for the FMF class of the
 * Master Universitario en Informatica Grafica, Juegos y Realidad Virtual.
 * Its purpose is to be didactic and easy to understand, not hard optimized.
 *
 * //TODO: Fill-in your name and email
 * Name of alumn: Víctor Orrios Barón
 * Email of alumn: v.orrios.2026@alumnos.urjc.es
 * Year: 2026
 *
 */

// This file is another solution to the meshTexture exercise.
// Use a sparse matrix, and a sparse solver to solve the system. It can
// successfully compute huge mesh efficiently.

#include <utility>
#define _CRT_NONSTDC_NO_DEPRECATE
#include <chrono>
#include <iostream>
#include <unordered_map>
#include <unordered_set>
using namespace std::chrono;

#include "SimpleMesh.hpp"
#include "TextureMesh.hpp"

// Check if Eigen  (a standar Matrix Library) is included. If not,
// you can get it at http://eigen.tuxfamily.org

// <Eigen/Dense> is the module for dense (traditional) matrix and vector.
// You can get a quick reference for using Eigen dense objects at:
// http://eigen.tuxfamily.org/dox/group__QuickRefPage.html
#include "Eigen/Dense"

// <Eigen/Sparse> is the module for sparse matrix and vectors, which
// are used when most of elements of the matrix will store a 0.0 value.
// You can get a quick reference for using sparse objects at:
// http://eigen.tuxfamily.org/dox/group__SparseQuickRefPage.html
#include "Eigen/Sparse"
using Eigen::SparseMatrix;

// Solvers
#include "Eigen/SparseCholesky"
#include "Eigen/SparseLU"

#include <algorithm>
#include <atomic>
#include <omp.h>

// Mesh with the most vertex is Venus3000 with 1.5M vertex
// 22 bits to store an index => 44 bits per edge => 11bit radix sort is optimal
inline uint64_t edge64bit(uint32_t a, uint32_t b) {
  return a < b ? (uint64_t(a) << 22) | b : (uint64_t(b) << 22) | a;
}

inline SimpleEdge bitToEdge(uint64_t bits) {
  bits &= 0x7FFFFFFFFFFFULL; // Clean unique edge bit flag
  return SimpleEdge(bits >> 22, bits & 0x3FFFFF);
}

inline void parallelRadixSort11bit(std::vector<uint64_t> &data, int startBit,
                                   int numPasses) {
  constexpr int RADIX_BITS = 11;
  constexpr int RADIX = 1 << RADIX_BITS; // 2048
  constexpr int RADIX_MASK = RADIX - 1;

  const std::size_t n = data.size();
  if (n == 0 || numPasses <= 0)
    return;

  const int maxThreads = omp_get_max_threads();

  // Scratch buffer. We use a stable MSD radix sort implementation
  std::vector<uint64_t> buffer(n);
  // Frecuency histogram
  std::vector<std::array<uint32_t, RADIX>> hist(maxThreads);
  // Pass copy global offsets
  std::vector<std::array<uint32_t, RADIX>> offsets(maxThreads);

#pragma omp parallel shared(hist, offsets)
  {
    const int tid = omp_get_thread_num();
    const int numThreads = omp_get_num_threads();

    const std::size_t begin = n * tid / numThreads;
    const std::size_t end = n * (tid + 1) / numThreads;

    // Ping pong buffer flag
    bool isPassEven = true;

    for (int pass = 0; pass < numPasses; ++pass) {
      uint64_t *src = isPassEven ? data.data() : buffer.data();
      uint64_t *dst = isPassEven ? buffer.data() : data.data();
      const int shift = startBit + pass * RADIX_BITS;

      // (a) Frecuency counting
      auto &localHist = hist[tid];
      localHist.fill(0);
      for (std::size_t i = begin; i < end; ++i)
        localHist[(src[i] >> shift) & RADIX_MASK]++;

      // (b) Compute offsets (exclusive scan)
#pragma omp barrier
#pragma omp single
      {
        std::size_t bucketBase = 0;
        for (int bucket = 0; bucket < RADIX; ++bucket) {
          std::size_t pos = bucketBase;
          for (int t = 0; t < numThreads; ++t) {
            offsets[t][bucket] = pos;
            pos += hist[t][bucket];
          }
          bucketBase = pos;
        }
      } // Implicit barrier

      // (c) Stable scatter
      auto next = offsets[tid];
      for (std::size_t i = begin; i < end; ++i)
        dst[next[(src[i] >> shift) & RADIX_MASK]++] = src[i];

      // Flip ping pong
      isPassEven = !isPassEven;
#pragma omp barrier
    }

    // Just in case num of passes is odd
    if (!isPassEven) {
      for (std::size_t i = begin; i < end; ++i)
        data[i] = buffer[i];
    }
  }
}

void sortAndClassifyEdges(std::vector<uint64_t> &edges,
                          std::vector<SimpleEdge> &externalEdges,
                          std::vector<SimpleEdge> &internalEdges) {
  const std::size_t n = edges.size();

  // 1. Radix sort by bits 22 to 43 (lower index in an edge)
  parallelRadixSort11bit(edges, 22, 2);

  // Unique edge bit flag
  // Used instead of bool/bit vector for less memory and better cache
  constexpr uint64_t UNIQUE_EDGE_BIT_MASK = 0x8000000000000000ULL;

  const int maxThreads = omp_get_max_threads();

  // Count to resize output vectors. [0] = internal [1] = external
  std::vector<std::array<uint32_t, 2>> intExtCounts(maxThreads, {0, 0});

#pragma omp parallel shared(externalEdges, internalEdges, intExtCounts)
  {
    const int tid = omp_get_thread_num();
    const int numThreads = omp_get_num_threads();

    const std::size_t begin = n * tid / numThreads;
    const std::size_t end = n * (tid + 1) / numThreads;

    // 2. Skip first group in thread block range
    std::size_t skipStart = begin;
    if (skipStart > 0 &&
        (edges[skipStart] >> 22) == (edges[skipStart - 1] >> 22)) {
      const uint64_t a = edges[skipStart] >> 22;
      while (skipStart < end && (edges[skipStart] >> 22) == a)
        skipStart++;
    }

    // 3. Sort + Mark + Count. One pass
    {
      std::size_t it = skipStart;
      uint32_t localExt = 0;

      while (it < end) {
        // 3.a. Get edge group range
        const uint64_t a = edges[it] >> 22;
        std::size_t gEnd = it + 1;
        while (gEnd < n && (edges[gEnd] >> 22) == a)
          gEnd++;

        // 3.b. Insertion sort the lower 22 bit
        // Adapted from https://en.wikipedia.org/wiki/Insertion_sort pseudocode
        // Note: Since the higher 22 bit are already sorted direct comparison
        // w/o bitwise AND is optimal
        for (std::size_t k = it + 1; k < gEnd; k++) {
          const uint64_t v = edges[k];
          std::size_t j = k;
          while (j > it && edges[j - 1] > v) {
            edges[j] = edges[j - 1];
            j--;
          }
          edges[j] = v;
        }

        // 3.c. Mark unique edges with bit flag and count the number of edges
        // marked
        std::size_t p = it;
        while (p < gEnd) {
          if (p + 1 < gEnd && edges[p] == edges[p + 1]) {
            p += 2;
          } else {
            edges[p] |= UNIQUE_EDGE_BIT_MASK;
            localExt++;
            p++;
          }
        }

        it = gEnd;
      }

      // Record counts, compute internal count using total visited - external
      // marked
      intExtCounts[tid][0] = (it - skipStart - localExt) / 2;
      intExtCounts[tid][1] = localExt;
    }

    // 4. Compute output vector sizes and acummulate thread write offsets
#pragma omp barrier
#pragma omp single
    {
      std::size_t internalAcumm = 0;
      std::size_t externalAcumm = 0;
      // Acummulate and write
      for (int i = 0; i < numThreads; i++) {
        const auto auxInt = intExtCounts[i][0];
        const auto auxExt = intExtCounts[i][1];

        intExtCounts[i][0] = internalAcumm;
        intExtCounts[i][1] = externalAcumm;
        internalAcumm += auxInt;
        externalAcumm += auxExt;
      }
      internalEdges.clear();
      externalEdges.clear();
      internalEdges.resize(internalAcumm);
      externalEdges.resize(externalAcumm);
    } // Implicit omp barrier given by omp single

    // 5. Fill output vectors using write offsets
    {
      std::size_t it = skipStart;
      uint32_t extOff = intExtCounts[tid][1];
      uint32_t intOff = intExtCounts[tid][0];

      // Use bit flag to now next edge to visit
      while (it < end) {
        if (edges[it] & UNIQUE_EDGE_BIT_MASK) {
          externalEdges[extOff++] = bitToEdge(edges[it]);
          it++;
        } else {
          internalEdges[intOff++] = bitToEdge(edges[it]);
          it += 2;
        }
      }
    }
  }
}

/// Write an Eigen Matrix to a matlab file
void exportDenseToMatlab(const Eigen::MatrixXd &m, const std::string &filename,
                         const std::string matrixName = "A") {
  cout << "Export matrix " << matrixName << " to file: " << filename
       << std::endl;

  // Open the file as a stream
  ofstream os(filename.c_str());
  if (!os.is_open())
    throw filename + string(": Error creating the file");

  os << "# name: " << matrixName << std::endl
     << "# type: matrix" << std::endl
     << "# rows: " << m.rows() << std::endl
     << "# columns: " << m.cols() << std::endl;
  Eigen::IOFormat fmt(Eigen::StreamPrecision, Eigen::DontAlignCols, " ", "\n",
                      "", "", "", "\n");
  os << m.format(fmt);

  os.close();
  std::cout << "To import " << matrixName
            << " into matlab use the command: load(\"" << filename << "\")"
            << std::endl;

} // void exportDenseToMatlab (const &MatrixXd m, const std::string &filename)

/// Write a sparse Eigen Matrix to a matlab file
void exportSparseToMatlab(const Eigen::SparseMatrix<double> &m,
                          const std::string &filename,
                          const std::string matrixName = "A") {
  cout << "Export sparse matrix " << matrixName << " to file: " << filename
       << std::endl;

  // Open the file as a stream
  ofstream os(filename.c_str());
  if (!os.is_open())
    throw filename + string(": Error creating the file");

  os << "# name: " << matrixName << std::endl
     << "# type: sparse matrix" << std::endl
     << "# nnz: " << m.nonZeros() << std::endl
     << "# rows: " << m.rows() << std::endl
     << "# columns: " << m.cols() << std::endl;

  for (int k = 0; k < m.outerSize(); ++k) {
    for (Eigen::SparseMatrix<double>::InnerIterator it(m, k); it; ++it) {
      os << 1 + it.row() << " " << 1 + it.col() << " " << it.value()
         << std::endl;
    }
  } // for
  os.close();
  std::cout << "To import " << matrixName
            << " into matlab use the command: load(\"" << filename << "\")"
            << std::endl;
} // void exportSparseToMatlab (const Eigen::SparseMatrix<double> &m, ...

/// Update the contents of externalEdges and internalEdges
void updateEdgeLists(const SimpleMesh &mesh,
                     std::vector<SimpleEdge> &externalEdges,
                     std::vector<SimpleEdge> &internalEdges) {
  externalEdges.clear();
  internalEdges.clear();

  // TODO 2.1: Implement the body of the updateEdgeLists() method
  const std::size_t maxEdgeOccurrences = 3 * mesh.triangles.size();

  // Get all edges
  std::vector<uint64_t> edges(maxEdgeOccurrences);
#pragma omp parallel for schedule(static)
  for (std::size_t i = 0; i < mesh.triangles.size(); ++i) {
    const auto &tri = mesh.triangles[i];

    const std::size_t j = 3 * i;

    edges[j] = edge64bit(tri.a, tri.b);
    edges[j + 1] = edge64bit(tri.b, tri.c);
    edges[j + 2] = edge64bit(tri.c, tri.a);
  }

  // Radix sort
  sortAndClassifyEdges(edges, externalEdges, internalEdges);

  // Find the smallest vertex on the boundary.
  vec3 minPos(INFINITY, INFINITY, INFINITY);
  unsigned minVertex = static_cast<unsigned>(-1);

  auto updateMin = [&](unsigned v) {
    const vec3 &pos = mesh.coordinates[v];

    if ((pos.X < minPos.X) || (pos.X == minPos.X && pos.Y < minPos.Y) ||
        (pos.X == minPos.X && pos.Y == minPos.Y && pos.Z < minPos.Z)) {
      minPos = pos;
      minVertex = v;
    }
  };

  auto eraseOne = [](std::vector<unsigned> &vec, unsigned value) {
    auto it = std::find(vec.begin(), vec.end(), value);
    if (it != vec.end())
      vec.erase(it);
  };

  // Build the boundary edge map.
  std::unordered_map<unsigned, std::vector<unsigned>> edgeMap(
      externalEdges.size());

  for (const auto &exEdge : externalEdges) {
    edgeMap[exEdge.a].push_back(exEdge.b);
    edgeMap[exEdge.b].push_back(exEdge.a);
    updateMin(exEdge.a);
    updateMin(exEdge.b);
  }

  // Order external edges into boundary loops.
  unsigned currVertex = minVertex;
  for (std::size_t i = 0; i < externalEdges.size(); i++) {
    auto itCur = edgeMap.find(currVertex);

    if (itCur == edgeMap.end() || itCur->second.empty()) {
      // New loop, find new minimum
      minPos = vec3(INFINITY, INFINITY, INFINITY);

      bool found = false;
      for (const auto &entry : edgeMap) {
        if (!entry.second.empty()) {
          updateMin(entry.first);
          found = true;
        }
      }
      if (!found)
        break; // Just in case

      currVertex = minVertex;
      itCur = edgeMap.find(currVertex);
    }

    // Get next vertex.
    unsigned next = itCur->second.back();
    itCur->second.pop_back();

    // Remove the reverse edge from next.
    auto itNext = edgeMap.find(next);
    if (itNext != edgeMap.end())
      eraseOne(itNext->second, currVertex);

    // Save current edge
    externalEdges[i] = SimpleEdge(currVertex, next);
    currVertex = next;
  }
  // END TODO 2.1
} // void updateEdgeLists()

int main(int argc, char *argv[]) {
  // Set dumpMatrix to true for debugging. Writting matrix to file can take some
  // time
  bool dumpMatrix = false;

  try {
    // Set default input mesh filename
    // std::string filename("mallas/16Triangles.off");  //Minimal case test
    // std::string filename("mallas/mask2.ply"); // Easy case test
    // std::string filename("mallas/mannequin2.ply"); //Medium case test
    // std::string filename("mallas/laurana50k.ply"); //Really hard for dense
    // std::string filename("mallas/MaxPlanck.45kv.ply");
    std::string filename("mallas/Nefertiti.990kv.ply");
    if (argc > 1)
      filename = std::string(argv[1]);

    ///////////////////////////////////////////////////////////////////////
    // Step 1.
    // Read an input mesh
    SimpleMesh mesh;
    cout << "Loading file " << filename << endl;
    mesh.readFile(filename, false);

    cout << "Num vertex: " << mesh.numVertex()
         << " Num triangles: " << mesh.numTriangles()
         << " Unreferenced vertex: " << mesh.checkUnreferencedVertex() << endl;

    // Set dumpMatrix to true for debugging. Writting matrix to file can take
    // some time
    dumpMatrix = dumpMatrix || (mesh.numVertex() < 40);

    // Time measure
    high_resolution_clock::time_point clock0 = high_resolution_clock::now();

    ///////////////////////////////////////////////////////////////////////
    // Step 2.
    // Compute edge list and show it. This step should work if you correctly
    // finished the meshBoundary exercise.
    std::vector<SimpleEdge> externalEdges;
    std::vector<SimpleEdge> internalEdges;
    updateEdgeLists(mesh, externalEdges, internalEdges);

    // Count num of internal and external vertex
    size_t numVertex = mesh.coordinates.size();
    size_t numOfExternalVertex = externalEdges.size();
    // size_t numOfInternalVertex = numVertex - numOfExternalVertex;

    // Dump external edges
    cout << numOfExternalVertex << " vertex in external boundary: " << endl;
    if (numOfExternalVertex < 80) {
      for (auto &e : externalEdges) {
        cout << "[" << e.a << "->" << e.b << "] ";
      }
      cout << endl;
    }

    // Build the list of external vertex
    std::vector<size_t> externalVertex;
    externalVertex.reserve(externalEdges.size());
    for (auto &e : externalEdges) {
      externalVertex.push_back(e.a);
    }

    // Optimization: Keep a lookup table to ask if one vertex is external
    std::vector<bool> isExternal(numVertex, false);
    for (size_t i = 0; i < numOfExternalVertex; i++) {
      isExternal[externalVertex[i]] = true;
    }

    cout << "Done compute Boundary. "
         << duration<float>(high_resolution_clock::now() - clock0).count()
         << " seconds" << endl;

    ///////////////////////////////////////////////////////////////////////
    // Step 3.
    // Build Laplace matrix (step 1)
    cout << "Computing Laplacian matrix (sparse):" << endl;

    // We will use a sparse  matrix. If you want to use sparse matrix
    // you have to check
    // http://eigen.tuxfamily.org/dox/group__TutorialSparse.html and declare a
    // sparse matrix instead. Eigen::MatrixXd meshMatrix(numVertex, numVertex);
    Eigen::SparseMatrix<double> meshMatrix(numVertex, numVertex);

    // TODO 3.1: Build the cotangent Laplacian matrix
    //
    //       /\         .  Lij = cot(α) + cot(β) , if vertex i neighbour of j
    //      /β \        .
    //     /    \       .  Lii = -Sum Lij , diagonal element is the sum of row
    //    /      \      .
    //  vi--------vj    .
    //    \      /      .
    //     \    /       .
    //      \α /        .
    //       \/         .
    //

    std::vector<char> isBoundary(numVertex, 0);
    for (const auto &edge : externalEdges) {
      isBoundary[edge.a] = 1;
      isBoundary[edge.b] = 1;
    }

    const std::size_t numTri = mesh.triangles.size();
    const int nThreads = omp_get_max_threads();

    std::vector<std::vector<Eigen::Triplet<double>>> localTriplets(nThreads);
    std::vector<std::vector<double>> localSumRow(
        nThreads, std::vector<double>(numVertex, 0.0));

    for (int t = 0; t < nThreads; ++t)
      localTriplets[t].reserve((numTri / nThreads + 1) * 6);

    #pragma omp parallel
    {
      const int tId = omp_get_thread_num();
      auto &myTriplets = localTriplets[tId];
      auto &mySumRow = localSumRow[tId];

      #pragma omp for schedule(static)
      for (std::size_t i = 0; i < numTri; ++i) {
        const auto &tri = mesh.triangles[i];

        const vec3 &A = mesh.coordinates[tri.a];
        const vec3 &B = mesh.coordinates[tri.b];
        const vec3 &C = mesh.coordinates[tri.c];

        const vec3 AB = A - B;
        const vec3 AC = A - C;
        const vec3 BC = B - C;

        const double Adot = AB.dot(AC);
        const double Bdot = (-AB).dot(BC);
        const double Cdot = (-AC).dot(-BC);

        const double crossMag = AB.cross(AC).module();
        const double Acot = Adot / crossMag;
        const double Bcot = Bdot / crossMag;
        const double Ccot = Cdot / crossMag;

        auto addTriplets = [&](unsigned a, unsigned b, double cot) {
          mySumRow[a] += cot;
          mySumRow[b] += cot;
          if (!isBoundary[a]) myTriplets.emplace_back(a, b, cot);
          if (!isBoundary[b]) myTriplets.emplace_back(b, a, cot);
        };
        addTriplets(tri.a, tri.b, Ccot);
        addTriplets(tri.a, tri.c, Bcot);
        addTriplets(tri.b, tri.c, Acot);
      }
    }

    std::vector<double> sumRow(numVertex, 0.0);
    #pragma omp parallel for schedule(static)
    for (std::size_t v = 0; v < numVertex; ++v) {
      double s = 0.0;
      for (int t = 0; t < nThreads; ++t) s += localSumRow[t][v];
      sumRow[v] = s;
    }

    std::size_t totalTriplets = 0;
    for (const auto &v : localTriplets) totalTriplets += v.size();

    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(totalTriplets + numVertex);
    for (auto &v : localTriplets) {
      triplets.insert(triplets.end(), v.begin(), v.end());
      v.clear();
      v.shrink_to_fit();
    }
    for (std::size_t i = 0; i < numVertex; ++i)
      if (!isBoundary[i])
        triplets.emplace_back(i, i, -sumRow[i]);


    // build triplets
    meshMatrix.setFromTriplets(triplets.begin(), triplets.end());
    meshMatrix.reserve(Eigen::VectorXi::Constant(numVertex, 8));

    // END TODO 3.1

    cout << "Done Laplacian matrix (sparse). "
         << duration<float>(high_resolution_clock::now() - clock0).count()
         << " seconds" << endl;

    if (dumpMatrix)
      exportSparseToMatlab(meshMatrix, "laplacian.mat", "L");

    ///////////////////////////////////////////////////////////////////////
    // Step 5.1
    // Build system matrix
    cout << "Computing System matrix:" << endl;

    // TODO 3.2: Patch Laplace matrix to generate a valid system of equations

    // Add frontier vertices
    for (auto &edge : externalEdges) {
      meshMatrix.coeffRef(edge.a, edge.a) = 1;
    }
    //meshMatrix.prune(0.0);
    meshMatrix.makeCompressed();

    // END TODO 3.2

    cout << "Done System matrix. "
         << duration<float>(high_resolution_clock::now() - clock0).count()
         << " seconds" << endl;

    if (dumpMatrix)
      exportSparseToMatlab(meshMatrix, "systemMatrix.mat", "A");

    ///////////////////////////////////////////////////////////////////////
    // Step 5.2
    // Build UV values for external vertex, using the boundary of a square.
    cout << "Computing contour conditions:" << endl;

    // Note that this is a matrix of numvertex rows and 2 columns
    Eigen::MatrixX2d UV_0(numVertex, 2);
    UV_0.setZero();

    // TODO 3.3: Build a set of valid UV values for external vertex, mapping
    // the vertex to the boundary of a square of side unit.

    // Precalculate all distances of the frontier
    std::vector<double> edgeDistances(externalEdges.size());
    double edgeDistanceSum = 0;
    for (size_t i = 0; i < externalEdges.size(); i++) {
      const auto &edge = externalEdges[i];
      double d = mesh.distance(edge.a, edge.b);
      edgeDistances[i] = d;
      edgeDistanceSum += d;
    }

    // Normalize frontier distance to [0,1] and then to ([0,1],[0,1])
    double normalizationSum = 0;
    for (size_t i = 0; i < edgeDistances.size(); i++) {
      auto &edge = externalEdges[i];
      auto &d = edgeDistances[i];
      normalizationSum += d;
      double norm =
          normalizationSum / edgeDistanceSum; // Edge length normalized
      // norm = double(i)/double(edgeDistances.size()); // Equal distance
      double u, v;
      if (norm < 0.25) {
        u = 4.0 * norm;
        v = 0;
      } else if (norm < 0.5) {
        u = 1;
        v = 4.0 * (norm - 0.25);
      } else if (norm < 0.75) {
        u = 1.0 - 4.0 * (norm - 0.5);
        v = 1;
      } else {
        u = 0;
        v = 1.0 - 4.0 * (norm - 0.75);
      }
      if (u < 0.0 || v < 0.0) {
        cout << "Negative UVs at frontier vertex " << edge.a << " UV " << u
             << "," << v << endl;
      }
      UV_0.coeffRef(edge.a, 0) = u;
      UV_0.coeffRef(edge.b, 1) = v;
    }

    // END TODO 3.3

    cout << "Done contour conditions. "
         << duration<float>(high_resolution_clock::now() - clock0).count()
         << " seconds" << endl;

    if (dumpMatrix)
      exportDenseToMatlab(UV_0, "boundaryUVO.mat", "UV0");

    ///////////////////////////////////////////////////////////////////////
    // Step 6

    // We will use Eigen to solve the matrix from here.
    cout << "Solving the system using Eigen (sparse):" << endl;

    // Solve Ax = b; where A = meshMatrix, x=UV,  and b = UV_0
    Eigen::MatrixX2d UV;
    // TODO 3.4: Solve the system using the sparse matrix. Leave solution in UV
    // Check http://eigen.tuxfamily.org/dox/group__TopicSparseSystems.html

    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    //Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;

    solver.compute(meshMatrix);

    if (solver.info() != Eigen::Success) {
      throw std::runtime_error("Error factorizing matrix");
    }

    UV = solver.solve(UV_0);

    if (solver.info() != Eigen::Success) {
      throw std::runtime_error("Error solving matrix system");
    }

    // END TODO 3.4
    cout << "Done solving the system. "
         << duration<float>(high_resolution_clock::now() - clock0).count()
         << " seconds" << endl;

    // Dump the computed solution to the system
    if (dumpMatrix)
      exportDenseToMatlab(UV, "solutionUV.mat", "UV");

    ///////////////////////////////////////////////////////////////////////
    // Step 6.4

    // Create a planar mesh (reverse UV mesh)
    SimpleMesh planarMesh;
    // Use UV values as geometrical coordinates and same triangles that input
    // mesh
    planarMesh.coordinates.resize(numVertex);
    for (size_t i = 0; i < numVertex; i++) {
      planarMesh.coordinates[i].set(UV(i, 0), UV(i, 1), 0);
    }
    planarMesh.triangles = mesh.triangles;

    string output_UVMesh1 = "output_UVMesh1.ply";
    cout << "Saving parameterization mesh to " << output_UVMesh1 << endl;
    // planarMesh.writeFileOBJ(output_UVMesh1);
    planarMesh.writeFilePLY(output_UVMesh1);

    // Create a TextureMesh mesh with UV coordinates
    TextureMesh textureMesh;
    textureMesh.coordinates = mesh.coordinates;
    textureMesh.triangles = mesh.triangles;

    // Set image filename to be used as texture
    textureMesh.textureFile = "UVchecker.jpg";

    // Set UV as texture-per-vertex coordinates
    textureMesh.UV.resize(numVertex);
    for (size_t i = 0; i < numVertex; i++)
      textureMesh.UV[i].set(float(UV(i, 0)), float(UV(i, 1)));

    // Dump textureMesh to file (.obj or .ply)
    string output_UVMesh2 = "output_UVMesh2.ply";
    cout << "Saving texture mesh to " << output_UVMesh2 << endl;
    // textureMesh.writeFileOBJ(output_UVMesh2);
    textureMesh.writeFilePLY(output_UVMesh2);

    // Visualize the file with an external viewer
#ifdef WIN32
    // string viewcmd = "\"C:\\Program Files
    // (x86)\\VCG\\MeshLab\\meshlab.exe\"";
    string viewcmd = "C:/meshlab/meshlab_32.exe";
#else
    string viewcmd =
        "/usr/bin/flatpak run --branch=stable --arch=x86_64 --command=meshlab "
        "--file-forwarding net.meshlab.MeshLab >/dev/null 2>&1 ";
#endif
    string cmd = viewcmd + " " + output_UVMesh2;
    cout << "Executing external command: " << cmd << endl;
    return system(cmd.c_str());

  } // try
  catch (const string &str) {
    std::cerr << "EXCEPTION: " << str << std::endl;
  } catch (const char *str) {
    std::cerr << "EXCEPTION: " << str << std::endl;
  } catch (std::exception &e) {
    std::cerr << "EXCEPTION: " << e.what() << std::endl;
  } catch (...) {
    std::cerr << "EXCEPTION (unknow)" << std::endl;
  }

#ifdef WIN32
  cout << "Press Return to end the program" << endl;
  cin.get();
#else
#endif

  return 0;
}
