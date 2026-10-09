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
 * Name of alumn:
 * Email of alumn:
 * Year: 2026
 *
 */

// This file is another solution to the meshTexture exercise.
// Use a sparse matrix, and a sparse solver to solve the system. It can
// successfully compute huge mesh efficiently.

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
  // TODO 2.1: Copy the your body of the updateEdgeLists() method from
  // meshBoundary.cpp
  //   Hash function for an edge. Only works with edges contaning vertices in
  //   [0,mesh.numvertex())
  struct SimpleEdgeHash {
    std::size_t numVertex;
    std::size_t operator()(const SimpleEdge &e) const {
      return min(e.a, e.b) + max(e.a, e.b) * numVertex;
    }
  };

  struct SimpleEdgeEqual {
    bool operator()(const SimpleEdge &lhs, const SimpleEdge &rhs) const {
      return (lhs.a == rhs.a && lhs.b == rhs.b) ||
             (lhs.a == rhs.b && lhs.b == rhs.a);
    }
  };

  std::unordered_set<SimpleEdge, SimpleEdgeHash, SimpleEdgeEqual> externalSet(
      0, SimpleEdgeHash{mesh.numVertex()}),
      internalSet(0, SimpleEdgeHash{mesh.numVertex()});

  for (auto &tri : mesh.triangles) {
    for (auto &edge : tri.edges()) {
      auto internalIt = internalSet.find(edge);
      if (internalIt != internalSet.end()) {
        throw("TODO 2.1 IMPOSIBLE: Input mesh is not a manifold");
      }

      auto externalIt = externalSet.find(edge);
      if (externalIt != externalSet.end()) {
        externalSet.erase(externalIt);
        internalSet.insert(internalIt, edge);
      } else {
        externalSet.insert(externalIt, edge);
      }
    }
  }

  externalEdges.clear();
  internalEdges.clear();
  externalEdges.reserve(externalSet.size());
  internalEdges.reserve(internalSet.size());

  vec3 minPos(INFINITY, INFINITY, INFINITY);
  unsigned minVertex = -1;

  auto updateMin = [&](unsigned v) {
    vec3 pos = mesh.coordinates[v];
    if ((pos.X < minPos.X) || (pos.X == minPos.X && pos.Y < minPos.Y) ||
        (pos.X == minPos.X && pos.Y == minPos.Y && pos.Z < minPos.Z)) {
      minPos = pos;
      minVertex = v;
    }
  };

  // Untangle the external edges
  std::unordered_map<unsigned, unsigned> edgeMap(externalSet.size());
  for (auto &exEdge : externalSet) {
    if (edgeMap.count(exEdge.a)) {
      edgeMap[exEdge.b] = exEdge.a;
      updateMin(exEdge.b);
    } else {
      edgeMap[exEdge.a] = exEdge.b;
      updateMin(exEdge.a);
    }
  }

  unsigned currVertex = minVertex;
  for (size_t i = 0; i < externalSet.size(); i++) {
    if (!edgeMap.count(currVertex)) {
      // New loop => Find new staring point
      minPos.set(INFINITY, INFINITY, INFINITY);
      minVertex = -1;
      for (auto &it : edgeMap)
        updateMin(it.first);
      currVertex = minVertex;
    }

    unsigned nextVertex = edgeMap.at(currVertex);
    externalEdges.push_back(SimpleEdge(currVertex, nextVertex));
    edgeMap.erase(currVertex);
    currVertex = nextVertex;
  }

  for (auto &edge : internalSet)
    internalEdges.push_back(edge);
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
    std::string filename("mallas/MaxPlanck.45kv.ply");
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

    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(mesh.triangles.size() * 6);

    // Out of order cotanget sum
    for (auto &tri : mesh.triangles) {
      auto addTriplets = [&](unsigned a, unsigned b, double cot) {
        triplets.emplace_back(a, b, cot);
        triplets.emplace_back(b, a, cot);
      };

      vec3 A = mesh.coordinates[tri.a];
      vec3 B = mesh.coordinates[tri.b];
      vec3 C = mesh.coordinates[tri.c];

      vec3 AB = A - B;
      vec3 AC = A - C;
      vec3 BC = B - C;

      double Adot = AB.dot(AC);
      double Bdot = (-AB).dot(BC);
      double Cdot = (-AC).dot(-BC);

      double Across = AB.cross(AC).module();
      double Bcross = AB.cross(BC).module();
      double Ccross = AC.cross(BC).module();

      double Acot = Adot / Across;
      double Bcot = Bdot / Bcross;
      double Ccot = Cdot / Ccross;

      addTriplets(tri.a, tri.b, Ccot);
      addTriplets(tri.a, tri.c, Bcot);
      addTriplets(tri.b, tri.c, Acot);
    }

    cout << "Done triplets. "
         << duration<float>(high_resolution_clock::now() - clock0).count()
         << " seconds" << endl;

    // Initialize matrix
    meshMatrix.setFromTriplets(triplets.begin(), triplets.end());
    // VERY IMPORTANT!!! Reserve 8 edges per vertex, estimate of max connectivity, speeds up insertion a lot
    meshMatrix.reserve(Eigen::VectorXi::Constant(numVertex, 8));

    cout << "Done matrix init. "
         << duration<float>(high_resolution_clock::now() - clock0).count()
         << " seconds" << endl;

    // Diagonal correction
    std::vector<double> sumRow(numVertex, 0.0);

    for (int k = 0; k < meshMatrix.outerSize(); ++k) {
      for (Eigen::SparseMatrix<double>::InnerIterator it(meshMatrix, k); it;
           ++it) {
        sumRow[it.row()] += it.value();
      }
    }

    for (size_t i = 0; i < numVertex; ++i) {
      meshMatrix.coeffRef(i, i) = -sumRow[i];
    }

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

    // Add frontier vertices
    for (auto &edge : externalEdges) {
      for (size_t i = 0; i < mesh.numVertex(); i++) {
        if (meshMatrix.coeff(edge.a, i) != 0) {
          meshMatrix.coeffRef(edge.a, i) = 0;
        }
      }
      meshMatrix.coeffRef(edge.a, edge.a) = 1;
    }

    // TODO 3.2: Patch Laplace matrix to generate a valid system of equations

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

    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver; // General solver, slow

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
