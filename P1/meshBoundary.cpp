/*
 * meshBoundary.cpp
 *
 * Written by Jose Miguel Espadero <josemiguel.espadero@urjc.es>
 *
 * This code is written as material for the FMF class of the
 * Master Universitario en Informatica Grafica, Juegos y Realidad Virtual.
 * Its purpose is to be didactic and easy to understand, not hard optimized.
 *
 * This file computes the list of internal and external edges of a mesh
 * and writes the boundary (list of external edges). Then create a copy
 * of the mesh over a colorMesh and assign the red color to the vertex
 * in the boundary.
 *
 * //TODO: Fill-in your name and email
 * Name of alumn: Víctor Orrios Barón
 * Email of alumn: v.orrios.2026@alumnos.urjc.es
 * Year: 2026
 *
 */

#include "ColorMesh.hpp"
#include <cstddef>
#include <cstdint>
#include <queue>
#include <sys/types.h>
#include <unordered_map>
#include <unordered_set>
#ifdef _MSC_VER
#pragma warning(error : 4101)
#endif

#define _CRT_NONSTDC_NO_DEPRECATE
#include "SimpleMesh.hpp"
#include "TextureMesh.hpp"
#include <chrono>
#include <cmath>
#include <iostream>
using namespace std::chrono;

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

void sortAndClassifyEdges(std::vector<uint64_t> &edges,
                          std::vector<SimpleEdge> &externalEdges,
                          std::vector<SimpleEdge> &internalEdges) {
  const std::size_t n = edges.size();

  // Radix sort params, 11bit 2 pases
  constexpr int RADIX_BITS = 11;
  constexpr int RADIX = 1 << RADIX_BITS;
  constexpr int RADIX_MASK = RADIX - 1;
  constexpr int PASSES = 2;

  // Unique edge bit flag
  // Used instead of bool/bit vector for less memory and better cache
  constexpr uint64_t UNIQUE_EDGE_BIT_MASK = 0x8000000000000000ULL;

  const int maxThreads = omp_get_max_threads();

  // Scratch buffer. We use a stable MSD radix sort implementation
  std::vector<uint64_t> buffer(n);
  // Frecuency histogram
  std::vector<std::array<uint32_t, RADIX>> hist(maxThreads);
  // Pass copy global offsets
  std::vector<std::array<uint32_t, RADIX>> offsets(maxThreads);
  // Count to resize output vectors. [0] = internal [1] = external
  std::vector<std::array<uint32_t, 2>> intExtCounts(maxThreads, {0, 0});

#pragma omp parallel shared(hist, offsets, externalEdges, internalEdges,       \
                                intExtCounts)
  {
    const int tid = omp_get_thread_num();
    const int numThreads = omp_get_num_threads();

    const std::size_t begin = n * tid / numThreads;
    const std::size_t end = n * (tid + 1) / numThreads;

    // Ping pong buffer flag
    bool isPassEven = true;

    // 1. Radix sort by bits 22 to 43 (lower index in an edge)
    for (int pass = 0; pass < PASSES; ++pass) {
      uint64_t *src = isPassEven ? edges.data() : buffer.data();
      uint64_t *dst = isPassEven ? buffer.data() : edges.data();
      const int shift = 22 + pass * RADIX_BITS; // Start at bit 22

      auto &localHist = hist[tid];
      localHist.fill(0);

      // 1.a. Frecuency counting
      for (std::size_t i = begin; i < end; ++i)
        localHist[(src[i] >> shift) & RADIX_MASK]++;

      // 1.b. Compute offsets
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
      } // Implicit omp barrier given by omp single

      // 1.c. Copy from ping to pong buffer
      auto next = offsets[tid];
      for (std::size_t i = begin; i < end; ++i)
        dst[next[(src[i] >> shift) & RADIX_MASK]++] = src[i];

      // Ping pong swap
      isPassEven = !isPassEven;
#pragma omp barrier
    }

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
          if (edges[p] != edges[p + 1]) {
            edges[p] |= UNIQUE_EDGE_BIT_MASK;
            localExt++;
            p++;
          } else {
            p += 2;
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
        cout << "It " << i << " internal " << auxInt << " external " << auxExt
             << " intAcumm " << internalAcumm << " extAcumm " << externalAcumm
             << endl;
      }
      internalEdges.clear();
      externalEdges.clear();
      internalEdges.resize(internalAcumm);
      externalEdges.resize(externalAcumm);
      cout << "internalEdges.size " << internalEdges.size()
           << " externalEdges.size " << externalEdges.size() << endl;
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

void classifyEdges(const std::vector<uint64_t> &edges,
                   std::vector<SimpleEdge> &externalEdges,
                   std::vector<SimpleEdge> &internalEdges) {

  const std::size_t n = edges.size();
  const int maxThreads = omp_get_max_threads();

  std::vector<uint8_t> isUnique(n, 0); // True if is an internal edge
  std::vector<std::pair<uint, uint>> offsets(
      maxThreads,
      std::pair<uint, uint>(0, 0)); // First = internal, second = external

#pragma omp parallel shared(edges, externalEdges, internalEdges, isUnique,     \
                                offsets)
  {
    const int tid = omp_get_thread_num();
    const int numThreads = omp_get_num_threads();

    // Start and finish index, thread secuential
    const std::size_t begin = n * tid / numThreads;
    const std::size_t end = n * (tid + 1) / numThreads;
    std::size_t it = begin;

    // Pair counting
    while (it + 1 < end) {
      if (edges[it] == edges[it + 1]) {
        it += 2;
      } else {
        isUnique[it] = true;
        offsets[tid].second++;
        it++;
      }
    }

    // Start block ammend
    if (begin > 0 && isUnique[begin] && (edges[begin] == edges[begin - 1])) {
      isUnique[begin] = false;
      offsets[tid].second--;
    }

    // End block ammend
    if (it < end && (tid == numThreads - 1 || edges[it] != edges[end])) {
      isUnique[it] = true;
      offsets[tid].second++;
    }

    // Get internal offset using number of external edges in section
    offsets[tid].first = end - begin - offsets[tid].second;

#pragma omp barrier

    // Compute global offsets
#pragma omp single
    {
      std::size_t internalAcumm = 0;
      std::size_t externalAcumm = 0;
      for (int i = 0; i < omp_get_num_threads(); i++) {
        auto aux = offsets[i];
        offsets[i].first = internalAcumm;
        offsets[i].second = externalAcumm;
        internalAcumm += aux.first;
        externalAcumm += aux.second;
      }

      // Resize to total edge count
      internalEdges.resize(internalAcumm);
      externalEdges.resize(externalAcumm);
    }

#pragma omp barrier

    // Fill vectors with edges
    {
      const int tid = omp_get_thread_num();
      const int numThreads = omp_get_num_threads();

      // Start and finish index, thread secuential
      const std::size_t begin = n * tid / numThreads;
      const std::size_t end = n * (tid + 1) / numThreads;
      std::size_t it = begin;

      while (it < end) {
        if (isUnique[it]) {
          externalEdges[offsets[tid].second++] = bitToEdge(edges[it]);
          it++;
        } else {
          internalEdges[offsets[tid].first++] = bitToEdge(edges[it]);
          it += 2;
        }
      }
    }
  }
}

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

  // Count non duplicates for boundary
  // classifyEdges(edges, externalEdges, internalEdges);

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
}

int main(int argc, char *argv[]) {
  try {
    // Set default input mesh filename
    // std::string filename("mallas/16Triangles.off");
    // std::string filename("mallas/mannequin.ply");
    // std::string filename("mallas/mask2.ply");
    // std::string filename("mallas/knot-hole.ply");
    std::string filename("mallas/Nefertiti.990kv.ply");
    // std::string filename("mallas/angel_kneeling.150kv.ply");
    // std::string filename("mallas/bunny.ply");

    if (argc > 1)
      filename = std::string(argv[1]);

    ///////////////////////////////////////////////////////////////////////
    // Read a mesh and write a mesh
    SimpleMesh mesh;
    cout << "Loading file " << filename << endl;
    mesh.readFile(filename, false);

    cout << "Num vertex: " << mesh.numVertex()
         << " Num triangles: " << mesh.numTriangles()
         << " Unreferenced vertex: " << mesh.checkUnreferencedVertex() << endl;

    // Init time measures used for profiling
    high_resolution_clock::time_point clock0 = high_resolution_clock::now();

    ///////////////////////////////////////////////////////////////////////

    // Compute the list of external and internal edges
    // Implement the body of the updateEdgeLists() method .
    std::vector<SimpleEdge> externalEdges;
    std::vector<SimpleEdge> internalEdges;
    updateEdgeLists(mesh, externalEdges, internalEdges);

    cout << "Done updateEdgeLists() "
         << duration<float>(high_resolution_clock::now() - clock0).count()
         << " seconds" << endl;
    cout << externalEdges.size() << " boundary edges" << endl;
    cout << internalEdges.size() << " internal edges" << endl;

    // Write external boundary and compute boundary length
    double boundaryLength = 0.0;
    cout << "Edges in the boundary:" << endl;
    for (const SimpleEdge &e : externalEdges) {
      cout << "[" << e.a << "->" << e.b << "] ";
      boundaryLength += mesh.edgeLength(e);
    }
    cout << endl;
    cout << "Boundary length: " << boundaryLength << endl;

    // Compute the Euler Characteristic for this mesh
    // For a conected mesh, the value means:
    // 2 -> The mesh is a closed surface with no holes (sphere-like topology)
    // 1 -> The mesh has one hole (disk-like topology)
    // 0 -> The has one handle (torus-like topology)
    //-2 -> The mesh has two holes (tube-like topology)
    int eulerCharacteristic = mesh.numVertex() + mesh.numTriangles() -
                              externalEdges.size() - internalEdges.size();
    cout << "Euler Characteristic: " << eulerCharacteristic << endl;

    // Vector to store the distance from a vertex to the nearest vertex in the
    // boundary
    std::vector<double> boundDist;
    // Index of the vertex with max distance to boundary
    unsigned deepestVertex = 0xDEADC0DE;

    // TODO 2.2:
    // Compute the distance from each vertex to the nearest vertex in the
    // boundary (euclidean distance measured along edges) and store it in the
    // boundDist vector. Store in deepestVertex the index of the vertex with
    // the maximum distance to boundary

    boundDist.assign(mesh.numVertex(), INFINITY);

    // Build adjacency matrix to construct graph

    // Adjacency matrix: adj[vertexA][conexionIdx] = std::pair(vertexB,
    // distance(vertexA,vertexB)) Conexion vector is sparse and not ordered =>
    // must iterate fully for query
    std::vector<std::vector<std::pair<unsigned, double>>> adj(mesh.numVertex());

    auto addEdge = [&](const SimpleEdge &e) {
      double w = mesh.edgeLength(e);

      adj[e.a].push_back({e.b, w});
      adj[e.b].push_back({e.a, w});
    };

    for (const auto &e : externalEdges)
      addEdge(e);

    for (const auto &e : internalEdges)
      addEdge(e);

    using qItem = std::pair<unsigned, double>;

    struct ComparePairInQueue {
      bool operator()(const qItem &lhs, const qItem &rhs) const {
        return lhs.second > rhs.second;
      }
    };

    std::priority_queue<qItem, std::vector<qItem>, ComparePairInQueue> queue;

    // Add source nodes
    for (auto &edge : externalEdges) {
      if (boundDist[edge.a] > 0.0) {
        boundDist[edge.a] = 0.0;
        queue.push({edge.a, 0.0});
      }
    }

    // Dijkstras
    while (!queue.empty()) {
      const qItem item = queue.top();
      queue.pop();

      if (boundDist[item.first] < item.second)
        continue;

      for (auto &conexion : adj[item.first]) {
        double newDist = conexion.second + item.second;
        if (boundDist[conexion.first] > newDist) {
          boundDist[conexion.first] = newDist;
          queue.push({conexion.first, newDist});
        }
      }
    }

    double maxDist = boundDist[0];
    for (size_t i = 0; i < boundDist.size(); i++) {
      double d = boundDist[i];
      if (d > maxDist) {
        maxDist = d;
        deepestVertex = i;
      }
    }

    if (maxDist == INFINITY) {
      cout << "meshBoundary IMPOSIBLE: No frontier found, fatal error" << endl;
      return 1;
    }

    // END TODO 2.2
    cout << "Done boundDist() "
         << duration<float>(high_resolution_clock::now() - clock0).count()
         << " seconds" << endl;

    // Dump the index and distance for the deepestVertex
    cout << "maxDistance to boundary: " << boundDist[deepestVertex]
         << " at vertex: " << deepestVertex << endl;

    // Dump distances to boundary (for small meshes) for debugging
    if (boundDist.size() < 40) {
      cout << "Distances to boundary: " << endl;
      for (size_t i = 0; i < boundDist.size(); i++)
        cout << "vertex: " << i << " : " << boundDist[i] << endl;
      cout << endl;
    }

    /* Dump path from deepestVertex to nearest boundary vertex
    cout << "Shortest path from " << deepestVertex<< " to boundary:" << endl;
    int i = deepestVertex;
    while (i != -1)
    {
        cout << i << " -> ";
        i = parent[i];
    }
    cout << endl;
    //Dump path */

    std::string outputFilename = "output_boundary.ply";
    cout << "Saving output to " << outputFilename << endl;

    // TODO 2.3:
    // Create a color mesh where the color of each vertex shows its distance
    // to boundary Save the color mesh to a PLY file named
    // "output_boundary.ply" See meshColor.cpp or meshColor2.cpp to see an
    // how-to example

    ColorMesh outputMesh;
    outputMesh.coordinates = mesh.coordinates;
    outputMesh.triangles = mesh.triangles;
    outputMesh.colors.resize(mesh.numVertex());

    // Iterate over the vertex of the mesh and set the color
    // of each vertex using its boundDistance
    for (size_t i = 0; i < mesh.numVertex(); i++) {
      outputMesh.colors[i].setTemperature(float(maxDist - boundDist[i]), 0,
                                          float(maxDist));
    }

    // Save result to a file in .ply format
    outputMesh.writeFilePLY(outputFilename);

    // END TODO 2.3
    cout << "Done colorize by distance to boundary "
         << duration<float>(high_resolution_clock::now() - clock0).count()
         << " seconds" << endl;

    // Visualize the file with an external viewer
#ifdef WIN32
    // string viewcmd = "\"C:\\Program Files
    // (x86)\\VCG\\MeshLab\\meshlab.exe\"";
    string viewcmd = "C:/meshlab/meshlab_32.exe";
#else
    string viewcmd = "/usr/bin/flatpak run --branch=stable --arch=x86_64 "
                     "--command=meshlab "
                     "--file-forwarding net.meshlab.MeshLab >/dev/null 2>&1 ";
#endif
    string cmd = viewcmd + " " + outputFilename;
    cout << "Executing external command: " << cmd << endl;
    return system(cmd.c_str());
  }

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
