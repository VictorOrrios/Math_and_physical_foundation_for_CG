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

/// Update the contents of externalEdges and internalEdges
void updateEdgeLists(const SimpleMesh &mesh,
                     std::vector<SimpleEdge> &externalEdges,
                     std::vector<SimpleEdge> &internalEdges) {
  // TODO 2.1: Implement the body of the updateEdgeLists() method

  struct SimpleEdgeHash {
    std::size_t operator()(const SimpleEdge &e) const noexcept {
      const std::uint32_t lo = std::min(e.a, e.b);
      const std::uint32_t hi = std::max(e.a, e.b);

      return static_cast<std::size_t>((static_cast<std::uint64_t>(lo) << 32) |
                                      hi);
    }
  };

  struct SimpleEdgeEqual {
    bool operator()(const SimpleEdge &lhs,
                    const SimpleEdge &rhs) const noexcept {
      return (lhs.a == rhs.a && lhs.b == rhs.b) ||
             (lhs.a == rhs.b && lhs.b == rhs.a);
    }
  };

  // Count how many times each undirected edge appears:
  // 1 -> external edge
  // 2 -> internal edge
  // >2 -> non-manifold mesh
  std::unordered_map<SimpleEdge, unsigned char, SimpleEdgeHash, SimpleEdgeEqual>
      edgeCounts;

  const std::size_t maxEdgeOccurrences = 3 * mesh.triangles.size();
  edgeCounts.reserve(maxEdgeOccurrences);

  uint32_t numInternalEdges = 0;

  for (const auto &tri : mesh.triangles) {
    for (const auto &edge : tri.edges()) {
      auto it = edgeCounts.find(edge);

      if (it == edgeCounts.end()) {
        edgeCounts.emplace(edge, 1);
      } else {
        if (it->second == 2) {
          throw "TODO 2.1 IMPOSIBLE: Input mesh is not a manifold";
        }

        it->second++;
        numInternalEdges++;
      }
    }
  }

  externalEdges.clear();
  internalEdges.clear();

  externalEdges.reserve(maxEdgeOccurrences - numInternalEdges);
  internalEdges.reserve(numInternalEdges);

  // Separate external and internal edges.
  for (const auto &entry : edgeCounts) {
    const SimpleEdge &edge = entry.first;
    const unsigned char count = entry.second;

    if (count == 1) {
      externalEdges.push_back(edge);
    } else {
      internalEdges.push_back(edge);
    }
  }

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

  // Build the boundary edge map.
  std::unordered_map<unsigned, unsigned> edgeMap;
  edgeMap.reserve(externalEdges.size());

  for (const auto &exEdge : externalEdges) {
    auto it = edgeMap.find(exEdge.a);

    if (it != edgeMap.end()) {
      edgeMap.emplace(exEdge.b, exEdge.a);
      updateMin(exEdge.b);
    } else {
      edgeMap.emplace(exEdge.a, exEdge.b);
      updateMin(exEdge.a);
    }
  }

  // Order external edges into boundary loops.
  unsigned currVertex = minVertex;

  while (!edgeMap.empty()) {
    auto it = edgeMap.find(currVertex);

    if (it == edgeMap.end()) {
      // Start a new boundary loop.
      minPos.set(INFINITY, INFINITY, INFINITY);
      minVertex = static_cast<unsigned>(-1);

      for (const auto &entry : edgeMap) {
        updateMin(entry.first);
      }

      if (minVertex == static_cast<unsigned>(-1)) {
        break;
      }

      currVertex = minVertex;
      it = edgeMap.find(currVertex);

      if (it == edgeMap.end()) {
        throw "Invalid boundary edge map";
      }
    }

    const unsigned nextVertex = it->second;

    externalEdges.push_back(SimpleEdge(currVertex, nextVertex));

    edgeMap.erase(it);
    currVertex = nextVertex;
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
    // boundDist vector. Store in deepestVertex the index of the vertex with the
    // maximum distance to boundary

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
    // Create a color mesh where the color of each vertex shows its distance to
    // boundary Save the color mesh to a PLY file named "output_boundary.ply"
    // See meshColor.cpp or meshColor2.cpp to see an how-to example

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
    string viewcmd =
        "/usr/bin/flatpak run --branch=stable --arch=x86_64 --command=meshlab "
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
