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
#include <algorithm>
#include <array>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <queue>
#include <sys/types.h>
#include <unordered_map>
#include <unordered_set>
#include <vector>
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

#include <algorithm>
#include <atomic>
#include <omp.h>

constexpr float DELTA = 0.3f;
constexpr uint MAX_VTX_VALENCE = 20;

struct WeightedEdge {
  uint32_t index;
  float weight;
};
struct VertexAdj {
  uint8_t valence = 0;
  std::array<WeightedEdge, MAX_VTX_VALENCE> edges;
};

// Mesh with the most vertex is Venus3000 with 1.5M vertex
// 22 bits to store an index => 44 bits per edge => 11bit radix sort is optimal
inline uint64_t edge64bit(uint32_t a, uint32_t b) {
  return a < b ? (uint64_t(a) << 22) | b : (uint64_t(b) << 22) | a;
}

inline uint64_t edge64bitWIdxEncoding(uint32_t idx, uint32_t a, uint32_t b) {
  return ((uint64_t(idx) << 44)) | (uint64_t(a) << 22) | b;
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

std::vector<VertexAdj> getAdjMatrix(SimpleMesh &mesh,
                                    std::vector<SimpleEdge> &internalEdges) {

  const unsigned N = mesh.numVertex();
  const unsigned M = internalEdges.size();

  std::vector<VertexAdj> adj(N);

  // 1. Compute edge weights
  std::vector<float> weight(M);
#pragma omp parallel for schedule(static)
  for (size_t i = 0; i < M; i++) {
    weight[i] = mesh.distance(internalEdges[i].a, internalEdges[i].b);
  }

  // 2. Build adjacency matrix

  std::vector<uint64_t> inverseEdge(M);

  // 2.a. First pass: Direct fill wirh internal edges
#pragma omp parallel shared(internalEdges, adj)
  {
    const int tid = omp_get_thread_num();
    const int numThreads = omp_get_num_threads();

    // Start and finish index, thread secuential
    const std::size_t begin = M * tid / numThreads;
    const std::size_t end = M * (tid + 1) / numThreads;

    {
      // 2.a.1 Skip start the first group
      std::size_t skipStart = begin;
      if (skipStart > 0 &&
          internalEdges[skipStart].a == internalEdges[skipStart - 1].a) {
        const uint64_t a = internalEdges[skipStart].a;
        while (skipStart < end && internalEdges[skipStart].a == a)
          skipStart++;
      }

      std::size_t it = skipStart;
      while (it < end) {
        // 2.a.2. Find group size
        const uint64_t a = internalEdges[it].a;
        std::size_t gEnd = it + 1;
        while (gEnd < M && internalEdges[gEnd].a == a)
          gEnd++;

        // 2.a.3 Fill direct edge conexxions
        for (size_t k = it; k < gEnd; ++k) {
          unsigned b = internalEdges[k].b;
          // Update weight
          assert(adj[a].valence < MAX_VTX_VALENCE &&
                 "======== Mesh execeds MAX_VTX_VALENCE ========");
          adj[a].edges[adj[a].valence] = {.index = uint32_t(b),
                                          .weight = weight[k]};
          // Cretae 64bit inverse edge with adj index encoding
          inverseEdge[k] = edge64bitWIdxEncoding(adj[a].valence, b, a);
          // Update valence on both vertex
          adj[a].valence++;
        }

        it = gEnd;
      }
    }
  }

  // 2.b. Radix sort the inverse edges
  parallelRadixSort11bit(inverseEdge, 22, 2);

// 2.c. Second pass: Indirect fill with inverse internal edges
#pragma omp parallel shared(inverseEdge, adj)
  {
    constexpr uint64_t BIT_MASK_22 = 0x3FFFFF;

    const int tid = omp_get_thread_num();
    const int numThreads = omp_get_num_threads();

    // Start and finish index, thread secuential
    const std::size_t begin = M * tid / numThreads;
    const std::size_t end = M * (tid + 1) / numThreads;

    {
      // 2.c.1 Skip start the first group
      std::size_t skipStart = begin;
      if (skipStart > 0 &&
          ((inverseEdge[skipStart] >> 22) & BIT_MASK_22) ==
              ((inverseEdge[skipStart - 1] >> 22) & BIT_MASK_22)) {
        const uint64_t a = ((inverseEdge[skipStart] >> 22) & BIT_MASK_22);
        while (skipStart < end &&
               ((inverseEdge[skipStart] >> 22) & BIT_MASK_22) == a)
          skipStart++;
      }

      std::size_t it = skipStart;
      while (it < end) {
        // 2.c.2. Find group size
        const uint64_t a = (inverseEdge[it] >> 22) & BIT_MASK_22;
        std::size_t gEnd = it + 1;
        while (gEnd < M && ((inverseEdge[gEnd] >> 22) & BIT_MASK_22) == a)
          gEnd++;

        // 2.c.3 Fill direct edge conexions, remove encoded index and add
        // light edge counter
        for (size_t k = it; k < gEnd; k++) {
          uint32_t index = uint32_t(inverseEdge[k] >> 44);
          uint32_t b = uint32_t(inverseEdge[k] & 0x3FFFFF);
          float weight = adj[b].edges[index].weight;
          assert(adj[a].valence < MAX_VTX_VALENCE &&
                 "======== Mesh execeds MAX_VTX_VALENCE ========");
          adj[a].edges[adj[a].valence++] = {.index = b, .weight = weight};
        }

        it = gEnd;
      }
    }
  }

  return adj;
}

// Compare and exchange macro
inline bool cas(std::atomic<float> &x, float old_val, float new_val) {
  return x.compare_exchange_strong(old_val, new_val);
}

// Adapted from sbeamer/gapbs/sssp.c
inline void RelaxEdges(const std::vector<VertexAdj> &adj, uint32_t u,
                       std::vector<std::atomic<float>> &dist,
                       std::vector<std::vector<uint32_t>> &local_bins) {
  const auto &a = adj[u];
  const float dist_u = dist[u];
  for (int i = 0; i < a.valence; ++i) {
    const auto &c = a.edges[i];
    const uint32_t v = c.index;
    float old_dist = dist[v];
    float new_dist = dist_u + c.weight;
    while (new_dist < old_dist) {
      if (cas(dist[v], old_dist, new_dist)) {
        size_t dest_bin = static_cast<size_t>(new_dist / DELTA);
        if (dest_bin >= local_bins.size())
          local_bins.resize(dest_bin + 1);
        local_bins[dest_bin].push_back(v);
        break;
      }
      old_dist = dist[v]; // CAS fail, retry
    }
  }
}

// Adapted from sbeamer/gapbs/sssp.c
void delta_stepping_mssp(const std::vector<VertexAdj> &adj,
                         const std::vector<SimpleEdge> &externalEdges,
                         unsigned N,
                         std::vector<std::atomic<float>> &boundDist) {
  const size_t kMaxBin = std::numeric_limits<size_t>::max() / 2;
  const size_t kBinSizeThreshold = 1000; // Gapbs default value

  // 1. Get total mesh valence for maximum frontier size
  // Uses reduction(+ : E) to avoid race
  size_t tot_val = 0;
#pragma omp parallel for schedule(static) reduction(+ : tot_val)
  for (std::size_t i = 0; i < adj.size(); ++i)
    tot_val += static_cast<size_t>(adj[i].valence);

  std::vector<uint32_t> frontier(tot_val, 0);

  // 2. Multisource initialization
  size_t init_tail = 0;
  for (const auto &e : externalEdges) {
    const uint32_t v = e.a;
    if (boundDist[v] > 0.0f) {
      boundDist[v] = 0.0f;
      frontier[init_tail++] = v;
    }
  }

  // Two element arrays for double buffering curr=iter&1, next=(iter+1)&1
  size_t shared_indexes[2] = {0, kMaxBin};
  size_t frontier_tails[2] = {init_tail, 0};

  // 3. Delta stepping
#pragma omp parallel
  {
    std::vector<std::vector<uint32_t>> local_bins(0);
    size_t iter = 0;

    // While there is still nodes in the bins to process
    while (shared_indexes[iter & 1] != kMaxBin) {
      // Ping pong
      size_t &curr_bin_index = shared_indexes[iter & 1];
      size_t &next_bin_index = shared_indexes[(iter + 1) & 1];
      size_t &curr_frontier_tail = frontier_tails[iter & 1];
      size_t &next_frontier_tail = frontier_tails[(iter + 1) & 1];

      // Phase A: Process current frontier
#pragma omp for nowait schedule(dynamic, 64)
      for (size_t i = 0; i < curr_frontier_tail; ++i) {
        const uint32_t u = frontier[i];
        // Lazy filter
        if (boundDist[u] >= DELTA * static_cast<float>(curr_bin_index))
          RelaxEdges(adj, u, boundDist, local_bins);
      }

      // Phase B: Process local small bins
      while (curr_bin_index < local_bins.size() &&
             !local_bins[curr_bin_index].empty() &&
             local_bins[curr_bin_index].size() < kBinSizeThreshold) {
        auto curr_bin_copy = std::move(local_bins[curr_bin_index]);
        local_bins[curr_bin_index].resize(0);
        for (uint32_t u : curr_bin_copy)
          RelaxEdges(adj, u, boundDist, local_bins);
      }

      // Phase C: Get next min non empty bin
      for (size_t i = curr_bin_index; i < local_bins.size(); ++i) {
        if (!local_bins[i].empty()) {
#pragma omp critical
          next_bin_index = std::min(next_bin_index, i);
          break;
        }
      }

#pragma omp barrier

      // Phase D: Reset current bin
#pragma omp single nowait
      {
        curr_bin_index = kMaxBin;
        curr_frontier_tail = 0;
      }

      // Phase E: Promote next bin to gobal frontier
      if (next_bin_index < local_bins.size()) {
        const size_t sz = local_bins[next_bin_index].size();
        const size_t copy_start = __sync_fetch_and_add(&next_frontier_tail, sz);
        std::copy(local_bins[next_bin_index].begin(),
                  local_bins[next_bin_index].end(),
                  frontier.begin() + copy_start);
        local_bins[next_bin_index].clear();
      }

      ++iter;
#pragma omp barrier
    }
  }
}

float findMax(unsigned int &deepestVertex,
              const std::vector<double> &boundDist) {
  float maxDist = -1.0f;
  size_t deepestIndex = 0;

#pragma omp parallel for schedule(static) reduction(max : maxDist, deepestIndex)
  for (size_t i = 0; i < boundDist.size(); ++i) {
    const float d = static_cast<float>(boundDist[i]);
    if (d > maxDist) {
      maxDist = d;
      deepestIndex = i;
    }
  }

  deepestVertex = static_cast<unsigned int>(deepestIndex);
  return maxDist;
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
    using Clock = std::chrono::high_resolution_clock;
    using TP = Clock::time_point;

    auto now = []() { return Clock::now(); };
    auto us = [](TP a, TP b) -> double {
      return std::chrono::duration_cast<
                 std::chrono::duration<double, std::micro>>(b - a)
          .count();
    };

    // -----------------------------------------------------------------------
    TP t0 = now();
    const unsigned N = mesh.numVertex();
    boundDist.assign(N, INFINITY);
    TP t1 = now();

    // -----------------------------------------------------------------------
    std::vector<VertexAdj> adj = getAdjMatrix(mesh, internalEdges);
    TP t2 = now();

    // -----------------------------------------------------------------------
    std::vector<std::atomic<float>> atomicBoundDist(boundDist.size());
#pragma omp parallel for schedule(static)
    for (std::size_t i = 0; i < boundDist.size(); ++i)
      atomicBoundDist[i].store(boundDist[i], std::memory_order_relaxed);
    TP t3 = now();

    // -----------------------------------------------------------------------
    delta_stepping_mssp(adj, externalEdges, N, atomicBoundDist);
    TP t4 = now();

// -----------------------------------------------------------------------
#pragma omp parallel for schedule(static)
    for (std::size_t i = 0; i < boundDist.size(); ++i)
      boundDist[i] = atomicBoundDist[i].load(std::memory_order_relaxed);
    TP t5 = now();

    // -----------------------------------------------------------------------
    float maxDist = findMax(deepestVertex, boundDist);
    TP t6 = now();

    // -----------------------------------------------------------------------
    const double t_init = us(t0, t1);
    const double t_adj = us(t1, t2);
    const double t_copyIn = us(t2, t3);
    const double t_sssp = us(t3, t4);
    const double t_copyOut = us(t4, t5);
    const double t_max = us(t5, t6);
    const double t_total =
        t_init + t_adj + t_copyIn + t_sssp + t_copyOut + t_max;

    auto pct = [t_total](double v) { return 100.0 * v / t_total; };

    std::cout << "\n=== Profiling meshBoundary ===\n";
    std::cout << "  N vertex         : " << N << "\n";
    std::cout << "  threads          : " << omp_get_max_threads() << "\n";
    std::cout << "  ------------------------------------------------\n";
    std::cout << "  init dist        : " << t_init << " us  (" << pct(t_init)
              << " %)\n";
    std::cout << "  getAdjMatrix     : " << t_adj << " us  (" << pct(t_adj)
              << " %)\n";
    std::cout << "  copy to atomic   : " << t_copyIn << " us  ("
              << pct(t_copyIn) << " %)\n";
    std::cout << "  delta-stepping   : " << t_sssp << " us  (" << pct(t_sssp)
              << " %)\n";
    std::cout << "  copy from atomic : " << t_copyOut << " us  ("
              << pct(t_copyOut) << " %)\n";
    std::cout << "  findMax          : " << t_max << " us  (" << pct(t_max)
              << " %)\n";
    std::cout << "  ------------------------------------------------\n";
    std::cout << "  TOTAL            : " << t_total << " us\n";
    std::cout << std::flush;

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
