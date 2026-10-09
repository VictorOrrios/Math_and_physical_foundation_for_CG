#include <glm/common.hpp>
#include <glm/ext/vector_float3.hpp>
#include <glm/geometric.hpp>
#include <glm/glm.hpp>
#include <vector>

using namespace std;

struct Vertex {
  glm::vec3 pos;
  glm::vec3 normal = glm::vec3(0);
};

struct Tri {
  uint v0;
  uint v1;
  uint v2;
  glm::vec3 normal = glm::vec3(0);
};

struct Mesh {
  Tri *tris;
  Vertex *vertices;
  uint triCount;
  uint vertexCount;
};

float triArea(const Mesh &m, const Tri &t) {
  glm::vec3 e1 = m.vertices[t.v1].pos - m.vertices[t.v0].pos;
  glm::vec3 e2 = m.vertices[t.v2].pos - m.vertices[t.v0].pos;

  return 0.5f * glm::length(glm::cross(e1, e2));
}

float computeAngle(glm::vec3 &a, glm::vec3 &b, glm::vec3 &c) {
  // From Mallas 2, 12
  glm::vec3 ab = b - a;
  glm::vec3 ac = c - a;
  // Compute |ab|·|ac|
  float m = glm::sqrt((glm::dot(ab, ab)) * glm::dot(ac, ac));
  float cosine = glm::dot(ab, ac) / m;
  float sine = glm::length(glm::cross(ab, ac)) / m;
  return atan2(sine, cosine);
}

glm::vec3 calcCentroide(const Mesh &m) {
  glm::vec3 centroide = glm::vec3(0);
  for (uint i = 0; i < m.vertexCount; i++) {
    const Vertex &v = m.vertices[i];
    centroide += v.pos;
  }
  return centroide / glm::vec3(m.vertexCount);
}

std::vector<glm::vec3> calcTriNormalsFacet(Mesh &m) {
  std::vector<glm::vec3> normals;
  normals.reserve(m.triCount);
  for (uint i = 0; i < m.triCount; i++) {
    Tri &t = m.tris[i];
    t.normal =
        glm::normalize(glm::cross(m.vertices[t.v1].pos - m.vertices[t.v0].pos,
                                  m.vertices[t.v2].pos - m.vertices[t.v0].pos));
    normals.push_back(t.normal);
  }
  return normals;
}

std::vector<glm::vec3> calcVertexNormalsAreaWeighted(Mesh &m) {
  std::vector<glm::vec3> normals;
  normals.resize(m.vertexCount);
  std::vector<float> areas(m.vertexCount, 0.0f);
  for (uint i = 0; i < m.triCount; i++) {
    const Tri &t = m.tris[i];
    const float area = triArea(m, t);
    const glm::vec3 weightedNormal = t.normal * area;
    m.vertices[t.v0].normal += weightedNormal;
    m.vertices[t.v1].normal += weightedNormal;
    m.vertices[t.v2].normal += weightedNormal;
    areas[t.v0] += area;
    areas[t.v1] += area;
    areas[t.v2] += area;
  }
  for (uint i = 0; i < m.vertexCount; i++) {
    Vertex &v = m.vertices[i];
    v.normal /= areas[i];
    normals[i] = v.normal;
  }
  return normals;
}

void addWeighted(uint a, uint b, uint c, const glm::vec3 &normal,
                 Vertex *vertices, std::vector<float> &angles) {
  float angle = computeAngle(vertices[a].pos, vertices[b].pos, vertices[c].pos);
  vertices[a].normal += normal * angle;
  angles[a] += angle;
}

std::vector<glm::vec3> calcVertexNormalsAngleWeighted(Mesh &m) {
  std::vector<glm::vec3> normals;
  normals.resize(m.vertexCount);
  std::vector<float> angles(m.vertexCount, 0.0f);
  for (uint i = 0; i < m.triCount; i++) {
    const Tri &t = m.tris[i];
    addWeighted(t.v0, t.v1, t.v2, t.normal, m.vertices, angles);
    addWeighted(t.v1, t.v2, t.v0, t.normal, m.vertices, angles);
    addWeighted(t.v2, t.v0, t.v1, t.normal, m.vertices, angles);
  }
  for (uint i = 0; i < m.vertexCount; i++) {
    Vertex &v = m.vertices[i];
    v.normal /= angles[i];
    normals[i] = v.normal;
  }
  return normals;
}

int main() { return 0; }
