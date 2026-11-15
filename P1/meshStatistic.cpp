/*
 * meshStatistic.cpp
 *
 * Written by Jose Miguel Espadero <josemiguel.espadero@urjc.es>
 *
 * This code is written as material for the FMF class of the
 * Master Universitario en Informatica Grafica, Juegos y Realidad Virtual.
 * Its purpose is to be didactic and easy to understand, not hard optimized.
 *
 * This file compute some statistics about a mesh and dump then to the console.
 * Also produce an output mesh with degenerate triangles remarked.
 *
 * //TODO: Fill-in your name and email
 * Name of alumn: Víctor Orrios Barón
 * Email of alumn: v.orrios.2026@alumnos.urjc.es
 * Year: 2026
 *
 */

#include <cmath>
#include <limits>
#include <unordered_set>
#include <vector>
#ifdef _MSC_VER
#pragma warning(error : 4101)
#endif

#define _CRT_NONSTDC_NO_DEPRECATE
#include "ColorMesh.hpp"
#include "SimpleMesh.hpp"
#include <chrono>
#include <cstdlib>
#include <iomanip>
#include <iostream>
using namespace std::chrono;

#include <omp.h>

int main(int argc, char *argv[]) {
  try {
    // Set default input mesh filename
    std::string filename("mallas/Nefertiti.990kv.ply");
    if (argc > 1)
      filename = std::string(argv[1]);

    // Set default degenerate triangle shapeFactor threshold
    double shapeFactorTh = 1 / (4 * sqrt(3));
    if (argc > 2)
      shapeFactorTh = strtod(argv[2], nullptr);

    ///////////////////////////////////////////////////////////////////////
    // Read a mesh from given filename
    SimpleMesh mesh;
    cout << "Loading file " << filename << endl;
    mesh.readFile(filename);

    cout << "\nNum vertex: " << mesh.numVertex()
         << " Num triangles: " << mesh.numTriangles()
         << " Unreferenced vertex: " << mesh.checkUnreferencedVertex() << endl;

    // Compute minimum, maximum and average area per triangle, and total area
    // for the mesh
    double minArea, maxArea, averageArea = 0, totalArea;

    // Compute minimum and maximum angle for the mesh
    double minAngle, maxAngle;

    // Compute minimum, maximum and average edge length
    double minEdgeLen, maxEdgeLen, averageEdgeLen;

    // Compute minimum, maximum and average shape factor
    double minShapeFactor, maxShapeFactor, averageShapeFactor;

    //////////////////////////////////////////////////////////////////////////////////
    // TODO 1.1:
    //  Compute the values for minArea, maxArea, averageArea, totalArea
    //  Compute the values for minAngle, maxAngle;
    //  Compute the values for minEdgeLen, maxEdgeLen, averageEdgeLen;
    //  Compute the values for minShapeFactor, maxShapeFactor,
    //  averageShapeFactor;

    high_resolution_clock::time_point clock0 = high_resolution_clock::now();

    // Ensure correct initialization
    const std::size_t numTriangles = mesh.numTriangles();
    const std::size_t numVertices = mesh.numVertex();

    minArea = std::numeric_limits<double>::max();
    maxArea = std::numeric_limits<double>::lowest();
    totalArea = 0;

    minAngle = std::numeric_limits<double>::max();
    maxAngle = std::numeric_limits<double>::lowest();

    minEdgeLen = std::numeric_limits<double>::max();
    maxEdgeLen = std::numeric_limits<double>::lowest();
    averageEdgeLen = 0;

    minShapeFactor = std::numeric_limits<double>::max();
    maxShapeFactor = std::numeric_limits<double>::lowest();
    averageShapeFactor = 0;

    std::vector<double> vertexAreas(numVertices, 0.0);
    double sumVertexAreas = 0;
    double minVertexArea = std::numeric_limits<double>::max();
    double maxVertexArea = std::numeric_limits<double>::lowest();

    std::vector<double> triShapesFactors(numTriangles);

    // A whole lot of reduce to not have per thread vars
#pragma omp parallel for schedule(static)                                      \
    reduction(min : minArea, minAngle, minEdgeLen, minShapeFactor)             \
    reduction(max : maxArea, maxAngle, maxEdgeLen, maxShapeFactor)             \
    reduction(+ : totalArea, averageEdgeLen, averageShapeFactor)
    for (std::size_t i = 0; i < numTriangles; i++) {
      const SimpleTriangle &tri = mesh.triangles[i];

      const vec3 &A = mesh.coordinates[tri.a];
      const vec3 &B = mesh.coordinates[tri.b];
      const vec3 &C = mesh.coordinates[tri.c];

      const vec3 AB = A - B;
      const vec3 AC = A - C;
      const vec3 BC = B - C;

      const double lenAB = AB.module();
      const double lenAC = AC.module();
      const double lenBC = BC.module();

      const double Adot = AB.dot(AC);
      const double Bdot = (-AB).dot(BC);
      const double Cdot = (-AC).dot(-BC);

      const double cosA = Adot / (lenAB * lenAC);
      const double cosB = Bdot / (lenAB * lenBC);
      const double cosC = Cdot / (lenAC * lenBC);

      minAngle = min(minAngle, min({1.0 - cosA, 1.0 - cosB, 1.0 - cosC}));
      maxAngle = max(maxAngle, max({-cosA, -cosB, -cosC}));

      minEdgeLen = min({minEdgeLen, lenAB, lenAC, lenBC});
      maxEdgeLen = max({maxEdgeLen, lenAB, lenAC, lenBC});
      averageEdgeLen += lenAB + lenAC + lenBC;

      // Area calc
      const double crossMag = AB.cross(AC).module();
      const double area = 0.5 * crossMag;

      minArea = min(minArea, area);
      maxArea = max(maxArea, area);
      totalArea += area;

      // Shape factor calc
      const double minL = min({lenAB, lenAC, lenBC});
      const double abc = lenAB * lenAC * lenBC;
      const double radius = abc / (2.0 * crossMag);
      const double shapeF =
          (radius == 0.0) ? 0.0 : minL / (radius * std::sqrt(3.0));

      minShapeFactor = min(minShapeFactor, shapeF);
      maxShapeFactor = max(maxShapeFactor, shapeF);
      averageShapeFactor += shapeF;
      triShapesFactors[i] = shapeF;

      double Aarea, Barea, Carea;
      if (Adot < 0.0 || Bdot < 0.0 || Cdot < 0.0) {
        const double triAreax2 = area * 2.0;
        const double triAreax4 = triAreax2 * 2.0;
        Aarea = Adot > 0.0 ? triAreax2 : triAreax4;
        Barea = Bdot > 0.0 ? triAreax2 : triAreax4;
        Carea = Cdot > 0.0 ? triAreax2 : triAreax4;
      } else {
        const double a = lenBC;
        const double b = lenAC;
        const double c = lenAB;
        const double Across = crossMag; // Reuse
        const double Bcross = AB.cross(BC).module();
        const double Ccross = AC.cross(BC).module();
        const double Aterm = a * a * (Adot / Across);
        const double Bterm = b * b * (Bdot / Bcross);
        const double Cterm = c * c * (Cdot / Ccross);
        Aarea = Bterm + Cterm;
        Barea = Aterm + Cterm;
        Carea = Aterm + Bterm;
      }

#pragma omp atomic
      vertexAreas[tri.a] += Aarea;
#pragma omp atomic
      vertexAreas[tri.b] += Barea;
#pragma omp atomic
      vertexAreas[tri.c] += Carea;
    }

#pragma omp parallel for schedule(static) reduction(min : minVertexArea)       \
    reduction(max : maxVertexArea) reduction(+ : sumVertexAreas)
    for (std::size_t v = 0; v < numVertices; ++v) {
      vertexAreas[v] /= 8.0;
      minVertexArea = min(minVertexArea, vertexAreas[v]);
      maxVertexArea = max(maxVertexArea, vertexAreas[v]);
      sumVertexAreas += vertexAreas[v];
    }

    // End average calculation
    averageArea = totalArea / mesh.numTriangles();
    averageShapeFactor /= mesh.numTriangles();
    averageEdgeLen /= mesh.numTriangles() * 3;
    // Angle correction, cheaper to store in cos(angle) format
    minAngle = acos(1 - minAngle);
    maxAngle = acos(-maxAngle);

    // END TODO 1.1

    // Dump statistics to console output
    cout << std::fixed << std::setprecision(4) << "\nArea    min: " << minArea
         << " max: " << maxArea << " average: " << averageArea
         << " total: " << totalArea << "\nAngle   min: " << minAngle
         << " max: " << maxAngle << "\nEdgeLen min: " << minEdgeLen
         << " max: " << maxEdgeLen << " average: " << averageEdgeLen
         << "\nShapeF  min: " << minShapeFactor << " max: " << maxShapeFactor
         << " average: " << averageShapeFactor << endl;

    //////////////////////////////////////////////////////////////////////////////////
    // TODO OPTATIVE 1:
    // Compute Vertex area for each vertex and store in vertexAreas vector
    // Compute minimun and maximun vertex area in minVertexArea, maxVertexArea

    cout << "Done statistic() "
         << duration<float>(high_resolution_clock::now() - clock0).count()
         << " seconds" << endl;

    // END TODO OPTATIVE 1

    // Check values for vertex areas and compute sum of areas:
    if (vertexAreas.size() != mesh.numVertex()) {
      cout << "VertexAreas OPTATIVE PART NOT DONE" << endl;
    } else {
      cout << "VertexA min: " << minVertexArea << " max: " << maxVertexArea
           << " average: " << sumVertexAreas / mesh.numVertex()
           << " total: " << sumVertexAreas << endl;
    }

    //////////////////////////////////////////////////////////////////////////////////
    // TODO 1.2:
    // Write triangles with ShapeFactor (radius / minEdge) greater than
    // shapeFactorTh Use messages formated as: "Triangle nnn has ShapeFactor
    // xxx"

    std::vector<uint> shapeTris;
    shapeTris.reserve(10);

    for (size_t i = 0; i < mesh.numTriangles(); i++) {
      auto shapeF = triShapesFactors[i];
      if (shapeF < shapeFactorTh)
        shapeTris.push_back(i);
    }

    cout << endl
         << shapeTris.size() << " triangles with ShapeFactor < "
         << shapeFactorTh << endl;
    for (auto &t : shapeTris) {
      cout << "Triangle " << t << " has ShapeFactor " << triShapesFactors[t]
           << endl;
    }
    cout << endl;

    // END TODO 1.2

    //////////////////////////////////////////////////////////////////////////////////
    // TODO 1.3:
    // Create a colorMesh where faces with ShapeFactor greater than
    // shapeFactorTh have their vertex colored in red. Save it to file named
    // output_statistic.ply and visualize it with meshlab or another external
    // viewer.

    ColorMesh outputMesh;
    outputMesh.coordinates = mesh.coordinates;
    outputMesh.triangles = mesh.triangles;

    // Set evry color to white
    RGBColor white;
    white.set(1, 1, 1);
    outputMesh.colors.resize(mesh.numVertex(), white);

    for (size_t i = 0; i < mesh.numTriangles(); i++) {
      const SimpleTriangle &tri = mesh.triangles[i];
      auto shapeF = triShapesFactors[i];
      // Recolor if degenerated
      if (shapeF < shapeFactorTh) {
        outputMesh.colors[tri.a].set(1, 0, 0);
        outputMesh.colors[tri.b].set(1, 0, 0);
        outputMesh.colors[tri.c].set(1, 0, 0);
      }
    }

    // Save to output_statistic.ply and open meshlab

    string outputFilename = "output_statistic.ply";
    cout << "Saving output to " << outputFilename << endl;
    outputMesh.writeFilePLY(outputFilename);

    // Visualize the .ply file with an external viewer
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

    // END TODO 1.3
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
