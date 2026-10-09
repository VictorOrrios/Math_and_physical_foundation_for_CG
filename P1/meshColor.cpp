/*
 * MeshColor.cpp
 *
 * Written by Jose Miguel Espadero <josemiguel.espadero@urjc.es>
 *
 * This code is written as material for the FMF class of the
 * Master Universitario en Informatica Grafica, Juegos y Realidad Virtual.
 * Its purpose is to be didactic and easy to understand, not hard optimized.
 *
 * Read a triangle mesh and compute a color for each vertex based on the
 * Euclidean distance to the first vertex of the mesh.
 *
 * //TODO: Nothing. Just read the code and learn how to read, modify and
 * save a mesh of type ColorMesh with RGB colors per vertex.
 */

#define _CRT_NONSTDC_NO_DEPRECATE
#include "ColorMesh.hpp"
#include <iostream>

#include <chrono> // import C++11 high_resolution_clock for profiling
using namespace std::chrono;

int main(int argc, char *argv[]) {
  try {
    // Set default input mesh filename
    std::string filename("./mallas/bunny.off");
    if (argc > 1)
      filename = std::string(argv[1]);

    // Set default coordinate index to 0 (first vertex)
    unsigned long baseIndex = 0;
    if (argc > 2)
      baseIndex = strtoul(argv[2], nullptr, 10);

    ////////////////////////////////////////////////////////////////////////
    // Read the input mesh
    SimpleMesh mesh;
    cout << "Loading file " << filename << endl;
    mesh.readFile(filename);

    // Show some basic statistic of the mesh
    cout << "Num vertex: " << mesh.numVertex()
         << " Num triangles: " << mesh.numTriangles()
         << " Unreferenced vertex: " << mesh.checkUnreferencedVertex() << endl;

    ///////////////////////////////////////////////////////////////////////////////
    // Compute the Euclidean distance from each vertex to a given point
    // and colorize the mesh to visually display the distance

    // Initialice profiler time.
    high_resolution_clock::time_point clock = high_resolution_clock::now();

    // Get the base point to measure all distances
    // Use .at() instead of [] because index may not be safe (greater than
    // numvertex)
    const vec3 basePoint = mesh.coordinates.at(baseIndex);
    cout << "Euclidean Distances to vertex " << baseIndex << " : " << basePoint
         << endl;

    // Compute the Euclidean distance from each vertex to basePoint
    // also the compute the maximun of all distances.
    size_t numVertex = mesh.numVertex();
    std::vector<double> distances;
    double maxDistance = 0;
    distances.resize(numVertex);
    for (size_t i = 0; i < mesh.numVertex(); i++) {
      // Compute the distance from basepoint to vertex
      distances[i] = basePoint.distance(mesh.coordinates[i]);

      // Search the max of all distances
      if (distances[i] > maxDistance)
        maxDistance = distances[i];
    }

    // Dump maxDistance
    cout << "Max distance to vertex " << baseIndex << " is " << maxDistance
         << endl;

    // Prepare a copy of the input mesh on a ColorMesh
    ColorMesh outputMesh;
    outputMesh.coordinates = mesh.coordinates;
    outputMesh.triangles = mesh.triangles;
    outputMesh.colors.resize(mesh.numVertex());

    // Iterate over the vertex of the mesh and set the color
    // of each vertex using its distance to basePoint.
    for (size_t i = 0; i < mesh.numVertex(); i++) {
      outputMesh.colors[i].setTemperature(float(distances[i]), 0,
                                          float(maxDistance));
    }

    // Measure time since profile timer was initializated
    cout << "\nDone in "
         << duration<float>(high_resolution_clock::now() - clock).count()
         << " seconds" << endl;

    ///////////////////////////////////////////////////////////////////////
    // Save result to a file in .ply format
    string outputFilename = "output_meshColor.ply";
    cout << "Saving output to " << outputFilename << endl;
    outputMesh.writeFilePLY(outputFilename);

    // Visualize the .ply file with an external viewer
#ifdef WIN32
    string viewcmd = "\"C:/Archivos\ de\ programa/VCG/MeshLab/meshlab.exe\"";
    // string viewcmd = "C:/meshlab/meshlab_32.exe";
#else
    string viewcmd = "meshlab >/dev/null 2>&1 ";
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
