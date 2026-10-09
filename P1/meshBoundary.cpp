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
 * Name of alumn: 
 * Email of alumn:
 * Year: 2026
 * 
 */

#ifdef _MSC_VER
#pragma warning(error: 4101)
#endif

#define _CRT_NONSTDC_NO_DEPRECATE
#include <iostream>
#include <cmath>
#include <SimpleMesh.hpp>
#include <ColorMesh.hpp>
#include <chrono>
using namespace std::chrono;

/// Update the contents of externalEdges and internalEdges 
void updateEdgeLists(const SimpleMesh &mesh, 
                     std::vector<SimpleEdge> &externalEdges,
                     std::vector<SimpleEdge> &internalEdges )
{
    //TODO 2.1: Implement the body of the updateEdgeLists() method
    throw ("updateEdgeLists has to be implemented as exercise");
    //END TODO 2.1

}//void updateEdgeLists()


int main (int argc, char *argv[])
{
    try
    {
        // Set default input mesh filename
        std::string filename("mallas/16Triangles.off");
        //std::string filename("mallas/mannequin.ply");
        //std::string filename("mallas/knot-hole.ply");
        //std::string filename("mallas/Nefertiti.990kv.ply");

        if (argc > 1)
            filename = std::string(argv[1]);

        ///////////////////////////////////////////////////////////////////////
        //Read a mesh and write a mesh
        SimpleMesh mesh;        
        cout << "Loading file " << filename << endl;
        mesh.readFile(filename, false);

        cout << "Num vertex: " << mesh.numVertex() << " Num triangles: " << mesh.numTriangles()
             << " Unreferenced vertex: " << mesh.checkUnreferencedVertex() << endl;

        //Init time measures used for profiling
        high_resolution_clock::time_point clock0 = high_resolution_clock::now();

        ///////////////////////////////////////////////////////////////////////

        //Compute the list of external and internal edges
        //Implement the body of the updateEdgeLists() method .
        std::vector<SimpleEdge> externalEdges;
        std::vector<SimpleEdge> internalEdges;
        updateEdgeLists(mesh, externalEdges, internalEdges);

 
        cout << "Done updateEdgeLists() " << duration<float>(high_resolution_clock::now() - clock0).count() << " seconds" << endl;
        cout << externalEdges.size() << " boundary edges" << endl;
        cout << internalEdges.size() << " internal edges" << endl;

        //Write external boundary and compute boundary length
        double boundaryLength = 0.0;
        cout << "Edges in the boundary:" << endl;
        for (const SimpleEdge &e : externalEdges)
        {
            cout << "[" << e.a << "->" << e.b << "] ";
            boundaryLength += mesh.edgeLength(e);
        }
        cout << endl;
        cout << "Boundary length: " << boundaryLength << endl;

        //Compute the Euler Characteristic for this mesh
        //For a conected mesh, the value means:
        //2 -> The mesh is a closed surface with no holes (sphere-like topology)
        //1 -> The mesh has one hole (disk-like topology)
        //0 -> The has one handle (torus-like topology)
        //-2 -> The mesh has two holes (tube-like topology)
        int eulerCharacteristic = mesh.numVertex() + mesh.numTriangles()
                                - externalEdges.size() - internalEdges.size();
        cout << "Euler Characteristic: " << eulerCharacteristic << endl;


        //Vector to store the distance from a vertex to the nearest vertex in the boundary
        std::vector<double> boundDist;
        //Index of the vertex with max distance to boundary
        unsigned deepestVertex = 0xDEADC0DE;

        //TODO 2.2:
        //Compute the distance from each vertex to the nearest vertex in the boundary (euclidean distance measured along edges)
        //and store it in the boundDist vector.
        //Store in deepestVertex the index of the vertex with the maximum distance to boundary

        //END TODO 2.2
        cout << "Done boundDist() " << duration<float>(high_resolution_clock::now() - clock0).count() << " seconds" << endl;

        //Dump the index and distance for the deepestVertex
        cout << "maxDistance to boundary: " << boundDist[deepestVertex] << " at vertex: " << deepestVertex << endl;

        //Dump distances to boundary (for small meshes) for debugging
        if (boundDist.size() < 40)
        {
            cout << "Distances to boundary: " << endl;
            for (size_t i=0; i < boundDist.size(); i++)
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

        std::string outputFilename="output_boundary.ply";
        cout << "Saving output to " << outputFilename << endl;

        //TODO 2.3:
        //Create a color mesh where the color of each vertex shows its distance to boundary
        //Save the color mesh to a PLY file named "output_boundary.ply"
        //See meshColor.cpp or meshColor2.cpp to see an how-to example

        //END TODO 2.3
        cout << "Done colorize by distance to boundary " << duration<float>(high_resolution_clock::now() - clock0).count() << " seconds" << endl;

        //Visualize the file with an external viewer
#ifdef WIN32
        //string viewcmd = "\"C:\\Program Files (x86)\\VCG\\MeshLab\\meshlab.exe\"";
        string viewcmd = "C:/meshlab/meshlab_32.exe";
#else
        string viewcmd = "meshlab >/dev/null 2>&1 ";
#endif
        string cmd = viewcmd+" "+outputFilename;
        cout << "Executing external command: " << cmd << endl;
        return system(cmd.c_str());
    }

    catch (const string &str) { std::cerr << "EXCEPTION: " << str << std::endl; }
    catch (const char *str) { std::cerr << "EXCEPTION: " << str << std::endl; }
    catch (std::exception& e)    { std::cerr << "EXCEPTION: " << e.what() << std::endl;  }
    catch (...) { std::cerr << "EXCEPTION (unknow)" << std::endl; }

#ifdef WIN32
    cout << "Press Return to end the program" <<endl;
    cin.get();
#else
#endif

    return 0;
}

