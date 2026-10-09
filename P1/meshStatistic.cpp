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
#pragma warning(error: 4101)
#endif

#define _CRT_NONSTDC_NO_DEPRECATE
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include "SimpleMesh.hpp"
#include "ColorMesh.hpp"

int main (int argc, char *argv[])
{
    try
    {
        //Set default input mesh filename
        std::string filename("mallas/mask2.ply");
        if (argc >1)
            filename = std::string(argv[1]);

        //Set default degenerate triangle shapeFactor threshold
        double shapeFactorTh = 1/(4*sqrt(3));
        if (argc > 2)
            shapeFactorTh = strtod(argv[2], nullptr);

        ///////////////////////////////////////////////////////////////////////
        //Read a mesh from given filename
        SimpleMesh mesh;
        cout << "Loading file " << filename << endl;
        mesh.readFile(filename);

        cout << "\nNum vertex: " << mesh.numVertex() << " Num triangles: " << mesh.numTriangles()
             << " Unreferenced vertex: " << mesh.checkUnreferencedVertex() << endl;


        //Compute minimum, maximum and average area per triangle, and total area for the mesh
        double minArea, maxArea, averageArea=0, totalArea;

        //Compute minimum and maximum angle for the mesh
        double minAngle, maxAngle;

        //Compute minimum, maximum and average edge length
        double minEdgeLen, maxEdgeLen, averageEdgeLen;

        //Compute minimum, maximum and average shape factor
        double minShapeFactor, maxShapeFactor, averageShapeFactor;

        //////////////////////////////////////////////////////////////////////////////////
        //TODO 1.1:
        // Compute the values for minArea, maxArea, averageArea, totalArea
        // Compute the values for minAngle, maxAngle;
        // Compute the values for minEdgeLen, maxEdgeLen, averageEdgeLen;
        // Compute the values for minShapeFactor, maxShapeFactor, averageShapeFactor;

        // Ensure correct initialization
        minArea = std::numeric_limits<double>::max();
        maxArea = std::numeric_limits<double>::min();
        averageArea = 0;
        totalArea = 0;

        minAngle = std::numeric_limits<double>::max();
        maxAngle = std::numeric_limits<double>::min();

        minEdgeLen = std::numeric_limits<double>::max();
        maxEdgeLen = std::numeric_limits<double>::min();
        averageEdgeLen = 0;

        minShapeFactor = std::numeric_limits<double>::max();
        maxShapeFactor = std::numeric_limits<double>::min();
        averageShapeFactor = 0;

        // Used for not visiting the same edge twice
        std::unordered_set<unsigned> vistedEdges;

        // Used for later 'TODO's
        std::vector<double> triShapesFactors(mesh.numTriangles());

        // Note, auxiliary functions are made in lamba form to fit inside the two 'TODO' comments
        for (size_t i = 0; i < mesh.numTriangles(); i++) {
            // Returns angle of vertex a over vertex b to c in format cos(alpha)
            auto angleCos = [](const vec3& a, const vec3& b, const vec3& c){
                vec3 ab = a-b;
                vec3 cb = c-b;
                return ab.dot(cb)/(ab.module()*cb.module());
            };

            // Adds the angle to max and min counters
            auto addAngle = [&](const unsigned a, const unsigned b, const unsigned c){
                double alpha = angleCos(
                    mesh.coordinates.at(a),
                    mesh.coordinates.at(b),
                    mesh.coordinates.at(c));
                minAngle = min(minAngle,1-alpha);
                maxAngle = max(maxAngle,-alpha);
            };

            // Transforms a SimpleEdge into a linear index
            auto edgeHash = [&](const SimpleEdge& edge){
                unsigned minV = min(edge.a,edge.b);
                unsigned maxV = max(edge.a,edge.b);
                return minV+maxV*mesh.numVertex();
            };

            // Adds edge to the max, min and average counters
            auto addEdge = [&](const SimpleEdge& edge){
                auto hash = edgeHash(edge);
                if(vistedEdges.count(hash)) return;
                vistedEdges.insert(hash);
                double eLenght = mesh.distance(edge.a, edge.b);
                minEdgeLen = min(minEdgeLen,eLenght);
                maxEdgeLen = max(maxEdgeLen,eLenght);
                averageEdgeLen += eLenght;
            };

            const SimpleTriangle& tri = mesh.triangles.at(i);
            double area = mesh.triangleArea(tri);
            double shapeF = mesh.triangleShapeFactor(tri);
            triShapesFactors[i] = shapeF;

            minArea = min(minArea,area);
            maxArea = max(maxArea,area);
            totalArea += area;

            minShapeFactor = min(minShapeFactor,shapeF);
            maxShapeFactor = max(maxShapeFactor,shapeF);
            averageShapeFactor += shapeF;

            // Add the three angles of a triangle to the counters
            addAngle(tri.a, tri.b, tri.c);
            addAngle(tri.b, tri.c, tri.a);
            addAngle(tri.c, tri.a, tri.b);

            // Add the three edges of a triangle to the counters
            for(const auto& edge: tri.edges())
                addEdge(edge);

        }

        // End average calculation
        averageArea = totalArea/mesh.numTriangles();
        averageShapeFactor /= mesh.numTriangles();
        averageEdgeLen /= vistedEdges.size();
        // Angle correction, cheaper to store in cos(angle) format
        minAngle = acos(1-minAngle);
        maxAngle = acos(-maxAngle);

        //END TODO 1.1

        // Dump statistics to console output
        cout << std::fixed << std::setprecision(4) <<
                "\nArea    min: " << minArea <<
                " max: " << maxArea  <<
                " average: " << averageArea  <<
                " total: " << totalArea <<
                "\nAngle   min: " << minAngle <<
                " max: " << maxAngle  <<
                "\nEdgeLen min: " << minEdgeLen <<
                " max: " << maxEdgeLen  <<
                " average: " << averageEdgeLen  <<
                "\nShapeF  min: " << minShapeFactor <<
                " max: " << maxShapeFactor  <<
                " average: " << averageShapeFactor  << endl;

        //////////////////////////////////////////////////////////////////////////////////
        //TODO OPTATIVE 1:
        //Compute Vertex area for each vertex and store in vertexAreas vector
        //Compute minimun and maximun vertex area in minVertexArea, maxVertexArea
        std::vector<double>vertexAreas(mesh.numVertex(),0.0);
        double sumVertexAreas = 0;
        double minVertexArea = INFINITY;
        double maxVertexArea = 0;

        for (size_t i = 0; i < mesh.numTriangles(); i++) {
            const SimpleTriangle& tri = mesh.triangles.at(i);

            vec3 A = mesh.coordinates[tri.a];
            vec3 B = mesh.coordinates[tri.b];
            vec3 C = mesh.coordinates[tri.c];

            vec3 AB = A-B; double c = (AB).module();
            vec3 AC = A-C; double b = (AC).module();
            vec3 BC = B-C; double a = (BC).module();

            double Aarea, Barea, Carea;

            double Adot = AB.dot(AC);
            double Bdot = (-AB).dot(BC);
            double Cdot = (-AC).dot(-BC);

            if(Adot < 0.0 || Bdot < 0.0 || Cdot < 0.0){
                // Obtuse triangle
                double triAreax2 = mesh.triangleArea(tri)*2;
                double triAreax4 = triAreax2*2;
                // Has to be times 2 and times 4 to cancel the later division by 8

                Aarea = Adot>0.0 ? triAreax2 : triAreax4;
                Barea = Bdot>0.0 ? triAreax2 : triAreax4;
                Carea = Cdot>0.0 ? triAreax2 : triAreax4;
            }else{
                // Acute triangle
                double Across = AB.cross(AC).module();
                double Bcross = AB.cross(BC).module();
                double Ccross = AC.cross(BC).module();

                double Aterm = a*a*(Adot/Across); // |BC|²cot(Â)
                double Bterm = b*b*(Bdot/Bcross); // |AC|²cot(B)
                double Cterm = c*c*(Cdot/Ccross); // |AB|²cot(Ĉ)

                Aarea = Bterm+Cterm;
                Barea = Aterm+Cterm;
                Carea = Aterm+Bterm;
            }

            vertexAreas[tri.a] += Aarea;
            vertexAreas[tri.b] += Barea;
            vertexAreas[tri.c] += Carea;
        }

        for (double& area: vertexAreas) {
            area /= 8.0;
            minVertexArea = min(minVertexArea,area);
            maxVertexArea = max(maxVertexArea,area);
            sumVertexAreas += area;
        }

        //END TODO OPTATIVE 1

        //Check values for vertex areas and compute sum of areas:
        if (vertexAreas.size() != mesh.numVertex() )
        {
          cout << "VertexAreas OPTATIVE PART NOT DONE" << endl;
        }
        else
        {
          cout << "VertexA min: " << minVertexArea <<
                  " max: " << maxVertexArea  <<
                  " average: " << sumVertexAreas / mesh.numVertex()  <<
                  " total: " << sumVertexAreas << endl;
        }


        //////////////////////////////////////////////////////////////////////////////////
        //TODO 1.2:
        //Write triangles with ShapeFactor (radius / minEdge) greater than shapeFactorTh
        //Use messages formated as: "Triangle nnn has ShapeFactor xxx"

        for (size_t i = 0; i < mesh.numTriangles(); i++) {
            auto shapeF = triShapesFactors[i];
            if(shapeF<shapeFactorTh)
                cout << "Triangle "<<i<<" has ShapeFactor < "<<shapeF<<endl;
        }

        //END TODO 1.2

        //////////////////////////////////////////////////////////////////////////////////
        //TODO 1.3:
        //Create a colorMesh where faces with ShapeFactor greater than shapeFactorTh
        //have their vertex colored in red. Save it to file named output_statistic.ply
        //and visualize it with meshlab or another external viewer.

        ColorMesh outputMesh;
        outputMesh.coordinates = mesh.coordinates;
        outputMesh.triangles = mesh.triangles;

        // Set evry color to white
        RGBColor white; white.set(1,1,1);
        outputMesh.colors.resize(mesh.numVertex(),white);

        for (size_t i = 0; i < mesh.numTriangles(); i++) {
            const SimpleTriangle& tri = mesh.triangles[i];
            auto shapeF = triShapesFactors[i];
            // Recolor if degenerated
            if(shapeF<shapeFactorTh){
                outputMesh.colors[tri.a].set(1,0,0);
                outputMesh.colors[tri.b].set(1,0,0);
                outputMesh.colors[tri.c].set(1,0,0);
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
        string viewcmd = "/usr/bin/flatpak run --branch=stable --arch=x86_64 --command=meshlab --file-forwarding net.meshlab.MeshLab >/dev/null 2>&1 ";
#endif
        string cmd = viewcmd + " " + outputFilename;
        cout << "Executing external command: " << cmd << endl;

        return system(cmd.c_str());

        //END TODO 1.3
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
