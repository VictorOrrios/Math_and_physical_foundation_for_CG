/*
 * meshTexture-Dense.cpp
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

//This file is the simplest solution to the meshTexture exercise.
//Use a dense matrix and does not reduce the problem, so take long time to execute
//for files as mannequin2.ply

#define _CRT_NONSTDC_NO_DEPRECATE
#include <iostream>
#include <iomanip>
#include <chrono>
using namespace std::chrono;

#include <SimpleMesh.hpp>
#include <TextureMesh.hpp>

//Check if Eigen  (a standar Matrix Library) is included. If not,
//you can get it at http://eigen.tuxfamily.org

// <Eigen/Dense> is the module for dense (traditional) matrix and vector.
// You can get a quick reference for using Eigen dense objects at:
// http://eigen.tuxfamily.org/dox/group__QuickRefPage.html
#include <Eigen/Core>
#include <Eigen/LU>

/// Update the contents of externalEdges and internalEdges 
void updateEdgeLists(const SimpleMesh &mesh, 
                     std::vector<SimpleEdge> &externalEdges,
                     std::vector<SimpleEdge> &internalEdges )
{
    //TODO 2.1: Copy the body of the updateEdgeLists() method from meshBoundary
    throw ("updateEdgeLists has to be implemented as exercise");
    //END TODO 2.1

}//void updateEdgeLists()


/// Write an Eigen Matrix to a matlab file
void exportDenseToMatlab (const Eigen::MatrixXd &m, const std::string &filename, const std::string matrixName ="A")
{
    cout << "Export matrix "<< matrixName << " to file: " << filename << std::endl;

    //Open the file as a stream
    ofstream os(filename.c_str());
    if (!os.is_open())
        throw filename + string(": Error creating the file");


    os << "# name: "<<matrixName << std::endl
       << "# type: matrix" << std::endl
       << "# rows: " << m.rows() << std::endl
       << "# columns: " << m.cols() << std::endl;
    Eigen::IOFormat fmt(Eigen::StreamPrecision, Eigen::DontAlignCols, " ", "\n", "", "", "", "\n");
    os << m.format(fmt);

    os.close();
    std::cout << "To import "<< matrixName <<" into matlab use the command: load(\""<<filename<<"\")"<<std::endl;

}//void exportDenseToMatlab (const &MatrixXd m, const std::string &filename)


int main (int argc, char *argv[])
{
    //Set dumpMatrix to true for debugging. Writting matrix to file can take some time
    bool dumpMatrix = false;

    try
    {
        // Set default input mesh filename
        //std::string filename("mallas/16Triangles.off");  //Minimal case test
        std::string filename("mallas/mask2.ply");      //Easy case test
        //std::string filename("mallas/mannequin2.ply"); //Medium case test
        //std::string filename("mallas/laurana50k.ply"); //Really hard for dense matrix
        if (argc > 1)
            filename = std::string(argv[1]);

        ///////////////////////////////////////////////////////////////////////
        //Step 1.
        //Read an input mesh
        SimpleMesh mesh;
        cout << "Loading file " << filename << endl;
        mesh.readFile(filename, false);

        cout << "Num vertex: " << mesh.numVertex() << " Num triangles: " << mesh.numTriangles()
             << " Unreferenced vertex: " << mesh.checkUnreferencedVertex() << endl;

        //Set dumpMatrix to true for debugging. Writting matrix to file can take some time
        dumpMatrix = dumpMatrix || (mesh.numVertex() < 40);

        //Time measure
        high_resolution_clock::time_point clock0 = high_resolution_clock::now();

        ///////////////////////////////////////////////////////////////////////
        //Step 2.
        //Compute edge list and show it. This step should work if you correctly
        //finished the meshBoundary exercise.
        std::vector<SimpleEdge> externalEdges;
        std::vector<SimpleEdge> internalEdges;
        updateEdgeLists(mesh, externalEdges, internalEdges);

        //Count num of internal and external vertex
        size_t numVertex = mesh.coordinates.size();
        size_t numOfExternalVertex = externalEdges.size();
        //size_t numOfInternalVertex = numVertex - numOfExternalVertex;


        //Dump external edges
        cout << numOfExternalVertex << " vertex in external boundary: " << endl;
        if (numOfExternalVertex < 80)
        {
            for (auto &e : externalEdges)
            {
                cout << "[" << e.a << "->" << e.b << "] ";
            }
            cout << endl;
        }

        //Build the list of external vertex
        std::vector<size_t>externalVertex;
        externalVertex.reserve(externalEdges.size());
        for (auto &e : externalEdges)
        {
            externalVertex.push_back(e.a);
        }

        //Optimization: Keep a lookup table to ask if one vertex is external
        std::vector<bool>isExternal(numVertex, false);
        for(size_t i=0; i< numOfExternalVertex; i++)
        {
            isExternal[externalVertex[i]] = true;
        }

        cout << "Done compute Boundary. " << duration<float>(high_resolution_clock::now() - clock0).count() << " seconds" << endl;

        ///////////////////////////////////////////////////////////////////////
        //Step 3.
        //Build Laplace matrix (step 1)
        cout << "Computing Laplacian matrix (dense):" << endl;

        //We will use a dense matrix. If you want to use sparse matrix
        //you have to check http://eigen.tuxfamily.org/dox/group__TutorialSparse.html
        //and declare a sparse matrix instead.
        Eigen::MatrixXd meshMatrix(numVertex, numVertex);

        //TODO 3.1: Build the Laplacian matrix using the weights
        //
        //      /\         .  Lij = cot(α) + cot(β) , if vertex i neighbour of j
        //     /β \        .
        //    /    \       .  Lii = -Sum Lij , diagonal element is the sum of row
        //   /      \      .
        // vi--------vj    .
        //   \      /      .
        //    \    /       .
        //     \α /        .
        //      \/         .
        //

        //END TODO 3.1

        cout << "Done Laplacian matrix (dense). " << duration<float>(high_resolution_clock::now() - clock0).count() << " seconds" << endl;

        if (dumpMatrix)
            exportDenseToMatlab(meshMatrix, "laplacian.mat", "L");

        ///////////////////////////////////////////////////////////////////////
        //Step 5.1
        //Build system matrix
        cout << "Computing System matrix:" << endl;

        //TODO 3.2: Patch Laplace matrix to generate a valid system of equations

        //END TODO 3.2

        cout << "Done System matrix. " << duration<float>(high_resolution_clock::now() - clock0).count() << " seconds" << endl;

        if (dumpMatrix)
            exportDenseToMatlab(meshMatrix, "systemMatrix.mat", "A");

        ///////////////////////////////////////////////////////////////////////
        //Step 5.2
        //Build UV values for external vertex, using the boundary of a square.
        cout << "Computing contour conditions:" << endl;

        //Note that this is a matrix of numvertex rows and 2 columns 
        Eigen::MatrixX2d UV_0(numVertex,2);
        UV_0.setZero();
        
        //TODO 3.3: Build a set of valid UV values for external vertex, mapping
        //the vertex to the boundary of a square of side unit.

        //END TODO 3.3

        cout << "Done contour conditions. " << duration<float>(high_resolution_clock::now() - clock0).count() << " seconds" << endl;

        if (dumpMatrix)
            exportDenseToMatlab(UV_0, "boundaryUVO.mat", "UV0");

        ///////////////////////////////////////////////////////////////////////
        //Step 6

        //We will use Eigen to solve the matrix from here. 
        cout << "Solving the system using Eigen (dense):" << endl;

        //Solve Ax = b; where A = meshMatrix and b = UV_0
        //This code is valid only for solving a dense matrix. If you want to use sparse matrix
        //you have to check http://eigen.tuxfamily.org/dox/group__TopicSparseSystems.html
        Eigen::PartialPivLU<Eigen::MatrixXd>meshMatrixLU = meshMatrix.lu();
        Eigen::MatrixX2d UV = meshMatrixLU.solve(UV_0);

        cout << "Done solving the system. " << duration<float>(high_resolution_clock::now() - clock0).count() << " seconds" << endl;

        //Dump the computed solution to the system
        if (dumpMatrix)
            exportDenseToMatlab(UV, "solutionUV.mat", "UV");

        ///////////////////////////////////////////////////////////////////////
        //Step 6.4 

        //Create a planar mesh (reverse UV mesh)
        SimpleMesh planarMesh;
        //Use UV values as geometrical coordinates and same triangles that input mesh
        planarMesh.coordinates.resize(numVertex);
        for (size_t i=0; i< numVertex; i++)
        {
            planarMesh.coordinates[i].set(UV(i,0), UV(i,1), 0);
        }
        planarMesh.triangles = mesh.triangles;

        string output_UVMesh1="output_UVMesh1.ply";
        cout << "Saving parameterization mesh to " << output_UVMesh1 << endl;
        //planarMesh.writeFileOBJ(output_UVMesh1);
        planarMesh.writeFilePLY(output_UVMesh1);


        //Create a TextureMesh mesh with UV coordinates
        TextureMesh textureMesh;
        textureMesh.coordinates = mesh.coordinates;
        textureMesh.triangles= mesh.triangles;

        //Set image filename to be used as texture
        textureMesh.textureFile= "UVchecker.jpg";

        //Set UV as texture-per-vertex coordinates
        textureMesh.UV.resize(numVertex);
        for (size_t i=0; i< numVertex; i++)
            textureMesh.UV[i].set(float(UV(i,0)), float(UV(i,1)));

        //Dump textureMesh to file (.obj or .ply)
        string output_UVMesh2="output_UVMesh2.ply";
        cout << "Saving texture mesh to " << output_UVMesh2 << endl;
        //textureMesh.writeFileOBJ(output_UVMesh2);
        textureMesh.writeFilePLY(output_UVMesh2);

        //Visualize the file with an external viewer
#ifdef WIN32
        //string viewcmd = "\"C:\\Program Files (x86)\\VCG\\MeshLab\\meshlab.exe\"";
        string viewcmd = "C:/meshlab/meshlab_32.exe";
#else
        string viewcmd = "meshlab >/dev/null 2>&1 ";
#endif
        string cmd = viewcmd+" "+output_UVMesh2;
        cout << "Executing external command: " << cmd << endl;
        return system(cmd.c_str());

    }//try
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

