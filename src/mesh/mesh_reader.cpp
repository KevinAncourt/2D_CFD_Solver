#include "mesh/mesh_reader.hpp"
#include <iostream>
#include <fstream>
#include <sstream>



MeshReader::MeshReader(const std::string& filename)
{
    filename_ = filename;
}

void MeshReader::read_node(std::ifstream& file)
{
    std::string line;
    while (std::getline(file, line))
    {
        if (line == "$Nodes")
        {
            std::getline(file, line);

            std::istringstream iss(line);

            int numEntityBlocks;
            int numNodes;
            int minNodeTag;
            int maxNodeTag;

            iss >> numEntityBlocks >> numNodes >> minNodeTag >> maxNodeTag;

            std::cout << "numEntityBlocks = " << numEntityBlocks << std::endl;
            std::cout << "numNodes        = " << numNodes << std::endl;
            std::cout << "minNodeTag      = " << minNodeTag << std::endl;
            std::cout << "maxNodeTag      = " << maxNodeTag << std::endl;

            CoordX.resize(numNodes);
            CoordY.resize(numNodes);
            CoordZ.resize(numNodes);
            lntogn.resize(numNodes);
            gntoln.resize(maxNodeTag + 1, -1);

            int localIndex = 0;

            for (int b = 0; b < numEntityBlocks; b++)
            {
                std::getline(file, line);
                std::istringstream blockHeader(line); // Convert the line into an input stream to extract values easily
 
                int entityDim;
                int entityTag;
                int parametric;
                int numNodesInBlock;

                blockHeader >> entityDim >> entityTag >> parametric >> numNodesInBlock;

                std::vector<int> nodeTags(numNodesInBlock);

                for (int i = 0; i < numNodesInBlock; i++)
                {
                    std::getline(file, line);
                    nodeTags[i] = std::stoi(line); // Convert the string line into an integer
                }

                for (int i = 0; i < numNodesInBlock; i++)
                {
                    std::getline(file, line);
                    std::istringstream coordStream(line);

                    double x, y, z;
                    coordStream >> x >> y >> z;

                    int nodeTag = nodeTags[i];

                    CoordX[localIndex] = x;
                    CoordY[localIndex] = y;
                    CoordZ[localIndex] = z;

                    lntogn[localIndex] = nodeTag;
                    gntoln[nodeTag] = localIndex;

                    localIndex++;
                }
            }

            std::getline(file, line); // lit $EndNodes

            std::cout << "Nodes finished" << std::endl;
                        
        }         
    }
}


void MeshReader::read()
{
    std::ifstream file(filename_);  

    if (!file)
    {
        std::cerr << "Error : Can't open Mesh " << filename_ << std::endl;
        return;
    }

    std::cout << "Mesh Open : " << filename_ << std::endl;
    MeshReader::read_node(file);

    for (int i = 0; i < 10; i++)
    {
        std::cout << "local = " << i
              << " gmsh = " << lntogn[i]
              << " x = " << CoordX[i]
              << " y = " << CoordY[i]
              << " z = " << CoordZ[i]
              << std::endl;
    }  


    
}
