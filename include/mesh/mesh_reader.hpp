#pragma once

#include <string>
#include <vector>
#include <array>

class MeshReader
{
public:
    MeshReader(const std::string& filename);
    void read_node(std::ifstream& file);
    void read_element(std::ifstream& file);
    void read();


private : 
    std::string filename_;
    std::vector<double> CoordX;
    std::vector<double> CoordY;
    std::vector<double> CoordZ;
    std::vector<int> lntogn;
    std::vector<int> gntoln;
    std::vector<std::array<int, 2>> edgebcwall;
    std::vector<std::array<int, 2>> edgebcfarfield;
    std::vector<std::array<int, 3>> triangles;
};