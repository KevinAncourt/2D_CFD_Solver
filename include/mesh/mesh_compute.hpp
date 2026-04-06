#pragma once

#include "mesh/mesh_reader.hpp"

class MeshCompute {

public:
    MeshCompute(const MeshReader& reader); 
    void compute_center_cells();
    void compute_cells_area();
    void compute_faces();
    void orient_faces();



private:
    const MeshReader& reader_;
    std::vector<double> centercells_X;
    std::vector<double> centercells_Y;
    std::vector<double> centercells_Z;
    std::vector<double> cellsarea;
    std::vector<std::array<int,2>> faces_nodes;   // Faces nodes
    std::vector<std::array<int,2>> faces_cells;   // Left and right cells 
    std::vector<double> faces_center_x;
    std::vector<double> faces_center_y;
    std::vector<double> faces_normal_x;
    std::vector<double> faces_normal_y;
    std::vector<double> faces_length;

    


};