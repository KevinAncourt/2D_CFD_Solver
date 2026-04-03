#pragma once

#include "mesh/mesh_reader.hpp"

class MeshCompute {

public:
    MeshCompute(const MeshReader& reader); 
    void compute_center_cells();
    void compute_cells_area();


private:
    const MeshReader& reader_;
    std::vector<double> centercells_X;
    std::vector<double> centercells_Y;
    std::vector<double> centercells_Z;
    std::vector<double> cellsarea;

    


};