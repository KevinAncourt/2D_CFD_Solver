#include "mesh/mesh_compute.hpp"
#include <cmath>


MeshCompute::MeshCompute(const MeshReader& reader)
    : reader_(reader)
{
    std::cout << "MeshCompute created " << std::endl;
}

void MeshCompute::compute_center_cells()
{
    const auto& X = reader_.get_CoordX();
    const auto& Y = reader_.get_CoordY();
    const auto& Z = reader_.get_CoordZ();
    const auto& T = reader_.get_Triangles();

    centercells_X.resize(T.size());
    centercells_Y.resize(T.size());
    centercells_Z.resize(T.size());

    for(size_t i=0; i< T.size(); i++)
    {
        int n0 = T[i][0];
        int n1 = T[i][1];
        int n2 = T[i][2];
        centercells_X[i] =(X[n0] + X[n1] + X[n2])/3.0;
        centercells_Y[i] =(Y[n0] + Y[n1] + Y[n2])/3.0;
        centercells_Z[i] =(Z[n0] + Z[n1] + Z[n2])/3.0;
    }

    std::cout << "Center Cells Calculated " << std::endl;
}

void MeshCompute::compute_cells_area()
{
    const auto& X = reader_.get_CoordX();
    const auto& Y = reader_.get_CoordY();
    const auto& Z = reader_.get_CoordZ();
    const auto& T = reader_.get_Triangles();

    cellsarea.resize(T.size());

    for(size_t i=0; i<T.size(); i++)
    {
        int n0 = T[i][0];
        int n1 = T[i][1];
        int n2 = T[i][2];

        double ux = X[n1] - X[n0];
        double uy = Y[n1] - Y[n0];
        double uz = Z[n1] - Z[n0];

        double vx = X[n2] - X[n0];
        double vy = Y[n2] - Y[n0];
        double vz = Z[n2] - Z[n0];

        double cx = uy * vz - uz * vy;
        double cy = uz * vx - ux * vz;
        double cz = ux * vy - uy * vx;

        cellsarea[i] = 0.5 * std::sqrt(cx*cx + cy*cy + cz*cz);
    }

    std::cout << "Cells area Calculated " << std::endl;
}