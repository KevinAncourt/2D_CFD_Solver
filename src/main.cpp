#include "mesh/mesh_reader.hpp"
#include "mesh/mesh_compute.hpp"



int main()
{
    MeshReader mesh("../mesh_generation/msh/naca0012_farfield150.msh");
    mesh.read();
    MeshCompute meshcompute(mesh);
    mesh.write_vtk("../output/mesh.vtk");
    mesh.write_vtk_wall("../output/wall.vtk");
    mesh.write_vtk_farfield("../output/farfield.vtk");
    meshcompute.compute_center_cells();
    meshcompute.compute_cells_area();
    return 0;

}