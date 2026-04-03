#include "mesh/mesh_reader.hpp"
#include <iostream>
#include <fstream>


int main()
{
    MeshReader mesh("../mesh_generation/msh/naca0012_farfield150.msh");
    mesh.read();
    mesh.write_vtk("../output/mesh.vtk");
    mesh.write_vtk_wall("../output/wall.vtk");
    mesh.write_vtk_farfield("../output/farfield.vtk");

    return 0;

}