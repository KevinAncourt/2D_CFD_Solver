#include "mesh/mesh_reader.hpp"
#include <iostream>
#include <fstream>


int main()
{
       MeshReader mesh("../mesh_generation/msh/naca0012_farfield150.msh");
    mesh.read();
    return 0;

}