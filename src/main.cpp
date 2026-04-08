#include "mesh/mesh_reader.hpp"
#include "mesh/mesh_compute.hpp"
#include "physics/euler_physics.hpp"
#include "solver/solver.hpp"
#include <cmath>
#include <iostream>
#include <chrono>

int main()
{
    MeshReader mesh("../mesh_generation/msh/naca0012_farfield150_coarse.msh");
    mesh.read();

    MeshCompute meshcompute(mesh);
    mesh.write_vtk("../output/mesh.vtk");
    mesh.write_vtk_wall("../output/wall.vtk");
    mesh.write_vtk_farfield("../output/farfield.vtk");

    meshcompute.compute_center_cells();
    meshcompute.compute_cells_area();
    meshcompute.compute_faces();
    meshcompute.orient_faces();
    meshcompute.compute_faces_bc_type();

    EulerPhysics Euler;
    FlowSolver solver(meshcompute, Euler);

    // Adimensional freestream state
    double gamma = 1.4;
    double mach  = 0.85;
    double alpha_deg = 2.0;
    double alpha = alpha_deg * M_PI / 180.0;
    double rho_inf = 1.0;
    double u_inf   = std::cos(alpha);
    double v_inf   = std::sin(alpha);
    double p_inf   = 1.0 / (gamma * mach * mach);

    auto t0 = std::chrono::high_resolution_clock::now();
    solver.initialize_freestream(rho_inf, u_inf, v_inf, p_inf);
    solver.run_explicit(500000, 0.5, 1e-15);
    // solver.run_rk4(5000, 0.3, 1e-10);
    solver.write_residual_history("../output/residual_history.dat");
    solver.write_solution_vtk("../output/solution.vtk");
    std::cout << "Explicit run finished." << std::endl;

    auto t1 = std::chrono::high_resolution_clock::now();
    double elapsed = std::chrono::duration<double>(t1 - t0).count();
    std::cout << "Elapsed time = " << elapsed << " s" << std::endl;

    return 0;
}
