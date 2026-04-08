#pragma once

#include <array>
#include <vector>
#include "mesh/mesh_compute.hpp"
#include "physics/euler_physics.hpp"
#include "bc/boundary_conditions.hpp"
#include "flux/rusanov_flux.hpp"
#include <string>

class FlowSolver
{
public:
    FlowSolver(const MeshCompute& mesh,
               const EulerPhysics& physics);

    void initialize_solution(const std::array<double,4>& W0);
    void initialize_freestream(double rho_inf,
                               double u_inf,
                               double v_inf,
                               double p_inf);

    void reset_residual();
    void compute_residual();
    void compute_local_dt(double cfl);
    void update_solution_explicit();
    void set_uniform_dt(double dt_value);
    void run_explicit(int niter, double cfl, double tol);

    double compute_residual_norm() const;
    double compute_residual_max_norm() const;
    double compute_residual_l2_raw() const;

    int get_ncells() const;
    const std::vector<std::array<double,4>>& get_W() const;
    const std::vector<std::array<double,4>>& get_residual() const;
    const std::vector<double>& get_dt() const;

    void write_solution_vtk(const std::string& filename) const;
    void write_residual_history(const std::string& filename) const;
    void run_rk4(int niter, double cfl, double tol);

private:
    const MeshCompute& mesh_;
    const EulerPhysics& physics_;
    BoundaryConditions boundary_conditions_;
    RusanovFlux rusanov_flux_;

    int ncells_;
    std::vector<std::array<double,4>> W_;
    std::vector<std::array<double,4>> residual_;
    std::vector<double> dt_;
    std::array<double,4> W_inf_;
    std::vector<double> residual_history_;
    double initial_residual_ = -1.0;
};
