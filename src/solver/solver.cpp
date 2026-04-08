#include "solver/solver.hpp"

#include <cmath>
#include <iostream>
#include <limits>
#include <fstream>
#include <omp.h>

FlowSolver::FlowSolver(const MeshCompute& mesh,
                       const EulerPhysics& physics)
    : mesh_(mesh),
      physics_(physics),
      boundary_conditions_(physics),
      rusanov_flux_(physics),
      W_inf_{1.0, 0.0, 0.0, 1.0}
{
    ncells_ = static_cast<int>(mesh_.get_centercells_X().size());

    W_.resize(ncells_);
    residual_.resize(ncells_);
    dt_.resize(ncells_, 0.0);

    reset_residual();
}

void FlowSolver::initialize_solution(const std::array<double,4>& W0)
{
    #pragma omp parallel for
    for (int i = 0; i < ncells_; i++)
    {
        W_[i] = W0;
    }
}

void FlowSolver::initialize_freestream(double rho_inf,
                                       double u_inf,
                                       double v_inf,
                                       double p_inf)
{
    std::array<double,4> V_inf = {rho_inf, u_inf, v_inf, p_inf};
    W_inf_ = physics_.prim_to_cons(V_inf);

    initialize_solution(W_inf_);
    boundary_conditions_.set_farfield_state(W_inf_);
}

void FlowSolver::reset_residual()
{
    #pragma omp parallel for
    for (int i = 0; i < ncells_; i++)
    {
        residual_[i][0] = 0.0;
        residual_[i][1] = 0.0;
        residual_[i][2] = 0.0;
        residual_[i][3] = 0.0;
    }
}

void FlowSolver::compute_residual()
{
    reset_residual();

    const auto& faces_cells    = mesh_.get_faces_cells();
    const auto& faces_normal_x = mesh_.get_faces_normal_x();
    const auto& faces_normal_y = mesh_.get_faces_normal_y();
    const auto& faces_length   = mesh_.get_faces_length();
    const auto& faces_bc_type  = mesh_.get_faces_bc_type();

    const int nfaces = static_cast<int>(faces_cells.size());
    #pragma omp parallel for
    for (int f = 0; f < nfaces; f++)
    {
        int left  = faces_cells[f][0];
        int right = faces_cells[f][1];

        double nx = faces_normal_x[f];
        double ny = faces_normal_y[f];
        double Sf = faces_length[f];

        std::array<double,4> flux = {0.0, 0.0, 0.0, 0.0};

        if (right != -1)
        {
            flux = rusanov_flux_.compute(W_[left], W_[right], nx, ny, Sf);

            for (int k = 0; k < 4; k++)
            {
                #pragma omp atomic
                residual_[left][k]  += flux[k];
                #pragma omp atomic
                residual_[right][k] -= flux[k];
            }
        }
        else
        {
            int bc = faces_bc_type[f];

            if (bc == 1)
            {
                flux = boundary_conditions_.wall_flux(W_[left], nx, ny, Sf);
            }
            else if (bc == 2)
            {
                flux = boundary_conditions_.farfield_flux(W_[left], nx, ny, Sf);
            }
            else
            {
                continue;
            }

            for (int k = 0; k < 4; k++)
            {
                #pragma omp atomic
                residual_[left][k] += flux[k];
            }
        }
    }
}

void FlowSolver::compute_local_dt(double cfl)
{
    const auto& cells_area     = mesh_.get_cells_area();
    const auto& faces_cells    = mesh_.get_faces_cells();
    const auto& faces_normal_x = mesh_.get_faces_normal_x();
    const auto& faces_normal_y = mesh_.get_faces_normal_y();
    const auto& faces_length   = mesh_.get_faces_length();
    const auto& faces_bc_type  = mesh_.get_faces_bc_type();

    std::vector<double> spectral_sum(ncells_, 0.0);

    const int nfaces = static_cast<int>(faces_cells.size());
    #pragma omp parallel for
    for (int f = 0; f < nfaces; f++)
    {
        int left  = faces_cells[f][0];
        int right = faces_cells[f][1];

        double nx = faces_normal_x[f];
        double ny = faces_normal_y[f];
        double Sf = faces_length[f];

        const std::array<double,4> VL = physics_.cons_to_prim(W_[left]);
        double unL = std::abs(VL[1] * nx + VL[2] * ny);
        double aL  = physics_.sound_speed(W_[left]);

        if (right != -1)
        {
            const std::array<double,4> VR = physics_.cons_to_prim(W_[right]);
            double unR = std::abs(VR[1] * nx + VR[2] * ny);
            double aR  = physics_.sound_speed(W_[right]);

            double lambda = std::max(unL + aL, unR + aR) * Sf;
            #pragma omp atomic
            spectral_sum[left]  += lambda;
            #pragma omp atomic
            spectral_sum[right] += lambda;
        }
        else
        {
            int bc = faces_bc_type[f];

            if (bc == 1)
            {
                #pragma omp atomic
                spectral_sum[left] += (unL + aL) * Sf;
            }
            else if (bc == 2)
            {
                const std::array<double,4> Vinf = physics_.cons_to_prim(W_inf_);
                double unInf = std::abs(Vinf[1] * nx + Vinf[2] * ny);
                double aInf  = physics_.sound_speed(W_inf_);

                double lambda = std::max(unL + aL, unInf + aInf) * Sf;
                #pragma omp atomic
                spectral_sum[left] += lambda;
            }
        }
    }
    #pragma omp parallel for
    for (int i = 0; i < ncells_; i++)
    {
        if (spectral_sum[i] > 1.0e-14)
        {
            dt_[i] = cfl * cells_area[i] / spectral_sum[i];
        }
        else
        {
            dt_[i] = std::numeric_limits<double>::max();
        }
    }
}

void FlowSolver::update_solution_explicit()
{
    const auto& cells_area = mesh_.get_cells_area();
    #pragma omp parallel for
    for (int i = 0; i < ncells_; i++)
    {
        double cell_area = cells_area[i];

        for (int k = 0; k < 4; k++)
        {
            W_[i][k] -= (dt_[i] / cell_area) * residual_[i][k];
        }
    }
}

void FlowSolver::set_uniform_dt(double dt_value)
{
    #pragma omp parallel for
    for (int i = 0; i < ncells_; i++)
    {
        dt_[i] = dt_value;
    }
}

double FlowSolver::compute_residual_norm() const
{
    const auto& area = mesh_.get_cells_area();

    double sum = 0.0;
    double total_volume = 0.0;
    #pragma omp parallel for
    for (int i = 0; i < ncells_; i++)
    {
        sum += residual_[i][0] * residual_[i][0] * area[i];
        total_volume += area[i];
    }

    if (total_volume == 0.0)
        return 0.0;

    return std::sqrt(sum / total_volume);
}

double FlowSolver::compute_residual_max_norm() const
{
    double max_res = 0.0;
    #pragma omp parallel for
    for (int i = 0; i < ncells_; i++)
    {
        max_res = std::max(max_res, std::abs(residual_[i][0]));
    }

    return max_res;
}

double FlowSolver::compute_residual_l2_raw() const
{
    double sum = 0.0;
    #pragma omp parallel for
    for (int i = 0; i < ncells_; i++)
    {
        sum += residual_[i][0] * residual_[i][0];
    }

    return std::sqrt(sum);
}


void FlowSolver::run_explicit(int niter, double cfl, double tol)
{
    for (int iter = 0; iter < niter; iter++)
    {
        compute_local_dt(cfl);
        compute_residual();

        
        double res_norm = compute_residual_norm();

        if (iter == 0)
        {
            initial_residual_ = res_norm;
        
            if (initial_residual_ < 1e-14)
                initial_residual_ = 1e-14;
        }
        double res_rel = res_norm / initial_residual_;

        
        double res_max = compute_residual_max_norm();
        double res_l2  = compute_residual_l2_raw();


        residual_history_.push_back(res_norm);

        if (iter % 100 == 0)
        {
            std::cout << "Iter " << iter
                      << "  Residual  Relative L2(rho) = " << res_rel
                      << "  Residual L2(rho) = " << res_norm
                      << "  L2raw(rho) = " << res_l2
                      << "  Linf(rho) = " << res_max
                      << std::endl;
        }

        if (res_norm < tol)
        {
            std::cout << "Convergence reached at iteration " << iter
                      << " with Residual L2(rho) = " << res_norm
                      << std::endl;
            break;
        }

        update_solution_explicit();
    }
}

int FlowSolver::get_ncells() const
{
    return ncells_;
}

const std::vector<std::array<double,4>>& FlowSolver::get_W() const
{
    return W_;
}

const std::vector<std::array<double,4>>& FlowSolver::get_residual() const
{
    return residual_;
}

const std::vector<double>& FlowSolver::get_dt() const
{
    return dt_;
}

void FlowSolver::write_solution_vtk(const std::string& filename) const
{
    std::ofstream out(filename);

    if (!out)
    {
        std::cerr << "Error: cannot open " << filename << std::endl;
        return;
    }

    const auto& X = mesh_.get_CoordX();
    const auto& Y = mesh_.get_CoordY();
    const auto& Z = mesh_.get_CoordZ();
    const auto& T = mesh_.get_Triangles();

    const int numPoints = static_cast<int>(X.size());
    const int numCells  = static_cast<int>(T.size());

    out << "# vtk DataFile Version 2.0\n";
    out << "CFD solution\n";
    out << "ASCII\n";
    out << "DATASET UNSTRUCTURED_GRID\n";

    // --------------------------------------------------
    // Points
    // --------------------------------------------------
    out << "POINTS " << numPoints << " float\n";
    for (int i = 0; i < numPoints; i++)
    {
        out << X[i] << " " << Y[i] << " " << Z[i] << "\n";
    }

    // --------------------------------------------------
    // Cells
    // For triangles: each cell = 3 node ids
    // VTK format line = "3 n0 n1 n2"
    // Total size = numCells * (1 + 3) = 4 * numCells
    // --------------------------------------------------
    out << "\nCELLS " << numCells << " " << 4 * numCells << "\n";
    for (int i = 0; i < numCells; i++)
    {
        out << "3 "
            << T[i][0] << " "
            << T[i][1] << " "
            << T[i][2] << "\n";
    }

    // --------------------------------------------------
    // Cell types
    // VTK_TRIANGLE = 5
    // --------------------------------------------------
    out << "\nCELL_TYPES " << numCells << "\n";
    for (int i = 0; i < numCells; i++)
    {
        out << "5\n";
    }

    // --------------------------------------------------
    // Cell data
    // --------------------------------------------------
    out << "\nCELL_DATA " << numCells << "\n";

    // Density
    out << "SCALARS density float 1\n";
    out << "LOOKUP_TABLE default\n";
    for (int i = 0; i < numCells; i++)
    {
        double rho = W_[i][0];
        out << rho << "\n";
    }

    // Pressure
    out << "\nSCALARS pressure float 1\n";
    out << "LOOKUP_TABLE default\n";
    for (int i = 0; i < numCells; i++)
    {
        double p = physics_.pressure(W_[i]);
        out << p << "\n";
    }

    // Velocity vector
    out << "\nVECTORS velocity float\n";
    for (int i = 0; i < numCells; i++)
    {
        double rho  = W_[i][0];
        double rhou = W_[i][1];
        double rhov = W_[i][2];

        double u = rhou / rho;
        double v = rhov / rho;

        out << u << " " << v << " " << 0.0 << "\n";
    }

    // Velocity magnitude
    out << "\nSCALARS velocity_magnitude float 1\n";
    out << "LOOKUP_TABLE default\n";
    for (int i = 0; i < numCells; i++)
    {
        double rho  = W_[i][0];
        double rhou = W_[i][1];
        double rhov = W_[i][2];

        double u = rhou / rho;
        double v = rhov / rho;

        double vel = std::sqrt(u * u + v * v);
        out << vel << "\n";
    }

    // Mach
    out << "\nSCALARS mach float 1\n";
    out << "LOOKUP_TABLE default\n";
    for (int i = 0; i < numCells; i++)
    {
        double rho  = W_[i][0];
        double rhou = W_[i][1];
        double rhov = W_[i][2];

        double u = rhou / rho;
        double v = rhov / rho;

        double vel = std::sqrt(u * u + v * v);
        double a   = physics_.sound_speed(W_[i]);

        double M = vel / a;
        out << M << "\n";
    }

    out.close();

    std::cout << "VTK solution written to: " << filename << std::endl;
}

void FlowSolver::write_residual_history(const std::string& filename) const
{
    std::ofstream out(filename);

    if (!out)
    {
        std::cerr << "Error: cannot open " << filename << std::endl;
        return;
    }

    out << "# iteration residual\n";
    for (size_t i = 0; i < residual_history_.size(); i++)
    {
        out << i << " " << residual_history_[i] << "\n";
    }

    out.close();
}

void FlowSolver::run_rk4(int niter, double cfl, double tol)
{
    const auto& cells_area = mesh_.get_cells_area();

    for (int iter = 0; iter < niter; iter++)
    {
        
        compute_local_dt(cfl);

        
        std::vector<std::array<double,4>> W0 = W_;

        
        double res_norm = compute_residual_norm();
        residual_history_.push_back(res_norm);

        if (iter % 100 == 0)
            std::cout << "Iter " << iter
                      << "  Residual L2(rho) = " << res_norm << std::endl;

        if (res_norm < tol) {
            std::cout << "Convergence at iter " << iter << std::endl;
            break;
        }

        // --- k1 ---
        compute_residual();   // R(W^n)
        std::vector<std::array<double,4>> k1(ncells_);
        for (int i = 0; i < ncells_; i++)
            for (int k = 0; k < 4; k++)
                k1[i][k] = -(dt_[i] / cells_area[i]) * residual_[i][k];

        // W = W^n + 0.5*k1
        for (int i = 0; i < ncells_; i++)
            for (int k = 0; k < 4; k++)
                W_[i][k] = W0[i][k] + 0.5 * k1[i][k];

        // --- k2 ---
        compute_residual();   // R(W^n + 0.5*k1)
        std::vector<std::array<double,4>> k2(ncells_);
        for (int i = 0; i < ncells_; i++)
            for (int k = 0; k < 4; k++)
                k2[i][k] = -(dt_[i] / cells_area[i]) * residual_[i][k];

        // W = W^n + 0.5*k2
        for (int i = 0; i < ncells_; i++)
            for (int k = 0; k < 4; k++)
                W_[i][k] = W0[i][k] + 0.5 * k2[i][k];

        // --- k3 ---
        compute_residual();   // R(W^n + 0.5*k2)
        std::vector<std::array<double,4>> k3(ncells_);
        for (int i = 0; i < ncells_; i++)
            for (int k = 0; k < 4; k++)
                k3[i][k] = -(dt_[i] / cells_area[i]) * residual_[i][k];

        // W = W^n + k3
        for (int i = 0; i < ncells_; i++)
            for (int k = 0; k < 4; k++)
                W_[i][k] = W0[i][k] + k3[i][k];

        // --- k4 ---
        compute_residual();   // R(W^n + k3)
        std::vector<std::array<double,4>> k4(ncells_);
        for (int i = 0; i < ncells_; i++)
            for (int k = 0; k < 4; k++)
                k4[i][k] = -(dt_[i] / cells_area[i]) * residual_[i][k];

        
        for (int i = 0; i < ncells_; i++)
            for (int k = 0; k < 4; k++)
                W_[i][k] = W0[i][k]
                         + (1.0/6.0) * (k1[i][k]
                         + 2.0 * k2[i][k]
                         + 2.0 * k3[i][k]
                         + k4[i][k]);

        
        
    }
}