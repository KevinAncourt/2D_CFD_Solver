#include "bc/boundary_conditions.hpp"
#include <cmath>


BoundaryConditions::BoundaryConditions(const EulerPhysics& physics)
    : physics_(physics), rusanov_flux_(physics), W_farfield_{1.0, 0.0, 0.0, 1.0}
{
}

void BoundaryConditions::set_farfield_state(const std::array<double,4>& W_farfield)
{
    W_farfield_ = W_farfield;
}

std::array<double,4> BoundaryConditions::wall_flux(const std::array<double,4>& W_inside,
                                                   double nx,
                                                   double ny,
                                                   double Sf) const
{
    double p = physics_.pressure(W_inside);

    std::array<double,4> flux;
    flux[0] = 0.0;
    flux[1] = p * nx * Sf;
    flux[2] = p * ny * Sf;
    flux[3] = 0.0;

    return flux;
}

std::array<double,4> BoundaryConditions::farfield_flux(const std::array<double,4>& W_inside,
                                                       double nx,
                                                       double ny,
                                                       double Sf) const
{
    static int count = 0;
    if (count < 10)
    {
        std::cout << "farfield_flux called" << std::endl;
        count++;
    }
    return rusanov_flux_.compute(W_inside, W_farfield_, nx, ny, Sf);
    
}

// Work in Progres : Better Farfield Conditions 
// std::array<double,4> BoundaryConditions::farfield_flux(const std::array<double,4>& W_inside,
//                                                        double nx,
//                                                        double ny,
//                                                        double Sf) const
// {
        
//     const std::array<double,4> V_in  = physics_.cons_to_prim(W_inside);
//     const std::array<double,4> V_inf = physics_.cons_to_prim(W_farfield_);

//     double rho_in = V_in[0],  u_in = V_in[1],  v_in = V_in[2],  p_in = V_in[3];
//     double rho_f  = V_inf[0], u_f  = V_inf[1], v_f  = V_inf[2], p_f  = V_inf[3];

//     double a_in = physics_.sound_speed(W_inside);
//     double a_f  = std::sqrt(physics_.get_gamma() * p_f / rho_f);


//     double Rm = (u_in * nx + v_in * ny) - 2.0 * a_in / (physics_.get_gamma() - 1.0); // sortant
//     double Rp = (u_f  * nx + v_f  * ny) + 2.0 * a_f  / (physics_.get_gamma() - 1.0); // entrant

//     double un_bc = 0.5 * (Rp + Rm);
//     double a_bc  = 0.25 * (physics_.get_gamma() - 1.0) * (Rp - Rm);

//     double un_in = u_in * nx + v_in * ny;
//     double un_f  = u_f  * nx + v_f  * ny;


//     double tx = -ny, ty = nx;
//     double ut_bc;
//     if (un_in >= 0.0)  // face sortante
//         ut_bc = u_in * tx + v_in * ty;
//     else               // face entrante
//         ut_bc = u_f  * tx + v_f  * ty;

//     double u_bc = un_bc * nx + ut_bc * tx;
//     double v_bc = un_bc * ny + ut_bc * ty;


//     double s_bc;
//     if (un_in >= 0.0)
//         s_bc = p_in / std::pow(rho_in, physics_.get_gamma());
//     else
//         s_bc = p_f  / std::pow(rho_f,  physics_.get_gamma());

//     double rho_bc = std::pow(a_bc * a_bc / (physics_.get_gamma() * s_bc),
//                              1.0 / (physics_.get_gamma() - 1.0));
//     double p_bc   = s_bc * std::pow(rho_bc, physics_.get_gamma());

//     std::array<double,4> V_bc = {rho_bc, u_bc, v_bc, p_bc};
//     std::array<double,4> W_bc = physics_.prim_to_cons(V_bc);

//     return physics_.physical_flux(W_bc, nx, ny, Sf);
    
// }