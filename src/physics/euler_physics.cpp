#include "physics/euler_physics.hpp"
#include <cmath>

EulerPhysics::EulerPhysics(double gamma)
    : gamma_(gamma)
{
}

std::array<double,4> EulerPhysics::cons_to_prim(const std::array<double,4>& W) const
{
    double rho  = W[0];
    double rhou = W[1];
    double rhov = W[2];
    double rhoE = W[3];

    double u = rhou / rho;
    double v = rhov / rho;
    double p = pressure(W);

    return {rho, u, v, p};
}

std::array<double,4> EulerPhysics::prim_to_cons(const std::array<double,4>& V) const
{
    double rho = V[0];
    double u   = V[1];
    double v   = V[2];
    double p   = V[3];

    double rhou = rho * u;
    double rhov = rho * v;
    double rhoE = p / (gamma_ - 1.0) + 0.5 * rho * (u * u + v * v);

    return {rho, rhou, rhov, rhoE};
}

double EulerPhysics::pressure(const std::array<double,4>& W) const
{
    double rho  = W[0];
    double rhou = W[1];
    double rhov = W[2];
    double rhoE = W[3];

    double u = rhou / rho;
    double v = rhov / rho;

    double kinetic = 0.5 * rho * (u * u + v * v);
    double p = (gamma_ - 1.0) * (rhoE - kinetic);

    return p;
}

double EulerPhysics::sound_speed(const std::array<double,4>& W) const
{
    double rho = W[0];
    double p   = pressure(W);

    double a = std::sqrt(gamma_ * p / rho);

    return a;
}

std::array<double,4> EulerPhysics::physical_flux(const std::array<double,4>& W,
                                                 double nx,
                                                 double ny,
                                                 double Sf) const
{
    double rho  = W[0];
    double rhou = W[1];
    double rhov = W[2];
    double rhoE = W[3];

    double u = rhou / rho;
    double v = rhov / rho;
    double p = pressure(W);

    double un = u * nx + v * ny;

    std::array<double,4> flux;
    // Normals are unit vectors → scale by face area (Sf) to obtain the flux through the face
    flux[0] = Sf * (rho * un);
    flux[1] = Sf * (rho * u * un + p * nx);
    flux[2] = Sf * (rho * v * un + p * ny);
    flux[3] = Sf * ((rhoE + p) * un);

    return flux;
}

double EulerPhysics::get_gamma() const
{
    return gamma_;
}