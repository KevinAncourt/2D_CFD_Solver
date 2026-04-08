#include "flux/rusanov_flux.hpp"

#include <algorithm>
#include <cmath>

RusanovFlux::RusanovFlux(const EulerPhysics& physics)
    : physics_(physics)
{
}

std::array<double,4> RusanovFlux::compute(const std::array<double,4>& WL,
                                          const std::array<double,4>& WR,
                                          double nx,
                                          double ny,
                                          double Sf) const
{
    const std::array<double,4> VL = physics_.cons_to_prim(WL);
    const std::array<double,4> VR = physics_.cons_to_prim(WR);

    double uL = VL[1];
    double vL = VL[2];
    double pL = VL[3];
    double rhoL = VL[0];

    double uR = VR[1];
    double vR = VR[2];
    double pR = VR[3];
    double rhoR = VR[0];

    double aL = std::sqrt(physics_.get_gamma() * pL / rhoL);
    double aR = std::sqrt(physics_.get_gamma() * pR / rhoR);

    double unL = uL * nx + vL * ny;
    double unR = uR * nx + vR * ny;

    double smax = std::max(std::abs(unL) + aL,
                           std::abs(unR) + aR);

    std::array<double,4> FL = physics_.physical_flux(WL, nx, ny, Sf);
    std::array<double,4> FR = physics_.physical_flux(WR, nx, ny, Sf);

    std::array<double,4> flux;

    for (int k = 0; k < 4; k++)
    {
        flux[k] = 0.5 * (FL[k] + FR[k])
                - 0.5 * smax * Sf * (WR[k] - WL[k]); 
    }

    return flux;
}

