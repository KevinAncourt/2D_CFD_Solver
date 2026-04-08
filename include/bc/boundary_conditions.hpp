#pragma once

#include <array>
#include "physics/euler_physics.hpp"
#include "flux/rusanov_flux.hpp"

class BoundaryConditions
{
public:
    explicit BoundaryConditions(const EulerPhysics& physics);

    void set_farfield_state(const std::array<double,4>& W_farfield);

    std::array<double,4> wall_flux(const std::array<double,4>& W_inside,
                                   double nx,
                                   double ny,
                                   double Sf) const;

    std::array<double,4> farfield_flux(const std::array<double,4>& W_inside,
                                       double nx,
                                       double ny,
                                       double Sf) const;

private:
    const EulerPhysics& physics_;
    RusanovFlux rusanov_flux_;
    std::array<double,4> W_farfield_;
};
