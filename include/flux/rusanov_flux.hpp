#pragma once

#include <array>
#include "physics/euler_physics.hpp"

class RusanovFlux
{
public:
    explicit RusanovFlux(const EulerPhysics& physics);

    std::array<double,4> compute(const std::array<double,4>& WL,
                                 const std::array<double,4>& WR,
                                 double nx,
                                 double ny,
                                 double Sf) const;


private:
    const EulerPhysics& physics_;
};
