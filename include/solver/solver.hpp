#pragma once

#include <vector>
#include "mesh/mesh_compute.hpp"
#include "physics/euler_physics.hpp"

class FlowSolver
{
public:
    FlowSolver(const MeshCompute& mesh,
               const EulerPhysics& physics);

private:
    const MeshCompute& mesh_;
    const EulerPhysics& physics_;
    int ncells_;
    std::vector<std::array<double,4>> W_;        // conservative variables
    std::vector<std::array<double,4>> residual_; 
    std::vector<double> dt_;                    
};