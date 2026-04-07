#include "solver/solver.hpp"

FlowSolver::FlowSolver(const MeshCompute& mesh,
                       const EulerPhysics& physics)
    : mesh_(mesh), physics_(physics)
{
    ncells_ = static_cast<int>(mesh_.get_centercells_X().size());

    W_.resize(ncells_);
    residual_.resize(ncells_);
    dt_.resize(ncells_, 0.0);
}


