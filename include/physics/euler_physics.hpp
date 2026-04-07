#pragma once

#include "mesh/mesh_compute.hpp"

class EulerPhysics {

public:
    EulerPhysics(double gamma = 1.4);
    std::array<double,4> cons_to_prim(const std::array<double,4>& W) const;
    std::array<double,4> prim_to_cons(const std::array<double,4>& V) const;
    double pressure(const std::array<double,4>& W) const;
    double sound_speed(const std::array<double,4>& W) const;
    std::array<double,4> physical_flux(const std::array<double,4>& W, double nx, double ny, double Sf) const;
    double get_gamma() const;
    



private:
    double gamma_;

    


};