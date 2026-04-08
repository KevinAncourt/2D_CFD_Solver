# 2D Euler CFD Solver (Prototype)

## 📌 Overview

This project is a **2D finite volume solver for compressible Euler
equations** on unstructured triangular meshes.

It is developed as a **personal prototype** to better understand and
implement: - Numerical methods for CFD - Unstructured mesh handling -
Flux-based finite volume schemes

> ⚠️ This is a **prototype / research-oriented code**, not a
> production-ready solver.

------------------------------------------------------------------------

## ✨ Features

-   Finite Volume Method (cell-centered)
-   Unstructured triangular mesh (Gmsh)
-   Euler equations (compressible flow)
-   Rusanov (local Lax-Friedrichs) flux
-   Local time stepping (CFL-based)
-   Explicit pseudo-time integration
-   Boundary conditions:
    -   Slip wall (Euler wall)
    -   Farfield (Rusanov-based, freestream state)
-   Residual monitoring (L2, relative, Linf)
-   VTK output (ParaView visualization)

------------------------------------------------------------------------

## 🧮 Governing Equations

The solver uses the **2D compressible Euler equations** in conservative
form:

    W = [ rho, rho*u, rho*v, rho*E ]

### Pressure relation

![pressure](https://latex.codecogs.com/png.latex?p%20%3D%20(%5Cgamma%20-%201)(%5Crho%20E%20-%200.5%5Crho(u%5E2%2Bv%5E2)))

### Flux projected along a normal

![flux](https://latex.codecogs.com/png.latex?F(W)%20%5Ccdot%20n%20%3D%20%5Cbegin%7Bpmatrix%7D%20%5Crho%20u_n%20%5C%5C%20%5Crho%20u%20u_n%20%2B%20p%20n_x%20%5C%5C%20%5Crho%20v%20u_n%20%2B%20p%20n_y%20%5C%5C%20(%5Crho%20E%20%2B%20p)%20u_n%20%5Cend%7Bpmatrix%7D)

------------------------------------------------------------------------

## 🔢 Numerical Method

### Finite Volume formulation

    W_i^{n+1} = W_i^n - (dt_i / |Omega_i|) * sum(F_faces)

### Numerical flux (Rusanov)

![rusanov](https://latex.codecogs.com/png.latex?F%5E%7B*%7D%20%3D%200.5(F_L%20%2B%20F_R)%20-%200.5%20s_%7Bmax%7D%20(W_R%20-%20W_L))

-   Robust and stable
-   Introduces numerical diffusion

------------------------------------------------------------------------

## 🌊 Boundary Conditions

### Wall (slip condition)

-   No mass flux
-   No energy flux
-   Only pressure contribution

![wall](https://latex.codecogs.com/png.latex?F_%7Bwall%7D%20%3D%20%5B0%2C%20p%20n_x%2C%20p%20n_y%2C%200%5D)

------------------------------------------------------------------------

### Farfield (current implementation)

The farfield boundary is implemented using a **Rusanov flux with
freestream state**:

    F_bc = F_rusanov(W_inside, W_inf)

-   Simple and robust
-   Does not use characteristic decomposition
-   Slightly diffusive

> ⚠️ A characteristic-based boundary condition would improve convergence
> behavior.

------------------------------------------------------------------------

## 🧱 Mesh Handling

-   Mesh generated with **Gmsh**
-   Read from `.msh` (v4.1 format)
-   Internal storage optimized for performance:
    -   `CoordX`, `CoordY`, `CoordZ`
    -   Triangle connectivity
    -   Face reconstruction
-   Face data:
    -   Normals
    -   Surface
    -   Left/right cells

------------------------------------------------------------------------

## ▶️ How to Run

### 1. Generate mesh

``` bash
cd mesh_generation/script
g++ naca0012.cpp -o naca0012 -lgmsh
./naca0012
```

### 2. Build solver

``` bash
mkdir build && cd build
cmake ..
make
```

### 3. Run simulation

``` bash
./solver
```

------------------------------------------------------------------------

## 📊 Output

-   `.vtk` files (ParaView)
-   Residual history file

Typical fields: - Density - Pressure - Velocity - Mach number

------------------------------------------------------------------------

## 📈 Example Results (Minf = 1.5, CFL = 0.1)

### Mesh

<p align="center">
  <img src="M0_5_CFL0_001/mesh.png" width="300"/>
</p>


### Solution (Mach number)

<p align="center">
  <img src="M1_5_AOA_1_CFL0_1/Naca12.png" width="300"/>
</p>

### Residual history

<p align="center">
  <img src="M1_5_AOA_1_CFL0_1/residus.png" width="300"/>
</p>

------------------------------------------------------------------------

## ⚠️ Numerical Behavior

For the supersonic test case ($M_\infty = 1.5$, AoA = $1^\circ$), the solver produces a qualitatively reasonable compressible flow field around the airfoil.

From the present results, the following observations can be made:

- A strong perturbation is visible near the leading edge, consistent with the expected supersonic compression effects.
- The downstream flow structure is captured qualitatively, but compressible features remain broadly smeared.
- No sharp shock is clearly resolved in the Mach contours, which indicates significant numerical diffusion.
- This behavior is expected from the current **first-order finite-volume formulation** combined with a **Rusanov flux**, which is robust but strongly diffusive.
- The current mesh is sufficient to recover the global flow pattern, but it is likely too coarse to capture thin shock structures accurately with the present scheme.
- The residual decreases overall, which indicates a stable evolution toward a steady state, but convergence remains relatively slow.

These results should therefore be interpreted as **prototype-level numerical experiments**:
the solver captures the main qualitative behavior of the flow, but the solution is not yet suitable for high-fidelity shock resolution or quantitative aerodynamic prediction.

Typical sources of limitation at this stage are:
- first-order spatial accuracy,
- diffusive Rusanov flux,
- moderate mesh resolution near compressible structures,
- simplified boundary-condition treatment.


------------------------------------------------------------------------

## 📈 Example Results (Minf = 0.5, CFL = 0.001)



### Solution (Mach number)

<p align="center">
  <img src="M0_5_CFL0_001/Naca12.png" width="300"/>
</p>

### Residual history

<p align="center">
  <img src="M0_5_CFL0_001/residus.png" width="300"/>
</p>

------------------------------------------------------------------------

## ⚠️ Numerical Behavior

-   Residuals decrease initially
-   A minimum is reached
-   Then a slight increase appears

This is typical for: - First-order scheme - Rusanov flux (diffusive) -
Simplified farfield BC



------------------------------------------------------------------------

## ⚠️ Limitations

-   First-order scheme (no MUSCL)
-   Diffusive flux (Rusanov)
-   Simplified farfield boundary
-   No viscous terms (Euler only)
-   Explicit solver only (no implicit yet)
-   No parallelization

------------------------------------------------------------------------

## 🚀 Future Improvements

-   MUSCL reconstruction (2nd order)
-   Roe / HLLC flux
-   Characteristic farfield BC
-   Implicit solver
-   Mesh adaptation (goal-oriented :eyes:)
-   MPI parallelization

------------------------------------------------------------------------

## 👨‍🔬 Author

Kevin Ancourt\
PhD in Computational Fluid Dynamics 

Specialized in: - Adjoint methods - Goal-oriented mesh adaptation -
High-performance scientific computing

📄 **PhD Thesis (French)**:  
[Analyse et mise en œuvre d'une méthode de raffinement de maillage ciblé pour les écoulements compressibles](https://theses.hal.science/tel-05480104)

📄 **Google Scholar**:  
[Link to GS](https://scholar.google.com/citations?user=QssUs3YAAAAJ&hl=fr)

------------------------------------------------------------------------

## 📜 License

This project is released for **educational and research purposes**.
