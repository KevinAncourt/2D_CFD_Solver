# 2D CFD Solver

> A personal **prototype** 2D CFD solver in C++ for the **compressible Euler equations** on **unstructured triangular meshes**, using a **cell-centered finite-volume** formulation.

![Status](https://img.shields.io/badge/status-prototype-orange)
![Language](https://img.shields.io/badge/language-C%2B%2B17-blue)
![Build](https://img.shields.io/badge/build-CMake-brightgreen)
![Mesh](https://img.shields.io/badge/mesh-Gmsh%204.1-informational)

---

## Overview

This repository is a **personal development project** whose goal is to progressively build the core ingredients of a 2D CFD solver from scratch:

- Gmsh mesh reading
- compact local indexing
- face and connectivity reconstruction
- geometric preprocessing
- Euler physics utilities
- numerical flux computation
- boundary-condition handling
- local CFL time stepping
- explicit pseudo-time marching
- VTK output for ParaView post-processing

The current codebase is intentionally simple and readable. It should be viewed as a **prototype for learning, validation, and progressive extension**, not as a production-ready or industrial CFD code.

---

## Current Features

The solver already supports:

- reading **Gmsh 4.1** `.msh` files,
- converting Gmsh node tags to **compact local indices**,
- extracting **fluid triangles** and selected **boundary edges**,
- reconstructing **unique mesh faces** from cell connectivity,
- computing:
  - cell centers,
  - cell areas,
  - face centers,
  - face lengths,
  - face normals,
- orienting face normals consistently from **left cell** to **right cell**,
- solving the **2D compressible Euler equations**,
- using a **Rusanov** numerical flux,
- applying:
  - a **slip wall** boundary condition,
  - a simple **farfield** boundary condition,
- advancing the solution with an **explicit local-CFL scheme**,
- exporting mesh and solution fields to **VTK** for visualization in **ParaView**.

---

## Project Status

This project is currently at the stage of a **working prototype**.

### What it is

- a clean starting point for a custom CFD solver,
- a testbed for numerical methods,
- a personal codebase to better understand solver architecture and implementation details.

### What it is not

- a production solver,
- a fully verified and validated research code,
- an optimized HPC implementation,
- a complete Navier–Stokes framework.

---

## Numerical Model

The solver works with the conservative state vector

\[
W =
\begin{pmatrix}
\rho \\
\rho u \\
\rho v \\
\rho E
\end{pmatrix}
\]

where:

- \(\rho\) is the density,
- \(u, v\) are the velocity components,
- \(E\) is the total specific energy.

The pressure is recovered from the perfect-gas equation of state:

\[
p = (\gamma - 1)\left(\rho E - \frac{1}{2}\rho(u^2+v^2)\right)
\]

The physical flux projected onto a unit normal \(\mathbf{n} = (n_x,n_y)\) is:

\[
F(W)\cdot \mathbf{n} =
\begin{pmatrix}
\rho u_n \\
\rho u u_n + p n_x \\
\rho v u_n + p n_y \\
(\rho E + p)u_n
\end{pmatrix}
\]

with

\[
u_n = u n_x + v n_y
\]

The internal numerical flux currently used is the **Rusanov flux**:

\[
\hat{F}(W_L, W_R)
= \frac{1}{2}\left(F(W_L)+F(W_R)\right)\cdot \mathbf{n}
- \frac{1}{2}s_{\max}(W_R-W_L)
\]

where

\[
s_{\max} = \max\left(|u_{n,L}| + a_L,\; |u_{n,R}| + a_R\right)
\]

and \(a = \sqrt{\gamma p/\rho}\) is the sound speed.

The explicit update is written as

\[
W_i^{n+1} = W_i^n - \frac{\Delta t_i}{|\Omega_i|} R_i
\]

with \(|\Omega_i|\) the cell area and \(R_i\) the sum of fluxes over the cell faces.

The local pseudo-time step is computed from a CFL estimate:

\[
\Delta t_i = \mathrm{CFL}\,\frac{|\Omega_i|}{\sum_f \lambda_f}
\]

where each \(\lambda_f\) is based on a local estimate of \(|u_n| + a\).

---

## Boundary Conditions

### Slip Wall

The wall boundary condition currently implemented is an **inviscid slip wall**.

At the wall, there is:

- **no mass flux** through the boundary,
- **no energy flux** through the boundary,
- only the **pressure force normal to the wall**.

This leads to the wall flux

\[
\Phi_{\text{wall}} =
\begin{pmatrix}
0 \\
p n_x S_f \\
p n_y S_f \\
0
\end{pmatrix}
\]

where \(S_f\) is the face length.

### Farfield

The farfield boundary condition is currently implemented in a simple way by computing a numerical flux between:

- the internal cell state,
- a prescribed freestream state.

This is appropriate for a prototype, although it is not yet a full characteristic farfield treatment.

---

## Repository Structure

```text
2D_CFD_Solver/
├── CMakeLists.txt
├── README.md
├── include/
│   ├── bc/
│   │   └── boundary_conditions.hpp
│   ├── flux/
│   │   └── rusanov_flux.hpp
│   ├── mesh/
│   │   ├── mesh_compute.hpp
│   │   └── mesh_reader.hpp
│   ├── physics/
│   │   └── euler_physics.hpp
│   └── solver/
│       └── solver.hpp
├── src/
│   ├── bc/
│   │   └── boundary_conditions.cpp
│   ├── flux/
│   │   └── rusanov_flux.cpp
│   ├── mesh/
│   │   ├── mesh_compute.cpp
│   │   └── mesh_reader.cpp
│   ├── physics/
│   │   └── euler_physics.cpp
│   ├── solver/
│   │   └── solver.cpp
│   └── main.cpp
├── mesh_generation/
│   ├── msh/
│   │   └── naca0012_farfield150.msh
│   └── script/
│       └── naca0012.cpp
├── output/
│   ├── mesh.vtk
│   ├── wall.vtk
│   ├── farfield.vtk
│   ├── normals.vtk
│   └── normals_boundary.vtk
└── build/
```

---

## Main Components

### `MeshReader`
Reads the Gmsh mesh and stores data in contiguous arrays:

- nodal coordinates,
- fluid triangles,
- selected boundary edges,
- Gmsh-to-local index mappings.

### `MeshCompute`
Builds all geometric quantities required by the solver:

- cell centers,
- cell areas,
- unique faces,
- left/right cell connectivity,
- oriented normals,
- boundary face types.

### `EulerPhysics`
Provides the main physical building blocks:

- conservative/primitive conversions,
- pressure,
- sound speed,
- physical flux.

### `RusanovFlux`
Computes the numerical flux at internal interfaces and for the simple farfield treatment.

### `BoundaryConditions`
Applies the current boundary-condition models:

- inviscid slip wall,
- simple farfield state.

### `FlowSolver`
Handles the solution loop:

- initialization,
- residual assembly,
- local time-step computation,
- explicit update,
- residual monitoring,
- VTK export.

---

## Build

From the project root:

```bash
mkdir -p build
cd build
cmake ..
make -j
```

The executable is generated in the `build/` directory.

---

## Run

From `build/`:

```bash
./solver
```

The current case configured in `main.cpp` corresponds to a **NACA0012** test case with freestream conditions such as:

- \(\gamma = 1.4\)
- \(M_\infty = 1.5\)
- \(\alpha = 1^\circ\)
- \(\rho_\infty = 1\)
- \(p_\infty = 1/(\gamma M_\infty^2)\)

The computed solution is exported to a `solution.vtk` file.

---

## Post-processing

The generated VTK files can be visualized in **ParaView**.

Typical exported quantities include:

- mesh geometry,
- wall and farfield boundaries,
- density,
- pressure,
- velocity vector,
- velocity magnitude,
- Mach number.

---

## Mesh Generation

The `mesh_generation/script/` directory contains a C++ script based on the **Gmsh API** to generate a test mesh around a **NACA0012** airfoil with a circular farfield of radius `150`.

At the moment, the mesh reading and boundary identification are still partly tied to this specific case. Future improvements should make:

- boundary-group handling more generic,
- boundary recognition independent of hard-coded entity tags,
- multi-case usage cleaner and more modular.

---

## Current Limitations

This repository should clearly be viewed as a **prototype**. Current limitations include:

- only a **Rusanov** flux is implemented,
- the spatial discretization is currently **first-order**,
- no **MUSCL** reconstruction,
- no limiter,
- no implicit solver,
- no viscous / Navier–Stokes terms,
- boundary handling is still partly **case-specific**,
- robustness near strong shocks is still limited,
- no automated testing framework,
- no claim of industrial-grade robustness or HPC optimization.

---

## Possible Roadmap

Natural next steps for the project include:

1. making boundary-condition handling more generic,
2. adding **MUSCL** reconstruction,
3. introducing limiters,
4. improving shock robustness,
5. adding Roe / HLLC for comparison,
6. implementing an implicit solver,
7. extending toward Navier–Stokes,
8. preparing a cleaner architecture for performance and parallelism.

---

## Disclaimer

This project is a **personal prototype** intended for solver development, learning, and progressive validation. It is useful as a base for understanding and extending the core pieces of a CFD solver, but it still requires substantial work in robustness, verification, and validation before being considered a reliable research or production tool.