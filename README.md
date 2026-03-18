# SMA FEM

A nonlinear finite element (FEM) solver for thermomechanical behavior of shape memory alloys (SMA), including phase transformations.

## Overview

This project implements a 3D finite element solver for shape memory alloys based on a constitutive model with phase transformations developed by A.A. Movchan and co-authors.

The implementation includes:

* phase transformation between austenite and martensite
* internal variables evolution (martensite fraction, transformation strain)
* temperature-dependent material response

---

## Model Description

The implemented constitutive model is based on a micromechanical and thermodynamically consistent framework for shape memory alloys, developed in Russian-language literature by Movchan et al.

The model describes nonlinear deformation of shape memory alloys undergoing phase and structural transformations between austenite and martensite.

Key features of the model:

- **Micromechanical foundation**  
  The material is treated as a representative volume with distributed internal microstresses.

- **Phase transformation modeling**  
  Forward (A → M) and reverse (M → A) transformations are governed by temperature and stress state.

- **Internal variables**  
  - martensite volume fraction `q`  
  - transformation strain `ω`  

- **Nonlinear constitutive behavior**  
  The model captures:
  - variation of elastic moduli during transformation  
  - accumulation of transformation strain  
  - hysteresis effects  
  - pseudoelastic and shape memory behavior  


- **Stress update scheme**  
  The constitutive equations are integrated in time using a tangent operator:

  σₙ₊₁ = σₙ + Cₙ : (εₙ₊₁ − εₙ − ΔT · Aₙ)

The model combines elements of micromechanics, phenomenology, and thermodynamics, allowing it to reproduce key features of shape memory alloy behavior, including phase-dependent stiffness, transformation-induced strain.

---

## Notes

- This is a research-oriented implementation and not a production-ready FEM code  
- The current implementation is known to contain bugs and is not fully validated  
- Numerical stability issues may be present, particularly in the time integration of the constitutive model  
- The solver should be treated as an experimental prototype rather than a reliable simulation tool  
- Some parts of the implementation are simplified for clarity and may require further refinement  

---

## Dependencies

* C++17 compatible compiler
* CMake ≥ 3.16
* Eigen3

---

## Build

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
make -j4
```

---

## Usage

1. Prepare a mesh (see below)
2. Adjust parameters in `consts.h`
3. Run:

```bash
./build/shape_memory_fem
```

Results are exported in VTK format and can be visualized in ParaView.

---

## Mesh Generation

Meshes should be generated using Gmsh.

The solver expects a pre-generated mesh (e.g. `.vtk` format).
You can convert meshes from Gmsh if needed.


## References

1. Мовчан А. А., Мовчан И. А., Сильченко Л. Г.
   Микромеханическая модель нелинейного деформирования сплавов с памятью формы при фазовых и структурных превращениях 
   // Механика твердого тела. 2010. № 3. 
   (Micromechanical model of nonlinear deformation of shape memory alloys under phase and structural transformations)

2. Мовчан А. А., Сильченко Л. Г., Казарина С. А., Тант Зин Аунг  
   Определяющие соотношения для сплавов с памятью формы — микромеханика, феноменология, термодинамика // 
   Ученые записки Казанского университета. Серия: Физико-математические науки. 2010. Т. 152, кн. 4. 
   (Constitutive relations for shape memory alloys: micromechanics, phenomenology, and thermodynamics)
   
3. Master thesis, included in the repository (in Russian).

