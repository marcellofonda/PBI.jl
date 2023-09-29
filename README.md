# X-ray Propagation-Based Phase-Contrast Imaging and Retrieval

This repository contains some code about Propagation-Based Phase-Contrast imaging of homogeneous samples using X-rays. More about it at [marcellofonda.github.io/Tesi](https://marcellofonda.github.io/Tesi)

## How to use it
1. Clone the repository using git, or copy-paste source code into your own files.

2. Download and install the latest version of Julia;

3. In Julia, type `cd("path\\to\\your\\directory")`. Beware to write double `\` signs, since they are used as special characters in Julia strings

4. If these packages are not already installed, type
```Julia
] add LinearAlgebra
] add Images
] add FFTW
] add FileIO
] add ImageMagick
] add ImageTransformations
] add Interpolations
] add https://github.com/JuliaPhysics/PhysicalOptics.jl
```

5. Copy the three `.jl` files from the `src` folder to your directory

6. You can now use the functions of this project your Julia script by adding the lines
```Julia
include("simulation.jl")
include("retrieval.jl")
```

## Main features

In the `simulation.jl` file, there are four main functions of interest:

```Julia
AbsoprtionRadiography(thickness, β, k, I_0)
```
which simulates a traditional absorption technique through the Beer-Lambert law:
\[
    I_T = I_0 e^{-2k\beta \mathrm{thickness}}
\]
* `thickness` is expected to be a matrix of thickness values representing the sample; 
* `β` is the imaginary part of the refractive index;
* `k` is the wave number;
* `I_0` the initial intensity of radiation
Note: `thickness` and `k` must be expressed in the same unit system.

```Julia
LaplacianRadiography(thickness, δ, β, k, I_0, R, pixelsize)
```
which simulates a Propagation-Based Phase-Contrast image through the analytical expression derived from the Transport of Intensity Equation. Namely:
\[
    I_R = \left(1-\frac{R\delta}{\mu}\nabla^2\right)I_0 e^{-2k\beta \mathrm{thickness}
\]
* `thickness` is expected to be a matrix of thickness values representing the sample;
* `δ` is the real decrement of the refractive index;
* `β` is the imaginary part of the refractive index;
* `k` is the wave number;
* `I_0` the initial intensity of radiation
* `R` is the propagation distance
* `pixelsize` is the size of the pixel in the imaging system
Note: `thickness`, `k`, `R`, and `pixelsize` must be expressed in the same unit system.



