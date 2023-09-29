# X-ray Propagation-Based Phase-Contrast Imaging and Retrieval

This repository contains some code about phase-contrast imaging using X-rays. More about it at [marcellofonda.github.io/Tesi](https://marcellofonda.github.io/Tesi)

## How to use it
1. Clone the repository using git, or copy-paste source code into your own files.

2. Download and install the latest version of Julia;

3. In Julia, type `cd("path\\to\\your\\directory")`. Beware to double `\` signs, since they are used as special characters in Julia strings

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

5. Copy the three `.jl` files of this project in your directory

6. You can now use the functions of this project your Julia script by adding the lines
```Julia
include("simulation.jl")
include("retrieval.jl")
```
