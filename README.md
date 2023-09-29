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


### `simulation.jl`
In the `simulation.jl` file, there are four main functions of interest:

```Julia
AbsoprtionRadiography(thickness, β, k, I_0)
```
which simulates a traditional absorption radiography image through the Beer-Lambert law:
$$I_T = I_0 e^{-2k\beta\cdot \mathrm{thickness}}$$
<details>
<summary>Arguments description</summary>

* `thickness` is expected to be a matrix of thickness values representing the sample; 
* `β` is the imaginary part of the refractive index;
* `k` is the wave number;
* `I_0` the initial intensity of radiation.
Note: `thickness` and `k` must be expressed in the same unit system.
</details>

```Julia
LaplacianRadiography(thickness, δ, β, k, I_0, R, pixelsize)
```
which simulates a Propagation-Based Phase-Contrast image through the analytical expression derived from the Transport of Intensity Equation. Namely:
$$I_R = \left(1-\frac{R\delta}{2k\beta}\nabla^2\right)I_0 e^{-2k\beta\cdot \mathrm{thickness}}$$
<details>
<summary>Arguments description</summary>

* `thickness` is expected to be a matrix of thickness values representing the sample;
* `δ` is the real decrement of the refractive index;
* `β` is the imaginary part of the refractive index;
* `k` is the wave number;
* `I_0` the initial intensity of radiation;
* `R` is the propagation distance;
* `pixelsize` is the size of the pixel in the imaging system.
Note: `thickness`, `k`, `R`, and `pixelsize` must be expressed in the same unit system.
</details>

```Julia
FresnelRadiography(thickness, δ, β, k, I_0, R, pixelsize)
```
which simulates a Propagation-Based Phase-Contrast image by propagating the absorption image through the Fresnel Propagation Integral:
$$\psi_R(x,y)=\frac{e^{ikR}}{i\lambda R}\int \exp \left( \frac{i\pi}{\lambda R}[(x-x')^2+(y-y')^2] \right) \psi_0(x',y')dx'dy'$$
<details>
<summary>Arguments description</summary>

* `thickness` is expected to be a matrix of thickness values representing the sample;
* `δ` is the real decrement of the refractive index;
* `β` is the imaginary part of the refractive index;
* `k` is the wave number;
* `I_0` the initial intensity of radiation;
* `R` is the propagation distance;
* `pixelsize` is the size of the pixel in the imaging system.
Note: `thickness`, `k`, `R`, and `pixelsize` must be expressed in the same unit system.
</details>



```Julia
montecarlo_radiography(size_x, size_y, photons_per_pixel, t, grad_t, δ, β, k, R, pixelsize; multithreading = true)
```
which simulates a Propagation-Based Phase-Contrast image through a MC simulation.

<details>
<summary>Arguments description</summary>

* `size_x` and `size_y` are the dimensions of the desired output image;
* `photons_per_pixel` is the average number of photons per pixel to be simulated;
* `t` is a function representing the thickness distribution of the sample, e.g.:
`thickness::Float64 = t(x::Tuple{Float64,Float64})`;
* `grad_t` is a function representing the gradient of `t`, e.g.: `grad_thickness::Tuple{Float64,Float64} = grad_t(x::Tuple{Float64,Float64})`;
* `δ` is the real decrement of the refractive index;
* `β` is the imaginary part of the refractive index;
* `k` is the wave number;
* `I_0` the initial intensity of radiation;
* `R` is the propagation distance;
* `pixelsize` is the size of the pixel in the imaging system.
* `multithreading` is an optional argument which defaults to true and controls whether the code should be executed on all available threads. You can check the number of available threads by:
```Julia
julia> Threads.nthreads()
```
Note: `thickness`, `k`, `R`, and `pixelsize` must be expressed in the same unit system.
</details>

### `retrieval.jl`

In the `retrieval.jl` file, there are two main functions of interest:


```Julia
PhaseRetrieve(input_image, δ, β, k, R, pixelsize)
```
which executes the numerical phase retrieval on a Propagation-Based Phase-Contrast image through the Paganin Method:
$$t = -\frac1{2k\beta}\ln\left(\mathcal{F}^{-1}\left(\frac{\mathcal{F}[I_R/I_0] (u,v)}{1+\frac{R\delta}{2k\beta}(u^2+v^2)}\right)\right),$$
assuming $I_0=1$.
<details>
<summary>Arguments description</summary>

* `input_image` is the measured image;
* `δ` is the real decrement of the refractive index;
* `β` is the imaginary part of the refractive index;
* `k` is the wave number;
* `R` is the propagation distance;
* `pixelsize` is the size of the pixel in the imaging system.

Note: `k`, `R`, and `pixelsize` must be expressed in the same unit system. The output will be in that unit system.
</details>


```Julia
DarkFieldRetrieve(input_image, δ, β, R, λ, pixelsize, ε; scale = 1, save_images = true)
```
which executes the numerical retrieval of dark-field signal on a Propagation-Based Phase-Contrast image through the Gureyev Method:
$$A_{\mathrm{Born}} = \frac1{I_{\mathrm{TIE}}}\mathcal{F}^{-1}\left(\frac{\hat I_R(u,v)-\hat I_{R,\mathrm{TIE}}(u,v)}{-2\sqrt{1+\gamma^2}\cos{[\pi\lambda R(u^2+v^2)-\arctan(\gamma)]}}\right),$$
with $\gamma=\delta/\beta$ and assuming $I_0=1$.
<details>
<summary>Arguments description</summary>

* `input_image` is the measured image;
* `δ` is the real decrement of the refractive index;
* `β` is the imaginary part of the refractive index;
* `R` is the propagation distance;
* `λ` is the wavelength;
* `pixelsize` is the size of the pixel in the imaging system;
* `ε` is the dimensionless parameter for Tykhonoff regularization;
* `scale` is the optional upscaling factor (default value 1) to use in forward propagations to avoid artifacts due to sampling;
* `save_images` is the optional parameter to save all intermediate images and defaults to `false`.

Note: `λ`, `R`, and `pixelsize` must be expressed in the same unit system.

The output is a `Tuple` containing `A_Born`, `I_TIE_PR` (TIE Phase-Retrieval), and, optionally (if `save_images` is set to `true`), `images_paths`, an array containing the paths to the saved images.
</details>

