# IonSim.jl
============================

<p align="center">
  <img src="https://github.com/HaeffnerLab/IonSim.jl/blob/media/logo3_SM.svg?raw=true", width="450px">
</p>

[![test status](https://github.com/HaeffnerLab/IonSim.jl/actions/workflows/test.yml/badge.svg)](https://github.com/HaeffnerLab/IonSim.jl/actions/workflows/test.yml)
[![codecov][codecov-badge]][codecov-url]
[![License: MIT][license-badge]][license-url]

A lightweight Julia package for simulating the dynamics of a configuration of trapped ions interacting with laser light. 

A simple tool, built on top of [QuantumOptics.jl](https://qojulia.org/), for simulating the dynamics of a configuration of
trapped ions interacting with laser light.

**IonSim.jl** primarily performs two jobs:
1. Keeps track of the physical parameters necessary for describing the system.
2. Using these parameters, constructs a function that quickly computes the system's Hamiltonian as a function of time.

The functional form of the Hamiltonian can then be used as input to any of the solvers implemented in
[`QuantumOptics.timeevolution`](https://qojulia.org/documentation/timeevolution/timeevolution/).

Also **Ionsim.jl** is:
+ Fast: runtimes comparable to QuTiP (Cython)
+ Intuitive: you set up your simulation the same way that you set up your experiments
+ Flexible: full control over RWA cutoff frequencies, Lamb-Dicke order approximations, Hilbert space truncation, and methods
+ Open Source: all source code is freely available and built with extensibility in mind.


 <!-- For more information see:

+ Main code: [https://github.com/HaeffnerLab/IonSim.jl/tree/master/src](https://github.com/HaeffnerLab/IonSim.jl/tree/master/src)
+ Documentation: [https://docs.ionsim.org](https://docs.ionsim.org)
+ Examples: [https://examples.ionsim.org](https://examples.ionsim.org) -->

## Getting Started 

- Installation of IonSim can be done by following instructions [here](https://examples.ionsim.org/tutorial-notebooks/Installation.html).
- This page contains a collection of interactive [Jupyter notebooks](https://examples.ionsim.org/tutorial-notebooks/getting_started.html) intended to familiarize you with using `IonSim.jl`.
- Documentation of IonSim is presented [here](https://docs.ionsim.org)


```{note}
Where available, interactive notebooks can be launched by clicking on the <i class="fas fa-rocket"></i> icon. This allows you to launch a Jupyter session in the cloud using [BinderHub](https://mybinder.org/). Because Binder first needs to build the appropriate environment (inlcuding all dependencies), it can take a while to launch initially. Computational resources are also limited in the free service and we are not thoroughly maintaining compatability. But perhaps you'll still find this feature useful.
```

<!-- ```{note}
IonSim is maintained by Hartmut Haeffner's trapped ion group at UC Berkeley.
``` -->

## Citation

+ IonSim is maintained by Hartmut Haeffner's trapped ion group at UC Berkeley.
[Paper Link]()

## To Be Implemented

If you have an idea for how to improve IonSim, need some help getting things working or have any other IonSim-related questions feel free to open a GitHub issue.


If you have any questions, please make a GitHub issue.

```{note}
This webpage is built dynamically using [Jupyter Book](https://jupyterbook.org/). At build time, notebooks are executed remotely on a GitHub server. Julia's `@time` macro is used in various places to give an idea of the runtime for snippets of code. However, because of our lack of control over the execution environment, these values should be taken with a grain of salt. For example, solving for the evolution corresponding to the [VAET Hamiltonian](vaet-simulation-solve) runs as much as six times slower as compared to a standard MacBook Pro.
```

[license-url]: https://github.com/HaeffnerLab/IonSim.jl/blob/master/LICENSE.md
[license-badge]: https://img.shields.io/badge/License-MIT-green.svg

[codecov-url]: https://codecov.io/gh/HaeffnerLab/IonSim.jl
[codecov-badge]: https://codecov.io/gh/HaeffnerLab/IonSim.jl/branch/master/graph/badge.svg

[twitter-url]: https://twitter.com/Berkeley_ions
[twitter-badge]: https://img.shields.io/twitter/follow/Berkeley_ions.svg?style=social&label=@Berkeley_ions

[logo-url]: https://github.com/HaeffnerLab/IonSim.jl/blob/media/smallest_logo.png?raw=true




