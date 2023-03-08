# IonSim.jl
============================

<p align="center">
  <img src="https://github.com/HaeffnerLab/IonSim.jl/blob/media/logo3_SM.svg?raw=true", width="450px">
</p>

[![test status](https://github.com/HaeffnerLab/IonSim.jl/actions/workflows/test.yml/badge.svg)](https://github.com/HaeffnerLab/IonSim.jl/actions/workflows/test.yml)
[![codecov][codecov-badge]][codecov-url]
[![License: MIT][license-badge]][license-url]

## Getting Started 

- Installation of IonSim can be done by following instructions [here](https://examples.ionsim.org/tutorial-notebooks/Installation.html).
- This page contains a collection of interactive [Jupyter notebooks](https://examples.ionsim.org/tutorial-notebooks/getting_started.html) intended to familiarize you with using `IonSim.jl`.
- Documentation of IonSim can be found at [https://docs.ionsim.org](https://docs.ionsim.org)
- The source code of IonSim.jl can be found at [https://github.com/HaeffnerLab/IonSim.jl](https://github.com/HaeffnerLab/IonSim.jl) 

 <!-- For more information see:

+ Main code: [https://github.com/HaeffnerLab/IonSim.jl/tree/master/src](https://github.com/HaeffnerLab/IonSim.jl/tree/master/src)
+ Documentation: [https://docs.ionsim.org](https://docs.ionsim.org)
+ Examples: [https://examples.ionsim.org](https://examples.ionsim.org) -->

## Overview

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


```{note}
Where available, interactive notebooks can be launched by clicking on the <i class="fas fa-rocket"></i> icon. This allows you to launch a Jupyter session in the cloud using [BinderHub](https://mybinder.org/). Because Binder first needs to build the appropriate environment (inlcuding all dependencies), it can take a while to launch initially. Computational resources are also limited in the free service and we are not thoroughly maintaining compatability. But perhaps you'll still find this feature useful.
```

<!-- ```{note}
IonSim is maintained by Hartmut Haeffner's trapped ion group at UC Berkeley.
``` -->

## Citation

+ IonSim is maintained by [Hartmut Haeffner's trapped ion group](https://ions.berkeley.edu/) at UC Berkeley.

If you employ IonSim in your research, please support its continued development and maintenance. Use of scqubits in research publications is appropriately acknowledged by citing: (Rewrite sentence later)

[Paper Link (To Be Implemented)](https://www.youtube.com/watch?v=dQw4w9WgXcQ)

## To Be Implemented

+ Heterogeneous ion chains
+ Magnetic energy transitions
+ Raman transitions
+ Integrated tools for simulating noisy systems
+ Performance improvements
+ "Pre-canned" incorporation of some topical tools, e.g. optimized pulse-shaping for fast, high-fidelity entanglement gates 
+ A Bloch-Redfield equation solver that is compatible with time-dependent Hamiltonians
+ Incorporation of non-linear couplings between vibrational modes in a linear chain
+ Visualization tools

If you’d like to contribute to IonSim.jl, head over to our GitHub page.

+ If you have a good idea of what you’d like to do and how to do it, the preferred method is to submit a pull request on GitHub.
+ If you’re less sure about your ideas, would like some feedback, or want to open up a discussion, then feel free to either open an issue on GitHub or get in contact with any of us directly at [here (To Be Implemented)](https://www.youtube.com/watch?v=dQw4w9WgXcQ).
+ If you have an idea for how to improve IonSim, need some help getting things working or have any other IonSim-related questions feel free to open a GitHub issue.


[license-url]: https://github.com/HaeffnerLab/IonSim.jl/blob/master/LICENSE.md
[license-badge]: https://img.shields.io/badge/License-MIT-green.svg

[codecov-url]: https://codecov.io/gh/HaeffnerLab/IonSim.jl
[codecov-badge]: https://codecov.io/gh/HaeffnerLab/IonSim.jl/branch/master/graph/badge.svg

[twitter-url]: https://twitter.com/Berkeley_ions
[twitter-badge]: https://img.shields.io/twitter/follow/Berkeley_ions.svg?style=social&label=@Berkeley_ions

[logo-url]: https://github.com/HaeffnerLab/IonSim.jl/blob/media/smallest_logo.png?raw=true




