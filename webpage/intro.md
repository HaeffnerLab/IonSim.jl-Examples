Introduction
============================

This page contains a collection of interactive [Jupyter notebooks](https://jupyter.org/) intended to familiarize you with using `IonSim.jl`. These notebooks are organized into several sections:
<br> 

- **Tutorials** and **Examples**: A collection of walkthroughs and examples designed to demonstrate a particular concept. 
<br><br>
- **Learn**: A pedagogical introduction to the core concepts in Trapped ions (state preparation, single-ion operations, entanglement-operations).

You can return to the main page [here](https://www.ionsim.org) (or by clicking on the logo in the sidebar).

```{note}
Where available, interactive notebooks can be launched by clicking on the <i class="fas fa-rocket"></i> icon. This allows you to launch a Jupyter session in the cloud using [BinderHub](https://mybinder.org/). Because Binder first needs to build the appropriate environment (inlcuding all dependencies), it can take a while to launch initially. Computational resources are also limited in the free service and we are not thoroughly maintaining compatability. But perhaps you'll still find this feature useful.
```

```{note}
This webpage is built dynamically using [Jupyter Book](https://jupyterbook.org/). At build time, notebooks are executed remotely on a GitHub server. Julia's `@time` macro is used in various places to give an idea of the runtime for snippets of code. However, because of our lack of control over the execution environment, these values should be taken with a grain of salt. For example, solving for the evolution corresponding to the [VAET Hamiltonian](vaet-simulation-solve) runs as much as six times slower as compared to a standard MacBook Pro. 
```


