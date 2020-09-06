# IonSim.jl Examples

A webpage for organizing IonSim.jl tutorials and examples. The webpage is built using [jupyter-book](https://www.jupyterbook.org) and you can refer to that link for details on how to edit or add content. But the gist of it is, you add Jupyter notebooks and/or markdown files to this repo and then update the `_toc.yaml` file to let jupyter-book know where to display this content.

A GitHub action has been created (workflow description located at `.github/workflows/book.yml`) such that anytime something is pushed to the master branch, the webpage will be rebuilt and pushed to the gh-pages branch (this branch should not be edited directly). During this rebuild any new or edited Jupyter notebooks will be run from scratch. This process is explained in the documentation for jupyter-book and detailed GitHub actions documentation can be found [here](https://docs.github.com/en/actions). 

*Note: you may need to edit the `book.yml` file to make sure the proper Julia packages and/or versions are loaded into the remote execution environment.*

*Note: after building, also please double check that the `custom domain` is set to `examples.ionsim.org` in the repo settings (and `enforce HTTPS` is selected)*

*Note: when adding a notebook, in order for interactivity to work, you'll need to make sure that any dependencies are reflected in the current state of the `Project.toml` file.*

The webpage is published to [https://examples.ionsim.org/intro.html](https://examples.ionsim.org).
