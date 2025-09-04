# multivariateEGPD
This repository contains the source code and replication files for ‘Joint modeling of low and high extremes using the multivariate eGPD,’ including simulation studies, simulation-based inference with the neural ‘NeuralEstimators’ package, and an application to Dutch precipitation extremes.


### Software dependencies

We suggest that users set up a [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html) environment, so that the dependencies of this repository do not affect the user's current installation. In your environment, install `Python`, `R`, and `Julia`. An environment instruction that worked in our case was the following:

```
conda create -n mEPGD -c conda-forge r-base nlopt
```

Then activate the conda environment with:

```
conda activate mEPGD
```

If you do not wish to use a conda environment, then install the software directly from the following websites:

- Install [R >= 4.4.0](https://www.r-project.org/).

Once `R` is setup, install their package dependencies by running the following commands from the top-level of the repository: 

```
Rscript dependencies_install.R
```

Finally, install [Julia >= 1.11.0](https://julialang.org/downloads/), for example: 

```
curl -fsSL https://install.julialang.org | sh
juliaup add 1.11.3
juliaup default 1.11.3
```

Then install the Julia package dependencies by running:

```
julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'
```

### Hardware requirements

In general, the fast training of neural networks requires GPUs. However, the code in this repository also runs on CPUs in a moderate amount of time. Therefore, there are no major hardware requirements for running this code. 

### Reproducing the results

First, download this repository and navigate to its top-level directory within terminal.

The repository is organised into folders containing source code (`src`), intermediate objects generated from the source code (`intermediates`), and figures (`img`). 

The replication script is `run.sh`, invoked using `bash run.sh` from the top level of this repository. The replication script will automatically train the neural networks, generate estimates/samples from both the neural and likelihood-based estimators/samplers, and populate the `img` folder with the figures and results of the manuscript.

Note that the nature of our experiments means that the run time for reproducing the results of the manuscript can be moderate (on the order of several hours). 

### Minor reproducibility difficulties

When training neural networks, there is often unavoidable non-determinism: see, for example, [here](https://discourse.julialang.org/t/flux-reproducibility-of-gpu-experiments/62092). In our reproducible code, this does not significantly affect the "story" of the final results in the sense that each method performs similarly well in each run, but there may be some slight numerical differences each time the code is executed.

