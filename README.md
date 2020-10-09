# Visualizing Spectral Bundle Adjustment Uncertainty
Kyle Wilson and Scott Wehrwein, 3DV 2020.

## Abstract
Bundle adjustment is the gold standard for refining solutions to geometric computer vision problems. This paper develops an uncertainty visualization technique for bundle adjustment solutions to Structure from Motion problems.

Propagating uncertainty through an optimization–from measurement uncertainties to uncertainties in the resulting parameter estimates–is well understood. However, the calculations involved fail numerically for real problems. Often we cope by considering only individual variances, but this ignores the important mutual dependencies between parameters. The dominant modes of uncertainty in most models are large motions involving nearly all parameters at once. These frequently look like flexions, stretchings, and bendings in the overall scene structure.

In this paper we present a numerically tractable method for computing dominant eigenvectors of the covariance of a Bundle Adjustment solution. We pay careful attention to the mismatched scales of rotational and translational parameters. Finally, we animate this spectral information. The resulting interactive visualizations (included in the supplemental) give insight into the quality and failure modes of a model. We hope that this work is a step towards broader uncertainty-aware computation for Structure from Motion.

![Teaser Figure](teaser_fig.png)

## Overview
This repository contains:
- Interactive web visualizations that accompany our 3DV20 paper
- The research code that produces these visualizations

## Visualizations
These visualizations are included in the supplemental to our paper. We also post them here. Please see the supplemental commentary for interpretive remarks. 

### How to run the demos
Please clone this repository. There are HTML files for each scene in the `visualizations` directory. Open these in any web browser (tested: Firefox, Chrome, Safari). 

### Demo controls

Action | Control
-------|------------
Rotate the model | `click + drag`
Translate the model | `altclick + drag`
Next/previous eigenvector | `j/k`
Increase/decrease translation scale | `h/l`
Increase/decrease rotation scale | `g/;`
Help popup | `?`

## Research Code
The code pipeline works as follows:
- Given a problem instance in the format of [Bundle Adjustment in the Large](https://grail.cs.washington.edu/projects/bal/)
- Run a stock bundle adjuster to guarantee that we are at a local minimum (otherwise the Gauss-Newton approximation will be terrible!)
- Extract the problem parameters and its derivatives and save them out to an HDF5 file.
- Read the HDF5 file, compute eigenvectors, and save the eigenvectors and scene parameters to a JSON file.
- Load the JSON file and display an interactive visualization.

The preprocessing stages are based on [ceres-solver](http://ceres-solver.org/) and its included bundle adjuster. This is in C++. The eigenvector computations are performed in Julia, and the final interactive visualization is based on ThreeJS in Javascript.

### Setup and Building
The `C++` targets depend on `Eigen`, `GFlags`, `Glog`, and `hdf5`. They also require `ceres-solver`, which should be [built from source](http://ceres-solver.org/installation.html). Be sure to build ceres with `suite-sparse`.

For convenience on Ubuntu: `apt-get install libeigen3-dev libgoogle-glog-dev libgflags-dev libhdf5-dev`

To build the `C++` targets:
```
mkdir build
cd build
cmake ..
make
```

The Julia code also has many dependencies. Install all of these by running `julia deps.jl` from the `code/` directory.

### Scripts and bins
Run all of these scripts from the `code/` directory.

**Preprocessing scripts:**
+ `build/bin/compute_jacobian <bal_problem> <jac_file>` : compute the Jacobian of a bundle adjustment problem.
+ `build/bin/bundle_adjuster` : an all-purpose bundle adjuster, lightly tweaked from the ceres-solver [example](https://github.com/ceres-solver/ceres-solver/blob/master/examples/bundle_adjuster.cc). There are many command line options.

**Dataset scripts:**
+ `julia scripts/download_bal.jl` : scape, download, and decompress all of the BAL datasets. Store them at `dataset/bal/`.
+ `julia scripts/run_ba_on_bal.jl` : run a bundle adjuster on all BAL problems found. This takes a long time.