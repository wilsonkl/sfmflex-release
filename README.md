# Visualizing Spectral Bundle Adjustment Uncertainty
Kyle Wilson and Scott Wehrwein, 3DV 2020.

## Abstract
Bundle adjustment is the gold standard for refining so- lutions to geometric computer vision problems. This paper develops an uncertainty visualization technique for bundle adjustment solutions to Structure from Motion problems.

Propagating uncertainty through an optimization—from measurement uncertainties to uncertainties in the resulting parameter estimates—is well understood. However, the cal- culations involved fail numerically for real problems. Often we cope by considering only individual variances, but this ignores the important mutual dependencies between param- eters. The dominant modes of uncertainty in most models are large motions involving nearly all parameters at once. These frequently look like flexions, stretchings, and bend- ings in the overall scene structure.

![Teaser Figure](teaser_fig/pdf)

In this paper we present a numerically tractable method for computing dominant eigenvectors of the covariance of a Bundle Adjustment solution. We pay careful attention to the mismatched scales of rotational and translational param- eters. Finally, we animate this spectral information. The resulting interactive visualizations (included in the supple- mental) give insight into the quality and failure modes of a model. We hope that this work is a step towards broader uncertainty-aware computation for Structure from Motion.


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
