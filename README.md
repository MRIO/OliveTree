# OliveTree

A biophysically detailed neuronal network model of the inferior olivary nucleus of the cerebellar complex in parallel MATLAB with gpu support via CUDA. It is designed to test interactions between excitatory and inhibitory inputs on ongoing oscillatory activity.

## Model Features:

- Multiple possible connectivity schemes including clusterization and based on reconstruction data. Particularly, a cluster bridge architecture is available.
- A cell property randomizer (conductances and any other parameters) and analysis (```profile_sim.m```, ```NDscatter.m```)
- An analysis suite including spike analysis, spectral clustering, cross correlations, raster plots, etc.

# Reproducing Manuscript Figures

- Negrello et al (Plos Comp Bio, 2019)

  - ```Warnaar_figures.m```

- Loyola et al (eLife 2023)

  - Figure 7:

    - Main script – ```Loyola_PRC.m```
    - The types of stimulation must be selected by toggling the types of experiment.
  
  - Figure 8 and 9
  
    - Exc_Inh_Reconstruction.m
  
  - Supplementary Movies:
  
    - Phase differences: ```Exc_Inh_Reconstruction.m``` (set variable ``plotvolume = 1``)
    - Membrane potential: ```Exc_Inh_Reconstruction.m```
  
  - Supplementary on Tuning and analysis of cell population: 
    - ```cells_pspace.m```
  



# Main Functions and Scripts

```IOcell.m``` - simulates a single cell

```IOnet.m``` - simulates a network

```singlesim.m```- creates and simulates an example network

```CreateDefaultNeurons.m``` randomizes neuronal parameters

```profile_sim.m``` - summarizes activity and dynamics

```hilbert_of_membrane_potential.m``` - computes phase transforms of subthreshold oscillations using Damoco toolbox



# References based on this model:

[1]	M. Negrello, P. Warnaar, V. Romano, C. B. Owens, S. Lindeman, E. Iavarone, J. K. Spanke, L. W. J. Bosman, and C. I. De Zeeuw, “Quasiperiodic rhythms of the inferior olive.,” PLoS Computational Biololgy, vol. 15, no. 5, p. e1006475, May 2019.


[2] Sebastián Loyola, Tycho M. Hoogland, Mario Negrello, Hugo Hoedemaker, VincenzoRomano, Chris I. De Zeeuw. "How inhibitory and excitatory inputs gate output of the inferior olive". bioRxiv 2022.09.04.506491; doi: https://doi.org/10.1101/2022.09.04.506491  (This article is a preprint and has not been certified by peer review)
TAG: 06Apr2023

## Foundations:

[1] Y. Manor, J. Rinzel, I. Segev, and Y. Yarom. Low-amplitude oscillations in the inferior olive: a model based on electrical coupling of neurons with heterogeneous channel densities. Journal of Neurophysiology, 77(5):2736 – 2752, 1997.

[2] N. Schweighofer, K. Doya, and M. Kawato. Electrophysiological properties of inferior olive neurons: A compartmental model. Journal of Neurophysiology, 82(2):804 – 817, 1999.

[3] P. Bazzigaluppi, J. R. D. Gruijl, R. S. v. d. Giessen, S. Khosrovani, C. I. D. Zeeuw, M. T. G. d. Jeu, and bazzigali. Olivary subthreshold oscillations and burst activity revisited. Frontiers in Neural Circuits, 6:91, 2012.

[4] J. D. Gruijl, P. Sokol, M. Negrello, and C. I. D. Zeeuw. Calcium Dependent Gap Junction Plasticity: Modulation of Electrotonic Coupling in the Inferior Olive Glomerulus. bioRxiv, 17:072041, 2016.

[5] J. R. D. Gruijl, P. Bazzigaluppi, M. T. G. d. Jeu, and C. I. D. Zeeuw. Climbing fiber burst size and olivary sub-threshold oscillations in a network setting. PLoS Computational Biololgy, 8(12):e1002814, 2012.



# Data Note:

Cell coordinates for the reconstruction of the Medial Accessory Olive originates from Marylka Uusisaari, Nora Vrieler and Nicholas Medvedev.

All rights to the reconstruction reserved to Marylka Uusisaari, Nora Vrieler.


# CONTRIBUTORS
 - Mario Negrello (main author)
 - Sungho Hong (spectral clustering)
 - Marylka Uusisaari, Nora Vrieler, Nick Medvedev - Anatomical Data
 - Sebastian Loyola -- electrophysiological data


# LICENCE

<a rel="license" href="http://creativecommons.org/licenses/by-sa/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by-sa/4.0/88x31.png" /></a>

<span xmlns:dct="http://purl.org/dc/terms/" property="dct:title">Inferior Olivary Model</span> by <a xmlns:cc="http://creativecommons.org/ns#" href="https://github.com/MRIO/OliveTree" property="cc:attributionName" rel="cc:attributionURL">Mario Negrello</a> and collaborators is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by-sa/4.0/">Creative Commons Attribution-ShareAlike 4.0 International License</a>.Based on work at <a xmlns:dct="http://purl.org/dc/terms/" href="https://github.com/MRIO/OliveTree" rel="dct:source">https://github.com/MRIO/OliveTree</a>.



Damoco toolbox (licenses and publications: http://www.stat.physik.uni-potsdam.de/~mros/damoco.html)



