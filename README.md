# OliveTree

A biophysically detailed neuronal network model of the inferior olivary nucleusof the cerebellar complex in parallel MATLAB with gpu support via CUDA.

## Model Features:

- Multiple possible connectivity schemes including clusterization and based on reconstruction data. Particularly, a cluster bridge architecture is available.
- A cell property randomizer (conductances and any other parameters).
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
  
    - Phase differences: Exc_Inh_Reconstruction.m (set variable ``plotvolume = 1``)
    - Membrane potential: Exc_Inh_Reconstruction.m
  
  - Supplementary on Tuning and analysis of cell population: cells_pspace.m
  
    - 
  
  - 
  
  - 
  
    



# Functions

```IOcell.m``` - simulates a single cell

```IOnet.m``` - simulates a network

```singlesim.m```- creates and simulates an example network

```CreateDefaultNeurons.m``` randomizes neuronal parameters

```profile_sim.m``` - summarizes activity and dynamics

```hilbert_of_membrane_potential.m``` -



# Reference:

[1]	M. Negrello, P. Warnaar, V. Romano, C. B. Owens, S. Lindeman, E. Iavarone, J. K. Spanke, L. W. J. Bosman, and C. I. De Zeeuw, “Quasiperiodic rhythms of the inferior olive.,” PLoS Computational Biololgy, vol. 15, no. 5, p. e1006475, May 2019.



[2] Sebastián Loyola, Tycho M. Hoogland, Mario Negrello, Hugo Hoedemaker, VincenzoRomano, Chris I. De Zeeuw. "How inhibitory and excitatory inputs gate output of the inferior olive". bioRxiv 2022.09.04.506491; doi: https://doi.org/10.1101/2022.09.04.506491  (This article is a preprint and has not been certified by peer review)

# Data Note:

Cell coordinates for the reconstruction of the Medial Accessory Olive originates from Marylka Uusisaari and Nora Vrieler.

All rights to the reconstruction reserved to Marylka Uusisaari, Nora Vrieler and Nick ….



# LICENCE

<a rel="license" href="http://creativecommons.org/licenses/by-sa/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by-sa/4.0/88x31.png" /></a>

<span xmlns:dct="http://purl.org/dc/terms/" property="dct:title">Inferior Olivary Model</span> by <a xmlns:cc="http://creativecommons.org/ns#" href="https://github.com/MRIO/OliveTree" property="cc:attributionName" rel="cc:attributionURL">Mario Negrello</a> and collaborators is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by-sa/4.0/">Creative Commons Attribution-ShareAlike 4.0 International License</a>.Based on work at <a xmlns:dct="http://purl.org/dc/terms/" href="https://github.com/MRIO/OliveTree" rel="dct:source">https://github.com/MRIO/OliveTree</a>.



Damoco toolbox (licenses and publications: http://www.stat.physik.uni-potsdam.de/~mros/damoco.html)



