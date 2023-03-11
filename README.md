# OliveTree

A biophysically detailed neuronal network model of the inferior olivary nucleusof the cerebellar complex in parallel MATLAB with gpu support via CUDA.

## Model Features:

- Multiple possible connectivity schemes including clusterization and based on reconstruction data. Particularly, a cluster bridge architecture is available.
- A cell property randomizer (conductances and any other parameters).
- An analysis suite including spike analysis, spectral clustering, cross correlations, raster plots, etc.

# Reproducing Manuscripts

- Negrello et al (Plos Comp Bio, 2019)

  - Warnaar_

- Loyola et al (eLife 2023)

  - Figure 7:

    - Run "Loyola_PRC.m"

  - Figure 8: Supplementary Materials

    - Exc_Inh_Reconstruction.m

  - Tuning: cells_pspace.m

  - 
  
  - 
  
    



# Functions

IOcell.m

IOnet.m

CreateDefaultNeurons.m

profile_sim.m





# Reference:

[1]	M. Negrello, P. Warnaar, V. Romano, C. B. Owens, S. Lindeman, E. Iavarone, J. K. Spanke, L. W. J. Bosman, and C. I. De Zeeuw, “Quasiperiodic rhythms of the inferior olive.,” PLoS Computational Biololgy, vol. 15, no. 5, p. e1006475, May 2019.
