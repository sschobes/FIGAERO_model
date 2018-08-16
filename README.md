

## FIGAERO model

### Purpose

This is a model framework for simulating:
1. the molecular-level evaporation of (single compounds of) a spherical aerosol particle in a heated clean flow of nitrogen, optionally including
   - the dissociation and/or formation of physical entities, themselves assumed non-volatile (e.g. low-volatility oligomers), but able to release and/or take up molecules able to evaporate, and
2. the subsequent ad- and desorption of the evaporated molecules to surfaces.

The underlying idea is that modeling these processes achieves a reproduction of the signal measured for individual compositions during the thermal desorption of aerosol particles using a Filter Inlet for Gases and AEROsols coupled to a Time-of-Flight Chemical Ionization Mass Spectrometer (FIGAERO-CIMS, see [Lopez-Hilfiker et al., 2014](https://www.atmos-meas-tech.net/7/983/2014/amt-7-983-2014.html)).
Particle measurements using the FIGAERO-CIMS consist of first collecting aerosol particles on a PTFE filter, followed by their thermal desorption using the controlled heat of a clean flow of nitrogen. Molecules evaporating from the particles travel through the filter and into the CIMS detector, where they are ionized and their exact mass measured, which determines their molecular formulas. Typically, the desorption temperature is ramped linearly with time, yielding so-called thermograms for each composition, i.e. signal vs. desorption temperature. Features of these thermograms have previously been related to molecular and hence aerosol properties. The aim of our model here is to improve the scientific interpretation of these data.
The model and its application to organic aerosol data are thoroughly described in [Schobesberger et al., (2018)](https://www.atmos-chem-phys-discuss.net/acp-2018-398/).

### Commitment history

The initially committed version of this model is based on the internal (pre-release) version v2j3, which was the one actually used for preparing [our paper](https://www.atmos-chem-phys-discuss.net/acp-2018-398/). For commitment, that version has been:
- cleaned up (removed derelict code that has not been in use anymore)
- a bit simplified (removed a few orphaned model options that have been of limited usefulness and stopped being maintained during model development)
- better commented

However, the model has not been thoroughly tested in this shape, so bugs are currently likely. Apologies in advance!

### Model content

The model consists of .m-files containing scripts to be run using [MATLAB](https://www.mathworks.com/products/matlab.html). The code was developed using version R2016b. MATLAB a proprietary software, but the model is open-source, obviously, and fairly simple, so adaptation to other high-level programming languages such as R or Python could be feasible.

To run the model, open **_FIGAEROmodel_MAIN.m_**.
For an initial demonstration (using included example experimental data), run its script from the beginning to the end. A description of the example data is provided in the corresponding [data folder](../example_exp_data).

#### Description of individual functions:

**_FIGAEROmodel_MAIN.m_** is the main script, where: (1) experimental data is selected, loaded, and prepared for comparison to simulation results, (2) general model parameters are chosen, (3) compound-specific parameters are chosen, and (4) the model is run by calling _thermogram_simul.m_

_prepare_paras_for_simul.m_ is a collection of scripts executed by _FIGAEROmodel_MAIN.m_ for preparing some parameters for the eventual model runs, getting additional variables, and putting variables into the structure _params_ (if not included yet), which will be passed on to _thermogram_simul.m_ (the function running the simulation).

_bundle_solutions.m_ is a script called by _FIGAEROmodel_MAIN.m_ to bundle all parameters that are typically being optimized for reproducing the observations into the structure _sols_.

_get_goal.m_ is a function called by _FIGAEROmodel_MAIN.m_, outputting the variables _goal_ (= experiental data, synchronized to the expected model output) and _shift_ (= shift in seconds needed for that sychronization).

**_thermogram_simul.m_** is the major function called by _FIGAEROmodel_MAIN.m_. It prepares the initial conditions plus some parameters for the ODE solver, calls it, and gathers and cleans up its output. This function's output is the "model results".

**_thermogram_simul_ODE.m_** is the major function called by _thermogram_simul.m_. It defines the ODE set and calculates temperature- or time-dependent parameters.

_diffusivity_H2O.m_ is a function that calculates the diffusivity of water vapor in air (m^2/s).

_make_rough_plots.m_ is a function containing scripts to create 5 to 8 plots to quickly check out the model results. It also calls the scripts in _N_vs_t_plots.m_. (Both these files are currently poorly documented.)

### Note about tofTools
This model package includes functions from tofTools (version 6.08), developed by Heikki Junninen, Gustaf Lönn, and others (University of Helsinki, see for example [Junninen et al., 2010](https://doi.org/10.5194/amt-3-1039-2010)). TofTools is licensed under [GNU GPL v.3](../model_scripts/from_tofTools_608/license.txt), permitting open-source use such as here. The full & latest tofTools package (not needed for this model) is available for [free download](http://junninen.net/tofTools/).
