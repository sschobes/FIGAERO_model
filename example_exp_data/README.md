## Data contents


### PSD.mat

MATLAB data file

**Content:**
- _sizeDistFrac_
SOA particle number size distribution
  - col. 1: size bin (nm)
  - col. 2: number size distribution dN/dlopDp (cm^-3)
- _sizeDistOrgVolMed_
Median diameter of SOA particle volume size distribution (nm)

**Source:**
Processed from SMPS measurements, representing particle size distribution in the PNNL Environmental Chamber during SOAFFEE campaign, for time 11 Jul 2015, 17:04 UTC


### FIGAERO_exp.mat

MATLAB data file

**Content:**
_FIGAERO_expdat_
Structure containing:
- _expName_
Name of experiment
- _C8H12O5_
Thermogram data for ion C8H12O5.I-
  - _measdat_
    - col. 1: desorption temperature (°C)
    - col. 2: normalized and corrected signal (s^-1)
  - _frac_
    - col. 1: fraction of total signal from organics
    - col. 2: fraction of signal remaining since start of desorption
- _allCHO_
Same as C8H12O5, but for all observed CxHyOz compounds ("bulk")
- _timres_
Time resolution of measdat (i.e. time passed between each row) (s)
- _wait_time_
Period of isothermal evaporation ("wait") prior to start of temperature ramp, typically zero (min)
  
**Source:**
Processed from FIGAERO (University of Washington FIGAERO/TOF1) measurements during SOAFFEE, 11 Jul 2015, 17:04 UTC