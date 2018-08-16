%% FIGAERO model, MAIN SCRIPT

% FIGAEROmodel_MAIN.m

%% Version notes:

% v1.00:
% based on pre-release version v2j3, but
% - cleaned up
% - a bit simplified
% - better commented

% Siegfried Schobesberger
% 16 Aug 2018
% Kuopio, Finland

%% Note about tofTools:

% This model package includes functions from tofTools (version 6.08), 
% created by Heikki Junninen, Gustaf L?nn, et al. (University of Helsinki, 
% see for instance https://doi.org/10.5194/amt-3-1039-2010)

% tofTools is licensed under GNU GPL v.3 (see 
% /model scripts/from_tofTools_608/license.txt), permitting open-source 
% use, as is the case here.

% The full & latest tofTools package is available at 
% http://junninen.net/tofTools/

%% Initializing (paths, tofTools)
% Hit Enter in pop-up to confirm tofTools folder

% Get and go to path this script is in
editor_service = com.mathworks.mlservices.MLEditorServices;
editor_app = editor_service.getEditorApplication;
active_editor = editor_app.getActiveEditor;
storage_location = active_editor.getStorageLocation;
file = char(storage_location.getFile);
homepath = [fileparts(file),'/']; clear active_* editor_* file storage_*
cd(homepath);

% Add home directory paths
addpath(genpath(homepath));

%% Load experimental data into workspace
% Using size distribution data to define particle size(s)
% Using FIGAERO data if option Raoult==2 is set

% example exp data:
% (for explanation, see 'example_exp_data/README.md')

% size distribution data (processed from SMPS measurements):
load([homepath,'example_exp_data/PSD.mat']); % loads 2 variables

% FIGAERO data (processed selection)
load([homepath,'example_exp_data/FIGAERO_exp.mat']); % loads 1 variable

%% Select experimental FIGAERO data used in model run
% Needed if Rest==1 is going to be used, for initial amount of selected composition
% Needed even more if Raoult==2 is going to be used (for obtaining Raoult terms)

% Compound selection
c = 'C8H12O5';

% Extract experimental FIGAERO data
    % measdat col. 1 = desorption T
    % measdat col. 2 = thermogram signal
    % fracs col. 1 = fraction of bulk
    % fracs col. 2 = fraction remaining since start
measdat = FIGAERO_expdat.(c).measdat; % measurement data for compound to be simulated ("compound i")
measdat = [nan(size(measdat,1),1) measdat]; % add a column of NaNs because time trace not available (if available, needs to be in s)
chi_exp = FIGAERO_expdat.(c).fracs(:,1); % chi (Raoult term)
rem_exp = FIGAERO_expdat.(c).fracs(:,2); % fraction remaining
measdatA = FIGAERO_expdat.allCHO.measdat; % measurement data for all CxHyOz ("bulk")
measdatA = [nan(size(measdatA,1),1) measdatA]; % add a column of NaNs because time trace not available (if available, needs to be in s)
remA_exp = FIGAERO_expdat.allCHO.fracs(:,2); % phi (fraction remaining of bulk)

% Experiment info:
expName = FIGAERO_expdat.expName; % experiment name
timres = FIGAERO_expdat.timres; % time resolution of measdat and measdatA (s)
wait_time = FIGAERO_expdat.wait_time; % isothermal evaporation wait time (min)
    
%% Choose model parameters

% General model choices
WaitInModel = 1;
    % 1 ... simulating isothermal ("wait") period prior to start of T ramp (if that time is > 0)
    % 0 ... ignoring such a period (if existing)
params.Soak = 1;
    % 1 ... for also simulating "soak period" at end of T ramp (default)
    % 0 ... T ramp only
TempTreat = 'Gauss'; % 'Simple' or 'Gauss'
    % 'Simple' ... single temperature (much faster, but leading to too narrow peaks)
    % 'Gauss' ... non-ideal heating via Gauss-like temperature distribution
SizeTreat = 1;
    % 0 ... set a single size manually below
    % 1 ... for using single diameter, i.e. sizeDistOrgVolMed (default)
    % 2 ... for using size distribution, i.e. sizeDistFrac (SLOWER)

% Desorption temperature profile
% ... should correspond to measurement data, obviously
params.T0 = 273.15+25; %273.15+25; % T at start of desorption (K)
% !!! model not tested for T0 other than 298.15 K !!!
params.rampRate = 0.14; % linear T ramp rate (K/s)
params.durRamp = 2000; % duration of T ramp (s) ... max. des. T in model will be = T0+rampRate*durRamp
params.durSoak = 2940; % duration of soak period at max. des. T (s)
params.maxT = params.T0+params.rampRate*params.durRamp;

% Miscellaneous parameters
params.R = 8.31446; % universal gas constant (J mol^-1 K^-1)
params.pressure = 1013.25; % ambient pressure (mbar)

% Model options
Raoult = 1;
    % 1 ... chi and phi are calculated explicitly during model run (default and faster)
    %   ... (=1 if only one compound modeled)
    % 2 ... chi and phi are interpolated from experimental data, using chi_exp
    %   ... (implicitly assumes that thermogram signal immediately relates to particle evaporation, i.e. with negligible transport delay)
    %   ... (if chi_exp is giving mole fraction, then that is used instead of mass fraction in this case!)
Walls = 1;
    % 1 ... simulating a single set of walls (default)
    % 2 ... two sets of walls (e.g. filter + IMR?)
Glu = 1;
    % 1 ... simulate reversible oligomerization ("gluing"), in equilibrium with monomers at start!
    % 0 ... no oligomers
GluOrder = 3; % ...
    % 1 ... pseudo-1st-order gluing
    % 2 ... "real" 2nd-order gluing (requires "all" compounds; anyway currently NOT working!)
    % 3 ... default: 2nd-order but simplified (using OA-remaining, derived depending on params.Raoult)
GluTimDep = 0; % time-dependent gluing parameters!
    % !!! still experimental !!!
Rock = 0; % 1, 2 or 3 ... include other variant of oligomer ("rock")
    % 0 ... off
    % 1 ... formed from glued but dissociating to monomer (currently requires GluOrder == 1 or 3)
    % 2 ... working identically (i.e. in parallel) to gluing (but with possibly different values) (currently requires GluOrder == 1 or 3)
    % 3 ... non-volatile oligomers with pre-set amount (not in equilibrium at start, set by RockFrac), and thermally decomposing into monomer
Rest = 1;
    % 1 ... including a bulk component (to be determined via measdat or measdatA)
    % 0 ... or not
Reff = 1;
    % 1 ... including refractory core  (non-evaporating, and non-interacting)
    % 0 ... or not
nrComps = 1; % number of compounds to simulate, each with own set of parameters, excluding bulk (if Rest==1) and core (if Reff==1)
UseNdist = 0;
    % 0 ... default. IF nrComps>1, but initial conditions defined by chi_exp are only for 1 compound, will EQUALLY divide initial molecule number to nrComps compounds
    % 1 ... allows for using a user-defined distribution for the starting amounts of the nrComps compounds, i.e. initial molecule numbers will be distributed according to Ndist
    % Ndist (which needs size=[nrComps 1]) will be defined below when choosing the compound-specific parameters

% Wall parameters
params.k_won = 1; % wall stick parameter, wall 1 (unitless (~adsorption probability))
params.k_wonI = 0.81; % wall stick parameter, wall 2 (if in use)
params.Tfrac_after = 1; % fraction of desorption T to apply for wall 1
params.Tfrac_afterI = 1; % fraction of desorption T to apply for wall 2
params.CW = 0.877*1e4; % equivalent sorbing aerosol conc, times adsorption time-scale, wall 1 (ug m^-3 s^-1)
    % CW = 8770 is from tuning for UW-FIGAERO (TOF1). Adjust as needed to fit calibrations. (And ideally use a value yielding sensible fits on isothermal blanks.)
    % set much smaller (e.g. /1e5) to effectively turn those surface interactions "off"
params.CWI = 1.49*1e5; % equivalent sorbing aerosol conc, times adsorption time-scale, wall 2 (ug m^-3 s^-1)

% Aerosol particle properties
params.dp = 100e-9; % diameter of simulated particle (m)
    % !!! ... will be overridden if SizeTreat>0 !!!
params.dens = 1200; % density of particle (kg m^-3)
params.ReffDia = 50e-9; % diameter of refractory core (m) (only used if Reff==1; should be smaller than total particle size(s)!)
params.densReff = 1770; % kg m^-3 ... density for refractory part (only used if Reff==1)

%% Choose compound-specific parameters
% (For fitting more-or-less manually at this stage)

% initializing parameters for gluing and rocking (in case they would not get set):
% gluing:
    KG = zeros(nrComps,1);  % rate constant for gluing k_g (s^-1)
    KD = KG; % rate constant for dissociation from glued k_d (s^-1)
    EAD = KG; % activation energy for dissociation E_A,d (kJ/mol)
    EAG = KG; % activation energy for gluing E_A,g (kJ/mol)
    KGRe = KG; KDRe = KG; % Bulk ("rest") will use mean EAD and EAG for now
% rocking:
    KR = KG; % analogue to KG (s^-1)
    KRD = KG; % analogue to KD (s^-1)
    EARD = KG; % analogue to EAD (kJ/mol)
    EAR = KG; % analogue to EAG (kJ/mol)
    RockFrac = zeros(nrComps,1);  % fraction initially in low-volatility state ("Rock" in this case) ... used if Rock==3 (each element <=1)
    RockFracRe = 0; % same, but for bulk ("rest") ... used if Rock==3
    % Bulk will use mean parameters for now
% time-dependent gluing (see extra_info.txt for more information):
    KGmax = KG; KGgainT = KG; KDmin = KG; KDdecT = KG;
    
% initializing Ndist (in case it is used but will not be set)
Ndist = ones(nrComps,1);

% some example value sets more-or-less suitable for fitting the example SOA data for composition C8H12O5
if strcmp(expName,'SOAFFEE-2015-0711-1704')
    if strcmp(c,'C8H12O5')
        mff = 1; % default = 1, increase to use smoothed data for determining maximum thermogram value
        MOLW = tof_exact_mass(c,1) * ones(nrComps,1); % molecular weight (g/mol)   
        if nrComps == 1
           if Glu==1 && Rock == 0
               CSTAR = 0.8; % saturation concentration C* at T0 (ug/m3)
               DH = 88; % vaporization enthalpy DeltaH (kJ/mol)
               MA = 1; % evaporation coefficient alpha
               KG = 1.7e-3; % rate constant for gluing k_g (s^-1)
               KD = 5e-4; % rate constant for dissociation from glued k_d (s^-1)
               EAG = 0; % activation energy for gluing E_A,g (kJ/mol)
               EAD = 6; % activation energy for dissociation E_A,d (kJ/mol)
               % for getting the bulk thermogram roughly right
                   DHRe = 40;
                   CSTARRe = 6;
                   KGRe = KG; KDRe = KD;
                   MWRe = 297;
            elseif Glu==0 && Rock == 3
               MOLW = tof_exact_mass(c,1);
               CSTAR = 0.8; % ug/m3
               DH = 88; % kJ/mol
               MA = 1 * ones(nrComps,1);
                   EARD = 6; % seems with pre-exponential factor A=3e10, I'd need a range from 85 to 125
                   KRD = 5.6e-3*exp(-EARD*1000/8.314/298.15); % A now adjusted to work with EARD = 6
                   RockFrac = 0.85; % fraction initially in low-volatility state ("Rock" in this case)
           else
               MOLW = tof_exact_mass(c,1);
               CSTAR = 0.8; % ug/m3
               DH = 88; % kJ/mol
               MA = 1 * ones(nrComps,1);
           end    
        elseif nrComps==8 && Glu==0 && Rock==3
                MOLW = tof_exact_mass(c,1) *ones(nrComps,1);
                CSTAR = 2 *ones(nrComps,1); % ug/m3
                DH = 80 *ones(nrComps,1); % kJ/mol
                MA = 1 * ones(nrComps,1);
                RockFrac = ones(nrComps,1);
                    Ndist = [1 0.6 0.5 0.35 0.3 0.3 0.3 0.3];
                    EAR = zeros(nrComps,1);
                    EARD = [80; 85; 90; 95; 100; 105; 110; 115];
                    KR = zeros(nrComps,1);
                    KRD = 3e10*exp(-EARD*1000/8.31446/298.15);
                DHRe = 15;
                CSTARRe = 3;
                MWRe = 297;
        else
            disp('Lacking properly set parameters (probably)! STOPPING.')
            return
        end
    end
end

%% Adjust parameters, initialize variables for model run, get chi and phi

% Get size distribution if desired
if SizeTreat==2
    SrunNos = 2:8; % which size bins to use
    params.dp = sizeDistFrac(SrunNos,1)'*1e-9; % particle diameter(s) to use (m)
    params.sizeDistFrac = sizeDistFrac;
elseif SizeTreat==1
    SrunNos = 0; % i.e. off (will use single size)
    params.dp = sizeDistOrgVolMed*1e-9; % particle diameter to use (m)
end
if Reff && sum(params.dp<=params.ReffDia)>0
   disp('Please do not include particle sizes < or = size of refractory core!')
end

% T distribution for non-ideal heating
    % Gaussian, sort of...
    Tx = [0:0.05:1]; % tune for different FIGAERO! (UW-FIGAERO (TOF1): Tx = [0:0.05:1])
    Tsigm = 0.28; % tune for different FIGAERO! (UW-FIGAERO (TOF1): Tsigm = 0.28)
    Tnorm = normpdf(Tx,1,Tsigm); % x, mu, sigma
    params.Ty = Tnorm/sum(Tnorm); % fraction 
    params.DTrange = Tx;
% which parts of T distributions to actually use:
switch TempTreat
    case 'Gauss'
        params.TrunNos = 7:21; % (Ty very small for x<7 (UW-FIGAERO case))
    otherwise
        params.TrunNos = 21:21;
end

% a few other things (incl. the chi and phi parts)
prepare_paras_for_simul

%% Setting up optimizable parameters for model runs

% collect solutions into single structure (idea of this is to eventually set up for an optimization routine, but WIP)
bundle_solutions; % creates strucuture "sols"

% For optimizing: get experimental data to fit to into "correct" format
% -> new variable "goal":
% - col 1: time (s)
% - col 2: temperature (K)
% - col 3: signal
[goal,shift] = get_goal(measdat,timres,params);

if shift<0 && shift>-100 % second is to avoid the next line going massively wrong if wait is not modelled, but data has wait
    params.TotDur = floor(params.TotDur + shift); % shorten model run according to shift between measurement data and goal
end

%% RUN MODEL

[simulation,simNi,simNg,simNr,simNw,dy] = thermogram_simul(sols,params);

% OUTPUT:
% simulation
    % - row 1: time (s)
    % - row 2: temperature (K)
    % - row 3+: signal (each individial compound, then bulk (if Rest==1), then core (if Reff==1)
% simNi: number of free monomer molecules in particle (one row for each compound, as in simulation, rows 3+)
% simNg: number of monomers in "glued" state
% simNr: number of monomers in "rock" state
% simNw: number of monomers on walls/surfaces
% dy: output from all runs of ODE solver

%% Several plots for checking results

make_rough_plots(simulation,simNi,simNg,simNr,simNw,mff,goal,measdat,measdatA,timres,params);

