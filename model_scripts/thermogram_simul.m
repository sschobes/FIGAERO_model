function [simulation,simNi,simNg,simNr,simNw,dy_sv2] = thermogram_simul(sols,params)

% This function prepares the initial conditions and some parameters for the 
% ODE solver, calls it, and gathers and cleans up its output.
% 
% INPUT:
% sols: structure containingn variables that are typically optimized for
%       reproducing observations
% params: structure containing model parameters
% 
% OUTPUT:
% simulation
%   - row 1: time (s)
%   - row 2: temperature (K)
%   - row 3+: signal (each individial compound, then bulk (if Rest==1), then core (if Reff==1)
% simNi: number of free monomer molecules in particle (one row for each compound, as in simulation, rows 3+)
% simNg: number of monomers in "glued" state
% simNr: number of monomers in "rock" state
% simNw: number of monomers on walls/surfaces
% dy_sv2: output from all runs of ODE solver

%% "Extract" all those variable and some parameters

% compound-specific parameters
CSTAR = sols.CSTAR;
MA = sols.MA;
DH = sols.DH;
EAD = sols.EAD;
EAG = sols.EAG;
KG = sols.KG;
KD = sols.KD;
KR = sols.KR;
KRD = sols.KRD;
EARD = sols.EARD;
EAR = sols.EAR;
KGmax = sols.KGmax;
KGgainT = sols.KGgainT;
KDmin = sols.KDmin;
KDdecT = sols.KDdecT;
KGRe = sols.KGRe;
KDRe = sols.KDRe;
params.dHvap(1:params.nrComps) = DH;
params.ma(1:params.nrComps) = MA;

% other model parameters
T0 = params.T0;
maxT = params.maxT;
rampRate = params.rampRate;
DTrange = params.DTrange;
Ty = params.Ty;
TrunNos = params.TrunNos;
SrunNos = params.SrunNos;
dp = params.dp;
Rest = params.Rest;
Reff = params.Reff;
ReffDia = params.ReffDia;
if ~Reff; ReffDia=0; end
nrComps = params.nrComps;
Raoult = params.Raoult;
chi = params.chi;
Glu = params.Glu;
GluOrder = params.GluOrder;
Rock = params.Rock;
Soak = params.Soak;
durRamp = params.durRamp;
durSoak = params.durSoak;
Wait = params.Wait;
wait_time = params.wait_time;
Walls = params.Walls;
CW = params.CW;
CWI = params.CWI;
MW = params.MW;
dHvap = params.dHvap;
ma = params.ma;

%% MAIN LOOP for running the model
% in "simple" case run only once

    TR = 0;
for DT = rampRate*DTrange(TrunNos)
    disp(['Running w/ T ramp rate = ' num2str(DT,3) ' K/min (' num2str(DT/rampRate*100) '%) ...']);
    TR = TR+1; % desorption T loop number
        SR = 0;
    for DP = dp
        SR = SR+1; % particle size loop number
        params.dp = DP;
        
        % converting C*s to working units
            Csat0_ug(1:params.nrComps,1) = CSTAR; % ug m^-3
            Csat0_ug(nrComps+1,1) = params.CSTARRe;
                % for refractory part, I set dHvap to 0. Csat should then just be really low:
                Csat0_ug(dHvap==0) = 1e-10; % ug m^-3
            % and in working units (molec cm^-3)
            Csat0 = Csat0_ug ./MW *6.022e23 *1e-12; % molec cm^-3
                         
        % collecting values for "gluing" and dissociating variables
            k_glu0 = [...
                KG; ... % (molec^-1 s^-1)
                KGRe; ... % "rest"
                1e-20; ... % refractive part
                ];
            K_glu0 = sqrt(k_glu0*k_glu0'); % use geometric mean for i gluing with j if k_{g,i} is different from k_{g,j}
            k_dis0 = [...
                KD; ... % (s^-1)
                KDRe; ... % "rest"
                1e-10; ... % refractive part
                ];
            params.Ea_dis = [...
                EAD; ...  % (kJ mol^-1)
                mean(EAD); ... % "rest"
                mean(EAD); ... % refractive part
                ];
                Ea_glu = [EAG; mean(EAG)*ones(2,1)];
            params.EA_glu = sqrt(Ea_glu*Ea_glu'); % same idea as for K_glu0
            params.EA_glu = sign(Ea_glu(1))*params.EA_glu; % force positive
        % do the same for "rockifying":
            k_rock0 = [KR; mean(KR); 0];
            k_rdis0 = [KRD; mean(KRD); 1e-10];
            params.Ea_rdis = [EARD; mean(EARD)*ones(2,1)];
            params.Ea_rock = [EAR; mean(EAR)*ones(2,1)];
        % get time-dependent parameters, if wanted:
            params.K_glu0 = K_glu0;
                k_glu0_max = [KGmax; mean(KGmax)*ones(2,1)];
            params.K_glu0_max = sqrt(k_glu0_max*k_glu0_max');
                k_glu0_gainT = [KGgainT; mean(KGgainT); 1e10];
            params.K_glu0_gainT = sqrt(k_glu0_gainT*k_glu0_gainT');
            params.k_dis0 = k_dis0;
            params.k_dis0_min = [KDmin; mean(KDmin)*ones(2,1)];
            params.k_dis0_decT = [KDdecT; mean(KDdecT); 1e10];

        % Get initial molecule numbers:       
            Totc = nrComps + Rest + Reff; % Total number of compounds that will be simulated
            if Reff
                Nrefr0 = ReffDia^3*pi/6 * params.densReff*1e3 / MW(end) * 6.022e23; % (molec) ... #molecules in refractory core
            else
                Nrefr0 = 0;
            end
            if nrComps ~= size(chi,2)-1
                disp('params.chi not the right size!')
                stop_this;
                return
            end
            if Rest
                params.chi = [chi 1-sum(chi(:,2:end),2)];
            else % re-adjust chi so that sum is 1
                params.chi = [chi(1,1) chi(1,2:end)./sum(chi(1,2:end),2)];
                if size(chi,1)>1
                    for i=2:size(chi,1)
                        params.chi(i,:) = [chi(i,1) chi(i,2:end)./sum(chi(i,2:end),2)];
                    end
                end
            end
            Nnonrefr0 = (DP.^3*pi/6-ReffDia^3*pi/6) * params.dens*1e3 / (params.chi(1,2:end)*MW(1:nrComps+Rest))*sum(params.chi(1,2:end)) * 6.022e23; % (molec) ... #molecules that can evaporate
            Ni0 = params.chi(1,2:end)' .* Nnonrefr0; % #molecules of specific composition
            if Rest
                NR0 = Ni0(end); % #molecules in remaining bulk ("Rest")
                Ni0 = Ni0(1:end-1);
            else
                NR0 = 0;
            end
            Ntot0 = sum(Ni0) + NR0 + Nrefr0; % total #molecules in particles
            % for checking that that worked:
            %test_dp = ( ([sum(Ni0) NR0 Nrefr0]*(params.MW./[params.dens; params.dens; params.densReff])) * 6/pi/1e3/6.022e23 )^(1/3); % ... should be = params.dp
            % add another column to params.chi in case of Reff==1 (zeros because shouldn't evaporate anyway)
            if Reff
                params.chi = [params.chi zeros(size(params.chi,1),1)];
            end
            % adjusting parameters according to Rest and Reff
                if Rest && Reff
                    N = [Ni0;NR0;Nrefr0];
                elseif Rest
                    N = [Ni0;NR0];
                elseif Reff
                    N = [Ni0;Nrefr0];
                else
                    N = [Ni0];
                end
                if Reff && ~Rest
                    I = [1:Totc-1 Totc+1];
                else
                    I = 1:Totc;
                end
                params.MW = MW(I);
                params.dHvap = dHvap(I);
                params.ma = ma(I);
                k_glu0 = k_glu0(I);
                K_glu0 = K_glu0(I,I);
                params.K_glu0 = params.K_glu0(I,I);
                params.K_glu0_max = params.K_glu0_max(I,I);
                params.K_glu0_gainT = params.K_glu0_gainT(I,I);
                k_dis0 = k_dis0(I);
                params.k_dis0 = params.k_dis0(I);
                params.k_dis0_min = params.k_dis0_min(I);
                params.k_dis0_decT = params.k_dis0_decT(I);
                k_rock0 = k_rock0(I);
                k_rdis0 = k_rdis0(I);
                params.Ea_dis = params.Ea_dis(I);
                params.EA_glu = params.EA_glu(I,I);
                params.Ea_rdis = params.Ea_rdis(I);
                params.Ea_rock = params.Ea_rock(I);
                Csat0 = Csat0(I);
                Csat0_ug = Csat0_ug(I);
                params.Csat0 = Csat0;
            
            % calculate starting equilibrum "glued" amounts
            if Glu
                if GluOrder == 2
                    if Rock == 1 || Rock == 2
                        disp('Rock == 1 or 2 only works for GluOrder == 1 or 3 so far!')
                        disp('STOPPING.')
                        stop_this;
                        return
                    end
                    Ng0 = N .* (K_glu0 * N) ./ k_dis0;
                    % re-adjust to original total numbers for initial equlibrium estimate:
                    Ng0 = Ng0./(N+Ng0).*N;
                    % ... thought I had this figured out, but didn't always work. Will try again later maybe
                        disp('GluOrder == 2 is actually not working at all currently, sorry!')
                        disp('STOPPING.')
                        stop_this;
                        return
                elseif GluOrder == 1
                    if Rock==0 || Rock==3
                        Ng0 = diag(K_glu0) .* N ./ (k_dis0 + diag(K_glu0));
                        Nr0 = zeros(length(Ni0)+Rest+Reff,1);
                    elseif Rock==1
                        % sol = solve([Ng == kg/(kd+kr) * Ni, Nr == kr*kg/(krd*(kd+kr)) * Ni, Ni + Ng + Nr == Ntot],[Ni, Ng, Nr]);
                        Ng0 = (N .* diag(K_glu0) .* k_rdis0) ./ (k_dis0.*k_rdis0 + diag(K_glu0).*k_rock0 + diag(K_glu0).*k_rdis0 + k_rock0.*k_rdis0);
                        Nr0 = (N .* diag(K_glu0) .* k_rock0) ./ (k_dis0.*k_rdis0 + diag(K_glu0).*k_rock0 + diag(K_glu0).*k_rdis0 + k_rock0.*k_rdis0);
                    elseif Rock==2
        				% sol = solve([Ng == kg/kd * Ni, Nr == kr/krd * Ni, Ni + Ng + Nr == Ntot],[Ni, Ng, Nr]);
        				Ng0 = (N.* diag(K_glu0) .* k_rdis0) ./ (k_dis0.*k_rock0 +  k_dis0.*k_rdis0 + diag(K_glu0).*k_rdis0);
        				Nr0 = (N.* k_dis0 .* k_rock0) ./ (k_dis0.*k_rock0 +  k_dis0.*k_rdis0 + diag(K_glu0).*k_rdis0);
                    end
                elseif GluOrder == 3
                    if Rock == 0 || Rock==3
                        Ng0 = diag(K_glu0) .* N ./ (k_dis0 + diag(K_glu0));
                        Nr0 = zeros(length(Ni0)+Rest+Reff,1);
                    elseif Rock == 1
                        Ng0 = (N .* diag(K_glu0) .* k_rdis0) ./ (k_dis0.*k_rdis0 + diag(K_glu0).*k_rock0 + diag(K_glu0).*k_rdis0 + k_rock0.*k_rdis0);
                        Nr0 = (N .* diag(K_glu0) .* k_rock0) ./ (k_dis0.*k_rdis0 + diag(K_glu0).*k_rock0 + diag(K_glu0).*k_rdis0 + k_rock0.*k_rdis0);
                    elseif Rock == 2
                        Ng0 = (N.* diag(K_glu0) .* k_rdis0) ./ (k_dis0.*k_rock0 +  k_dis0.*k_rdis0 + diag(K_glu0).*k_rdis0);
        				Nr0 = (N.* k_dis0 .* k_rock0) ./ (k_dis0.*k_rock0 +  k_dis0.*k_rdis0 + diag(K_glu0).*k_rdis0);
                    end
                end
                if Rest && Reff
                    Ni0 = Ni0 - Ng0(1:end-2) - Nr0(1:end-2);
                    NR0 = NR0 - Ng0(end-1) - Nr0(end-1);
                    Nrefr0 = Nrefr0 - Ng0(end) - Nr0(end);
                    N = [Ni0;NR0;Nrefr0];
                elseif Rest
                    Ni0 = Ni0 - Ng0(1:end-1) - Nr0(1:end-1);
                    NR0 = NR0 - Ng0(end) - Nr0(end);
                    N = [Ni0;NR0];
                elseif Reff
                    Ni0 = Ni0 - Ng0(1:end-1) - Nr0(1:end-1);
                    Nrefr0 = Nrefr0 - Ng0(end) - Nr0(end);
                    N = [Ni0;Nrefr0];
                else
                    Ni0 = Ni0 - Ng0 - Nr0;
                    N = [Ni0];
                end
                % for 2nd-order stuff: employ ODE solver to get it right (because I'm too stupid otherwise):
                y0 = [Ng0' N'];
                if GluOrder == 1 || GluOrder == 3
                    Ng0 = y0(1:length(y0)/2);
                    if Rest && Reff
                        Ni0 = y0(length(y0)/2+1:end-2);
                        NR0 = y0(end-1);
                        Nrefr0 = y0(end);
                    elseif Rest
                        Ni0 = y0(length(y0)/2+1:end-1);
                        NR0 = y0(end);
                        Nrefr0 = 0;
                    elseif Reff
                        Ni0 = y0(length(y0)/2+1:end-1);
                        NR0 = 0;
                        Nrefr0 = y0(end);
                    else
                        Ni0 = y0(length(y0)/2+1:end);
                        NR0 = 0;
                        Nrefr0 = 0;
                    end
                elseif GluOrder == 2
                    disp('Hmm, GluOrder==2 ... shouldn''t this have stopped already?')
                    disp('Trying stopping again.')
                    stop_this;
                    return
                end
            else
                Ng0 = zeros(length(Ni0)+Rest+Reff,1);
                Nr0 = zeros(length(Ni0)+Rest+Reff,1);
            end
            if size(Ng0,2) > size(Ng0,1)
                Ng0=Ng0';
            end
            if size(Ni0,2) > size(Ni0,1)
                Ni0=Ni0';
            end
            if Rock==3
                    Itmp = length(params.RockFrac);
                    if length(Ni0)~=Itmp
                        disp('Need RockFrac for each i!');
                        disp('STOPPING.');
                        stop_this;
                    end
                Nr0(1:Itmp) = (N(1:Itmp)+Ng0(1:Itmp)).*params.RockFrac;
                    Ni0tmp = Ni0(1:Itmp);
                Ni0 = Ni0(1:Itmp) - Nr0(1:Itmp) .* Ni0(1:Itmp)./(Ni0(1:Itmp)+Ng0(1:Itmp));
                Ng0(1:Itmp) = Ng0(1:Itmp) - Nr0(1:Itmp) .* Ng0(1:Itmp)./(Ni0tmp+Ng0(1:Itmp));
                if params.Rest
                    Nr0(nrComps+1) = (NR0+Ng0(nrComps+1)).*params.RockFracRe;
                        NR0tmp = NR0;
                    NR0 = NR0 - Nr0(nrComps+1) .* NR0./(NR0+Ng0(nrComps+1));
                    Ng0(nrComps+1) = Ng0(nrComps+1) - Nr0(nrComps+1) .* Ng0(nrComps+1)./(NR0tmp+Ng0(nrComps+1));
                end
            end

        % Let's also get some wall effects:
            % paramaters for ODE solver: 
                if ( Reff && ~Rest )||( ~Reff && Rest )
                    params.C_w = CW ./mean(params.MW(1:end-1)) *6.022e23 *1e-12; % molec cm^-3
                elseif Reff && Rest
                    params.C_w = CW ./mean(params.MW(1:end-2)) *6.022e23 *1e-12; % molec cm^-3
                else
                    params.C_w = CW ./mean(params.MW) *6.022e23 *1e-12; % molec cm^-3
                end
                %params.k_woff = k_won ./ (params.C_w./Csat_p); % s^-1
                % initial conditions:
                kwoff0 = params.k_won ./ (params.C_w./Csat0); % unitless
                Nw0 = zeros(Totc,1);
            % now also 2nd wall (Walls==2), i.e. IMR wall:
                if ( Reff && ~Rest )||( ~Reff && Rest )
                    params.C_wI = CWI ./mean(params.MW(1:end-1)) *6.022e23 *1e-12;
                elseif Reff && Rest
                    params.C_wI = CWI ./mean(params.MW(1:end-2)) *6.022e23 *1e-12;
                else
                    params.C_wI = CWI ./mean(params.MW) *6.022e23 *1e-12;
                end
                kwoffI0 = params.k_wonI ./ (params.C_wI./Csat0);
                NwI0 = zeros(Totc,1);

        % ---
        % Run
        % ---
        
        for rr=1:(Soak+1+Wait)

            % Setting up desorption temperature as a function of time
            if Wait && rr==1 % isothermal wait period (no heating)
                params.T = struct(...
                    'equation','@(x) T0',...
                    'type','linear',...
                    'a',0);
                y0(1) = T0;
                t1 = 1:(60*wait_time);
            elseif rr==1+Wait % linear desorption T ramp
                params.T = struct(...
                    'equation','@(x) T0+rampRate*x',...
                    'type','linear',...
                    'a',DT);
                t1 = 1:durRamp;
                if Wait
                    t1 = 60*wait_time+t1;
                end
                y0(1) = T0;
            elseif rr==2+Wait % soak period at maximum desorption T
                params.T = struct(...
                    'equation','@(x) T0+(maxT-T0)*DT/0.14',...
                    'type','linear',...
                    'a',0);
                y0(1) = T0 + (maxT-T0)*DT/0.14;
                t1 = durRamp:(durRamp+durSoak);
                if t1(end) > params.TotDur
                    t1 = t1(1):floor(params.TotDur);
                end
                if Wait
                    t1 = 60*wait_time+t1;
                end
            elseif rr==3+Wait % This if for cool-down phase after soak. Currently not used.
                params.T = struct(...
                    'equation','@(x) ((maxT-T0)*exp((durRamp+durSoak)*b))*exp(-b*x)+c',...
                    'type','exponential decay',...
                    'b',0.005117,...
                    'c',298.15);
                    b = params.T.b; c = params.T.c;
                    tmp = eval(params.T.equation);
                y0(1) = tmp(durSoak+durRamp);
                t1 = (durSoak+durRamp):params.TotDur;
                if Wait
                    t1 = 60*wait_time+t1;
                    params.T.equation = ['@(x) ((maxT-T0)*exp((durRamp+durSoak+' num2str(60*wait_time) ')*b))*exp(-b*x)+c'];
                end
            end

            % Running the ODE solver:

            % forcing some steps (in addition to discontinuities in dT/dt) if time here is > 30 min (because long-term (> hours) behavior of ODE solver solution is weird)
                nSteps = ceil(length(t1)/3600/2);
                stepszMax = 2*3600;
            for r=1:nSteps
                t = [t1(1+(r-1)*stepszMax) min(t1(end),t1(1)+r*stepszMax-1)];
                
                % Setting up vector with initial conditions for ODE solver
                if r==1 && rr==1
                    y0 = [y0(1) ...
                        Csat0' ...
                        reshape(K_glu0',1,Totc*Totc) ...
                        k_dis0' ...
                        k_rock0' ...
                        k_rdis0' ...
                        Ng0' ...
                        Nr0' ...
                        Ni0' ...
                        ];
                    if Rest
                        y0 = [y0 NR0];
                    end
                    if Reff % if params.chi is adjusted above and used for getting Ni0 initially, Ni0 will already include Nrefr0
                        y0 = [y0 Nrefr0];
                    end
                    if Walls == 1
                        y0 = [y0 ...
                            y0(1)' ...
                            Csat0' ...
                            kwoff0' ...
                            Nw0' ...
                            ];
                    elseif Walls == 2
                        y0 = [y0 ...
                            y0(1)' ...
                            Csat0' ...
                            kwoff0' ...
                            Nw0' ...
                            y0(1)' ...
                            Csat0' ... 
                            kwoffI0' ...
                            NwI0' ...
                            ];
                    end
                else
                    y0 = dy.y(:,end)';
                end

                % FINALLY RUN THAT ODE SOLVER!
                    %options = odeset('NonNegative',logical(ones(size(y0))));
                    %options = odeset('AbsTol',[1e-6 1e-6*ones(1,length(y0)-1)]);
                    options = odeset('RelTol',1e-5);
                dy = ode15s(@thermogram_simul_ODE,t,y0,options,params);

                % collect (save) output for each step in this loop into dy_sv
                if r==1 && rr==1
                    dy_sv = dy;
                else
                    dy_sv.x = [dy_sv.x dy.x(:,2:end)];
                    dy_sv.y = [dy_sv.y dy.y(:,2:end)];
                end
            end

        end

        % collect output for each run of TR and SR loops in dy_sv2
        % ... and also some info on the correspondingly used desorption temperatures and particles sizes
        dy_sv2{TR,1} = T0+DT*(1:durRamp); % save separately the intended desorption T ramp
        dy_sv2{TR,2}{SR,2} = dy_sv;
            Tcontr = flipud(Ty');
            Tcontr = Ty(TrunNos);
        dy_sv2{TR,3} = Tcontr(TR);
        dy_sv2{TR,2}{SR,1} = DP;
        if SrunNos == 0
            dy_sv2{TR,2}{1,3} = 1;
        else
            dy_sv2{TR,2}{SR,3} = params.sizeDistFrac(SrunNos(SR),2);
        end

    end
end

%% Getting out rate of molecules entering CIMS
% (i.e. what should be proportional to observed count rate)

% ... need to use dy_sv2, in case SR or TR > 1

% get time base for end result
TotTimes = 0;
for t = 1:size(dy_sv2,1)
    for s = 1:size(dy_sv2{t,2},1)
        TotTimes = TotTimes + length(dy_sv2{t,2}{s,2}.x);
    end
end
tibase = nan(TotTimes,1);
rti=0;
for t = 1:size(dy_sv2,1)
    for s = 1:size(dy_sv2{t,2},1)
        tibase(rti+1:rti+length(dy_sv2{t,2}{s,2}.x)) = dy_sv2{t,2}{s,2}.x;
        rti=rti+length(dy_sv2{t,2}{s,2}.x);
    end
end
tibase = unique(sortrows(tibase));

% initiate main output file (and put there time and temperature)
simulation = nan(2+Totc,length(tibase)); % time, temperature, signal
simulation(1,:) = tibase;
    % for bugfixing:
    %disp(['x: ' num2str(size(dy_sv2{end,2}{1,2}.x)) ', y1: ' num2str(size(dy_sv2{end,2}{1,2}.y(1,:))) ', ti: ' num2str(size(tibase))]);
    % if something went seriously wrong (early on!), "stop"
    try
        simulation(2,:) = interp1(dy_sv2{end,2}{1,2}.x,dy_sv2{end,2}{1,2}.y(1,:),tibase); % assuming the last model run is actually for nominal desorption temperature!
    catch
        simulation = zeros(3,1);
        return
    end
% initiate other output files
simNi = nan(Totc,size(simulation,2));
simNg = simNi;
simNr = simNi;
simNw = simNi;
   
% get individual results and put that into output files
Start = length(params.MW)^2+6*length(params.MW)+2;
for i=Start:Start+Totc-1
    r = i-Start+1;
    simMat = cell(size(dy_sv2,1) * size(dy_sv2{1,2},1),2);
    simNiMat = cell(size(dy_sv2,1) * size(dy_sv2{1,2},1),1);
    simNgMat = simNiMat;
    simNrMat = simNiMat;
    simNwMat = simNiMat;
    rsi = 0;
    for t = 1:size(dy_sv2,1)
        for s = 1:size(dy_sv2{t,2},1)
            rsi=rsi+1;
            if Walls==1
                simMat{rsi,1} = [dy_sv2{t,2}{s,2}.x; dy_sv2{t,2}{s,2}.y(Start+2*length(params.MW)+r,:) .* dy_sv2{t,2}{s,2}.y(Start+3*length(params.MW)+r,:)]; % time; kwoff * Nw
            elseif Walls==2
                simMat{rsi,1} = [dy_sv2{t,2}{s,2}.x; dy_sv2{t,2}{s,2}.y(Start+1+5*length(params.MW)+r,:) .* dy_sv2{t,2}{s,2}.y(Start+1+6*length(params.MW)+r,:)]; % time; kwoff * Nw ("IMR")
            end
            simMat{rsi,2} = dy_sv2{t,2}{s,3} * dy_sv2{t,3}/sum(cell2mat(dy_sv2(:,3)));
                simNgMat{rsi,1} = dy_sv2{t,2}{s,2}.y(i-2*length(params.MW),:);
                simNrMat{rsi,1} = dy_sv2{t,2}{s,2}.y(i-length(params.MW),:);
                if Walls==1
                    simNwMat{rsi,1} = dy_sv2{t,2}{s,2}.y(Start+3*length(params.MW)+r,:);
                elseif Walls==2
                    simNwMat{rsi,1} = dy_sv2{t,2}{s,2}.y(Start+1+6*length(params.MW)+r,:);
                end
            simNiMat{rsi,1} = dy_sv2{t,2}{s,2}.y(i,:);
        end
    end
    simulation(2+r,:) = zeros(1,size(simulation,2));
    simNi(r,:) = zeros(1,size(simulation,2));
    simNg(r,:) = simNi(r,:);
    simNr(r,:) = simNi(r,:);
    simNw(r,:) = simNi(r,:);
    for j=1:size(simMat,1)
        simulation(2+r,:) = simulation(2+r,:) + interp1(simMat{j,1}(1,:),simMat{j,1}(2,:),tibase')*simMat{j,2};
        simNi(r,:) = simNi(r,:) + interp1(simMat{j,1}(1,:),simNiMat{j,1},tibase')*simMat{j,2};
        simNg(r,:) = simNg(r,:) + interp1(simMat{j,1}(1,:),simNgMat{j,1},tibase')*simMat{j,2};
        simNr(r,:) = simNr(r,:) + interp1(simMat{j,1}(1,:),simNrMat{j,1},tibase')*simMat{j,2};
        simNw(r,:) = simNw(r,:) + interp1(simMat{j,1}(1,:),simNwMat{j,1},tibase')*simMat{j,2};
    end
end

end

function [simulation,simNi,simNg,simNr,simNw,dy_sv2] = stop_this
    simulation = NaN;
    simNi = NaN;
    simNg = NaN;
    simNr = NaN;
    simNw = NaN;
    dy_sv2 = NaN;
end

    