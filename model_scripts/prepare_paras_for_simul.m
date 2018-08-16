
%%
% These scripts prepare some parameters for the eventual model runs, get
% additional ones, and then puts them into the structure "params" (if not
% included yet), which will be passed on down to the function running the
% simulation.

%%

% Limit model length to larger of T ramp+soak or length of measdat (to hopefully avoid problems if Raoult==2)
    tmp = interp1(1:timres:timres*length(measdat(:,1)),measdat(:,1),1:timres*length(measdat(:,1))); % length/duration of measdat in s
    T = 273.15 + [(params.T0-273.15):0.14:200 200*ones(1,4190-length((params.T0-273.15):0.14:200)) tmp(4191:end)]; % -> length(T) = larger of ramp+soak or tmp
params.TotDur = length(T);
params.TotDur = params.TotDur - timres*5; % might be a problem those last steps, at least if Raoult==2

% Decide if actually doing a wait period
if ~exist('wait_time') || ~wait_time>0
    WaitInModel = 0;
end

% Revert to all ones if UseNdist==0
if ~UseNdist
    Ndist = ones(nrComps,1);
end

% get either:
% - time series (in model time) for the Raoult term (chi) and OA fraction remaining (phi), based on experimental data, if needed (i.e. if Raoult==2)
% - or only initial fractions for each compound (i.e. if Raoult==1)
if Raoult==2
    % get model time, and relation to experimental time (make time vector as per model parameters and use the point of passing 50 deg C as reference)
    if WaitInModel
        T_vs_t_theo = [[(1:(60*wait_time))'; 60*wait_time+(1:params.durRamp)'; 60*wait_time+((params.durRamp+1):(params.durSoak+params.durRamp))'] [params.T0*ones(60*wait_time,1); params.T0+params.rampRate*(0:(params.durRamp-1))'; (273.15+params.maxT)*ones((params.durSoak+params.durRamp)-params.durRamp,1)]];
    else
        T_vs_t_theo = [[(1:params.durRamp)'; ((params.durRamp+1):(params.durSoak+params.durRamp))'] [params.T0+params.rampRate*(0:(params.durRamp-1))'; (273.15+params.maxT)*ones((params.durSoak+params.durRamp)-params.durRamp,1)]];
    end
    Tm50_theo = find(T_vs_t_theo(:,2)>273.15+50,1,'first');
    tm50_theoExact = interp1(T_vs_t_theo(Tm50_theo-1:Tm50_theo,2),T_vs_t_theo(Tm50_theo-1:Tm50_theo,1),273.15+50);
    T_vs_t_data = [(0:timres:timres*(size(measdat,1)-1))' measdat(:,2)+273.15];
    Tm50_data = find(T_vs_t_data(:,2)>273.15+50,1,'first');
    tm50_dataExact = interp1(T_vs_t_data(Tm50_data-1:Tm50_data,2),T_vs_t_data(Tm50_data-1:Tm50_data,1),273.15+50);
    shift = tm50_theoExact - tm50_dataExact; % s
    % get Raoult term (chi)
    for i=1:size(chi_exp,2)
            params.chi = [T_vs_t_theo(:,1) interp1(T_vs_t_data(:,1)+shift,chi_exp(:,i),T_vs_t_theo(:,1),'linear')];
            tmp = interp1(T_vs_t_data(:,1)+shift,chi_exp(:,i),T_vs_t_theo(:,1),'nearest','extrap');
        params.chi(isnan(params.chi(:,1+i)),2) = tmp(isnan(params.chi(:,1+i)));
        if nrComps>1 && size(params.chi,2)==2 % i.e. chi_exp is for only 1 compount, but nrComps>1
            params.chi = [params.chi(:,1) repmat(params.chi(:,2),1,nrComps)];
            for i=1:nrComps
                params.chi(:,i+1) = params.chi(:,i+1)*Ndist(i)/sum(Ndist);
            end
        end
    end
    % get OA remaining (phi)
    params.OArem = params.chi(:,1);
    params.OArem(:,2) = interp1(T_vs_t_data(:,1)+shift,remA_exp,T_vs_t_theo(:,1),'linear');
    tmp = interp1(T_vs_t_data(:,1)+shift,remA_exp,T_vs_t_theo(:,1),'nearest','extrap');
    params.OArem(isnan(params.OArem(:,2)),2) = tmp(isnan(params.OArem(:,2)));
elseif Raoult==1 % need chi_exp only for initial conditions
    params.chi = [1 chi_exp(1,:)];
    if nrComps>1 && size(params.chi,2)==2 % i.e. chi_exp is for only 1 compount, but nrComps>1
        params.chi = [params.chi(:,1) repmat(params.chi(:,2),1,nrComps)];
        for i=1:nrComps
            params.chi(:,i+1) = params.chi(:,i+1)*Ndist(i)/sum(Ndist);
        end
    end
end

% put parameters into params for passing down to other functions:
params.Wait = WaitInModel;
params.wait_time = wait_time;
params.Walls = Walls;
params.Glu = Glu;
params.GluOrder = GluOrder;
params.GluTimDep = GluTimDep;
params.Rock = Rock;
params.RockFrac = RockFrac;
params.RockFracRe = RockFracRe;
params.Raoult = Raoult;
params.Rest = Rest;
params.Reff = Reff;
params.c = c;
params.nrComps = nrComps;
params.SrunNos = SrunNos;
params.MW = [...
    MOLW; ...
    MWRe; ... % Rest
    230; ... % Reff (just something, shouldn't matter)
    ];
params.ma = [...
    MA; ...
    MA; ... % Rest (i.e. assuming single evaporation coefficient for all compounds)
    1; ... % Reff (just something, shouldn't matter)
    ];
params.dHvap = [...
    DH; ...
    DHRe; ... % Rest
    0; ... % Reff (just something, shouldn't matter)
    ];
params.CSTARRe = CSTARRe;
