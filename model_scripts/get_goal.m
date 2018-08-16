function [goal,shift] = get_goal(measdat,timres,params)

% Get the variables "goal" and "shift".
% 
% INPUT:
% measdat: processed experimental data in 3+ columns: time (s), T (C), 
%       composition-specific signals
% timres: time-resolution of measdat (s)
% params: structure containing model parameters
% 
% OUTPUT:
% goal: experiental data, synchronized to the expected model output.
% shift: shift (s) needed for that sychronization.

% break out parameters from params as needed
WaitInModel = params.Wait;
wait_time = params.wait_time;

% make time vector as per model parameters
if WaitInModel
    T_vs_t_theo = [[(1:(60*wait_time))'; 60*wait_time+(1:params.durRamp)'; 60*wait_time+((params.durRamp+1):(params.durSoak+params.durRamp))'] [params.T0*ones(60*wait_time,1); params.T0+params.rampRate*(0:(params.durRamp-1))'; (273.15+params.maxT)*ones((params.durSoak+params.durRamp)-params.durRamp,1)]];
else
    T_vs_t_theo = [[(1:params.durRamp)'; ((params.durRamp+1):(params.durSoak+params.durRamp))'] [params.T0+params.rampRate*(0:(params.durRamp-1))'; (273.15+params.maxT)*ones((params.durSoak+params.durRamp)-params.durRamp,1)]];
end
% use the point of passing 50 deg C as reference.
Tm50_theo = find(T_vs_t_theo(:,2)>273.15+50,1,'first');
tm50_theoExact = interp1(T_vs_t_theo(Tm50_theo-1:Tm50_theo,2),T_vs_t_theo(Tm50_theo-1:Tm50_theo,1),273.15+50);
if size(measdat,2)>=3 && diff(measdat(1:2,1))==timres % concluding that first column in measdat is time series (and in s),
    T_vs_t_data = measdat(:,1:2);
    T_vs_t_data(:,2) = T_vs_t_data(:,2) + 273.15;
elseif isnan(measdat(1,1)) && size(measdat,2)>1 && ~isnan(measdat(1,2)) % assume time is not available, but that measdat is otherwise OK (i.e. containing data
    T_vs_t_data = [(0:timres:timres*(size(measdat,1)-1))' measdat(:,2)+273.15];
    Tm50_data = find(T_vs_t_data(:,2)>273.15+50,1,'first');
    tm50_dataExact = interp1(T_vs_t_data(Tm50_data-1:Tm50_data,2),T_vs_t_data(Tm50_data-1:Tm50_data,1),273.15+50);
    shift = tm50_theoExact - tm50_dataExact; % s
    [xu,XU] = unique(T_vs_t_data(:,1)+shift);
    goal = [T_vs_t_theo(:,1) interp1(xu,measdat(XU,2),T_vs_t_theo(:,1),'linear')+273.15 interp1(xu,measdat(XU,3),T_vs_t_theo(:,1),'linear')];
    tmp = interp1(xu,measdat(XU,2)+273.15,T_vs_t_theo(:,1),'nearest','extrap');
    goal(isnan(goal(:,2)),2) = tmp(isnan(goal(:,2)));
    tmp = interp1(xu,measdat(XU,3),T_vs_t_theo(:,1),'nearest','extrap');
    goal(isnan(goal(:,3)),3) = tmp(isnan(goal(:,3)));
else
    goal = nan(1,3); shift = NaN;
end

end