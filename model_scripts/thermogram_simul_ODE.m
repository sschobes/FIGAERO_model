function dy = thermogram_simul_ODE(t,y,params)

% INPUT:
% t: time range (s) for which to solve the set of ODE: dy/dt = f(t,y).
% y: initial conditions for y(t), i.e. y(t=0).
% params: structure containing model parameters
% 
% OUTPUT:
% dy: structure containing the solution:
%   - x: time t (s)
%   - y: solution y(t) for each y

%% Get, derive, set and calculate some parameters

% Extract model parameters again
Wait = params.Wait;
Walls = params.Walls;
Glu = params.Glu;
GluOrder = params.GluOrder;
GluTimDep = params.GluTimDep;
Rock = params.Rock;
Raoult = params.Raoult;
R = params.R; % gas constant
dp = params.dp; % particle diameter (m)
dens = params.dens; % particle density (kg m^-3)
densReff = params.densReff; % density of refractory part (kg m^-3)
MW = params.MW; % g mol^-1 
dHvap = params.dHvap; % kJ mol^-1
Csat0 = params.Csat0; % molec cm^-3
R = params.R; % J mol^-1 K^-1
pressure = params.pressure; % mbar
ma = params.ma; % mass accommodation coefficient
Ea_dis = params.Ea_dis; % kJ mol^-1
EA_glu = params.EA_glu; % kJ mol^-1
Ea_rdis = params.Ea_rdis; % kJ mol^-1
Ea_rock = params.Ea_rock; % kJ mol^-1
K_glu0 = params.K_glu0; % s^-1
K_glu0_max = params.K_glu0_max; % s^-1
K_glu0_gainT = params.K_glu0_gainT; % s
k_dis0 = params.k_dis0; % s^-1
k_dis0_min = params.k_dis0_min; % s^-1
k_dis0_decT = params.k_dis0_decT; % % s
k_won = params.k_won; % s^-1
C_w = params.C_w; % molec cm^-3
Tfrac_after = params.Tfrac_after;
k_wonI = params.k_wonI;
C_wI = params.C_wI;
Tfrac_afterI = params.Tfrac_afterI;
T0 = params.T0; % reference temperature for thermodynamics (K)
maxT = params.maxT; % maximum desorption temperature (at end of ramp and during soak) (K)

% Derive, or adjust units:
%Csat0 = 10.^((131-dHvap)/11); % ug m^-3 ("simple" Epstein et al. (EST, 2009))
Csat0 = Csat0 ./MW *6.022e23 *1e-12; % molec cm^-3
dHvap = dHvap*1e6; % g m^2 s^-2 mol^-1
R = R*1e3; % g m^2 s-^2 mol^-1 K^-1
kB = 1.38064852e-23; % m^2 kg s^-2 K^-1 (= J/K)
    NA = 6.022140857e23; % molec mol^?1
mm = MW / NA *1e-3; % kg molec^-1
Ea_dis = Ea_dis*1e6; % g m^2 s^-2 mol^-1
EA_glu = EA_glu*1e6; % g m^2 s^-2 mol^-1
Ea_rdis = Ea_rdis*1e6; % g m^2 s^-2 mol^-1
Ea_rock = Ea_rock*1e6; % g m^2 s^-2 mol^-1

%% Define ODE set (and calculate T- or time-dependent parameters)

% initialize output variable
dy=zeros(length(y),1);

% Temperature (K):
switch params.T.type
    case 'linear'
        dy(1) = params.T.a;
    case 'exponential decay'
        if Wait
            dy(1) = -(maxT-T0)*params.T.b*exp(params.durSoak*params.T.b).*exp(-params.T.b*t);
        else
            dy(1) = -(maxT-T0)*params.T.b*exp((params.durSoak+60*params.wait_time)*params.T.b).*exp(-params.T.b*t);
        end
    otherwise
        dy(1) = 0;
end
r=1;

% Csat
for i=1:length(MW)
    r=r+1;
    dy(r) = y(r) * dy(1) * dHvap(i)/R./y(1).^2;
end

% k_glu
    I_glu = nan(length(MW));
for i=1:length(MW)
    for j=1:length(MW)
        r=r+1;
        I_glu(i,j) = r;
        if Glu % kg = kg0 * exp(-E/R*(1/T-1/T0))
            dy(r) = y(r) * dy(1) * EA_glu(i,j)/R./y(1).^2;
            if GluTimDep % kg = (kg0 + kg0maxInc*(1-exp(-t/gain)) * exp(-E/R*(1/T-1/T0))
                K_glu0_maxInc = K_glu0_max(i,j) - K_glu0(i,j);
                %dy(r) = dy(r) + dy(r)/K_glu0(i,j)*K_glu0_max(i,j) - dy(r)/K_glu0(i,j)*K_glu0_max(i,j)*exp(-t/K_glu0_gainT(i,j)) + K_glu0_max(i,j)/K_glu0_gainT(i,j)*exp(EA_glu(i,j)/R*(1/y(1)-1/T0))*exp(-t/K_glu0_gainT(i,j));
                % = 
                dy(r) = dy(r) + dy(r)/K_glu0(i,j)*K_glu0_maxInc - exp(-t/K_glu0_gainT(i,j))*K_glu0_maxInc * (dy(r)/K_glu0(i,j) - 1/K_glu0_gainT(i,j)*exp(EA_glu(i,j)/R*(1/y(1)-1/T0)) );
            end
        else
            dy(r) = 0;
        end
    end
end

% k_dis
    I_dis = nan(1,length(MW));
for i=1:length(MW)
    r=r+1;
    I_dis(i) = r;
    if Glu % kd = kd0 * exp(-E/R*(1/T-1/T0))
        dy(r) = y(r) * dy(1) * Ea_dis(i)/R./y(1).^2;
        if GluTimDep % kg = (kd0 - kd0maxDec*(1-exp(-t/decay)) * exp(-E/R*(1/T-1/T0))
            k_dis0_maxDec = k_dis0(i) - k_dis0_min(i);
            dy(r) = dy(r) - dy(r)/k_dis0(i)*k_dis0_maxDec + exp(-t/k_dis0_decT(i))*k_dis0_maxDec * (dy(r)/k_dis0(i) - 1/k_dis0_decT(i)*exp(Ea_dis(i)/R*(1/y(1)-1/T0)) );
        end
    else
        dy(r) = 0;
    end
end

% k_rock
    I_rock = nan(length(MW));
for i=1:length(MW)
    r=r+1;
    I_rock(i) = r;
    if Rock
        dy(r) = y(r) * dy(1) * Ea_rock(i)/R./y(1).^2;
    else
        dy(r) = 0;
    end
end

% k_rdis
    I_rdis = nan(1,length(MW));
for i=1:length(MW)
    r=r+1;
    I_rdis(i) = r;
    if Rock
        dy(r) = y(r) * dy(1) * Ea_rdis(i)/R./y(1).^2;
    else
        dy(r) = 0;
    end
end

% Nig
    I_Ng = nan(1,length(MW));
for i=1:length(MW)
    r=r+1;
    if i==1; r0 = r; end
    I_Ng(i) = r;
    if Glu
        if GluOrder == 2
            % 2nd order:
            glu_term = y(r+2*length(MW))*(y(I_glu(i,:))'*y(r0+2*length(MW):r0+3*length(MW)-1));
        elseif GluOrder == 1
            % pseudo-1st order:
            glu_term = y(r+2*length(MW))*y(I_glu(i,i));
        elseif GluOrder == 3
            % 2nd order, but simplified (requires OArem)
            if Raoult==1
                dp0 = dp;
                if params.Reff
                    dp1 = ( ((y([r0:(r0+length(MW)-1)]+2*length(MW))+y([r0:(r0+length(MW)-1)])+y([r0:(r0+length(MW)-1)]+length(MW)))'*(MW./[dens*ones(length(MW)-1,1); densReff])) * 6/pi/1e3/6.022e23 )^(1/3);
                    OArem = (dp1^3 - params.ReffDia^3) / (dp0^3 - params.ReffDia^3);
                else
                    dp1 = ( ((y([r0:(r0+length(MW)-1)]+2*length(MW))+y([r0:(r0+length(MW)-1)])+y([r0:(r0+length(MW)-1)]+length(MW)))'*(MW./(dens*ones(length(MW),1)))) * 6/pi/1e3/6.022e23 )^(1/3);
                    OArem = dp1^3 / dp0^3;
                end
            else % i.e. Raoult==2
                OArem = interp1(params.OArem(:,1),params.OArem(:,2),t);
            end
            glu_term = y(r+2*length(MW))*y(I_glu(i,i))*OArem;
        end
        dis_term = - y(I_dis(i))*y(r)';
        if Rock == 1
            rock_term = - y(I_rock(i))*y(r)';
            dy(r) = glu_term + dis_term + rock_term;
        else
            dy(r) = glu_term + dis_term;
        end
    else
        dy(r) = 0;
    end
end

% Nr
    I_Nr =nan(1,length(MW));
for i=1:length(MW)
    r=r+1;
    if i==1; r0 = r; end
    I_Nr(i) = r;
    if Rock == 1
        rock_term = y(I_rock(i))*y(I_Ng(i));
        rdis_term = - y(I_rdis(i))*y(r)';
        dy(r) = rock_term + rdis_term;
    elseif Rock == 2
        if GluOrder == 1
            rock_term = y(I_rock(i))*y(r+length(MW));
        elseif GluOrder == 3
            % 2nd order, but simplified (requires OArem)
            if Raoult==1
                dp0 = dp;
                if params.Reff
                    dp1 = ( ((y([r0:(r0+length(MW)-1)]+length(MW))+y(I_Ng)+y([r0:(r0+length(MW)-1)]))'*(MW./[dens*ones(length(MW)-1,1); densReff])) * 6/pi/1e3/6.022e23 )^(1/3);
                    OArem = (dp1^3 - params.ReffDia^3) / (dp0^3 - params.ReffDia^3);
                else
                    dp1 = ( ((y([r0:(r0+length(MW)-1)]+length(MW))+y(I_Ng)+y([r0:(r0+length(MW)-1)]))'*(MW./dens*ones(length(MW),1))) * 6/pi/1e3/6.022e23 )^(1/3);
                    OArem = dp1^3 / dp0^3;
                end
            else % i.e. Raoult==2
                OArem = interp1(params.OArem(:,1),params.OArem(:,2),t);
            end
            rock_term = y(I_rock(i))*y(r+length(MW))*OArem;
        end
    elseif Rock == 3
        rock_term = 0;
    end
    if Rock >= 1
    	rdis_term = - y(I_rdis(i))*y(r)';
    	dy(r) = rock_term + rdis_term;
    else
        dy(r) = 0;
    end
end

% Ni
    I_Ni = nan(1,length(MW));
    I_Ni = r+1:r+length(MW);
for i=1:length(MW)
        
    r=r+1;
    if i==1; r0 = r; end

    % time-/temperature-dependent parameters in Hertz-Knudsen equation:
    D = sqrt(18./MW) .* diffusivity_H2O(y(1),pressure); % m^2/s (~1.5e-5)
    v_rms = sqrt(8*R*y(1)/pi./MW); % m/s
    lambda = 3*D./v_rms;
    if Raoult==1
        if params.Reff
            dp = ( ((y(I_Ni)+y(I_Ng)+y(I_Nr))'*(MW./[dens*ones(length(MW)-1,1); densReff])) * 6/pi/1e3/6.022e23 )^(1/3);
        else
            dp = ( ((y(I_Ni)+y(I_Ng)+y(I_Nr))'*(MW./dens*ones(length(MW),1))) * 6/pi/1e3/6.022e23 )^(1/3);
        end
    else % i.e. Raoult==2
        if GluOrder~=3 || Glu~=1 % otherwise would have OArem already
            OArem = interp1(params.OArem(:,1),params.OArem(:,2),t);
        end
        if ~params.Reff
            dp = dp * OArem^(1/3);
        else
            dp = (params.ReffDia^3 + OArem*(dp^3-params.ReffDia^3))^(1/3);
        end
    end
    Kn = 2*lambda/dp; % Knudsen number
    % correction for non-free molecular regime:
    Gamma = ma .* Kn.*(1+Kn)./(Kn.^2+Kn+0.283*Kn*1+0.75*1); % Trump et al., AST, 2016 (10.1080/02786826.2016.1232858)
    
    % Raoult factor
    if Raoult==2
        if length(MW)~=size(params.chi,2)-1
            disp('params.chi not the right size!')
            return
        end
        chi = interp1(params.chi(:,1),params.chi(:,2),t);
        % including max-statement to try to avoid getting negative molecule numbers
        chi = chi*max([y(r)-2 0])/(y(r)+y(I_Ng(i))+y(I_Nr(i)));
    elseif Raoult==1
        %chi = y(r)*MW(i) / (y(I_Ni(1:end-params.Reff))'*MW(1:end-params.Reff) + y(I_Ng(1:end-params.Reff))'*MW(1:end-params.Reff) + y(I_Nr(1:end-params.Reff))'*MW(1:end-params.Reff));
        % including max-statement and also a "cushion" of 100 molecules to try to avoid getting negative molecule numbers
        chi = max([y(r)*MW(i)-2 0])/(100+ (y(I_Ni(1:end-params.Reff))'*MW(1:end-params.Reff) + y(I_Ng(1:end-params.Reff))'*MW(1:end-params.Reff) + y(I_Nr(1:end-params.Reff))'*MW(1:end-params.Reff)) );
    end
    % next one kosher? (once more trying to avoid negative molecule numbers!)
    if chi<0; chi=0; end
    
    % get gluing terms
    if Glu
        if GluOrder == 2
            glu_term = - y(r)*(y(I_glu(i,:))'*y(r0:r0+length(MW)-1)) + y(I_dis(i))*y(r-2*length(MW));
        elseif GluOrder == 1
            glu_term = - y(r)*y(I_glu(i,i)) + y(I_dis(i))*y(r-2*length(MW));
        elseif GluOrder == 3
            glu_term = - y(r)*y(I_glu(i,i))*OArem + y(I_dis(i))*y(r-2*length(MW));
        end
    else
        glu_term = 0;
    end

    % get rocking terms
    if Rock == 1
        rock_term = y(I_rdis(i))*y(r-length(MW));
    elseif Rock == 2
        if GluOrder == 1
            rock_term = y(I_rdis(i))*y(r-length(MW)) - y(r)*y(I_rock(i));
        elseif GluOrder == 3
            rock_term = y(I_rdis(i))*y(r-length(MW)) - y(r)*y(I_rock(i))*OArem;
        end
    elseif Rock == 3
        rock_term = y(I_rdis(i))*y(r-length(MW));
    else
        rock_term = 0;
    end
    
    % get evaporation based on Hertz-Knudsen equation (except surface area and chi)
    evap_base = - 1/sqrt(2*pi*kB*mm(i)*y(1)) * Gamma(i) * y(1+i)/6.022e23*1e3*(pressure*100)*(22.4*(y(1)/273.15));
   
        
    % include surface area and gluing and rocking
        SA = pi*dp^2; % surface area (m^2) ... i.e. spherical particle
    dy(r) = evap_base*chi*SA + glu_term + rock_term; % molec s^-1 (per collected particle)
    
end

if Walls == 1 || Walls == 2

    % T after filter (i.e. at walls)
    r=r+1;
    I_Tw = r;
    dy(r) = dy(1) * Tfrac_after;

    % Csat for walls
        I_Csw = nan(1,length(MW));
    for i=1:length(MW)
        r=r+1;
        I_Csw(i) = r;
        dy(r) = y(r) * dy(I_Tw) * dHvap(i)/R./y(I_Tw).^2;
    end

    % k_woff
        I_kw = nan(1,length(MW));
    for i=1:length(MW)
        r=r+1;
        I_kw(i) = r;
        dy(r) = k_won / C_w * dy(I_Csw(i));
    end

    % Nwall
        I_Nw = nan(1,length(MW));
    for i=1:length(MW)
        r=r+1;
        if i==1; r0 = r; end
        I_Nw(i) = r;
        dy(r) = k_won * (-dy(I_Ni(i))-dy(I_Ng(i))+y(I_rdis(i))*y(I_Nr(i))) - y(I_kw(i)) * y(r);
    end
    
    if Walls == 2
        
        % T after filter (i.e. at 2nd wall set)
        r=r+1;
        I_TwI = r;
        dy(r) = dy(1) * Tfrac_afterI;

        % Csat for 2nd wall set
            I_CswI = nan(1,length(MW));
        for i=1:length(MW)
            r=r+1;
            I_CswI(i) = r;
            dy(r) = y(r) * dy(I_TwI) * dHvap(i)/R./y(I_TwI).^2;
        end

        % k_woff
            I_kwI = nan(1,length(MW));
        for i=1:length(MW)
            r=r+1;
            I_kwI(i) = r;
            dy(r) = k_wonI / C_wI * dy(I_CswI(i));
        end
        
        % NwallI
            I_NwI = nan(1,length(MW));
        for i=1:length(MW)
            r=r+1;
            if i==1; r0 = r; end
            I_NwI(i) = r;
            dy(r) = k_wonI * y(I_kw(i)) * y(I_Nw(i)) - y(I_kwI(i)) * y(r);
        end
    
    end
    
end

if r ~= length(y)
    disp('Problem!')
end

