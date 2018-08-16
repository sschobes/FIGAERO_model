%%
% Scripts for plots of number of molecules vs. time

%%

if ShowInclRest && PlusRock3
    disp('Proper treatment of Rock required if Rest is used.');
    disp('... Changing to NoPlotNr=0 and PlusRock3=0');
    NoPlotNr=0; PlusRock3=0;
end

hold on

NrExists=0; try; if size(simNr,1)==size(simNi,1)  && size(simNr,2)==size(simNi,2); NrExists=1; end; end
if NrExists && ~NoPlotNr
    simTot = simNg+simNi+simNr;
    if params.Reff
        simTot = simTot-simNi(end,:)-simNg(end,:)-simNr(end,:);
    end
    if params.Rest
        simR = simNi(end-params.Reff,:)+simNg(end-params.Reff,:)+simNr(end-params.Reff,:);
        simTot = simTot-simR;
    end
    if ~NoPlotReff && ~NoPlotRest
        h3 = plot(simulation(1,:)/(1+(PlotMinutes*59)),simNr,'linewidth',2,'color',[0 0 0]);
    elseif ~NoPlotRest % && NoPlotReff
        h3 = plot(simulation(1,:)/(1+(PlotMinutes*59)),simNr(1:size(simNi,1)-params.Reff,:),'linewidth',2,'color',[0 0 0]);
    elseif ~NoPlotReff % && NoPlotRest
        if params.Reff && params.Rest
            h3 = plot(simulation(1,:)/(1+(PlotMinutes*59)),simNr([1:size(simNi,1)-2 size(simNi,1)],:),'linewidth',2,'color',[0 0 0]);
        elseif params.Reff % && ~params.Rest
            h3 = plot(simulation(1,:)/(1+(PlotMinutes*59)),simNr,'linewidth',2,'color',[0 0 0]);
        else % i.e. ~params.Reff
            h3 = plot(simulation(1,:)/(1+(PlotMinutes*59)),simNr(1:size(simNi,1)-params.Rest,:),'linewidth',2,'color',[0 0 0]);
        end
    else % i.e. NoPlot either
        h3 = plot(simulation(1,:)/(1+(PlotMinutes*59)),simNr(1:size(simNi,1)-params.Reff-params.Rest,:),'linewidth',2,'color',[0 0 0]);
    end
end
if PlusRock3
    simTot = sum(simNg+simNi+simNr,1);
    if params.Reff
        simTot = simTot-simNi(end,:)-simNg(end,:)-simNr(end,:);
    end
    if params.Rest
        simR = simNi(end-params.Reff,:)+simNg(end-params.Reff,:)+simNr(end-params.Reff,:);
        simTot = simTot-simR;
    end
    rock3s = rock3*simTot(1)*ones(size(simulation(1,:)))/(1-rock3); rock3s(simulation(1,:)>wait_time*60) = NaN;
end
if WaitOnly
    simNi(:,simulation(1,:)>wait_time*60) = NaN;
    simNg(:,simulation(1,:)>wait_time*60) = NaN;
    simNr(:,simulation(1,:)>wait_time*60) = NaN;
    simNw(:,simulation(1,:)>wait_time*60) = NaN;
end
if NrExists && ~PlusRock3
    simTot = sum(simNg+simNi+simNr,1);
    simR = simNi(end-params.Reff,:)+simNg(end-params.Reff,:)+simNr(end-params.Reff,:);
elseif ~NrExists && ~PlusRock3
    simTot = sum(simNg+simNi,1);
    simR = simNi(end-params.Reff,:)+simNg(end-params.Reff,:);
elseif NrExists && PlusRock3
    simTot = sum(simNg+simNi+simNr,1)+rock3s;
    simR = simNi(end-params.Reff,:)+simNg(end-params.Reff,:)+simNr(end-params.Reff,:);
else
    simTot = sum(simNg+simNi,1)+rock3s;
    simR = simNi(end-params.Reff,:)+simNg(end-params.Reff,:);
end
if ~params.Rest
    simR = zeros(size(simR));
end
if params.Reff
    simTot = simTot-simNi(end,:)-simNg(end,:)-simNr(end,:);
end
if params.Rest
    simTot = simTot-simR;
end
if ShowInclRest
    if ~NoPlotNr && NrExists
        h73 = plot(simulation(1,:)/(1+(PlotMinutes*59)),sum(simNr(1:size(simNr,1)-params.Reff,:),1),'linewidth',2,'color',[0.7 0.7 0.7]);
    end
    if ~NoPlotNg
        h72 = plot(simulation(1,:)/(1+(PlotMinutes*59)),sum(simNg(1:size(simNg,1)-params.Reff,:),1),'linewidth',2,'color',[0.7 0.5 1]);
    end
    h71 = plot(simulation(1,:)/(1+(PlotMinutes*59)),sum(simNi(1:size(simNi,1)-params.Reff,:),1),'linewidth',2,'color',[0.5 0.7 0.5]);
    if Extended>0
        if Extended==1
            h74 = plot(simulation(1,:)/(1+(PlotMinutes*59)),simTot+simR,'linewidth',1,'color',[1 0.4 0.4]);
        end
                apparentRemInclRest = simTot+simR+sum(simNw(1:size(simNw,1)-params.Reff,:),1);
                apparentRemInclRest(simulation(1,:)>60*wait_time) = NaN;
        if Extended==1
            h75 = plot(simulation(1,:)/(1+(PlotMinutes*59)),sum(simNw(1:size(simNi,1)-params.Reff,:),1),'--','linewidth',1,'color',[0.7 0.7 0.7]);  
        end
        h76 = plot(simulation(1,:)/(1+(PlotMinutes*59)),apparentRemInclRest,'linewidth',2,'color',[1 0.4 0.4]);
    end
end
if ~NoPlotNg
    if ~NoPlotReff && ~NoPlotRest
        h1 = plot(simulation(1,:)/(1+(PlotMinutes*59)),simNg,'linewidth',2,'color',[0 0 1]);
    elseif ~NoPlotRest % && NoPlotReff
        h1 = plot(simulation(1,:)/(1+(PlotMinutes*59)),simNg(1:size(simNi,1)-params.Reff,:),'linewidth',2,'color',[0 0 1]);
    elseif ~NoPlotReff % && NoPlotRest
        if params.Reff && params.Rest
            h1 = plot(simulation(1,:)/(1+(PlotMinutes*59)),simNg([1:size(simNi,1)-2 size(simNi,1)],:),'linewidth',2,'color',[0 0 1]);
        elseif params.Reff % && ~params.Rest
            h1 = plot(simulation(1,:)/(1+(PlotMinutes*59)),simNg,'linewidth',2,'color',[0 0 1]);
        else % i.e. ~params.Reff
            h1 = plot(simulation(1,:)/(1+(PlotMinutes*59)),simNg(1:size(simNi,1)-params.Rest,:),'linewidth',2,'color',[0 0 1]);
        end
    else % i.e. NoPlot either
        h1 = plot(simulation(1,:)/(1+(PlotMinutes*59)),simNg(1:size(simNi,1)-params.Reff-params.Rest,:),'linewidth',2,'color',[0 0 1]);
    end
end
if ~SumHigh
    if ~NoPlotReff && ~NoPlotRest
        h2 = plot(simulation(1,:)/(1+(PlotMinutes*59)),simNi,'linewidth',2,'color',[0 0.7 0.7]);
    elseif ~NoPlotRest % && NoPlotReff
        h2 = plot(simulation(1,:)/(1+(PlotMinutes*59)),simNi(1:size(simNi,1)-params.Reff,:),'linewidth',2,'color',[0 0.7 0.7]);
    elseif ~NoPlotReff % && NoPlotRest
        if params.Reff && params.Rest
            h2 = plot(simulation(1,:)/(1+(PlotMinutes*59)),simNi([1:size(simNi,1)-2 size(simNi,1)],:),'linewidth',2,'color',[0 0.7 0.7]);
        elseif params.Reff % && ~params.Rest
            h2 = plot(simulation(1,:)/(1+(PlotMinutes*59)),simNi,'linewidth',2,'color',[0 0.7 0.7]);
        else % i.e. ~params.Reff
            h2 = plot(simulation(1,:)/(1+(PlotMinutes*59)),simNi(1:size(simNi,1)-params.Rest,:),'linewidth',2,'color',[0 0.7 0.7]);
        end
    else % i.e. NoPlot either
        h2 = plot(simulation(1,:)/(1+(PlotMinutes*59)),simNi(1:size(simNi,1)-params.Reff-params.Rest,:),'linewidth',2,'color',[0 0.7 0.7]);
    end
else
    if ~NoPlotReff && ~NoPlotRest
        h2 = plot(simulation(1,:)/(1+(PlotMinutes*59)),sum(simNi,1),'linewidth',2,'color',[0 0.7 0.7]);
    elseif ~NoPlotRest % && NoPlotReff
        h2 = plot(simulation(1,:)/(1+(PlotMinutes*59)),sum(simNi(1:size(simNi,1)-params.Reff,:),1),'linewidth',2,'color',[0 0.7 0.7]);
    elseif ~NoPlotReff % && NoPlotRest
        if params.Reff && params.Rest
            h2 = plot(simulation(1,:)/(1+(PlotMinutes*59)),sum(simNi([1:size(simNi,1)-2 size(simNi,1)],:),1),'linewidth',2,'color',[0 0.7 0.7]);
        elseif params.Reff % && ~params.Rest
            h2 = plot(simulation(1,:)/(1+(PlotMinutes*59)),sum(simNi,1),'linewidth',2,'color',[0 0.7 0.7]);
        else % i.e. ~params.Reff
            h2 = plot(simulation(1,:)/(1+(PlotMinutes*59)),sum(simNi(1:size(simNi,1)-params.Rest,:),1),'linewidth',2,'color',[0 0.7 0.7]);
        end
    else % i.e. NoPlot either
        h2 = plot(simulation(1,:)/(1+(PlotMinutes*59)),sum(simNi(1:size(simNi,1)-params.Reff-params.Rest,:),1),'linewidth',2,'color',[0 0.7 0.7]);
    end
end
if PlusRock3; plot(simulation(1,:)/(1+(PlotMinutes*59)),rock3s,'linewidth',2,'color',[0 0 0]); end
if Extended>0
    if Extended==1
        h4 = plot(simulation(1,:)/(1+(PlotMinutes*59)),simTot,'linewidth',1,'color',[1 0 1]);
        h5 = plot(simulation(1,:)/(1+(PlotMinutes*59)),simNw(1:size(simNi,1)-params.Reff-params.Rest,:),'--','linewidth',1,'color',[0 0 0]);  
    end
        apparentRem = simTot+sum(simNw(1:size(simNw,1)-params.Reff-params.Rest,:),1);
        apparentRem(simulation(1,:)>60*wait_time) = NaN;
    h6 = plot(simulation(1,:)/(1+(PlotMinutes*59)),apparentRem,'linewidth',2,'color',[1 0 1]);
end
if PlotMinutes; xlabel('Time (min)'); else xlabel('Time (s)'); end

hold off

yl = get(gca,'ylim'); ylim([0 yl(2)]); set(gca,'fontsize',16);
