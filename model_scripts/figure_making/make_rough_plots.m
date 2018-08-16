function make_rough_plots(simulation,simNi,simNg,simNr,simNw,mff,goal,measdat,measdatA,timres,params);

%%
% Scripts for making 5-8 plots for quickly checking out the model resultss

%%
c = params.c;
wait_time = params.wait_time;
Glu = params.Glu;
Rock = params.Rock;

%% signal vs temperature ("thermograms")

% absolute
    warning('off','all');
figure(1);
    scaler = 1;
hold on;
plot(simulation(2,:)-273.15,simulation(3:end-params.Reff-params.Rest,:)*scaler,':r');
plot(simulation(2,:)-273.15,sum(simulation(3:end-params.Reff-params.Rest,:),1)*scaler);
plot(goal(:,2)-273.15,goal(:,3),'linewidth',1,'color','k');
title(c)
hold off
    if params.Rest
        figure(191);
        hold on;
        plot(simulation(2,:)-273.15,simulation(end-params.Reff,:)*scaler,':r');
        plot(simulation(2,:)-273.15,sum(simulation(end-params.Reff,:),1)*scaler);
            [xu,U] = unique(measdat(:,2)); U = find(~isnan(xu)); xu = xu(U);
            y = measdat(U,3);
            x2 = measdatA(:,size(measdatA,2)-1); X2 = find(~isnan(x2)); x2 = x2(~isnan(x2));
        plot(measdatA(X2,2),(measdatA(X2,3)-interp1(xu,y,x2))*scaler,'linewidth',1,'color','k');
        hold off
        title('"Rest"')
    end

% normalized
figure(11)
    PlotRes = 0; % show some sort of residual
hold on
h1 = plot(goal(:,2)-273.15,goal(:,3)/max(medfilt1(goal(:,3),mff)),'linewidth',1,'color','k');
    xlim([25 200]); ylim([0 1]); set(gca,'fontsize',16); ylabel('Normalized signal'); xlabel('Desorption temperature (C)'); grid on;
    for i=3:size(simulation,1)-params.Reff-params.Rest
        plot(simulation(2,:)-273.15,simulation(i,:)/max(sum(simulation(3:end-params.Reff-params.Rest,:),1)),'-r','color',[0.8 0.7 0]);
    end
h2 = plot(simulation(2,:)-273.15,sum(simulation(3:end-params.Reff-params.Rest,:),1)/max(sum(simulation(3:end-params.Reff-params.Rest,:),1)));
if PlotRes && isreal(simulation)
    grn_y = sum(simulation(3:end-params.Reff-params.Rest,:),1)/max(sum(simulation(3:end-params.Reff-params.Rest,:),1));
    blk_y = goal(:,3)/max(medfilt1(goal(:,3),mff));
    [grn_xU,GU] = unique(simulation(2,:)-273.15);
    grn_xU = grn_xU(~isnan(grn_xU));
    y = grn_y(GU); y = y(~isnan(grn_xU));
    blk_x = goal(:,2)-273.15;
    resid_y = blk_y-interp1(grn_xU,y,blk_x);
    h3 = plot(blk_x,resid_y/max(resid_y),'--');
end
    ylim([0 1]); set(h2,'linewidth',1); set(h2,'color',[0 0.8 0]); legend([h1 h2],'Exp. data','Model')
    if PlotRes; set(h3,'color',[1 0 1]); legend([h1 h2 h3],'Exp. data','Model','"Residual"'); end
hold off
title(c)
    if params.Rest
        figure(1191)
        hold on
                [xu,U] = unique(measdat(:,2)); U = find(~isnan(xu)); xu = xu(U);
                y = measdat(U,3);
                x2 = measdatA(:,2); X2 = find(~isnan(x2)); x2 = x2(~isnan(x2));
            y = measdatA(X2,3)-interp1(xu,y,x2);
        h1 = plot(measdatA(X2,2),y/max(medfilt1(y,mff)),'linewidth',1,'color','k');
            xlim([25 200]); ylim([0 1]); set(gca,'fontsize',16); ylabel('Normalized signal'); xlabel('Desorption temperature (C)'); grid on;
            for i=size(simulation,1)-params.Reff
                plot(simulation(2,:)-273.15,simulation(i,:)/max(sum(simulation(end-params.Reff,:),1)),'-r','color',[0.8 0.7 0]);
            end
        h2 = plot(simulation(2,:)-273.15,sum(simulation(end-params.Reff,:),1)/max(sum(simulation(end-params.Reff,:),1)));
        if PlotRes && isreal(simulation)
            grn_y = sum(simulation(end-params.Reff,:),1)/max(sum(simulation(end-params.Reff,:),1));
            blk_y = goal(:,3)/max(medfilt1(goal(:,3),mff));
            [grn_xU,GU] = unique(simulation(2,:)-273.15);
            grn_xU = grn_xU(~isnan(grn_xU));
            y = grn_y(GU); y = y(~isnan(grn_xU));
            blk_x = goal(:,2)-273.15;
            resid_y = blk_y-interp1(grn_xU,y,blk_x);
            h3 = plot(blk_x,resid_y/max(resid_y),'--');
        end
            ylim([0 1]); set(h2,'linewidth',1); set(h2,'color',[0 0.8 0]); legend([h1 h2],'Exp. data','Model')
            if PlotRes; set(h3,'color',[1 0 1]); legend([h1 h2 h3],'Exp. data','Model','"Residual"'); end;    
        hold off
        title('"Rest"')
    end
    warning('on','all');

%% signals vs. time

% absolute
    warning('off','all');
figure(2);
    scaler = 1;
hold on;
plot(goal(:,1),goal(:,3),'k');
plot(simulation(1,:),simulation(3:end-params.Rest-params.Reff,:)*scaler,':r');
plot(simulation(1,:),sum(simulation(3:end-params.Rest-params.Reff,:),1)*scaler);
hold off
title(c)
if params.Rest
    figure(291)
    hold on;
        x = 0:timres:timres*(size(measdatA,1)-1);
    plot(x,measdatA(:,3)*scaler,'k');
    plot(simulation(1,:),simulation(end-params.Reff,:)*scaler,':r');
    plot(simulation(1,:),sum(simulation(end-params.Reff,:),1)*scaler);
    title('"Rest"');
end

% normalized
figure(22)
hold on
h1 = plot(goal(:,1),goal(:,3)/max(medfilt1(goal(:,3),mff)),'color','k');
    for i=3:size(simulation,1)-params.Rest-params.Reff
        plot(simulation(1,:),simulation(i,:)/max(sum(simulation(3:end-params.Rest-params.Reff,:),1)),'-r','color',[0.8 0.7 0]);
    end
h2 = plot(simulation(1,:),sum(simulation(3:end-params.Rest-params.Reff,:),1)/max(sum(simulation(3:end-params.Rest-params.Reff,:),1)));
ylim([0 1]); set(gca,'fontsize',16); ylabel('Normalized signal'); xlabel('Time (s)'); grid on;
hold off
title(c)
    if params.Rest
        figure(2291)
        hold on
        h1 = plot(0:timres:timres*(size(measdatA,1)-1),measdatA(:,end)/max(medfilt1(measdatA(:,end),mff)),'color','k');
            for i=size(simulation,1)-params.Reff
                plot(simulation(1,:),simulation(i,:)/max(simulation(i,:)),'-r','color',[0.8 0.7 0]);
            end
            ylim([0 1]); set(gca,'fontsize',16); ylabel('Normalized signal'); xlabel('Time (s)'); grid on;
        hold off
        title('"Rest"')
    end
warning('on','all');

%% Number of molecules vs. time

    figure(201);
    close; % always make new plot
if params.Rest
    figure(201)
        PlotMinutes = 1; % else s
        Extended = 1; % 1 ... show Nw, Tot, apparentRem, 2 ... same but not Nw
        PlusRock3 = 0; % 2 ... fine-tuning tests (tried for avgs2bLow only so far)
        if wait_time==0; WaitOnly=0; else; WaitOnly = 1; end
        SumHigh = 0; % sum up simNi (for a paper plot)
        if ~Glu; NoPlotNg = 1; % disregard Ng (for a paper plot, where it's just 0)
        else; NoPlotNg = 0; end;
        if Rock>0; NoPlotNr = 0;
        else; NoPlotNr = 1; end; % disregard Nr (for a paper plot, where it's just 0)
        NoPlotReff = 1; % disregard Nreff if used
        NoPlotRest = 1; % disregard Nrest if used
        ShowInclRest = 1; % plot Nrest separately (irrespective of NoPlotRest
    N_vs_t_plots
end
figure(2011);
close;figure(2011) % always make new plot
    PlotMinutes = 1; % else s
    Extended = 1; % 1 ... show Nw, Tot, apparentRem, 2 ... same but not Nw
    PlusRock3 = 0; % 2 ... fine-tuning tests (tried for avgs2bLow only so far)
    if wait_time==0; WaitOnly=0; else; WaitOnly = 1; end
    SumHigh = 0; % sum up simNi (for a paper plot)
    if ~Glu; NoPlotNg = 1; % disregard Ng (for a paper plot, where it's just 0)
    else; NoPlotNg = 0; end;
    if Rock>0; NoPlotNr = 0;
    else; NoPlotNr = 1; end; % disregard Nr (for a paper plot, where it's just 0)
    NoPlotReff = 1; % disregard Nreff if used
    NoPlotRest = 1; % disregard Nrest if used
    ShowInclRest = 0; % plot Nrest separately (irrespective of NoPlotRest
N_vs_t_plots

end