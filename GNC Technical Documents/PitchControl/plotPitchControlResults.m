function files = plotPitchControlResults(r,outputDirectory)
%PLOTPITCHCONTROLRESULTS Plot in chronological order, including Mach axes.
folder = fullfile(outputDirectory,'plots');
mkdir(folder);
files = strings(0,1);
quantities = ["M_vehicle_Nm","M_trim_required_Nm","M_maneuver_required_Nm", ...
    "DeltaCm_trim","DeltaCm_maneuver"];
titles = ["Vehicle aerodynamic pitching moment","Required static-trim control moment", ...
    "Required dynamic control moment","Required static moment coefficient","Required dynamic moment coefficient"];
notice = 'Provisional reference inputs - see Inputs_Verification';
if ~startsWith(r.status,"PROVISIONAL"), notice = 'References user-confirmed - review interpolation diagnostics'; end
for q = 1:numel(quantities)
    for start = 1:r.config.maxCasesPerPlot:height(r.cases)
        last = min(start+r.config.maxCasesPerPlot-1,height(r.cases));
        fig = newFigure();
        closer = onCleanup(@() close(fig));
        ax = axes('Parent',fig); hold(ax,'on'); grid(ax,'on');
        colors = lines(last-start+1);
        for k = start:last
            H = r.caseHistory(r.caseHistory.CaseID == r.cases.CaseID(k),:);
            x = H.Time_s;
            if q >= 4, x = H.Mach; end
            label = sprintf('%s: xCG %.6g m, Iyy %.6g kg m^2', ...
                r.cases.CaseLabel(k),r.cases.xCG_m(k),r.cases.Iyy_kg_m2(k));
            plot(ax,x,H.(quantities(q)),'LineWidth',1.5,'Color',colors(k-start+1,:), ...
                'DisplayName',label);
            [~,peak] = max(abs(H.(quantities(q))));
            plot(ax,x(peak),H.(quantities(q))(peak),'o','Color',colors(k-start+1,:), ...
                'HandleVisibility','off');
        end
        title(ax,{char(titles(q)),notice},'Interpreter','none');
        if q >= 4
            xlabel(ax,'Mach (trajectory order retained)'); ylabel(ax,'Required moment coefficient');
        else
            xlabel(ax,'Time (s)'); ylabel(ax,'Pitch moment (N m)');
        end
        legend(ax,'Location','southoutside','Interpreter','none','FontSize',9);
        padY(ax);
        files = [files;saveBoth(fig,folder,sprintf('%s_cases_%d_%d',quantities(q),start,last))]; %#ok<AGROW>
        clear closer
    end
end
fig = newFigure(); closer = onCleanup(@() close(fig));
layout = tiledlayout(fig,3,1,'TileSpacing','compact');
motion = {'PitchAttitude_deg','PitchRate_rad_s','PitchAccel_rad_s2'};
labels = {'Pitch attitude (deg)','Pitch rate (rad/s)','Pitch acceleration (rad/s^2)'};
for k = 1:3
    ax = nexttile(layout); plot(ax,r.common.Time_s,r.common.(motion{k}),'LineWidth',1.5);
    grid(ax,'on'); ylabel(ax,labels{k}); xlabel(ax,'Time (s)');
    padY(ax);
end
title(layout,{'Spline-derived pitch motion',notice},'Interpreter','none');
files = [files;saveBoth(fig,folder,'pitch_motion')]; clear closer

fig = newFigure(); closer = onCleanup(@() close(fig));
layout = tiledlayout(fig,2,1,'TileSpacing','compact');
names = {'Alpha_deg','FlightPath_deg'};
labels = {'Angle of attack (deg)','Flight-path angle (deg)'};
for k = 1:2
    ax = nexttile(layout);
    plot(ax,r.common.Time_s,r.common.(names{k}),'LineWidth',1.5,'DisplayName','Smoothing spline');
    hold(ax,'on');
    plot(ax,r.nodeDiagnostics.Time_s,r.nodeDiagnostics.(names{k}),'o', ...
        'MarkerSize',4,'DisplayName','Original collocation nodes');
    grid(ax,'on'); ylabel(ax,labels{k}); xlabel(ax,'Time (s)'); legend(ax,'Location','best');
    padY(ax);
end
title(layout,{'Interpolation fit at original trajectory nodes',notice},'Interpreter','none');
files = [files;saveBoth(fig,folder,'spline_fit')]; clear closer

if ~isempty(r.sensitivity)
    variants = unique(r.sensitivity.Variant,'stable');
    first = r.sensitivity(r.sensitivity.CaseID == r.cases.CaseID(1),:);
    for start = 1:r.config.maxCasesPerPlot:height(r.cases)
        last = min(start+r.config.maxCasesPerPlot-1,height(r.cases));
        fig = newFigure(); closer = onCleanup(@() close(fig));
        layout = tiledlayout(fig,2,1,'TileSpacing','compact');
        ax = nexttile(layout);
        plot(ax,1:height(first),first.MaxAbsPitchAccel_rad_s2,'o-','LineWidth',1.5);
        xticks(ax,1:numel(variants)); xticklabels(ax,variants); grid(ax,'on'); padY(ax);
        ylabel(ax,'Max |pitch acceleration| (rad/s^2)');
        ax = nexttile(layout); hold(ax,'on');
        for k = start:last
            S = r.sensitivity(r.sensitivity.CaseID == r.cases.CaseID(k),:);
            plot(ax,1:height(S),S.PeakDynamicMagnitude_Nm,'o-','DisplayName',r.cases.CaseID(k));
        end
        xticks(ax,1:numel(variants)); xticklabels(ax,variants); grid(ax,'on'); padY(ax);
        ylabel(ax,'Max |dynamic control moment| (N m)');
        legend(ax,'Location','eastoutside','Interpreter','none');
        title(layout,{'Fixed interpolation sensitivity variants; baseline unchanged',notice},'Interpreter','none');
        name = 'interpolation_sensitivity';
        if height(r.cases) > r.config.maxCasesPerPlot
            name = sprintf('%s_cases_%d_%d',name,start,last);
        end
        files = [files;saveBoth(fig,folder,name)]; %#ok<AGROW>
        clear closer
    end
end
end

function padY(ax)
limits = ylim(ax);
padding = 0.04*diff(limits);
ylim(ax,limits+[-padding padding]);
end

function fig = newFigure()
fig = figure('Visible','off','Color','w','Position',[100 100 1100 750]);
end

function paths = saveBoth(fig,folder,name)
png = fullfile(folder,[name '.png']);
figFile = fullfile(folder,[name '.fig']);
exportgraphics(fig,png,'Resolution',150);
savefig(fig,figFile);
paths = [string(png);string(figFile)];
end
