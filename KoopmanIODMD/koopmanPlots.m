function varargout = koopmanPlots(task, varargin)
%KOOPMANPLOTS Plots and tables of the Koopman model identification
%   varargout = koopmanPlots(task, ...) plots the identification and
%   validation results of the Koopman models (koopmanEstimation).
%
%   Tasks (original function names in brackets):
%   - 'effWind' (plotEDMDinputsEffWind): inputs, effective wind speeds or
%     farm power of WFSim and the Koopman model
%   - 'effWindUpdate' (plotEDMDinputsEffWindUpdate): as 'effWind' with the
%     online updated model
%   - 'powerStatesDiag' (plotPowerStatesDiag): diagnostic plots of the
%     input-only power models (thesis figures 03aFig_WRC_*)
%   - 'aicNoV1Sigma' (plotAICnoV1Sigma): singular values of the online
%     updated models of the adaptive AIC
%   - 'aicNoV1IdVal' (plotAICnoV1IdVal): adaptive AIC model on its
%     identification and validation set
%   - 'aicNoV1Validation' (validateAICnoV1): open-loop validation of the
%     non-adaptive and the online updated model (wind ramp 7 to 9 m/s)
%   - 'aicNoV1Test9' (testAICnoV1Vinf9): offline and learned models on
%     unseen data at 9 m/s
%   - 'wrcModelTable' (writeWRCModelTable): LaTeX table of the WRC-AIC
%     model comparison
switch task
    case 'effWind',           [varargout{1:nargout}] = plotEDMDinputsEffWind(varargin{:});
    case 'effWindUpdate',     [varargout{1:nargout}] = plotEDMDinputsEffWindUpdate(varargin{:});
    case 'powerStatesDiag',   [varargout{1:nargout}] = plotPowerStatesDiag(varargin{:});
    case 'aicNoV1Sigma',      [varargout{1:nargout}] = plotAICnoV1Sigma(varargin{:});
    case 'aicNoV1IdVal',      [varargout{1:nargout}] = plotAICnoV1IdVal(varargin{:});
    case 'aicNoV1Validation', [varargout{1:nargout}] = validateAICnoV1(varargin{:});
    case 'aicNoV1Test9',      [varargout{1:nargout}] = testAICnoV1Vinf9(varargin{:});
    case 'wrcModelTable',     [varargout{1:nargout}] = writeWRCModelTable(varargin{:});
    otherwise, error('koopmanPlots:task', 'Unknown task %s', task);
end
end

%% ======================================================================
%% plotEDMDinputsEffWind  (task 'effWind')
function fig1 = plotEDMDinputsEffWind(ysim_val,Inputs_val,Deterministic_val,FITje_val,dirFig,Vin,nx,strVal,displayOff,Vinf,strOut,tVec,interp,ysimRef,refName,estName,FITjeRef)
% plotEDMDinputsEffWind plots the simulated outputs of the eDMD estimations
% and prints them into the directory given in dirFig.
% Optional tVec: time (s) of the samples, x axis 'Time (s)' (default: sample index)
% Optional interp: 'latex' (default) or 'tex' (labels converted: no $, units
% in (), no legends: the first title names the measured (WT1, WT2) and the
% estimated signals in their colors)
% Optional ysimRef: estimates of a reference model (same size as ysim_val),
% plotted as grey dashed lines; refName: its name in the legend/title key
% Optional estName, FITjeRef (tex): names of the model and the reference
% model and the VAF of the reference in the title key ('Estimate P_{WF}:' or
% 'Estimate V_r:' with the names and the VAF in the line colors)

if nargin < 6
    Vin = 0;
end
if nargin < 7
    nx = 0;
end

if nargin < 8
    strVal = 'Val.';
end

if nargin < 9
    displayOff = 0;
end

if nargin < 10
    Vinf = [];
end

if nargin < 11
    strOut = 'U';
end

if nargin < 12
    tVec = [];
end

if nargin < 13
    interp = 'latex';
end

if nargin < 14
    ysimRef = [];
end

if nargin < 15 || isempty(refName)
    refName = 'Ref.';
end

if nargin < 16
    estName = '';
end

if nargin < 17
    FITjeRef = [];
end

if ~isempty(ysimRef) && size(ysimRef,2) < size(ysimRef,1)
    ysimRef = ysimRef';
end

if nx > 0
    strK =  sprintf(', $n_K$ = %d', nx);
else
    strK = '';
end


if size(ysim_val,2) < size(ysim_val,1)
    ysim_val = ysim_val';
end


% (path: see WFSim_addpaths) % setFigureSize, printFigure
noSubIn = size(Inputs_val,1);
noSubOut =  min(size(ysim_val));
noSub = noSubOut + noSubIn -1;
fig1 = figure;
tl = tiledlayout(fig1, noSub, 1, 'TileSpacing', 'compact', 'Padding', 'compact'); % less white space than subplot
hFac = 1; % height: 5 subplots 1.6 (as the closed-loop plots), 4 subplots 1.3, 3 subplots 1.15 times the standard figure
if noSub == 5, hFac = 1.6; elseif noSub == 4, hFac = 1.3; elseif noSub == 3, hFac = 1.15; end
if hFac > 1
    posDefault = get(gcf, 'position');
    set(gcf, 'position', [posDefault(1),posDefault(2)-150,posDefault(3),posDefault(4)*hFac]);
end
setFigureSize(fig1, [1 hFac]); % export size independent of the screen
cl = lines;
lw = 1;
greyRef = [0.6 0.6 0.6]; % reference model estimates (dashed)
lwRef = 1.7; % thicker than the estimate, so that the grey line stays visible
fs = getTimeSeriesFontSize(); % default 13 (runGenerateAllThesisFigures)
kplot = 60: length(ysim_val);
if isempty(tVec)
    xk = 1:length(kplot);
else % time axis: the data start after the start-up (t0), show all (as the VAF)
    kplot = 1:length(ysim_val);
    xk = tVec(kplot);
end

nexttile(tl,1)
plot(xk,Inputs_val(1,kplot),'LineWidth',lw,'Color', cl(1,:)); %1
hold on;
plot(xk,Inputs_val(2,kplot),'LineWidth',lw,'Color', cl(2,:)) %3
grid on; axis tight; ylim([0 2+0.2])
ylabel('Thrust $C_T^\prime$ [-]','interpreter','latex','Fontsize',fs); % name and symbol as in chapter 4

if displayOff && length(FITje_val) >1
    if ~isappdata(groot, 'ShowThesisTitle') || getappdata(groot, 'ShowThesisTitle') % title switch of runGenerateAllThesisFigures
        title(sprintf('%s:  VAF %sr1: %d, r2: %d',...
            strVal,strOut,round(FITje_val(1:2))),'interpreter','latex','Fontsize',fs);
    end
elseif length(FITje_val)>1
    if ~isappdata(groot, 'ShowThesisTitle') || getappdata(groot, 'ShowThesisTitle') % title switch of runGenerateAllThesisFigures
        title(sprintf('%s Set%s:  VAF %sr1: %d, r2: %d',...
            strVal,strK,strOut,round(FITje_val(1:2))),'interpreter','latex','Fontsize',fs);
    end
elseif displayOff
    if ~isappdata(groot, 'ShowThesisTitle') || getappdata(groot, 'ShowThesisTitle') % title switch of runGenerateAllThesisFigures
        title(sprintf('%s:  VAF %sT: %d',...
            strVal,strOut,round(FITje_val)),'interpreter','latex','Fontsize',fs);
    end
else
    if ~isappdata(groot, 'ShowThesisTitle') || getappdata(groot, 'ShowThesisTitle') % title switch of runGenerateAllThesisFigures
        title(sprintf('%s Set%s:  VAF %sT: %d',...
            strVal,strK,strOut,round(FITje_val)),'interpreter','latex','Fontsize',fs);
    end
end


strLeg = sprintf('WT%d,',1:2);
legend(regexp(strLeg(1:end-1), ',', 'split'),'interpreter','latex','Fontsize',fs, 'Orientation','horizontal', 'location','southeast');
%legend('boxoff')
set(gca,'fontsize', fs)

if noSubIn > 2
    nexttile(tl,2)
    if isempty(Vinf)
        plot(xk,Inputs_val(3,kplot),'LineWidth',lw,'Color', cl(1,:)); %1
    else
        plot(xk,Vinf(kplot),'k','LineWidth',lw); hold on; plot(xk,Inputs_val(3,kplot),'LineWidth',lw,'Color', cl(1,:));
    end
    grid on; axis tight
    if ~isempty(Vinf) % V_inf line not on the axis limit
        ylim(ylim + [-0.2 0.2]);
    end
    if Vin && ~isempty(Vinf)
        ylabel('Wind $V$ [m/s]','interpreter','latex','Fontsize',fs);
        legend('$V_{\inf}$','$V_1$','interpreter','latex','Fontsize',fs,'orientation', 'horizontal','Location','NorthWest');
        %legend('boxoff')
    elseif (Vin && noSub == 4) || (~Vin && noSub == 3) 
         ylabel('Yaw $\gamma$ [deg]','interpreter','latex','Fontsize',fs);       
    else
        ylabel('Wind $V_1$ [m/s]','interpreter','latex','Fontsize',fs);
    end
    set(gca,'fontsize', fs)
    
    if noSubIn > 3 %&& noSub ~= 4
        nexttile(tl,3)
        plot(xk,Inputs_val(4,kplot),'LineWidth',lw,'Color', cl(1,:)); %1
        ylabel('Wind $V_1$ [m/s]','interpreter','latex','Fontsize',fs);
        grid on; axis tight
    end
    
end

if noSubOut == 1
    strAdd = {'_{WF}$ [MW]'};
else
    strAdd = {'_{r1}$ [m/s]','_{r2}$ [m/s]'};
end

for idx = 1: noSubOut
    idxT = 2-mod(idx,2);
    nexttile(tl,noSub - (noSubOut-idx))
    plot(xk,Deterministic_val(idx,kplot),'LineWidth',lw,'Color', cl(idxT,:)); %1
    hold on;
    if ~isempty(ysimRef) % reference model (grey) below the estimate
        plot(xk,ysimRef(idx,kplot),'Color',greyRef,'LineWidth',lwRef,'LineStyle','--')
    end
    plot(xk,ysim_val(idx,kplot),'Color',[0,0,0],'LineWidth',lw,'LineStyle','--') %3
    grid on; axis tight;
    pos = axis; axis([pos(1:3), pos(4)+0.3]);

    %xlabel('k [-]','interpreter','latex')
    sigName = 'Wind';
    if noSubOut == 1, sigName = 'Power'; end
    ylabel([sigName,' $',strOut,strAdd{idxT}],'interpreter','latex')
    if isempty(ysimRef)
        legend({'Real','Est.'},'Location','northwest','interpreter','latex','fontsize',fs,'orientation','horizontal')
    else
        legend({'Real',refName,'Est.'},'Location','northwest','interpreter','latex','fontsize',fs,'orientation','horizontal')
    end
    %legend('boxoff')
    set(gca,'fontsize', fs)
    
end
if ~isempty(tVec)
    xlabel('Time (s)','interpreter','latex','Fontsize',fs);
end

if strcmpi(interp,'tex') % convert all LaTeX strings to TeX
    for hTxt = findall(fig1,'Type','text')'
        set(hTxt,'String',latex2tex(get(hTxt,'String')),'Interpreter','tex');
    end
    delete(findall(fig1,'Type','legend'));
    axs = findall(fig1,'Type','axes');
    set(axs,'TickLabelInterpreter','tex');
    axTop = axs(end); % first subplot
    if noSubOut == 1 % estimated output(s), black dashed
        estSig = [strOut,'_{WF}'];
    else
        estSig = 'V_r';
    end
    measNames = {'WT1','WT2'}; measColors = [cl(1,:); cl(2,:)];
    if Vin && ~isempty(Vinf) % free stream wind speed (black) in the wind subplot
        measNames{end+1} = 'V_\infty'; measColors(end+1,:) = [0 0 0];
    end
    strMeas = ['Measured sig.: ', colorKey(measNames, measColors)];
    if isempty(ysimRef)
        axTop.Title.String = {axTop.Title.String, ...
            [strMeas, '; ', colorKey({['Estimate ',estSig]}, [0 0 0])]};
    else % model (black) and reference model (grey) with their VAF
        vafStr = @(v) ['(VAF ', strjoin(arrayfun(@(x) sprintf('%d',round(x)), v, 'UniformOutput', false), '/'), ')'];
        estKey = {strtrim([estName,' ',vafStr(FITje_val)]), refName};
        if ~isempty(FITjeRef), estKey{2} = [refName,' ',vafStr(FITjeRef)]; end
        axTop.Title.String = {[latex2tex(strVal), ': ', strMeas], ...
            ['Estimate ', estSig, ': ', colorKey(estKey, [0 0 0; greyRef])]};
    end
    if isappdata(groot, 'ShowThesisTitle') && ~getappdata(groot, 'ShowThesisTitle') % title switch of runGenerateAllThesisFigures
        axTop.Title.String = '';
    end
end
setFigureFont(fig1, fs); % all text 12 pt, bold titles
setTimeLabelFontSize(fig1); % time axis label 14 (runGenerateAllThesisFigures)
for ax = findall(fig1,'Type','axes')' % time in s without x10^4 (long identification sets): no extra row
    ax.XAxis.Exponent = 0;
end

%sys_red{1,1}.Xint = xo;
strVal1 = matlab.lang.makeValidName(strVal);
if strcmpi(interp,'tex') % thesis figures: y axes not completely tight
    padYLim(fig1);
end
% x tick labels only at the bottom, same small gaps between the tiles
% (findall: last created axes first)
tightenTiles(flipud(findall(fig1,'Type','axes')));
printFigure(fig1, fullfile(dirFig,['Val',strVal1])); % fixed size (setFigureSize)
end

function str = latex2tex(str)
% latex2tex LaTeX label to TeX: no $, \inf -> \infty, units [..] -> (..),
% U_r -> V_r, C_T^\prime -> C_T'
if iscell(str), str = cellfun(@latex2tex, str, 'UniformOutput', false); return; end
str = strrep(str,'$','');
str = strrep(str,'^\prime','''');
str = strrep(str,'\inf}','\infty}');
str = regexprep(str,'\[([^\]]*)\]','($1)');
str = strrep(str,'(deg)',['(',char(176),')']); % degree sign as in the qLMPC plots
str = regexprep(str,'\<Ur(\d)','V_{r$1}'); % VAF Ur1 -> V_{r1}
str = regexprep(str,'VAF PWF:','VAF P_{WF}:'); % farm power as in the y label
str = regexprep(str,', r(\d):',', V_{r$1}:'); % VAF ..., r2 -> V_{r2}
str = regexprep(str,'\<U(_?\{?r)','V$1'); % effective wind speed V_r as in the thesis
end

function str = colorKey(names, colors)
% colorKey signal names in their line colors (TeX), used instead of a legend
parts = cell(size(names));
for idx = 1:numel(names)
    parts{idx} = sprintf('{\\color[rgb]{%.3f,%.3f,%.3f}%s}', colors(idx,:), names{idx});
end
str = strjoin(parts, ', ');
end

%% ======================================================================
%% plotEDMDinputsEffWindUpdate  (task 'effWindUpdate')
function fig1 = plotEDMDinputsEffWindUpdate(ysim_val,Inputs_val,Deterministic_val,FITje_val,dirFig,Vin,nx,strVal,Pvec,displayOff,Vinf,tVec,interp,marks,figName)
% plotEDMDinputsEffWind plots the simulated outputs of the eDMD estimations
% and prints them into the directory given in dirFig.
% Optional tVec: time (s) of the samples, x axis 'Time (s)' (default: sample index)
% Optional interp: 'latex' (default) or 'tex' (labels converted: no $, units
% in (), no legends: the first title names the measured signals WT1, WT2,
% V_inf and the estimated signals V_r1,est, V_r2,est in their colors)
% Optional marks: table with Time_s and Color (plotAICnoV1Sigma), colored
% lines at these times in the update subplot (tVec required)
% Optional figName: file name without folder and extension (default
% ValP<strVal as valid name>)

if nargin < 6
    Vin = 0;
end
if nargin < 7
    nx = 0;
end

if nargin < 8
    strVal = 'Val.';
end

if nargin < 9
    displayOff = 0;
end

if nargin < 11
    Vinf = [];
end

if nargin < 12
    tVec = [];
end

if nargin < 13
    interp = 'latex';
end

if nargin < 14
    marks = [];
end

if nargin < 15
    figName = '';
end

if nx > 0
    strK =  sprintf(', $n_K$ = %d', nx);
else
    strK = '';
end

% (path: see WFSim_addpaths) % setFigureSize, printFigure
noSub = size(Inputs_val,1) + 2;
fig1 = figure;
tl = tiledlayout(fig1, noSub, 1, 'TileSpacing', 'compact', 'Padding', 'compact'); % less white space than subplot
posDefault = get(gcf, 'position');
hFac = 1.3; % 4 subplots: 1.3 times the standard figure
if noSub == 5 % as the closed-loop plots (plotWFSimTimeseriesTex)
    hFac = 1.6;
end
set(gcf, 'position', [posDefault(1),posDefault(2)-150,posDefault(3),posDefault(4)*hFac]);
setFigureSize(fig1, [1 hFac]); % export size independent of the screen
cl = colormap('lines');
lw = 1;
fs = getTimeSeriesFontSize(); % default 13 (runGenerateAllThesisFigures)
if isempty(tVec)
    kplot = 60: length(ysim_val);
    xk = 1:length(kplot);
else % time axis: the data start after the start-up (t0), show all
    kplot = 1: length(ysim_val);
    xk = tVec(kplot);
end

nexttile(tl,1)
plot(xk,Inputs_val(1,kplot),'LineWidth',lw,'Color', cl(1,:)); %1
hold on;
plot(xk,Inputs_val(2,kplot),'LineWidth',lw,'Color', cl(2,:)) %3
grid on; axis tight; ylim([0 2+0.2])
ylabel('Thrust $C_T^\prime$ [-]','interpreter','latex','Fontsize',fs); % name and symbol as in chapter 4

if displayOff
    if ~isappdata(groot, 'ShowThesisTitle') || getappdata(groot, 'ShowThesisTitle') % title switch of runGenerateAllThesisFigures
        title(sprintf('%s: VAF Ur1: %d, Ur2: %d',...
        strVal,round(FITje_val)),'interpreter','latex','Fontsize',fs);
    end
else
if ~isappdata(groot, 'ShowThesisTitle') || getappdata(groot, 'ShowThesisTitle') % title switch of runGenerateAllThesisFigures
    title(sprintf('%s Set%s:  VAF Ur1: %d, Ur2: %d',...
        strVal,strK,round(FITje_val)),'interpreter','latex','Fontsize',fs);
end
end


strLeg = sprintf('WT%d,',1:2);
legend(regexp(strLeg(1:end-1), ',', 'split'),'interpreter','latex','Fontsize',fs, 'Orientation','horizontal', 'location','southeast');
legend('boxoff')
set(gca,'fontsize', fs)

if noSub > 3
    nexttile(tl,2)
    if isempty(Vinf)
    plot(xk,Inputs_val(3,kplot),'LineWidth',lw,'Color', cl(1,:)); %1
    else
         plot(xk,Vinf(kplot),'k','LineWidth',lw); hold on; plot(xk,Inputs_val(3,kplot),'LineWidth',lw,'Color', cl(1,:));
    end
    grid on; axis tight
    if Vin && ~isempty(Vinf)
        ylabel('Wind $V$ [m/s]','interpreter','latex','Fontsize',fs);
        legend('$V_{\inf}$','$V_1$','interpreter','latex','Fontsize',fs,'orientation', 'horizontal','Location','NorthWest');
        legend('boxoff')
    elseif Vin
        ylabel('Wind $V_1$ [m/s]','interpreter','latex','Fontsize',fs);
    else
        ylabel('Yaw $\gamma$ [deg]','interpreter','latex','Fontsize',fs);
    end
    set(gca,'fontsize', fs)
end

nexttile(tl,noSub-2)
plot(xk,Deterministic_val(1,kplot),'LineWidth',lw,'Color', cl(1,:)); %1
hold on;
plot(xk,ysim_val(1,kplot),'Color',[0,0,0.6],'LineWidth',lw,'LineStyle','--') %3
grid on; axis tight;
pos = axis; axis([pos(1:3), pos(4)+0.3]);
%xlabel('k [-]','interpreter','latex')
ylabel('Wind $U_{r1}$ [m/s]','interpreter','latex')
lgh = legend({'Real','Est.'},'Location','southeast','interpreter','latex','fontsize',fs,'orientation','horizontal');
%  0.6155    0.4802    0.2702    0.0397
lgh.Position = [lgh.Position(1) + 0.02,lgh.Position(2:4)];
legend('boxoff')
set(gca,'fontsize', fs)

nexttile(tl,noSub-1)
plot(xk,Deterministic_val(2,kplot),'LineWidth',lw,'Color', cl(2,:)); %1
hold on;
plot(xk,ysim_val(2,kplot),'Color',[0.6, 0, 0.0],'LineWidth',lw,'LineStyle','--') %3
grid on; axis tight
pos = axis; axis([pos(1:3), pos(4)+0.3]);
ylabel('Wind $U_{r2}$ [m/s]','interpreter','latex','fontsize', fs)
lgh = legend({'Real','Est.'},'Location','southeast','interpreter','latex','fontsize',fs,'orientation','horizontal');
lgh.Position = [0.2595    0.3585    0.2702    0.0397];
legend('boxoff')
set(gca,'fontsize', fs)

nexttile(tl,noSub)
if isempty(tVec)
    plot(min((Pvec-1),1),'k','LineWidth',lw); %1
else
    plot(xk,min((Pvec(kplot)-1),1),'k','LineWidth',lw);
end
if ~isempty(marks) % updates after which the model is shown elsewhere
    hold on;
    for idx = find(marks.Time_s(:)' >= xk(1) & marks.Time_s(:)' <= xk(end))
        % first sample (initial model): shifted to be visible next to the axis
        tMark = max(marks.Time_s(idx), xk(1) + 0.004*(xk(end)-xk(1)));
        plot(tMark*[1 1], [0 1], 'Color', marks.Color(idx,:), 'LineWidth', 2, 'LineStyle', '-'); % markers, not signals: 2
        % time stamp at the bottom, right of the line (left if TextSide < 0,
        % with more space, the updates before are dense), the magnitude of
        % TextSide scales the distance
        side = 1;
        if ismember('TextSide', marks.Properties.VariableNames), side = marks.TextSide(idx); end
        hAlign = {'right','left'};
        gap = 0.008 + 0.012*(side < 0);
        text(tMark + side*gap*(xk(end)-xk(1)), 0.05, sprintf('%d s', marks.Time_s(idx)), ...
            'Color', marks.Color(idx,:), 'FontSize', 12, 'FontWeight', 'bold', ...
            'VerticalAlignment', 'bottom', 'HorizontalAlignment', hAlign{(side > 0) + 1});
    end
end

grid on; axis tight
pos = axis; axis([pos(1:3), pos(4)+0.3]);
if isempty(tVec)
    xlabel('k [-]','interpreter','latex')
else
    xlabel('Time (s)','interpreter','latex')
end
ylabel('$Update$ ','interpreter','latex')
axis tight;
if ~isempty(tVec) % zero line visible above the axis, as the closed loop
    ylim([-0.1 1.1]);
end

%legend('P','Location','northeast','interpreter','latex','fontsize',fs,'orientation','horizontal');
%legend('boxoff')
set(gca,'fontsize', fs)

if strcmpi(interp,'tex') % convert all LaTeX strings to TeX, color key instead of legends
    for hTxt = findall(fig1,'Type','text')'
        set(hTxt,'String',latex2tex(get(hTxt,'String')),'Interpreter','tex');
    end
    delete(findall(fig1,'Type','legend'));
    axs = findall(fig1,'Type','axes');
    set(axs,'TickLabelInterpreter','tex');
    axTop = axs(end); % first subplot
    axTop.Title.String = {axTop.Title.String, ...
        ['Measured sig.: ', colorKey({'WT1','WT2','V_\infty'}, [cl(1,:); cl(2,:); 0 0 0]), ...
        '; Estimated sig.: ', colorKey({'V_{r1,est}','V_{r2,est}'}, [0 0 0.6; 0.6 0 0])]};
    if isappdata(groot, 'ShowThesisTitle') && ~getappdata(groot, 'ShowThesisTitle') % title switch of runGenerateAllThesisFigures
        axTop.Title.String = '';
    end
    axs(1).YLabel.String = 'Update (-)';
end
setFigureFont(fig1, fs); % all text 12 pt, bold titles
setTimeLabelFontSize(fig1); % time axis label 14 (runGenerateAllThesisFigures)

%sys_red{1,1}.Xint = xo;
strVal1 = matlab.lang.makeValidName(strVal);
if strcmpi(interp,'tex') % thesis figures: y axes not completely tight
    padYLim(fig1);
end
% x tick labels only at the bottom, same small gaps between the tiles
% (findall: last created axes first)
tightenTiles(flipud(findall(fig1,'Type','axes')));
if isempty(figName), figName = ['ValP',strVal1]; end
printFigure(fig1, fullfile(dirFig,figName)); % fixed size (setFigureSize)
end



%% ======================================================================
%% plotPowerStatesDiag  (task 'powerStatesDiag')
function R = plotPowerStatesDiag(modelFile, filenameVal, t0, refModelFile, refName, filenameIdPlot, estName, setsPlot)
% plotPowerStatesDiag diagnostic plots of the open-loop identification of
% the input-only Koopman power model (KoopFeatModel<k>.mat, runKoopStatesCL):
% per identification set (model.filenameId) and validation set one figure
% with the rotor effective wind speeds Ur1,
% Ur2 and one with the farm power, WFSim and the linear model simulated
% from the measured initial state (as testPowerStates), VAF in the legend.
% Plots with plotEDMDinputsEffWind (as MainWfSimVinYawPowerOut), saved in
% DataInOutWfSim/PowerInputsOnly/Diag_<model>_Ur and Diag_<model>_P
%
% Inputs
% - modelFile: default KoopFeatModel3.mat ('Ur n=6')
% - filenameVal: validation sets (cell)
% - t0: first sample (s) as in the identification
% - refModelFile: optional reference model (KoopFeatModel<k>.mat), whose
%   estimates are shown as grey lines (simulated in the same way)
% - refName: name of the reference model in the plots (default 'Ref. model')
% - filenameIdPlot: identification sets (default: all of model.filenameId)
% - estName: name of the model in the title key (with a reference model)
% - setsPlot: sets that are plotted (default: all; {}: none, VAF only)

mainDir = fileparts(fileparts(mfilename('fullpath')));
dirOut = fullfile(mainDir,'DataInOutWfSim','PowerInputsOnly');
if nargin < 1 || isempty(modelFile), modelFile = 'KoopFeatModel3.mat'; end
if nargin < 2 || isempty(filenameVal)
    filenameVal = {'Vinf8dot0_sowfa_2turb_yaw_ct_levels_val','Vinf8dot0_sowfa_2turb_yaw_steps', ...
        'Vinf8dot0_sowfa_2turb_yaw_steps_Ct_comb'};
end
if nargin < 3, t0 = 240; end
if nargin < 4, refModelFile = ''; end
if nargin < 5 || isempty(refName), refName = 'Ref. model'; end
if nargin < 6, filenameIdPlot = []; end
if nargin < 7, estName = ''; end
if nargin < 8, setsPlot = []; end
% (path: see WFSim_addpaths) % powerFeatures, koopFeatState

tmp = load(fullfile(dirOut,modelFile),'model');
model = tmp.model;
sys = ss(model.A, model.B, model.C, model.D, 1);
if ~isempty(refModelFile)
    tmp = load(fullfile(dirOut,refModelFile),'model');
    modelRef = tmp.model;
end
[~, modelName] = fileparts(modelFile);
dirUr = fullfile(dirOut, sprintf('Diag_%s_Ur', modelName));
dirP = fullfile(dirOut, sprintf('Diag_%s_P', modelName));
if ~isfolder(dirUr), mkdir(dirUr); end
if ~isfolder(dirP), mkdir(dirP); end
% Identification sets (each separately, as in the identification from t0
% on) and validation sets
filenameId = model.filenameId;
if ischar(filenameId), filenameId = {filenameId}; end
if ~isempty(filenameIdPlot)
    if ischar(filenameIdPlot), filenameIdPlot = {filenameIdPlot}; end
    filenameId = filenameIdPlot;
end
setNames = [filenameId(:); filenameVal(:)];
isId = [true(numel(filenameId),1); false(numel(filenameVal),1)];
R = table();
for idxV = 1:numel(setNames)
    if isId(idxV) % one identification set: D.id with all samples
        evalc('[~, D] = koopmanEstimation(''powerData'', setNames{idxV}, filenameVal{1}, 8, t0);');
        S = D.id;
    else
        evalc('[~, D] = koopmanEstimation(''powerData'', model.filenameId, setNames{idxV}, 8, t0);');
        S = D.val;
    end
    close(findobj('Type','figure','Name','PowerInputsOnly'));

    % States from the second sample on (initial state: measured Ur, history)
    urS = struct('M',S.lifted(1:2,1),'UrPrev',S.lifted(1:2,1),'k',1);
    X = zeros(model.n, size(S.lifted,2));
    for k = 1:size(S.lifted,2)
        [X(:,k), urS] = koopmanLiftedStates('featState', model, S.lifted(1:2,k), urS);
    end
    W = koopmanLiftedStates('powerFeatures', S.lifted(end-2:end,:), model.featureSet, model.tau);
    if endsWith(model.featureSet, ',lin') % linear model: feature states, the controls as inputs
        [F, W] = koopmanLiftedStates('featStatesLin', W);
        X = [X; F];
    end
    W = W(:,2:end); X = X(:,2:end);
    P = S.P(2:end)/1e6;
    Ur = S.lifted(1:2,2:end);
    [y, ~, xSim] = lsim(sys, W', [], X(:,1));
    y = y'; UrSim = xSim(:,1:2)';
    UrRef = []; yRef = [];
    if ~isempty(refModelFile) % reference model, simulated as the model
        [UrRef, yRef] = simModel(modelRef, S);
    end

    vafUr = [vaf(Ur(1,:), UrSim(1,:)), vaf(Ur(2,:), UrSim(2,:))];
    vafP = vaf(P, y);
    vafRef = [];
    if ~isempty(refModelFile)
        vafRef = [vaf(Ur(1,:), UrRef(1,:)), vaf(Ur(2,:), UrRef(2,:)), vaf(P, yRef)];
    end
    setName = strrep(setNames{idxV},'Vinf8dot0_sowfa_2turb_','');
    setType = {'Val.','Id.'};
    R = [R; table(string(setType{isId(idxV)+1}), string(setName), vafUr(1), vafUr(2), vafP, ...
        'VariableNames', {'Type','Set','VAF_Ur1','VAF_Ur2','VAF_P'})]; %#ok<AGROW>

    if iscell(setsPlot) && ~any(strcmp(setsPlot, setNames{idxV}))
        continue % VAF only
    end
    % Diagnostic plots as MainWfSimVinYawPowerOut (plotEDMDinputsEffWind):
    % inputs ct1, ct2, phi1 and WFSim (Real) vs. model (Est.)
    % Title: set, states and (added by plotEDMDinputsEffWind) the VAF
    strVal = sprintf('%s %s, $n_K$ = %d (%s)', setType{isId(idxV)+1}, strrep(setName,'_','-'), ...
        model.n, stateStr(model.stateSet));
    if ~isempty(refModelFile) % states of both models in the title key (estName, refName)
        strVal = sprintf('%s %s', setType{isId(idxV)+1}, strrep(setName,'_','-'));
    end
    vafRefUr = []; vafRefP = [];
    if ~isempty(vafRef), vafRefUr = vafRef(1:2); vafRefP = vafRef(3); end
    inp = S.lifted(end-2:end,2:end);
    t = t0 + (1:numel(P)); % time (s), samples from t0 + 1 on (dt = 1 s)
    fig = plotEDMDinputsEffWind(UrSim, inp, Ur, vafUr, dirUr, 1, model.n, strVal, 1, [], 'U', t, 'tex', UrRef, refName, estName, vafRefUr);
    set(fig,'Name',sprintf('Ur %s',setName));
    fig = plotEDMDinputsEffWind(y, inp, P, vafP, dirP, 0, model.n, strVal, 1, [], 'P', t, 'tex', yRef, refName, estName, vafRefP); % Vin = 0: 3rd input labelled gamma
    set(fig,'Name',sprintf('P %s',setName));
end
disp(R)
end

function [UrSim, y] = simModel(model, S)
% simModel simulates a power model from the measured initial state (from
% the second sample on, as above): estimates of Ur1, Ur2 and the farm power
urS = struct('M',S.lifted(1:2,1),'UrPrev',S.lifted(1:2,1),'k',1);
X = zeros(model.n, size(S.lifted,2));
for k = 1:size(S.lifted,2)
    [X(:,k), urS] = koopmanLiftedStates('featState', model, S.lifted(1:2,k), urS);
end
W = koopmanLiftedStates('powerFeatures', S.lifted(end-2:end,:), model.featureSet, model.tau);
if endsWith(model.featureSet, ',lin') % linear model: feature states, the controls as inputs
    [F, W] = koopmanLiftedStates('featStatesLin', W);
    X = [X; F];
end
[y, ~, xSim] = lsim(ss(model.A, model.B, model.C, model.D, 1), W(:,2:end)', [], X(:,2));
y = y'; UrSim = xSim(:,1:2)';
end

function str = stateStr(stateSet)
% Koopman states of the state set (testPowerStates) as LaTeX string
nL = sscanf(stateSet,'Ur n=%d');
names = {'$U_r$','$U_r^2$','$U_r^3$','$\Delta U_r$'};
str = strjoin(names(1:nL/2), ', ');
if endsWith(stateSet,'+d')
    str = [str, ', $U_r(k-1)$'];
end
end

%% ======================================================================
%% plotAICnoV1Sigma  (task 'aicNoV1Sigma')
function [fig, sel] = plotAICnoV1Sigma(dirFig, filenameVal, t0, ad, tSel, figSuffix)
% plotAICnoV1Sigma singular values of the Koopman effective wind speed model
% without V1 (identifyAICnoV1: inputs CT1, CT2, outputs Vr1, Vr2) during the
% online update, as Fig. 3 of the CDC22 paper: the model is updated in the
% open-loop validation with the wind speed ramp from 7 to 9 m/s
% (simAdaptiveKoopman, settings as validateAICnoV1) and compared with the
% models identified offline at constant 7 m/s (initial model) and 9 m/s.
% The singular values are those of the lifted linear model
% phi(k+1) = A phi(k) + B w(k), [Vr1; Vr2] = [I 0] phi(k), sample time 1 s.
%
% Inputs
% - dirFig: folder of AdaptiveAIC_OL_Sigma (.eps, .png), default
%   DataInOutWfSim/AdaptiveAIC
% - filenameVal, t0, ad: validation data and update settings (validateAICnoV1)
% - tSel: table with the times (Time_s) and colors (Color) of the shown
%   models, at most 4 (default aicNoV1ModelTimes: 240, 5100, 6000, 8000 s)
% - figSuffix: suffix of the file name, e.g. '_nu010_h5' (default ''); if
%   given, the title also shows nu and h
%
% Outputs
% - fig: figure handle
% - sel: table of the shown models (time, number of updates so far, line
%   color), input marks of validateAICnoV1

mainDir = fileparts(fileparts(mfilename('fullpath')));
% (path: see WFSim_addpaths) % koopmanAdaptStep
% (path: see WFSim_addpaths) % setFigureSize, printFigure
dirOut = fullfile(mainDir,'DataInOutWfSim','AdaptiveAIC');
if nargin < 1 || isempty(dirFig), dirFig = dirOut; end
if nargin < 2 || isempty(filenameVal), filenameVal = 'sowfa_2turb_ct_levels_val'; end
if nargin < 3 || isempty(t0), t0 = 240; end
if nargin < 4 || isempty(ad), ad = struct('h',5,'nu',0.15,'resetState',true,'trigChannels',1); end
if nargin < 5 || isempty(tSel), tSel = koopmanEstimation('aicNoV1ModelTimes'); end
if nargin < 6, figSuffix = ''; end

% Offline models at 7 m/s (initial model of the update) and 9 m/s
file7 = fullfile(dirOut,'KoopAICnoV1.mat');
if ~isfile(file7), koopmanEstimation('aicNoV1'); end
file9 = fullfile(dirOut,'KoopAICnoV1_Vinf9.mat');
if ~isfile(file9), koopmanEstimation('aicNoV1', 'Vinf9dot0_sowfa_2turb_ct_levels', [], 'KoopAICnoV1_Vinf9.mat'); end
tmp = load(file7,'model'); model7 = tmp.model;
tmp = load(file9,'model'); model9 = tmp.model;

% Validation data with the ramp (as validateAICnoV1)
f = dir(fullfile(mainDir,'DataT2OLWFSim','Vinf7dot0To9dot0_OL_Ct',['*',filenameVal,'.mat']));
S = load(fullfile(f(1).folder,f(1).name),'Ur1','Ur2','Ct1','Ct2');
k = t0:numel(S.Ur1);
Ur = [S.Ur1(k); S.Ur2(k)];
W = [S.Ct1(k); S.Ct2(k)];

% Shown models at the times tSel, colors as the marks in the validation
% plot (aicNoV1ModelTimes)
[~, upd] = koopmanEstimation('simAdaptive', model7, Ur, W, true, ad);
saveAt = tSel.Time_s' - k(1) + 1;
[~, ~, ~, ~, models] = koopmanEstimation('simAdaptive', model7, Ur, W, true, ad, saveAt);
cl = tSel.Color;
sel = table(k([models.k])', arrayfun(@(m) sum(upd(1:m.k)), models)', cl, ...
    'VariableNames', {'Time_s','NoUpdates','Color'});

% Singular values
omega = logspace(-5, log10(pi), 600);
sv = @(A,B) sigma(ss(A, B, eye(2,size(A,1)), zeros(2,size(B,2)), 1), omega);

fig = figure;
loglog(omega, sv(model9.A, model9.B)', 'k--', 'LineWidth', 2); hold on; % below the updated models, 2: visible next to them
ls = {'-','-.','--',':'};
lw = 1.7*ones(1,4);
names = cell(1,numel(models));
for idx = 1:numel(models)
    loglog(omega, sv(models(idx).A, models(idx).B)', 'Color', cl(idx,:), ...
        'LineStyle', ls{idx}, 'LineWidth', lw(idx));
    names{idx} = sprintf('%d s', k(models(idx).k));
end
grid on; axis tight;
xlabel('Frequency (rad/s)');
ylabel('Singular values (m/s)');
set(gca, 'TickLabelInterpreter', 'tex');
if ~isappdata(groot, 'ShowThesisTitle') || getappdata(groot, 'ShowThesisTitle') % title switch of runGenerateAllThesisFigures
    strTitle = {'Model C_{T1}, C_{T2} \rightarrow V_{r1}, V_{r2}, initial: offline ID 7 m/s', ...
        ['Model at t = ', colorKey(names, cl), '; offline ID ', colorKey({'9 m/s'}, [0 0 0])]};
    if ~isempty(figSuffix) % settings of the update
        adF = adaptSettings(ad);
        strTitle{end+1} = ['Update: ', sprintf('\\nu = %.2f m/s, h = %d', adF.nu, adF.h)];
    end
    title(strTitle, 'Interpreter', 'tex');
end
setFigureFont(fig, 12); % all text 12 pt, bold titles
setFigureSize(fig, [1 1]); % standard figure, independent of the screen
padYLim(fig); % y axis not completely tight
printFigure(fig, fullfile(dirFig, ['AdaptiveAIC_OL_Sigma', figSuffix]));
disp(sel)
end


%% ======================================================================
%% plotAICnoV1IdVal  (task 'aicNoV1IdVal')
function [R, files] = plotAICnoV1IdVal(t0, tEndVal)
% plotAICnoV1IdVal open-loop plots of the Koopman effective wind speed
% model without V1 (identifyAICnoV1, 7 m/s, non-adaptive) on its
% identification set (CT levels at V_inf = 7 m/s) and on the validation
% set up to the wind speed ramp (also 7 m/s, other CT sequence), plotted as
% the diagnostic plots (plotEDMDinputsEffWind: CT, V_inf and V1, Vr1, Vr2).
% The validation after the ramp is shown by validateAICnoV1.
%
% Inputs
% - t0: first sample (s), default 240 (wake established)
% - tEndVal: last sample (s) of the validation set, default 5000 (ramp)
%
% Outputs
% - R: VAF of Vr1, Vr2 on both sets
% - files: png/eps files without extension (DataInOutWfSim/AdaptiveAIC),
%   {identification, validation}

if nargin < 1 || isempty(t0), t0 = 240; end
if nargin < 2 || isempty(tEndVal), tEndVal = 5000; end
mainDir = fileparts(fileparts(mfilename('fullpath')));
% (path: see WFSim_addpaths) % simAdaptiveKoopman
koopmanDir = fileparts(mfilename('fullpath')); % vaf, plotEDMDinputsEffWind
% (path: see WFSim_addpaths)
dirOut = fullfile(mainDir,'DataInOutWfSim','AdaptiveAIC');
tmp = load(fullfile(dirOut,'KoopAICnoV1.mat'),'model');
model = tmp.model;

% {data file, last sample, title}
sets = {fullfile('Vinf7dot0_OL_Ct',model.filenameId), inf, 'Id. 7 m/s'
    fullfile('Vinf7dot0To9dot0_OL_Ct','Vinf9dot0_sowfa_2turb_ct_levels_val'), tEndVal, 'Val. 7 m/s'};
Wp.mesh = struct('Lx',1882.1,'Ly',800.0,'Nx',50,'Ny',25);
nTx = floor(0.4e3/Wp.mesh.Lx*Wp.mesh.Nx);
nTy = round(400/Wp.mesh.Ly*Wp.mesh.Ny);
R = table();
files = cell(1,size(sets,1));
for idx = 1:size(sets,1)
    S = load(fullfile(mainDir,'DataT2OLWFSim',sets{idx,1}),'Ur1','Ur2','Ct1','Ct2','u');
    k = t0:min(numel(S.Ur1), sets{idx,2});
    Ur = [S.Ur1(k); S.Ur2(k)];
    W = [S.Ct1(k); S.Ct2(k)];
    utemp = reshape(S.u, Wp.mesh.Nx, Wp.mesh.Ny, []); % V_inf and V1 for the plot only
    V1 = squeeze(utemp(nTx,nTy,k))';
    Vinf = squeeze(utemp(1,nTy,k))';

    UrEst = koopmanEstimation('simAdaptive', model, Ur, W, false); % non-adaptive
    vafUr = [vaf(Ur(1,:), UrEst(1,:)), vaf(Ur(2,:), UrEst(2,:))];
    R = [R; table(string(sets{idx,3}), vafUr(1), vafUr(2), ...
        'VariableNames', {'Set','VAF_Ur1','VAF_Ur2'})]; %#ok<AGROW>
    fig = plotEDMDinputsEffWind(UrEst, [W; V1], Ur, max(vafUr,0), dirOut, 1, model.n, ...
        sets{idx,3}, 1, Vinf, 'U', k, 'tex');
    set(fig,'Name',sets{idx,3});
    files{idx} = fullfile(dirOut,['Val',matlab.lang.makeValidName(sets{idx,3})]); % as plotEDMDinputsEffWind
end
disp(R)
end

%% ======================================================================
%% validateAICnoV1  (task 'aicNoV1Validation')
function R = validateAICnoV1(filenameVal, t0, ad, marks, figSuffix)
% validateAICnoV1 open-loop validation of the Koopman effective wind speed
% model without V1 (identifyAICnoV1, inputs CT1, CT2 only, no yaw) on the
% data set with the wind speed ramp from 7 to 9 m/s: non-adaptive and
% online updated model (simAdaptiveKoopman), plotted as in the thesis
% (plotEDMDinputsEffWindUpdate: CT, V_inf and V1, Ur1, Ur2, update signal).
% V1 is only shown, it is not used by the model.
%
% Inputs
% - filenameVal: open-loop data with the ramp (default: CT levels, ramp at 5000 s)
% - t0: first sample (s)
% - ad: settings of the online update (koopmanAdaptStep, simAdaptiveKoopman),
%   default h = 5, nu = 0.15 m/s, resetState = true, trigger on Vr1 only
%   (trigChannels = 1, as in the controller)
% - marks: times marked in color in the adaptive plot, table with Time_s
%   and Color (default aicNoV1ModelTimes: the models shown with their
%   singular values in plotAICnoV1Sigma), [] for none
% - figSuffix: suffix of the file name of the adaptive plot, e.g.
%   '_nu010_h5' (default ''); if given, the title also shows nu and h

mainDir = fileparts(fileparts(mfilename('fullpath')));
% (path: see WFSim_addpaths) % koopmanAdaptStep
koopmanDir = fileparts(mfilename('fullpath')); % vaf, plotEDMDinputsEffWindUpdate
% (path: see WFSim_addpaths)
if nargin < 1 || isempty(filenameVal)
    filenameVal = 'sowfa_2turb_ct_levels_val';
end
if nargin < 2 || isempty(t0), t0 = 240; end
if nargin < 3 || isempty(ad), ad = struct('h',5,'nu',0.15,'resetState',true,'trigChannels',1); end
if nargin < 4, marks = koopmanEstimation('aicNoV1ModelTimes'); end
if nargin < 5, figSuffix = ''; end
dirOut = fullfile(mainDir,'DataInOutWfSim','AdaptiveAIC');
modelFile = fullfile(dirOut,'KoopAICnoV1.mat');
if ~isfile(modelFile)
    koopmanEstimation('aicNoV1');
end
tmp = load(modelFile,'model');
model = tmp.model;

% Validation data (folder of WFSim_demo for the wind time series)
f = dir(fullfile(mainDir,'DataT2OLWFSim','Vinf7dot0To9dot0_OL_Ct',['*',filenameVal,'.mat']));
S = load(fullfile(f(1).folder,f(1).name),'Ur1','Ur2','Ct1','Ct2','u');
k = t0:numel(S.Ur1);
Ur = [S.Ur1(k); S.Ur2(k)];
W = [S.Ct1(k); S.Ct2(k)];

% V_inf and V1 (u in the grid cell in front of WT1) for the plot only
Wp.mesh = struct('Lx',1882.1,'Ly',800.0,'Nx',50,'Ny',25);
nTx = floor(0.4e3/Wp.mesh.Lx*Wp.mesh.Nx);
nTy = round(400/Wp.mesh.Ly*Wp.mesh.Ny);
utemp = reshape(S.u, Wp.mesh.Nx, Wp.mesh.Ny, []);
V1 = squeeze(utemp(nTx,nTy,k))';
Vinf = squeeze(utemp(1,nTy,k))';

R = table();
for adaptive = [false true]
    [UrEst, upd] = koopmanEstimation('simAdaptive', model, Ur, W, adaptive, ad);
    vafUr = [vaf(Ur(1,:), UrEst(1,:)), vaf(Ur(2,:), UrEst(2,:))];
    R = [R; table(adaptive, vafUr(1), vafUr(2), sum(upd), ...
        'VariableNames', {'Adaptive','VAF_Ur1','VAF_Ur2','NoUpdates'})]; %#ok<AGROW>
    if adaptive
        strVal = 'Val. adaptive';
        marksPlot = marks;
    else
        strVal = 'Val. non-adaptive';
        marksPlot = [];
    end
    figName = ['ValP', matlab.lang.makeValidName(strVal), repmat(figSuffix, 1, adaptive)];
    if adaptive && ~isempty(figSuffix) % settings of the update in the title
        adF = adaptSettings(ad);
        strVal = sprintf('%s ($\\nu$ = %.2f m/s, h = %d)', strVal, adF.nu, adF.h);
    end
    fig = plotEDMDinputsEffWindUpdate(UrEst, [W; V1], Ur, max(vafUr,0), dirOut, 1, model.n, ...
        strVal, 1 + upd, 1, Vinf, k, 'tex', marksPlot, figName); % time axis (s), TeX
    set(fig,'Name',strVal);
end
disp(R)
end

%% ======================================================================
%% testAICnoV1Vinf9  (task 'aicNoV1Test9')
function [R, fig] = testAICnoV1Vinf9(dirFig, t0, ad, figSuffix)
% testAICnoV1Vinf9 test of the Koopman effective wind speed models without
% V1 (inputs CT1, CT2, outputs Vr1, Vr2) on unseen data at constant
% V_inf = 9 m/s (validation CT levels, seed 4, not used for any model):
%   - offline ID 7 m/s: initial model of the online update
%   - offline ID 9 m/s: identified on the 9 m/s training data (seed of
%     sowfa_2turb_ct_levels)
%   - learned: 7 m/s model after the online update on the ramp from 7 to
%     9 m/s (validateAICnoV1, plotAICnoV1Sigma), frozen (an update on the
%     test data would not trigger: the Vr1 error stays below nu)
% Estimation as in the MPC (simAdaptiveKoopman, free run from the measured
% initial state).
%
% Inputs
% - dirFig: folder of AdaptiveAIC_Test9 (.eps, .png), default
%   DataInOutWfSim/AdaptiveAIC
% - t0: first sample (s), default 240 (wake established, as validateAICnoV1)
% - ad: settings of the online update (validateAICnoV1)
% - figSuffix: suffix of the file name, e.g. '_nu010_h5' (default ''); if
%   given, the name of the learned model also shows nu and h
%
% Outputs
% - R: table of the VAF and RMSE of Vr1, Vr2
% - fig: figure handle

mainDir = fileparts(fileparts(mfilename('fullpath')));
% (path: see WFSim_addpaths) % koopmanAdaptStep
% (path: see WFSim_addpaths) % setFigureSize, printFigure
koopmanDir = fileparts(mfilename('fullpath')); % vaf
% (path: see WFSim_addpaths)
dirOut = fullfile(mainDir,'DataInOutWfSim','AdaptiveAIC');
if nargin < 1 || isempty(dirFig), dirFig = dirOut; end
if nargin < 2 || isempty(t0), t0 = 240; end
if nargin < 3, ad = struct('h',5,'nu',0.15,'resetState',true,'trigChannels',1); end
if nargin < 4, figSuffix = ''; end

% Offline models
file7 = fullfile(dirOut,'KoopAICnoV1.mat');
if ~isfile(file7), koopmanEstimation('aicNoV1'); end
file9 = fullfile(dirOut,'KoopAICnoV1_Vinf9.mat');
if ~isfile(file9), koopmanEstimation('aicNoV1', 'Vinf9dot0_sowfa_2turb_ct_levels', [], 'KoopAICnoV1_Vinf9.mat'); end
tmp = load(file7,'model'); model7 = tmp.model;
tmp = load(file9,'model'); model9 = tmp.model;

% Learned model: online update on the ramp validation data (as validateAICnoV1)
f = dir(fullfile(mainDir,'DataT2OLWFSim','Vinf7dot0To9dot0_OL_Ct','*sowfa_2turb_ct_levels_val.mat'));
S = load(fullfile(f(1).folder,f(1).name),'Ur1','Ur2','Ct1','Ct2');
k = 240:numel(S.Ur1);
[~, ~, ~, modelL] = koopmanEstimation('simAdaptive', model7, [S.Ur1(k); S.Ur2(k)], [S.Ct1(k); S.Ct2(k)], true, ad);

% Test data at 9 m/s
S = load(fullfile(mainDir,'DataT2OLWFSim','Vinf9dot0_OL_Ct', ...
    'Vinf9dot0_sowfa_2turb_ct_levels_val.mat'),'Ur1','Ur2','Ct1','Ct2');
k = t0:numel(S.Ur1);
Ur = [S.Ur1(k); S.Ur2(k)];
W = [S.Ct1(k); S.Ct2(k)];

% The VAF ignores offsets (offline ID 7 m/s: about 1.5 m/s too low), hence
% also the RMSE
cases = {'Offline ID 7 m/s', model7
         'Offline ID 9 m/s', model9
         'Learned', modelL};
if ~isempty(figSuffix) % settings of the update
    adF = adaptSettings(ad);
    cases{3,1} = ['Learned (', sprintf('\\nu = %.2f m/s, h = %d', adF.nu, adF.h), ')'];
end
nC = size(cases,1);
UrEst = cell(nC,1);
R = table();
for idx = 1:nC
    UrEst{idx} = koopmanEstimation('simAdaptive', cases{idx,2}, Ur, W, false);
    rmse = sqrt(mean((Ur - UrEst{idx}).^2, 2));
    R = [R; table(string(cases{idx,1}), vaf(Ur(1,:), UrEst{idx}(1,:)), ...
        vaf(Ur(2,:), UrEst{idx}(2,:)), rmse(1), rmse(2), ...
        'VariableNames', {'Model','VAF_Ur1','VAF_Ur2','RMSE_Ur1','RMSE_Ur2'})]; %#ok<AGROW>
end

% Plot: measured and estimated Vr1, Vr2
fig = figure;
cl = [0 0 0; lines(nC)];
iPlot = 1:nC;
tl = tiledlayout(fig, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
axV = gobjects(2,1);
for idxV = 1:2
    axV(idxV) = nexttile(tl, idxV);
    plot(k, Ur(idxV,:), 'Color', cl(1,:), 'LineWidth', 1); hold on;
    for idx = iPlot
        plot(k, UrEst{idx}(idxV,:), 'Color', cl(idx+1,:), 'LineWidth', 1);
    end
    grid on; axis tight;
    ylabel(sprintf('Wind V_{r%d} (m/s)', idxV));
    if idxV == 2, xlabel('Time (s)'); end
    rmseStr = sprintf('%.2f/', R{iPlot, sprintf('RMSE_Ur%d', idxV)});
    rmseStr = ['RMSE = ', rmseStr(1:end-1), ' m/s'];
    if idxV == 1
        if ~isappdata(groot, 'ShowThesisTitle') || getappdata(groot, 'ShowThesisTitle') % title switch of runGenerateAllThesisFigures
            title({['Test data 9 m/s (unseen), V_{r1}: ', rmseStr], ...
                ['Measured, ', colorKey(cellstr(R.Model(iPlot)), cl(iPlot+1,:))]}, ...
                'Interpreter', 'tex');
        end
    else
        if ~isappdata(groot, 'ShowThesisTitle') || getappdata(groot, 'ShowThesisTitle') % title switch of runGenerateAllThesisFigures
            title(['V_{r2}: ', rmseStr], 'Interpreter', 'tex');
        end
    end
end
setFigureFont(fig, 12); % all text 12 pt, bold titles
setTimeLabelFontSize(fig); % time axis label 14 (runGenerateAllThesisFigures)
setFigureSize(fig, [1 1]); % standard figure, independent of the screen
padYLim(fig); % y axes not completely tight
tightenTiles(axV); % x tick labels only at the bottom, no y tick labels at the limits
printFigure(fig, fullfile(dirFig, ['AdaptiveAIC_Test9', figSuffix]));
disp(R)
end


%% ======================================================================
%% writeWRCModelTable  (task 'wrcModelTable')
function writeWRCModelTable(VAFall, resCL, featureSets, modelNames, maxEig, fileOut)
% writeWRCModelTable writes the comparison of the input-only Koopman power
% models of the WRC-AIC chapter as LaTeX table (booktabs): per model (row)
% the VAF of Ur1, Ur2 and the farm power on the diagnostic data sets, the
% largest eigenvalue magnitude of A and TE / AA of the three closed loops
%
% Inputs
% - VAFall: table of plotPowerStatesDiag with column Model (feature set)
% - resCL: closed loop results (columns Model, Case, TE_kW, AA from 500 s,
%   optional TE240_kW, AA240 from sample 240 on, in brackets)
% - featureSets: feature sets of the rows (cell, order of the rows)
% - modelNames: names of the rows (cell)
% - maxEig: max|eig(A)| per row
% - fileOut: file name of the table (.tex)

sets = unique(VAFall.Set(VAFall.Model == featureSets{1}), 'stable');
types = VAFall.Type(VAFall.Model == featureSets{1});
cases = unique(resCL.Case, 'stable');

fid = fopen(fileOut, 'w');
fprintf(fid, '%% Created by writeWRCModelTable (mainGeneratePlotThesis) on %s\n', char(datetime('today','Format','yyyy-MM-dd')));
fprintf(fid, '%s\n', '\begin{table}[h]');
fprintf(fid, '%s\n', '    \centering');
fprintf(fid, '%s\n', ['    \caption[Input-only Koopman power models of the WRC-AIC]{\label{tab:WRCModelComparison} ', ...
    'Input-only Koopman power models of the WRC-AIC: VAF of the effective wind speeds $V_{r1}$, $V_{r2}$ and ', ...
    'the wind farm power $P_\text{WF}$ (simulated from the measured initial state, VAF values below 0 are set to 0), ', ...
    'largest eigenvalue magnitude of $A$ and closed-loop tracking error TE and actuator activity AA ', ...
    'from 500\,s on (as in Figures~\ref{fig:CLbaseline} to~\ref{fig:CLMPC2}) and, in brackets, from sample 240 on.}']);
fprintf(fid, '%s\n', '    \resizebox{\textwidth}{!}{%');
fprintf(fid, '    \\begin{tabular}{@{}l*{%d}{c}@{}}\n', numel(sets) + 1 + numel(cases));
fprintf(fid, '%s\n', '        \toprule');
fprintf(fid, '         & \\multicolumn{%d}{c}{VAF $V_{r1}$ / $V_{r2}$ / $P_\\text{WF}$ (\\%%)} & & \\multicolumn{%d}{c}{TE (kW) / AA $\\cdot 10^{-3}$} \\\\\n', ...
    numel(sets), numel(cases));
fprintf(fid, '        \\cmidrule(lr){2-%d}\\cmidrule(lr){%d-%d}\n', numel(sets) + 1, numel(sets) + 3, numel(sets) + 2 + numel(cases));
head = 'Model';
for idx = 1:numel(sets)
    head = sprintf('%s & %s \\texttt{%s}', head, types(idx), strrep(sets(idx),'_','\_'));
end
head = [head, ' & $\max|\lambda(A)|$'];
for idx = 1:numel(cases)
    head = sprintf('%s & %s', head, cases(idx));
end
fprintf(fid, '        %s \\\\\n', head);
fprintf(fid, '%s\n', '        \midrule');
for idxR = 1:numel(featureSets)
    row = modelNames{idxR};
    for idx = 1:numel(sets)
        v = VAFall(VAFall.Model == featureSets{idxR} & VAFall.Set == sets(idx), :);
        row = sprintf('%s & %.1f / %.1f / %.1f', row, v.VAF_Ur1, v.VAF_Ur2, v.VAF_P);
    end
    row = sprintf('%s & %.4f', row, maxEig(idxR));
    for idx = 1:numel(cases)
        r = resCL(resCL.Model == featureSets{idxR} & resCL.Case == cases(idx), :);
        row = sprintf('%s & %.1f / %.1f', row, r.TE_kW, r.AA*1e3);
        if ismember('TE240_kW', r.Properties.VariableNames) && ~isnan(r.TE240_kW)
            row = sprintf('%s (%.1f / %.1f)', row, r.TE240_kW, r.AA240*1e3);
        end
    end
    fprintf(fid, '        %s \\\\\n', row);
end
fprintf(fid, '%s\n', '        \bottomrule');
fprintf(fid, '%s\n', '    \end{tabular}}');
fprintf(fid, '%s\n', '\end{table}');
fclose(fid);
fprintf('Table written: %s\n', fileOut);
end
