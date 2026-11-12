function testProfileViewLargerPlot(testOn)
close all;

% --- SWITCH FOR PAPER ---
useLatexForPaper = false; 
% ------------------------

if useLatexForPaper
    strInterpreter = 'latex';
    % Set colors to Gray for paper version
    colInit = [0.5, 0.5, 0.5]; 
    colFarm = [0.5, 0.5, 0.5]; 
    colCtrl = [0.5, 0.5, 0.5];
else
    strInterpreter = 'tex'; % or 'tex'
    % Original colors (Red/Blue)
    colInit = [0, 0, 0];
    colFarm = [0, 0, 1];
    colCtrl = [1, 0, 0];
end

if ~nargin
    testOn = 0;
end

if testOn
    profile on -history
    figure(1);
    plot(magic(4));
    p = profile('info');
    profile off;
    save p p;
else
    %load('profileInfoSMPC.mat','profileInfo');
    load(fullfile('OutDir','profileInfoCtrl1.mat'),'profileInfo');
    %load('profileInfo.mat', 'profileInfo');
    p = profileInfo;
end

xPosIndent = 0.05;
yPosIndent = 0.0325; %0.0475;
nrFig = 1;

strPwd = pwd;

% Passed strInterpreter to the subfunction
testProfilePlotSubFun(p, xPosIndent, yPosIndent, nrFig, strPwd, testOn, strInterpreter);

close(2);

figDir = 'figDir';

if ~isdir(figDir) %#ok<ISDIR>
    mkdir(figDir)
end

figure1 = gcf;
moveCorner = 0.25;
moveCorner2 = 0.2;


%Create textbox (around initilization)
shiftUpWF = 0.098;
leftCorner = 0.11;
scaleWF = 0.7;

annotation(figure1,'textbox',...
   [leftCorner 0.63+shiftUpWF 0.86-moveCorner 0.32*scaleWF],...
   'Color', colInit,...
   'FitBoxToText','off', ...
   'EdgeColor', colInit, ... % Set edge color
   'Interpreter', strInterpreter); 

annotation(figure1,'textbox',...
    [0.65-moveCorner2 0.64+shiftUpWF 0.07 0.06*scaleWF],...
    'String','Initialization',...
    'FontSize', 15,...
    'Color', colInit, ... % Text color
    'EdgeColor',[1,1,1],...
    'FitBoxToText','off', ...
    'Interpreter', strInterpreter);


% Create textbox Wind farm
shiftUpWF = 0.288;
scaleWF = 0.702;
annotation(figure1,'textbox',...
    [leftCorner 0.07+shiftUpWF  0.86-moveCorner 0.525*scaleWF],...
    'FitBoxToText','off',...
    'EdgeColor', colFarm, ...
    'Interpreter', strInterpreter); 

annotation(figure1,'textbox',...
    [0.65-moveCorner2 0.08+shiftUpWF  0.5 0.06*scaleWF],...
    'Color', colFarm,...
    'String','Farm Model',...
    'FontSize',15,...
    'FitBoxToText','off',...
    'EdgeColor',[1 1 1], ...
    'Interpreter', strInterpreter);

% Create textbox controller 
shiftUpCtl = -0.481;
scaleCtl = 0.675;
annotation(figure1,'textbox',...
    [leftCorner 0.595 + shiftUpCtl 0.86-moveCorner 0.36*scaleCtl],...
    'FitBoxToText','off',...
    'EdgeColor', colCtrl, ...
    'Interpreter', strInterpreter); 

annotation(figure1,'textbox',[0.65-moveCorner2 0.58 + shiftUpCtl 0.07 0.1*scaleCtl],...
    'String','Controller',...
    'FontSize',15,...
    'Color', colCtrl,...
    'FitBoxToText','off',...
    'EdgeColor',[1 1 1], ...
    'Interpreter', strInterpreter);

% Create textbox
strFig ='FlowchartInitCtrlPlant_1';
print(fullfile(figDir,strFig), '-dpng');
print(fullfile(figDir,strFig), '-depsc');

end

function testProfilePlotSubFun(p, xPosIndent, yPosIndent, nrFig, strPwd, testOn, strInterpreter)

xPosIndentMult = -1;
yPosIndentMult = -1;
cntAxAy = 1;

if strcmp(strInterpreter,'latex')
    strCellAxAy = {' (for $A_x$)',' (for $A_y$)'};
else
    strCellAxAy = {' (for A_x)',' (for A_y)'};
end

newFigListFun;
listFunNotAdded = {'vec'};
curFunName = '';

mainDir = fileparts(fileparts(fileparts(fileparts(pwd))));
YalmipStr = 'YALMIP-master'; %fullfile(mainDir,'YALMIP-master');

yalmipFunStr = fullfile(strPwd,'YALMIP-master','extras','yalmip.m');
yalmipCnt = 0; %only display the first call to yalmip

optimizeFunStr = fullfile(strPwd,'YALMIP-master','extras','optimize.m');
sdpsettingsFunStr = fullfile(strPwd,'YALMIP-master','extras','sdpsettings.m');

for n = 1:size(p.FunctionHistory,2)
    
    aFunction = p.FunctionTable(p.FunctionHistory(2,n));
    isScriptOrFunction = strcmp(aFunction.Type,'M-function') || ...
        strcmp(aFunction.Type,'M-script'); %only add m files here

    if testOn
        isInDir = 1;
    else

        strOrig = 'C:\Projects\201221_FastFarm';
        strNew = 'D:\ditt_aj\Projekte\201221_FastFarmWFSimVeeresh';
        aFunction.FileName = strrep(aFunction.FileName,strOrig, strNew);

        isInDir = ~isempty(strfind(aFunction.FileName,strPwd)); %check WFSim function
    end
    
    yalmipFun = ~isempty(strfind(aFunction.FileName,yalmipFunStr));
    optimizeFun = ~isempty(strfind(aFunction.FileName,optimizeFunStr));
    
    isNotInYALMIP = isempty(strfind(aFunction.FileName,YalmipStr)) || ...   %check WFSim function
        (yalmipFun && ~yalmipCnt) || optimizeFun;
    
    if yalmipFun
        yalmipCnt = 1;
    end
    
    addFun = ~contains(aFunction.FunctionName,listFunNotAdded);
    
    % Note: Latex handles \_ as underscore, so this replacement logic works for both
    tmp =  strrep(strrep( aFunction.FunctionName, '\','/'), '_','\_');
    
    if  isInDir && isScriptOrFunction && addFun &&  isNotInYALMIP %&&  ~strcmp(tmp, curFunName)
             
        isNewFun = ~strcmp(tmp, curFunName); %Don't add a function name if it was just called
        if isNewFun
            curFunName = tmp;
        end
        
        isEnterFun = ~p.FunctionHistory(1,n);
        if isEnterFun %&& isNewFun
            str = 'entering function: ';
            xPosIndentMult = xPosIndentMult + 1;
            yPosIndentMult = yPosIndentMult + 1;
            curXPos = 0 + xPosIndent*xPosIndentMult;
            curYPos = 1 - yPosIndent*yPosIndentMult;
            if curYPos < 0
                yPosIndentMult = 0;
                curYPos = 1;
                newFigListFun;
            end
            
            % curFunName = strrep(curFunName,'9turb\_apc\_alm\_turbl','2turb\_adm\_uniform');
            % curFunName = strrep(curFunName,'sowfa\_9turb\_apc\_alm\_turbl','sowfa\_2turb\_apc\_alm\_turbl'); %controlSet_sowfa_2turb_apc_alm_turbl(Wp)
            % 

            if strcmp(curFunName,'Actuator') || strcmp(curFunName,'BoundaryConditions') %|| contains(curFunName,'layout')
                text(curXPos,curYPos,curFunName,'Color','k', 'Interpreter', strInterpreter);
            elseif strcmp(curFunName,'MakingSparseMatrix') && cntAxAy <= 2
                text(curXPos,curYPos,[curFunName strCellAxAy{cntAxAy}], 'Interpreter', strInterpreter);
                cntAxAy = cntAxAy + 1;
            else
                text(curXPos,curYPos,curFunName, 'Interpreter', strInterpreter);
            end
            
        else
            str = ' exiting function: ';
            xPosIndentMult = max(xPosIndent,xPosIndentMult - 1);
        end
        %disp([str curFunName]);
    end
    %end
end

    function newFigListFun
        figure(nrFig); nrFig = nrFig +1;
        pos0 = get(0,'defaultFigurePosition');
        set(nrFig-1,'position',[pos0(1), pos0(2)-pos0(4)*0.4, pos0(3),pos0(4)*1.4]);
        set(gca,'XColor','white');
        set(gca,'YColor','white');
        set(gcf,'Color','white');
    end
end