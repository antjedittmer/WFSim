function testProfileView(ctrlNum,OutDir)
close all;

if ~nargin
    ctrlNum = 3;
end

if nargin <2
    OutDir = 'OutDir';
end

if ctrlNum == -1 % test based on plot magic
    profile on -history
    figure(1);
    plot(magic(4));
    p = profile('info');
    profile off;
    save p p;
    nrFig = 2;
else
    %load('profileInfoSMPC.mat','profileInfo');
    matStr = sprintf('profileInfoCtrl%d.mat',ctrlNum);
    fileStr = fullfile(OutDir,matStr);
    load(fileStr,'profileInfo');
    p = profileInfo;
    nrFig = 1;


end

xPosIndent = 0.05;
yPosIndent = 0.0475;
if ctrlNum == 3, yPosIndent = 0.045; end

strPwd = pwd;

testProfilePlotSubFun(p, xPosIndent, yPosIndent, nrFig,strPwd,ctrlNum);

end

function testProfilePlotSubFun(p, xPosIndent, yPosIndent, nrFig,strPwd,ctrlNum)

xPosIndentMult = -1;
yPosIndentMult = -1;

newFigListFun;
listFunNotAdded = {'vec'};
curFunName = '';

YalmipStr = fullfile(strPwd,'YALMIP-master');

yalmipFunStr = fullfile(strPwd,'YALMIP-master','extras','yalmip.m');
yalmipCnt = 0; % only display the first call to yalmip
optimizeFunStr = fullfile(strPwd,'YALMIP-master','extras','optimize.m');
%sdpsettingsFunStr = fullfile(strPwd,'YALMIP-master','extras','sdpsettings.m');

for n = 1:size(p.FunctionHistory,2)

    aFunction = p.FunctionTable(p.FunctionHistory(2,n));
    isScriptOrFunction = strcmp(aFunction.Type,'M-function') || ...
        strcmp(aFunction.Type,'M-script'); %only add m files here

    isInDir = ~isempty(strfind(aFunction.FileName,strPwd)); %check WFSim function

    yalmipFun = ~isempty(strfind(aFunction.FileName,yalmipFunStr));
    optimizeFun = ~isempty(strfind(aFunction.FileName,optimizeFunStr));

    isNotInYALMIP = isempty(strfind(aFunction.FileName,YalmipStr)) || ...   %check WFSim function
        (yalmipFun && ~yalmipCnt) || optimizeFun;

    if yalmipFun
        yalmipCnt = 1;
    end

    addFun = ~contains(aFunction.FunctionName,listFunNotAdded);

    tmp =  strrep(strrep( aFunction.FunctionName, '\','/'), '_','\_');

    if  (isInDir || ctrlNum == -1) && isScriptOrFunction && addFun &&  isNotInYALMIP %&&  ~strcmp(tmp, curFunName)


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
                if nrFig > 4
                    return
                end
                newFigListFun;
            end
            text(curXPos,curYPos,curFunName);

        else
            str = ' exiting function: ';
            xPosIndentMult = max(xPosIndent,xPosIndentMult - 1);
        end
        if ctrlNum < 0 % only display in Command Window for test cases
            disp([str curFunName]);
        end
    end
    %end


end

    function newFigListFun

        figure(nrFig); nrFig = nrFig +1;
        %set(gcf,'position',[360.0000  117.6667  560.0000  500.3333]);
        set(gca,'XColor','white');
        set(gca,'YColor','white');
        set(gcf,'Color','white');
    end
end