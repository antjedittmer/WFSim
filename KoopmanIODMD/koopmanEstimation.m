function varargout = koopmanEstimation(task, varargin)
%KOOPMANESTIMATION Identification of the Koopman models from open-loop WFSim data
%   varargout = koopmanEstimation(task, ...) identifies the Koopman models
%   of the wind farm (EDMD / least squares on the lifted states of
%   koopmanLiftedStates) from the open-loop data of WFSim_demo
%   (DataT2OLWFSim) and saves them in DataInOutWfSim.
%
%   Tasks (original function names in brackets):
%   - 'modelVin' (MainWfSimVin): effective wind speed model with the inputs
%     CT1, CT2 and the measured wind speed V1 (eDMDresults_UasOutput)
%   - 'modelVinYaw' (MainWfSimVinYaw): as 'modelVin' with the yaw angle as
%     input (IFAC23 paper)
%   - 'modelVinYawPower' (MainWfSimVinYawPowerOut): effective wind speed and
%     farm power model with yaw and V1 (eDMDresults_UasOutput_MIMO_PT)
%   - 'powerData' (testPowerInputsOnly): lifted identification and
%     validation data of the input-only power models
%   - 'powerStates' (testPowerStates): input-only Koopman power models for
%     state and feature sets (KoopFeatModel<k>.mat)
%   - 'linearPowerModel' (buildLinearPowerModel): linear power model with 21
%     states from the control-affine model
%   - 'aicNoV1' (identifyAICnoV1): effective wind speed model of the
%     adaptive AIC with the inputs CT1, CT2 only (KoopAICnoV1.mat)
%   - 'simAdaptive' (simAdaptiveKoopman): open-loop simulation of the
%     online updated model (koopmanAdaptStep)
%   - 'aicNoV1ModelTimes' (aicNoV1ModelTimes): times of the updated models
%     shown in the thesis figures
switch task
    case 'modelVin',          [varargout{1:nargout}] = MainWfSimVin(varargin{:});
    case 'modelVinYaw',       [varargout{1:nargout}] = MainWfSimVinYaw(varargin{:});
    case 'modelVinYawPower',  [varargout{1:nargout}] = MainWfSimVinYawPowerOut(varargin{:});
    case 'powerData',         [varargout{1:nargout}] = testPowerInputsOnly(varargin{:});
    case 'powerStates',       [varargout{1:nargout}] = testPowerStates(varargin{:});
    case 'linearPowerModel',  [varargout{1:nargout}] = buildLinearPowerModel(varargin{:});
    case 'aicNoV1',           [varargout{1:nargout}] = identifyAICnoV1(varargin{:});
    case 'simAdaptive',       [varargout{1:nargout}] = simAdaptiveKoopman(varargin{:});
    case 'aicNoV1ModelTimes', [varargout{1:nargout}] = aicNoV1ModelTimes(varargin{:});
    otherwise, error('koopmanEstimation:task', 'Unknown task %s', task);
end
end

%% ======================================================================
%% MainWfSimVin  (task 'modelVin')
function MainWfSimVin(yawmode,filenameId,filenameVal,noStates,PolyVec,useVal,percentTrain,NTMstr,t0)
% Code based on
% Model Predictive Control For Wake Steering: a Koopman Dynamic Mode
% Decomposition Approach
% Master Thesis Dissertation
% Author: Nassir Rodrigues Cassamo
% Supervisors: Professor Jan-Willem Van Wingerden and Professor Jo?o Sousa
% First official commit (1.0): December 2020

%% RELEVAT INFORMATION
%This script has the following goals:
% (1) Assess Simulation data, both qualitatively (animations) and
% quantitatively (graphics)
% (2) Derives a low dimensional model using Dynamic Mode Decomposition
% (variations included to take into account input-output data and other
% known states - deterministic states)
% (3) Validates the models with a set of validaiton data

%% (0) INITIALISE
% define type of simulation
%clc; clear; close all;
if ~nargin
    % directory for identification and validation data, yaw and pitch
    yawmode = 0; %0 for ct (pitch/torque) control, 1 for additional yaw control
    

    filenameId = 'Vinf8dot0_sowfa_2turb_alm_turbl_AllComb.mat'; NTMstr = '';
    filenameVal = 'Vinf8dot0_sowfa_2turb_alm_turbl_AllComb.mat'; NTMstr = '';

    
    % User input
    noStates = 6; %number of Koopman states [,12,12,14,18,24]
    useVal = 1; % use extra data for validation
    percentTrain = .6;% 60% data for validation
    NTMstr = '';
    PolyVec = zeros(size(noStates)); % 1: use only polynomial, 0: otherwise
    idxPolyOn = find(noStates == 12, 1, 'first');
    PolyVec(idxPolyOn) = 1;
    t0 = 1;   
end


% Code based on
% Model Predictive Control For Wake Steering: a Koopman Dynamic Mode
% Decomposition Approach
% Master Thesis Dissertation
% Author: Nassir Rodrigues Cassamo
% Supervisors: Professor Jan-Willem Van Wingerden and Professor Jo?o Sousa
% First official commit (1.0): December 2020

%% RELEVAT INFORMATION
%This script has the following goals:
% (1) Assess Simulation data, both qualitatively (animations) and
% quantitatively (graphics)
% (2) Derives a low dimensional model using Dynamic Mode Decomposition
% (variations included to take into account input-output data and other
% known states - deterministic states)
% (3) Validates the models with a set of validaiton data


method = -2; %0: DMD ; 1:DMDc; 2:IODMD; 3:EXTIODMD -1: EDMD
koopman = 1; %to add deterministic states to flow field data

% Simulation characterisitc (resampling)
dt = 1; %time sampling (s)

KoopmanDir = fileparts(mfilename('fullpath')); % Get Koopman directory
pathFunctions = genpath(fullfile(KoopmanDir,'Functions'));

addpath(fullfile(KoopmanDir,'2.DYNAMIC_MODE_DECOMPOSITION'),...
    pathFunctions);
codedir = mfilename('fullpath');
parentdir = fileparts(fileparts(codedir));

% DEFINE MAIN output data DIRECTORY TO STORE ALL RESULTS
DataOut = fullfile(parentdir,'DataInOutWfSim');
if ~isfolder(DataOut)
    mkdir(DataOut)
end

DataIn = fullfile(parentdir,'DataT2OLWFSim'); %,...
sepStr = strfind(filenameId,'_');
if ~isfolder(DataIn)
    warning('DataIn directory does not exist')
    return
end

% load identification data
dirName = fullfile(DataIn, [filenameId(1:sepStr-1),'_OL_Ct'], filenameId);
tmpId = load(dirName,'Ct1', 'Ct2', 'FT1', 'FT2', 'PT1', 'PT2', ...
    'Ur1', 'Ur2', 'p', 'phi1', 'phi2', 'u', 'v');

if useVal == 1 % load different data set for validation
    sepStr = strfind(filenameVal,'_');
    filenameVal2 = strrep(filenameVal,'To','');
    dirNameVal = fullfile(DataIn, [filenameVal(1:sepStr-1),'_OL_Ct'], filenameVal2);
    tmpVal = load(dirNameVal,'Ct1', 'Ct2', 'FT1', 'FT2', 'PT1', 'PT2', ...
        'Ur1', 'Ur2', 'p', 'phi1', 'phi2', 'u', 'v');
end


%% (1) ASSESS DATA
% easy solution to augment u flow field data matricx with other flow
% field data
fieldCell = fieldnames(tmpId);

useVal = exist('tmpVal','var') && useVal;

for idx = 1:length(fieldCell)
    if useVal
        tmp = [tmpId.(fieldCell{idx}), tmpVal.(fieldCell{idx})]; %#ok<NASGU>
    else
        tmp = [tmpId.(fieldCell{idx})]; %#ok<NASGU>
    end
    eval([fieldCell{idx},'= tmp;']);
end

if useVal % get time
    tident = length(tmpId.Ct1);
else
    tident = floor(length(Ur1)*percentTrain);
end
tval = tident+t0;

% Get index of 
Wp.turbine = struct(... 
    'Crx',[0.4000    1.0321]*1e3,... % X-coordinates of turbines (m)
    'Cry',[400.0000  398.5747]... % Y-coordinates of turbines (m)
    );
Wp.mesh = struct(...
    'Lx',1882.1,... % Domain length in x-direction
    'Ly',800.0,... % Domain length in y-direction
    'Nx',50,... % Number of cells in x-direction
    'Ny',25 ... % Number of cells in y-direction
    );

nTx = round(Wp.turbine.Crx(1)/Wp.mesh.Lx * Wp.mesh.Nx)-1;
nTy = round(Wp.turbine.Cry(1)/Wp.mesh.Ly * Wp.mesh.Ny);

kk = length(u)/Wp.mesh.Ny;

utemp = reshape(u,Wp.mesh.Nx,Wp.mesh.Ny,kk);
utempT = squeeze(utemp(nTx,nTy,:));
vtemp = reshape(v,Wp.mesh.Nx,Wp.mesh.Ny,kk);
vtempT = squeeze(vtemp(nTx,nTy,:));

QQ_u1 = utempT(t0:tident)'; %#ok<*USENS>
QQ_v1 = vtempT(t0:tident)';
QQ_p1 = p(:,t0:tident);

valid.QQ_u1 = utempT(tval:end)';
valid.QQ_v1 = vtempT(tval:end)';
valid.QQ_p1 = u(:,tval:end);

VinfAll = squeeze(utemp(1,1,:));
Vinf = VinfAll(tval:end)';

%% (2) DYNAMIC MODE DECOMPOSITION
states_u = QQ_u1;%fluid flow as states, identification data set
statesvalid_u = valid.QQ_u1; %fluid flow as states, validation data set for comparison

if yawmode == 0
    Inputs = [Ct1(t0:tident);Ct2(t0:tident)];
    Inputs_val = [Ct1(tval:end);Ct2(tval:end)];
else
    Inputs = [Ct1(t0:tident);Ct2(t0:tident); phi1(t0:tident)];
    Inputs_val = [Ct1(tval:end);Ct2(tval:end); phi1(tval:end)];
end

Outputs = [PT1(t0:tident);PT2(t0:tident);FT1(t0:tident);FT2(t0:tident)];
Outputs_val = [PT1(tval:end);PT2(tval:end);FT1(tval:end);FT2(tval:end)];

% ToDo check

strVAF = cell(length(noStates),1);
for idx = 1: length(noStates)
    
    n = noStates(idx);
    poly = PolyVec(idx);
    
    % subdirectories in dataOut
    if poly == 1 && yawmode == 0
        dirdmdName = 'eDMDresults_UasOutput_poly';
    elseif poly == 0 && yawmode == 0
        dirdmdName = 'eDMDresults_UasOutput';
    elseif poly == 1
        dirdmdName = 'eDMDresults_UasOutput_MIMO_poly';
    else
        dirdmdName = 'eDMDresults_UasOutput_MIMO';
    end
    
    dirdmd = fullfile(DataOut,dirdmdName);
    if ~ isfolder(dirdmd)
        mkdir(dirdmd);
    end
    
    Deterministic = [Ur1(t0:tident); Ur2(t0:tident)];
    Deterministic_val = [Ur1(tval:end); Ur2(tval:end)];
    %include non linear observables - Koopman extensions to better recover non linear dynamics
    if koopman == 1
        
        structPrev.Ur1_prev1 = Deterministic(1,1);
        structPrev.Ur2_prev1 = Deterministic(2,1);
        structPrev.dUr1_prev1 = 0;
        structPrev.dUr2_prev1 = 0;
        structPrev.M1(1) = structPrev.Ur1_prev1; % Moving mean
        structPrev.M2(1) = structPrev.Ur2_prev1;
        structPrev.k = 1;
        
        statesUr = nan(n+2,size(Deterministic,2));
        for idxK = 1: length(Deterministic)
            [statesUr(:,idxK),stateNameUr,structPrev] = koopmanLiftedStates('urStates', [Deterministic(:,idxK);Inputs(:,idxK)],poly,n,structPrev); %koopmanstateextensionRT(Deterministic(:,idxK),poly,n,structPrev);
        end
        
        windAtTurbine = [QQ_u1; QQ_v1];
        statesUV = koopmanLiftedStates('windStates', windAtTurbine,poly,n);
        
        stateNameUV = regexprep(stateNameUr,{'Ur1','Ur2','M'},{'u','v','Muv'});
        states = [statesUr]; %;statesUV];
        stateName = [stateNameUr]; %,';',stateNameUV];
        
        statesvalidUr = nan(n+2,size(Deterministic_val,2));
        structPrev.Ur1_prev1 = Deterministic_val(1,1);
        structPrev.Ur2_prev1 = Deterministic_val(2,1);
        structPrev.dUr1_prev1 = 0;
        structPrev.dUr2_prev1 = 0;
        structPrev.M1(1) = Deterministic_val(1,1); % Moving mean
        structPrev.M2(1) = Deterministic_val(2,1);
        structPrev.k = 1;
        for idxK = 1: length(Deterministic_val)
            [statesvalidUr(:,idxK),st1,structPrev] = koopmanLiftedStates('urStates', [Deterministic_val(:,idxK);Inputs_val(:,idxK)],poly,n,structPrev); %koopmanstateextensionRT(Deterministic(:,idxK),poly,n,structPrev);
        end
        windAtTurbineVal = [valid.QQ_u1; valid.QQ_v1];
        statesUV_val = koopmanLiftedStates('windStates', windAtTurbineVal,poly,n);
        statesvalid = [statesvalidUr];%;statesUV_val];
    else
        states = states_u;
        statesvalid = statesvalid_u;
        temp = sprintf('u%d;',1: size(states,1));
        stateName = temp(1:end-1);
    end
    
    Vinf1 = utemp(1,1,t0);
    Vinf2 = utemp(1,1,tident);  Vstr = '';
    if abs(Vinf1 -Vinf2) >0.01 && isempty(NTMstr)
        Vstr = sprintf('To%2.1f',Vinf2);
    end
    folderName = strrep(sprintf('Vin_Vinf%2.1f%s%s_states%02d',Vinf1,Vstr,NTMstr,n),'.','dot');%'Vinf%d_diffComb_states%d'
    fileName = sprintf('stateName_K%02d_P%d.mat',n,poly); %stateName);
    
    dirFig = fullfile(dirdmd,folderName); % this depends on vinf
    if ~exist(dirFig,'dir')
        mkdir(dirFig);
    end
    
    %     cellStName = regexp(stateName,';','split'); cellStName(end+1) ={'U2p'};
    %     aTable = array2table([states(:,1:end-1)',states(1,2:end)'],...
    %         'VariableNames',cellStName);
    %     aMatrix = [states(:,1:end-1)',Inputs(:,1:end-1)',states(2,2:end)'];
    
    % Sys id with ct1,ct2 as inputs
%     [sys_red0,FITje0,Xd0,Xd_p0,xsim0,FITje_val0,xo_0,Koop0,...
%         ysim,Inputs,Deterministic,ysim_val,Inputs_val,Deterministic_val]= eDMD_RT_Vin(states,stateName,statesvalid,poly,n);
    
    %     koopmanPlots('effWind', ysim,Inputs,Deterministic, FITje0,dirFig,Vin,n,strVal);
    %     koopmanPlots('effWind', ysim_val,Inputs_val,Deterministic_val,FITje_val0,dirFig,Vin,n);
    %
    % Sys id with ct1,ct2,V1 as inputs
    strVal = 'Id.'; Vin = 1;
    statesVinf = [states;statesUV(1,:)];
    statesvalidVinf = [statesvalid; statesUV_val(1,:)];
    stateNameVinf = [stateName,';V1'];
    [sys_red,FITje,Xd,Xd_p,xsim,FITje_val,xo,Koop,...
        ysim,Inputs,Deterministic,ysim_val,Inputs_val,Deterministic_val] = eDMD_RT_Vin(statesVinf,stateNameVinf,statesvalidVinf,poly,n);
    koopmanPlots('effWind', ysim,Inputs,Deterministic, FITje,dirFig,Vin,n,strVal);
    set(gcf,'name',strrep([strVal,fileName],'.stateName', ['Vin',num2str(Vin)]));
    koopmanPlots('effWind', ysim_val,Inputs_val,Deterministic_val,FITje_val,dirFig,Vin,n);
    set(gcf,'name',strrep([strVal,fileName],'.stateName', ['Vin',num2str(Vin)]));
%     
    save(fullfile(dirFig,fileName),'sys_red','FITje','FITje_val','dirName','xo','Koop');
    
    strVAF{idx} = sprintf('%s \t\t %d\t\t\t %d\t\t\t %d \t\t%d\n',num2str(n,'%02d'),...
        round(FITje(1)), round(FITje_val(1)),round(FITje(2)), round(FITje_val(2)));
%     strVAF0{idx} = sprintf('%s \t\t %d\t\t\t %d\t\t\t %d \t\t%d\n',num2str(n,'%02d'),...
%         round(FITje0(1)), round(FITje_val0(1)),round(FITje0(2)), round(FITje_val0(2)));
    
end
end

%% ======================================================================
%% MainWfSimVinYaw  (task 'modelVinYaw')
function MainWfSimVinYaw(yawmode,filenameId,filenameVal,noStates,PolyVec,useVal,percentTrain,NTMstr,t0)% Code based on

% Code based on
% Model Predictive Control For Wake Steering: a Koopman Dynamic Mode
% Decomposition Approach
% Master Thesis Dissertation
% Author: Nassir Rodrigues Cassamo
% Supervisors: Professor Jan-Willem Van Wingerden and Professor Jo?o Sousa
% First official commit (1.0): December 2020

%% RELEVAT INFORMATION
%This script has the following goals:
% (1) Assess Simulation data, both qualitatively (animations) and
% quantitatively (graphics)
% (2) Derives a low dimensional model using Dynamic Mode Decomposition
% (variations included to take into account input-output data and other
% known states - deterministic states)
% (3) Validates the models with a set of validaiton data

%% (0) INITIALISE
% define type of simulation
%% (0) INITIALISE
% define type of simulation
if ~nargin
    
    % directory for identification and validation data, yaw and pitch
    yawmode = 1; %0 for ct (pitch/torque) control, 1 for additional yaw control
    % User input
    noStates = 6; %[6,12,12,14,16,18,24];
    useVal = 1; % use extra data for validation
    percentTrain = .9;
    t0 = 340;
    
    NTMstr = '';
    filenameId = 'Vinf8dot0_sowfa_2turb_yaw_noise_step'; %Vinf8dot0_sowfa_2turb_yaw_noise2';
    filenameVal = 'Vinf8dot0_sowfa_2turb_yaw_steps_Ct_comb';
       
    PolyVec = zeros(size(noStates)); % 1: use only polynomial, 0: otherwise
    % idxPolyOn = find(noStates == 12, 1, 'first');
    % PolyVec(idxPolyOn) = 1;
end

method = -2; %0: DMD ; 1:DMDc; 2:IODMD; 3:EXTIODMD -1: EDMD
koopman = 1; %to add deterministic states to flow field data

% Simulation characterisitc (resampling)
dt = 1; %time sampling (s)

KoopmanDir = fileparts(mfilename('fullpath')); % Get Koopman directory
pathFunctions = genpath(fullfile(KoopmanDir,'Functions'));

addpath(fullfile(KoopmanDir,'2.DYNAMIC_MODE_DECOMPOSITION'),...
    pathFunctions);
codedir = mfilename('fullpath');
parentdir = fileparts(fileparts(codedir));

% DEFINE MAIN output data DIRECTORY TO STORE ALL RESULTS
DataOut = fullfile(parentdir,'DataInOutWfSim');
if ~isfolder(DataOut)
    mkdir(DataOut)
end

if isempty(NTMstr)
    DataIn = fullfile(parentdir,'DataT2OLWFSim'); %,...
    sepStr = strfind(filenameId,'_');
else
    DataIn = fullfile(parentdir,'DataT2OLWFSim','Vinf8dot0NTM_OL_Ct');
end
if ~isfolder(DataIn)
    warning('DataIn directory does not exist')
    return
end

%% Loop for  'Vinf8dot5_sowfa_2turb_yaw_alm_combined2.mat'
if yawmode == 0 %pitch control
    dirName = fullfile(DataIn, [filenameId(1:sepStr-1),'_OL_Ct'], filenameId);
    
    %     sepStr = strfind(filenameId,'_');
    %     dirNameId = fullfile(DataIn, [filenameVal(1:sepStr-1),'_OL_Ct'], filenameVal);
    tmpId = load(dirName,'Ct1', 'Ct2', 'FT1', 'FT2', 'PT1', 'PT2', ...
        'Ur1', 'Ur2', 'p', 'phi1', 'phi2', 'u', 'v');
    
    if useVal == 1
        sepStr = strfind(filenameVal,'_');
        
        filenameVal2 = strrep(filenameVal,'To','');
        dirNameVal = fullfile(DataIn, [filenameVal(1:sepStr-1),'_OL_Ct'], filenameVal2);
        tmpVal = load(dirNameVal,'Ct1', 'Ct2', 'FT1', 'FT2', 'PT1', 'PT2', ...
            'Ur1', 'Ur2', 'p', 'phi1', 'phi2', 'u', 'v');
    end
    
else
    
    dirName = fullfile(DataIn,[filenameId(1:sepStr-1),'_OL_Ct'],filenameId);
    tmpId = load(dirName,'Ct1', 'Ct2', 'FT1', 'FT2', 'PT1', 'PT2', ...
        'Ur1', 'Ur2', 'p', 'phi1', 'phi2', 'u', 'v');
    if useVal == 1
        dirNameVal = fullfile(DataIn,[filenameVal(1:9),'_OL_Ct'],filenameVal);
        tmpVal = load(dirNameVal,'Ct1', 'Ct2', 'FT1', 'FT2', 'PT1', 'PT2', ...
            'Ur1', 'Ur2', 'p', 'phi1', 'phi2', 'u', 'v');
    end
end

%% (1) ASSESS DATA
% easy solution to augment u flow field data matricx with other flow
% field data
fieldCell = fieldnames(tmpId);

useVal = exist('tmpVal','var') && useVal;

for idx = 1:length(fieldCell)
    if useVal
        tmp = [tmpId.(fieldCell{idx}), tmpVal.(fieldCell{idx})]; %#ok<NASGU>
    else
        tmp = [tmpId.(fieldCell{idx})]; %#ok<NASGU>
    end
    eval([fieldCell{idx},'= tmp;']);
end

if useVal
    tident = length(tmpId.Ct1);
else
    tident = floor(length(Ur1)*percentTrain);
end
tval = tident+t0;

Wp.turbine = struct(...
    'Crx',[0.4000    1.0321]*1e3,... % X-coordinates of turbines (m)
    'Cry',[400.0000  398.5747],... % Y-coordinates of turbines (m)
    'Drotor',126.4,... % Rotor diameter (m), note that WFSim only supports a uniform Drotor for now
    'powerscale',0.95,... % Turbine power scaling
    'forcescale',1.40 ... % Turbine force scaling
    );
Wp.mesh = struct(...
    'Lx',1882.1,... % Domain length in x-direction
    'Ly',800.0,... % Domain length in y-direction
    'Nx',50,... % Number of cells in x-direction
    'Ny',25 ... % Number of cells in y-direction
    );

nTx = round(Wp.turbine.Crx(1)/Wp.mesh.Lx * Wp.mesh.Nx);
nTy = round(Wp.turbine.Cry(1)/Wp.mesh.Ly * Wp.mesh.Ny);

kk = length(u)/Wp.mesh.Ny;

utemp = reshape(u,Wp.mesh.Nx,Wp.mesh.Ny,kk);
utempT = squeeze(utemp(nTx,nTy,:));
vtemp = reshape(v,Wp.mesh.Nx,Wp.mesh.Ny,kk);
vtempT = squeeze(vtemp(nTx,nTy,:));

QQ_u1 = utempT(t0:tident)'; %#ok<*USENS>
QQ_v1 = vtempT(t0:tident)';
QQ_p1 = p(:,t0:tident);

valid.QQ_u1 = utempT(tval:end)';
valid.QQ_v1 = vtempT(tval:end)';
valid.QQ_p1 = u(:,tval:end);

%% (2) DYNAMIC MODE DECOMPOSITION
states_u = QQ_u1;%fluid flow as states, identification data set
statesvalid_u = valid.QQ_u1; %fluid flow as states, validation data set for comparison

if yawmode == 0
    Inputs = [Ct1(t0:tident);Ct2(t0:tident)];
    Inputs_val = [Ct1(tval:end);Ct2(tval:end)];
else
    Inputs = [Ct1(t0:tident);Ct2(t0:tident); phi1(t0:tident)];
    Inputs_val = [Ct1(tval:end);Ct2(tval:end); phi1(tval:end)];
end

Outputs = [PT1(t0:tident);PT2(t0:tident);FT1(t0:tident);FT2(t0:tident)];
Outputs_val = [PT1(tval:end);PT2(tval:end);FT1(tval:end);FT2(tval:end)];

% ToDo check

strVAF = cell(length(noStates),1);
for idx = 1: length(noStates)
    
    n = noStates(idx);
    poly = PolyVec(idx);
    
    % subdirectories in dataOut
    if poly == 1 && yawmode == 0
        dirdmdName = 'eDMDresults_UasOutput_poly';
    elseif poly == 0 && yawmode == 0
        dirdmdName = 'eDMDresults_UasOutput';
    elseif poly == 1
        dirdmdName = 'eDMDresults_UasOutput_MIMO_poly';
    else
        dirdmdName = 'eDMDresults_UasOutput_MIMO';
    end
    
    dirdmd = fullfile(DataOut,dirdmdName);
    if ~ isfolder(dirdmd)
        mkdir(dirdmd);
    end
    
    Deterministic = [Ur1(t0:tident); Ur2(t0:tident)];
    Deterministic_val = [Ur1(tval:end); Ur2(tval:end)];
    %include non linear observables - Koopman extensions to better recover non linear dynamics
    if koopman == 1
        
        structPrev.Ur1_prev1 = Deterministic(1,1);
        structPrev.Ur2_prev1 = Deterministic(2,1);
        structPrev.dUr1_prev1 = 0;
        structPrev.dUr2_prev1 = 0;
        structPrev.M1(1) = structPrev.Ur1_prev1; % Moving mean
        structPrev.M2(1) = structPrev.Ur2_prev1;
        structPrev.k = 1;
        
        statesUr = nan(n+2+yawmode,size(Deterministic,2));
        for idxK = 1: length(Deterministic)
            [statesUr(:,idxK),stateNameUr,structPrev] = koopmanLiftedStates('urStates', [Deterministic(:,idxK);Inputs(:,idxK)],poly,n,structPrev,yawmode); %koopmanstateextensionRT(Deterministic(:,idxK),poly,n,structPrev);
        end
        
        windAtTurbine = [QQ_u1; QQ_v1];
        statesUV = koopmanLiftedStates('windStates', windAtTurbine,poly,n);
        
        stateNameUV = regexprep(stateNameUr,{'Ur1','Ur2','M'},{'u','v','Muv'});
        states = [statesUr]; %;statesUV];
        stateName = [stateNameUr]; %,';',stateNameUV];
        
       
        structPrev.Ur1_prev1 = Deterministic_val(1,1);
        structPrev.Ur2_prev1 = Deterministic_val(2,1);
        structPrev.dUr1_prev1 = 0;
        structPrev.dUr2_prev1 = 0;
        structPrev.M1(1) = Deterministic_val(1,1); % Moving mean
        structPrev.M2(1) = Deterministic_val(2,1);
        structPrev.k = 1;
        
        statesvalidUr = nan(n+2+yawmode,size(Deterministic_val,2));
        for idxK = 1: length(Deterministic_val)
            [statesvalidUr(:,idxK),st1,structPrev] = koopmanLiftedStates('urStates', [Deterministic_val(:,idxK);Inputs_val(:,idxK)],poly,n,structPrev); %koopmanstateextensionRT(Deterministic(:,idxK),poly,n,structPrev);
        end
        
        windAtTurbineVal = [valid.QQ_u1; valid.QQ_v1];
        statesUV_val = koopmanLiftedStates('windStates', windAtTurbineVal,poly,n);
        statesvalid = [statesvalidUr];%;statesUV_val];
    else
        states = states_u;
        statesvalid = statesvalid_u;
        temp = sprintf('u%d;',1: size(states,1));
        stateName = temp(1:end-1);
    end
    
    Vinf = utemp(1,1,t0);
    Vinf2 = utemp(1,1,tident);  Vstr = '';
    if abs(Vinf -Vinf2) >0.01 && isempty(NTMstr)
        Vstr = sprintf('To%2.1f',Vinf2);
    end
    folderName = strrep(sprintf('Vin_Vinf%2.1f%s%s_states%02d',Vinf,Vstr,NTMstr,n),'.','dot');%'Vinf%d_diffComb_states%d'
    fileName = sprintf('stateName_K%02d_P%d.mat',n,poly); %stateName);
    
    dirFig = fullfile(dirdmd,folderName); % this depends on vinf
    if ~exist(dirFig,'dir')
        mkdir(dirFig);
    end
    
    %     cellStName = regexp(stateName,';','split'); cellStName(end+1) ={'U2p'};
    %     aTable = array2table([states(:,1:end-1)',states(1,2:end)'],...
    %         'VariableNames',cellStName);
    %     aMatrix = [states(:,1:end-1)',Inputs(:,1:end-1)',states(2,2:end)'];
    
    % Sys id with ct1,ct2 as inputs
    [sys_red0,FITje0,Xd0,Xd_p0,xsim0,FITje_val0,xo_0,Koop0,...
        ysim,Inputs0,Deterministic0,ysim_val,Inputs_val0,Deterministic_val0]= eDMD_RT_Vin(states,stateName,statesvalid,poly,n,yawmode);
    strVal = 'Id.'; Vin = 0;
    koopmanPlots('effWind', ysim,Inputs0,Deterministic0, FITje0,dirFig,Vin,n,strVal);
    set(gcf,'name',strrep([strVal,fileName],'.stateName', ['Vin',num2str(Vin)]));
    koopmanPlots('effWind', ysim_val,Inputs_val0,Deterministic_val0,FITje_val0,dirFig,Vin,n);
    set(gcf,'name',strrep([strVal,fileName],'.stateName', ['Vin',num2str(Vin)]));
    
    % Sys id with ct1,ct2,V1 as inputs
    statesVinf = [states;statesUV(1,:)];
    statesvalidVinf = [statesvalid; statesUV_val(1,:)];
    stateNameVinf = [stateName,';V1'];
    [sys_red,FITje,Xd,Xd_p,xsim,FITje_val,xo,Koop,...
        ysim,InputsVinf,DeterministicVinf,ysim_val,Inputs_valOut,DeterministicVinf_val] = eDMD_RT_Vin(statesVinf,stateNameVinf,statesvalidVinf,poly,n,yawmode);
    strVal = 'Id.'; Vin = 1;
    koopmanPlots('effWind', ysim,InputsVinf,DeterministicVinf, FITje,dirFig,Vin,n,strVal);
    set(gcf,'name',strrep([strVal,fileName],'.stateName', ['Vin',num2str(Vin)]));
    strVal = 'Val.';
    koopmanPlots('effWind', ysim_val,Inputs_valOut,DeterministicVinf_val,FITje_val,dirFig,Vin,n,strVal);
    set(gcf,'name',strrep([strVal,fileName],'.stateName', ['Vin',num2str(Vin)]));
    
    save(fullfile(dirFig,fileName),'sys_red','FITje','FITje_val','dirName','xo','Koop');
    
    strVAF{idx} = sprintf('%s \t\t %d\t\t\t %d\t\t\t %d \t\t%d\n',num2str(n,'%02d'),...
        round(FITje(1)), round(FITje_val(1)),round(FITje(2)), round(FITje_val(2)));
    strVAF0{idx} = sprintf('%s \t\t %d\t\t\t %d\t\t\t %d \t\t%d\n',num2str(n,'%02d'),...
        round(FITje0(1)), round(FITje_val0(1)),round(FITje0(2)), round(FITje_val0(2)));
    
end

%% Print output to file and save results of offline sys id
% print this to workspace and into file

fprintf('Three Inputs: ct1, ct2, V \nNo K.\tT1(Id)\t\tT1(Val)\t\t T2(Id)\t\tT2(Val)\n');
for idx = 1: length(noStates)
    fprintf('%s',strVAF{idx});
end

fprintf('\n\n');
fprintf('Two Inputs: ct1, ct2\nNo K.\tT1(Id)\t\tT1(Val)\t\t T2(Id)\t\tT2(Val)\n');
for idx = 1: length(noStates)
    fprintf('%s',strVAF0{idx});
end

fid = fopen(['VAF_Vin_',strrep(filenameId,'.mat',''),'.txt'],'w');
fprintf(fid,'No K.\tT1(Id)\t\tT1(Val)\t\t T2(Id)\t\tT2(Val)\n');
for idx = 1: length(noStates)
    fprintf(fid,'%s',strVAF{idx});
end

fprintf(fid,'\n\n');
fprintf(fid,'No K.\tT1(Id)\t\tT1(Val)\t\t T2(Id)\t\tT2(Val)\n');
for idx = 1: length(noStates)
    fprintf(fid,'%s',strVAF0{idx});
end
fclose(fid);

Koffline = Koop;

end

%% ======================================================================
%% MainWfSimVinYawPowerOut  (task 'modelVinYawPower')
function MainWfSimVinYawPowerOut(yawmode,filenameId,filenameVal,noStates,PolyVec,useVal,percentTrain,NTMstr,t0)% Code based on
% Model Predictive Control For Wake Steering: a Koopman Dynamic Mode
% Decomposition Approach
% Master Thesis Dissertation
% Author: Nassir Rodrigues Cassamo
% Supervisors: Professor Jan-Willem Van Wingerden and Professor Jo?o Sousa
% First official commit (1.0): December 2020

%% RELEVAT INFORMATION
%This script has the following goals:
% (1) Assess Simulation data, both qualitatively (animations) and
% quantitatively (graphics)
% (2) Derives a low dimensional model using Dynamic Mode Decomposition
% (variations included to take into account input-output data and other
% known states - deterministic states)
% (3) Validates the models with a set of validaiton data

%% (0) INITIALISE
% define type of simulation
if ~nargin
    
    % directory for identification and validation data, yaw and pitch
    yawmode = 1; %0 for ct (pitch/torque) control, 1 for additional yaw control
    % User input
    noStates = 6; %[6,12,12,14,16,18,24];
    useVal = 1; % use extra data for validation
    percentTrain = .9;
    t0 = 240;
    
    NTMstr = '';
%     filenameId = 'Vinf8dot0_sowfa_2turb_yaw_noise_step3'; %Vinf8dot0_sowfa_2turb_yaw_noise2';
%     filenameVal = 'Vinf8dot0_sowfa_2turb_yaw_steps_Ct_comb';
%     
    filenameId = 'Vinf8dot0_sowfa_2turb_yaw_noise_step';
    filenameVal = 'Vinf8dot0_sowfa_2turb_yaw_steps'; %'Vinf8dot0_sowfa_2turb_alm_turbl.mat'; %
       
    PolyVec = zeros(size(noStates)); % 1: use only polynomial, 0: otherwise
    % idxPolyOn = find(noStates == 12, 1, 'first');
    % PolyVec(idxPolyOn) = 1;
end

% method = -2; %0: DMD ; 1:DMDc; 2:IODMD; 3:EXTIODMD -1: EDMD
koopman = 1; %to add deterministic states to flow field data

% Simulation characterisitc (resampling)
dt = 1; %time sampling (s)

KoopmanDir = fileparts(mfilename('fullpath')); % Get Koopman directory
pathFunctions = genpath(fullfile(KoopmanDir,'Functions'));

addpath(fullfile(KoopmanDir,'2.DYNAMIC_MODE_DECOMPOSITION'),...
    pathFunctions);
codedir = mfilename('fullpath');
parentdir = fileparts(fileparts(codedir));

% DEFINE MAIN output data DIRECTORY TO STORE ALL RESULTS
DataOut = fullfile(parentdir,'DataInOutWfSim');
if ~isfolder(DataOut)
    mkdir(DataOut)
end

if isempty(NTMstr)
    DataIn = fullfile(parentdir,'DataT2OLWFSim'); %,...
    sepStr = strfind(filenameId,'_');
else
    DataIn = fullfile(parentdir,'DataT2OLWFSim','Vinf8dot0NTM_OL_Ct');
end
if ~isfolder(DataIn)
    warning('DataIn directory does not exist')
    return
end

%% Load data from corresponding subdirectory
dirName = fullfile(DataIn, [filenameId(1:sepStr-1),'_OL_Ct'], filenameId);
tmpId = load(dirName,'Ct1', 'Ct2', 'FT1', 'FT2', 'PT1', 'PT2', ...
    'Ur1', 'Ur2', 'p', 'phi1', 'phi2', 'u', 'v');

if useVal == 1 % load different data set for validation
    sepStr = strfind(filenameVal,'_');
    filenameVal2 = strrep(filenameVal,'To','');
    dirNameVal = fullfile(DataIn, [filenameVal(1:sepStr-1),'_OL_Ct'], filenameVal2);
    tmpVal = load(dirNameVal,'Ct1', 'Ct2', 'FT1', 'FT2', 'PT1', 'PT2', ...
        'Ur1', 'Ur2', 'p', 'phi1', 'phi2', 'u', 'v');
end


%% (1) ASSESS DATA
% easy solution to augment u flow field data matricx with other flow
% field data

% Define layout
Wp.turbine = struct(...
    'Crx',[0.4000    1.0321]*1e3,... % X-coordinates of turbines (m)
    'Cry',[400.0000  398.5747],... % Y-coordinates of turbines (m)
    'Drotor',126.4,... % Rotor diameter (m), note that WFSim only supports a uniform Drotor for now
    'powerscale',0.95,... % Turbine power scaling
    'forcescale',1.40 ... % Turbine force scaling
    );
Wp.mesh = struct(...
    'Lx',1882.1,... % Domain length in x-direction
    'Ly',800.0,... % Domain length in y-direction
    'Nx',50,... % Number of cells in x-direction
    'Ny',25 ... % Number of cells in y-direction
    );

nTx = floor(Wp.turbine.Crx(1)/Wp.mesh.Lx * Wp.mesh.Nx);
nTy = round(Wp.turbine.Cry(1)/Wp.mesh.Ly * Wp.mesh.Ny);

fieldCell = fieldnames(tmpId);

useVal = exist('tmpVal','var') && useVal;


for idx = 1:length(fieldCell)
    if useVal
        idxPhiVal = tmpVal.phi1 <= 30;
        if length(tmpVal.(fieldCell{idx})) == length(idxPhiVal)
            tmp = [tmpId.(fieldCell{idx}), tmpVal.(fieldCell{idx})(idxPhiVal)]; %#ok<NASGU>
        else
            tmp = [tmpId.(fieldCell{idx}), tmpVal.(fieldCell{idx})];
        end
    else
        tmp = [tmpId.(fieldCell{idx})]; %#ok<NASGU>
    end
    eval([fieldCell{idx},'= tmp;']);
end

% phiTemp = cos(phi1/180 *pi); %phi1/180 *pi; %
% phi1 = (phiTemp - mean(phiTemp))/std(phiTemp);

if useVal
    tident = length(tmpId.Ct1);
else
    tident = floor(length(Ur1)*percentTrain);
end
tval = tident+t0;

kk = length(u)/Wp.mesh.Ny;

utemp = reshape(u,Wp.mesh.Nx,Wp.mesh.Ny,kk);
utempT = squeeze(utemp(nTx,nTy,:));
vtemp = reshape(v,Wp.mesh.Nx,Wp.mesh.Ny,kk);
vtempT = squeeze(vtemp(nTx,nTy,:));
ptemp = reshape(p,Wp.mesh.Nx,Wp.mesh.Ny,kk);
ptempT = squeeze(ptemp(nTx,nTy,:));

QQ_u1 = utempT(t0:tident)'; %#ok<*USENS>
QQ_v1 = vtempT(t0:tident)';
QQ_p1 = ptempT(t0:tident);

valid.QQ_u1 = utempT(tident+1:end)';
valid.QQ_v1 = vtempT(tident+1:end)';
valid.QQ_p1 = ptempT(tident+1:end)';

if ~useVal
    idxPhiVal = phi1(tident+1:end) <= 40;
end
tempQ_u1 = valid.QQ_u1(idxPhiVal);
tempQ_v1 = valid.QQ_v1(idxPhiVal);
tempQ_p1 = valid.QQ_p1(idxPhiVal);
valid.QQ_u1 = tempQ_u1(t0:end);
valid.QQ_v1 = tempQ_v1(t0:end);
valid.QQ_w1 = tempQ_v1(t0:end);

%% (2) DYNAMIC MODE DECOMPOSITION
states_u = QQ_u1;%fluid flow as states, identification data set
statesvalid_u = valid.QQ_u1; %fluid flow as states, validation data set for comparison

if yawmode == 0
    Inputs = [Ct1(t0:tident);Ct2(t0:tident)];
    Inputs_val = [Ct1(tval:end);Ct2(tval:end)];
else
    Inputs = [Ct1(t0:tident);Ct2(t0:tident); phi1(t0:tident)];
    Inputs_val = [Ct1(tval:end);Ct2(tval:end); phi1(tval:end)];
end

Outputs = [PT1(t0:tident);PT2(t0:tident);FT1(t0:tident);FT2(t0:tident)];
Outputs_val = [PT1(tval:end);PT2(tval:end);FT1(tval:end);FT2(tval:end)];
Psum = PT1(t0:tident)+PT2(t0:tident);

% ToDo check

strVAF = cell(length(noStates),1);
for idx = 1: length(noStates)
    
    n = noStates(idx);
    poly = PolyVec(idx);
    
    % subdirectories in dataOut
    if poly == 1 && yawmode == 0
        dirdmdName = 'eDMDresults_UasOutput_poly_PT';
    elseif poly == 0 && yawmode == 0
        dirdmdName = 'eDMDresults_UasOutput_PT';
    elseif poly == 1
        dirdmdName = 'eDMDresults_UasOutput_MIMO_poly_PT';
    else
        dirdmdName = 'eDMDresults_UasOutput_MIMO_PT';
    end
    
    dirdmd = fullfile(DataOut,dirdmdName);
    if ~ isfolder(dirdmd)
        mkdir(dirdmd);
    end
    
    Deterministic = [Ur1(t0:tident); Ur2(t0:tident)];
    Deterministic_val = [Ur1(tval:end); Ur2(tval:end)];
    %include non linear observables - Koopman extensions to better recover non linear dynamics
    if koopman == 1
        
        structPrev.Ur1_prev1 = Deterministic(1,1);
        structPrev.Ur2_prev1 = Deterministic(2,1);
        structPrev.dUr1_prev1 = 0;
        structPrev.dUr2_prev1 = 0;
        structPrev.M1(1) = structPrev.Ur1_prev1; % Moving mean
        structPrev.M2(1) = structPrev.Ur2_prev1;
        structPrev.k = 1;
        
        statesUr = nan(n+2+yawmode,size(Deterministic,2));
        for idxK = 1: length(Deterministic)
            [statesUr(:,idxK),stateNameUr,structPrev] = koopmanLiftedStates('urStates', [Deterministic(:,idxK);Inputs(:,idxK)],poly,n,structPrev,yawmode); %koopmanstateextensionRT(Deterministic(:,idxK),poly,n,structPrev);
        end
        
        windAtTurbine = [QQ_u1; QQ_v1];
        statesUV = koopmanLiftedStates('windStates', windAtTurbine,poly,n);
        
        stateNameUV = regexprep(stateNameUr,{'Ur1','Ur2','M'},{'u','v','Muv'});
        states = [statesUr]; %;statesUV];
        stateName = [stateNameUr]; %,';',stateNameUV];
        
        structPrev.Ur1_prev1 = Deterministic_val(1,1);
        structPrev.Ur2_prev1 = Deterministic_val(2,1);
        structPrev.dUr1_prev1 = 0;
        structPrev.dUr2_prev1 = 0;
        structPrev.M1(1) = Deterministic_val(1,1); % Moving mean
        structPrev.M2(1) = Deterministic_val(2,1);
        structPrev.k = 1;
        
        statesvalidUr = nan(n+2+yawmode,size(Deterministic_val,2));
        for idxK = 1: length(Deterministic_val)
            [statesvalidUr(:,idxK),st1,structPrev] = koopmanLiftedStates('urStates', [Deterministic_val(:,idxK);Inputs_val(:,idxK)],poly,n,structPrev); %koopmanstateextensionRT(Deterministic(:,idxK),poly,n,structPrev);
        end
        
        windAtTurbineVal = [valid.QQ_u1; valid.QQ_v1];
        statesUV_val = koopmanLiftedStates('windStates', windAtTurbineVal,poly,n);
        statesvalid = [statesvalidUr];%;statesUV_val];
    else
        states = states_u;
        statesvalid = statesvalid_u;
        temp = sprintf('u%d;',1: size(states,1));
        stateName = temp(1:end-1);
    end
    
    Vinf = utemp(1,1,t0);
    Vinf2 = utemp(1,1,tident);  Vstr = '';
    if abs(Vinf -Vinf2) >0.01 && isempty(NTMstr)
        Vstr = sprintf('To%2.1f',Vinf2);
    end
    folderName = strrep(sprintf('Vin_Vinf%2.1f%s%s_states%02d',Vinf,Vstr,NTMstr,n),'.','dot');%'Vinf%d_diffComb_states%d'
    fileName = sprintf('stateName_K%02d_P%d.mat',n,poly); %stateName);
    
    dirFig = fullfile(dirdmd,folderName); % this depends on vinf
    if ~exist(dirFig,'dir')
        mkdir(dirFig);
    end
    
    %     cellStName = regexp(stateName,';','split'); cellStName(end+1) ={'U2p'};
    %     aTable = array2table([states(:,1:end-1)',states(1,2:end)'],...
    %         'VariableNames',cellStName);
    %     aMatrix = [states(:,1:end-1)',Inputs(:,1:end-1)',states(2,2:end)'];
    statesP = [Outputs(1:2,:);states];
    statesvalidP = [Outputs_val(1:2,:)/10^6;statesvalid];
    stateNameP = ['P1;P2;',stateName];
    
    %     aTable = array2table(statesP');
    %     tmp = strsplit(stateNameP,';');
    %     aTable.Properties.VariableNames = tmp';
    %
    % Sys id with ct1,ct2 as inputs
    %     [sys_red0,FITje0,Xd0,Xd_p0,xsim0,FITje_val0,xo_0,Koop0,...
    %         ysim,Inputs0,Deterministic0,ysim_val,Inputs_val0,Deterministic_val0] = ...
    %         eDMD_RT_Vin_P(statesP,stateNameP,statesvalidP,poly,n+2,yawmode);
    
    [outStatesAndInputsP] = [sum(Outputs(1:2,:)); states];
    [outStatesAndInputsPvalid] = [sum(Outputs_val(1:2,:)); statesvalid];
    outStateName = ['PT; ',stateName];
    
    [sys_red0,FITje0,Xd0,Xd_p0,xsim0,FITje_val0,xo_0,Koop0,...
        ysim,Inputs0,Deterministic0,ysim_val,Inputs_val0,Deterministic_val0] = ...
        eDMD_RT_Vin_P(outStatesAndInputsP,outStateName,outStatesAndInputsPvalid,poly,n,yawmode);
    
    
    %eDMD_RT_Vin_P(statesP,stateNameP,statesvalidP,poly,n+2,yawmode);
    strVal = 'Id.'; Vin = 0; displayOff = 0; Vinf = '';strOut = 'P';
    koopmanPlots('effWind', ysim,Inputs0,Deterministic0, FITje0,dirFig,Vin,n,strVal,displayOff,Vinf,strOut);
    set(gcf,'name',strrep([strVal,fileName],'.stateName', ['Vin',num2str(Vin)]));
    strVal = 'Val.';
    koopmanPlots('effWind', ysim_val,Inputs_val0,Deterministic_val0,FITje_val0,dirFig,Vin,n,strVal,displayOff,Vinf,strOut);
    set(gcf,'name',strrep([strVal,fileName],'.stateName', ['Vin',num2str(Vin)]));
    
    % Sys id with ct1,ct2,V1 as inputs
    statesVinfP = [outStatesAndInputsP;statesUV(1,:)]; %[statesP;statesUV(1,:)];
    statesvalidVinfP = [outStatesAndInputsPvalid; statesUV_val(1,:)]; %[statesvalidP; statesUV_val(1,:)];
    stateNameVinfP =  [outStateName,';V1']; %[stateNameP,';V1'];
    
    [sys_red,FITje,Xd,Xd_p,xsim,FITje_val,xo,Koop,...
        ysim,InputsVinf,DeterministicVinf,ysim_val,Inputs_valOut,DeterministicVinf_val] = ...
        eDMD_RT_Vin_P(statesVinfP,stateNameVinfP,statesvalidVinfP,poly,n,yawmode);
    strVal = 'Id.'; Vin = 1;
    koopmanPlots('effWind', ysim,InputsVinf,DeterministicVinf, FITje,dirFig,Vin,n,strVal,displayOff,Vinf,strOut);
    set(gcf,'name',strrep([strVal,fileName],'.stateName', ['Vin',num2str(Vin)]));
    set(gcf,'name', [strVal,fileName])
    strVal = 'Val.';
    koopmanPlots('effWind', ysim_val,Inputs_valOut,DeterministicVinf_val,FITje_val,dirFig,Vin,n,strVal,displayOff,Vinf,strOut);
    set(gcf,'name',strrep([strVal,fileName],'.stateName', ['Vin',num2str(Vin)]));   
    %     statesVinfAll = [statesVinfP, statesUV(1,:)]; %[statesP;statesUV(1,:)];
    %     statesvalidVinfP = [outStatesAndInputsPvalid; statesUV_val(1,:)]; %[statesvalidP; statesUV_val(1,:)];
    %     stateNameVinfP =  [outStateName,';V1']; %[stateNameP,';V1'];
     
    save(fullfile(dirFig,fileName),'sys_red','FITje','FITje_val','dirName','xo','Koop');
    
    %     strVAF{idx} = sprintf('%s \t\t %d\t\t\t %d\t\t\t %d \t\t%d\n',num2str(n,'%02d'),...
    %         round(FITje(1)), round(FITje_val(1)),round(FITje(2)), round(FITje_val(2)));
    %     strVAF0{idx} = sprintf('%s \t\t %d\t\t\t %d\t\t\t %d \t\t%d\n',num2str(n,'%02d'),...
    %         round(FITje0(1)), round(FITje_val0(1)),round(FITje0(2)), round(FITje_val0(2)));
    %
end

end

%% ======================================================================
%% testPowerInputsOnly  (task 'powerData')
function [R, D] = testPowerInputsOnly(filenameId,filenameVal,noStates,t0)
% testPowerInputsOnly compares Koopman power models for constant free
% stream wind: the model of MainWfSimVinYawPowerOut with the measured wind
% speed V1 (u in the grid cell of turbine 1) as known input, and models with
% the control inputs ct1, ct2, phi1 only. Data, lifting and identification
% as in MainWfSimVinYawPowerOut (eDMD_RT_Vin_P, yawmode 1).
%
% Models (output: farm power PT1 + PT2)
%   'V1':   inputs ct1, ct2, phi1, V1 (as in MainWfSimVinYawPowerOut)
%   'u':    inputs ct1, ct2, phi1 (no wind speed, no offset)
%   'u+1':  inputs ct1, ct2, phi1 and a constant one (offset); for constant
%           wind the free stream wind speed is part of this offset
%
% Each model is simulated over the validation set from its first sample
%   'nonlin': lifting of the simulated Ur1, Ur2 in each step (eDMD_RT_Vin_P)
%   'lin':    linear state space model sys_red (as used in the MPC)
%   'lin,ss': as 'lin', initial state the steady state for the first inputs
%             (no measurement at all, model 'u+1' only)

if ~nargin
    filenameId = 'Vinf8dot0_sowfa_2turb_yaw_noise_step';
    filenameVal = 'Vinf8dot0_sowfa_2turb_yaw_steps';
    noStates = 6;
    t0 = 240;
end
yawmode = 1;
poly = 0;
n = noStates;

KoopmanDir = fileparts(mfilename('fullpath'));
% (path: see WFSim_addpaths)
DataIn = fullfile(fileparts(KoopmanDir),'DataT2OLWFSim');
dirFig = fullfile(fileparts(KoopmanDir),'DataInOutWfSim','PowerInputsOnly');
if ~isfolder(dirFig)
    mkdir(dirFig);
end

%% Load data (as MainWfSimVinYawPowerOut, useVal = 1)
vars = {'Ct1','Ct2','PT1','PT2','Ur1','Ur2','phi1','u'};
if ischar(filenameId) % several identification sets (cell): prepared from t0 each and appended
    filenameId = {filenameId};
end
Ids = cell(size(filenameId));
for idxId = 1:numel(filenameId)
    sepStr = strfind(filenameId{idxId},'_');
    Ids{idxId} = load(fullfile(DataIn,[filenameId{idxId}(1:sepStr-1),'_OL_Ct'],filenameId{idxId}),vars{:});
end
sepStr = strfind(filenameVal,'_');
Val = load(fullfile(DataIn,[filenameVal(1:sepStr-1),'_OL_Ct'],filenameVal),vars{:});

Wp.mesh = struct('Lx',1882.1,'Ly',800.0,'Nx',50,'Ny',25);
Crx1 = 0.4000e3; Cry1 = 400;
nTx = floor(Crx1/Wp.mesh.Lx*Wp.mesh.Nx);
nTy = round(Cry1/Wp.mesh.Ly*Wp.mesh.Ny);

D = struct();
[D.id, Vinf] = prepareSet(Ids{1}, t0, nTx, nTy, Wp, true(size(Ids{1}.Ct1)), n);
for idxId = 2:numel(Ids)
    Sadd = prepareSet(Ids{idxId}, t0, nTx, nTy, Wp, true(size(Ids{idxId}.Ct1)), n);
    for fn = {'V1','P','P12','lifted'}
        D.id.(fn{1}) = [D.id.(fn{1}), Sadd.(fn{1})];
    end
end
D.val = prepareSet(Val, t0, nTx, nTy, Wp, Val.phi1 <= 30, n);
fprintf('Free stream wind: %.2f m/s (id.), V1 from %.2f to %.2f m/s\n', Vinf, min(D.id.V1), max(D.id.V1));

%% Identify and simulate the models
models = {'V1','u','u+1'};
R = table();
for idxM = 1:numel(models)
    S = {D.id, D.val};
    for k = 1:2
        switch models{idxM}
            case 'V1',  S{k}.extra = S{k}.V1;
            case 'u',   S{k}.extra = zeros(0,numel(S{k}.P));
            case 'u+1', S{k}.extra = ones(1,numel(S{k}.P));
        end
        S{k}.Z = [S{k}.P; S{k}.lifted; S{k}.extra];
    end
    stateName = ['PT;', S{1}.stateName, repmat(';extra', 1, size(S{1}.extra,1))];
    [sys_red,FITid,~,~,~,FITval,~,~,~,~,~,ysim_val] = ...
        eDMD_RT_Vin_P(S{1}.Z,stateName,S{2}.Z,poly,n,yawmode);

    % Linear model (MPC): inputs [ct1;ct2;phi1;extra], x0 from the data
    uVal = S{2}.Z(n+2:end,1:end-1);
    pVal = S{2}.P(1:end-1);
    x0 = S{2}.lifted(1:n,1);
    yLin = lsim(sys_red, uVal', [], x0)';
    row = {models{idxM}, FITid, FITval, vaf(pVal(1:end-1), ysim_val), rmse(pVal(1:end-1), ysim_val), ...
        vaf(pVal, yLin), rmse(pVal, yLin), NaN, NaN, max(abs(eig(sys_red.A)))};
    yLinSS = nan(size(yLin));
    if strcmp(models{idxM},'u+1')
        xss = (eye(n) - sys_red.A) \ (sys_red.B*uVal(:,1));
        yLinSS = lsim(sys_red, uVal', [], xss)';
        row(8:9) = {vaf(pVal, yLinSS), rmse(pVal, yLinSS)};
    end
    R = [R; cell2table(row, 'VariableNames', {'Model','VAFid','VAFval','VAFval_nonlin', ...
        'RMSEval_nonlin_MW','VAFval_lin','RMSEval_lin_MW','VAFval_linSS','RMSEval_linSS_MW','maxAbsEig'})]; %#ok<AGROW>
    Y.(matlab.lang.makeValidName(models{idxM})) = struct('nonlin', ysim_val, 'lin', yLin, 'linSS', yLinSS, 'sys', sys_red);
end
disp(R)

%% Plot validation
t = (0:numel(D.val.P)-2);
cl = lines(4);
fig = figure('Name','PowerInputsOnly');
ax = subplot(3,1,1);
plot(ax, t, D.val.lifted(n+1:n+2,1:end-1)'); hold(ax,'on');
yyaxis(ax,'right'); plot(ax, t, D.val.lifted(n+3,1:end-1),'k'); ylabel(ax,'\phi_1 (deg)');
yyaxis(ax,'left'); ylabel(ax,'C_T'''); legend(ax,'C_{T,1}''','C_{T,2}''','\phi_1'); grid(ax,'on'); axis(ax,'tight');
title(ax, sprintf('Validation %s, V_\\infty = %.1f m/s', strrep(filenameVal,'_','\_'), Vinf));
ax = subplot(3,1,2);
plot(ax, t, D.val.P(1:end-1)/1e6, 'k', 'LineWidth', 1.5); hold(ax,'on');
plot(ax, t, Y.V1.lin/1e6, '--', 'Color', cl(1,:));
plot(ax, t, Y.u.lin/1e6, ':', 'Color', cl(2,:), 'LineWidth', 1.2);
plot(ax, t, Y.u_1.lin/1e6, '-', 'Color', cl(3,:));
plot(ax, t, Y.u_1.linSS/1e6, '-.', 'Color', cl(4,:));
ylabel(ax,'P_{WF} (MW)'); grid(ax,'on'); axis(ax,'tight');
legend(ax,'WFSim', ...
    sprintf('V1 input (%.0f%%)', R.VAFval_lin(1)), sprintf('inputs only (%.0f%%)', R.VAFval_lin(2)), ...
    sprintf('inputs + offset (%.0f%%)', R.VAFval_lin(3)), sprintf('inputs + offset, x_0 = x_{ss} (%.0f%%)', R.VAFval_linSS(3)), ...
    'Location','best');
title(ax,'Linear model (as in the MPC), VAF in brackets');
ax = subplot(3,1,3);
plot(ax, t, (Y.V1.lin - D.val.P(1:end-1))/1e6, '--', 'Color', cl(1,:)); hold(ax,'on');
plot(ax, t, (Y.u_1.lin - D.val.P(1:end-1))/1e6, '-', 'Color', cl(3,:));
plot(ax, t, (Y.u_1.linSS - D.val.P(1:end-1))/1e6, '-.', 'Color', cl(4,:));
ylabel(ax,'error (MW)'); xlabel(ax,'time (s)'); grid(ax,'on'); axis(ax,'tight');
legend(ax,'V1 input','inputs + offset','inputs + offset, x_0 = x_{ss}','Location','best');
set(fig,'Position',[100 100 700 700]);
set(fig,'PaperUnits','centimeters','PaperPositionMode','manual','PaperSize',[17 17],'PaperPosition',[0 0 17 17]);
print(fig, fullfile(dirFig,'PowerInputsOnly_val'), '-dpng', '-r200');
print(fig, fullfile(dirFig,'PowerInputsOnly_val'), '-depsc', '-painters', '-loose');
save(fullfile(dirFig,'PowerInputsOnly.mat'), 'R', 'Y', 'D');
end

function [S, Vinf] = prepareSet(T, t0, nTx, nTy, Wp, idxKeep, n)
% prepareSet lifted states (K_psi) and power of one data set from t0 on
kk = length(T.u)/Wp.mesh.Ny;
utemp = reshape(T.u, Wp.mesh.Nx, Wp.mesh.Ny, kk);
V1 = squeeze(utemp(nTx,nTy,:))';
Vinf = utemp(1,1,t0);
f = @(x) x(idxKeep);
Ur = [f(T.Ur1); f(T.Ur2)];
inp = [f(T.Ct1); f(T.Ct2); f(T.phi1)];
V1 = f(V1);
idx = t0:numel(Ur(1,:));
Ur = Ur(:,idx); inp = inp(:,idx);
S.V1 = V1(idx);
S.P12 = [f(T.PT1); f(T.PT2)]; S.P12 = S.P12(:,idx); % turbine powers (W)
S.P = sum(S.P12, 1);
structPrev = struct('Ur1_prev1',Ur(1,1),'Ur2_prev1',Ur(2,1),'dUr1_prev1',0,'dUr2_prev1',0, ...
    'M1',Ur(1,1),'M2',Ur(2,1),'k',1);
S.lifted = nan(n+3, numel(idx));
for k = 1:numel(idx)
    [S.lifted(:,k),S.stateName,structPrev] = koopmanLiftedStates('urStates', [Ur(:,k); inp(:,k)], 0, n, structPrev, 1);
end
end

function e = rmse(y, yhat)
e = sqrt(mean((y(:) - yhat(:)).^2))/1e6;
end

%% ======================================================================
%% testPowerStates  (task 'powerStates')
function [T, Mdl] = testPowerStates(filenameId, filenameVal, featureSets, t0, tau, stateSets)
% testPowerStates compares Koopman states of the input-only farm power
% model: the lifted rotor effective wind speeds Ur (as K_psi) and liftings
% of the turbine powers P1, P2 of the previous time step (MW), which are
% measured directly. Inputs: features of powerFeatures (ct1, ct2, phi1).
%
% Model: x(k+1) = A x(k) + B w(k), P(k) = C x(k) + D w(k), identified by
% least squares as testPowerFeatures, simulated linearly over the
% validation set from the measured initial state (as in the MPC).
%
% Feature sets: see powerFeatures; '<set>+V1' appends the measured wind
% speed V1 (u in front of turbine 1) as input, 'u+V1' = ct1, ct2, phi1, V1
% (inputs of MainWfSimVinYawPowerOut, no constant).
%
% States x(k) (Pi = Pi(k-1) in MW, cti = cti(k-1)):
%   'Ur n=4'        Ur1, Ur2, Ur1^2, Ur2^2 (reference, as exportPowerFeatureModel)
%   'Ur n=6'        'Ur n=4', Ur1^3, Ur2^3 (K_psi)
%   'Ur n=8'        'Ur n=6', DUr1, DUr2 (K_psi: deviation from the moving
%                   mean of the past Ur over 25 s)
%   'Ur n=4+d', 'Ur n=6+d', 'Ur n=8+d'   plus the delayed Ur1(k-1), Ur2(k-1)
%   'Ur n=6+D<N>'   'Ur n=6' plus DUr1, DUr2 with the moving mean over N
%                   samples (recursion of K_psi, which uses N = 25)
%   'Ur n=6+D<N>,M' as 'Ur n=6+D<N>' plus the moving means M1, M2
%   'Ur n=6+M<N>fix' 'Ur n=6' plus M1, M2; the rows of M in A are not fitted
%                   but set to the known recursion M(k+1) = (1-1/N)*M(k) + Ur(k)/N
%   'Ur n=4+Ur1(25,50,100)'  'Ur n=4' plus the shift register Ur1(k-1..k-dmax)
%                   (fixed shift rows); only the listed delays enter the Ur rows
%                   and the output; Ur1 before the first sample held at Ur1(1)
%   'P'             P1, P2
%   'P,sqrtP'       P1, P2, sqrt(P1), sqrt(P2)
%   'P,cbrtP'       P1, P2, P1^(1/3), P2^(1/3)
%   'P,Ue'          P1, P2, (P1/ct1)^(1/3), (P2/ct2)^(1/3)  (~ Ur, P = c*ct'*Ue^3)
%   'P,sqrtP,cbrtP' P1, P2, sqrt(P1), sqrt(P2), P1^(1/3), P2^(1/3)

if nargin < 1 || isempty(filenameId)
    filenameId = 'Vinf8dot0_sowfa_2turb_yaw_noise_step';
end
if nargin < 2 || isempty(filenameVal)
    filenameVal = {'Vinf8dot0_sowfa_2turb_yaw_steps','Vinf8dot0_sowfa_2turb_yaw_steps_Ct_comb'};
end
if ischar(filenameVal), filenameVal = {filenameVal}; end
if nargin < 3 || isempty(featureSets), featureSets = {'u+1','NL+1+M+x'}; end
if nargin < 4, t0 = 240; end
if nargin < 5, tau = 100; end

mainDir = fileparts(fileparts(mfilename('fullpath')));
% (path: see WFSim_addpaths) % powerFeatures
Dall = cell(numel(filenameVal),1);
for idxV = 1:numel(filenameVal)
    evalc('[~, Dall{idxV}] = testPowerInputsOnly(filenameId, filenameVal{idxV}, 8, t0);');
    close all
end

if nargin < 6 || isempty(stateSets)
    stateSets = {'Ur n=4','P','P,sqrtP','P,cbrtP','P,Ue','P,sqrtP,cbrtP'};
end
T = table();
Mdl = struct('A',{},'B',{},'C',{},'D',{},'n',{},'stateSet',{},'featureSet',{},'tau',{},'maxAbsEig',{});
for idxF = 1:numel(featureSets)
    for idxS = 1:numel(stateSets)
        [X, W, P] = dataMatrices(Dall{1}.id, stateSets{idxS}, featureSets{idxF}, tau);
        n = size(X,1);
        Psi = [X(:,1:end-1); W(:,1:end-1)];
        Y = [X(:,2:end); P(1:end-1)];
        sc = sqrt(mean(Psi.^2, 2));
        K = (Y*(Psi./sc)')/((Psi./sc)*(Psi./sc)' + 1e-8*size(Psi,2)*eye(size(Psi,1)));
        K = K./sc';
        tokFix = regexp(stateSets{idxS},'\+M(\d+)fix$','tokens','once');
        if ~isempty(tokFix) % known moving mean dynamics (last two states)
            N = str2double(tokFix{1});
            K(n-1:n,:) = 0;
            K(n-1,[1,n-1]) = [1/N, 1-1/N];
            K(n,[2,n]) = [1/N, 1-1/N];
        end
        tokDel = regexp(stateSets{idxS},'^Ur n=(\d+)\+Ur1\(([\d,]+)\)$','tokens','once');
        if ~isempty(tokDel) % shift register Ur1(k-1..k-dmax): fixed shift rows,
            % only the delays dl enter the Ur rows and the output (refit)
            nL = str2double(tokDel{1});
            dl = str2double(split(tokDel{2},','))';
            iSel = [1:nL, nL + dl, n+1:size(Psi,1)];
            iFree = [1:nL, n+1];
            sc = sqrt(mean(Psi(iSel,:).^2, 2));
            Ps = Psi(iSel,:)./sc;
            K(iFree,:) = 0;
            K(iFree,iSel) = ((Y(iFree,:)*Ps')/(Ps*Ps' + 1e-8*size(Ps,2)*eye(numel(iSel))))./sc';
            K(nL+1:n,:) = 0;
            K(nL+1,1) = 1;
            K(sub2ind(size(K), nL+2:n, nL+1:n-1)) = 1;
        end
        sys = ss(K(1:n,1:n), K(1:n,n+1:end), K(n+1,1:n), K(n+1,n+1:end), 1);
        Mdl(end+1) = struct('A',sys.A,'B',sys.B,'C',sys.C,'D',sys.D,'n',n, ... %#ok<AGROW>
            'stateSet',stateSets{idxS},'featureSet',featureSets{idxF},'tau',tau,'maxAbsEig',max(abs(eig(sys.A))));
        for idxV = 1:numel(filenameVal)
            [Xv, Wv, Pv] = dataMatrices(Dall{idxV}.val, stateSets{idxS}, featureSets{idxF}, tau);
            y = lsim(sys, Wv', [], Xv(:,1))';
            T = [T; table(string(featureSets{idxF}), string(stateSets{idxS}), ...
                string(strrep(filenameVal{idxV},'Vinf8dot0_sowfa_2turb_','')), ...
                vaf(Pv, y), sqrt(mean((y - Pv).^2)), mean(y - Pv), max(abs(eig(sys.A))), ...
                'VariableNames', {'Features','States','ValSet','VAF','RMSE_MW','MeanErr_MW','maxAbsEig'})]; %#ok<AGROW>
        end
    end
end
disp(T)
dirOut = fullfile(mainDir,'DataInOutWfSim','PowerInputsOnly');
writetable(T, fullfile(dirOut,'PowerStates.csv'));
end

function [X, W, P] = dataMatrices(S, stateSet, featureSet, tau)
% States from k = 2 on (power states use the previous time step)
inp = S.lifted(end-2:end,:);
if endsWith(featureSet,'+V1') % measured wind speed V1 as additional input
    baseSet = extractBefore(featureSet, strlength(featureSet) - 2);
    if strcmp(baseSet,'u') % as the paper: ct1, ct2, phi1, V1 (no constant)
        W = [inp; S.V1];
    else
        W = [koopmanLiftedStates('powerFeatures', inp, baseSet, tau); S.V1];
    end
else
    W = koopmanLiftedStates('powerFeatures', inp, featureSet, tau);
end
P12 = S.P12/1e6;
Pp = P12(:,1:end-1);      % P(k-1)
ctp = inp(1:2,1:end-1);   % ct(k-1)
switch stateSet
    case {'Ur n=4','Ur n=6','Ur n=8'}
        X = S.lifted(1:sscanf(stateSet,'Ur n=%d'),2:end);
    case {'Ur n=4+d','Ur n=6+d','Ur n=8+d'}
        X = [S.lifted(1:sscanf(stateSet,'Ur n=%d'),2:end); S.lifted(1:2,1:end-1)];
    case 'P'
        X = Pp;
    case regexp(stateSet,'^Ur n=\d+\+Ur1\([\d,]+\)$','match','once')
        tok = regexp(stateSet,'^Ur n=(\d+)\+Ur1\(([\d,]+)\)$','tokens','once');
        dmax = max(str2double(split(tok{2},',')));
        Ur1 = [repmat(S.lifted(1,1), 1, dmax), S.lifted(1,:)]; % held before the start
        N = size(S.lifted,2);
        X = [S.lifted(1:str2double(tok{1}),:); zeros(dmax, N)];
        for j = 1:dmax
            X(str2double(tok{1}) + j,:) = Ur1((1:N) + dmax - j); % Ur1(k-j)
        end
        X = X(:,2:end);
    case regexp(stateSet,'^Ur n=\d+\+M\d+fix$','match','once')
        tok = regexp(stateSet,'^Ur n=(\d+)\+M(\d+)','tokens','once');
        [~, M] = movingMeanDev(S.lifted(1:2,:), str2double(tok{2}));
        X = [S.lifted(1:str2double(tok{1}),:); M];
        X = X(:,2:end);
    case regexp(stateSet,'^Ur n=\d+\+D\d+(,M)?$','match','once')
        tok = regexp(stateSet,'^Ur n=(\d+)\+D(\d+)','tokens','once');
        [DUr, M] = movingMeanDev(S.lifted(1:2,:), str2double(tok{2}));
        X = S.lifted(1:str2double(tok{1}),:);
        X = [X; DUr];
        if endsWith(stateSet,',M'), X = [X; M]; end
        X = X(:,2:end);
    case 'P,sqrtP'
        X = [Pp; sqrt(Pp)];
    case 'P,cbrtP'
        X = [Pp; nthroot(Pp,3)];
    case 'P,Ue'
        X = [Pp; nthroot(Pp./ctp,3)];
    case 'P,sqrtP,cbrtP'
        X = [Pp; sqrt(Pp); nthroot(Pp,3)];
end
W = W(:,2:end);
P = S.P(2:end)/1e6;
end

function [DUr, M] = movingMeanDev(Ur, N)
% Deviation of Ur from the moving mean of the past Ur over N samples, with
% the recursion of K_psi (initial mean Ur(:,1), k = 1 at the first sample)
M = zeros(size(Ur));
Mk = Ur(:,1);
UrPrev = Ur(:,1);
for k = 1:size(Ur,2)
    l = min(N,k);
    Mk = (l-1)/l*Mk + UrPrev/l;
    M(:,k) = Mk;
    UrPrev = Ur(:,k);
end
DUr = Ur - M;
end

%% ======================================================================
%% buildLinearPowerModel  (task 'linearPowerModel')
function model = buildLinearPowerModel(modelAff, filenameId, filenameVal, t0)
% buildLinearPowerModel linear Koopman power model of the WRC-AIC from the
% control-affine model (feature set 'NL+1+M+x,aff'): the features become
% states, so that the MPC needs no linearization (MatsysPfeat)
%
%   z(k+1) = A z(k) + B u(k),  P_WF(k) = C z(k) + D u(k)  (MW)
%
% States z (21): lifted Ur (Ur1, Ur2, Ur1^2, Ur2^2, Ur1^3, Ur2^3) and the
% feature states of featStatesLin: filtered inputs s = [ct1; ct2; phi1],
% offset 1, g1, g2, moving means Mc, Mgamma, Ma, Mac, products g2*Ma,
% g2*Mgamma, g2*Mac, induction a1 and ac = a1*cos^2(phi1) (all of s);
% inputs u (3): the controls before the filter of the turbine model.
% Rows of A, B, C, D:
% - lifted Ur and output: A_w, B_w, B_u, C_w, D_w, D_u of the control-affine
%   model (linear in its features and controls, i.e. unchanged)
% - filtered inputs: s(k+1) = (1-tauU) s(k) + tauU u(k) (inputFilterCoeff)
% - offset: 1(k+1) = 1(k)
% - moving means: M(k+1) = (1-1/tau) M(k) + 1/tau [s1; s3; a1; ac](k)
% - g1, g2, products, a1, ac: nonlinear functions of s, i.e. of the inputs;
%   identified by least squares from the identification data (EDMD) with
%   all states and inputs as regressors (the only approximated rows)
%
% Inputs
% - modelAff: control-affine model (KoopFeatModel<k>.mat, 'NL+1+M+x,aff')
% - filenameId, filenameVal: data sets as in testPowerStates (identification
%   from t0 on, the validation set is only needed by testPowerInputsOnly)
% - t0: first sample

assert(strcmp(modelAff.featureSet, 'NL+1+M+x,aff'), 'buildLinearPowerModel: control-affine model expected');
mainDir = fileparts(fileparts(mfilename('fullpath')));
% (path: see WFSim_addpaths)
if iscell(filenameVal), filenameVal = filenameVal{1}; end
evalc('[~, D] = testPowerInputsOnly(filenameId, filenameVal, 8, t0);');
close(findobj('Type','figure','Name','PowerInputsOnly'));
S = D.id;

% State data as in testPowerStates (from the second sample on)
featureSet = 'NL+1+M+x,lin';
tau = modelAff.tau;
W = koopmanLiftedStates('powerFeatures', S.lifted(end-2:end,:), featureSet, tau);
[F, U] = koopmanLiftedStates('featStatesLin', W);
nL = modelAff.n; % lifted Ur
Z = [S.lifted(1:nL,:); F];
Z = Z(:,2:end); U = U(:,2:end);
nz = size(Z,1); nu = size(U,1);

% Indices of the states
iS = nL + (1:3); iOne = nL + 4; iG = nL + [5 6]; iM = nL + (7:10);
iX = nL + (11:13); iA = nL + [14 15];
iNL = [iG, iX, iA]; % rows identified by least squares

A = zeros(nz); B = zeros(nz, nu);
A(1:nL, 1:nL + 13) = [modelAff.A, modelAff.B(:,1:13)]; % lifted Ur: A_w, B_w
B(1:nL,:) = modelAff.B(:,14:16); % B_u
C = [modelAff.C, modelAff.D(1:13), 0, 0]; % C_w, D_w
Dm = modelAff.D(14:16); % D_u
tauU = koopmanLiftedStates('filterCoeff');
A(iS,iS) = (1 - tauU)*eye(3);
B(iS,:) = tauU*eye(3);
A(iOne,iOne) = 1;
A(iM,iM) = (1 - 1/tau)*eye(4);
A(iM,[iS(1), iS(3), iA]) = 1/tau*eye(4); % h = [s1; s3; a1; ac]

% Nonlinear rows: least squares (scaled regressors, small ridge term as
% testPowerStates); regressors: the feature states and the inputs, not the
% lifted Ur, as the features are functions of the filtered inputs and the
% moving means only (the feature states are then driven by u only)
iF = nL+1:nz;
Psi = [Z(iF,1:end-1); U(:,1:end-1)];
Y = Z(iNL,2:end);
sc = sqrt(mean(Psi.^2, 2));
K = (Y*(Psi./sc)')/((Psi./sc)*(Psi./sc)' + 1e-8*size(Psi,2)*eye(size(Psi,1)));
K = K./sc';
A(iNL,iF) = K(:,1:numel(iF));
B(iNL,:) = K(:,numel(iF)+1:end);

% Check: one-step prediction of the states and the output on the data
Zp = A*Z(:,1:end-1) + B*U(:,1:end-1);
errNL = sqrt(mean((Zp(iNL,:) - Z(iNL,2:end)).^2, 2))';
errExact = max(max(abs(Zp(setdiff(nL+1:nz, iNL),:) - Z(setdiff(nL+1:nz, iNL),2:end))));
eigA = eig(A(setdiff(1:nz, iOne), setdiff(1:nz, iOne))); % without the offset (eigenvalue 1)
fprintf('Linear power model: %d states, %d inputs, max|eig(A)| without offset %.4f\n', nz, nu, max(abs(eigA)));
fprintf('  one-step RMSE of the identified rows (g1 g2 g2Ma g2Mg g2Mac a1 ac): %s\n', sprintf('%.2e ', errNL));
fprintf('  max one-step error of the known filter rows: %.2e\n', errExact);

model = struct('A',A,'B',B,'C',C,'D',Dm,'n',nL,'nz',nz,'stateSet',modelAff.stateSet, ...
    'featureSet',featureSet,'tau',tau,'maxAbsEig',max(abs(eigA)), ...
    't0',t0,'filenameId',{filenameId},'modelAff',modelAff);
end

%% ======================================================================
%% identifyAICnoV1  (task 'aicNoV1')
function model = identifyAICnoV1(filenameId, t0, fileOut)
% identifyAICnoV1 Koopman model of the effective wind speeds for the
% adaptive AIC without the measured wind speed V1 (thesis Section
% "Koopman Matrix Creation" with the thrust control signals only,
% w = [CT1; CT2] instead of [V1; CT1; CT2], no yaw):
% phi(k+1) = A phi(k) + B w(k), phi = [Ur1 Ur2 Ur1^2 Ur2^2 Ur1^3 Ur2^3].
% Saved as DataInOutWfSim/AdaptiveAIC/<fileOut> with the correlation
% matrices L = G+ Gw', Lu = Gw Gw' (sums over the samples, as in the thesis)
% as initial values of the online update.
%
% Inputs
% - filenameId: open-loop data (default: random CT levels at V_inf = 7 m/s)
% - t0: first sample (s), default 2: the start-up of the wake (Vr1 = Vr2)
%   is needed, the closed loop starts from the undisturbed flow
% - fileOut: file name of the model (default KoopAICnoV1.mat), e.g.
%   KoopAICnoV1_Vinf9.mat for the offline model at 9 m/s

if nargin < 1 || isempty(filenameId)
    filenameId = 'Vinf7dot0_sowfa_2turb_ct_levels';
end
if nargin < 2 || isempty(t0), t0 = 2; end
if nargin < 3 || isempty(fileOut), fileOut = 'KoopAICnoV1.mat'; end
mainDir = fileparts(fileparts(mfilename('fullpath')));
dirData = [extractBefore(filenameId,'_sowfa'),'_OL_Ct']; % folder of WFSim_demo, e.g. Vinf7dot0_OL_Ct
S = load(fullfile(mainDir,'DataT2OLWFSim',dirData,filenameId),'Ur1','Ur2','Ct1','Ct2');
k = t0:numel(S.Ur1);
Ur = [S.Ur1(k); S.Ur2(k)];
W = [S.Ct1(k); S.Ct2(k)];

G = [Ur; Ur.^2; Ur.^3];
Gw = [G(:,1:end-1); W(:,1:end-1)];
Gp = G(:,2:end);
N = size(Gw,2);
L = Gp*Gw';
Lu = Gw*Gw';
Kmat = L/Lu;
n = size(G,1);
model = struct('A',Kmat(:,1:n),'B',Kmat(:,n+1:end),'L',L,'Lu',Lu,'n',n, ...
    'inputs',{{'CT1','CT2'}},'filenameId',filenameId,'t0',t0);
model.maxAbsEig = max(abs(eig(model.A)));

dirOut = fullfile(mainDir,'DataInOutWfSim','AdaptiveAIC');
if ~isfolder(dirOut), mkdir(dirOut); end
save(fullfile(dirOut,fileOut),'model');
fprintf('%s: %d samples, max|eig(A)| = %.4f\n', fileOut, N, model.maxAbsEig);
end

%% ======================================================================
%% simAdaptiveKoopman  (task 'simAdaptive')
function [UrEst, upd, err, modelEnd, modelsAt] = simAdaptiveKoopman(model, Ur, W, adaptive, ad, saveAt)
% simAdaptiveKoopman open-loop estimation of the effective wind speeds with
% the Koopman model as in the MPC (wfmodelNL): the model is simulated from
% the measured initial state, the predicted Ur are lifted again in each
% step (K_psi). With adaptive = true, the model is updated online with the
% measured Ur (koopmanAdaptStep), as in the controller.
%
% Inputs
% - model: struct with A, B, L, Lu (identifyAICnoV1)
% - Ur: measured effective wind speeds (2 x N)
% - W: inputs (nw x N), here [CT1; CT2]
% - adaptive: online update on/off
% - ad: settings of the update (h, nu, maxEig), default h = 5, nu = 0.15 m/s;
%   ad.resetState = true: after an update the estimator continues from the
%   measured Ur (default false: free run as in wfmodelNL)
% - saveAt: samples after which the model is returned in modelsAt (optional)
%
% Outputs
% - UrEst: estimated effective wind speeds (2 x N)
% - upd: update flags (1 x N)
% - err: estimation error of the last h samples (1 x N)
% - modelEnd: model after the last sample
% - modelsAt: models (fields A, B, k) after the samples saveAt

if nargin < 5, ad = []; end
if nargin < 6, saveAt = []; end
modelsAt = struct('A',{},'B',{},'k',{});
if ismember(1, saveAt) % initial model
    modelsAt(1) = struct('A',model.A,'B',model.B,'k',1);
end
resetState = isstruct(ad) && isfield(ad,'resetState') && ad.resetState;
N = size(Ur,2);
UrEst = nan(2,N); UrEst(:,1) = Ur(:,1);
upd = zeros(1,N); err = zeros(1,N);
x = lift(Ur(:,1));
for k = 2:N
    xn = model.A*x + model.B*W(:,k-1);
    UrEst(:,k) = xn(1:2);
    x = lift(xn(1:2));
    if adaptive
        [model, ad, upd(k)] = koopmanAdaptStep(model, ad, Ur(:,k), W(:,k), UrEst(:,k));
        err(k) = ad.err;
        if upd(k) && resetState
            x = lift(Ur(:,k));
        end
    end
    if ismember(k, saveAt)
        modelsAt(end+1) = struct('A',model.A,'B',model.B,'k',k); %#ok<AGROW>
    end
end
modelEnd = model;
end

function phi = lift(Ur)
phi = [Ur; Ur.^2; Ur.^3];
end

%% ======================================================================
%% aicNoV1ModelTimes  (task 'aicNoV1ModelTimes')
function sel = aicNoV1ModelTimes()
% aicNoV1ModelTimes times of the models of the adaptive AIC without V1 that
% are shown with their singular values (plotAICnoV1Sigma) and marked in the
% update subplot of the adaptive open-loop validation (validateAICnoV1):
% the initial model at t0 = 240 s and the updated models at 5100, 6000 and
% 8000 s (ramp from 7 to 9 m/s at 5000 s, last update at 6235 s).
% Colors: initial model gray, the others in colors not used for WT1, WT2
% (lines 1, 2).
%
% Output
% - sel: table with Time_s, Color and TextSide (time stamp in the update
%   subplot right (> 0) or left (< 0) of the line, the magnitude scales the
%   distance; 5100 s left, next to the first updates, 6000 s right of
%   the last update at 6235 s)

cl = lines(5);
sel = table([240; 5100; 6000; 8000], [0.5 0.5 0.5; cl(3:5,:)], [1; -1; 4; 1], ...
    'VariableNames', {'Time_s','Color','TextSide'});
end

%% ======================================================================
%% eDMD_RT_Vin
function [sys_red,FITje,Xd,Xd_p,xsim,FITje_val,xo,Koop,ysim,Inputs,Deterministic,ysim_val,Inputs_val,Deterministic_val] ...
    = eDMD_RT_Vin(states,stateName,statesvalid,poly,n,yawmode)
%[sys_red,FITje,Xd,Xd_p,xsim,FITje_val,fig1,xo,Koop]
%[sys_red,FITje,U,S,V,X,X_p,x]
% eDMD_RTUpdate builds a state space model from the states
% and input gathered in the simulation 

if nargin < 6
    yawmode = 0;
end

% Size states, inputs, outputs
nPsi = size(states,1); % number of lifted states + inputs
nx = n; % number of lifted states: Ur1,Ur2 plus lifted states
nu = nPsi - n; %size(inp,1); % number inputs: ct1, ct2,v1 
ny = 2; % number outputs:  Ur1,Ur2 (no lifted states)

t0 = 1;
tend = round(size(states,2));
tendVal = size(statesvalid,2);

psi_xk_1 = states(:,t0:tend-1); % States delayed
psi_xk = states(:,t0+1:tend); % States  states(1:end-1,t0:tend-1); % %X_k = states(1:nx,t0:tend-1);

Inputs = psi_xk(end-(nu-1): end,:);
Deterministic = psi_xk(1:ny,:); 

Inputs_val = statesvalid((nPsi -(nu-1): nPsi), (t0:tendVal-1));
Deterministic_val = statesvalid(1:ny, (t0+1) : tendVal); 

% Psi_k = K * Psi_k_1 -> K = (Psi_k*Psi_k_1^T) * pinv(Psi_k_1*Psi_k_1^T)  =
%     K = A * pinv(G)
% Psi_k'  = Psi_k_1'* K' -> K' = pinv(psi_xk_1 * psi_xk_1') * psi_xk_1 * psi_xk'
%     K' = pinv(G) * A
Koop.G = psi_xk_1 * psi_xk_1'; % Update G (GT
Koop.A = psi_xk_1 * psi_xk';   % Update A (At)
Koop.K = pinv(Koop.G) * Koop.A;% Update Koopman operator
K = Koop.K'; Ky = K(1:ny,:);

approxA = K(1:nx,1:nx); % system matrix
approxB = K(1:nx,nx+1:end);
approxC = [eye(ny),zeros(ny,nx-ny)];
approxD = zeros(ny,nu);
sys_red = ss(approxA,approxB,approxC,approxD,1);

% sysC = d2c(sys_red,'tustin');
% max(real(eig(sysC)))

xo = states(1:nx,t0);% for States Ur1, Ur2 and lifted states
% [ysim,~,xsim] = lsim(sys_red, Inputs',[],xo);

structPrev.Ur1_prev1 = psi_xk_1(1,1);
structPrev.Ur2_prev1 = psi_xk_1(2,1);
structPrev.dUr1_prev1 = 0;
structPrev.dUr2_prev1 = 0;
structPrev.M1(1) = structPrev.Ur1_prev1; % Moving mean
structPrev.M2(1) = structPrev.Ur2_prev1;
structPrev.k = 1;

xsim = nan(size(psi_xk_1));
xsim(1:ny,1) = Ky * psi_xk_1(:,1); % states: 2 + 6 + 3: P1, P2, Ur1, Ur2, n
xsim(:,1) = koopmanLiftedStates('urStates', [xsim(1:ny,1);psi_xk(nx+1:end,1)],poly,n,structPrev,yawmode);
for idx = 1: length(psi_xk_1)-1
    xsim(1:ny,idx+1) = Ky * xsim(:,idx);
    [xsim(:,idx+1),~,structPrev] = koopmanLiftedStates('urStates', [xsim(1:ny,idx+1);psi_xk_1(nx+1:end,idx+1)], poly,n,structPrev,yawmode);
end
ysim = xsim(1:ny,:);

FITje = vaf(Deterministic,ysim);

% xo_val = statesvalid(1:nx,t0);%for States Ur1,Ur2 and lifted states
% ysim_val = lsim(sys_red, Inputs_val',[],xo_val);
structPrev.Ur1_prev1 = statesvalid(1,1);
structPrev.Ur2_prev1 = statesvalid(2,1);
structPrev.dUr1_prev1 = 0;
structPrev.dUr2_prev1 = 0;
structPrev.M1(1) = structPrev.Ur1_prev1; % Moving mean
structPrev.M2(1) = structPrev.Ur2_prev1;
structPrev.k = 1;

xsim_val = nan(size(statesvalid));
xsim_val(1:ny,1) = Ky * statesvalid(:,1);
xsim_val(:,1) = koopmanLiftedStates('urStates', [xsim_val(1:ny,1);statesvalid(nx+1:end,1)],poly,n,structPrev);
for idx = 1: length(xsim_val)-1
    xsim_val(1:ny,idx+1) = Ky * xsim_val(:,idx);
    [xsim_val(:,idx+1),~,structPrev]  = koopmanLiftedStates('urStates', [xsim_val(1:ny,idx+1);statesvalid(nx+1:end,idx+1)],poly,n,structPrev,yawmode);
end
ysim_val = xsim_val(1:ny,2:end);

FITje_val = vaf(Deterministic_val,ysim_val);

tmp = strsplit(stateName,';');
sys_red.StateName = tmp(1:nx); % state names
sys_red.OutputName = {'Ur1';'Ur2'};

if max(size(sys_red)) == 4
    sys_red.InputName  =  {'Ct1';'Ct2';'phi1';'Vin'};
elseif max(size(sys_red)) == 3 && yawmode
    sys_red.InputName  =  {'Ct1';'Ct2';'phi1'};
elseif max(size(sys_red)) == 3
    sys_red.InputName  =  {'Ct1';'Ct2';'Vin'};
else
    sys_red.InputName  =  {'Ct1';'Ct2'};
end

% Ouput names should be consistent
Xd = psi_xk_1; Xd_p = psi_xk;

% Plot for validation data
% Vin = 1;
% fig1 = koopmanPlots('effWind', ysim_val,Inputs_val,Deterministic_val, FITje_val,dirFig,Vin,nx);
% strVal = 'Id.';
% koopmanPlots('effWind', ysim,Inputs,Deterministic, FITje,dirFig,Vin,nx,strVal);


end

%% ======================================================================
%% eDMD_RT_Vin_P
function [sys_red,FITje,Xd,Xd_p,xsim,FITje_val,xo,Koop,ysim,Inputs,Deterministic,ysim_val,Inputs_val,Deterministic_val] ...
    = eDMD_RT_Vin_P(states,stateName,statesvalid,poly,n,yawmode)
%[sys_red,FITje,Xd,Xd_p,xsim,FITje_val,fig1,xo,Koop]
%[sys_red,FITje,U,S,V,X,X_p,x]
% eDMD_RTUpdate builds a state space model from the states
% and input gathered in the simulation 

if nargin < 6
    yawmode = 0;
end
 
% Size states, inputs, outputs
nx = 2; % number of original states: Ur1,Ur2 
ny = 1; % number outputs:  sum(P1,P2 (no lifted states)
nxl = n; % number of original and lifted states: Ur1,Ur2 plus lifted states
nxle = nxl + ny; % index last state in vector
nPsi = size(states,1); % number of lifted states + outputs + inputs
nu = nPsi - nxle; % number inputs: ct1, ct2,v1 

t0 = 1; % Start point; This is just for testing: 
tend = size(states,2);
tendVal = size(statesvalid,2);

% Identification set
psi_xk_1 = states(ny+1:end,t0:(tend-1)); % delayed states and inputs [x_k-1;u_k-1]
out_xk_1 = states(1:ny,t0:(tend-1)); % delayed outputs [x_k-1;u_k-1]
X_xk = states(ny+1:nxle,(t0+1):tend); % States x_k 
psi_xk = [X_xk; out_xk_1]; % [x_k+1;y_k] % States and delayed outputs 
Inputs = psi_xk_1(nxl+1:end,:); % delayed inputs [x_k-1;u_k-1]
Deterministic = out_xk_1; % outputs (ToDo Should this line be removed?)

% Validation set
psi_xk_1_val = statesvalid(ny+1:end,t0:(tendVal-1)); % delayed states and inputs [x_k-1;u_k-1]
out_xk_1_val = statesvalid(1:ny,t0:(tendVal-1)); % delayed outputs [x_k-1;u_k-1]
% X_xk_val = statesvalid(ny+1:nxle,(t0+1):tendVal); % States x_k 
% psi_xk_val = [X_xk_val; out_xk_1_val]; % [x_k+1;y_k] % States and delayed outputs 
Inputs_val = psi_xk_1_val(nxl+1:end,:); % delayed inputs [x_k-1;u_k-1]
Deterministic_val = out_xk_1_val; 

% Psi_k = K * Psi_k_1 -> K = (Psi_k*Psi_k_1^T) * pinv(Psi_k_1*Psi_k_1^T)  =
%     K = A * pinv(G)
% Psi_k'  = Psi_k_1'* K' -> K' = pinv(psi_xk_1 * psi_xk_1') * psi_xk_1 * psi_xk'
%     K' = pinv(G) * A
Koop.G = psi_xk_1 * psi_xk_1'; % Update G (GT
Koop.A = psi_xk_1 * psi_xk';   % Update A (At)
Koop.K = pinv(Koop.G) * Koop.A;% Update Koopman operator
K = Koop.K'; K_P = K((nxl+1):end,:);
K_V = K(1:nx,:);

approxA = K(1:nxl,1:nxl); % system matrix
approxB = K(1:nxl,nxl+1:end);
approxC = K(nxl+1:end,1:nxl);  %[eye(ny),zeros(ny,nx-ny)];
approxD = K(nxl+1:end,nxl+1:end); %zeros(ny,nu);
sys_red = ss(approxA,approxB,approxC,approxD,1);

% sysC = d2c(sys_red,'tustin');
% max(real(eig(sysC)))

xo = states(1:nxle,t0);% for States Ur1, Ur2 and lifted states
[ysim1,~,xsim1] = lsim(sys_red, Inputs',[],xo(ny+1:end)); %FITjeLin = vaf(Deterministic(1:end-1),ysim1(1:end-1));

structPrev.Ur1_prev1 = psi_xk_1(1,1);
structPrev.Ur2_prev1 = psi_xk_1(2,1);
structPrev.dUr1_prev1 = 0;
structPrev.dUr2_prev1 = 0;
structPrev.M1(1) = structPrev.Ur1_prev1; % Moving mean
structPrev.M2(1) = structPrev.Ur2_prev1;
structPrev.k = 1;

xsim = nan(size(psi_xk_1));
ysim = nan(ny,length(psi_xk_1)-1);
ysim(:,1) = K_P * psi_xk_1(:,1);
xsim(1:nx,1) = K_V * psi_xk_1(:,1);
xsim(:,1) = koopmanLiftedStates('urStates', [xsim(1:nx,1);Inputs(:,1)],poly,n,structPrev,yawmode);

for idx = 1: length(psi_xk_1)-1
    ysim(:,idx) = K_P * xsim(:,idx);
    xsim(1:nx,idx+1) = K_V * xsim(:,idx);
    [xsim(:,idx+1),~,structPrev] = koopmanLiftedStates('urStates', [xsim(1:nx,idx+1); Inputs(:,idx+1)], poly,n,structPrev,yawmode);
end
FITje = vaf(Deterministic(1:end-1),ysim);

% xo_val = statesvalid(1:nx,t0);%for States Ur1,Ur2 and lifted states
% ysim_val = lsim(sys_red, Inputs_val',[],xo_val);

structPrev.Ur1_prev1 = psi_xk_1_val(1,1);
structPrev.Ur2_prev1 = psi_xk_1_val(2,1);
structPrev.dUr1_prev1 = 0;
structPrev.dUr2_prev1 = 0;
structPrev.M1(1) = structPrev.Ur1_prev1; % Moving mean
structPrev.M2(1) = structPrev.Ur2_prev1;
structPrev.k = 1;

xsim_val = nan(size(psi_xk_1_val));
ysim_val = nan(ny,length(psi_xk_1_val)-1);
ysim_val(:,1) =  K_P * psi_xk_1_val(:,1);
xsim_val(1:nx,1) = K_V * psi_xk_1_val(:,1);
xsim_val(:,1) = koopmanLiftedStates('urStates', [xsim_val(1:nx,1);Inputs_val(:,1)],poly,n,structPrev,yawmode);

for idx = 1: length(xsim_val)-1
    ysim_val(:,idx) = K_P * xsim_val(:,idx);
    xsim_val(1:nx,idx+1) = K_V * xsim_val(:,idx);
    [xsim_val(:,idx+1),~,structPrev] = koopmanLiftedStates('urStates', [xsim_val(1:nx,idx+1); Inputs_val(:,idx+1)], poly,n,structPrev,yawmode);
end

FITje_val = vaf(Deterministic_val(:,1:end-1),ysim_val);

tmp = strsplit(stateName,';');
sys_red.StateName = tmp(ny+1:nxle); % state names
sys_red.OutputName = tmp(1:ny);
%sys_red.InputName = tmp(nxle+1: nxle + nu);

if nu == 4
    sys_red.InputName  =  {'Ct1';'Ct2';'phi1';'Vin'};
elseif nu == 3 && yawmode
    sys_red.InputName  =  {'Ct1';'Ct2';'phi1'};
elseif nu == 3
    sys_red.InputName  =  {'Ct1';'Ct2';'Vin'};
else
    sys_red.InputName  =  {'Ct1';'Ct2'};
end

% Ouput names should be consistent
Xd = psi_xk_1; Xd_p = psi_xk;

% Plot for validation data
% Vin = 1;
% fig1 = koopmanPlots('effWind', ysim_val,Inputs_val,Deterministic_val, FITje_val,dirFig,Vin,nx);
% strVal = 'Id.';
% koopmanPlots('effWind', ysim,Inputs,Deterministic, FITje,dirFig,Vin,nx,strVal);


end
