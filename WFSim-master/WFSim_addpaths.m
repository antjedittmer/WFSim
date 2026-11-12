function OutDir = WFSim_addpaths
% WFSim_addpaths adds necessary paths for WFSim
WFSimFolder = fileparts(mfilename('fullpath'));% [WFSimFolder, ~, ~] = fileparts(which([mfilename '.m']));         % Get WFSim directory

FastFarmFolder = fileparts(fileparts(fileparts(WFSimFolder)));
ProjectFolder = fileparts(FastFarmFolder);

addpath(genpath(fullfile(ProjectFolder,'YALMIP-master')))

addpath(fullfile(WFSimFolder,'layoutDefinitions')); % Folder with predefined wind farm layouts
addpath(fullfile(WFSimFolder,'controlDefinitions')); % Make use of a predefined timeseries of control inputs
addpath(fullfile(WFSimFolder,'solverDefinitions')); % Folder with model options, solver settings, etc.

addpath(fullfile(WFSimFolder,'bin' ,'core'));      % Add core files
addpath(fullfile(WFSimFolder,'bin', 'analysis'));  % Add analysis files
addpath(fullfile(WFSimFolder,'bin', 'mpc'));  % Add analysis files
addpath(fullfile(WFSimFolder,'bin', 'mpc','mpcFilesWF'));  % Add analysis files


addpath(genpath(fullfile(WFSimFolder,'libraries'))); % Add external libraries
% addpath(genpath(fullfile(WFSimFolder,'YALMIP-master')));

% Output folder
OutDir = fullfile(WFSimFolder,'OutDir'); % Create path to OutDir
if ~isfolder(OutDir), mkdir(OutDir), end % Create OutDir, if it doesn't exist
