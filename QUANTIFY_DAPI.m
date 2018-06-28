function QUANTIFY_DAPI(varargin)

%% Header
% This is a main function for running Nuclear Quantification Code on a single
% experimental group. The Exp_Params text file in the folder containing the data
% should describe all relevant information for a particular experimental design.
% Sample Exp_params text files can be found in the GitHub repository.

%% Parameter Processing
% Check to see if anything has been passed as parameter, if anything has
% been passed, make sure it is a folder.

if (not(isempty(varargin)))
    if (exist(varargin{1},'dir'))
        folder = varargin{1};
    else
        error('Expected first parameter to be a folder with images to process.');
    end
else
    %% Set up
    clear;
    close all;
    clc;
end

%% Read in Pre-processing parameters

if (exist('folder','var'))
    [~,params_file] = GetParamsFile(folder); %#ok<ASGLU>
else
    [folder,params_file] = GetParamsFile; %#ok<ASGLU>
end
ProcessParamsFile;

%% Manual Background Subtraction and Auto Contrasting
addpath(folder);
rehash
PreProcessDAPI_only

%% Segmentation
rehash
if strcmpi(segmentation,'y')
    if strcmpi(structure,'Nuclei')
        NucleiSeg_only
    end
end

%% Blob Analysis
rehash
if strcmpi(segmentation,'y') && strcmpi(banalyze,'y') 
    if strcmpi(structure,'Nuclei')
        NucleiAnalyze_only
    end
end