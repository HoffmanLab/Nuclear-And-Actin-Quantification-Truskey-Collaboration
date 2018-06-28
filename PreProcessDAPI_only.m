% This script performs pre-processing functions for Nuclei Quantification, 
% which only includes background subtraction for now

% If you want to subtract a particular value off of each image as
% background subtraction, make sure you have a variable called "bvals" that
% has your chosen background value for each imaging channel.


%% Preprocess images using PreParams.mat file in GoogleDrive (Protocols -> Analysis Protocols -> FRET)
rehash
if isempty(file_search('pre_\w+.TIF',folder))
    mkdir(fullfile(folder,'Preprocessed Images'))
    if exist('bvals','var')
        preprocess_BS_AutoContrast(exp_name,folder,bvals)
    else
        preprocess_BS_AutoContrast(exp_name,folder)
    end
end
addpath(fullfile(folder,'Preprocessed Images'))