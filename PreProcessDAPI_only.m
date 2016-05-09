% This script performs pre-processing functions for stitched DAPI image
% only, and includes scripts to unstitch images, filter based on clarity
% directional derivative, and background subtract them.

%% Preprocess Images
rehash
if isempty(file_search('pre_\w+.TIF',folder))
    mkdir(fullfile(folder,'Preprocessed Images'))
    preprocess_stitched(folder,exp_name)
end
addpath(fullfile(folder,'Preprocessed Images'))