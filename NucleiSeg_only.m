% This script deals with only segmentation protocols for focal adhesion
% structures. It optimizes parameters if desired (optimize_params = 'y')
% and then generates focal adhesion masks based on these parameters.

rehash
if isempty(file_search('fa_\w+.TIF',folder))
    mkdir(folder,'Nuclei Images')
    nuclei_gen_nowater([prefix exp_name '\w+.TIF'],blob_params,sizemin,sizemax,folder);
end
addpath(fullfile(folder,'Nuclei Images'))