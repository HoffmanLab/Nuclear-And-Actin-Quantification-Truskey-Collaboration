% This script quantifies actin intensity and orientations across images

rehash

if ~exist(fullfile(folder,'Sharpness Images'),'dir')
    mkdir(folder,'Sharpness Images')
    imgnames = file_search([prefix exp_name '\w+.TIF'],folder);
    for i = 1:length(imgnames)
        DDcall(imgnames{i},10,folder);
    end
end
addpath(fullfile(folder,'Sharpness Images'))