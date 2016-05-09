function preprocess_stitched(folder,exp_name)
% This function allows one to preprocess a set of images (for now this is
% only background subtraction (automated)

% Brief overview of steps
% (1) Calculation of blur factor (directional derivative)
% (2) Automatic background subtraction

imgNames = file_search([exp_name '\w+.TIF'],folder);

for i = 1:length(imgNames)    
    img = single(imread(fullfile(folder,imgNames{i})));
    params.bin = 1;
    params.nozero = 0;
    img = bs_ff(img,params);
    imwrite2tif(img,[],fullfile(folder,'Preprocessed Images',['pre_' imgNames{i}]),'single')
end

end