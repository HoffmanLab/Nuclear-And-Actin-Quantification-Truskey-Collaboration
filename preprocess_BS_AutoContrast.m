function preprocess_BS_AutoContrast(exp_name,folder,varargin)
% This function allows one to perform background subtraction on DAPI images

% Overview of steps
% (1) Background subtraction

n = nargin;
imgNames = file_search([exp_name '\w+.TIF'],folder);
for i = 1:length(imgNames)
    img = single(imread(imgNames{i}));
    %Background subtraction
    params.bin = 1;
    params.nozero = 0;
    if n ==2
        img = bs_ff_mod(img,params);
    else
        img = bs_ff_mod(img,varargin{1}(1),params);
    end
    img = img./(2^12);
    
   
    
    % Option 1
    adjustFun = @(block_struct) imadjust(block_struct.data); %% Use cell finder function rather than water (set boarder pixels to be zero, then 8 bit connectivity to generate bwlabel fa image)
    block_size = [150 150];
    img = blockproc(img,block_size,adjustFun,'BorderSize',[50,50]); 
    
%     Option 2
%     img = imadjust(img);


    img = img.*(2^12);
    %Write out as 32bit TIFs
    imwrite2tif(img,[],fullfile(folder,'Preprocessed Images',['pre_' imgNames{i}]),'single')
end

end