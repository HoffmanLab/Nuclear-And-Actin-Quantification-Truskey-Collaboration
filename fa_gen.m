function fa_gen(fname,params,fold)
% fa_gen('bsa_pre_VinTS\w+Venus.TIF,[25 500 50],'VinTS 093015')

% A simple program to generate and save masks using the water program.
% A typical set of params is [25,500,50].

files = file_search(fname,fold);

for i = 1:length(files)
    img = double(imread(files{i}));
    waterFun = @(block_struct) water(block_struct.data,params); %% Use cell finder function rather than water (set boarder pixels to be zero, then 8 bit connectivity to generate bwlabel fa image)
%     img(:,end-68:end) = [];
%     img(end-68:end,:) = [];
    block_size = [50 50];
    block_fas = blockproc(img,block_size,waterFun,'BorderSize',[10,10]);    
%     w = water(img,params);
    imwrite2tif(block_fas,[],fullfile(fold,'FA Images',['fac_' files{i}]),'single');
end

end