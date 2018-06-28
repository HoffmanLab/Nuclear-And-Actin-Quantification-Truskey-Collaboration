function nuclei_gen_nowater(img_files,parray,sizemin,sizemax,folder)
% Generates cell masks automatically based on user-adjusted threshold
% (through ThreshSelect.m). Always requires the user to use ThreshSelect.m
% to optimize threshold parameter.
imgs = file_search(img_files,folder);
high_pass_filt_width = parray(1);
thresh = parray(2);

%% Image Filtering and Prep. for Analysis
for i = 1:length(imgs)
    imOrig = imread(imgs{i});
    [r,c] = size(imOrig);
    AvgFilt = fspecial('average',high_pass_filt_width);
    pImg = padarray(imOrig,[high_pass_filt_width high_pass_filt_width],'symmetric');
    pSmImg = conv2(pImg,AvgFilt,'same');
    SmImg = pSmImg(1+high_pass_filt_width:r+high_pass_filt_width,1+high_pass_filt_width:c+high_pass_filt_width);
    FiltImg = double(imOrig) - SmImg;
    im = ones(size(imOrig));
    bw = im.*(FiltImg > thresh);
    
    % Dilate then erode once to smooth
    SE = strel('disk',3);
    bw = imerode(bw,SE);
    bw = imdilate(bw,SE);
    
    D = -bwdist(~bw);
    mask = imextendedmin(D,0.4);
    D2 = imimposemin(D,mask); % Watershed from only certain seed regions
    Ld2 = watershed(D2); % Find watershed boundaries with these seed constraints
    bw(Ld2 == 0) = 0; % Add boundaries between cells
    nuclei = imfill(bw,'holes'); % Fill holes after watershedfigure
%     nuclei_filt = bwareaopen(nuclei,sizemin);
%     nuclei_filt = bwpropclose(nuclei_filt,'Area',sizemax);
%     nuclei_filt = bwpropopen(nuclei_filt,'Solidity',0.9);
    % Label each blob with 8-connectivity, so we can make measurements of it
    % and get rid of small cells, save out masks
    nuclei = bwlabel(nuclei, 8);
%     nuclei_filt = bwlabel(nuclei_filt, 8);
    % Write out tif of cells (grayscale, not rgb)
    imwrite2tif(nuclei,[],fullfile(folder,'Nuclei Images',['nuclei_' imgs{i}]),'single')
%     imwrite2tif(nuclei_filt,[],fullfile(folder,'Nuclei Images',['nuclei_filt_' imgs{i}]),'single')
end
end

