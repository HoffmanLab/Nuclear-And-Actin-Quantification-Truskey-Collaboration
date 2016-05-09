function DDcall(imgname,gaussian_kernel_size,folder)

%Script designed to compute structure tensor over an image in a
%user-input region size
% SP.actin_optimize = 'y';
% save(fullfile(pwd,SP.folder,['SaveParams_' SP.folder '.mat']),'-struct','SP');

img = double(imread(fullfile(folder,'Preprocessed Images',imgname)));
[x,y] = size(img);

% % Initialize variables for speed
% e1 = cell(x,y);
% e2 = cell(x,y);
% l1 = zeros(x,y);
% l2 = zeros(x,y);
% theta1 = zeros(x,y);
% theta2 = zeros(x,y);

DoG = difference_of_gaussian_kernels(gaussian_kernel_size);
Ix = conv2(img, DoG.Gx,'same');
Iy = conv2(img, DoG.Gy,'same');
DD_mag = sqrt(Ix.^2 + Iy.^2);
DD_phase = atan(Ix./Iy);
DD_phase = DD_phase + pi/2;
DD_phase = DD_phase.*(180/pi);

% %% Brad's version (may need debugging)
% t = Ix.*Ix;
% u = Ix.*Iy;
% v = Iy.*Iy;
% 
% m = reshape([t(:)';u(:)';u(:)';v(:)'],2,2,x*y);
% c = num2cell(m,[1,2]);
% Ac2 = reshape(c,x,y);
% [e3,e4,l3,l4] = cellfun(@eigen_decomposition,Ac2,'un',0);
% l3 = cell2mat(l3);
% l4 = cell2mat(l4);
% 
% first_e4 = cellfun(@(v) v(1), e4);
% second_e4 = cellfun(@(v) v(2), e4);
% theta4 = atan(second_e4./first_e4);
% theta4 = (theta4+pi/2).*(180/pi);

imwrite2tif(DD_mag,[],fullfile(folder,'Actin Quantification Images',['DDmag_' imgname(1:end-4) '.TIF']),'single')
imwrite2tif(DD_phase,[],fullfile(folder,'Actin Quantification Images',['DDphase_' imgname(1:end-4) '.TIF']),'single')
% imwrite2tif(l3,[],fullfile(folder,'Actin Quantification Images',['STl1_' imgname(1:end-4) '.TIF']),'single')
% imwrite2tif(theta4,[],fullfile(folder,'Actin Quantification Images',['STtheta2_' imgname(1:end-4) '.TIF']),'single')


% %% My version
% for i = 1:x
%     for j = 1:y
%         IxI = Ix(i,j);
%         IyI = Iy(i,j);
%         MI = [IxI*IxI, IxI*IyI;...
%             IxI*IyI, IyI*IyI];
%         % From structureTensorDemo
%         [e1{i,j},e2{i,j},l1(i,j),l2(i,j)] = eigen_decomposition(MI);
%         theta1(i,j) = atan(e1{i,j}(2)/e1{i,j}(1));
%         theta2(i,j) = atan(e2{i,j}(2)/e2{i,j}(1));
%     end
% end
% 
% theta1 = (theta1+pi/2).*(180/pi);
% theta2 = (theta2+pi/2).*(180/pi);

% imwrite2tif(DD_mag,[],fullfile(folder,'Actin Quantification Images',['DDmag_' imgname(1:end-4) '.TIF']),'single')
% imwrite2tif(DD_phase,[],fullfile(folder,'Actin Quantification Images',['DDphase_' imgname(1:end-4) '.TIF']),'single')
% imwrite2tif(l1,[],fullfile(folder,'Actin Quantification Images',['STl1_' imgname(1:end-4) '.TIF']),'single')
% % imwrite2tif(l2,[],fullfile(folder,'Actin Quantification Images',['STl2_' imgname(1:end-4) '.TIF']),'single')
% % imwrite2tif(theta1,[],fullfile(folder,'Actin Quantification Images',['STtheta1_' imgname(1:end-4) '.TIF']),'single')
% imwrite2tif(theta2,[],fullfile(folder,'Actin Quantification Images',['STtheta2_' imgname(1:end-4) '.TIF']),'single')
function DoG = difference_of_gaussian_kernels(maskSize)

% difference_of_gaussian_kernels - output first to fourth order diff. of Gaussian kernels (used for partial derivative calc) %%%%%%%%%%%%%%%%%%%%%%%%%
%     DoG = difference_of_gaussian_kernels(maskSize)
% 
%     Given the size of the kernel output the appropriate difference of Gaussians from first to fourth order derivative
%     Based on Salden et al., 1992
% 
% Example:
%  DoG = difference_of_gaussian_kernels(11); %----- creates the DoG structure of scale 11x11
% 
%     Author: Shawn Arseneau
%     Created: September 15, 2006
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    gmap = double(fspecial('gaussian', maskSize));
   
    gLeft  = [gmap, zeros(maskSize,1)];
    gRight = [zeros(maskSize,1), gmap];
    gDiff = gLeft-gRight;
    DoG.Gx = gDiff(:,1:maskSize);

    gTop    = [gmap; zeros(1,maskSize)];
    gBottom = [zeros(1,maskSize); gmap];
    gDiff = gTop-gBottom;
    DoG.Gy = gDiff(1:maskSize,:);
    
    DoG.Gx2 = conv2(DoG.Gx,DoG.Gx,'same');
    DoG.Gy2 = conv2(DoG.Gy,DoG.Gy,'same');
    DoG.Gxy = conv2(DoG.Gx,DoG.Gy,'same');
    
    DoG.Gx3 = conv2(DoG.Gx2,DoG.Gx,'same');
    DoG.Gx2y = conv2(DoG.Gx2,DoG.Gy,'same');
    DoG.Gxy2 = conv2(DoG.Gx,DoG.Gy2,'same');
    DoG.Gy3 = conv2(DoG.Gy2,DoG.Gy,'same');
    
    DoG.Gx4 = conv2(DoG.Gx3,DoG.Gx,'same');
    DoG.Gx3y = conv2(DoG.Gx3,DoG.Gy,'same');
    DoG.Gx2y2 = conv2(DoG.Gx2,DoG.Gy2,'same');
    DoG.Gxy3 = conv2(DoG.Gx,DoG.Gy3,'same');
    DoG.Gy4 = conv2(DoG.Gy3,DoG.Gy,'same');
end

function [e1,e2,e3,l1,l2,l3] = eigen_decomposition(M)

% eigen_decomposition = uses eigs.m and appropriately labels output (for both 2D and 3D)  %%%%%%%%%%%%%%
%
%   [e1,e2,l1,l2] = eigen_decomposition(M)
%
%   Perform eigen-based decomposition of the input tensor M.
%   
%   INPUT:
%    M = structure tensor of the form [dxdx dxdy; dxdy dydy]
%
% Example:
% [e1,e2,l1,l2] = eigen_decomposition([1,0; 0,1]);
% 
%   Author: Shawn Arseneau
%   Created: August 22, 2005
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%     if (nargout~=4 && nargout~=6)
%        error('input to eigen_decomposition must be a 2x2 matrix and there must be four output variables'); 
%     end

%     if ndims(M)==2
        [V, D] = eig(M);
        D = abs(D);
        l2 = D(1,1);
        l1 = D(2,2);
        e2 = V(:,1);
        e1 = V(:,2);

        if nargout==6
           l3 = D(3,3); 
           e3 = V(:,3);
           return;
        else %------- form [e1,e2,l1,l2]
            e3 = l1;
            l1 = l2;
        end
%     elseif ndims(M)==4
%         [rows, cols, st1, st2] = size(M);
%         tMat = zeros(st1, st2);
%         
%         if nargout==6 %---- 3D info
%             for r=1:rows
%                 for c=1:cols
%                     tMat(:,:) = M(r,c,:,:);
%                     [V D] = eigs(tMat);
%                     D = abs(D);
%                     l1(r,c) = D(1,1);
%                     l2(r,c) = D(2,2);
%                     l3(r,c) = D(3,3); 
%                     e1(r,c,:) = V(:,1);
%                     e2(r,c,:) = V(:,2);
%                     e3(r,c,:) = V(:,3);
%                 end
%             end
%         else  %---- 2D info
%             for r=1:rows
%                 for c=1:cols
%                     tMat(:,:) = M(r,c,:,:);
%                     [V D] = eigs(tMat);
%                     D = abs(D);
%                     e1(r,c,:) = V(:,1);
%                     e2(r,c,:) = V(:,2);
%                     e3(r,c) = D(1,1);
%                     l1(r,c) = D(2,2);
%                 end
%             end
%         end
%     else
%         error('Input to eigen decomposition in wrong array amount (NxN) or (NxMx2x2)');
%     end
end
end