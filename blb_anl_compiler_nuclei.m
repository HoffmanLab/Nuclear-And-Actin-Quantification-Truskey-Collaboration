%% Description
% This function is used to automatically concatinate and add onto a set of
% blb_anl_exp.txt files generated from a large experiment where it may be
% cumbersome to do so manually.

function blb_anl_compiler_nuclei(CompParam)

global geo_x...
    geo_y...
    venus_x...
    venus_y...
    venus...
    std_venus...
    FA_size...
    FA_eccentricity...
    FA_orientation...
    FAID... %10
    ImageID...
    construct...
    stain_name...
    FA_relative_orientation_AvgFA...
    FA_mean_orientation...
    sd_FA_mean_orientation...
    FA_1_nearest_neighbor_dist...
    FA_4_nearest_neighbors_dist...
    FA_8_nearest_neighbors_dist...
    unique_marker...
    n_FAs... %20

%% Input parameters for naming and exclusion
CompParam.FAminsize = 50;
CompParam.num_expIDs = length(CompParam.cell);
for n = 1:length(CompParam.cell)
    CompParam.exp_file{n} = file_search(['blb_anl_' CompParam.cell{n} '_' CompParam.dose{n} '_nuclei.txt'],CompParam.folder);
    rehash
    if n == 1
        b = load(fullfile(CompParam.folder,CompParam.exp_file{n}{1}));
        b(:,end+1:end+2) = 1;
        a = b;
    else
        b = load(fullfile(CompParam.folder,CompParam.exp_file{n}{1}));
        b(:,end+1:end+2) = n;
        a = vertcat(a,b);
    end
    clear('b');
end
sz = 13;
construct_col =  12;
stain_col =  13;

%% Make column modifications
% Invert y centroids (1948-col)
a(:,geo_y) = 1948-a(:,geo_y);
a(:,venus_y) = 1948-a(:,venus_y);

% Delete FAs with the following constraints
a(a(:,FA_size) < CompParam.FAminsize,:) = [];

%% Add columns for
% Unique marker
a(:,unique_marker) = 1000000.*a(:,construct)+100.*a(:,stain_name)+a(:,ImageID);

%% Calculations that require cell-average measurements
unique_imgs = unique(a(:,unique_marker));

for i = 1:length(unique_imgs)
    % All-or-nothing exclusion criteria
    rows = a(:,unique_marker)==unique_imgs(i);
    a(rows,n_FAs) = length(nonzeros(rows));
    if ~isempty(nonzeros(rows))
        % Calculate FA orientation relative to mean(FA orientation) cell-by-cell
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % THIS SECTION OF CODE MAY PRODUCE ERRORS FOR CELLS WITH FAS PRIMARILY
        % ORIENTED VERTICALLY
        vec = -pi/2:pi/10:pi/2;
        %         hist(a(rows,FA_orientation),length(vec));
        %         pause()
        %         close all
        up = linspace(pi/4,pi/2,6);
        down = linspace(pi/2,pi/4,6);
        vecref = horzcat(up(1:6),down(2:6));
        temp = a(rows,FA_orientation);
        [counts,~] = histc(temp,vec);
        peak_loc = counts==max(counts);
        peak = vec(peak_loc);
        peak = peak(1);
        for k = 1:length(vec)
            if peak == vec(k) && peak <= 0
                temp1 = abs(peak - temp);
                temprows = temp1 > vecref(k);
                temp(temprows) = temp(temprows) - pi;
            elseif peak == vec(k) && peak >= 0
                temp1 = abs(peak - temp);
                temprows = temp1 > vecref(k);
                temp(temprows) = temp(temprows) + pi;
            end
        end
        meanFAori = mean(temp);
        sdFAori = std(temp);
        a(rows,FA_mean_orientation) = meanFAori;
        a(rows,sd_FA_mean_orientation) = sdFAori;
        if meanFAori > (pi/2)
            a(rows,FA_mean_orientation) = a(rows,FA_mean_orientation) - pi;
        elseif meanFAori < (-pi/2)
            a(rows,FA_mean_orientation) = a(rows,FA_mean_orientation) + pi;
        end
        angls1 = a(rows,FA_orientation)+pi/2;
        angls2 = a(rows,FA_mean_orientation) + pi/2;
        a(rows,FA_relative_orientation_AvgFA) = rangle(angls1,angls2);
        clear('angls1','angls2')
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        % Calculate k nearest neighbors cell-by-cell
        temp = a(rows,[geo_x,geo_y]);
        [~,D1] = knnsearch(temp,temp,'K',2);
        [~,D5] = knnsearch(temp,temp,'K',5);
        [~,D9] = knnsearch(temp,temp,'K',9);
        D1(:,1) = [];
        D5(:,1) = [];
        D9(:,1) = [];
        D5 = mean(D5,2);
        D9 = mean(D9,2);
        a(rows,FA_1_nearest_neighbor_dist) = D1; %Distance to closest neighboring FA
        a(rows,FA_4_nearest_neighbors_dist) = D5; %Distance to closest 4 neighbors
        a(rows,FA_8_nearest_neighbors_dist) = D9; %Distance to closest 8 neighbors
        clear('rows','temp','angls1','angls2');
    end
end

%% Additional cell parameter columns
rows = find(a(:,FA_relative_orientation_AvgFA)>(pi/2));
a(rows,FA_relative_orientation_AvgFA) = abs(a(rows,FA_relative_orientation_AvgFA)-pi);
a(:,FA_relative_orientation_AvgFA) = a(:,FA_relative_orientation_AvgFA).*(180/pi);

%% Convert a to a cell so that it can contain both numerical and text data
% After final exclusion of NaN rows
nanrows = isnan(a);
nanrows = sum(nanrows,2);
nanrows = logical(nanrows);
a(nanrows,:) = [];

a = num2cell(a);
% Substitute in construct and stain names in the correct columns (52 and 53)
for j = 1:CompParam.num_expIDs
    [row,col] = find(cellfun(@(x,y) isequal(x,j),a));
    tmp = horzcat(row,col);
    tmp_c = tmp(tmp(:,2)==construct_col,:);
    tmp_s = tmp(tmp(:,2)==stain_col,:);
    for i = 1:length(tmp_c)
        a{tmp_c(i,1),tmp_c(i,2)} = CompParam.cell{j};
        a{tmp_s(i,1),tmp_s(i,2)} = CompParam.dose{j};
    end
    clear('tmp','tmp_c','tmp_s','row','col')
end

% Substitute text unique for the unique column [sz+length(intens)+4] and "low"
% "high" or "middle" for the final four columns.
for i = 1:length(a)
    a{i,unique_marker}=[a{i,construct} '_' a{i,stain_name} '_img' num2str(a{i,ImageID})]; % exp_stain_image_cell
end

headers = CompParam.headers;

a = vertcat(headers,a);
xlswrite(fullfile(CompParam.folder,'nuclear_shape_compiled.xlsx'),a);
end

function angles = rangle(angles1,angles2)
angles1 = angles1 + pi/2;
angles2 = angles2 + pi/2;
angles(:,1) = abs(angles1-angles2);
angles(:,2) = angles2 - (angles1-pi);
angles(:,3) = abs(angles2 - (angles1+pi));
angles = min(angles,[],2);
end

