%% Script to call blb_anl_compiler_nuclei and feed it the correct inputs
clc
clear
close all

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
    unique_marker... %20
    n_FAs...
    
geo_x = 1;
geo_y = 2;
venus_x = 3;
venus_y = 4;
venus = 5;
std_venus = 6;
FA_size = 7;
FA_eccentricity = 8;
FA_orientation = 9;
FAID = 10; %10
ImageID = 11;
construct = 12;
stain_name = 13;
FA_relative_orientation_AvgFA = 14;
FA_mean_orientation = 15;
sd_FA_mean_orientation = 16;
FA_1_nearest_neighbor_dist = 17;
FA_4_nearest_neighbors_dist = 18;
FA_8_nearest_neighbors_dist = 19;
unique_marker = 20; %20a
n_FAs = 21;

CompParam.headers = {...
    'Geometric Centroid (X)',...
    'Geometric Centroid (Y)',...
    'DAPI Weighted Centroid (X)',...
    'DAPI Weighted Centroid (Y)',...
    'DAPI Intensity (a.u.)',...
    'std(DAPI Intensity) (a.u.)',...
    'Nuclear Area (microns^2)',...
    'Nuclear Eccentricity',...
    'Nuclear Orientation',...
    'Nuclei ID',... %10
    'Image ID',...
    'Cell Type',...
    'Dose',...
    'Nuclei Orientation (Rel. Avg. Nuclei)',...
    'Nuclei Orientation (Avg.)',...
    'Nuclei Orientation (Std.)',...
    'Nuclei 1 Nearest Neighbor Distance',...
    'Nuclei 4 Nearest Neighbors Distance',...
    'Nuclei 8 Nearest Neighbors Distance',...
    'Unique Image',...
    'n Nuclei'};

CompParam.folder = input('Folder containing text files: ','s');
CompParam.cell = {...
    'C2C12_1',...
    'C2C12_1',...
    'C2C12_1',...
    'C2C12_1',...
    'C2C12_2',...
    'C2C12_2',...
    'C2C12_2',...
    'C2C12_2',...
    'C2C12_3',...
    'C2C12_3',...
    'S101',...
    'S101',...
    'S101',...
    'S101',...
    'S103',...
    'S103',...
    'S103',...
    'S103',...
    'T100',...
    'T100',...
    'T100',...
    'T100',...
    'T101',...
    'T101',...
    'T101',...
    'T101'};
CompParam.dose = {...
    '0',...
    '100',...
    '1000',...
    '10000',...
    '0',...
    '100',...
    '1000',...
    '10000',...
    '100',...
    '1000',...
    '0',...
    '100',...
    '1000',...
    '10000',...
    '0',...
    '100',...
    '1000',...
    '10000',...
    '0',...
    '100',...
    '1000',...
    '10000',...
    '0',...
    '100',...
    '1000',...
    '10000'};

blb_anl_compiler_nuclei(CompParam);