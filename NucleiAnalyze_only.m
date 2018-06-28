% This script organizes relevant files into a certain logical order and
% passes them onto the blob analyze function to calculate average value of
% each channel within the structures defined by the mask image.

rehash
imageset1 = {[prefix exp_name '\w+.TIF'],['nuclei_' prefix exp_name '\w+.TIF']};
% imageset2 = {[prefix exp_name '\w+.TIF'],['nuclei_filt_' prefix exp_name '\w+.TIF']};
col_labels1 = nuclei_analyze(imageset1,sizemin,sizemax,strcat(exp_name,'_nuclei'),folder);
% col_labels2 = nuclei_analyze(imageset2,sizemin,sizemax,strcat(exp_name,'_nuclei_filtered'),folder);

