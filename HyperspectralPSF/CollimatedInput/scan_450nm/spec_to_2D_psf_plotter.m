close all; clear all; clc;
%% Read Lambda values
lamda_values_matrix = importdata('lambda_values.txt');
%% Read dark file
% lambda = 550;
%filepath = strcat('D:\Apratim\2022_12_09_Large_area_MDL_spectrometer_psf\Collimated_beam\8um_fiber\scan_450nm\');
filepath = [];
dark_filename = strcat(filepath,'spec_dark_10ms.txt');
dark0 = importdata(dark_filename);
% dark = dark0(:,2);
% avg_dark = sum(dark)./length(dark);
%% Read spectrum files and store data
samples = 101;
lambda_value = 450.231598;
% lambda_value = 750.219260;
% lambda_value = 645.072958;
% lambda_value = 655.147077;
% lambda_value = 650.113574;
% lambda_value = 550.133358;
% lambda_value = 555.302346
% lambda_value = 544.957666;
% lambda_value = 550.564363;
% lambda_value = 551.426234;
% lambda_value = 650.533304;
% lambda_value = 555.302346;
% lambda_value = 560.464611;
% lambda_value = 600.667353;
% lambda_value = 650.533304;
% lambda_value = 750.219260;
% lambda_value = 535.018933;
spec_matrix0 = zeros(samples,samples);
for i = 1:1:samples
    for j = 1:1:samples
        scan_name = strcat('spec_','x',num2str(j),'_y',num2str(i),'.txt');
        spec_filename = strcat(filepath,scan_name);
        spec0 = importdata(spec_filename);
%         spec1 = spec0(:,2);
        [M,I] = find(spec0==lambda_value);
        psf0 = spec0(M,2);
        dark = dark0(M,2);
%         psf1 = psf0;
        psf1 = psf0 - dark; 
        if psf1<0
                psf1 = 0;
        end
        spec_matrix0(j,i) = psf1;
        %fprintf(strcat('\n',scan_name,'\n'));     
    end  
end
%%
psf_norm = spec_matrix0./max(max(spec_matrix0));
figure(1)
imagesc(psf_norm); axis square; colorbar; caxis([0 1]);
title(strcat('lambda =',num2str(lambda_value),'nm'));
filename_save = strcat('lambda_',num2str(lambda_value),'nm','.png');
saveas(gcf, filename_save);