close all; clear all; clc;
%% At focus
filepath = []; %strcat('../scan_550nm/');
dark_filename = strcat(filepath,'spec_dark_10ms.txt');
dark0 = importdata(dark_filename);
 
%% Read spectrum files and store data
[M, N] = size(dark0);
wavelength_samples = M;
samples = 101;
spec_matrix0 = zeros(samples,samples, wavelength_samples);
for i = 1:1:samples
    for j = 1:1:samples
        scan_name = strcat('spec_','x',num2str(j),'_y',num2str(i),'.txt');
        spec_filename = strcat(filepath,scan_name);
        spec0 = importdata(spec_filename);
        psf1 = spec0(:,2) - dark0(:,2);
        index = find(psf1 < 0); psf1(index) = 0;
        spec_matrix0(j,i,:) = psf1';
   end
       
end

save('Coll_450nm.mat','spec_matrix0');

wavelengths = spec0(:,1); 
% Note that for wavelength index below 661, all values are zero.
% Corresponds to 486.2757nm.

% max_values = max(max(spec_matrix0,[],1),[],2);
% spec_matrix0 = spec_matrix0./max_values;

% zslice = wavelengths(661:50:end);
% figure; h = slice(x_vec, y_vec, wavelengths, spec_matrix0, 0, [], zslice);
% set(h,'EdgeColor','none','FaceColor','interp','FaceAlpha',0.25);
% figure; imagesc(x_vec, wavelengths, squeeze(spec_matrix0(:,round(samples/2),:))');
% xlabel('Wavelengths(nm)'); ylabel('X-slice'); axis([min(x_vec) max(x_vec) 450 850 ]); colormap('hot');
% figure; imagesc(y_vec, wavelengths, squeeze(spec_matrix0(round(samples/2),:,:))');
% xlabel('Wavelengths(nm)'); ylabel('Y-slice'); axis([min(x_vec) max(x_vec) 450 850]); colormap('hot');


temp = 466:635;
hsi_psf = zeros(61,61,length(temp));
for cnt=1:length(temp)
    temp = spec_matrix0(:,:,cnt+673);
    if max(max(temp)) > 0 & ~isnan(max(max(temp)))
    [i1 i2] = find(temp == max((max(temp))));
    hsi_psf(:,:,cnt) = temp(i1-30:i1+30,i2-30:i2+30);
    hsi_psf(:,:,cnt) = hsi_psf(:,:,cnt)/max(max(hsi_psf(:,:,cnt)));
    end
end

zslice = wavelengths(674:100:1600);
x_vec = (0:60)*10; x_vec = x_vec - mean(x_vec); % micrometers
y_vec = x_vec;
figure(3); h = slice(x_vec, y_vec, wavelengths(674:1600), hsi_psf, 0, [], zslice);
set(h,'EdgeColor','none','FaceColor','interp','FaceAlpha',0.25);
colormap('winter');
figure(4); h = slice(x_vec, y_vec, wavelengths(674:1600), hsi_psf, 0, 0, []);
colormap('winter');
set(h,'EdgeColor','none','FaceColor','interp','FaceAlpha',0.25);

% plot individual PSFs
design_wavelengths = 500:25:825; %492:515;% 590:613; %540:1:563; 
lens_radius = 50e-3;
focal_length = 200e-3;
dx = 1; % microns
x1 = min(x_vec):dx:max(x_vec); % 1um for interpolation.
y1 = x1;
for cnt=1:length(design_wavelengths)
    [minValue,closestIndex] = min(abs(wavelengths-design_wavelengths(cnt)));
    %figure(1); subplot(4,6,cnt); plot(x_vec, hsi_psf(:,30,closestIndex-673)/max(hsi_psf(:,30,closestIndex-673)), 'LineWidth',2); grid on; title(sprintf('%0.2f nm',wavelengths(closestIndex)));
    psf_interp = interp2(x_vec, y_vec', hsi_psf(:,:,closestIndex-673),x1,y1', 'linear');
    figure(1); subplot(4,6,cnt); plot(x1, psf_interp(:,301)/max(max(psf_interp)), 'LineWidth',2); grid on; title(sprintf('%0.2f nm',wavelengths(closestIndex))); set(gca,'FontSize',16);
    %[mtf_1D,fx_1D] = psf2mtf(hsi_psf(:,:,closestIndex-673),61);
    [mtf_1D,fx_1D] = psf2mtf(psf_interp,length(x1)*1e-3);
    figure(2); hold on; plot(fx_1D-min(fx_1D), mtf_1D, 'LineWidth',2); grid on; % title(sprintf('%0.2f nm',wavelengths(closestIndex))); axis([0 100 0 1]);
    % plot the diffraction-limited MTF
%     fx0 = lens_radius/(design_wavelengths(cnt)*focal_length)*1e6; % convert to 1/um.
%     mtf_DL = 2/pi*(acos(fx_1D/(2*fx0))-fx_1D/(2*fx0).*sqrt(1-(fx_1D/(2*fx0)).^2)); 
%     index = find(fx_1D > 2*fx0); mtf_DL(index) = 0;
%     hold on; subplot(4,6,cnt); plot(fx_1D, mtf_DL, '--','LineWidth',2);
    %figure(5); subplot(4,6,cnt); imagesc(x_vec, y_vec, hsi_psf(:,:,closestIndex-673)); axis equal; colormap('hot'); title(sprintf('%0.2f nm',wavelengths(closestIndex)));
    %hold on; 
    %axes('Position',[.7 .7 .2 .2])
    %box on;
    %plot(x_vec, hsi_psf(:,30,closestIndex-673)/max(hsi_psf(:,30,closestIndex-673)), 'y', 'LineWidth',2); grid on;
end
legend('500nm', '525nm', '550nm', '575nm', '600nm', '625nm', '650nm', '675nm', '700nm', '725nm', '750nm', '775nm', '800nm', '825nm');
set(gca,'FontSize',16);

% Diffraction-limited MTF
for cnt=1:length(design_wavelengths)
    % plot the diffraction-limited MTF
    fx0 = lens_radius/(design_wavelengths(cnt)*focal_length)*1e6; % convert to 1/um.
    mtf_DL = 2/pi*(acos(fx_1D/(2*fx0))-fx_1D/(2*fx0).*sqrt(1-(fx_1D/(2*fx0)).^2)); 
    index = find(fx_1D > 2*fx0); mtf_DL(index) = 0;
    hold on; plot(fx_1D, mtf_DL,'LineWidth',2);
end
legend('500nm', '525nm', '550nm', '575nm', '600nm', '625nm', '650nm', '675nm', '700nm', '725nm', '750nm', '775nm', '800nm', '825nm');
set(gca,'FontSize',16);
axis([0 50 0 1]); grid on;

% plot individual PSFs no sub-plots
% design_wavelengths = [499.82 542 601.94 656.08 757.95];
% %500:25:825;
% %[499.82 542 601.94 656.08 700.95 757.95 800 825 850] %492:515;% 590:613; %540:1:563; 
% for cnt=1:length(design_wavelengths)
%     [minValue,closestIndex] = min(abs(wavelengths-design_wavelengths(cnt)));
%     figure; imagesc(x_vec, y_vec, hsi_psf(:,:,closestIndex-673)); axis equal; colormap('hot'); title(sprintf('%0.2f nm',wavelengths(closestIndex)));
%     hold on; 
%     axes('Position',[0.272 0.6 .5 .3])
%     % box on;
%     plot(x_vec, hsi_psf(:,30,closestIndex-673)/max(hsi_psf(:,30,closestIndex-673)), 'y', 'LineWidth',1); grid on; axis off;
% end

 % plot efficiency. 
[M N P] = size(hsi_psf);
eff = zeros(1,P);
dark_filename = strcat('../Efficiency_measurements/before_mdl_bb_dark_2000ms.txt');
dark0 = importdata(dark_filename);
spec_filename = strcat('../Efficiency_measurements/before_mdl_bb_spec_2000ms.txt');
spec0 = importdata(spec_filename);
incident_power = spec0-dark0;
incident_power = incident_power/2000*10; % scale to 10ms integration time.
incident_power = incident_power*(100e-3/10e-6)^2; % scale to area of MDL.
for cnt=1:P
    eff(cnt) = sum(sum(hsi_psf(:,:,cnt)))./incident_power(cnt+673,2)*100;
end
figure; plot(wavelengths(674:1600),eff); xlabel('Wavelengths(nm'); ylabel('Efficiency(%)');
set(gca,'FontSize',16);

% transmission efficiency
dark_filename = '../Efficiency_measurements/after_mdl_bb_dark.txt';
dark1 = importdata(dark_filename);
spec_filename = '../Efficiency_measurements/after_mdl_bb_spec.txt';
spec1 = importdata(spec_filename);
transmitted_power = spec1-dark1;

dark_filename = strcat('../Efficiency_measurements/before_mdl_bb_dark_2000ms.txt');
dark0 = importdata(dark_filename);
spec_filename = strcat('../Efficiency_measurements/before_mdl_bb_spec_2000ms.txt');
spec0 = importdata(spec_filename);
incident_power = spec0-dark0;
incident_power = incident_power/2000*10; % scale to 10ms integration time.

transmitted_eff = transmitted_power./incident_power;
figure; plot(wavelengths,transmitted_eff*100); xlabel('Wavelengths(nm)'); ylabel('Efficiency(%)');
set(gca,'FontSize',16);

