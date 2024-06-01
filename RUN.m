clc; clear all; close all;
%--------------------------------------SWI------------------------------------------------------------------------------------
%------ SPECIFY DATA ------------------------------------------------------
path = 'C:\Users\Daniel\Documents\FHKärnten\Medical Engineering\2. Semester\Applied Medical Image Analysis\Project_Part1\';
helper_path = 'C:\Users\Daniel\Documents\FHKärnten\Medical Engineering\2. Semester\Applied Medical Image Analysis\Project_Part1\helper_functions\';
data_path = 'C:\Users\Daniel\Documents\FHKärnten\Medical Engineering\2. Semester\Applied Medical Image Analysis\Project_Part1\data\swi\';

%Load Helper Functions 
addpath(helper_path);
nifti_toolbox_path = 'C:\Users\Daniel\AppData\Roaming\MathWorks\MATLAB Add-Ons\Collections\Tools for NIfTI and ANALYZE image';
addpath(nifti_toolbox_path);

%Load NIfTI files
magnitudenii = load_nii([data_path, 'magnitude.nii']);
phasenii = load_nii([data_path, 'phase_unw.nii']);
masknii = load_nii([data_path, 'mask.nii']);

%extract img data of nifti structure
magnitude = magnitudenii.img;
phase = phasenii.img;
mask = masknii.img;

% Display the raw data
figure_handle = figure('Name', 'Raw MRI maps', 'NumberTitle', 'on');
%figure;
subplot(1,3,1);
imshow(magnitude(:,:,round(size(magnitude,3)/2)), []);
title('Magnitude Image');

subplot(1,3,2);
imshow(phase(:,:,round(size(phase,3)/2)), []);
title('Phase Image');

subplot(1,3,3);
imshow(mask(:,:,round(size(mask,3)/2)), []);
title('Mask');


%convert mask and magnitude from int16 to single values
mask = single(mask);
magnitude = single(magnitude);

% Apply the mask to the unwrapped phase image
masked_phase = phase .* mask;

% Create a new NIfTI structure for the masked phase image
masked_phasenii = make_nii(masked_phase);

% Save the masked phase image
save_nii(masked_phasenii, [data_path, 'masked_phase_unw.nii']);


% Display the magnitude and masked phase images
figure_handle = figure('Name', 'Masked Phase', 'NumberTitle', 'on');
subplot(1,2,1);
imshow(magnitude(:,:,round(size(magnitude,3)/2)), []);
title('Magnitude Image');

subplot(1,2,2);
imshow(masked_phase(:,:,round(size(masked_phase,3)/2)), []);
title('Masked Unwrapped Phase Image');

% Parameters for Gaussian low-pass filter ---> Unsharp Mask
sigma = 6;  %standard deviation of the Gaussian filter

% Apply the Gaussian low-pass filter and perform high-pass filtering slice-wise
high_pass_phase = zeros(size(masked_phase), 'like', masked_phase);
h = fspecial('gaussian', [6*sigma 6*sigma], sigma); % Create Gaussian filter

for i = 1:size(masked_phase, 3)
    low_pass_slice = imfilter(masked_phase(:,:,i), h, 'replicate');
    phase_hpf(:,:,i) = masked_phase(:,:,i) - low_pass_slice;
end

%make nii of phase high pass
phase_hpfnii = make_nii(phase_hpf);

% Save the enhanced phase image
save_nii(phase_hpfnii, [data_path, 'phase_highpassfiltered.nii']);

% Display the filtered phase image
figure_handle = figure('Name', 'High Pass Filtered Phase', 'NumberTitle', 'on');
subplot(1,2,1);
imshow(masked_phase(:,:,round(size(masked_phase,3)/2)), []);
title('Masked Unwrapped Phase Image');

subplot(1,2,2);
imshow(phase_hpf(:,:,round(size(phase_hpf,3)/2)), []);
title('High Pass Filtered Phase Image');


%Create a Phase Mask
swi_mask = zeros(size(phase_hpf), 'like', phase_hpf);

%apply the function elementwise
swi_mask(phase_hpf < 0) = 1;
swi_mask(phase_hpf > 1) = 0;
swi_mask(phase_hpf > 0 & phase_hpf < 1) = 1 - phase_hpf(phase_hpf > 0 & phase_hpf < 1);

%make nii
swi_masknii = make_nii(swi_mask);

%save the SWI mask
save_nii(swi_masknii, [data_path, 'swi_mask.nii']);

% Display the SWI mask
figure_handle = figure('Name', 'SWI Mask', 'NumberTitle', 'on');
imshow(swi_mask(:,:,round(size(swi_mask,3)/2)), []);
title('SWI Mask');

%Generate SWI apply mask 5 times
n = 5;
swi = magnitude .* (swi_mask * n);

%make nii
swinii = make_nii(swi);

%save the SWI image
save_nii(swinii, [data_path, 'swi.nii']);

%Display the SWI image
figure_handle = figure('Name', 'SWI', 'NumberTitle', 'on');
imshow(swi(:,:,round(size(swi,3)/2)), []);
title('SWI');


%Initialize MIP image as SWI image
mipswi = swi; 

% Iterate over each slice
for i = 3:size(swi, 3)-2
    for j = 1:size(swi,1)
        for k = 1:size(swi,2)
            %Search for minimum of current slice and its 4 neighbors for each pixel
            min_pixel = min([swi(j, k, i-2), swi(j, k, i-1), swi(j, k, i), swi(j, k, i+1), swi(j, k, i+2)]);
            mipswi(j, k, i) = min_pixel;
        end
    end
end



% Create a new NIfTI structure for the MIP SWI image
mipswinii = make_nii(mipswi);

% Save the MIP SWI image
save_nii(mipswinii, [data_path, 'mip_swi.nii']);

% Display the MIP SWI image
figure_handle = figure('Name', 'MIP SWI', 'NumberTitle', 'on');
imshow(mipswi(:,:,round(size(mipswi,3)/2)), []);
title('Minimum Intensity Projection SWI');


%Compare initial Magnitude to SWI and mipSWI
figure_handle = figure('Name', 'Mag/SWI/MIP Comparison', 'NumberTitle', 'on');
subplot(1,3,1);
imshow(magnitude(:,:,round(size(magnitude,3)/2)), []);
title('Magnitude');

subplot(1,3,2);
imshow(swi(:,:,round(size(swi,3)/2)), []);
title('SWI');
subplot(1,3,3);
imshow(mipswi(:,:,round(size(mipswi,3)/2)), []);
title('MinIP SWI');




%--------------------------------------QSM------------------------------------------------------------------------------------
qsm_path = 'C:\Users\Daniel\Documents\FHKärnten\Medical Engineering\2. Semester\Applied Medical Image Analysis\Project_Part1\data\qsm\';

%Load QSM nii files
magnitudeqsmnii = load_nii([qsm_path, 'magnitude.nii']);
demagnii = load_nii([qsm_path, 'demagnetization_field.nii']);
chi_modelnii = load_nii([qsm_path, 'chi_model.nii']);

%extract images
magnitudeqsm = magnitudeqsmnii.img;
demag_field = demagnii.img;
chi_model = chi_modelnii.img;

% Display the raw magnitude data
figure_handle = figure('Name', 'QSM Magnitude Raw', 'NumberTitle', 'on');
imshow(magnitudeqsm(:,:,round(size(magnitudeqsm,3)/2)), []);
title('Magnitude Image');
figure_handle = figure('Name', 'QSM Demag Field', 'NumberTitle', 'on');
imshow(demag_field(:,:,round(size(demag_field,3)/2)), []);
title('Demag Field');

% Initialize a 3D array to store the processed volume
brain_mask = zeros(size(magnitudeqsm));

% Process each slice
for i = 1:size(magnitudeqsm, 3)
    % Extract the current slice
    slice = magnitudeqsm(:, :, i);
    
    % Normalize and convert to grayscale if needed
    slice = mat2gray(slice);
    slice = imadjust(slice);

    % Detect edges using Canny edge detection
    edge_canny = edge(slice, 'canny', 0.05);

    % Close the disconnected edges using morphological closing
    edge_canny_closed = imclose(edge_canny, strel('disk', 1));

    % Automatically select the center point of the slice for filling
    [rows, cols] = size(edge_canny_closed);
    center_row = round(rows / 2);
    center_col = round(cols / 2);
    seed_point = false(rows, cols);
    seed_point(center_row, center_col) = true;

    % Fill the region starting from the center point
    edge_canny_filled = imfill(edge_canny_closed, 'holes');
    edge_canny_filled = imreconstruct(seed_point, edge_canny_filled);

    % Remove small, not-filled contours using morphological opening
    edge_canny_segm = imopen(edge_canny_filled, strel('disk', 5));

    % Store the processed slice back into the 3D array
    brain_mask(:, :, i) = edge_canny_segm;
end

brain_mask = single(brain_mask);
figure_handle = figure('Name', 'QSM Brain Mask', 'NumberTitle', 'on');
imshow(brain_mask(:,:,round(size(brain_mask,3)/2)), []);
title('Brain Mask');


% Save the brain mask as NIfTI
brain_masknii = make_nii(brain_mask);
save_nii(brain_masknii, [qsm_path, 'brain_mask.nii']);

% Apply the brain mask to the demagnetization field
masked_demag_field = demag_field .* brain_mask;

% Save the masked demagnetization field as a NIfTI file
masked_demag_fieldnii = make_nii(masked_demag_field);
save_nii(masked_demag_fieldnii, [qsm_path, 'masked_demag_field.nii']);

% Display the masked demagnetization field
figure_handle = figure('Name', 'Masked Demag Field', 'NumberTitle', 'on');
imshow(masked_demag_field(:,:,round(size(masked_demag_field,3)/2)), []);
title('Masked Demagnetization Field');

%------------------Create Dipole------------------------------
% Dipole Field init
voxel_size = [0.64, 0.64, 0.64]; % Voxel size in mm
B0_dir = [0, 0, 1]; % B0 direction

% Generate the dipole kernel using the provided function
kernel_in_FD = 1; % Set to 1 for Fourier domain
dipole_FD = create_dipole_kernel(B0_dir, voxel_size, size(magnitudeqsm), kernel_in_FD);

%permute kernel to correct rotation
% dipole_FD = permute(dipole_FD, [1, 3, 2]);

%In case the Kernel is of the wrong size, interpolate and resize to the
%correct dimensions
%target_size = size(demag_field);
%dipole_FD = imresize3(dipole_FD, target_size, 'linear');

%replace zeros with infinitesimally small values
%dipole_FD(dipole_FD == 0) = eps(10.0);

% Display the dipole kernel in Fourier domain
figure_handle = figure('Name', 'Dipole Kernel', 'NumberTitle', 'on');
imshow(abs(dipole_FD(:,:,round(size(dipole_FD,3)/2))), []);
title('Dipole Kernel k Domain');

%---------------------Create QSM-------------------------------
% FFT demag field
masked_demag_field_fft = fftshift(fftn(fftshift(masked_demag_field)));


figure_handle = figure('Name', 'Demag kSpace', 'NumberTitle', 'on');
imshow(masked_demag_field_fft(:,:,round(size(masked_demag_field_fft,3)/2)), []);
title('Demag kSpace');

%Dipole Division
qsm_freq_domain = masked_demag_field_fft ./ dipole_FD;
qsm_freq_domain(isinf(qsm_freq_domain)) = eps(10);
%qsm_freq_domain = single(qsm_freq_domain);


% IFFT
qsm = real(ifftshift(ifftn(ifftshift(qsm_freq_domain))));

% NaN treatment
%qsm(isnan(qsm)) = 0;

% Display QSM
figure_handle = figure('Name', 'QSM', 'NumberTitle', 'on');
imshow(qsm(:,:,round(size(qsm,3)/2)), []);
title('Quantitative Susceptibility Mapping (QSM)');

% Save the QSM result as a NIfTI file
qsmnii = make_nii(qsm);
save_nii(qsmnii, [qsm_path, 'qsm_result.nii']);








%--------------------------------------------------------------------------------
% %----------------Inverse filtering and Tikhonov regularization (Not working well)-------------------
% % % Shift the demagnetization field to center of k-space
% % masked_demag_field_shifted = fftshift(masked_demag_field);
% % 
% 
% % FFT
% masked_demag_field_fft = fftn(masked_demag_field);
% 
% %regularization factor
% lambda = 0.001;
% qsm_freq_domain = masked_demag_field_fft .* conj(dipole_FD) ./ (abs(dipole_FD).^2 + lambda^2);
% 
% 
% % Shift back the result to original position
% %qsm_freq_domain = ifftshift(qsm_freq_domain);
% 
% % IFFT
% qsm = real(ifftn(qsm_freq_domain));
% 
% % NaN treatment
% qsm(isnan(qsm)) = 0;
% 
% % Display QSM
% figure_handle = figure('Name', 'QSM Tikhonov', 'NumberTitle', 'on');
% imshow(qsm(:,:,round(size(qsm,3)/2)), []);
% title('Quantitative Susceptibility Mapping (QSM)');
% 
% % Save the QSM result as a NIfTI file
% qsmnii = make_nii(qsm);
% save_nii(qsmnii, [qsm_path, 'qsm_thikonov.nii']);
% 
% figure_handle = figure('Name', 'QSM Ground Truth', 'NumberTitle', 'on');
% imshow(chi_model(:,:,round(size(chi_model,3)/2)), []);
% title('QSM Ground Truth');



% Threshold values to test
threshold_values = [0.01, 0.05, 0.1, 0.2, 0.25];

% Store results for comparison
qsm_results = cell(length(threshold_values), 1);


% Generate QSMs for each threshold
for i = 1:length(threshold_values)
    % Generate QSM for the current threshold value
    qsm_results{i} = generate_qsm_thresholded_inverse(masked_demag_field_fft, dipole_FD, threshold_values(i));
end

% Display and save QSM results
for i = 1:length(threshold_values)
    % Display the QSM for each threshold
    figure_handle = figure('Name', ['QSM Thresh. = ', num2str(threshold_values(i))], 'NumberTitle', 'on');
    imshow(qsm_results{i}(:,:,round(size(qsm_results{i}, 3)/2)), []);
    title(['QSM with Threshold = ', num2str(threshold_values(i))]);
    
    % Save QSM result as a NIfTI file
    qsmnii = make_nii(single(qsm_results{i})); % Ensure saving in single format
    save_nii(qsmnii, [qsm_path, 'qsm_result_thresh_', num2str(threshold_values(i)), '.nii']);
end

% Initialize arrays to store error metrics
mse_values = zeros(length(threshold_values), 1);
ssim_values = zeros(length(threshold_values), 1);

% Apply the brain mask to the chi model
masked_chi_model = chi_model .* brain_mask;
figure_handle = figure('Name', 'QSM Ground Truth Enhanced', 'NumberTitle', 'on');
imshow(masked_chi_model(:,:,round(size(masked_chi_model,3)/2)), []);
title('QSM Ground Truth Enhanced');

% Compare each generated QSM with the ground truth (masked chi model)
for i = 1:length(qsm_results)
    
    % Calculate MSE
    mse_values(i) = immse(single(qsm_results{i}), single(masked_chi_model));
    
    % Calculate SSIM
    ssim_values(i) = ssim(single(qsm_results{i}), single(masked_chi_model));
    
    % Display results
    fprintf('Threshold = %f: MSE = %f, SSIM = %f\n', threshold_values(i), mse_values(i), ssim_values(i));
end

% Plot the comparison results
figure_handle = figure('Name', 'QSM Threshold Optimization', 'NumberTitle', 'on');
subplot(2,1,1);
plot(threshold_values, mse_values, '-o');
xlabel('Threshold');
ylabel('MSE');
title('MSE vs Threshold');

subplot(2,1,2);
plot(threshold_values, ssim_values, '-o');
xlabel('Threshold');
ylabel('SSIM');
title('SSIM vs Threshold');

% Select the best threshold based on the metrics (e.g., highest SSIM)
%[~, best_idx] = max(ssim_values);
[~, best_idx] = min(mse_values);
best_threshold = threshold_values(best_idx);
best_qsm = qsm_results{best_idx};

% Display the best QSM
figure_handle = figure('Name', 'Best QSM', 'NumberTitle', 'on');
imshow(best_qsm(:,:,round(size(best_qsm, 3)/2)), []);
title(['Best QSM with Threshold = ', num2str(best_threshold)]);

% Save the best QSM result as a NIfTI file
best_qsmnii = make_nii(best_qsm);
save_nii(best_qsmnii, [qsm_path, 'qsm_result_best_thresh_', num2str(best_threshold), '.nii']);


%----------------------------------SEGMENTATION------------------------------------------------------

addpath(genpath(path))
% best_qsm = qsm_results{3};
im = best_qsm(:,:,160);
roi_path1 = fullfile(qsm_path, 'roi_slice_160_seg1.mat');
roi_path2 = fullfile(qsm_path, 'roi_slice_160_seg2.mat');

%---------------------------- Interactive Snakes Segmentation------------------------

% % Display the slice
% figure_handle = figure('Name', 'Segmentation Slice', 'NumberTitle', 'on');
% imshow(im, []);
% title('Select Initial Contour Points');
% 
% % Initialize combined segmentation mask
% combined_segm = zeros(size(im));
% roi = combined_segm;
% 
% % Number of ROIs to segment
% num_rois = 2;
% 
% % Store brightness values
% brightness_values = zeros(1, num_rois);
% 
% for i = 1:num_rois
%     % Select initial contour points
%     [y, x] = getpts;
%     points = [x(:) y(:)];
% 
%     % Snake segmentation options
%     Options = struct;
%     Options.Verbose = true;
%     Options.Iterations = 20;
%     Options.nPoints = 50;
%     Options.Wedge = 5;
%     Options.Sigma1 = 5;
%     Options.Sigma2 = 6;
%     Options.Mu = 0.1;
%     Options.Sigma3 = 3;
%     Options.GIterations = 100;
% 
%     % Run the snake segmentation
%     figure_handle = figure('Name', 'Snake Segmentation', 'NumberTitle', 'on');
%     [curve, segm] = Snake2D(im, points, Options);
% 
%     % Save each segmented region separately
%     roi_path = fullfile(qsm_path, sprintf('roi_slice_160_seg%d.mat', i));
%     save(roi_path, 'segm');
%     disp(['ROI saved to ', roi_path]);
% 
%     % Combine the segmentation results
%     combined_segm = combined_segm | segm;
% 
%     % Calculate and store the average brightness value within the segmented region
%     brightness_values(i) = mean(im(segm));
% 
%     % Display the current segmentation result
%     figure_handle = figure('Name', 'Segmented Regions', 'NumberTitle', 'on');
%     imshow(im, []);
%     hold on;
% 
%     % Create a color mask for the overlay
%     red = cat(3, ones(size(im)), zeros(size(im)), zeros(size(im))); % Red color mask
%     h = imshow(red);
%     set(h, 'AlphaData', 0.5 * combined_segm); % Adjust transparency as needed
% 
%     title(['Segmented Regions (ROI ', num2str(i), ')']);
% end
% 
% % Display the final combined segmentation result
% figure_handle = figure('Name', 'Segmentation Result', 'NumberTitle', 'on');
% imshow(im, []);
% hold on;
% 
% % Create a color mask for the overlay
% red = cat(3, ones(size(im)), zeros(size(im)), zeros(size(im))); % Red color mask
% h = imshow(red);
% set(h, 'AlphaData', 0.5 * combined_segm); % Adjust transparency as needed
% 
% title('Final Combined Segmented Regions');
% 
% % Print the average brightness values for each ROI
% for i = 1:num_rois
%     fprintf('Average brightness value in ROI %d: %f\n', i, brightness_values(i));
% end
% 


%----------------------------------Load the Segmentation Result------------------------------------------------

% Define the paths to the saved ROI files
roi_path1 = fullfile(qsm_path, 'roi_slice_160_seg1.mat');
roi_path2 = fullfile(qsm_path, 'roi_slice_160_seg2.mat');

% Load each ROI from the saved MAT files
loaded_roi_data1 = load(roi_path1);
loaded_roi1 = loaded_roi_data1.segm;

loaded_roi_data2 = load(roi_path2);
loaded_roi2 = loaded_roi_data2.segm;

% Combine the loaded ROIs into one mask
combined_loaded_roi = loaded_roi1 | loaded_roi2;

% Display the loaded segmentation result
figure_handle = figure('Name', 'Loaded Segmentation Result', 'NumberTitle', 'on');
imshow(im, []);
hold on;

% Create a color mask for the overlay
red = cat(3, ones(size(im)), zeros(size(im)), zeros(size(im))); % Red color mask
h = imshow(red);
set(h, 'AlphaData', 0.5 * combined_loaded_roi); % Adjust transparency as needed

title('Loaded Segmented Regions');

% Calculate and display the average brightness value within each loaded segmentation region
mean_brightness1 = mean(im(loaded_roi1));
mean_brightness2 = mean(im(loaded_roi2));

fprintf('Average brightness value in region 1: %f\n', mean_brightness1);
fprintf('Average brightness value in region 2: %f\n', mean_brightness2);








%-----------------------------FUNCTIONS----------------------------------------------------
function qsm = generate_qsm_thresholded_inverse(masked_demag_field_fft, dipole_FD, threshold_value)

    % Truncate the dipole kernel in k-space
    dipole_thresholded_fft = dipole_FD;
    dipole_thresholded_fft(abs(dipole_FD) < threshold_value) = threshold_value;

    % Inverse filtering in the frequency domain
    qsm_freq_domain = masked_demag_field_fft ./ dipole_thresholded_fft;
    qsm_freq_domain(isinf(qsm_freq_domain)) = 0;

    % Inverse Fourier Transform to get the susceptibility map
    qsm = real(ifftshift(ifftn(ifftshift(qsm_freq_domain))));

    %NaN treatment
    %qsm(isnan(qsm)) = 0;
    
    %Permutation of dipole for plotting
    dipole_thresholded_fft = permute(dipole_thresholded_fft, [1,3,2]);
    
    figure_handle = figure('Name', ['Dipole = ', num2str(threshold_value)], 'NumberTitle', 'on');
    imshow(dipole_thresholded_fft(:,:,round(size(dipole_thresholded_fft, 3)/2)), []);
    title(['Dipole = ',num2str(threshold_value)]);
end


