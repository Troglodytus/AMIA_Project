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
figure;
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
figure;
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

% Display the enhanced phase image
figure;
subplot(1,2,1);
imshow(masked_phase(:,:,round(size(masked_phase,3)/2)), []);
title('Masked Unwrapped Phase Image');

subplot(1,2,2);
imshow(phase_hpf(:,:,round(size(phase_hpf,3)/2)), []);
title('High Pass Filtered Phase Image');


%Create a Phase Mask-------------
swi_mask = zeros(size(phase_hpf), 'like', phase_hpf);

% Applying the function element-wise
swi_mask(phase_hpf < 0) = 1;
swi_mask(phase_hpf > 1) = 0;
swi_mask(phase_hpf > 0 & phase_hpf < 1) = 1 - phase_hpf(phase_hpf > 0 & phase_hpf < 1);

%make nii
swi_masknii = make_nii(swi_mask);

%save the SWI mask
save_nii(swi_masknii, [data_path, 'swi_mask.nii']);

% Display the SWI mask
figure;
imshow(swi_mask(:,:,round(size(swi_mask,3)/2)), []);
title('SWI Mask');

%Generate SWI
n = 5;
swi = magnitude .* (swi_mask * n);

%make nii
swinii = make_nii(swi);

%save the SWI image
save_nii(swinii, [data_path, 'swi.nii']);

%Display the SWI image
figure;
imshow(swi(:,:,round(size(swi,3)/2)), []);
title('SWI');


%Initialize MIP image as SWI image
mipswi = swi; 

% Iterate over each slice
for i = 2:size(swi, 3)-1
    for j = 1:size(swi,1)
        for k = 1:size(swi,2)
            % Take the minimum of the current slice and its neighbors for each pixel
            prev_pixel = swi(j, k, i-1);
            current_pixel = swi(j, k, i);
            next_pixel = swi(j, k, i+1);
            mipswi(j, k, i) = min(current_pixel, min(prev_pixel, next_pixel));
        end
    end
end


%mipswi = min(swi, [], 3);

% Create a new NIfTI structure for the MIP SWI image
mipswinii = make_nii(mipswi);

% Save the MIP SWI image
save_nii(mipswinii, [data_path, 'mip_swi.nii']);

% Display the MIP SWI image
figure;
imshow(mipswi(:,:,round(size(mipswi,3)/2)), []);
title('Minimum Intensity Projection SWI');


%Compare initial Magnitude to SWI and mipSWI
figure;
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
figure;
imshow(magnitudeqsm(:,:,round(size(magnitudeqsm,3)/2)), []);
title('Magnitude Image');

%Brain mask with thresholding
threshold = 0.05 * max(magnitudeqsm(:));
brain_mask = magnitudeqsm > threshold;

%Apply morphological fix operations to mask
brain_mask = imfill(brain_mask, 'holes'); %fill holes
brain_mask = imopen(brain_mask, strel('disk', 5)); %Remove artefacts
brain_mask = single(brain_mask);

%show brain mask
figure;
imshow(brain_mask(:,:,round(size(brain_mask,3)/2)), []);
title('Brain Mask');

%sve the brain mask as nii
brain_masknii = make_nii(brain_mask);
save_nii(brain_masknii, [qsm_path, 'brain_mask.nii']);

% Apply the brain mask to the demagnetization field
masked_demag_field = demag_field .* brain_mask;

% Save the masked demagnetization field as a NIfTI file
masked_demag_fieldnii = make_nii(masked_demag_field);
save_nii(masked_demag_fieldnii, [qsm_path, 'masked_demag_field.nii']);

% Display the masked demagnetization field
figure;
imshow(masked_demag_field(:,:,round(size(masked_demag_field,3)/2)), []);
title('Masked Demagnetization Field');

%Create Dipole
% Dipole Field init
voxel_size = [0.64, 0.64, 0.64]; % Voxel size in mm
B0_dir = [0, 0, 1]; % B0 direction

% Generate the dipole kernel using the provided function
kernel_in_FD = 1; % Set to 1 for Fourier domain
dipole_FD = create_dipole_kernel(B0_dir, voxel_size, size(masked_demag_field), kernel_in_FD);

% Display the dipole kernel in Fourier domain
figure;
imshow(abs(dipole_FD(:,:,round(size(dipole_FD,3)/2))), []);
title('Dipole Kernel k Domain');

% Create QSM 
%FFT
masked_demag_field_fft = fftn(masked_demag_field);

%lambda regularization factor
lambda = 0.01;
%inverse filtering and Tikhonov regularization
qsm_freq_domain = masked_demag_field_fft .* conj(dipole_FD) ./ (abs(dipole_FD).^2 + lambda^2);

%ifft
qsm = real(ifftn(qsm_freq_domain));

%NaN treatment
qsm(isnan(qsm)) = 0;

% Display QSM
figure;
imshow(qsm(:,:,round(size(qsm,3)/2)), []);
title('Quantitative Susceptibility Mapping (QSM)');

% Save the QSM result as a NIfTI file
qsmnii = make_nii(qsm);
save_nii(qsmnii, [qsm_path, 'qsm_result.nii']);




% Threshold values to test
threshold_values = [0.01, 0.05, 0.1, 0.5, 1];

% Store results for comparison
qsm_results = cell(length(threshold_values), 1);

% Generate QSMs for each threshold
for i = 1:length(threshold_values)
    qsm_results{i} = single(generate_qsm_thresholded_inverse(masked_demag_field, dipole_FD, threshold_values(i)));
    
    % Display the QSM for each threshold
    figure;
    imshow(qsm_results{i}(:,:,round(size(qsm_results{i}, 3)/2)), []);
    %title(['QSM with Threshold = ', num2str(threshold_values[i]]);
end

% Save QSM results
for i = 1:length(threshold_values)
    qsmnii = make_nii(qsm_results{i});
    %save_nii(qsmnii, [qsm_path, 'qsm_result_thresh_', num2str(threshold_values[i]), '.nii']);
end


% Initialize arrays to store error metrics
mse_values = zeros(length(threshold_values), 1);
ssim_values = zeros(length(threshold_values), 1);

% Compare each generated QSM with the ground truth
for i = 1:length(qsm_results)
    % Calculate MSE
    mse_values(i) = immse(qsm_results{i}, chi_model);
    
    % Calculate SSIM
    ssim_values(i) = ssim(qsm_results{i}, chi_model);
    
    % Display results
    fprintf('Threshold = %f: MSE = %f, SSIM = %f\n', threshold_values(i), mse_values(i), ssim_values(i));
end

% Plot the comparison results
figure;
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
[~, best_idx] = max(ssim_values);
best_threshold = threshold_values(best_idx);
best_qsm = qsm_results{best_idx};

% Display the best QSM
figure;
imshow(best_qsm(:,:,round(size(best_qsm, 3)/2)), []);
title(['Best QSM with Threshold = ', num2str(best_threshold)]);

% Save the best QSM result as a NIfTI file
best_qsmnii = make_nii(best_qsm);
save_nii(best_qsmnii, [qsm_path, 'qsm_result_best_thresh_', num2str(best_threshold), '.nii']);


%----------------------------------SEGMENTATION------------------------------------------------------
addpath(genpath(path))
best_qsm = qsm_results{2};
im = best_qsm(:,:,160);

% Display the slice
figure;
imshow(im, []);
title('Select Initial Contour Points');

% Initialize combined segmentation mask
combined_segm = zeros(size(im));

% Number of ROIs to segment
num_rois = 2;

for i = 1:num_rois
    % Select initial contour points
    [y, x] = getpts;
    points = [x(:) y(:)];

    % Snake segmentation options
    Options = struct;
    Options.Verbose = true;
    Options.Iterations = 500;
    Options.nPoints = 50;
    Options.Wedge = 5;
    Options.Sigma1 = 5;
    Options.Sigma2 = 5;
    Options.Mu = 0.1;
    Options.Sigma3 = 3;
    Options.GIterations = 100;

    % Run the snake segmentation
    figure;
    [curve, segm] = Snake2D(im, points, Options);

    % Combine the segmentation results
    combined_segm = combined_segm | segm;

    % Display the current segmentation result
    figure;
    imshow(imoverlay(im, combined_segm, 'cyan'));
    title(['Segmented Regions (ROI ', num2str(i), ')']);
end

% Save the combined ROI
roi = combined_segm;
roi_path = fullfile(qsm_path, 'roi_slice_160.mat');
save(roi_path, 'roi');
disp(['ROI saved to ', roi_path]);

% Display the final combined segmentation result
figure;
imshow(imoverlay(im, combined_segm, 'cyan'));
title('Final Combined Segmented Regions');


















% Function to generate QSM with thresholded inverse filtering
function qsm_thresh = generate_qsm_thresholded_inverse(masked_demag_field, dipole_FD, threshold)
    % FFT
    masked_demag_field_fft = fftn(masked_demag_field);

    % Thresholded Inverse Filtering
    H = abs(dipole_FD) > threshold;
    qsm_freq_domain = zeros(size(masked_demag_field_fft));
    qsm_freq_domain(H) = masked_demag_field_fft(H) ./ dipole_FD(H);

    % Inverse Fourier Transform
    qsm_thresh = real(ifftn(qsm_freq_domain));

    % Replace NaNs with zeros in the resulting QSM
    qsm_thresh(isnan(qsm_thresh)) = 0;
end