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
threshold = 0.1 * max(magnitudeqsm(:)); % Adjust threshold value if necessary
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




