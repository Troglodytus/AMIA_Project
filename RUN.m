clc; clear all; close all;

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

% Apply the mask to the unwrapped phase image
masked_phase = phase .* mask;

% Create a new NIfTI structure for the masked phase image
masked_phasenii = make_nii(masked_phase);












