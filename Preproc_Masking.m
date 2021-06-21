% Preprocessing
% Network masking of resting-state functional connectivity maps
% Networks: Default Mode Network, Basal Ganglia, and Motor hubs

% Written by Kristia Pamungkas
% Schweinfurt, Germany
% 21.06.2021

%% Initialization
tic
clear all
clc

%% Setting Up
addpath('...'); % add path
atlas = load_nii('/Users/kristiapamungkas/Desktop/Master Thesis/AAL3/ROI_MNI_V7.nii'); % load Atlas

%% Prepare DMN mask and control masks

% Atlas = AAL 3 Brain Atlas

DMN = [];
BG = []; % control: basal ganglia
Motor = []; % control: motor network

DMN_index = [19, 20, 39, 40, 69, 70, 71, 72, 151, 152, 153, 154, 155, 156];
BG_index = [75, 76, 77, 78, 79, 80, 161, 162, 163, 164]; 
Motor_index = [1, 2, 15, 16];

mask_DMN = atlas;
mask_BG = atlas;
mask_Motor = atlas;

DMN_Voxel = [];
BG_Voxel = [];
Motor_Voxel = [];

temp = [];

for i = 1:numel(DMN_index)
    temp = mask_DMN.img;
    temp(temp ~= (DMN_index(i))) = 0;
    DMN_Voxel(:,:,:,i) = temp;
end

for i = 1:numel(BG_index)
    temp = mask_BG.img;
    temp(temp ~= (BG_index(i))) = 0;
    BG_Voxel(:,:,:,i) = temp;
end

for i = 1:numel(Motor_index)
    temp = mask_Motor.img;
    temp(temp ~= (Motor_index(i))) = 0;
    Motor_Voxel(:,:,:,i) = temp;
end

DMN_Voxel = squeeze(sum(DMN_Voxel,4));
BG_Voxel = squeeze(sum(BG_Voxel,4));
Motor_Voxel = squeeze(sum(Motor_Voxel,4));

DMN_Voxel(DMN_Voxel > 0) = 1;
BG_Voxel(BG_Voxel > 0) = 1;
Motor_Voxel(Motor_Voxel > 0) = 1;


%% Specify Subjects Directory
mainfolder = '...'; % add path to RSFC maps
cd(mainfolder)

sub_dir = {}; % add subjects' ID


%% Looping

% Load RSFCMaps from every subject
% Mask them with DMN mask as well as M1 mask
% save them as a 3d array/matrix for each subject and each mask

for isub = 1:numel(sub_dir)
    subNii = load_nii(fullfile(mainfolder, ['RSFCMaps_C1_S' sub_dir{isub} '.nii']));
    
    subVoxel = subNii.img;
    
    % get overlapping voxels in atlas structure
    subVoxel_DMN = subVoxel .* DMN_Voxel;
    subVoxel_BG = subVoxel .* BG_Voxel;
    subVoxel_Motor = subVoxel .* Motor_Voxel;
    
    % store back to the nifti format
    subNii_DMN = subNii;
    subNii_DMN.img = subVoxel_DMN;
    
    subNii_BG = subNii;
    subNii_BG.img = subVoxel_BG;
    
    subNii_Motor = subNii;
    subNii_Motor.img = subVoxel_Motor;
    
    % save to nifti
    save_nii(subNii_DMN, fullfile(mainfolder, ['DMN_RSFCMaps_C1_S' sub_dir{isub} '.nii']));
    save_nii(subNii_BG, fullfile(mainfolder, ['BG_RSFCMaps_C1_S' sub_dir{isub} '.nii']));
    save_nii(subNii_Motor, fullfile(mainfolder, ['Motor_RSFCMaps_C1_S' sub_dir{isub} '.nii']));
    
end

toc