% Preprocessing
% Features Extraction of Resting-State Functional Connectivity Map
% Network: Default Mode Network

% Written by Kristia Pamungkas
% Schweinfurt, Germany
% 21.06.2021

%% Initialization
tic
clear all
clc

%% Setting Up
addpath('...'); % add path
atlas = load_nii('/Users/kristiapamungkas/Desktop/Master Thesis/AAL3/ROI_MNI_V7.nii'); % load atlas

%% Specify Subjects Directory
mainfolder = '...'; % add path of RSFC data
cd(mainfolder)

sub_dir = {}; % add ID of subjects

% DMN Hubs
FSM = [];
ACC = [];
PCC = [];
AngularGyrus = [];
Precuneus = [];

FSM_index = [19, 20];
ACC_index = [151, 152, 153, 154, 155, 156];
PCC_index = [39, 40];
AngularGyrus_index = [69, 70];
Precuneus_index = [71, 72];


%% Looping
% break down atlas image into layers containing 1 atlas structure each
nClusters = max(atlas.img(:));
for i = 1:nClusters
    atlasstructures{i} = atlas.img==i;
end

for isub = 1:numel(sub_dir)
    subNii = load_nii(fullfile(mainfolder, ['RSFCMaps_C1_S' sub_dir{isub} '.nii']));
    
    subVoxel = subNii.img;
    
    % get overlapping voxels in atlas structure
    masked = {};
    voxel_sum = {};
    
    for i = 1:nClusters
        thisStructure = atlasstructures{i};
        voxel_sum{i} = sum(thisStructure(:));
        masked{i} = sum(subVoxel(:).*thisStructure(:));
    end
    
    FSM(isub) = sum(vertcat(masked{FSM_index}))/sum(vertcat(voxel_sum{FSM_index}));
    ACC(isub) = sum(vertcat(masked{ACC_index}))/sum(vertcat(voxel_sum{ACC_index}));
    PCC(isub) = sum(vertcat(masked{PCC_index}))/sum(vertcat(voxel_sum{PCC_index}));
    AngularGyrus(isub) = sum(vertcat(masked{AngularGyrus_index}))/sum(vertcat(voxel_sum{AngularGyrus_index}));
    Precuneus(isub) = sum(vertcat(masked{Precuneus_index}))/sum(vertcat(voxel_sum{Precuneus_index}));
    
end

roi_tab = table(FSM', ACC', PCC', AngularGyrus', Precuneus');
writetable(roi_tab, 'DMN_AvgRSFC_table_C1.csv'); % save table
toc