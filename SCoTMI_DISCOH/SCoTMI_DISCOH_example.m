
%% Loading data (can be done differently if needed)

% Defining the data locations
baseDir = fileparts(which('SCoTMI_DISCOH_example'));
funcDataPath = fullfile(baseDir, 'testData', 'filtered_func_data.nii.gz');
parcDataPath = fullfile(baseDir, 'testData', 'aparc+aseg.nii.gz');

% Which freesurfer regions do we use?
FS_numbers_toUse = [ 1, 2, 5, 11:20];  % just some subcortical regions
DO_LOADFUNCTIONAL = true;

% load full fMRI data and parcellation, 
% requires FSL and FreeSurfer matlab packages
[funcData, brainNumbers, u_brainNumbers, brainNames] ...
    = SCoTMI_DISCOH_dataLoad_FS(...
        funcDataPath, parcDataPath, FS_numbers_toUse,DO_LOADFUNCTIONAL);

% u_brainNumbers (and equivalent ROI names in brainNames) are possibly a
% subset of the ROI numbers we specified in FS_numbers_toUse; having
% removed those ROI numbers which did not exist in the parcellation volume.
    
% use the first principal component from the voxels in each node
Y = SCoTMI_DISCOH_dataPrep_fMRI(funcData, brainNumbers, u_brainNumbers);

%% Network estimation

% Y should be a matrix of size (# TimePoints x # nodes) without any empty
% or NaN elements.

% set some options (optional)
opt.vb = false;
opt.DO_PARALLEL = 0;

% run the calculation
[SCoTMI_raw,DISCOH_raw, NWvarTopo] = SCoTMI_DISCOH_calc(Y, opt);

% calculate an empirical null threshold using Bonferroni correction
pval = 0.001;
[nTimePoints,nNodes] = size(Y);
[nullThresh] = SCoTMI_H0calc(pval, nTimePoints, nNodes, [], [], [], []);

% calculate the thresholded SCoTMI links
SCoTMI_thresh = SCoTMI_raw;
SCoTMI_thresh(SCoTMI_thresh < nullThresh) = 0;

% plot some results
% figure();
% imagesc(SCoTMI_raw)
% title('Raw SCoTMI link strengths')
% figure(); imagesc(DISCOH_raw), title('Raw DISCOH link strengths')
% figure(); imagesc(SCoTMI_thresh);
% title('Thresholded SCoTMI link strengths')