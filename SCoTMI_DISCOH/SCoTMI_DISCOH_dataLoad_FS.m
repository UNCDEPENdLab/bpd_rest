%
% SCoTMI_DISCOH_dataLoad_FS
%
% Simple data loading script for functional MRI data, using standard
% freesurfer parcellation numbers/labels. 
% 
% Inputs: 
%       funcDataPath      : full path and file name for 4D fMRI data, 
%           .nii or  .nii.gz
%       parcDataPath      : full path and file name for 3D anatomical
%           parcellation template, where each voxel is a number associated
%           with the freesurfer lookup table.
%       FS_numbers_toUse  : a vector of numbers corresponding to freesurfer
%           labels, as found in $FREESURFER_HOME/FreeSurferColorLUT.txt
%           which we want to use (can be a small subset of the numbers in
%           the parcellation template volume file).
%       DO_LOADFUNCTIONAL : true OR false
%           default:true , skip loading functional data if we only want the
%           parcellation numbers etc.
%
% Outputs:
%       funcData          : 2D matrix, time x voxels
%       brainNumbers      : 1D vector of Freesurfer numbers actually being
%                           used, per voxel
%       u_brainNumbers    : The unique set from brainNumbers
%       brainNames        : Freesurfer names corresponding to
%                           u_brainNumbers
% Ben Cassidy, (UNSW, NeuRA) Dec 2014

function [funcData, brainNumbers, u_brainNumbers, brainNames] ...
    = SCoTMI_DISCOH_dataLoad_FS(funcDataPath, parcDataPath, FS_numbers_toUse,DO_LOADFUNCTIONAL)

assert(ischar(funcDataPath));
assert(ischar(parcDataPath));
assert(isnumeric(FS_numbers_toUse))
if isempty(DO_LOADFUNCTIONAL)
    DO_LOADFUNCTIONAL = true;
end

if ~exist('read_avw', 'file')
    path('/usr/local/fsl/etc/matlab', path);
end

brainNumbers = read_avw(parcDataPath);
brainNumbers = brainNumbers(:);

if DO_LOADFUNCTIONAL
    funcData = read_avw(funcDataPath);
    funcData = squeeze(reshape(funcData, [],1,1,size(funcData,4)));
else
    funcData = [];
end
%% Which voxels to keep?

voxToKeep = zeros(length(brainNumbers),1);
voxToKeep(ismember(brainNumbers, FS_numbers_toUse(:))) = 1;
if DO_LOADFUNCTIONAL, 
    funcData(voxToKeep==0,:) = [];
end
brainNumbers(voxToKeep==0) = [];
%%

%% Reshaping 
if DO_LOADFUNCTIONAL
funcData = funcData'; % now time x voxels
end
%%

if ~exist('read_fscolorlut', 'file')
    path(fullfile(getenv('FREESURFER_HOME'), 'matlab'), path)
end

u_brainNumbers = unique(brainNumbers(brainNumbers~=0));

[code, name, rgbv] = read_fscolorlut();
brainNames = name;
bufExistPos = ismember(u_brainNumbers, code);
assert(all(bufExistPos));
bufPos = ismember(code, u_brainNumbers);
brainNames(~bufPos,:) = [];



end