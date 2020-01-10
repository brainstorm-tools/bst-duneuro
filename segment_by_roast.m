function segment_by_roast(fullPathToT1,fullPathToT2)
% function segment_by_roast(fullPathToT1,fullPathToT2)
% This function uses roast implementation to segment the mri into 6 tissus
% {WM, GM, CSF, Skull, Scalp}
% and save the masks in the path of the T1 specified in the arguments.
% The roast toolbox is needed to run this function 

% File created on December, Takfarinas MEDANI
% Adapted from Roast toolbox.

%% Lookfor the roast toolbox and add it to the path
%  find the bst_duneuro_toolbox path
str = which('roast','-all'); % 
if ~isempty(str)
[filepath,~,~] = fileparts(str{1});    
% add to the path
addpath(genpath(filepath));
end
% read the input 
subjRSPD = fullPathToT1;
% you can either add a T2 image in order to improve the mesh
T2 = fullPathToT2;

[dirname,baseFilenameRSPD] = fileparts(subjRSPD);

%%  STEP 1 : SEGMENT THE MRI...
% step 1  : segmentation using spm ... goo deeper in order to have more
% details about this function
if (isempty(T2) && ~exist([dirname filesep 'c1' baseFilenameRSPD '_T1orT2.nii'],'file')) ||...
        (~isempty(T2) && ~exist([dirname filesep 'c1' baseFilenameRSPD '_T1andT2.nii'],'file'))
    disp('======================================================')
    disp('       STEP 1 : SEGMENT THE MRI...          ')
    disp('======================================================')
    start_seg(subjRSPD,T2);
else
    disp('======================================================')
    disp('          MRI ALREADY SEGMENTED, SKIP STEP 1          ')
    disp('======================================================')
end

%%  STEP 2 : SEGMENTATION TOUCHUP...
if (isempty(T2) && ~exist([dirname filesep baseFilenameRSPD '_T1orT2_masks.nii'],'file')) ||...
        (~isempty(T2) && ~exist([dirname filesep baseFilenameRSPD '_T1andT2_masks.nii'],'file'))
    disp('======================================================')
    disp('     STEP 2 : SEGMENTATION TOUCHUP...       ')
    disp('======================================================')
    segTouchup(subjRSPD,T2);
else
    disp('======================================================')
    disp('    SEGMENTATION TOUCHUP ALREADY DONE, SKIP STEP 2    ')
    disp('======================================================')
end

end