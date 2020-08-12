function  [node,elem]  = roast_GenerateHeadModel(varargin)
% function roast_GenerateHeadModel(fullPathToT1,fullPathToT2)
% function roast_GenerateHeadModel(fullPathToT1,fullPathToT2,opts)

% This function uses roast implementation to segment the mri into 6 tissus
% {WM, GM, CSF, Skull, Scalp}
% and save the masks in the path of the T1 specified in the arguments.
% The roast toolbox is needed to run this function

% File created on December 2019, Takfarinas MEDANI
%                          July 2020, combine all the functions in one function
% Adapted from Roast toolbox to brainstoem call.

% opts can be used to tune the roast mesh option

fullPathToT1 = varargin{1};
fullPathToT2 =varargin{2};

if nargin<3
    opts.maxvol = 10;
    opts.reratio = 3; 
    opts.radbound = 5; 
    opts.angbound = 30; 
    opts.distbound = 0.3;
    % see cgalv2m for more information
end

%% Lookfor the roast toolbox and add it to the path
subjRSPD = fullPathToT1;
% you can either add a T2 image in order to improve the mesh
T2 = fullPathToT2;
[dirname,baseFilenameRSPD] = fileparts(subjRSPD);

% The rest of the function is as it is from the original function
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

%%  STEP 3: MESH GENERATION...

% suplementary options // maybe not required to brainstorm ===> remove it 
saveMeshFormatMat = 0;
saveMeshFormatMsh= 0;
meshOpt = struct('radbound',opts.radbound,'angbound',opts.angbound,...
                             'distbound',opts.distbound,'reratio',opts.reratio,...
                             'maxvol',opts.maxvol,'saveMeshFormatMat',saveMeshFormatMat,...
                             'saveMeshFormatMsh',saveMeshFormatMsh);
                         
uniqueTag = ['MeshModel_', num2str(opts.maxvol),'_',num2str(opts.reratio)...
                                            '_',num2str(opts.radbound), '_',num2str(opts.angbound) , '_',num2str(opts.distbound)];

if ~exist([dirname filesep baseFilenameRSPD '_' uniqueTag '.mat'],'file')
    disp('======================================================')
    disp('        STEP 3: MESH GENERATION...         ')
    disp('======================================================')
    [node,elem,face,allMask] = meshByIso2meshWithoutElectrode(subjRSPD,subjRSPD,T2,meshOpt,[],uniqueTag);
else
    disp('======================================================')
    disp('          MESH ALREADY GENERATED, SKIP STEP 3         ')
    disp('======================================================')
    load([dirname filesep baseFilename '_' uniqueTag '.mat'],'node','elem','face');
end

%% close all the figures 
close all

disp('======================================================')
disp('     STEP 3bis : SAVE THE MESH ...       ')
disp('======================================================')
if meshOpt.saveMeshFormatMat == 1;
    disp('saving mesh mat fomat...')
    save([dirname filesep baseFilename '_' uniTag '.mat'],'node','elem','face','allMask');
end
if meshOpt.saveMeshFormatMsh == 1;
    disp('saving mesh msh fomat...')
    maskName(1:numOfTissue) = {'WHITE','GRAY','CSF','BONE','SKIN','AIR'};
    savemsh(node(:,1:3),elem,[dirname filesep baseFilename '_' uniTag '.msh'],maskName);
end

end