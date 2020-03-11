function  errMsg = generate_conductivity_tensor(iSubject, iAnatomy)


%% Check if brainsuite is available otherwise ask user to install it and add link to that
[status,cmdout] = system('bdp --version');
if status ~= 0
    errMsg = [('Brainsuite is not installed in you computer, please visite : http://forums.brainsuite.org/download/ and xxx brainstorm tutorial')];
    return
end

[status,cmdout] = system('headreco --h');
if status ~= 0
    errMsg = ['The FEM tensor generation works only with SimNibs mesh and mask :  xxx brainstorm tutorial'];
    return
end

% maybe you need to save also the mask from SimNibs and laod it to
% brainstorm .... 
% generate the mesh 
bst_progress('start', 'Generate Conductivity tensor', 'Generate FEM Conductivity tensor...');

iSubject

iAnatomy
sSubject
% Get Protocol information
ProtocolInfo     = bst_get('ProtocolInfo');
ProtocolSubjects = bst_get('ProtocolSubjects');
% Default subject
if (iSubject == 0)
	sSubject = ProtocolSubjects.DefaultSubject;
% Normal subject 
else
    sSubject = ProtocolSubjects.Subject(iSubject);
end

%%  Ask the user for the DTI data path : The asociated bvec and bval should be in the same folder 
%% ===== SELECT MRI FILE =====
% If MRI file to load was not defined : open a dialog box to select it
% Get last used directories
LastUsedDirs = bst_get('LastUsedDirs');
% Get last used format
DefaultFormats = bst_get('DefaultFormats');
% Get MRI file
[DtiFile, FileFormat, FileFilter] = java_getfile( 'open', ... % DialogType
    'Import DWI data...', ...              % Window title
    LastUsedDirs.ImportAnat, ...      % Default directory
    'single', 'files_and_dirs', ... % Selection mode   ===>  TODO : could be multiple by including the bval and bvec 
    bst_get('FileFilters', 'mri'), ... % FilesOrDir
    DefaultFormats.MriIn... % Filters
    );
% If no file was selected: exit
if isempty(DtiFile)
    return
end
% Save default import directory
LastUsedDirs.ImportAnat = bst_fileparts(DtiFile);
bst_set('LastUsedDirs', LastUsedDirs);
% Save default import format
DefaultFormats.MriIn = FileFormat;
bst_set('DefaultFormats',  DefaultFormats);

%% Check the DTI data 
bst_progress('text', 'Checking the DWI data...');
% check the format of the input
[filepath,name,ext]  =  fileparts(DtiFile);
if ~strcmp(ext,'.nii')
        error('The selected file is not a NifTi file')   
end
% Load the DTI and check the size
[DTI_MRI, DTI_vox2ras] = in_mri_nii(DtiFile, 1, 1);
if length(size(DTI_MRI.Cube)) ~= 4
    errMsg = ['The selected file is not a 4D format'];
    return
end
if size(DTI_MRI.Cube) <= 6
    errMsg = ['The DWI should have at least 6 values on the 4th dimention (6 gradient directions)'];
    return
end
%% Assuming that the name of the DTI is similar to the name of the bval and bvec 
% TODO adapte the code for other namings
% for now it assumes that the data will have the name toto.nii  toto.bval
% and toto.bvec
% Check if the bval and bvec are in the same folder
% MriFile = MriFile{:};

bvalFileData = dir(fullfile(filepath,'*.bval'));
if isempty(bvalFileData)
    errMsg = ['There is no "*.bval" file in the selected folder'];
    return
end

bvecFileData = dir(fullfile(filepath,'*.bvec'));
if isempty(bvecFileData)
    errMsg = ['There is no "*.bvec" file in the selected folder'];
    return
end

% check the name of the data ==> shoul match to be arranged by subject
if ~contains(bvalFileData.name,name(1:4))
    errMsg = [('We can not match the  "*.bval" file to the MRI, please associate the same name for all the data belonging to the same subject')];
    return
end

if ~contains(bvecFileData.name,name(1:4))
    errMsg = [('We can not match the  "*.bvec" file to the MRI, please associate the same name for all the data belonging to the same subject')];
    return
end

% At this level we assume that all the data are present on the folder
%% Start here 
% === SAVE MRI AS NII ===
bst_progress('text', 'Exporting MRI...');
% Empty temporary folder, otherwise it may reuse previous files in the folder
gui_brainstorm('EmptyTempFolder');
% Create temporary folder for segmentation files
brainsuiteDir = bst_fullfile(bst_get('BrainstormTmpDir'), 'brainsuite');
mkdir(brainsuiteDir);
% Save MRI in .nii format
T1File = file_fullpath(sSubject.Anatomy(1).FileName);

subjid = strrep(sSubject.Name, '@', '');
T1Nii = bst_fullfile(brainsuiteDir, [subjid 'T1.nii']);
sMriT1 = in_mri_bst(T1File);
out_mri_nii(sMriT1, T1Nii);

%% Call the  brainsuite pipline
bst_progress('text', 'Calling Brainsuite ...');
OutPutFolder = brainsuiteDir;
T1Filename = T1Nii;
dtiFilename = DtiFile;
bvecFilename = fullfile(bvecFileData.folder,bvecFileData.name);
bvalFilename = fullfile(bvalFileData.folder,bvalFileData.name);
% 1-  Brain Surface Extractor or BSE tool 
bst_progress('text', 'Calling Brainsuite "bse" ...');
% This version take account of the weired files names and paths
% T1Filename = '1001_T1.nii'
[status_bse,cmdout_bse]  = system(['bse -i ' T1Filename  ' --auto '  ...
                                                                                            ' -o ' fullfile(OutPutFolder,'skull_stripped_mri.nii.gz')...
                                                                                            ' --mask ' fullfile(OutPutFolder,'bse_smooth_brain.mask.nii.gz')...
                                                                                            ' --hires ' fullfile(OutPutFolder, 'bse_detailled_brain.mask.nii.gz')]);%...
%                                                                                             ' --cortex ' fullfile(OutPutFolder, 'bse_cortex_file.nii.gz')]);
if status_bse ~= 0
    brainsuite_bdp_logfile = fullfile(OutPutFolder,'brainsuite_bdp_logfile.txt');
    fid = fopen(brainsuite_bdp_logfile , 'w');
    fprintf(fid, '%st', cmdout_bse);
    fprintf(fid, '%st', ' Load the file brainsuite_bdp_logfile for manual debuging/data checking');
    fclose(fid);
    errMsg = ['Brainsuite Process failled at the "bse" step 01, please check the output logfile: %s ', brainsuite_bdp_logfile];
    return    
end
% [status_bse,cmdout_bse]  = system(['bse -i ' pathToT1 ' -o skull_stripped_mri.nii.gz --mask bse_mask.nii.gz  --hires bse_brainDetailled.nii.gz --cortex bse_cortex.nii.gz ' ]);

% 2- Non-uniformity correction:  using BrainSuite’s bias field correction
bst_progress('text', 'Calling Brainsuite "bfc" ...');
% software (BFC tool) :  Important: BFC should only be applied to skull-stripped brain images.
% Use bse data
[status_bfc,cmdout_bfc]  = system(['bfc -i ' fullfile(OutPutFolder,'skull_stripped_mri.nii.gz')...
                                                                                                                ' -o ' fullfile(OutPutFolder,'output_mri.bfc.nii.gz')...
                                                                                                                ' -L 0.5 -U 1.5']);
if status_bfc ~= 0
    brainsuite_bdp_logfile = fullfile(OutPutFolder,'brainsuite_bdp_logfile.txt');
    fid = fopen(brainsuite_bdp_logfile , 'w');
    fprintf(fid, '%st', cmdout_bse);
    fprintf(fid, '%st', ' Load the file brainsuite_bdp_logfile for manual debuging/data checking');
    fclose(fid);
    errMsg = ['Brainsuite Process failled at the "bse" step 01, please check the output logfile: %s '  brainsuite_bdp_logfile];
    return
end

% 3-BDP :
bst_progress('text', 'Calling Brainsuite "bdp" ...');
% Check bval and bvec
bval = load(bvalFilename);
bvec = load(bvecFilename);
if length(bval) ~= length(bvec)
    error('bval and bvec have diffrent length, check these files before calling bdp');
end
[status_bdp,cmdout_bdp]  = system(['bdp ' fullfile(OutPutFolder,'output_mri.bfc.nii.gz') ' --tensor --nii ' dtiFilename...
                                                                                           '  --t1-mask '  fullfile(OutPutFolder,'bse_smooth_brain.mask.nii.gz')...
                                                                                           '  -g ' bvecFilename ' -b ' bvalFilename]);

if status_bdp ~= 0
    brainsuite_bdp_logfile = fullfile(OutPutFolder,'brainsuite_bdp_logfile.txt');
    fid = fopen(brainsuite_bdp_logfile , 'w');
    fprintf(fid, '%st', cmdout_bse);
    fprintf(fid, '%st', ' Load the file brainsuite_bdp_logfile for manual debuging/data checking');
    fclose(fid);
    errMsg = ['Brainsuite Process failled at the "bdp" step 01, please check the output logfile: %s ' brainsuite_bdp_logfile ];
    return
end
                                                                                       
% the main needed output is this : output_mri.dwi.RAS.correct.T1_coord.eig.nii.gz

%% Convert the eigen file to Vi and Li files
bst_progress('text', 'Calling Brainsuite "convert the eig.nii.gz to Vi and Li" ...');

% Check of the .eig.nii.gz files are generated by bdp
listing = dir(fullfile(OutPutFolder,'*.eig.nii.gz')) ;% filename can be specified with full path
if isempty(listing)
    brainsuite_bdp_logfile = fullfile(OutPutFolder,'brainsuite_bdp_logfile.txt');
    fid = fopen(brainsuite_bdp_logfile , 'w');
    fprintf(fid, '%st', cmdout_bdp);
    fprintf(fid, '%st', ' Load the file brainsuite_bdp_logfile for manual debuging/data checking');
    fclose(fid);    
    errMsg = ['It seems that there  are no files " *.eig.nii.gz " generated by the BDP .... please check the output logfile:  ' brainsuite_bdp_logfile]; return
end

eig_filename = listing.name;
output_base =  fullfile(OutPutFolder,[subjid '.tensor']); % file prefix string
%% add eig2nifti to matlab path 
eig2nifti( fullfile(OutPutFolder,eig_filename), output_base);
% ASK ANAND if we can have the source code of eig2nifti unstead of the
% binary or this p code.

%% Check the SimBio team files and code for the rest of this 
% The files generated by brainsuite:  unzip all the tensors files
gunzip( fullfile(OutPutFolder, '*.tensor.*.gz'))
% Load all the files : keep the fieldtrip vesion as it's originally used
dir(fullfile(OutPutFolder, '*.tensor.*.nii')); %subjid = 'A2016'
% ft_defaults
V1File = fullfile(OutPutFolder,[ subjid '.tensor.V1.nii']);
% V1Data = ft_read_mri(V1File); DTI{1}  = V1Data;
[V1_data, vox2ras] = in_mri_nii(V1File,1,[]); DTI{1}  = bstMri2ftMri(V1_data, vox2ras);
V2File = fullfile(OutPutFolder,[ subjid '.tensor.V2.nii']);
% V2Data = ft_read_mri(V2File); DTI{2}  = V2Data;
[V2_data, vox2ras] = in_mri_nii(V2File,1,[]); DTI{2}  = bstMri2ftMri(V2_data, vox2ras);
V3File = fullfile(OutPutFolder,[ subjid '.tensor.V3.nii']);
% V3Data = ft_read_mri(V3File); DTI{3}  = V3Data;
[V3_data, vox2ras] = in_mri_nii(V3File,1,[]);  DTI{3}  = bstMri2ftMri(V3_data, vox2ras);
L1File = fullfile(OutPutFolder,[ subjid '.tensor.L1.nii']);
% L1Data = ft_read_mri(L1File); DTI{4}   = L1Data;
[L1_data, vox2ras] = in_mri_nii(L1File,1,[]);  DTI{4}  = bstMri2ftMri(L1_data, vox2ras);
L2File = fullfile(OutPutFolder,[ subjid '.tensor.L2.nii']);
% L2Data = ft_read_mri(L2File); DTI{5}   = L2Data;
[L2_data, vox2ras] = in_mri_nii(L2File,1,[]);   DTI{5}  = bstMri2ftMri(L2_data, vox2ras);
L3File = fullfile(OutPutFolder,[ subjid '.tensor.L3.nii']);
% L3Data = ft_read_mri(L3File); DTI{6}   = L3Data;
[L3_data, vox2ras] = in_mri_nii(L3File,1,[]);   DTI{6}  = bstMri2ftMri(L3_data, vox2ras);

%% Load the mask and compute the anisotropy
% find the associated mask ==> need to have the mask
%% load the final mask / segmentation
% TodoUnzip and save the final masks to :  bst_fullfile(bst_get('BrainstormTmpDir'), 'simnibs_files',[subjid '_final_contr.nii']);
% Load the saved files on the User directory from SimNibs
simNibsMask = bst_fullfile(bst_get('BrainstormUserDir'), 'simnibs_files',subjid,[subjid '_final_contr.nii']);
if exist(simNibsMask) ~= 2 % the file exist 
        errMsg = ['It seems that there  are no mask for the subject ' subjid ', you need to run FEM mesh with SimNibs before calling the tensor']; return
end
simNibsFolder = bst_fullfile(bst_get('BrainstormUserDir'), 'simnibs_files');
if exist(simNibsFolder) ~= 7 
            errMsg = ['It seems that there  are no SimNibs folder, you need to run FEM mesh with SimNibs before calling the tensor']; return
end

% Load the final mask from SimNibs outputs 
meshGenerationMethod = 'simnibs';
switch meshGenerationMethod
    case  'simnibs'
% mri_segmented_final = ft_read_mri(simNibsMask);    
[sMri, vox2ras]  = in_mri_nii(simNibsMask);
% Re-write to the fieldtrip format
mri_segmented = bstMri2ftMri(sMri, vox2ras);
% Assuming that the data are comming from SimNibs
mri_segmented.anatomylabel = {'white' 'gray' 'csf' 'skull' 'scalp'};
% Replace the eyes by the scalp
mri_segmented.anatomy(mri_segmented.anatomy == 6) = 5;
    case 'fieldtrip'
        % todo : the finals masjks are needed
    case 'brain2mesh'
        % todo 
    case 'userDefined' % <== work on this, in the case where the user imported his mesh and masks

    otherwise
          bst_error(['FEM Mesh method is not recognized .' 10 'Check the Matlab command window for additional information.' 10], 'Generate FEM tensor', 0);       
end

%% Find the associated mesh from the data base
 if isempty(sSubject) || isempty(sSubject.iFEM)
     errMsg = 'No FEM head model available for this subject.';
 end
 % Path to the head fem model
 FemFiles = file_fullpath(sSubject.Surface(sSubject.iFEM).FileName);
 % Get the information from the head model
 fields = whos('-file',FemFiles);
 % Get number of  layers
 ivar = find(strcmpi({fields.name}, 'TissueLabels'));
 numberOfLayer =  max(fields(ivar).size);
  
if ~length(mri_segmented.anatomylabel) == numberOfLayer
          Msg = 'It seems that the mesh do not overlay the mask';
          bst_error(['The mask does not overlay the selected mesh.' 10 'Check the Matlab command window for additional information.' 10], 'Generate FEM tensor', 0);
end     
     
 % Get mesh element type
 ivar = find(strcmpi({fields.name}, 'Elements'));
 numberOfEdges =  (fields(ivar).size(2));
 if numberOfEdges == 4
     elementType = 'tetrahedron';
     bst_error(['Tetrahedron element are not supported for now.' 10 'Check the Matlab command window for additional information.' 10], 'Generate FEM tensor', 0);    
 elseif numberOfEdges == 8
     elementType = 'hexahedron';
 else
     bst_error(['Mesh element are not supported (not hexa and not tetra).' 10 'Check the Matlab command window for additional information.' 10], 'Generate FEM tensor', 0);    
 end
 
 % Get the default conductivity values : check the
 % units S/m or S/mm
 default_conductivity  = get_standard_conductivity((numberOfLayer));
%  % Get the names of the tissues
%  load(FemFiles, 'TissueLabels');
%  FemNames = TissueLabels;

%% create conductivity tensor
fprintf('...create conductivity tensor\n')
%only white matter
%owm = 0;
owm = 1;    % this is not implemented and not adapted for brainstorm for now .... only for wm 
wmIndex = 1;
% Convert the DTW to donductivity 
% This step will estimate the conductivity tensor only on the white matter
% and assighe zero every where. 
if exist('DTI','var')    
    cfg =[];
    cfg.compartments = mri_segmented.anatomylabel;
    cfg.conductivity = default_conductivity;
    if owm == 0
        [condcell, s, fail] = sb_calcTensorCond_tuch(cfg,mri_segmented,DTI{1},DTI{2},DTI{3},DTI{4},DTI{5},DTI{6});
    else
        [condcell, s, fail] = sb_calcTensorCond_tuch_wm(cfg,mri_segmented,DTI{1},DTI{2},DTI{3},DTI{4},DTI{5},DTI{6});
    end
end
% Check the value of the the tensor 
wmVoxels = find(mri_segmented.anatomy == wmIndex); 
OtherVoxels = find(mri_segmented.anatomy ~= wmIndex);
ind = randperm(length(wmTensor),1);
wmTensor = condcell{wmVoxels(ind)};
OtherTensor = condcell{OtherVoxels(ind)};

%% Load the hexa hedra mesh 

if numberOfEdges == 8
   femHead=  load(FemFiles);
else
  bst_error(['Tetrahedron element are not supported for now.' 10 'Check the Matlab command window for additional information.' 10], 'Generate FEM tensor', 0);    
end
% fpor testing generate hexa mesh :
ft_defaults
cfg        = [];
cfg.shift  = 0.3;
cfg.method = 'hexahedral';
cfg.spmversion = 'spm12';
cfg.downsample = 4;
mri_segmented.seg = mri_segmented.anatomy;
% mri_segmented_simNibs.seg = mri_segmented.anatomy;
mesh = ft_prepare_mesh(cfg,(mri_segmented));
[tetraElem,tetraNode,tetraLabel] = hex2tet(double(mesh.hex), mesh.pos, double(mesh.tissue), 3);

figure; plotmesh(tetraNode, [tetraElem tetraLabel],'x>0'); hold on; 

tetraNode2 = ft_warp_apply(inv(mri_segmented.transform), tetraNode, 'homogeneous');

[tetraNode3, Transf] = cs_convert(sMri, src, dest, P)

figure; plotmesh(tetraNode2, [tetraElem tetraLabel],'x>0'); hold on; 


% Use the hexa without fieldtrip
cfg.background = 0; % ===> it seems that the transfomation is not apllied here ... MRI coordinate system
mesh2 = build_mesh_hexahedral(cfg,(mri_segmented));
% Correct the orientation 
mesh2.pos_alignedToMri  = mesh2.pos;
mesh2.pos_alignedToMri(:,[1 2])  = mesh2.pos(:,[2 1]);

[tetraElem2,tetraNode2,tetraLabel2] = hex2tet(double(mesh2.hex), mesh2.pos_alignedToMri, double(mesh2.labels), 3);
% figure; plotmesh(tetraNode2, [tetraElem2 tetraLabel2],'x>0'); hold on; 

hFig = figure;
hAx         = axes('Parent',hFig);
fixedVolume = mri_segmented.seg;
centerFixed = size(fixedVolume)/2;
slice(hAx,double(fixedVolume),centerFixed(1),centerFixed(2),centerFixed(3));
shading(hAx,'interp');
set(hAx,'Xgrid','on','YGrid','on','ZGrid','on');
set(hFig,'Colormap',colormap('gray'));
hold on;
% plotmesh([tetraNode2(:,2) tetraNode2(:,1) tetraNode2(:,3)], [tetraElem2 tetraLabel2],'x>0'); hold on; 
% plotmesh(tetraNode2, [tetraElem2 tetraLabel2],'x>0'); hold on; 
plotmesh(tetraNode2, [tetraElem tetraLabel],'x>0'); hold on; 
xlabel('X'); ylabel('Y');  zlabel('Z'); 

 

% assigne for all element of the model
if 0
    % modifications for white matter anisotropy
    idx_wma         = find(sum(gm_wm_tensors,1)~=0 & labels' == wmIndex ); % wmIndex is the label of white matter
    nonzero_tensors = gm_wm_tensors(:,idx_wma);
    tensors         = zeros(9, numel(conductivity)+size(nonzero_tensors,2));
    tensors(:, 1:numel(conductivity))     = (conductivity' * [1 0 0 0 1 0 0 0 1])';
    tensors(:, numel(conductivity)+1:end) = nonzero_tensors;
    labels(idx_wma) = (1:numel(idx_wma)) + size(conductivity,2);
end
%% Load the associated mesh and assigne the conductivity
%%%% ---- PROBLEM when we assigne the tensors need to remove the negative
%%%% index or translate the MRI to new coordinate system
if 0
                        % Get current subject
                        iSubject
                    iAnatomy
                    sSubject
                    
                       
    
    %% Assigne the conducivity to the element of the mesh
    if exist('condcell','var')
        %create brain mask
        mask = mri_segmented;
        if owm == 0
            mask.anatomy((mask.anatomy~=5) & (mask.anatomy~=6))=0;
            mask.anatomy((mask.anatomy==5) | (mask.anatomy==6))=1;
        else
            mask.anatomy(mask.anatomy~=wmIndex)=0;
            mask.anatomy(mask.anatomy==wmIndex)=1;
        end
        mask.dim=size(mask.anatomy);
        
        if owm == 0
            tensors = sb_assiTensorCond(mask,mesh.nodes,mesh.elements,mesh.labels,condcell);
            %         [tensor_mcp, maxCond] = sb_assiTensorCond_mcp(mask,mesh.nodes,mesh.elements,condcell);
        else
            tensors = sb_assiTensorCond_wm(mask,mesh.HexaMesh.pos,mesh.HexaMesh.hex,mesh.HexaMesh.tissue,condcell);
            %         [tensor_mcp, maxCond] = sb_assiTensorCond_mcp(mask,mesh.nodes,mesh.elements,condcell);
        end
        fprintf('mn = %.4f and mx = %.4f\n',min(min(tensors)),max(max(tensors)));
    end
end

%% Save the tenso as a file with 1:nb element, coordinate of the centroide aand the the valu of the 9 tensors in the ctf coordinate system
% elemIndex  Centroide coordinates Tensor value 
% to display ... we need just to overlay the tensor file to the mesh file
% in the data base 

% Display tensors 


end