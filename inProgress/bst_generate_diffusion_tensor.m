function  [DTI, errMsg] = bst_generate_diffusion_tensor(iSubject)

% tree_callbacks.m  10227  gui_component('MenuItem', jPopup, [], 'Generate FEM conductivity tensor', IconLoader.ICON_FEM, [], @(h,ev)generate_conductivity_tensor(iSubject, iAnatomy));

%% Check if brainsuite is available otherwise ask user to install it and add link to that
disp('Check brainsite instalation ...')
[status,cmdout] = system('bdp --version');
if status ~= 0
    errMsg = [('Brainsuite is not installed in you computer, please visite : http://forums.brainsuite.org/download/ and xxx brainstorm tutorial')];
    return
end
disp('Check simnibs instalation ...')
[status,cmdout] = system('headreco --h');
if status ~= 0
    errMsg = ['The FEM tensor generation works only with SimNibs mesh and mask :  xxx brainstorm tutorial'];
    return
end

bst_progress('start', 'Generate Conductivity tensor', 'Generate FEM Conductivity tensor...');

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

%% Check the presence of the data on the subject
% Check the presence of the anatomy
if isempty(sSubject.Anatomy)
    errMsg = ['There no anatomy data on this subject'];
    return
end

if length(sSubject.Anatomy) < 2
    errMsg = ['At leadt the T1 and the segmentation data (tissue) should be available'];
    return
end

% TODO : Francois ==> add an exact ID to identify the tissus anatomy maybe : iTissu on the sSubject sturct
if length(sSubject.Anatomy) == 2
    disp('It seems that only T1 is used for the mesh and segmentation ... assuming the second file as tissue');
    iT1Mri = sSubject.iAnatomy;
    iTissues = 2;
elseif length(sSubject.Anatomy) == 3
    disp('It seems that T1 and T2 are used for the mesh and segmentation ... assuming the third file as tissue');
    iT1Mri = sSubject.iAnatomy;
    iTissues = 3;
end

% Check the tissu file ==> the comment should contains tissues
if strfind(sSubject.Anatomy(iTissues).Comment,  'tissues') >0
    disp(['High chance that the index  ' num2str(iTissues) ' is the tissues'] )
else
    disp(['High chance that the index  ' num2str(iTissues) ' is not the tissues'] )
end

% Test the presence of the MESH
if isempty(sSubject.iFEM)
    errMsg = ['There no FEM mesh on this subject'];
    return
end
if isempty(sSubject.Surface(sSubject.iFEM))
    errMsg = ['There no FEM mesh on this subject'];
    return
end

%%  Ask the user for the DTI data path : The asociated bvec and bval should be in the same folder
% ODO : ftadel : find a better way to link thee files to the subject
%% ===== SELECT DTI FILE =====
% If MRI file to load was not defined : open a dialog box to select it
% Get last used directories
LastUsedDirs = bst_get('LastUsedDirs');
% Get last used format
DefaultFormats = bst_get('DefaultFormats');
% Get DTI file
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
[filepath,filename,ext]  =  fileparts(DtiFile);
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
% check the name of the data ==> should match to be arranged by subject
if ~strfind(bvalFileData.name,filename(1:4))
    errMsg = [('We can not match the  "*.bval" file to the MRI, please associate the same name for all the data belonging to the same subject')];
    return
end

if ~strfind(bvecFileData.name,filename(1:4))
    errMsg = [('We can not match the  "*.bvec" file to the MRI, please associate the same name for all the data belonging to the same subject')];
    return
end
% At this level we assume that all the data are present on the folder

%% If all the tests are OK .... We can start
% ==== Empty temporary folder ====
gui_brainstorm('EmptyTempFolder');
% Create temporary folder for segmentation files
brainsuiteDir = bst_fullfile(bst_get('BrainstormTmpDir'), 'brainsuite');
mkdir(brainsuiteDir);
% === SAVE MRI AS NII ===
bst_progress('text', 'Exporting MRI : T1...');
% Save MRI in .nii format <==== be sure that the selected mri is the T1
T1File = file_fullpath(sSubject.Anatomy(iT1Mri).FileName);
subjid = strrep(sSubject.Name, '@', '');
T1Nii = bst_fullfile(brainsuiteDir, [subjid 'T1.nii']);
sMriT1 = in_mri_bst(T1File);
out_mri_nii(sMriT1, T1Nii);
% === SAVE TISSUES AS NII ===
bst_progress('text', 'Exporting MRI : Tissues/masks...');
MaskFile = file_fullpath(sSubject.Anatomy(iTissues).FileName);
MaskNii = bst_fullfile(brainsuiteDir, [subjid 'Mask.nii']);
sMriMask = in_mri_bst(MaskFile);
out_mri_nii(sMriMask, MaskNii);

%% Call the  brainsuite pipline
bst_progress('text', 'Calling Brainsuite ...');
OutPutFolder = brainsuiteDir;
T1Filename = T1Nii;
dtiFilename = DtiFile;
bvecFilename = fullfile(bvecFileData.folder,bvecFileData.name);
bvalFilename = fullfile(bvalFileData.folder,bvalFileData.name);
% DTI = bst_brainsuite_generate_dti_tensor(T1Filename, dtiFilename, bvecFilename, bvalFilename, OutPutFolder);

DTI = bst_brainsuite_generate_dti_tensor(T1Filename, dtiFilename, bvecFilename, bvalFilename, subjid, OutPutFolder);


end
%
% %% Find the associated mesh from the data base
% if isempty(sSubject) || isempty(sSubject.iFEM)
%     errMsg = 'No FEM head model available for this subject.';
% end
% % Path to the head fem model
% FemFiles = file_fullpath(sSubject.Surface(sSubject.iFEM).FileName);
% % Get the information from the head model
% fields = whos('-file',FemFiles);
% % Get number of  layers
% ivar = find(strcmpi({fields.name}, 'TissueLabels'));
% numberOfLayer =  max(fields(ivar).size);
% if ~length(mri_segmented.anatomylabel) == numberOfLayer
%     Msg = 'It seems that the mesh do not overlay the mask';
%     bst_error(['The mask does not overlay the selected mesh.' 10 'Check the Matlab command window for additional information.' 10], 'Generate FEM tensor', 0);
% end
% % Get mesh element type
% ivar = find(strcmpi({fields.name}, 'Elements'));
% numberOfEdges =  (fields(ivar).size(2));
%
% %% Load the mesh
% femHead=  load(FemFiles);
%
% if numberOfEdges == 8
%     %    quick check the mesh
%     [tetraElem,tetraNode,tetraLabel] = hex2tet(double(femHead.Elements), femHead.Vertices, double(femHead.Tissue), 3);
%     figure; plotmesh(tetraNode, [tetraElem tetraLabel],'x>0'); hold on;
%     title('SCS/CTF coordinates systems')
% elseif numberOfEdges == 4
%     figure; plotmesh(femHead.Vertices,[femHead.Elements femHead.Tissue] ,'x>0');
%     title('SCS/CTF coordinates systems')
% else
%     bst_error(['Tetrahedron element are not supported for now.' 10 'Check the Matlab command window for additional information.' 10], 'Generate FEM tensor', 0);
% end
%
%
%
% %% SECTION 2 : Convert DWI to Conductivity
% % Compute the tensor for the whole model
% %% Build the tensor for the isotropic tissues
% cfg = [];
% default_iso_conductivity = get_standard_conductivity((numberOfLayer)) ;
% % default_iso_conductivity = (default_iso_conductivity/max(default_iso_conductivity))/1000;
% conversion_m2mm = 1000;
% default_iso_conductivity = default_iso_conductivity/conversion_m2mm;
%
% cfg.elem = [femHead.Elements femHead.Tissue];
% cfg.node = femHead.Vertices;
% cfg.conductivity = default_iso_conductivity;
% cfg = bst_generate_tensor_on_elem(cfg);
% % Tensors isotropic
% cfg.IsotropicTensor = cfg.eigen;
% IsotropicTensor = cfg.eigen;
% cfg.eigen = IsotropicTensor;
% % Display the tensors
% cfg.indElem = 1:50000:length(cfg.elem);
% cfg.ellipse = 1;
% cfg.arrow = 0;
% bst_display_tensor_as_ellipse(cfg)
% view([90 0 0])
%
% figure;
% plotmesh(cfg.node, cfg.elem)
% %% Assigne anisotropy ===> convert from DTI to conductivity
% % There are many methds discussed on the litterature, but most of them are
% % related to the Tuch approach.
% % Here we implement these methods :
% % 1- Direc approach (Tuch formulation) [without the isotropic value]
% % 2- Direct approach with volume constraint/normalized [normalized with the isotropic value]
% % 3- Volume normalized  [don't use the Tuch relation]
% % 4- Wang constrain and volume constrain
% % 5- Artificial conductivity with  ratio
%
% %% Normalized volume approach
%
% % This maybe done during the process of FEM computation since the value of
% % the conductivity may change according to the users...Discuss with  ftadel
% %% Call the mapping from voxel to mesh and coordinates changes
% [V1rot,V2rot,V3rot, L1a, L2a, L3a ] = convert3Dtensors(sMriMask, femHead, DTI);
% % the output of this function is the eigen value and vector interpolated to
% % the centroide of each fem element and the vector are reoriented to the
% % SCS coordinate system
% % aniso_conductivity= [];
% % iso_conductivity = default_iso_conductivity(1);  % get this valu from the input
% % factor = 1;
% % L1a = L1a * factor;
% % L2a = L2a * factor;
% % L3a = L3a * factor;
% % ind = 1;
% % fail = [];
% % for globalIndex = 1 : length(femHead.Tissue)
% %     if femHead.Tissue(globalIndex) == 1 % 1 is the label of the wm
% %         L=diag([L1a(ind), L2a(ind), L3a(ind)]);
% %         V=[V1rot(ind,:)',V2rot(ind,:)',V3rot(ind,:)'];
% %         T = V*L*V';
% %         if sum(sum(L)) == 0 % useful in the case where BDP fails or it's not part of the wm mask
% %             fail = [fail ind];
% %             aniso_conductivity(:,:,ind) = diag([iso_conductivity,iso_conductivity,iso_conductivity]);
% %             eigen.eigen_vector{ind} = IsotropicTensor.eigen_vector{globalIndex};
% %             eigen.eigen_value{ind} =  IsotropicTensor.eigen_value{globalIndex};
% %         else % apply the volume approcah
% %                         aniso_conductivity(:,:,ind) = iso_conductivity * [(L1a(ind)*L2a(ind)*L3a(ind))^(-1/3)]*T;
% %                         [v,l]=eig(aniso_conductivity(:,:,ind));
% %                         eigen.eigen_vector{ind} = v;
% %                         eigen.eigen_value{ind} = l;
% % %             eigen.eigen_vector{ind} = V;
% % %             eigen.eigen_value{ind} = L;
% %         end
% %     anisotropicTensor.eigen_vector{globalIndex} = eigen.eigen_vector{ind} ;
% %     anisotropicTensor.eigen_value{globalIndex} = eigen.eigen_value{ind};
% %     ind = ind +1;
% %     else
% %     anisotropicTensor.eigen_vector{globalIndex} = IsotropicTensor.eigen_vector{globalIndex}  ;
% %     anisotropicTensor.eigen_value{globalIndex} = IsotropicTensor.eigen_value{globalIndex} ;
% %     end
% % end
%
% [aniso_conductivity, anisotropicTensor ] = bst_compute_anisotropy_tensors(femHead, default_iso_conductivity, IsotropicTensor, L1a,L2a,L3a,V1rot,V2rot,V3rot);
%
% % check and display
% % cfg.conductivity_tensor3x3 = aniso_conductivity/display_factor;
% cfg.indElem = find(cfg.elem(cfg.indElem ,end)==1);
% cfg.indElem = 1:500:length(cfg.elem);
% % cfg.indElem = cfg.indElem (1:10000:end);
% cfg.eigen = anisotropicTensor;
% cfg.noCsf  = 1; % do not display CSF
% cfg.ellipse = 0;
% cfg.arrow = 0;
% cfg.plotMesh = 1;
% bst_display_tensor_as_ellipse(cfg)
% view([0 90 0])
%
%
% %% Interpolate the tensor to the face for displaying
% z0 = mean(femHead.Vertices(:,3));
% z0 = 0.038;
% plane=[min(femHead.Vertices(:,1)) min(femHead.Vertices(:,2)) (z0 + z0/2)
%            min(femHead.Vertices(:,1)) max(femHead.Vertices(:,2)) (z0 + z0/2)
%            max(femHead.Vertices(:,1)) min(femHead.Vertices(:,2)) (z0 + z0/2)];
%
% y0 = 0.022;
% plane=[min(femHead.Vertices(:,1)) y0 min(femHead.Vertices(:,3))
%            min(femHead.Vertices(:,1)) y0 max(femHead.Vertices(:,3))
%            max(femHead.Vertices(:,1)) y0 min(femHead.Vertices(:,3))];
%
%
% x0 = 0.0;
% plane=[x0  min(femHead.Vertices(:,2))  min(femHead.Vertices(:,3))
%            x0  max(femHead.Vertices(:,2)) max(femHead.Vertices(:,3))
%            x0  min(femHead.Vertices(:,2)) min(femHead.Vertices(:,3))];
%
% % run qmeshcut to get the cross-section information at z=mean(node(:,1))
% % use the x-coordinates as the nodal values
%
% cfg.elem = cfg.elem(cfg.elem(:,5) == 1, :);
%
% [cutpos,cutvalue,facedata,elemid]= ...
%                     qmeshcut(cfg.elem(:,1:4),cfg.node,zeros(length(cfg.node),1),plane);
%
% figure;
% plotmesh(femHead.Vertices,[femHead.Elements femHead.Tissue], 'facealpha', 0.3,'edgecolor','none');
% hold on
% plotmesh(femHead.Vertices,[femHead.Elements(elemid,:) femHead.Tissue(elemid,:)]);
% view([0 90 0])
%
% % cfg.elem = cfg.elem(cfg.elem(:,5) == 1, :);
%
% cfg.indElem = elemid(1:1:end);
% % cfg.indElem = find(cfg.elem(cfg.indElem ,end)==1);
% cfg.eigen = anisotropicTensor;
% cfg.noCsf  = 1; % do not display CSF
% cfg.ellipse = 1;
% cfg.arrow = 0;
% cfg.plotMesh = 0;
% bst_display_tensor_as_ellipse(cfg)
% view([0 -90 0])
% axis([ -0.063 0.095 ...
%          min(cfg.node(:,2))  max(cfg.node(:,2)) ...
%          -0.02  0.105])
% hold on
% plotmesh(femHead.Vertices,femHead.Elements(elemid,:),'facealpha',0.2,'edgecolor','none');
% view([0 90 0])
%
% hold on
% plotmesh(femHead.Vertices,femHead.Elements(elemid,:),'facealpha',0.2);
%
% view([0 0 90])
%
%
% figure;
% plotmesh(femHead.Vertices,femHead.Elements,'facealpha',0.2,'edgecolor','none');
% view([0 0 90])
%
%
% figure;
% plotmesh(femHead.Vertices,femHead.Elements(elemid,:),'facealpha',0.2);
% view([0 0 90])
%
% %% TODO : PBLM with the conversion ... check the value of the eigen values & vectors
%
% % get the index of the iso and aniso
% wmAniso = ones(length(V1rot),1);
% wmAniso(fail) = 0;
% allTensorIndex = zeros(length(femHead.Tissue),1);
% allTensorIndex(femHead.Tissue ==1) = wmAniso; % 1 aniso, 0 iso
%
% size(allTensorIndex)
% size(femHead.Tissue )
%
%
% % replace the isotropic conductivity by the computed values from DWI
%
% %% --------------------------------------------------------------------------------------------------------
% %To be tested later
% if 0
%
%     %% Load the mask and compute the tensors
%     % load the final mask / segmentation
%     bst_progress('text', 'Exporting MRI : Tissues/masks...');
%     meshGenerationMethod = 'simnibs';
%     switch meshGenerationMethod
%         case  'simnibs'
%             % mri_segmented_final = ft_read_mri(simNibsMask);
%             [sMri, vox2ras]  = in_mri(MaskNii);
%             % Re-write to the fieldtrip format
%             %         sMriMask ;         mri_segmented = bstMri2ftMri(sMriMask, sMriMask.InitTransf{2});
%             mri_segmented = bstMri2ftMri(sMri, vox2ras);
%             % Assuming that the data are comming from SimNibs
%             mri_segmented.anatomylabel = {'white' 'gray' 'csf' 'skull' 'scalp'};
%             % Replace the eyes by the scalp
%             mri_segmented.anatomy(mri_segmented.anatomy == 6) = 5;
%         case 'fieldtrip'
%             % todo : the finals masjks are needed
%         case 'brain2mesh'
%             % todo ==> get the final masks from the segmentations and aligned
%             % to the right coordinates systemes
%         case 'userDefined' % <== work on this, in the case where the user imported his mesh and masks
%
%         otherwise
%             bst_error(['FEM Mesh method is not recognized .' 10 'Check the Matlab command window for additional information.' 10], 'Generate FEM tensor', 0);
%     end
%     % TODO : Get the default conductivity values : check the units S/m or S/mm
%     % TODO : add a panel to control these values  : Specify the anisotropy
%     % values and ask user to validate these values
%     default_conductivity = get_standard_conductivity((numberOfLayer));
%     %  % Get the names of the tissues
%     %% Create conductivity tensor
%     bst_progress('text', 'Calling Brainsuite  create conductivity tensor...');
%
%     fprintf('...create conductivity tensor\n')
%     %only white matter
%     owm = 1;    % this is not implemented and not adapted for brainstorm for now .... only for wm
%     wmIndex = 1;
%     % Convert the DTW to donductivity
%     % This step will estimate the conductivity tensor only on the white matter
%     % and assighe zero every where.
%     if exist('DTI','var')
%         cfg =[];
%         cfg.compartments = mri_segmented.anatomylabel;
%         cfg.conductivity = default_conductivity;
%         if owm == 0
%             [condcell, s, fail] = sb_calcTensorCond_tuch(cfg,mri_segmented,DTI{1},DTI{2},DTI{3},DTI{4},DTI{5},DTI{6});
%         else
%             [condcell, s, fail] = sb_calcTensorCond_tuch_wm(cfg,mri_segmented,DTI{1},DTI{2},DTI{3},DTI{4},DTI{5},DTI{6});
%         end
%     end
%     % Check the value of the the tensor
%     wmVoxels = find(mri_segmented.anatomy == wmIndex);
%     OtherVoxels = find(mri_segmented.anatomy ~= wmIndex);
%     ind = randperm(length(wmVoxels),1);
%     wmTensor = condcell{wmVoxels(ind)};
%     OtherTensor = condcell{OtherVoxels(ind)};
%
%
%
%     %% Convert the coordinate of the mesh from SCS to the VOXEL/or RAS
%     % The nodes should be in the same coordinates as the output of the brainsuite
%
%     % Apply the transoformation from scs to vox
%     [femHead.Vertices_voxel, Transf] = cs_convert(sMriMask,  'scs',  'voxel', femHead.Vertices);
%     [femHead.Vertices_ras, Transf] = cs_convert(sMriMask,  'scs',  'world', femHead.Vertices);
%     [femHead.Vertices_mri, Transf] = cs_convert(sMriMask,  'scs',  'mri', femHead.Vertices);
%
%
%     figure;
%     subplot(2,2,1)
%     plotmesh(femHead.Vertices, [femHead.Elements femHead.Tissue],'edgecolor','none'); title('SCS ');xlabel('x');ylabel('y');zlabel('z');grid on; grid minor;
%     subplot(2,2,2)
%     plotmesh(femHead.Vertices_voxel, [femHead.Elements femHead.Tissue],'edgecolor','none'); title('VOXEL');xlabel('x');ylabel('y');zlabel('z');grid on; grid minor;
%     subplot(2,2,3)
%     plotmesh(femHead.Vertices_ras, [femHead.Elements femHead.Tissue],'edgecolor','none'); title('RAS ');xlabel('x');ylabel('y');zlabel('z');grid on; grid minor;
%     subplot(2,2,4)
%     plotmesh(femHead.Vertices_mri, [femHead.Elements femHead.Tissue],'edgecolor','none'); title('MRI ');xlabel('x');ylabel('y');zlabel('z');grid on; grid minor;
%
%     % save('sessionGenerateTensors21Avril')
%     [V1rot,V2rot,V3rot, L1a, L2a, L3a ]=convert3Dtensors(sMriMask, femHead, DTI);
%     [V1rot,V2rot,V3rot, L1a, L2a, L3a ];
%
%     V=[V1rot,V2rot,V3rot];
%
%     % L=diag([L1,L2,L3]);
%     %  T = VLV'
%     aniso_conductivity= [];
%     iso_conductivity = 0.33;  % get this valu from the input
%     fail = [];
%     for ind = 1 : length(V1rot)
%         L=diag([L1a(ind), L2a(ind), L3a(ind)]);
%         V=[V1rot(ind,:)',V2rot(ind,:)',V3rot(ind,:)'];
%         T = V*L*V';
%         %     tensor(ind) = T;
%         if sum(sum(L)) == 0
%             fail = [fail ind];
%             aniso_conductivity(:,:,ind) = diag([iso_conductivity,iso_conductivity,iso_conductivity]);
%         else
%             aniso_conductivity(:,:,ind) = iso_conductivity * [(L1a(ind)*L2a(ind)*L3a(ind))^(-1/3)]*T;
%         end
%     end
%     cfg.node = femHead.Vertices;
%     cfg.elem = [femHead.Elements(femHead.Tissue==1,:) femHead.Tissue(femHead.Tissue==1)];
%     cfg.elem_centroide = bst_generate_centroide_on_elem(cfg.node,cfg.elem);
%     % figure; plotmesh(cfg.node, cfg.elem(1,:),'facealpha',0.3); hold on;plotmesh(cfg.elem_centroide(1,:),'r.')
%     cfg.indElem = fail(1:1:end);
%     cfg.ellipse = 0;
%     cfg.arrow =1;
%     cfg.conductivity_tensor3x3 = aniso_conductivity/300;
%     bst_display_tensor_as_ellipse(cfg)
%     view([0 90 0])
%
%     %% Load the associated mesh and assigne the conductivity
%     %%%% ---- PROBLEM when we assigne the tensors need to remove the negative
%     %%%% index or translate the MRI to new coordinate system
%     if 0
%
%
%         checkMRI =0;
%         if checkMRI == 1
%             hFig = figure;
%             hAx         = axes('Parent',hFig);
%             fixedVolume = mri_segmented.anatomy;
%             centerFixed = size(fixedVolume)/2;
%             slice(hAx,double(fixedVolume),centerFixed(1),centerFixed(2),centerFixed(3));
%             shading(hAx,'interp');
%             set(hAx,'Xgrid','on','YGrid','on','ZGrid','on');
%             set(hFig,'Colormap',colormap('gray'));
%             hold on;
%             plotmesh(femHead.Vertices_voxel,  [femHead.Elements femHead.Tissue],'x>0'); hold on;
%             xlabel('X'); ylabel('Y');  zlabel('Z');
%         end
%
%
%         %% Assigne the conducivity to the element of the mesh
%         if exist('condcell','var')
%             %create brain mask
%             mask = mri_segmented;
%             if owm == 0
%                 mask.anatomy((mask.anatomy~=5) & (mask.anatomy~=6))=0;
%                 mask.anatomy((mask.anatomy==5) | (mask.anatomy==6))=1;
%             else
%                 mask.anatomy(mask.anatomy~=wmIndex)=0;
%                 mask.anatomy(mask.anatomy==wmIndex)=1;
%             end
%             mask.dim=size(mask.anatomy);
%
%             if owm == 0
%                 gm_wm_tensors = sb_assiTensorCond(mask,mesh.nodes,mesh.elements,mesh.labels,condcell);
%                 %         [tensor_mcp, maxCond] = sb_assiTensorCond_mcp(mask,mesh.nodes,mesh.elements,condcell);
%             else % femHead.Vertices_voxel
%                 [gm_wm_tensors,  indexPos] = sb_assiTensorCond_wm(mask,femHead.Vertices_voxel,femHead.Elements,femHead.Tissue,condcell);
%                 %         [tensor_mcp, maxCond] = sb_assiTensorCond_mcp(mask,mesh.nodes,mesh.elements,condcell);
%             end
%             fprintf('mn = %.4f and mx = %.4f\n',min(min(gm_wm_tensors)),max(max(gm_wm_tensors)));
%         end
%     end
%
%     % assigne for all element of the model
%     if 0
%         conductivity = default_conductivity;
%         labels = double(femHead.Tissue');
%         % modifications for white matter anisotropy
%         idx_wma  = find(sum(gm_wm_tensors,1)~=0 &  labels == wmIndex ); % wmIndex is the label of white matter
%         nonzero_tensors = gm_wm_tensors(:,idx_wma);
%         %     tensors         =   zeros(9, numel(conductivity) + size(nonzero_tensors,2));
%         tensors         =   zeros(numel(labels),9) ;
%
%         iso_conductivity =  (conductivity' * [1 0 0 0 1 0 0 0 1])
%         aniso_conductivity = nonzero_tensors';
%         tensors = iso_conductivity(labels,:);
%         tensors(idx_wma,:) = aniso_conductivity;
%         %
%         %     tensors(:, 1:numel(conductivity))     = (conductivity' * [1 0 0 0 1 0 0 0 1])';
%         %     tensors(:, numel(conductivity) + 1 : end) = nonzero_tensors;
%         %     labels(idx_wma) = (1:numel(idx_wma)) + size(conductivity,2);
%     end
%
%     %% Dsiplay the mesh and the tensors :
%     cfg = [];
%     % tensors = gm_wm_tensors';
%     cfg.node = double(femHead.Vertices_voxel);
%     cfg.elem = double([femHead.Elements femHead.Tissue]);
%
%     % conductivity_tensor3x3: [3󫢬7822 double]
%     for ind = 1 : length(cfg.elem)
%         conductivity_tensor3x3(1:3,1:3,ind) =  [tensors(ind,1:3); tensors(ind,4:6); tensors(ind,7:9)];
%     end
%     cfg.conductivity_tensor3x3 = conductivity_tensor3x3;
%
%     % comput centroide of the elemts
%     cfg.elem_centroide = bst_generate_centroide_on_elem(cfg.node,cfg.elem);
%     cfg.plotMesh =0;
%     el = find(cfg.elem(:,end) == 1);
%     % cfg.elem = cfg.elem(el,:) ;
%     cfg.indElem = el(1:3:end),% 1:50:length(cfg.elem) ;%el(1:50:end) ;% 4500:6000;
%
%     bst_display_tensor_as_ellipse(cfg)
%
%     cfg = bst_display_tensor_as_mainEigenVector(cfg)
%
%     figure;
%     plotmesh(cfg.node,cfg.elem,'facealpha', 0.2);
%     hold on;
%     plotmesh(cfg.elem_centroide(1,:),'ro');
%
%
%
%     %% Reshape the Tensor and apply the transformation on the coordinates
%     tensors_matrix = reshape(tensors',3,[]);
%     tensors_matrix = tensors_matrix';
%     % clc % checking
%     % tensors_matrix(1:6,:)
%     % tensors(1:3,:)
%     %% Apply the transfert coordinates
%     [tensors_matrix_scs, Transf] = cs_convert(sMriMask,  'voxel',  'scs', tensors_matrix);
%
%
%     %% Save the tenso as a file with 1:nb element, coordinate of the centroide and the the valu of the 9 tensors in the ctf coordinate system
%     % elemIndex  Centroide coordinates Tensor value
%     % to display ... we need just to overlay the tensor file to the mesh file
%     % in the data base
%
%     % Display tensors
% end
%
% end