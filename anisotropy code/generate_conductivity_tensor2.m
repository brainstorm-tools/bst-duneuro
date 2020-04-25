function  errMsg = generate_conductivity_tensor2(iSubject, iAnatomy)

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

% maybe you need to save also the mask from SimNibs and laod it to
bst_progress('start', 'Generate Conductivity tensor', 'Generate FEM Conductivity tensor...');
iSubject
iAnatomy
% sSubject
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
   disp('It seems that only T1 is used for the mesh and segmentation ... assuming the second file as tissu');    
   iT1Mri = sSubject.iAnatomy;
   iTissues = 2;
elseif length(sSubject.Anatomy) == 3
   disp('It seems that T1 and T2 are used for the mesh and segmentation ... assuming the third file as tissu');       
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
%% Load the mask and compute the anisotropy
% load the final mask / segmentation
bst_progress('text', 'Exporting MRI : Tissues/masks...');
% MaskFile = file_fullpath(sSubject.Anatomy(iTissues).FileName); % could be
% useful in the case were we can assigne the conductivity in the right
% ctf/scs coordinate systems coordinate Load the final mask from SimNibs outputs
meshGenerationMethod = 'simnibs';
switch meshGenerationMethod
    case  'simnibs'
        % mri_segmented_final = ft_read_mri(simNibsMask);
        [sMri, vox2ras]  = in_mri(MaskNii);
        % Re-write to the fieldtrip format
%         sMriMask ;         mri_segmented = bstMri2ftMri(sMriMask, sMriMask.InitTransf{2});
        mri_segmented = bstMri2ftMri(sMri, vox2ras);
        % Assuming that the data are comming from SimNibs
        mri_segmented.anatomylabel = {'white' 'gray' 'csf' 'skull' 'scalp'};
        % Replace the eyes by the scalp
        mri_segmented.anatomy(mri_segmented.anatomy == 6) = 5;
    case 'fieldtrip'
        % todo : the finals masjks are needed
    case 'brain2mesh'
        % todo ==> get the final masks from the segmentations and aligned
        % to the right coordinates systemes
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
%      bst_error(['Tetrahedron element are not supported for now.' 10 'Check the Matlab command window for additional information.' 10], 'Generate FEM tensor', 0);    
 elseif numberOfEdges == 8
     elementType = 'hexahedron';
 else
     bst_error(['Mesh element are not supported (not hexa and not tetra).' 10 'Check the Matlab command window for additional information.' 10], 'Generate FEM tensor', 0);    
 end
 
%% Load the mesh 
if numberOfEdges == 8
    femHead=  load(FemFiles);
    %    quick check the mesh
    [tetraElem,tetraNode,tetraLabel] = hex2tet(double(femHead.Elements), femHead.Vertices, double(femHead.Tissue), 3);
    figure; plotmesh(tetraNode, [tetraElem tetraLabel],'x>0'); hold on;
    title('model from bst')
elseif numberOfEdges == 4
    femHead=  load(FemFiles);
    figure; plotmesh(femHead.Vertices,[femHead.Elements femHead.Tissue] ,'x>0');
    title('SCS/CTF coordinates systems')
else
    bst_error(['Tetrahedron element are not supported for now.' 10 'Check the Matlab command window for additional information.' 10], 'Generate FEM tensor', 0);
end

%% Convert the coordinate of the mesh from SCS to the VOXEL/or RAS ... need to check       
% The nodes should be in the same coordinates as the output of the brainsuite 

% Apply the transoformation from scs to vox
[femHead.Vertices_voxel, Transf] = cs_convert(sMriMask,  'scs',  'voxel', femHead.Vertices);
[femHead.Vertices_ras, Transf] = cs_convert(sMriMask,  'scs',  'world', femHead.Vertices);
[femHead.Vertices_mri, Transf] = cs_convert(sMriMask,  'scs',  'mri', femHead.Vertices);


figure;
subplot(2,2,1)
plotmesh(femHead.Vertices, [femHead.Elements femHead.Tissue],'edgecolor','none'); title('SCS ');xlabel('x');ylabel('y');zlabel('z');grid on; grid minor;
subplot(2,2,2)
plotmesh(femHead.Vertices_voxel, [femHead.Elements femHead.Tissue],'edgecolor','none'); title('VOXEL');xlabel('x');ylabel('y');zlabel('z');grid on; grid minor;
subplot(2,2,3)
plotmesh(femHead.Vertices_ras, [femHead.Elements femHead.Tissue],'edgecolor','none'); title('RAS ');xlabel('x');ylabel('y');zlabel('z');grid on; grid minor;
subplot(2,2,4)
plotmesh(femHead.Vertices_mri, [femHead.Elements femHead.Tissue],'edgecolor','none'); title('MRI ');xlabel('x');ylabel('y');zlabel('z');grid on; grid minor;

% femHead.Vertices_voxel(find(femHead.Vertices_voxel(:,3)<0),3) = (femHead.Vertices_voxel(find(femHead.Vertices_voxel(:,3)<0),3)) + 1
% femHead.Vertices_voxel(find(femHead.Vertices_voxel(:,3)==0),3) = (femHead.Vertices_voxel(find(femHead.Vertices_voxel(:,3)==0),3)) + 1
% 
% find(femHead.Vertices_voxel(:,3)==0)
% find(femHead.Vertices_voxel(:,2)==0)
% find(femHead.Vertices_voxel(:,1)==0)


%% Extract the EigenValue and vector for each voxel
% apply the anand model
% eport head model from brainstorm a206 ds 2 noshift
[Vertices_ras, Transf] = cs_convert(sMriMask,  'scs',  'voxel', femHead.Vertices);

cfg2.elem_wm = [femHead.Elements(femHead.Tissue == 1, :) femHead.Tissue(femHead.Tissue==1)];
cfg2.node = Vertices_ras;
elem_centroide_wm = bst_generate_centroide_on_elem(cfg2.node,cfg2.elem_wm);
cfg2.elem_centroide_wm = elem_centroide_wm;
if femHead.Elements >5
% display the WM
[tetraElem,tetraNode,tetraLabel] = hex2tet(cfg2.elem_wm(:,1:end-1), cfg2.node, cfg2.elem_wm(:,end), 3);
figure; plotmesh(tetraNode,  [tetraElem tetraLabel],'x>100','facealpha',0.3); hold on;
else
    figure; plotmesh(cfg2.node, cfg2.elem_wm,'x>100','facealpha',0.3); hold on;
end
hold on
plotmesh(cfg2.elem_centroide_wm,'r.')

SZ=size(DTI{1}.anatomy);
res=DTI{1}.hdr.dim.pixdim(2:4); 

Vertices_ras = elem_centroide_wm;
clear V1a;
V1a(:,1) = interpn(DTI{1}.anatomy(:,:,:,1), Vertices_ras(:,1)/res(1),Vertices_ras(:,2)/res(2),Vertices_ras(:,3)/res(3));
V1a(:,2) = interpn(DTI{1}.anatomy(:,:,:,2), Vertices_ras(:,1)/res(1),Vertices_ras(:,2)/res(2),Vertices_ras(:,3)/res(3));
V1a(:,3) = interpn(DTI{1}.anatomy(:,:,:,3), Vertices_ras(:,1)/res(1),Vertices_ras(:,2)/res(2),Vertices_ras(:,3)/res(3));
clear V2a;
V2a(:,1) = interpn(DTI{2}.anatomy(:,:,:,1), Vertices_ras(:,1)/res(1),Vertices_ras(:,2)/res(2),Vertices_ras(:,3)/res(3));
V2a(:,2) = interpn(DTI{2}.anatomy(:,:,:,2), Vertices_ras(:,1)/res(1),Vertices_ras(:,2)/res(2),Vertices_ras(:,3)/res(3));
V2a(:,3) = interpn(DTI{2}.anatomy(:,:,:,3), Vertices_ras(:,1)/res(1),Vertices_ras(:,2)/res(2),Vertices_ras(:,3)/res(3));
clear V3a;
V3a(:,1) = interpn(DTI{3}.anatomy(:,:,:,1), Vertices_ras(:,1)/res(1),Vertices_ras(:,2)/res(2),Vertices_ras(:,3)/res(3));
V3a(:,2) = interpn(DTI{3}.anatomy(:,:,:,2), Vertices_ras(:,1)/res(1),Vertices_ras(:,2)/res(2),Vertices_ras(:,3)/res(3));
V3a(:,3) = interpn(DTI{3}.anatomy(:,:,:,3), Vertices_ras(:,1)/res(1),Vertices_ras(:,2)/res(2),Vertices_ras(:,3)/res(3));
clear L1a L2a L3a;
L1a = interpn(DTI{4}.anatomy, Vertices_ras(:,1)/res(1),Vertices_ras(:,2)/res(2),Vertices_ras(:,3)/res(3));
L2a = interpn(DTI{5}.anatomy, Vertices_ras(:,1)/res(1),Vertices_ras(:,2)/res(2),Vertices_ras(:,3)/res(3));
L3a = interpn(DTI{6}.anatomy, Vertices_ras(:,1)/res(1),Vertices_ras(:,2)/res(2),Vertices_ras(:,3)/res(3));

La = inv(Transf(1:3,1:3)/Transf(4,4));

[V1rot,V2rot,V3rot] = PPD_linear(V1a,V2a,V3a,La);


%% check the results on the SCS coordinates
cfg3.node = femHead.Vertices;
cfg3.elem = [femHead.Elements(femHead.Tissue==1,:) femHead.Tissue(femHead.Tissue==1)];
cfg3.elem_centroide = bst_generate_centroide_on_elem(cfg3.node,cfg3.elem);


for i = 1 : length(V1rot)
eigen_vector{i} = [V1rot(i,:)' V2rot(i,:)' V3rot(i,:)'];
eigen_value{i} = [L1a(i) 0 0; 0 L2a(i) 0 ; 0 0 L3a(i)];
end

cfg3.eigen.eigen_vector = eigen_vector;
cfg3.eigen.eigen_value = eigen_value;


% comput centroide of the elemts 
cfg3.plotMesh = 1;
% cfg.elem = cfg.elem(el,:) ;
cfg3.indElem = 1: 2000:length(cfg3.elem);% 1:50:length(cfg.elem) ;%el(1:50:end) ;% 4500:6000;
cfg3.arrow = 1;
cfg3.ellipse = 0;


bst_display_tensor_as_ellipse(cfg3)

end