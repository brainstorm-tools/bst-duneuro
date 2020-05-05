function  [aniso_conductivity_tensor, eigenData, errMsg ] = bst_generate_conductivity_tensor(iSubject, DTI)

% tree_callbacks.m  10227  gui_component('MenuItem', jPopup, [], 'Generate FEM conductivity tensor', IconLoader.ICON_FEM, [], @(h,ev)generate_conductivity_tensor(iSubject, iAnatomy));

%% SECTION 1 : Get the data
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

% Get the Masks
bst_progress('text', 'Exporting MRI : Tissues/masks...');
MaskFile = file_fullpath(sSubject.Anatomy(iTissues).FileName);
MaskNii = bst_fullfile(brainsuiteDir, [subjid 'Mask.nii']);
sMriMask = in_mri_bst(MaskFile);


% Find the associated mesh from the data base
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

% Load the mesh
femHead=  load(FemFiles);
if numberOfEdges == 8
    %    quick check the mesh
    [tetraElem,tetraNode,tetraLabel] = hex2tet(double(femHead.Elements), femHead.Vertices, double(femHead.Tissue), 3);
    figure; plotmesh(tetraNode, [tetraElem tetraLabel],'x>0'); hold on;
    title('SCS/CTF coordinates systems')
elseif numberOfEdges == 4
    figure; plotmesh(femHead.Vertices,[femHead.Elements femHead.Tissue] ,'x>0');
    title('SCS/CTF coordinates systems')
else
    bst_error(['Tetrahedron element are not supported for now.' 10 'Check the Matlab command window for additional information.' 10], 'Generate FEM tensor', 0);
end

%% SECTION 2: Compute the isotropic tensor for the whole model
cfg = [];
default_iso_conductivity = get_standard_conductivity((numberOfLayer)) ;
conversion_m2mm = 1000;
default_iso_conductivity = default_iso_conductivity/conversion_m2mm;
cfg.elem = [femHead.Elements femHead.Tissue];
cfg.node = femHead.Vertices;
cfg.conductivity = default_iso_conductivity;
% Call the function that generate the tensors
cfg = bst_generate_tensor_on_elem(cfg);
% Save the isotropic tensor
IsotropicTensor = cfg.eigen;

check = 0;
if check = 1
    % Display the tensors
    cfg = [];
    cfg.IsotropicTensor = IsotropicTensor;
    cfg.eigen = IsotropicTensor;
    cfg.indElem = 1:50000:length(cfg.elem);
    cfg.ellipse = 1;
    cfg.arrow = 0;
    bst_display_tensor_as_ellipse(cfg)
    view([90 0 0])
    figure;
    plotmesh(cfg.node, cfg.elem)
end

%% SECTION 3 : Convert the tensors from the Voxel to SCS and interpolate to mesh
% DTI = load ('saved DTI from brainstorm');
[V1rot,V2rot,V3rot, L1a, L2a, L3a ] = bst_tensor_cs_convert(sMriMask, femHead, DTI);

%% SECTION 4 : Convert DWI to Conductivity
[aniso_conductivity_tensor, eigenData ] = bst_compute_anisotropy_tensors(femHead, default_iso_conductivity, IsotropicTensor, L1a,L2a,L3a,V1rot,V2rot,V3rot);



% check and display
check = 1;
if check
    % cfg.conductivity_tensor3x3 = aniso_conductivity/display_factor;
    cfg.indElem = find(cfg.elem(cfg.indElem ,end)==1);
    cfg.indElem = 1:500:length(cfg.elem);
    % cfg.indElem = cfg.indElem (1:10000:end);
    cfg.eigen = eigenData;
    cfg.noCsf  = 1; % do not display CSF
    cfg.ellipse = 0;
    cfg.arrow = 0;
    cfg.plotMesh = 1;
    bst_display_tensor_as_ellipse(cfg)
    view([0 90 0])
end
%
% %% Display the tensore within a mesh slice : defined by a plan
% z0 = mean(femHead.Vertices(:,3));
% z0 = 0.038;
% plane=[min(femHead.Vertices(:,1)) min(femHead.Vertices(:,2)) (z0 + z0/2)
%     min(femHead.Vertices(:,1)) max(femHead.Vertices(:,2)) (z0 + z0/2)
%     max(femHead.Vertices(:,1)) min(femHead.Vertices(:,2)) (z0 + z0/2)];
%
% y0 = 0.022;
% plane=[min(femHead.Vertices(:,1)) y0 min(femHead.Vertices(:,3))
%     min(femHead.Vertices(:,1)) y0 max(femHead.Vertices(:,3))
%     max(femHead.Vertices(:,1)) y0 min(femHead.Vertices(:,3))];
%
%
% x0 = 0.0;
% plane=[x0  min(femHead.Vertices(:,2))  min(femHead.Vertices(:,3))
%     x0  max(femHead.Vertices(:,2)) max(femHead.Vertices(:,3))
%     x0  min(femHead.Vertices(:,2)) min(femHead.Vertices(:,3))];
%
% % run qmeshcut to get the cross-section information at z=mean(node(:,1))
% % use the x-coordinates as the nodal values
%
% cfg.elem = cfg.elem(cfg.elem(:,5) == 1, :);
%
% [cutpos,cutvalue,facedata,elemid]= ...
%     qmeshcut(cfg.elem(:,1:4),cfg.node,zeros(length(cfg.node),1),plane);
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
%     min(cfg.node(:,2))  max(cfg.node(:,2)) ...
%     -0.02  0.105])
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
end