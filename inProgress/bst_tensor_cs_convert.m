function  [V1rot,V2rot,V3rot, L1a, L2a, L3a ] = bst_tensor_cs_convert(sMriMask, femHead, DTI)

% This function convert the tensors from the voxel space to the scs
% coordinate system and interpolate theire value to the centroide of the
% mesh element

% Takfarinas & Anand

% input 
wmIndex = 1;

% much faster to use only the anisotrop tissues
% apply the anand model
[Vertices_vox, Transf] = cs_convert(sMriMask,  'scs',  'voxel', femHead.Vertices);

% extract only the wm ... faster and only wm is needed
elem_wm = [femHead.Elements(femHead.Tissue == 1, :) femHead.Tissue(femHead.Tissue==wmIndex)];
elem_centroide_wm = bst_generate_centroide_on_elem(Vertices_vox,elem_wm);

%% Interplate the tensors to the centroides of the anisotropic tissue 
dti_resolution = DTI{1}.hdr.dim.pixdim(2:4);
% Interpolate the results on the centroides of the WM elements
Vertices_query = elem_centroide_wm;

% Interpolate the eigen vectors
clear V1a;
V1a(:,1) = interpn(DTI{1}.anatomy(:,:,:,1), Vertices_query(:,1)/dti_resolution(1),Vertices_query(:,2)/dti_resolution(2),Vertices_query(:,3)/dti_resolution(3));
V1a(:,2) = interpn(DTI{1}.anatomy(:,:,:,2), Vertices_query(:,1)/dti_resolution(1),Vertices_query(:,2)/dti_resolution(2),Vertices_query(:,3)/dti_resolution(3));
V1a(:,3) = interpn(DTI{1}.anatomy(:,:,:,3), Vertices_query(:,1)/dti_resolution(1),Vertices_query(:,2)/dti_resolution(2),Vertices_query(:,3)/dti_resolution(3));
clear V2a;
V2a(:,1) = interpn(DTI{2}.anatomy(:,:,:,1), Vertices_query(:,1)/dti_resolution(1),Vertices_query(:,2)/dti_resolution(2),Vertices_query(:,3)/dti_resolution(3));
V2a(:,2) = interpn(DTI{2}.anatomy(:,:,:,2), Vertices_query(:,1)/dti_resolution(1),Vertices_query(:,2)/dti_resolution(2),Vertices_query(:,3)/dti_resolution(3));
V2a(:,3) = interpn(DTI{2}.anatomy(:,:,:,3), Vertices_query(:,1)/dti_resolution(1),Vertices_query(:,2)/dti_resolution(2),Vertices_query(:,3)/dti_resolution(3));
clear V3a;
V3a(:,1) = interpn(DTI{3}.anatomy(:,:,:,1), Vertices_query(:,1)/dti_resolution(1),Vertices_query(:,2)/dti_resolution(2),Vertices_query(:,3)/dti_resolution(3));
V3a(:,2) = interpn(DTI{3}.anatomy(:,:,:,2), Vertices_query(:,1)/dti_resolution(1),Vertices_query(:,2)/dti_resolution(2),Vertices_query(:,3)/dti_resolution(3));
V3a(:,3) = interpn(DTI{3}.anatomy(:,:,:,3), Vertices_query(:,1)/dti_resolution(1),Vertices_query(:,2)/dti_resolution(2),Vertices_query(:,3)/dti_resolution(3));
clear L1a L2a L3a;
L1a = interpn(DTI{4}.anatomy, Vertices_query(:,1)/dti_resolution(1),Vertices_query(:,2)/dti_resolution(2),Vertices_query(:,3)/dti_resolution(3));
L2a = interpn(DTI{5}.anatomy, Vertices_query(:,1)/dti_resolution(1),Vertices_query(:,2)/dti_resolution(2),Vertices_query(:,3)/dti_resolution(3));
L3a = interpn(DTI{6}.anatomy, Vertices_query(:,1)/dti_resolution(1),Vertices_query(:,2)/dti_resolution(2),Vertices_query(:,3)/dti_resolution(3));

%% Anand's code
La = inv(Transf(1:3,1:3)/Transf(4,4));
% Rotate the eigen vectors
[V1rot,V2rot,V3rot] = PPD_linear(V1a,V2a,V3a,La);

%% check the results on the SCS coordinates
check = 0;
if check
    cfg3 = [];
    cfg3.node = femHead.Vertices;
    cfg3.elem = [femHead.Elements(femHead.Tissue==1,:) femHead.Tissue(femHead.Tissue==1)];
    cfg3.elem_centroide = bst_generate_centroide_on_elem(cfg3.node,cfg3.elem);
    figure; plotmesh(cfg3.node, cfg3.elem(1,:),'facealpha',0.3); hold on;plotmesh(cfg3.elem_centroide(1,:),'r.')
    
    clear eigen_vector eigen_value
    for i = 1 : length(V1rot)
        eigen_vector{i} = [V1rot(i,:)' V2rot(i,:)' V3rot(i,:)'];
        eigen_value{i} = [L1a(i) 0 0; 0 L2a(i) 0 ; 0 0 L3a(i)];
    end    
    cfg3.eigen.eigen_vector = eigen_vector;
    cfg3.eigen.eigen_value = eigen_value;    
    % comput centroide of the elemts
    cfg3.plotMesh = 1;
    % cfg.elem = cfg.elem(el,:) ;
    cfg3.indElem = 1: 1000:length(cfg3.elem);% 1:50:length(cfg.elem) ;%el(1:50:end) ;% 4500:6000;
    cfg3.arrow = 0;
    cfg3.ellipse = 1;
    bst_display_tensor_as_ellipse(cfg3)
end
end