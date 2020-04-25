function  [V1rot,V2rot,V3rot, L1a, L2a, L3a ]=convert3Dtensors(sMriMask, femHead, DTI)
% Takfarinas & Anand

%% Extract the EigenValue and vector for each voxel
% apply the anand model
[Vertices_ras, Transf] = cs_convert(sMriMask,  'scs',  'voxel', femHead.Vertices);
% extract only the wm ... faster and only wm is needed
cfg2.elem_wm = [femHead.Elements(femHead.Tissue == 1, :) femHead.Tissue(femHead.Tissue==1)];
cfg2.node = Vertices_ras;
elem_centroide_wm = bst_generate_centroide_on_elem(cfg2.node,cfg2.elem_wm);
cfg2.elem_centroide_wm = elem_centroide_wm;
% if femHead.Elements >5
% % display the WM
% [tetraElem,tetraNode,tetraLabel] = hex2tet(cfg2.elem_wm(:,1:end-1), cfg2.node, cfg2.elem_wm(:,end), 3);
% figure; plotmesh(tetraNode,  [tetraElem tetraLabel],'x>100','facealpha',0.3); hold on;
% else
%     figure; plotmesh(cfg2.node, cfg2.elem_wm,'x>100','facealpha',0.3); hold on;
%     figure; plotmesh(cfg2.node, cfg2.elem_wm(1,:),'facealpha',0.3); hold on;
% end
% hold on
% plotmesh(cfg2.elem_centroide_wm(1,:),'r.')

% SZ=size(DTI{1}.anatomy);
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
check = 0;
if check
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
    cfg3.indElem = 1: 500:length(cfg3.elem);% 1:50:length(cfg.elem) ;%el(1:50:end) ;% 4500:6000;
    cfg3.arrow = 0;
    cfg3.ellipse = 1;
    bst_display_tensor_as_ellipse(cfg3)
end
end