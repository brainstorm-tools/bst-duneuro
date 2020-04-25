%% Display the DTI tensor on the mesh:

%% 1 - import the tissu from matlab or from the function 

%% Import the Masks / final results from the segmentation
tissues;
mri_segmented = bstMri2ftMri(tissues, tissues.InitTransf{2});
mri_segmented.anatomy(mri_segmented.anatomy==6) = 5;
mri_segmented.seg = mri_segmented.anatomy;
% keep only the White matter voxel to 1 and the rest to zero
wmLabel = 1;
wmAnatomy = mri_segmented;
wmAnatomy.anatomy(wmAnatomy.anatomy~=wmLabel)=0;
wmAnatomy.anatomy(wmAnatomy.anatomy==wmLabel)=1;
wmAnatomy.seg = wmAnatomy.anatomy;

% Generate the mesh 
cfg = [];
cfg.method = 'hexahedral';
cfg.spmversion = 'spm12';
cfg.downsample = 2;
cfg.shift = 0;
cfg.background = 0;
% ft part advanced version
mesh = ft_prepare_mesh(cfg, mri_segmented);
[tetraElem,tetraNode,tetraLabel] = hex2tet(double(mesh.hex), mesh.pos, double(mesh.tissue), 3);
figure; plotmesh(tetraNode,  [tetraElem tetraLabel],'x>0'); hold on;
xlabel('X'); ylabel('Y');  zlabel('Z');

mesh.pos_voxel = ft_warp_apply(inv(mri_segmented.transform), mesh.pos, 'homogeneous');
[tetraElem,tetraNode,tetraLabel] = hex2tet(double(mesh.hex), mesh.pos_voxel, double(mesh.tissue), 3);
figure; plotmesh(tetraNode,  [tetraElem tetraLabel],'x>100'); hold on;
xlabel('X'); ylabel('Y');  zlabel('Z');

%% Find the index of the WM 
cfg = [];
cfg.elem = [mesh.hex mesh.tissue];
cfg.node = mesh.pos_voxel;

cfg.elem_wm = cfg.elem((cfg.elem(:,end) == 1),:);
cfg.elem_centroide_wm = bst_generate_centroide_on_elem(cfg.node,cfg.elem_wm);

% display the WM
[tetraElem,tetraNode,tetraLabel] = hex2tet(cfg.elem_wm(:,1:end-1), cfg.node, cfg.elem_wm(:,end), 3);
figure; plotmesh(tetraNode,  [tetraElem tetraLabel],'x>100','facealpha',0.3); hold on;
plotmesh(cfg.elem_centroide_wm,'r.')

isWm = 1;
isNotWm = 1;
for i = 1 : size(cfg.elem_wm,1)       
    if cfg.elem_wm(i,end) == 1 % only white matter
        pos = round(cfg.elem_centroide_wm(i,:));        
        %     pos = cfg.node(cfg.elem(i,1),:);
        if sum(pos == 0) == 1; pos((pos == 0)) = 1;end
        if sum(pos <= 0); pos((pos <=  0)) = 1;end
        
        if (~(wmAnatomy.seg(pos(1),pos(2),pos(3))) == 0)
            disp([num2str(i) ' : is WM OK'])            
            V1{i} = squeeze(DTI{1}.anatomy(pos(1),pos(2),pos(3),:));
            V2{i} = squeeze(DTI{2}.anatomy(pos(1),pos(2),pos(3),:));
            V3{i} = squeeze(DTI{3}.anatomy(pos(1),pos(2),pos(3),:));
            L1{i} = squeeze(DTI{4}.anatomy(pos(1),pos(2),pos(3),:));
            L2{i} = squeeze(DTI{5}.anatomy(pos(1),pos(2),pos(3),:));
            L3{i} = squeeze(DTI{6}.anatomy(pos(1),pos(2),pos(3),:));

            isWm = isWm +1;
            positionWmTenso{isWm} = pos;
            isWmLocation(i) = 1;
        else
            disp(['index ' num2str(i) ' or elem ' num2str(cfg.elem_wm(i)) ' is not WM :KO '])
            V1{i} = squeeze(DTI{1}.anatomy(pos(1),pos(2),pos(3),:));
            V2{i} = squeeze(DTI{2}.anatomy(pos(1),pos(2),pos(3),:));
            V3{i} = squeeze(DTI{3}.anatomy(pos(1),pos(2),pos(3),:));
            L1{i} = squeeze(DTI{4}.anatomy(pos(1),pos(2),pos(3),:));
            L2{i} = squeeze(DTI{5}.anatomy(pos(1),pos(2),pos(3),:));
            L3{i} = squeeze(DTI{6}.anatomy(pos(1),pos(2),pos(3),:));
            isNotWm = isNotWm + 1;
            isWmLocation(i) = 0;
        end
    end
%     tensor{i} = [L1{i}*V1{i} L2{i}*V2{i}  L3{i}*V3{i}];
% conductivity_tensor3x3(1:3,1:3,ind) = [L1{i}*V1{i} L2{i}*V2{i}  L3{i}*V3{i}];
eigen_vector{i} = [V1{i} V2{i} V3{i}];
eigen_value{i} = [L1{i} 0 0; 0 L2{i} 0 ; 0 0 L3{i}];
end

i=75026
isWmLocation(i)
tensor{i} = [L1{i}*V1{i} L2{i}*V2{i}  L3{i}*V3{i}] ;

%% disply the tensors
% conductivity_tensor3x3: [3󫢬7822 double]

cfg.conductivity_tensor3x3 = conductivity_tensor3x3;

cfg.eigen.eigen_vector = eigen_vector;
cfg.eigen.eigen_value = eigen_value;

% comput centroide of the elemts 
cfg.plotMesh = 1;
% cfg.elem = cfg.elem(el,:) ;
cfg.indElem = 1: 100:length(cfg.elem_wm);% 1:50:length(cfg.elem) ;%el(1:50:end) ;% 4500:6000;
cfg.elem = cfg.elem_wm;
cfg.elem_centroide = cfg.elem_centroide_wm;

cfg.arrow = 0;
cfg.ellipse = 1;


bst_display_tensor_as_ellipse(cfg)


%% Extract the EigenValue and vector for each voxel
% apply the anand model
% eport head model from brainstorm a206 ds 2 noshift
[Vertices_ras, Transf] = cs_convert(tissues,  'scs',  'voxel', bstModel.Vertices);

cfg2.elem_wm = [bstModel.Elements(bstModel.Tissue == 1, :) bstModel.Tissue(bstModel.Tissue==1)];
cfg2.node = Vertices_ras;
elem_centroide_wm = bst_generate_centroide_on_elem(cfg2.node,cfg2.elem_wm);
cfg2.elem_centroide_wm = elem_centroide_wm;
% display the WM
[tetraElem,tetraNode,tetraLabel] = hex2tet(cfg2.elem_wm(:,1:end-1), cfg2.node, cfg2.elem_wm(:,end), 3);
figure; plotmesh(tetraNode,  [tetraElem tetraLabel],'x>100','facealpha',0.3); hold on;
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
cfg3.node = bstModel.Vertices;
cfg3.elem = [bstModel.Elements(bstModel.Tissue==1,:) bstModel.Tissue(bstModel.Tissue==1)];
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
cfg3.indElem = 1: 100:length(cfg3.elem);% 1:50:length(cfg.elem) ;%el(1:50:end) ;% 4500:6000;
cfg3.arrow = 0;
cfg3.ellipse = 1;


bst_display_tensor_as_ellipse(cfg3)



