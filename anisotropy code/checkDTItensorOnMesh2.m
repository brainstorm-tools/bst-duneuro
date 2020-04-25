%% Display the DTI tensor on the mesh:

%% 1 - import the mesh 
% load the mesh within bst scs coordinates systemes
% Apply the transoformation from scs to vox

% Convert to Voxel coordinates systems
[femHead.Vertices_voxel, Transf] = cs_convert(sMriMask,  'scs',  'voxel', femHead.Vertices);
[femHead.Vertices_ras, Transf] = cs_convert(sMriMask,  'scs',  'world', femHead.Vertices);
[femHead.Vertices_mri, Transf] = cs_convert(sMriMask,  'scs',  'mri', femHead.Vertices);

% display
figure;
subplot(2,2,1)
plotmesh(femHead.Vertices, [femHead.Elements femHead.Tissue],'edgecolor','none'); title('SCS ');xlabel('x');ylabel('y');zlabel('z');grid on; grid minor;
subplot(2,2,2)
plotmesh(femHead.Vertices_voxel, [femHead.Elements femHead.Tissue],'edgecolor','none'); title('VOXEL');xlabel('x');ylabel('y');zlabel('z');grid on; grid minor;
subplot(2,2,3)
plotmesh(femHead.Vertices_ras, [femHead.Elements femHead.Tissue],'edgecolor','none'); title('RAS ');xlabel('x');ylabel('y');zlabel('z');grid on; grid minor;
subplot(2,2,4)
plotmesh(femHead.Vertices_mri, [femHead.Elements femHead.Tissue],'edgecolor','none'); title('MRI ');xlabel('x');ylabel('y');zlabel('z');grid on; grid minor;


%% 2 - compute the centroid of each element of the mesh 
cfg = [];
cfg.node = double(femHead.Vertices_voxel);
cfg.elem = double([femHead.Elements femHead.Tissue]);
elem_centroide = bst_generate_centroide_on_elem(cfg.node,cfg.elem);

figure;
if size(cfg.elem,2) == 9
    % convert to tetra
    [tetraElem,tetraNode,tetraLabel] = hex2tet(double(cfg.elem(:,1:end-1)), cfg.node, double(cfg.elem(:,end)), 3);
else
    % hold on; plotmesh(cfg.node,cfg.elem(indElem,:),'facealpha',0.2); % hold on;plotmesh(cfg.elem_centroide(indElem,:),'k.')
    tetraElem = cfg.elem(:,1:end-1);
    tetraNode = cfg.node;
    tetraLabel = cfg.elem(:,end);
end
plotmesh(tetraNode, [tetraElem tetraLabel],'facealpha',0.4); title('VOXEL');xlabel('x');ylabel('y');zlabel('z');grid on; grid minor;
hold on;
plotmesh(cfg.elem_centroide,'k.');


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
figure; plotmesh(tetraNode,  [tetraElem tetraLabel],'x>0'); hold on;
xlabel('X'); ylabel('Y');  zlabel('Z');

% light version 
% tmpcfg.downsample = cfg.downsample;
% tmpcfg.showcallinfo = 'yes';
% mri = ft_volumedownsample(tmpcfg, mri_segmented);
% seg_build.seg = mri.seg;
% seg_build.dim = mri.dim;
% % build the mesh
% mesh = build_mesh_hexahedral(cfg, seg_build);
% % convert the mesh to tetr in order to use plotmesh
% [tetraElem,tetraNode,tetraLabel] = hex2tet(double(mesh.hex), mesh.pos, double(mesh.labels), 3);
% figure; plotmesh(tetraNode,  [tetraElem tetraLabel],'x>0'); hold on;
% xlabel('X'); ylabel('Y');  zlabel('Z');

% converting position of meshpoints to the head coordinate system
% mesh.pos = ft_warp_apply(mri.transform, mesh.pos, 'homogeneous');

%% Find the index of the WM 
wmIndex = find(wmAnatomy.seg);

% find the position associated with the mesh 
cfg = [];
cfg.elem = [mesh.hex mesh.tissue];
cfg.node = mesh.pos_voxel;
cfg.elem_centroide = bst_generate_centroide_on_elem(cfg.node,cfg.elem);

cfg.elem_wm = cfg.elem((cfg.elem(:,end) == 1),:);

for i = 1 : size(cfg.elem_wm,1)   
    
    if cfg.elem_wm(i,end) == 1 % only white matter
        pos = round(cfg.elem_centroide(cfg.elem_wm(i),:));        
        %     pos = cfg.node(cfg.elem(i,1),:);
        pos = round(pos);
        if sum(pos == 0) == 1; pos((pos == 0)) = 1;end
        if sum(pos <= 0); pos((pos <=  0)) = 1;end
        
        if (~(wmAnatomy.seg(pos(1),pos(2),pos(3))) == 0)
            disp([num2str(i) ' : is WM OK'])            
            V1{i} = squeeze(DTI{1}.anatomy(pos(1),pos(2),pos(3),:));
        else
            disp(['index ' num2str(i) ' or elem ' num2str(cfg.elem_wm(i)) ' is not WM :KO '])
            V1{i} = squeeze(DTI{1}.anatomy(pos(1),pos(2),pos(3),:));
        end
    end    
end

%% Extract the EigenValue and vector for each voxel
% 













