%% Display the DTI tensor on the mesh:

%% 1 - import the mesh 
% load the mesh within bst scs coordinates systemes
% Apply the transoformation from scs to vox

% Convert to Voxel coordinates systems
[femHead.Vertices_voxel, Transf] = cs_convert(sMriMask,  'scs',  'voxel', femHead.Vertices);
[femHead.Vertices_ras, Transf] = cs_convert(sMriMask,  'scs',  'world', femHead.Vertices);
[femHead.Vertices_mri, Transf] = cs_convert(sMriMask,  'scs',  'mri', femHead.Vertices);

%  view_mri
%  view_mri(MaskFile);
%  script_view_mri_overlay(MaskFile,MaskFile)
% %  view_mri_3d(MriFile, OverlayFile, SurfAlpha, hFig)
%   view_mri_3d(MaskFile)

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
cfg.elem_centroide = bst_generate_centroide_on_elem(cfg.node,cfg.elem);

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
[femMesh.Vertices_voxel, Transf] = cs_convert(tissues,  'scs',  'voxel', femMesh.Vertices);

cs_convert(tissues, 'world', 'voxel')

mri_segmented = bstMri2ftMri(tissues, cs_convert(tissues, 'world', 'voxel'));
mri_segmented.anatomy(mri_segmented.anatomy==6) = 5;
mri_segmented.seg = mri_segmented.anatomy;

mri_segmented.transform =cs_convert(sMri, 'world', 'mri');


cfg = [];
cfg.method = 'hexahedral';
cfg.spmversion = 'spm12';
cfg.downsample = 5;
cfg.shift = 0;
mesh = ft_prepare_mesh(cfg, mri_segmented);
[tetraElem,tetraNode,tetraLabel] = hex2tet(double(mesh.hex), mesh.pos, double(mesh.tissue), 3);

hFig = figure;
hAx  = axes('Parent',hFig);
fixedVolume = mri_segmented.anatomy;
centerFixed = size(fixedVolume)/2;
slice(hAx,double(fixedVolume),centerFixed(1),centerFixed(2),centerFixed(3));
shading(hAx,'interp');
set(hAx,'Xgrid','on','YGrid','on','ZGrid','on');
set(hFig,'Colormap',colormap('gray'));
hold on;
plotmesh(tetraNode,  [tetraElem tetraLabel],'x>0'); hold on;
xlabel('X'); ylabel('Y');  zlabel('Z');


mri_segmented = tissues;
% keep only the White matter voxel to 1 and the rest to zero
wmLabel = 1;
wmAnatomy = mri_segmented;
wmAnatomy.Cube(wmAnatomy.Cube~=wmLabel)=0;
wmAnatomy.Cube(wmAnatomy.Cube==wmLabel)=1;
% 
hFig = figure;
hAx         = axes('Parent',hFig);
fixedVolume = mri_segmented.Cube;
centerFixed = size(fixedVolume)/2;
slice(hAx,double(fixedVolume),centerFixed(1),centerFixed(2),centerFixed(3));
shading(hAx,'interp');
set(hAx,'Xgrid','on','YGrid','on','ZGrid','on');
set(hFig,'Colormap',colormap('gray'));
hold on;
plotmesh(femMesh.Vertices_voxel,  [femMesh.Elements femMesh.Tissue],'x>0'); hold on;
xlabel('X'); ylabel('Y');  zlabel('Z');

figure;
plotmesh(femMesh.Vertices,  [femMesh.Elements femMesh.Tissue],'x>0'); hold on;
xlabel('X'); ylabel('Y');  zlabel('Z');


for i = 1 : size(cfg.elem,1)
    pos = cfg.node(cfg.elem(i,1),:);
    pos = round(pos); 
    if sum(pos == 0) == 1; pos((pos == 0)) = 1;end
    if sum(pos <= 0); pos((pos <=  0)) = 1;end
 
    indexPos{i} = pos;
    
    if (~(mask.anatomy(pos(1),pos(2),pos(3)) == 0))
    end
end

%% Extract the EigenValue and vector for each voxel














