function bst_display_tensors(iSubject, eigenData)

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


%% Display the tensore within a mesh slice : defined by a plan
z0 = mean(femHead.Vertices(:,3));
z0 = 0.038;
plane=[min(femHead.Vertices(:,1)) min(femHead.Vertices(:,2)) (z0 + z0/2)
    min(femHead.Vertices(:,1)) max(femHead.Vertices(:,2)) (z0 + z0/2)
    max(femHead.Vertices(:,1)) min(femHead.Vertices(:,2)) (z0 + z0/2)];

y0 = 0.022;
plane=[min(femHead.Vertices(:,1)) y0 min(femHead.Vertices(:,3))
    min(femHead.Vertices(:,1)) y0 max(femHead.Vertices(:,3))
    max(femHead.Vertices(:,1)) y0 min(femHead.Vertices(:,3))];


x0 = 0.0;
plane=[x0  min(femHead.Vertices(:,2))  min(femHead.Vertices(:,3))
    x0  max(femHead.Vertices(:,2)) max(femHead.Vertices(:,3))
    x0  min(femHead.Vertices(:,2)) min(femHead.Vertices(:,3))];

% run qmeshcut to get the cross-section information at z=mean(node(:,1))
% use the x-coordinates as the nodal values

cfg.elem = cfg.elem(cfg.elem(:,5) == 1, :);

[cutpos,cutvalue,facedata,elemid]= ...
    qmeshcut(cfg.elem(:,1:4),cfg.node,zeros(length(cfg.node),1),plane);

figure;
plotmesh(femHead.Vertices,[femHead.Elements femHead.Tissue], 'facealpha', 0.3,'edgecolor','none');
hold on
plotmesh(femHead.Vertices,[femHead.Elements(elemid,:) femHead.Tissue(elemid,:)]);
view([0 90 0])

% cfg.elem = cfg.elem(cfg.elem(:,5) == 1, :);

cfg.indElem = elemid(1:1:end);
% cfg.indElem = find(cfg.elem(cfg.indElem ,end)==1);
cfg.eigen = eigenData;
cfg.noCsf  = 1; % do not display CSF
cfg.ellipse = 1;
cfg.arrow = 0;
cfg.plotMesh = 0;
bst_display_tensor_as_ellipse(cfg)
view([0 -90 0])
axis([ -0.063 0.095 ...
    min(cfg.node(:,2))  max(cfg.node(:,2)) ...
    -0.02  0.105])
hold on
plotmesh(femHead.Vertices,femHead.Elements(elemid,:),'facealpha',0.2,'edgecolor','none');
view([0 90 0])

hold on
plotmesh(femHead.Vertices,femHead.Elements(elemid,:),'facealpha',0.2);

view([0 0 90])


figure;
plotmesh(femHead.Vertices,femHead.Elements,'facealpha',0.2,'edgecolor','none');
view([0 0 90])


figure;
plotmesh(femHead.Vertices,femHead.Elements(elemid,:),'facealpha',0.2);
view([0 0 90])

end