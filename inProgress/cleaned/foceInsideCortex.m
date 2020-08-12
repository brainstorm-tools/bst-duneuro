clear all
FemMat = load('C:\Users\33649\Documents\MATLAB\brainstorm\brainstorm_db\brainstorm_db\bst-duneuro\anat\A0206\tess_fem200609_1411.mat');
cortex = load('C:\Users\33649\Documents\MATLAB\brainstorm\brainstorm_db\brainstorm_db\bst-duneuro\anat\A0206\tess_cortex_200521_1847.mat');

% extract GM
gmTissueID = 2;
% gm_tetra = FemMat.Elements(FemMat.Tissue <= gmTissueID,:);
[ngmt,egmt] = removeisolatednode(FemMat.Vertices,FemMat.Elements(FemMat.Tissue <= gmTissueID,:));
% figure; plotmesh(nn,ee(1:50,:),'facealpha',0.5)

gm_face = volface(FemMat.Elements(FemMat.Tissue <= gmTissueID,:));
[ngmf, egmf] = removeisolatednode(FemMat.Vertices,gm_face);
gMfv.vertices = ngmf;
gMfv.faces = egmf;
tic
gMout = ~inpolyhedron(gMfv, cortex.Vertices); sum(gMout)
gMindex_out = find(gMout);
twm = toc;
%% Compute the centroide of the gM elements
elem_centroide = zeros(size(egmt, 1), 3);
tic;
for l = 1:3
    elem_centroide(:, l) = sum(reshape(ngmt(egmt(:, :), l), size(egmt, 1), size(egmt, 2))')'/size(egmt, 2);
end
t2 = toc;
figure; hold on ; plotmesh(elem_centroide(1:50,:),'r.')

% For gM
k = dsearchn(elem_centroide,cortex.Vertices(gMindex_out,:));
NewVertices = cortex.Vertices;
NewVertices(gMindex_out ,:) = elem_centroide(k,:);
gMoutFinal = ~inpolyhedron(gMfv, NewVertices); sum(gMoutFinal)

figure; 
hold on; plotmesh(cortex.Vertices,cortex.Faces,'edgecolor','none')
hold on; plotmesh(cortex.Vertices,'k.')
hold on; plotmesh(NewVertices,'r*')

figure; 
hold on; plotmesh(cortex.Vertices(gMindex_out,:),'k.')
hold on; plotmesh(NewVertices(gMindex_out,:),'r.')


% extract WM
wmTissueID = 1;
wm_tetra = FemMat.Elements(FemMat.Tissue <= wmTissueID,:);

wm_face = volface(wm_tetra);
[wMno, wMel] = removeisolatednode(FemMat.Vertices,wm_face);
wMfv.vertices = wMno;
wMfv.faces = wMel;
% figure; plotmesh(wMfv.vertices,wMfv.faces ,'facealpha',0.5)

tic
wMin = inpolyhedron(wMfv, cortex.Vertices); sum(wMin)
wMindex_in = find(wMin);
twm = toc;
figure; plotmesh(wMfv.vertices,wMfv.faces ,'facealpha',0.5)
hold on; plotmesh(cortex.Vertices(wMindex_in,:),'r.')

% For WM
% NewVertices = cortex.Vertices;
GMcentroide = 0;
if GMcentroide == 1
    k = dsearchn(elem_centroide,cortex.Vertices(wMindex_in,:));
    NewVertices(wMindex_in ,:) = elem_centroide(k,:);
    
else
    k = dsearchn(gMfv.vertices,cortex.Vertices(wMindex_in,:));
    NewVertices(wMindex_in ,:) = gMfv.vertices(k,:);    
end
gMoutFinal = inpolyhedron(wMfv, NewVertices); sum(gMoutFinal)
% move again  to the centoide 
k = dsearchn(elem_centroide,NewVertices(wMindex_in,:));
NewVertices(wMindex_in ,:) = elem_centroide(k,:);
gMoutFinal = inpolyhedron(wMfv, NewVertices); sum(gMoutFinal)



% figure; plotmesh(NewVertices,'k.')
% hold on; plotmesh(cortex.Vertices,'ro')
% 
% 
% NewCortex = cortex;
% NewCortex.Vertices = NewVertices;
% NewCortex.Comment = [cortex.Comment ' corrected'];