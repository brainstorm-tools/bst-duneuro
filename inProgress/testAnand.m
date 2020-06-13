%% test FEM
headFEM ;

surf.vertices = headFEM.Vertices;
surf.faces = headFEM.Elements;
surf.tissue = headFEM.Tissue;
% set a threshold point :
maxHeadPointZ = max(headFEM.Vertices(:,3));
minHeadPointZ = min(headFEM.Vertices(:,3));
delta = (maxHeadPointZ - minHeadPointZ);
ratio = 1/3; % par to remove from the mesh
threshold = minHeadPointZ+ratio*delta;
z0 = threshold;
% % find the outer skull point
% IndexC = strfind(headFEM.TissueLabels, 'csf');
% IndexSkull = find(not(cellfun('isempty',IndexC)));
% % extract the skull
% skullElem = headFEM.Elements(headFEM.Tissue ==IndexSkull,:);
% figure;
% plotmesh(headFEM.Vertices, skullElem,'facealpha',0.2);
% skullElem = unique(skullElem(:));
% minSkullPoint = min(headFEM.Vertices(skullElem,3));
% 
% figure;
% plotmesh(headFEM.Vertices, headFEM.Elements,'facealpha',0.2,'edgecolor','none');
% hold on
% plotmesh([mean(headFEM.Vertices(:,1)), ...
%                  mean(headFEM.Vertices(:,2)),...
%                  minSkullPoint],'ko','markersize',25);

% find all the faces that used these points:
tic,
m = 0;
if m == 0
    newFace = headFEM.Elements;
    newTissue = headFEM.Tissue;

    zCoordinate =  headFEM.Vertices(:,3);
    x = zCoordinate(newFace) <= z0;
    z = find(sum(x,2));
    toRemove = z;
end

if m==1
    belowPoint = find(surf.vertices(:,3)<=z0);
    length(belowPoint);
    length(surf.vertices);
    newFace = surf.faces;
    newTissue = surf.tissue;
    toRemove = [];
    for ind = 1 : length(belowPoint)
        temp = find(sum(surf.faces == belowPoint(ind),2));
        toRemove = [toRemove; temp];
    end
    toRemove = unique(toRemove);
end

newFace(toRemove,:) = [];
newTissue(toRemove,:)  = [];
%
% figure;
% plotmesh(surf.vertices, [surf.faces, surf.tissue],'x>0');
% hold on;
% plotmesh(surf.vertices, [newFace, newTissue],'x>0');

surfIN.vertices = headFEM.Vertices;
surfIN.faces = newFace;
surfIN.tissue = newTissue;
[NewSurf,UsedV]=delete_unused_vertices(surfIN);
timeOfCuting = toc;
figure;
plotmesh(surfIN.vertices, [surfIN.faces surfIN.tissue],'x>0');
