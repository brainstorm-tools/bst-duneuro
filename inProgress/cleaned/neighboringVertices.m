% Generate points
N = 10;
NC = N^2;
[X,Y] = meshgrid(linspace(0,1,N),linspace(0,1,N));
x = reshape(X,N^2,1);  
y = reshape(Y,N^2,1);

% Generate triangulation
dtr = delaunayTriangulation(x,y);

x = dtr.Points(:,1);
y = dtr.Points(:,2);

% 1. Get all the triangles attached to a particular vertex in the
% triangulation.  
attachedTriangles = vertexAttachments(dtr);


for i = 1: size(x,1)
    % 2. Use the connectivity list to get the vertex indices of all these
    % triangles
    verticesOfTI         = dtr.ConnectivityList(attachedTriangles{i},:);
    % 3. Find all the unique vertices and remove the current vertex
    neighboursOfInternal{i} = setdiff(unique(verticesOfTI), i);
end


% Lets Check this
for j = 1: 5: size(x,1)
    triplot(dtr);
    hold on;
    plot(x(j),y(j),'rx')
    plot(x(neighboursOfInternal{j}),y(neighboursOfInternal{j}),'ro')
    hold off;
    pause(0.5)
end
close
