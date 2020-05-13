function [cfg] =  bst_generate_triedre_on_elem(cfg)

% This function will calculate the centeroid of each element of the mesh
% and the generate a triedre (or 3 basis vector). The first vector is
% is the normal pointing on the direction defined  from the center of the mesh to the centroid of the element.
% Then two other unit vectors are generated and perpenducular to the normal vector a basis on each element
% centroid.
% input : cfg.node and cfg.elem
% Takfarinas MEDANI, November 22, 2019

%Compute centroide:
elem_centroide = bst_generate_centroide_on_elem(cfg.node,cfg.elem);

% generate tiedre on each point
Startpoint = mean(cfg.node);
Endpoint = elem_centroide;
vector_centroid = Endpoint - Startpoint;
vector_norm_centroid = vector_centroid./sqrt(vector_centroid(:,1).^2+vector_centroid(:,2).^2+vector_centroid(:,3).^2); % vector with length 1

% compute the tangential vector
vector_norm_centroid_t1 = zeros(length(vector_norm_centroid),3);
vector_norm_centroid_t2 = zeros(length(vector_norm_centroid),3);

for ind = 1 : length(vector_norm_centroid)
    r1= vector_norm_centroid(ind,:);
    r = null(r1(:).');
    vector_norm_centroid_t1(ind,:) = r(:,1)';
    vector_norm_centroid_t2(ind,:) =r(:,2)';
end

cfg.elem_centroide = elem_centroide;
cfg.vector_norm_centroid = vector_norm_centroid;
cfg.vector_norm_centroid_t1 = vector_norm_centroid_t1;
cfg.vector_norm_centroid_t2 = vector_norm_centroid_t2;

%% Checking ... not activated
view = 0; check = 0;
if check ==1
    % check the length should be 1, unit vectors
    sqrt(vector_norm_centroid(:,1).^2+vector_norm_centroid(:,2).^2+vector_norm_centroid(:,3).^2);
    sqrt(vector_norm_centroid_t1(:,1).^2+vector_norm_centroid_t1(:,2).^2+vector_norm_centroid_t1(:,3).^2);
    sqrt(vector_norm_centroid_t2(:,1).^2+vector_norm_centroid_t2(:,2).^2+vector_norm_centroid_t2(:,3).^2);
    
    % check scallar products : should be zero
    A = vector_norm_centroid;
    B = vector_norm_centroid_t1;
    C = vector_norm_centroid_t2;
    dotP = dot(A,B,2);max(dotP)
    dotP = dot(A,C,2);max(dotP)
    dotP = dot(B,C,2);max(dotP)
end

if view ==1
    elem_centroide = cfg.elem_centroide;
    vector_norm_centroid = cfg.vector_norm_centroid;
    figure;plotmesh(cfg.node,cfg.elem,'x>0','facealpha',0.5);  hold on;plotmesh(elem_centroide,'x>0','r.')
    hold on; quiver3(elem_centroide(:,1),elem_centroide(:,2),elem_centroide(:,3),vector_norm_centroid(:,1),vector_norm_centroid(:,2),vector_norm_centroid(:,3))
    
    num_view = 1:1000;
    figure;plotmesh(cfg.node,cfg.elem(num_view,:),'facealpha',0.5);  hold on;plotmesh(elem_centroide(num_view,:),'rx')
    hold on; quiver3(elem_centroide(num_view,1),elem_centroide(num_view,2),elem_centroide(num_view,3),vector_norm_centroid(num_view,1),vector_norm_centroid(num_view,2),vector_norm_centroid(num_view,3))
    
    % check the view
    num_view = 1500:1800;
    figure;plotmesh(cfg.node,cfg.elem(num_view,:),'facealpha',0.5);  hold on;plotmesh(elem_centroide(num_view,:),'rx')
    hold on; quiver3(elem_centroide(num_view,1),elem_centroide(num_view,2),elem_centroide(num_view,3),vector_norm_centroid(num_view,1),vector_norm_centroid(num_view,2),vector_norm_centroid(num_view,3))
    hold on; quiver3(elem_centroide(num_view,1),elem_centroide(num_view,2),elem_centroide(num_view,3),vector_norm_centroid_t1(num_view,1),vector_norm_centroid_t1(num_view,2),vector_norm_centroid_t1(num_view,3))
    hold on; quiver3(elem_centroide(num_view,1),elem_centroide(num_view,2),elem_centroide(num_view,3),vector_norm_centroid_t2(num_view,1),vector_norm_centroid_t2(num_view,2),vector_norm_centroid_t2(num_view,3))
    legend('elem','centr','Vr','Vt1','Vt2');
end
end