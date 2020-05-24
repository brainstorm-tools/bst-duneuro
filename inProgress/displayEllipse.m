
% c0 = [0,0,0];
% r = 1;
% rr = [1,1,1];
% tsize = max(rr);
% 
%  [node,face,elem]=meshanellip(c0,rr,tsize);
% [node,face,elem]=meshasphere(c0,r,tsize);
% 
% figure; plotmesh(node,face)

[node,face] = mesh_sphere(30);
figure; plotmesh(node,face)

% centers = [1 2 3; 4 5 6; 7 8 9];
centers = femHeah.tensors.position;
% radii = [1 1 1; 1 2 1; 3 4 2];
factor = 100;
bemMerge = {};  temp = zeros(size(node));
for ind = 1 : 10000%length(centers)
    disp([num2str(ind) ' / '  num2str(ind(end))])
    radii = diag(femHeah.tensors.eigen_value{ind});
    vector = femHeah.tensors.eigen_vector{ind};
    % scaling
    temp(:,1) = (node(:,1)*radii(1))/factor;
    temp(:,2) = (node(:,2)*radii(2))/factor;
    temp(:,3) = (node(:,3)*radii(3))/factor;    
    % rotation
    for jnd = 1 : length(node)
        temp(jnd,:) = vector * temp(jnd,:)';
    end    
    % translation
    temp = temp+ centers(ind,:);
    
    bemMerge = cat(2, bemMerge, temp, face);
end
disp('Merging...')
[newnode, newelem] = mergemesh(bemMerge{:});
disp('Plotting...')
figure; h = plotmesh(newnode,newelem(:,1:3),'edgecolor','none','facecolor','k');


[newnode,newelem]=mergemesh(node,face,node1,face,node2,face);
plotmesh(newnode,newelem)



