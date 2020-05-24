
% centers = [1 2 3; 4 5 6; 7 8 9];
centers = femHeah.tensors.position;
% radii = [1 1 1; 1 2 1; 3 4 2];
factor = 100;
lengthElem = 1000;
L1xV1 = [];
figure

for ind = 1 : lengthElem
    disp([num2str(ind) ' / '  num2str(lengthElem)])
    l1 = femHeah.tensors.eigen_value{ind};
    v1= femHeah.tensors.eigen_vector{ind}(:,1);
    L1xV1(ind,:)=l1*v1;   
%     hold on;
%     quiver3(centers(ind,1), centers(ind,2), centers(ind,3), ...
%                 L1xV1(ind,1), L1xV1(ind,2), L1xV1(ind,3),'ShowArrowHead','off','color',abs([v1(2) v1(1) v1(3)]));
end

figure
 quiver3(centers(1:lengthElem,1), centers(1:lengthElem,2), centers(1:lengthElem,3), ...
                L1xV1(:,1), L1xV1(:,2), L1xV1(:,3),'ShowArrowHead','off');

axis('equal')

