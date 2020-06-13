function [surf,UsedV]=delete_unused_vertices(surf)
%Author Anand A. Joshi ajoshi@sipi.usc.edu
%UsedV=reshape(surf.faces,size(surf.faces,1)*3,1);
UsedV=surf.faces(:);
UsedV=unique(UsedV);

%usedmarker=zeros(size(surf.vertices,1),1); usedmarker(UsedV)=1;
new_indx=zeros(size(surf.vertices,1),1);
%num_used=sum(usedmarker);
num_used=size(UsedV,1);
new_indx(UsedV)=[1:num_used];

surf.vertices=surf.vertices(UsedV,:);
surf.faces=new_indx(surf.faces);
% surf.tissue=new_indx(surf.tissue);
end