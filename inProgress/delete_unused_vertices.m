% SVReg: Surface-Constrained Volumetric Registration
% Copyright (C) 2019 The Regents of the University of California and the University of Southern California
% Created by Anand A. Joshi, Chitresh Bhushan, David W. Shattuck, Richard M. Leahy 
% Update and adapted to Tetrahedral mesh by Takfarinas MEDANI

% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; version 2.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301,
% USA.



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
surf.tissue=new_indx(surf.tissue);
