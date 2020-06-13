function bst_cut_femMesh(iSubject)
% BST_CUT_FEMMESH: Cut the FEM mesh and remove the bottom part.

% @=============================================================================
% This function is part of the Brainstorm software:
% https://neuroimage.usc.edu/brainstorm
% 
% Copyright (c)2000-2020 University of Southern California & McGill University
% This software is distributed under the terms of the GNU General Public License
% as published by the Free Software Foundation. Further details on the GPLv3
% license can be found at http://www.gnu.org/copyleft/gpl.html.
% 
% FOR RESEARCH PURPOSES ONLY. THE SOFTWARE IS PROVIDED "AS IS," AND THE
% UNIVERSITY OF SOUTHERN CALIFORNIA AND ITS COLLABORATORS DO NOT MAKE ANY
% WARRANTY, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO WARRANTIES OF
% MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE, NOR DO THEY ASSUME ANY
% LIABILITY OR RESPONSIBILITY FOR THE USE OF THIS SOFTWARE.
%
% For more information type "brainstorm license" at command prompt.
% =============================================================================@
%
% Authors: Takfarinas Medani, Anand Joshi, 2020

% This function is 
%% SECTION 1 : Get the data
disp('cut fem mesh')
bst_progress('start', 'Cutting Mesh','Load mesh...')

% Get the input data
ProtocolInfo     = bst_get('ProtocolInfo');
ProtocolSubjects = bst_get('ProtocolSubjects');
% Default subject
if (iSubject == 0)
    sSubject = ProtocolSubjects.DefaultSubject;
else
    sSubject = ProtocolSubjects.Subject(iSubject);
end

% Get the mesh file
FemFiles = file_fullpath(sSubject.Surface(sSubject.iFEM).FileName);
[filepath,name,ext]  = fileparts(FemFiles);
% Load the mesh
headFEM=  load(FemFiles);

%% SECTION 2 : Cut the mesh

% ask the user for the threshold ratio to removefrtom the bottom part
ratio = 0.33; % par to remove from the mesh
[res, isCancel] = java_dialog(  'input' , 'Ratio of the bottom part to remove  =  ' ,...
    'Cut the FEM mesh' ,[],num2str(ratio) );
if isCancel;         return;    end
ratio = str2double( res);

% set a threshold point :
maxHeadPointZ = max(headFEM.Vertices(:,3));
minHeadPointZ = min(headFEM.Vertices(:,3));
delta = (maxHeadPointZ - minHeadPointZ);
threshold = minHeadPointZ+ratio*delta;
z0 = threshold;

bst_progress('start', 'Cutting Mesh','Cutting mesh...')

tic,
% this is the fastest method
newFace = headFEM.Elements;
newTissue = headFEM.Tissue;
zCoordinate =  headFEM.Vertices(:,3);
temp = zCoordinate(newFace) <= z0;
toRemove = find(sum(temp,2));
newFace(toRemove,:) = [];
newTissue(toRemove,:)  = [];
surfIN.vertices = headFEM.Vertices;
surfIN.faces = newFace;
surfIN.tissue = newTissue;
[NewSurf,UsedV]=delete_unused_vertices(surfIN);
timeOfCuting = toc;

bst_progress('start', 'Cutting Mesh','save mesh to database ...')
% Create output structure
FemMat = db_template('femmat');
FemMat = headFEM;
FemMat.Vertices = NewSurf.vertices;
FemMat.Elements =NewSurf.faces;
FemMat.Tissue =NewSurf.tissue;
FemMat.Comment = sprintf('FEM cut %dV (ratio : %d%s, %d layers)', length(FemMat.Vertices), round(ratio,2)*100,'%', length(FemMat.TissueLabels));
% Add history
FemMat = bst_history('add', FemMat, 'process_cut_fem');
% Save to database
FemFile = file_unique(bst_fullfile((filepath), sprintf('tess_fem_%s_%dV.mat', 'cut_mesh', length(FemMat.Vertices))));
bst_save(FemFile, FemMat, 'v7');
db_add_surface(iSubject, FemFile, FemMat.Comment);
bst_progress('stop')
end



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
end