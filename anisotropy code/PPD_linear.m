% SVReg: Surface-Constrained Volumetric Registration
% Copyright (C) 2019 The Regents of the University of California and the University of Southern California
% Created by Anand A. Joshi, Chitresh Bhushan, David W. Shattuck, Richard M. Leahy
%
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

function [W1 W2 W3] = PPD_linear(V1, V2, V3, L, ~)
%PPD Rotates the eigenvectors of the diffusion matrix with Preservation of
%Principle component(PPD) algorithm & generates new eigenvectors. Saves the
%output in .nii.gz file.
%
%   V1 - Principal eigenvector (corresponding to eigenvalue L1)
%   V2 - Second eigenvector (corresponding to eigenvalue L2)
%   V3 - Third eigenvector (corresponding to eigenvalue L3)
%       such that L1>L2>L3
%
%   F1, F2, F3 - Rows of gradient matrix (of the map)
%
%   W1 - Rotated principal eigenvector
%   W2 - Rotated second eigenvector
%   W3 - Rotated third eigenvector
%

%fprintf('\n==== Applying PPD ====\n');

[len bre dep dim4] = size(V1);

% Invert the gradient matrix
% detG =   G1(:,:,:,1) .* G3(:,:,:,3) .* G2(:,:,:,2) ...
%        - G1(:,:,:,1) .* G3(:,:,:,2) .* G2(:,:,:,3) ...
%        - G2(:,:,:,1) .* G3(:,:,:,3) .* G1(:,:,:,2) ...
%        + G2(:,:,:,1) .* G3(:,:,:,2) .* G1(:,:,:,3) ...
%        + G3(:,:,:,1) .* G2(:,:,:,3) .* G1(:,:,:,2) ...
%        + G3(:,:,:,1) .* G2(:,:,:,2) .* G1(:,:,:,3);
%
% F1(:,:,:,1) = (G3(:,:,:,3).*G2(:,:,:,2) - G3(:,:,:,2).*G2(:,:,:,3))./detG ;
% F1(:,:,:,2) = (G3(:,:,:,2).*G1(:,:,:,3) - G3(:,:,:,3).*G1(:,:,:,2))./detG ;
% F1(:,:,:,3) = (G2(:,:,:,3).*G1(:,:,:,2) - G2(:,:,:,2).*G1(:,:,:,3))./detG ;
%
% F2(:,:,:,1) = (G3(:,:,:,1).*G2(:,:,:,3) - G3(:,:,:,3).*G2(:,:,:,1))./detG ;
% F2(:,:,:,2) = (G3(:,:,:,3).*G1(:,:,:,1) - G3(:,:,:,1).*G1(:,:,:,3))./detG ;
% F2(:,:,:,3) = (G2(:,:,:,1).*G1(:,:,:,3) - G2(:,:,:,3).*G1(:,:,:,1))./detG ;
%
% F3(:,:,:,1) = (G3(:,:,:,2).*G2(:,:,:,1) - G3(:,:,:,1).*G2(:,:,:,2))./detG ;
% F3(:,:,:,2) = (G3(:,:,:,1).*G1(:,:,:,2) - G3(:,:,:,2).*G1(:,:,:,1))./detG ;
% F3(:,:,:,3) = (G2(:,:,:,2).*G1(:,:,:,1) - G2(:,:,:,1).*G1(:,:,:,2))./detG ;

%F1=L(1,:); F2=L(2,:); F3=L(3,:);
%clear L

% convert to double for higher precision
V1 = double(V1);
V2 = double(V2);
V3 = double(V3);


% Normalizing V1, V2, V3
V1 = norm_vector(V1);
V2 = norm_vector(V2);
V3 = norm_vector(V3);


% applying Jacobian to V1
%n1 = zeros([len bre dep dim4]);
n1 = (L*V1')';
n1 = norm_vector(n1);

% applying Jacobian to V2
n2 = (L*V2')';
n2 = norm_vector(n2);

cosTheta1 = sum(V1.*n1, 2);
axis1 = cross(V1', n1')';

% Projection of n2 perpendicular to n1
temp = zeros(size(n1));
n1_dot_n2 = sum(n1.*n2, 2);
temp(:,1) = n1(:,1).*n1_dot_n2;
temp(:,2) = n1(:,2).*n1_dot_n2;
temp(:,3) = n1(:,3).*n1_dot_n2;
proj_n2 = n2 - temp;
proj_n2 = norm_vector(proj_n2);
clear temp n2 n1_dot_n2

V2rot = rot_vector(V2, axis1, cosTheta1);

cosTheta2 = sum(V2rot.*proj_n2, 2);
axis2 = cross(V2rot, proj_n2);


W1 = V1;
W1 = n1;

W2 = V2;
W2 = rot_vector(V2rot, axis2, cosTheta2);

W3 = V3;
V3rot = rot_vector(V3, axis1, cosTheta1);
W3 = rot_vector(V3rot, axis2, cosTheta2);

end


function [ v_out ] = norm_vector( v )
%norm_vectorTOR Returns normalized vectors,
v_norm = sqrt(sum(v.^2, 2));

v_out = v ./ v_norm;

v_out(isnan(v_out)) = 0;

end



function [ v_rot ] = rot_vector( v, axis, cosTheta )
% Rotates vector v about axis by angle theta. Function takes cosine of
% theta as argument. axis may not be a unit vector. Here vectors are stored
% in 4th dimesion of matrix a & b (usually the eigenvectors)
%
%  Uses Rodrigues' rotation formula

axis = norm_vector(axis);

term1 = zeros(size(v));
term1(:,1) = v(:,1).*cosTheta;
term1(:,2) = v(:,2).*cosTheta;
term1(:,3) = v(:,3).*cosTheta;


term2 = zeros(size(v));
sinTheta = sqrt(1 - (cosTheta.^2));
temp = cross(axis, v);
term2(:,1) = temp(:,1).*sinTheta;
term2(:,2) = temp(:,2).*sinTheta;
term2(:,3) = temp(:,3).*sinTheta;
clear temp sinTheta


term3 = zeros(size(v));
temp = sum(axis.*v, 2).*(1-cosTheta);
term3(:,1) = axis(:,1).*temp;
term3(:,2) = axis(:,2).*temp;
term3(:,3) = axis(:,3).*temp;

v_rot = term1 + term2 + term3;

end

