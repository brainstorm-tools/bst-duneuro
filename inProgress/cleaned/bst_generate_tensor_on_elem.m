function cfg = bst_generate_tensor_on_elem(cfg)
% This function  tensor_elem = generate_tensor_on_elem(conductivity)
% will generate a tensor value on each element of a mesh.
% It maps the isotropy on the associated direction of the elem.
% The value are projected into the basis vectors generated on each element
% centroid.
% cfg.elem The input is the list of elem of the mesh,
% cfg.vectors : Th 3 vector basis defined on the cenroid of each elements.
% cfg.conductivity : conductivity value
% output :
% cfg.conductivity_tensor = tensor  [temp(1) temp(5) temp(9) temp(4) temp(7) temp(8)];
%                                                        x11       x22           x33       x12          x23        x33
% generate tensor from conductivity :

% Takfarinas MEDANI, december 2019


% conductivité isotropic
if ~isfield(cfg,'isotropic'); cfg.isotropic = 1; end

if cfg.isotropic == 1
    % Define tensor for isotropic model ()
    disp('Define tensor for isotropic model')
    conductivity=[cfg.conductivity;
        cfg.conductivity];
    if ~isfield(cfg,'conductivity_radial');     cfg.conductivity_radial = cfg.conductivity; end
    if ~isfield(cfg,'conductivity_tangential');     cfg.conductivity_tangential = cfg.conductivity;end
else
    disp('Define tensor for anisotropic model')
    conductivity = [cfg.conductivity_radial; cfg.conductivity_tangential];
end

if size(conductivity,2) ~= length(unique(cfg.elem(:,5)))
    error ('The size of the conductivity file is different to the number of layers')
end

% generate the local tensor according to the golobal coordinates X,Y,Z
A = zeros(3,3,size(conductivity,2)) ;
for ind = 1: size(conductivity,2)
    A(:,:,ind) = diag([conductivity(1,ind),conductivity(2,ind), conductivity(2,ind)]); % radial tangential tenagitial
end

%% Compute the basis vector on each elem
cfg =  bst_generate_triedre_on_elem(cfg);

%% Transformation matrix  and tensor mapping on each direction
tensor = zeros(length(cfg.elem),6) ;
conductivity_tensor3x3 = zeros(3,3,length(cfg.elem)) ;
for ind =1 : length(cfg.elem)
    T1 = [cfg.vector_norm_centroid(ind,:)' ...
        cfg.vector_norm_centroid_t1(ind,:)' ...
        cfg.vector_norm_centroid_t2(ind,:)'];
    
    temp = T1 * A(:,:,cfg.elem(ind,5)) * T1';
    tensor(ind,:) = [temp(1) temp(5) temp(9) temp(4) temp(7) temp(8)];
    conductivity_tensor3x3(:,:,ind) = (temp);   
    
    cfg.eigen.eigen_vector{ind} = T1;
    cfg.eigen.eigen_value{ind} = A(:,:,cfg.elem(ind,5));    
end
cfg.conductivity_tensor = tensor;
cfg.conductivity_tensor3x3 = conductivity_tensor3x3;
end
% Hi,
%
% it T should be T = [Vr(:) Vt1(:) Vt2(:)], so that T * T^t is and identity matrix. ^t indicates the transpose, yes.
%
% Unfortunately, I am not aware of any thesis that writes this down in detail, since it is usually only used for evaluation…
%
% In the end, you can also always visualize your tensors with, e.g., SCIRun to check if they make sense.
%
% Best,
% Johannes
%
% On 22. Oct 2019, at 19:57, Takfarinas Medani <takfarinesmedani@live.fr> wrote:
%
% Super Johannes,
%
% I follow your recommendation.
%
% To be sure about the computations, I have some basic questions/confirmations :
%
% >> Then you write those three vectors into a matrix to have a transformation matrix T
%
% if my vectors are :
%
% Vr = [Vrx Vry Vrz],
% Vt1 = [Vt1x Vt1y Vt1z],
% Vt2 = [Vt2x Vt2y Vt2z],
%
% the T should be :
%
% T = [Vr; Vt1; Vt2]      == > each line is vector
% or
% T = [Vr(:) Vt1(:) Vt2(:)]     == > each row is a vector
%
% >> By computing T * A * T^t you get your anisotropy tensor.
% What do you mean by  T^t  ?
% transpose of T?
%
% Just to have a more and clear idea, do you have any recommendations for a paper or thesis for this process?
%
%
% Thanks again and have a nice evening
%
%
% Amicalement
% Takfarinas
%
% From: Johannes <j.vorw01@gmail.com>
% Sent: Tuesday, October 22, 2019 06:35
% To: Takfarinas Medani <takfarinesmedani@live.fr>
% Subject: Re: Need help
%
% Hi Takfarinas,
%
% you're already halfway done ;-)
%
% As next step, you need to compute two tangential vector (i.e., normal to the radial vector you already have) so that you have a set of three orthonormal vectors. Then you write those three vectors into a matrix to have a transformation matrix T. Additionally, you set up a matrix with the (anisotropic) conductivities on the diagonal, i.e., something like A = diag(0.33,0.033, 0.033) if you the first vector you wrote into the matrix is the radial one. By computing T * A * T^t you get your anisotropy tensor.
%
% Best,
%     Johannes
%
% Am Di., 22. Okt. 2019 um 16:01 Uhr schrieb Takfarinas Medani <takfarinesmedani@live.fr>:
% Hi Johannes,
%
% Again I need your precious advice.
%
% In order to show the performance of the Duenuero/SimBio, I want to work in a simple scenario.
%
% I started from a spherical head model, I created a spherical mesh 4 sphere (wm,csf, skull and scalp),
% I compute the analytical than the FEM then comparison ... simple, in the case of the isotropy, there is no problem.
%
% In the case of the anisotropy, I'm facing a problem.
% I have an available analytical solution (Zhang and al) and the script is available on FT under the external/openmeeg folder.
%
% My question here is how to define the anisotropy for the elements (for WM inner sphere).
%
% My process is I computed the centroid of each element, then I compute the radial vector from the centroid.
% but after that I'm stuck, I don't really understand how I should define the conductivity format xx, yy, zz, xy, yz, zx ?
% for the FEM in this case.
%
% Any recommendations ...
%
% Thanks in advance for all your help.
%
% Kind regards
% Takfarinas
