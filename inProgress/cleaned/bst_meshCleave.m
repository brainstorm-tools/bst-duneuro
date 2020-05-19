function [varargout] = bst_meshCleave(varargin)

% function [logicAt,logicAbove,logicBelow]=meshCleave(E,V,P,n,inclusiveSwitch)
% ------------------------------------------------------------------------
%
% 
% Kevin Mattheus Moerman
% 2020/05/18: Created
% Takfarinas MEDANI 
% 2020/05/2020 : adapted  to brainstorm head model 
%------------------------------------------------------------------------

%% Parse input
switch nargin
    case 2
        E=varargin{1};
        V=varargin{2};        
        P=mean(V,1);
        n=[0 0 1]; 
        inclusiveSwitch=[0 0];
    case 3
        E=varargin{1};
        V=varargin{2};        
        P=varargin{3};
        n=[0 0 1];        
        inclusiveSwitch=[0 0];
    case 4
        E=varargin{1};
        V=varargin{2};
        P=varargin{3};
        n=varargin{4};        
        inclusiveSwitch=[0 0];
    case 5
        E=varargin{1};
        V=varargin{2};
        P=varargin{3};
        n=varargin{4};
        inclusiveSwitch=varargin{5};
end

%Check P
if isempty(P)
    P=mean(V,1); %Use default as mean of coordinates
end

%Check normal direction
if isempty(n)
    n=[0 0 1]; %Default z-direction
end
n=vecnormalize(n); %Normalize n vector

%% Construct rotation matrix
nz=[0 0 1]; %The z-vector
dn=dot(n,nz);
rotAngle=acos(dot(n,nz)); %Rotation angle
%Determine rotation axis and formulate rotation matrix
if isapprox(dn,1,eps(1)) %n=nz
    R=eye(3,3);
else
    if isapprox(dn,-1,eps(1)) %n=-nz
        rotAxis=cross(n,[0 1 0]); %Use x-axis as rotation vector
    else % n~=nz and n~=-nz
        rotAxis=vecnormalize(cross(n,nz));
    end
    R=vecAngle2Rot(rotAngle,rotAxis);
end

%%

%Rotate coordinates
Vr=V-P(ones(size(V,1),1),:);
Vr=Vr*R';
Z=Vr(:,3);

%%

if inclusiveSwitch(1)==1
    logicVerticesBelow=Z<=0; %Logic for nodes below slice
else
    logicVerticesBelow=Z<0; %Logic for nodes below slice
end

if inclusiveSwitch(2)==1
    logicVerticesAbove=Z>=0; %Logic for nodes above slice
else
    logicVerticesAbove=Z>0; %Logic for nodes above slice
end

logicAt=any(logicVerticesBelow(E),2) & any(logicVerticesAbove(E),2); %Logic for slice elements
logicBelow=all(logicVerticesBelow(E),2); %Logic for elements with nodes all below
logicAbove=all(logicVerticesAbove(E),2); %Logic for elements with nodes all above

%%

varargout{1}=logicAt;
varargout{2}=logicAbove;
varargout{3}=logicBelow;

end


%% included function from Gibbon

function [V_norm]=vecnormalize(V)

if isvector(V)
    H=sqrt(sum(V.^2));
else
    H=sqrt(sum(V.^2,2));    
end
logicInvalid=H==0;

V_norm=V./H(:,ones(size(V,2),1));
V_norm(logicInvalid,:)=0; %Set invalid lengths to 0
end

function A=crossProdMat(a)

A=[ 0    -a(3)  a(2);...
    a(3)  0    -a(1);...
   -a(2)  a(1) 0];

end

function [R]=vecAngle2Rot(theta,w)

W=crossProdMat(w);
R=eye(3,3)+W*sin(theta)+W^2*(1-cos(theta));

end

function [L]=isapprox(varargin)
% function [L]=isapprox(A,B,tolLevel)
% ------------------------------------------------------------------------
% Check if A is approximately equal to B (to within tolLevel). 
%
% ------------------------------------------------------------------------
%%

switch nargin
    case 1
        A=varargin{1};
        B=zeros(size(A));
        tolLevel=eps(1);
    case 2
        A=varargin{1};
        B=varargin{2};
        tolLevel=eps(1);
    case 3
        A=varargin{1};
        B=varargin{2};
        tolLevel=varargin{3};
end

if numel(B)==1
    B=B.*ones(size(A));
end

L=abs(A-B)<tolLevel;
end