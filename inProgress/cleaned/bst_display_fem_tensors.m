function bst_display_fem_tensors(cfg)


if ~isfield(cfg,'plotMesh')
    cfg.plotMesh =1;
end

if ~isfield(cfg,'elemid')
    cfg.elemid = randperm(length(cfg.elem),50);
end

if ~isfield(cfg,'ellipse') && ~isfield(cfg,'arrow')
    cfg.ellipse = 1;
end

if  ~isfield(cfg,'conversion_m2mm') 
 cfg.conversion_m2mm = 1000;
 conversion_m2mm = cfg.conversion_m2mm;
end

tissueColor = {'y','k','b','g','r'};


% display
figure
for i = 1 : length(cfg.elemid)
    
    % get the center of the element
    xc = cfg.elem_centroide(cfg.elemid(i),1);
    yc = cfg.elem_centroide(cfg.elemid(i),2);
    zc = cfg.elem_centroide(cfg.elemid(i),3);
    
    % get the eigen data
    v = cfg.eigen.eigen_vector{cfg.elemid(i)} ; % vector 
    l = cfg.eigen.eigen_value{cfg.elemid(i)};   % value
    
    disp( [ num2str(i) ' / ' num2str(length(cfg.elemid)) ...
        ' | tissue : '  num2str(cfg.elem(cfg.elemid(i),5))...
        ' | volume : '  num2str(prod([l(1,1),l(2,2),l(3,3)]))]);
    
%     conversion_m2mm = 1; % 1000;
    factor = 4; % 4 is the optimal value for the SCS coordinates
    factor = factor/cfg.conversion_m2mm; % this is done because the conductivity is on S/meter, otherwise the size of the ellipse is bigger than the head
    
    % generate the grid
    if cfg.ellipse == 1
        hold on
        meshResolution = 10;
 
        [X,Y,Z] = ellipsoid(0,0,0,factor*norm(l(1,1)),factor*norm(l(2,2)),factor*norm(l(3,3)),meshResolution);
        % figure; surf(X,Y,Z); xlabel('X'); ylabel('Y'); zlabel('Z');
        % rotation
        sz=size(X);
        for x=1:sz(1)
            for y=1:sz(2)
                A=[X(x,y) Y(x,y) Z(x,y)]';
                A= v*A;
                X(x,y)=A(1);Y(x,y)=A(2);Z(x,y)=A(3);
            end
        end
        % translation
        X=X+xc; Y=Y+yc; Z=Z+zc;
        %         [i cfg.elem(i,end)]
        h= surf(X,Y,Z);
        %         	surf(sx*S(j)+X(j), sy*S(j)+Y(j), sz*S(j)+Z(j),...
        % 		'LineStyle','none',...
        % 		'FaceColor',C(j,:),...
        % 		'FaceAlpha',transp);
        %        shading flat
        if cfg.elem(cfg.elemid(i),5) == 1
            h.FaceColor = abs([v(2,1) v(1,1)  v(3,1)]);
        else
            h.FaceColor = tissueColor{cfg.elem(cfg.elemid(i),5)};
        end
        h.LineStyle = 'none';
    end
    
    if cfg.arrow == 1
        hold on
        quiver3(xc,yc,zc,factor*v(1,1)*l(1,1),factor*v(2,1)*l(1,1),factor*v(3,1)*l(1,1),1,'LineWidth', 1, 'Color',abs([v(2,1) v(1,1)  v(3,1)]));
        plot3(xc,yc,zc,'k.')
    end
    hold on;
end

disp('Plotting ...');

if cfg.ellipse == 1
    %     shading interp
    %     colormap([0.8 0.8 0.8])
    lighting phong
    light('Position',[0 0 1],'Style','infinite')
    %     light('Position',[0 0 1],'Style','infinite','Color',[ 1.000 0.584 0.000]);
    % from scatter
    %     daspect([ratios(1), ratios(2), ratios(3)]);
    % light('Position',[1 1 1],'Style','infinit','Color',[1 1 1]);
    % lighting gouraud;
    % view(30,30)
end
axis equal
xlabel('x');ylabel('y');zlabel('z');


if  cfg.plotMesh == 1
    shg; rotate3d on;
    if size(cfg.elem,2) == 9
        % convert to tetra
        [tetraElem,tetraNode,tetraLabel] = hex2tet(double(cfg.elem(:,1:end-1)), cfg.node, double(cfg.elem(:,end)), 3);
    else
        % hold on; plotmesh(cfg.node,cfg.elem(elemid,:),'facealpha',0.2); % hold on;plotmesh(cfg.elem_centroide(elemid,:),'k.')
        tetraElem = cfg.elem(:,1:end-1);
        tetraNode = cfg.node;
        tetraLabel = cfg.elem(:,end);
    end    
    %     hold on; plotmesh(tetraNode,[tetraElem, tetraLabel],'y>0','facealpha',0.1,'edgecolor','none','facecolor',[0.9 0.9 0.9]);
    disp('Press keyboard to dsiplay the head model')
    pause
    hold on; plotmesh(tetraNode,[tetraElem, tetraLabel],'y>0','facealpha',0.1,'edgecolor','none');
    disp('Press keyboard to dsiplay the mesh of the aniso tissue')
%     pause
%     hold on; plotmesh(tetraNode,[tetraElem(cfg.elemid,:) tetraElem(cfg.elemid,:)],'facealpha',0.4)
%     view([0 90 0])
    % hold on;plotmesh(cfg.elem_centroide(elemid,:),'k.')
    %     hold on; plotmesh(tetraNode,[tetraElem, tetraLabel],'x>50'); % hold on;plotmesh(cfg.elem_centroide(elemid,:),'k.')
end
end
