function bst_VolMesh2SurfMesh(TessFiles)
% This function will be used to extract the outer surfaces mesh
% from the fem volume mesh file of the bst (brainstrom format)
global GlobalData % ToDO : is there a better way to reach the MRI ? ==> without calling GlobalData
data = GlobalData.DataBase.ProtocolInfo(GlobalData.DataBase.iProtocol);
% out_mri_nii(bst_fullfile(data.SUBJECTS, MriFile), NiiFile);

% Extract surface from brainstorm FEM mesh
% femheadfile = 'tess_fem_simnibs_717852V.mat';
femheadfile = fullfile(data.SUBJECTS ,TessFiles);
femhead  = load(fullfile(data.SUBJECTS ,TessFiles));

% Open progress bar
bst_progress('start', 'FEM mesh', 'Surface Extraction From Mesh mesh (iso2mesh)...');

% load the data
% femhead = load(femheadfile);
elemID = unique(femhead.Tissue);
for ind = 1 : length(elemID)
    % assuming that the tissu are ordered from inner to outer
    index(1:ind) =  elemID(1:ind);
    combined = femhead.Tissue == index;
    combined =   sum(combined,2);
    
    % extract the surface enclosed the combined tissu
    [openface,~]=volface(femhead.Elements(find(combined), :));
    figure; plotmesh(femhead.Vertices,openface,'z>0'), title(femhead.TissueLabels{ind})
    [Vertices,Faces]=removeisolatednode(femhead.Vertices,openface);
    VertConn = []; VertNormals = []; Curvature = []; SulciMap = []; Atlas = []; iAtlas =[]; tess2mri_interp = []; Reg = [];
    History = ['surface initially generated from file "' femheadfile '"'];
    Comment = femhead.TissueLabels{ind};
    % save to brainstrom surface file
    eval([(([genvarname(['tess_' femhead.TissueLabels{ind}]), '.Comment' ])) '= femhead.TissueLabels{ind}']);
    eval([(([genvarname(['tess_' femhead.TissueLabels{ind}]), '.Vertices' ])) '= Vertices']);
    eval([(([genvarname(['tess_' femhead.TissueLabels{ind}]), '.Faces' ])) '= Faces']);
    eval([(([genvarname(['tess_' femhead.TissueLabels{ind}]), '.VertConn' ])) '= VertConn']);
    eval([(([genvarname(['tess_' femhead.TissueLabels{ind}]), '.VertNormals' ])) '= VertNormals']);
    eval([(([genvarname(['tess_' femhead.TissueLabels{ind}]), '.Curvature' ])) '= Curvature']);
    eval([(([genvarname(['tess_' femhead.TissueLabels{ind}]), '.SulciMap' ])) '= SulciMap']);
    eval([(([genvarname(['tess_' femhead.TissueLabels{ind}]), '.Atlas' ])) '= Atlas']);
    eval([(([genvarname(['tess_' femhead.TissueLabels{ind}]), '.iAtlas' ])) '= iAtlas']); %%%% Check if this should be imlemented
    eval([(([genvarname(['tess_' femhead.TissueLabels{ind}]), '.tess2mri_interp' ])) '= tess2mri_interp']);
    eval([(([genvarname(['tess_' femhead.TissueLabels{ind}]), '.Reg' ])) '= Reg']);
    eval([(([genvarname(['tess_' femhead.TissueLabels{ind}]), '.History' ])) '= History']);
    
    % save to the data base
    save(fullfile(data.SUBJECTS, genvarname(['tess_' femhead.TissueLabels{ind}])), ...
        'Comment' , ...
        'Vertices', ...
        'Faces' , ...
        'VertConn' , ...
        'VertNormals' , ...
        'Curvature' , ...
        'SulciMap' , ...
        'Atlas' , ...
        'iAtlas', ...
        'tess2mri_interp' , ...
        'Reg' , ...
        'History' )
    
    % % Load the bst data base
% ===== SAVE IN DATABASE =====
isSave = 1;
fileTag = 'VolMesh2SurfMesh';
if isSave
    % Create new filename
    NewTessFile = bst_fullfile( ['tess_' femhead.TissueLabels{ind} '.mat']);
    FemFile = file_unique(NewTessFile);
    % Save file
%     bst_save(FemFile, ['tess_' femhead.TissueLabels{ind}], 'v7');    
    % Get subject
    [sSubject, iSubject] = bst_get('SurfaceFile', TessFiles);   
    % Make output filename relative
    NewTessFile = file_short(NewTessFile);    
    % Register this file in Brainstorm database
    NewComment = [Comment '_VolMesh2SurfMesh'];         fileType = 'other';
    iSurface = db_add_surface(iSubject, FemFile, NewComment);
else
    NewTessFile = NewTess;
    iSurface = [];
end


end


% Close progress bar
bst_progress('stop');

end

