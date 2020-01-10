function femhead = bst_msh2bst(mshfilename)
% convert the mat file to bst fem format 
% dependecies :  mesh_load_gmsh4
% Takfarinas MEDANI 9/19/2019
 
 m = mesh_load_gmsh4(mshfilename);
% covert to the bst format
femhead.Comment =mshfilename;
femhead.Vertices =   m.nodes(:,1:3);
femhead.Elements =    m.tetrahedra(:,1:4);
femhead.Tissue =     m.tetrahedron_regions;
%femhead.TissueLabels = {'1-WM'; '2-GM'; '3-CSF'; '4-Skull'; '5-Scalp'; '6-Eyes'};
tissuID = length(unique(femhead.Tissue));
for ind = 1 : tissuID
     femhead.TissueLabels{ind}  = [ 'tissu_' num2str(ind)];
end
femhead.History =     [];
femhead.Faces = m.triangles; 
femhead.Faces_tissue = m.triangle_regions; 
end