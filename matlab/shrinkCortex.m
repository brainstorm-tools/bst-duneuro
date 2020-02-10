% shrink the cortex frm
    cortex_surface.Vertices = Vertices_gm;
    cortex_surface.Faces = Faces_gm;
    opts.decimate_cortex = 0;
    opts.source_depth = 1.5; % set the depth to 1.5 mm from the pial surface
    opts.under_node_or_face = 1; % 1 under node, 0 under face
    [dip_pos, dip_ori] = bst_generate_source_space(cortex_surface,opts);

    figure;
    plotmesh(dip_pos,Faces,'x>0');
    hold on;
    quiver3(dip_pos(:,1),dip_pos(:,2),dip_pos(:,3),dip_ori(:,1),dip_ori(:,2),dip_ori(:,3))