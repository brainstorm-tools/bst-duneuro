function [node,elem,face,allMask] = mesh_by_iso2mesh(fullPathToT1,fullPathToT2)

subjRSPD =fullPathToT1;
T2 = fullPathToT2;
[dirname,baseFilename] = fileparts(subjRSPD);

%%  STEP 3: MESH GENERATION...
% see cgalv2m for more information

model = 2;
if model == 0,    maxvol = 50; reratio = 3; radbound = 5; angbound = 30; distbound = 0.4; end
if model == 1,    maxvol = 10; reratio = 3; radbound = 5; angbound = 30; distbound = 0.4; end
if model == 2,    maxvol =   5; reratio = 3; radbound = 5; angbound = 30; distbound = 0.4; end
if model == 3,    maxvol =   1; reratio = 3; radbound = 5; angbound = 30; distbound = 0.4; end
 
saveMeshFormatMat = 1;
saveMeshFormatMsh= 0;

meshOpt = struct('radbound',radbound,'angbound',angbound,...
                             'distbound',distbound,'reratio',reratio,...
                             'maxvol',maxvol,'saveMeshFormatMat',saveMeshFormatMat,...
                             'saveMeshFormatMsh',saveMeshFormatMsh);
% uniqueTag = ['reratio_', num2str(reratio),'_maxvol_',num2str(maxvol)];
uniqueTag = ['MeshModel_', num2str(maxvol),'_',num2str(reratio)...
                                            '_',num2str(radbound), '_',num2str(angbound) , '_',num2str(distbound)];

if ~exist([dirname filesep baseFilename '_' uniqueTag '.mat'],'file')
    disp('======================================================')
    disp('        STEP 3: MESH GENERATION...         ')
    disp('======================================================')
    [node,elem,face,allMask] = meshByIso2meshWithoutElectrode(subjRSPD,subjRSPD,T2,meshOpt,[],uniqueTag);
else
    disp('======================================================')
    disp('          MESH ALREADY GENERATED, SKIP STEP 3         ')
    disp('======================================================')
    load([dirname filesep baseFilename '_' uniqueTag '.mat'],'node','elem','face');
end

end
