function bst_help_duneuro_meeg()
% This function will display the help the of the of the duneuro binaries
% usage : bst_help_duneuro_meeg


%% find the toolbox path
disp(['============================== '])
str = which('bst_unique_readme.txt','-all');
[filepath,~,~] = fileparts(str{1});

% option 2
currentFolder = pwd;
cd(fullfile(filepath,'bin'))
if ispc
    system(['bst_duneuro_meeg.exe ' '--help']);
else
    system(['/.bst_duneuro_meeg ' '--help']);
end
cd(currentFolder)
end