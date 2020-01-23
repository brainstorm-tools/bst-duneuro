function bst_help_duneuro()
% This function will display the help the of the of the duneuro binaries
% usage : bst_help_duneuro


%% find the toolbox path
disp(['============================== '])
str = which('bst_dueneuro_readme.txt','-all');
[filepath,~,~] = fileparts(str{1});

% option 1
% if ispc
% system([fullfile(filepath,'bst_duneuro.exe ') '--help']);
% else
% system([fullfile(filepath,'./bst_duneuro ') '--help']);
% end
% ==> fails if the full path contains a space

% option 2
currentFolder = pwd;
cd(fullfile(filepath,'bst-duneuro','bin'))
if ispc
    system(['bst_duneuro.exe ' '--help']);
else
    system(['/.bst_duneuro ' '--help']);
end
cd(currentFolder)
end