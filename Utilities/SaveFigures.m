% Save all current figures into certain folder
%
% Inputs
%   n/a
%
% Outputs
%   n/a
%
% Derived from code written by:
% - https://www.mathworks.com/matlabcentral/answers/385220-create-a-subfolder-with-datetime-as-its-name-using-mkdir
% - Dr. Dryan J. Caverly
%
% Written by:   Keegan Bunker 
%               PhD student, University of Minnesota 
%               bunke029@umn.edu
%
% Last major modifications: May 31, 2023
%-------------------------------------------------------------------------%
% TODO:
%   -   Define all settings and options in a seperate seciton -   Use said
%   seciton to create an optional input to change those settings and change
%   this to an actual function that can be ran in a workspace without input
%   or with these optional inputs defining the settings like directory,
%   file names, size, etc.
%   - Add the ability to change font sizes
%-------------------------------------------------------------------------%

% Pick where to save figures
% currentFolder = pwd;
% dname = uigetdir(currentFolder);
dname = pwd;

currDate = strrep(datestr(datetime), ':', '_');
FolderName =  ['Figures_',currDate];
mkdir(dname,FolderName)

% Save figures name just in order
FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
for iFig = 1:length(FigList)
    FigHandle = FigList(iFig);
    FigName   = num2str(get(FigHandle, 'Number'));
    set(0, 'CurrentFigure', FigHandle);
    
    % Set size and position
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % width = 20;
    % height = 10;
    % 
    % 
    % set(gcf, 'PaperPositionMode', 'manual');
    % set(gcf,'paperunits','centimeters')
    % set(gcf,'papersize',[width,height])
    % set(gcf,'paperposition',[0,0,width,height])
    % set(gcf,'renderer','Painters')
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    %print(FigHandle, fullfile(dname, FolderName,strcat(FigName)),'-depsc')
    saveas(FigHandle, fullfile(dname, FolderName,strcat(FigName)),'epsc')
    % print(FigHandle, strcat(FigName),'-depsc')
    % fullfile(dname, FolderName,strcat(FigName))

    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    saveas(FigHandle, fullfile(dname, FolderName,strcat(FigName, '.fig'))); % specify the full path
    saveas(FigHandle, fullfile(dname, FolderName,strcat(FigName, '.png'))); % specify the full path
end