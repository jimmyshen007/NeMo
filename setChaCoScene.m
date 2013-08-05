function [startupStruct, reschoice] = setChaCoScene(handles, fileIn)
% 
resultsFile = fileIn.ChaCoResultsFile;

load(resultsFile,'ChaCoResults'); 
[saveFolder, ~, ~] = fileparts(resultsFile);    

atlassize = size(ChaCoResults(1).Regions,2);

rescheck = [ 0 0 ]; %[ 86y/n 116y/n]

% Check for parallel file
switch atlassize
    case 86
        rescheck(1) = 1;
        parallelcase = [saveFolder filesep 'ChaCo116_MNI.mat'];
        set(handles.radiobutton1, 'Value', 1);
%         set(handles.popupmenu5, 'String', {'Positive/Negative';'Standard Lobes';'Functional 7-atlas'});
        set(handles.popupmenu5, 'String', {'Standard Lobes';'Functional 7-atlas'});
    case 116
        rescheck(2) = 1;
        parallelcase = [saveFolder filesep 'ChaCo86_MNI.mat'];
        set(handles.radiobutton2, 'Value', 1);
%         set(handles.popupmenu5, 'String', {'Positive/Negative';'Standard Lobes'});
        set(handles.popupmenu5, 'String', {'Standard Lobes'});
end

% Modify GUI
if exist(parallelcase, 'file')
    rescheck = [1 1];
    set(handles.radiobutton1, 'Enable', 'On');
    set(handles.radiobutton2, 'Enable', 'On');
else
    set(handles.radiobutton1, 'Enable', 'Off');
    set(handles.radiobutton2, 'Enable', 'Off');
end

% Create renderStructs
newStruct = importChaCoStruct(ChaCoResults(1), fileIn.plotOutputOptions);

set(handles.popupmenu6, 'Value', 1);
set(handles.edit2, 'String', newStruct(1).figstr);

if length(find(rescheck)) == 2
    load(parallelcase, 'ChaCoResults');
    parallelStruct = importChaCoStruct(handles, ChaCoResults(1), fileIn.plotOutputOptions);
    if atlassize == 86
        reschoice = 1;
        startupStruct = [newStruct; parallelStruct];
    else
        reschoice = 2;
        startupStruct = [parallelStruct; newStruct];
    end
else
    startupStruct = newStruct;
    reschoice = 1;
end


set(handles.pushbutton14,'BackgroundColor',newStruct(1).singleColor);
% 
% 
% % Set image properties
% 
% if plotOutputOptions.SurfPlot.flag
%     newStruct(1).opacity = 1;
% elseif plotOutputOptions.GBPlot.flag
%     newStruct(1).opacity = 0.08;
% end
% 
% 
% % Set render settings
% 
% 
% 
% 

% plotOutputOptions.SurfPlot.flag = 0;
% plotOutputOptions.SurfPlot.PlotHemi = 'left';
% plotOutputOptions.SurfPlot.MAP = [];
% plotOutputOptions.SurfPlot.movie = 0;
% plotOutputOptions.GBPlot.flag = 0;
% plotOutputOptions.GBPlot.movie = 0;
% plotOutputOptions.BoxPlot.flag =0;
% plotOutputOptions.GraphPlot.flag = 1;
% plotOutputOptions.GraphPlot.Global = 0;
% plotOutputOptions.GraphPlot.Local.flag = 1;
% plotOutputOptions.GraphPlot.Local.EF = 0;
% plotOutputOptions.GraphPlot.Local.BC = 0;
% plotOutputOptions.GraphPlot.Local.CP = 0;
% plotOutputOptions.GraphPlot.Local.MD = 1;
% plotOutputOptions.GraphPlot.Local.EC = 0;
% plotOutputOptions.GraphPlot.Local.CC = 0;
% plotOutputOptions.figsave = '_116_AD';
