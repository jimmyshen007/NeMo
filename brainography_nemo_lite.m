function varargout = brainography_nemo_lite(varargin)
dbstop if error
% BRAINOGRAPHY_NEMO_LITE MATLAB code for brainography_nemo_lite.fig
%      BRAINOGRAPHY_NEMO_LITE, by itself, creates a new BRAINOGRAPHY_NEMO_LITE or raises the existing
%      singleton*.
%
%      H = BRAINOGRAPHY_NEMO_LITE returns the handle to a new BRAINOGRAPHY_NEMO_LITE or the handle to
%      the existing singleton*.
%
%      BRAINOGRAPHY_NEMO_LITE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in BRAINOGRAPHY_NEMO_LITE.M with the given input arguments.
%
%      BRAINOGRAPHY_NEMO_LITE('Property','Value',...) creates a new BRAINOGRAPHY_NEMO_LITE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before brainography_nemo_lite_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to brainography_nemo_lite_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help brainography_nemo_lite

% Last Modified by GUIDE v2.5 22-May-2013 13:21:31

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @brainography_nemo_lite_OpeningFcn, ...
    'gui_OutputFcn',  @brainography_nemo_lite_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);


if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before brainography_nemo_lite is made visible.
function brainography_nemo_lite_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to brainography_nemo_lite (see VARARGIN)

% Choose default command line output for brainography_nemo_lite
handles.output = hObject;
% setappdata(handles.output,'mainHandle',handles.output);
fileIn = varargin{1};

colorMapOpts = {'autumn'; 'bone'; 'cool'; 'hot'; 'jet'; 'spring'; 'summer'; 'winter'};
set(handles.popupmenu3,'String',colorMapOpts);

[startupStruct, res] = setChaCoScene(handles, fileIn);
startupStruct(1,2).mainHandle = hObject;
if size(startupStruct,1) > 1
    for i = 1:2
        atlassize = startupStruct(i,1).numberROI;
        savestr = ['renderStruct' num2str(atlassize)];
        setappdata(hObject,savestr,startupStruct(i,:));
    end
end

renderStruct = startupStruct(res,:);

populateGUI(handles, renderStruct);

guidata(hObject, renderStruct);

function populateGUI(H, popVal)

atlassize = popVal(1).numberROI;

if atlassize == 86
    set(H.radiobutton1,'Value',1);
    set(H.radiobutton2,'Value',0);
else
    set(H.radiobutton1,'Value',0);
    set(H.radiobutton2,'Value',1);
end

set(H.pushbutton14,'BackgroundColor',popVal(1).singleColor);

if popVal(1).singleColorFlag
    enableGlassBrainOpts(atlassize, H);
    set(H.radiobutton3,'Value',1);
    set(H.radiobutton4,'Value',0);
    set(H.popupmenu4,'Value',popVal(1).nodeStyle); % Choice of statistic EF/EC/etc
else
    set(H.radiobutton3,'Value',0);
    set(H.radiobutton4,'Value',1); %Gummi option
    enableGummiBrainOpts(atlassize, H);
end

set(H.popupmenu6,'Value',popVal(2).custom_colormap); % hemisphere
set(H.popupmenu5,'Value',popVal(2).nodeSchema); % Lobe schema menu
set(H.popupmenu3,'Value',popVal(1).brain_colormapidx); % Gummi colormap option

set(H.pushbutton14,'BackgroundColor',popVal(1).singleColor);

set(H.edit2,'String',popVal(1).figstr);
set(H.checkbox2,'Value',popVal(2).saveImages);
set(H.checkbox3,'Value',popVal(2).saveMovie);

% --- Outputs from this function are returned to the command line.
function varargout = brainography_nemo_lite_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version ostore struct in guif MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
% varargout{1} = handles.output;
varargout{1}=hObject;
% varargout{2}=handles; %If I can find some way to output volume struct to
% command line
% disp(hObject);


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
disp(guihandles(hObject));
disp(handles);
% EXECUTE FULL-SCALE IMAGE GENERATION USING GUIDATA(1:end-1)
figure;
handles(1).nodeStyle = 1;
H = guihandles(hObject);
disp(handles(1).singleColor);
if get(H.radiobutton4,'Value')
    regionvals = cell2mat(handles(1).regionvalues(:,2)).*handles(1).custom_colormap;
    handles(1).regionvalues(:,2) = num2cell(regionvals);
end
BrainographyRender(handles,gca,1);

% --- Executes on button press in checkbox2.  Save Images Option
function checkbox2_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox2
handles(end).saveImages=get(hObject,'Value');
guidata(hObject,handles);

% --- Executes on button press in checkbox3. Save Movie Option
function checkbox3_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox3
handles(end).saveMovie=get(hObject,'Value');
guidata(hObject,handles);

function resetMyAxes(H)
ac = allchild(H.axes1);
for i = 1:size(ac,1)
    delete(ac(i));
end
view([1 0 0]);
clmo(handlem('light'));
if get(H.checkbox5,'Value')
    set(H.checkbox5,'Value',0);
end

% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

H = guihandles(hObject);
handles(end).saveMovie = 0;
handles(end).saveImages = 0;
%cla reset;
resetMyAxes(H);

if length(handles) > 1
    handles(1).nodeStyle = 1;
    if get(H.radiobutton4,'Value')
        regionvals = cell2mat(handles(1).regionvalues(:,2)).*handles(1).custom_colormap;
        handles(1).regionvalues(:,2) = num2cell(regionvals);
    end
    BrainographyRender(handles,H.axes1,3);
%     view([1 1 0]);
end

% --- Executes on button press in pushbutton11.
function pushbutton11_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
view([1 0 0]);

% --- Executes on button press in pushbutton12.
function pushbutton12_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
view([0 0 1]);

% --- Executes on button press in pushbutton13.
function pushbutton13_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
view([0 1 0]);


% --- Executes on button press in checkbox5.
function checkbox5_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox5
clval = get(hObject, 'Value');
clmo(handlem('light'))
if clval
    camlight left;
end
    

function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double
handles(1).figstr = get(hObject, 'String');

guidata(hObject,handles);




% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in pushbutton14.
function pushbutton14_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
V = uisetcolor(get(hObject,'BackgroundColor'));
% disp(V);
set(hObject,'BackgroundColor',V);
handles(1).singleColor = V;
handles(1).regionvalues(:,3:5) = num2cell(repmat(V,handles(1).numberROI,1));
guidata(hObject, handles);

% 
% currentVol = handles(end).currentVol;
% handles(currentVol).singleColor = V;
% regionvalues = handles(currentVol).regionvalues;
% numberROI = size(regionvalues,1);
% regionvalues(:,3:5) = num2cell(repmat(V,numberROI,1));
% handles(currentVol).regionvalues = regionvalues;
% guidata(hObject,handles);


% --- Executes on selection change in popupmenu3.
function popupmenu3_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu3 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu3
H = guihandles(hObject);
colorMapStr = get(H.popupmenu3,'String');
colorMapIdx = get(H.popupmenu3,'Value');

colorMapSel = colorMapStr{colorMapIdx};
eval(['CMap = ' colorMapSel '(200);']);

entriez = cell2mat(handles(1).regionvalues(:,2));
entriezRGB = getRGBTriple(CMap, min(entriez), max(entriez), entriez);
handles(1).regionvalues(:,3:5) = num2cell(entriezRGB);
handles(1).brain_colormapidx = colorMapIdx;
handles(1).brain_colormap = colorMapSel;
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function popupmenu3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu4.
function popupmenu4_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu4 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu4
H = guihandles(hObject);
GB= get(H.radiobutton3,'Value');

if GB
    statChoice = get(H.popupmenu4,'Value');
    stat = num2cell(handles(1).pipeColorHyperCube(statChoice,:)');
    handles(1).nodeProps(:,2) = stat;
    handles(1).nodeStyle = statChoice;
    guidata(hObject, handles);
end




% --- Executes during object creation, after setting all properties.
function popupmenu4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu5.
function popupmenu5_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu5 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu5
H = guihandles(hObject);
schemaVal = get(H.popupmenu5,'Value');
schemaStr = get(H.popupmenu5,'String');
atlassize = handles(1).numberROI;

[nodeSchema, colorSchema] = setLobeSchema(atlassize, schemaVal);
handles(1).nodeProps(:,3) = num2cell(nodeSchema);
handles(1).nodeSchema = colorSchema;
handles(2).nodeSchema = schemaVal;
guidata(hObject,handles);
% 
% function setLobeSchema(schemaChoice)
% H = getappdata(handles.output,'renderStruct');
% 
% if 
% schema = get

% --- Executes during object creation, after setting all properties.
function popupmenu5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in pushbutton15.
function pushbutton15_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
resetMyAxes(guihandles(hObject));

% --- Executes when selected object is changed in uipanel6.
function uipanel6_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel6 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
H = guihandles(hObject);
switch get(eventdata.NewValue, 'Tag') % Get Tag of selected object.
    case 'radiobutton1' %86 regions
        % Code for when radiobutton1 is selected.
        LobeStrings = {'Positive/Negative';'Standard Lobes';'Functional 7-atlas'};
        set(handles.popupmenu5, 'Enable', 'On');
        set(handles.popupmenu5, 'Value', 1);
        set(handles.popupmenu5, 'String', LobeStrings);
        setappdata(H.figure1, 'renderStruct116', handles);
        handles = getappdata(H.figure1, 'renderStruct86');
        populateGUI(H, handles);
        guidata(hObject, handles);
    case 'radiobutton2'
        % Code for when radiobutton2 is selected.
        set(handles.popupmenu5, 'Value', 1);
        LobeStrings = {'Positive/Negative'; 'Standard Lobes'};
        set(handles.popupmenu5, 'String', LobeStrings);
        set(handles.popupmenu5, 'Enable', 'Off');
        setappdata(H.figure1,'renderStruct86', handles);
        handles = getappdata(H.figure1, 'renderStruct116');
        populateGUI(H, handles);
        guidata(hObject, handles);
    otherwise
        disp('Oh no.');
        % Code for when there is no match.
end

% 
% function nodeSchema = getNodeSchema(lobeChoice)
% if lobeChoice == 3
% 	nodeSchema =   [0 0 1; %blue: Visual           
%                 1 0 1; %magenta: Somatomotor
%                 0 1 0; %green: Dorsal Attn
%                 1 0 0; %red: Ventral Attn
%                 0 1 1; %cyan: Limbic  
%                 1 1 0; %yellow: Fronto-pareital 
%                 0 0 0; %black: Default Mode   
%                 1 0.5 0]; %orange Cerebellar/SubCort
% elseif lobeChoice == 2
%     nodeSchema = [0 0 1; %blue: Frontal           
%                 1 0 1; %magenta: Parietal
%                 0 1 0; %green: Temporal
%                 1 0 0; %red: Occipital
%                 0 1 1; %cyan: Cingulate 
%                 0 0 0];%black: everyewhere else
% elseif lobeChoice == 1
%     nodeSchema = [0 0 1; 1 0 0];
% end

% --- Executes when selected object is changed in uipanel7.
function uipanel7_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel7 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
H = guihandles(hObject);
atlassize = handles(1).numberROI;

switch get(eventdata.NewValue,'Tag') % Get Tag of selected object.
    case 'radiobutton3'
        % Code for when radiobutton3 is selected. 'Glass Brain'
        enableGlassBrainOpts(atlassize, H);
        handles(1).opacity = 0.08;
        handles(1).singleColorFlag = 1;
        handles(1).nodes = 1;
    case 'radiobutton4'
        % Code for when radiobutton4 is selected. 'Gummi Brain'
        enableGummiBrainOpts(atlassize, H);
        handles(1).opacity = 1.0;
        handles(1).singleColorFlag = 0;
        handles(1).nodes = 0;
    otherwise
        disp('Oh no.');
        % Code for when there is no match.
end
guidata(hObject, handles);

function enableGlassBrainOpts(atlassize, H)

set(H.popupmenu5,'Enable','On');

set(H.pushbutton14,'Enable','On');

set(H.popupmenu3,'Enable','Off');

statOpts = {'Efficiency'; 'Betweenness Centrality'; 'Avg Shortest Path Length'; 'Modularity'; 'Eccentricity'; 'Clustering Coefficient'};
set(H.popupmenu4,'String',statOpts);
set(H.popupmenu4,'Value',1);
set(H.popupmenu6,'Enable','Off');
% set


function enableGummiBrainOpts(atlassize, H)

set(H.popupmenu5,'Enable','Off'); %Node schema off

set(H.pushbutton14,'Enable','Off');

set(H.popupmenu3,'Enable','On'); %Colormap choice

set(H.popupmenu4,'String',{'ChaCo Score'});
set(H.popupmenu4,'Value',1);

set(H.popupmenu6,'Enable','On');


% --- Executes on selection change in popupmenu6.
function popupmenu6_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu6 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from
%        popupmenu6
H = guihandles(hObject);

atval = get(H.radiobutton1,'Value');
if atval
    atlassize = 86;
else
    atlassize = 116;
end
hemival = get(H.popupmenu6,'Value');

switch hemival
    case 1
        hemidot = ones(atlassize,1);
    case 2
        hemidot = setHemiDot(atlassize,'left');
    case 3
        hemidot = setHemiDot(atlassize,'right');
    case 4
        hemidot = setHemiDot(atlassize,'subcort');
end
handles(1).custom_colormap = hemidot;
handles(2).custom_colormap = hemival;
guidata(hObject, handles);
% setappdata(H.figure1, 'hemidot', hemidot);


% --- Executes during object creation, after setting all properties.
function popupmenu6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
