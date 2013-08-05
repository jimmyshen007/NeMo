function varargout = NeMo_GUI(varargin)
% NEMO_GUI MATLAB code for NEMO_GUI.fig
%      NeMo_GUI, by itself, creates a new NEMO_GUI or raises the existing
%      singleton*.
%
%      H = NEMO_GUI returns the handle to a new NEMO_GUI or the handle to
%      the existing singleton*.
%
%      NEMO_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in NEMO_GUI.M with the given input arguments.
%
%      NEMO_GUI('Property','Value',...) creates a new NEMO_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before NEMO_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to NEMO_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help NeMo_GUI

% Last Modified by GUIDE v2.5 31-May-2013 00:36:41

% Begin initialization code - DO NOT EDIT

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @NeMo_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @NeMo_GUI_OutputFcn, ...
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


% --- Executes just before NeMo_GUI is made visible.
function NeMo_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to NeMo_GUI (see VARARGIN)
% if( ismcc || isdeployed)
%     addpath(genpath('mymfiles'));
% end
global dispvar disptype sMRItype atlas outputtype damagefile saveFolder referenceFile
damagefile = 0;
saveFolder = 0;
referenceFile = 0;
if ~(ismcc() || isdeployed())
    addpath(genpath(['mymfiles' filesep 'spm8']));
end
% Choose default command line output for NeMo_GUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes NeMo_GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = NeMo_GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

        
% --- Executes during object creation, after setting all properties.
function uipanel1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uipanel1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
        set(hObject,'Visible','off');
        
% --- Executes during object creation, after setting all properties.
function uipanel2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uipanel2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
        set(hObject,'Visible','off');
% --------------------------------------------------------------------

    
% --- Check if the file is a valid file type %    
function [valid_flg] = isValidFileType(filePath)
    valid_flg = 0;
    [fp, fn, fext] = fileparts(filePath);
    switch lower(fext)
        case '.nii'
            valid_flg = 1;
        case '.hdr'
            valid_flg = 1;
        case '.mat'
            valid_flg = 1;
    end

% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
        %global fullFilePath save_flg saveFolder conventional_flg spatial_flg graph_flg muT1 muT2 muS1 muS2...
        %        spa_autoparam_flg grf_autoparam_flg autoMCSS autob1cor defaultSaveFolder referenceFile TE_ui isEchoTimeSet noIter autoSS
        global dispvar disptype sMRItype barplotoption atlas outputtype damagefile saveFolder referenceFile
        if ~isnumeric(damagefile)
            if exist(damagefile, 'file')
                if isValidFileType(damagefile)
                    if isnumeric(saveFolder)
                            if get(handles.checkbox_saveResultDefault, 'Value')
                                [fp, fn, ext] = fileparts(damagefile);
                                saveFolder = [fp filesep 'default_output'];
                            end
                            if ~isdir(saveFolder)
                                mkdir(saveFolder);
                            end
                    else
                        if ~isdir(saveFolder)
                            mkdir(saveFolder);
                        end
                    end
                    fileDim = getFileStat(damagefile);
                    if (length(fileDim) == 3) 
                        if (fileDim(3) > 1)
                                if isnumeric(referenceFile)
                                    coregMNI = [];
                                else
                                    switch sMRItype
                                        case 1
                                            coregMNI.StructImageType = 't1';
                                        case 2
                                            coregMNI.StructImageType = 't2';
                                        case 3
                                            coregMNI.StructImageType = 'epi';
                                        otherwise
                                    end
                                    coregMNI.ImageFileName = referenceFile;
                                end         
                                ChaCoCalc(damagefile, coregMNI, 'MNI', atlas, saveFolder, 1, 0, 0);
                                if exist([saveFolder filesep 'ChaCo' num2str(atlas) '_MNI.mat'], 'file')
                                    set(handles.pushbuttonV, 'Enable', 'on');
                                    commPlot.flag = 1;
                                    switch disptype(1)
                                        case 1 
                                            commPlot.PlotHemi = 'both';
                                        case 2
                                            commPlot.PlotHemi = 'left';
                                        case 3
                                            commPlot.PlotHemi = 'right';
                                        case 4
                                            commPlot.PlotHemi = 'subcortical';
                                        otherwise

                                    end
                                    commPlot.movie = disptype(2);
                                    switch dispvar(1)
                                        case 1 
                                            commPlot.Global = 1;
                                            commPlot.Local.flag = 0;
                                        case 0
                                            commPlot.Global = 0;
                                            commPlot.Local.flag = 1;
                                            if dispvar(2)
                                                commPlot.Local.CP = 1;
                                            else
                                                commPlot.Local.CP = 0;
                                            end
                                            if dispvar(3)
                                                commPlot.Local.EF = 1;
                                            else
                                                commPlot.Local.EF = 0;
                                            end
                                            if dispvar(4)
                                                commPlot.Local.DG = 1;
                                            else
                                                commPlot.Local.DG = 0;
                                            end
                                            if dispvar(5)
                                                commPlot.Local.EC = 1;
                                            else
                                                commPlot.Local.EC = 0;
                                            end
                                            if dispvar(6)
                                                commPlot.Local.MD = 1;
                                            else
                                                commPlot.Local.MD = 0;
                                            end
                                            if dispvar(7)
                                                commPlot.Local.BC = 1;
                                            else
                                                commPlot.Local.BC = 0;
                                            end
                                            if dispvar(8)
                                                commPlot.Local.CC = 1;
                                            else
                                                commPlot.Local.CC = 0;
                                            end
                                        otherwise

                                    end
                                    commPlot.MAP = [];
                                    switch outputtype
                                        case 1
                                            GBPlot = commPlot;
                                            SurfPlot.flag = 0; 
                                            BoxPlot.flag = 0;
                                            GraphPlot = commPlot;
                                            if GraphPlot.Global == 1
                                                GraphPlot.Global = 0;
                                            end
                                        case 2
                                            GBPlot.flag = 0;
                                            SurfPlot = commPlot; 
                                            BoxPlot.flag = 0;
                                            GraphPlot = commPlot;
                                            if GraphPlot.Global == 1
                                                GraphPlot.Global = 0;
                                            end
                                        case 3
                                            GBPlot.flag = 0;
                                            SurfPlot.flag = 0; 
                                            GraphPlot.flag = 0;       
                                            if barplotoption(2)
                                                GraphPlot.flag = 1;
                                                GraphPlot.Global = 1;
                                                GraphPlot.Local.flag = 0;
                                            end
                                            if barplotoption(1)
                                                BoxPlot = commPlot;
                                            else
                                                BoxPlot.flag = 0;
                                            end
                                        case 4
                                            GBPlot.flag = 0;
                                            SurfPlot.flag = 0;
                                            BoxPlot.flag = 0;
                                            GraphPlot.flag = 0;
                                            writeExcel(saveFolder, atlas, [saveFolder filesep 'ChaCo' num2str(atlas) '_MNI']);
                                        otherwise
                                    end
                                    figsave = ['_' num2str(atlas) '_AD'];
                                    %PlotChaCoResults([saveFolder filesep 'ChaCo' num2str(atlas) '_MNI'], GBPlot,SurfPlot,BoxPlot,GraphPlot, figsave, 1);
                                else
                                    set(handles.pushbuttonV, 'Enable', 'off');
                                end
                        else
                            modaldlg('Title','Warning', 'String', 'Can not run, 3-dimensional mult-slice image needed');
                        end
                    else
                        modaldlg('Title','Warning', 'String', 'Can not run, 3-dimensional image needed');
                    end
                else
                     modaldlg('Title','Warning', 'String', 'Invalid file selected');
                end
            else
                modaldlg('Title','Warning', 'String', 'File selected not found');
            end
        else
            modaldlg('Title','Warning', 'String', 'No file selected');
        end


% --------------------------------------------------------------------
function file_menu_Callback(hObject, eventdata, handles)
% hObject    handle to file_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function view_menu_Callback(hObject, eventdata, handles)
% hObject    handle to view_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function tools_menu_Callback(hObject, eventdata, handles)
% hObject    handle to tools_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function help_menu_Callback(hObject, eventdata, handles)
% hObject    handle to help_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function about_item_Callback(hObject, eventdata, handles)
% hObject    handle to about_item (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



    
% --------------------------------------------------------------------
function imageView_item_Callback(hObject, eventdata, handles)
% hObject    handle to imageView_item (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
   spm_image('init', '');
              

    
function [valid_flg] = isValidSaveType(filePath)
    valid_flg = 0;
    [fp, fn, fext] = fileparts(filePath);
    switch lower(fext)
        case '.nii'
            valid_flg = 1;
        case '.hdr'
            valid_flg = 1;
    end


%  ---- get file extension -----------
function [ext] = getExtension(filename)
    [fp, fn, ext] = fileparts(filename);

%  ----  get file statistics ------    
function [sz] = getFileStat(filename)
    try
        ext = getExtension(filename);
        switch lower(ext)
            case '.mat'
                v = load(filename);
                sz = size(v);
            case '.nii'
                v = spm_vol(filename);
                sz = [v(1).dim length(v)];
            case '.hdr'
                v = spm_vol(filename);
                sz = v(1).dim;
        end        
    catch
        sz = 0;
    end
    

    

% --------------------------------------------------------------------
function exit_item_Callback(hObject, eventdata, handles)
% hObject    handle to exit_item (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    delete(handles.figure1);

    

% --- Executes on button press in checkbox_saveResultDefault.
function checkbox_saveResultDefault_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_saveResultDefault (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_saveResultDefault
    

% --- Executes during object creation, after setting all properties.
function checkbox_saveResultDefault_CreateFcn(hObject, eventdata, handles)
% hObject    handle to checkbox_saveResultDefault (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
    
        
        
function [isEV] = checkFileExistValid(fullFilePath)
    isEV = 0;
    if ~isnumeric(fullFilePath)
        if exist(fullFilePath, 'file')
            if isValidFileType(fullFilePath)
                isEV = 1;
            else
                 modaldlg('Title','Warning', 'String', 'Invalid file selected');
            end
        else
            modaldlg('Title','Warning', 'String', 'File selected not found');
        end
    else
        modaldlg('Title','Warning', 'String', 'No file selected');
    end        




% --- Executes when figure1 is resized.
function figure1_ResizeFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    currentFFPos = get(handles.editwm, 'Position');
    windowPos = get(hObject, 'Position');
    currentFFPos(3) = (84.5/140.66) * windowPos(3);
    if currentFFPos(3) >= 84.5
        set(handles.editwm, 'Position', currentFFPos);
    else
        currentFFPos(3) = 84.5;
        set(handles.editwm, 'Position', currentFFPos);
    end
    currentFFPos = get(handles.editsMRI, 'Position');
    windowPos = get(hObject, 'Position');
    currentFFPos(3) = (52.5/140.66) * windowPos(3);
    if currentFFPos(3) >= 52.5
        currentFFPos(3) = currentFFPos(3)+(currentFFPos(3)*(32/140.66));
        set(handles.editsMRI, 'Position', currentFFPos);
    else
        currentFFPos(3) = 52.5;
        set(handles.editsMRI, 'Position', currentFFPos);
    end


% --------------------------------------------------------------------
function OpenImage_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to OpenImage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function editwm_Callback(hObject, eventdata, handles)
% hObject    handle to editwm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editwm as text
%        str2double(get(hObject,'String')) returns contents of editwm as a double



% --- Executes during object creation, after setting all properties.
function editwm_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editwm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editsMRI_Callback(hObject, eventdata, handles)
% hObject    handle to editsMRI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editsMRI as text
%        str2double(get(hObject,'String')) returns contents of editsMRI as a double


% --- Executes during object creation, after setting all properties.
function editsMRI_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editsMRI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbuttonbrwm.
function pushbuttonbrwm_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonbrwm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    [filename, pathname] = uigetfile( ...
    {'*.hdr','header Files (*.hdr)';
    '*.nii',  'NIFTI files (*.nii)'}, ...
    'Choose a T1 reference image file');
    global damagefile
    old_damagefile = damagefile;
    damagefile = [pathname filename];
    if ~isnumeric(damagefile)
        if ~checkFileExistValid(damagefile)
            damagefile = 0;
            set(handles.editwm, 'String', '');
        else
            set(handles.editwm, 'String', damagefile);
            set(handles.pushbuttonViewWM, 'Enable', 'on');
        end
    else
        damagefile = old_damagefile;
    end
    


% --- Executes on button press in pushbuttonbrsMRI.
function pushbuttonbrsMRI_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonbrsMRI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    global referenceFile
    [filename, pathname] = uigetfile( ...
    {'*.nii',  'NIFTI files (*.nii)';
    '*.hdr','header Files (*.hdr)'}, ...
    'Choose a T1 reference image file');
    old_referenceFile = referenceFile; 
    referenceFile = [pathname filename];
    if ~isnumeric(referenceFile)
        if ~checkFileExistValid(referenceFile)
            referenceFile = 0;
            set(handles.editsMRI, 'String', 'None');
        else
            set(handles.editsMRI, 'String', referenceFile);
            set(handles.pushbuttonViewsMRI, 'Enable', 'on');
        end
    else
        referenceFile = old_referenceFile;
    end

% --- Executes on button press in T1radio.
function T1radio_Callback(hObject, eventdata, handles)
% hObject    handle to T1radio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of T1radio


% --- Executes on button press in T2radio.
function T2radio_Callback(hObject, eventdata, handles)
% hObject    handle to T2radio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of T2radio


% --- Executes on button press in EPIradio.
function EPIradio_Callback(hObject, eventdata, handles)
% hObject    handle to EPIradio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of EPIradio


% --- Executes on button press in radioAAL.
function radioAAL_Callback(hObject, eventdata, handles)
% hObject    handle to radioAAL (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radioAAL


% --- Executes on button press in radiobuttonFS.
function radiobuttonFS_Callback(hObject, eventdata, handles)
% hObject    handle to radiobuttonFS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobuttonFS


% --- Executes on button press in radiobuttonGlass.
function radiobuttonGlass_Callback(hObject, eventdata, handles)
% hObject    handle to radiobuttonGlass (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobuttonGlass


% --- Executes on button press in radiobuttonGummi.
function radiobuttonGummi_Callback(hObject, eventdata, handles)
% hObject    handle to radiobuttonGummi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobuttonGummi


% --- Executes on button press in radiobuttonBar.
function radiobuttonBar_Callback(hObject, eventdata, handles)
% hObject    handle to radiobuttonBar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobuttonBar


% --- Executes on button press in radiobuttonExcel.
function radiobuttonExcel_Callback(hObject, eventdata, handles)
% hObject    handle to radiobuttonExcel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobuttonExcel


% --- Executes on button press in radiobuttonLH.
function radiobuttonLH_Callback(hObject, eventdata, handles)
% hObject    handle to radiobuttonLH (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobuttonLH
    


% --- Executes on button press in radiobuttonMV.
function radiobuttonMV_Callback(hObject, eventdata, handles)
% hObject    handle to radiobuttonMV (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobuttonMV
    

% --- Executes on button press in radiobuttonWB.
function radiobuttonWB_Callback(hObject, eventdata, handles)
% hObject    handle to radiobuttonWB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobuttonWB
    

% --- Executes on button press in radiobuttonRH.
function radiobuttonRH_Callback(hObject, eventdata, handles)
% hObject    handle to radiobuttonRH (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobuttonRH
    


% --- Executes when selected object is changed in uipanelDV.
function uipanelDV_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanelDV 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
%get(eventdata.NewValue, 'Tag')
global dispvar
tmptag = get(eventdata.NewValue, 'Tag');
switch tmptag
    case 'radiobuttonCIC'
        set(handles.uipanelCha, 'visible', 'off');
        dispvar(1) = 1;
    otherwise
        set(handles.uipanelCha, 'visible', 'on');
        dispvar(1) = 0;
        
end

% --- Executes when selected object is changed in uipanelCha.
function uipanelCha_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanelCha 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
% global dispvar
% tmptag = get(eventdata.NewValue, 'Tag');
% switch tmptag
%     case 'checkboxASPL'
%         dispvar(2) = get(handles.checkboxASPL, 'Value');
%     case 'checkboxEff'
%         dispvar(3) = get(handles.checkboxEff, 'Value');
%     case 'checkboxDeg'
%         dispvar(4) = get(handles.checkboxDeg, 'Value');  
%     case 'checkboxEcc'
%         dispvar(5) = get(handles.checkboxEcc, 'Value');          
%     case 'checkboxMod'
%         dispvar(6) = get(handles.checkboxMod, 'Value');  
%     case 'checkboxBet'
%         dispvar(7) = get(handles.checkboxBet, 'Value');
%     otherwise
% end


% --- Executes when selected object is changed in uipanelsMRI.
function uipanelsMRI_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanelsMRI 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
global sMRItype
tmptag = get(eventdata.NewValue, 'Tag');
switch tmptag
    case 'T1radio'
        sMRItype = 1;
    case 'T2radio'
        sMRItype = 2;
    case 'EPIradio'
        sMRItype = 3;        
    otherwise
end


% --- Executes when selected object is changed in uipanelCA.
function uipanelCA_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanelCA 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
global atlas
tmptag = get(eventdata.NewValue, 'Tag');
switch tmptag
    case 'radioAAL'
        atlas = 116;
    case 'radiobuttonFS'
        atlas = 86;
    otherwise
end

% --- Executes when selected object is changed in uipanelCOT.
function uipanelCOT_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanelCOT 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
global outputtype
tmptag = get(eventdata.NewValue, 'Tag');
switch tmptag
    case 'radiobuttonGlass'
        set(handles.uipanelPar, 'visible', 'on');
        set(handles.uipanelBarplot, 'visible', 'off');
        set(handles.radiobuttonLNM, 'Enable', 'on');
        outputtype = 1;
    case 'radiobuttonGummi'
        set(handles.uipanelPar, 'visible', 'on');
        set(handles.uipanelBarplot, 'visible', 'off');
        set(handles.radiobuttonLNM, 'Enable', 'off');
        outputtype = 2;
    otherwise    
        set(handles.uipanelPar, 'visible', 'off');
        if strcmp(tmptag, 'radiobuttonBar')
            set(handles.uipanelBarplot, 'visible', 'on');
            outputtype = 3;
        else
            outputtype = 4;
        end    
end

% --- Executes when selected object is changed in uipanelDT.
function uipanelDT_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanelDT 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
global disptype
tmptag = get(eventdata.NewValue, 'Tag');
switch tmptag
    case 'radiobuttonWB'
        disptype(1) = 1;
    case 'radiobuttonLH'
        disptype(1) = 2;
    case 'radiobuttonRH'
        disptype(1) = 3;
    case 'radiobuttonSub'
        disptype(1) = 4;
    case 'checkboxMov'
        disptype(2) = get(handles.checkboxMov, 'Value');
    otherwise
end

% --- Executes during object creation, after setting all properties.
function uipanelDV_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uipanelDV (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
global dispvar
dispvar = [1 0 0 0 0 0 0 0];


% --- Executes during object creation, after setting all properties.
function radiobuttonLNM_CreateFcn(hObject, eventdata, handles)
% hObject    handle to radiobuttonLNM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
set(hObject, 'Value', 1);


% --- Executes during object creation, after setting all properties.
function uipanelCA_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uipanelCA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
global atlas
atlas = 116;


% --- Executes during object creation, after setting all properties.
function uipanelCOT_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uipanelCOT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
global outputtype
outputtype = 1;


% --- Executes during object creation, after setting all properties.
function uipanelDT_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uipanelDT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
global disptype
disptype = [1 0];


% --- Executes during object creation, after setting all properties.
function uipanelsMRI_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uipanelsMRI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
global sMRItype
sMRItype = 1;


% --------------------------------------------------------------------
function uipushtoolSave_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to uipushtoolSave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    global saveFolder
    saveFolder_old = saveFolder;
    [saveFolder pathname] = uigetdir('', 'Create a file to save processed data');
    if isnumeric(saveFolder)
        saveFolder = saveFolder_old;
    else
        saveFolder = [pathname filesep saveFolder];
        
    end


% --- Executes on button press in pushbuttonViewWM.
function pushbuttonViewWM_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonViewWM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    global damagefile
    if ~isnumeric(damagefile)
        if exist(damagefile, 'file')
            if isValidSaveType(damagefile)
                %try
                    spm_image('init', damagefile);
                %catch
                %    modaldlg('Title','Warning', 'String', 'Image file damaged');
                %end
            else
                modaldlg('Title','Warning', 'String', 'Can not view this type of image file');
            end
        else
            modaldlg('Title','Warning', 'String', 'Invalid image file');
        end
    else
        modaldlg('Title','Warning', 'String', 'No image file selected');
    end


% --- Executes on button press in pushbuttonViewsMRI.
function pushbuttonViewsMRI_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonViewsMRI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    global referenceFile
    if ~isnumeric(referenceFile)
        if exist(referenceFile, 'file')
            if isValidSaveType(referenceFile)
                %try
                    spm_image('init', referenceFile);
                %catch
                %    modaldlg('Title','Warning', 'String', 'Image file damaged');
                %end
            else
                modaldlg('Title','Warning', 'String', 'Can not view this type of image file');
            end
        else
            modaldlg('Title','Warning', 'String', 'Invalid image file');
        end
    else
        modaldlg('Title','Warning', 'String', 'No image file selected');
    end


% --- Executes when selected object is changed in uipanelBarplot.
function uipanelBarplot_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanelBarplot 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
% global barplotoption
% tmptag = get(eventdata.NewValue, 'Tag');
% switch tmptag
%     case 'radiobuttonDS'
%         barplotoption = 1;
%     case 'radiobuttonGNM'
%         barplotoption = 2;
%     otherwise      
% end

% --- Executes during object creation, after setting all properties.
function uipanelBarplot_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uipanelBarplot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
global barplotoption
barplotoption = [0 0];


% --- Executes on button press in checkboxMov.
function checkboxMov_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxMov (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxMov
        global disptype
        disptype(2) = get(handles.checkboxMov, 'Value');


% --- Executes on button press in checkboxASPL.
function checkboxASPL_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxASPL (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxASPL
        global dispvar
        dispvar(2) = get(handles.checkboxASPL, 'Value');
        
    

% --- Executes on button press in checkboxEff.
function checkboxEff_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxEff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxEff
        global dispvar
        dispvar(3) = get(handles.checkboxEff, 'Value');
        

% --- Executes on button press in checkboxEcc.
function checkboxEcc_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxEcc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxEcc
        global dispvar
        dispvar(5) = get(handles.checkboxEcc, 'Value');          
        

% --- Executes on button press in checkboxMod.
function checkboxMod_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxMod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxMod
        global dispvar
        dispvar(6) = get(handles.checkboxMod, 'Value');  

% --- Executes on button press in checkboxBet.
function checkboxBet_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxBet (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxBet
        global dispvar
        dispvar(7) = get(handles.checkboxBet, 'Value');


% --- Executes on button press in checkboxDeg.
function checkboxDeg_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxDeg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxDeg
        global dispvar
        dispvar(4) = get(handles.checkboxDeg, 'Value');  


% --- Executes on button press in checkboxDS.
function checkboxDS_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxDS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxDS
global barplotoption
barplotoption(1) = get(hObject,'Value');

% --- Executes on button press in checkboxGNM.
function checkboxGNM_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxGNM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxGNM
global barplotoption
barplotoption(2) = get(hObject,'Value');


% --- Executes on button press in pushbuttonV.
function pushbuttonV_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonV (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    global dispvar disptype barplotoption atlas outputtype saveFolder damagefile
    if isnumeric(saveFolder)
        if get(handles.checkbox_saveResultDefault, 'Value')
            [fp, fn, ext] = fileparts(damagefile);
            saveFolder = [fp filesep 'default_output'];
        end
    end
    if exist([saveFolder filesep 'ChaCo' num2str(atlas) '_MNI.mat'], 'file')
        set(handles.pushbuttonV, 'Enable', 'on');
        commPlot.flag = 1;
        switch disptype(1)
            case 1 
                commPlot.PlotHemi = 'both';
            case 2
                commPlot.PlotHemi = 'left';
            case 3
                commPlot.PlotHemi = 'right';
            case 4
                commPlot.PlotHemi = 'subcortical';
            otherwise

        end
        commPlot.movie = disptype(2);
        switch dispvar(1)
            case 1 
                commPlot.Global = 1;
                commPlot.Local.flag = 0;
            case 0
                commPlot.Global = 0;
                commPlot.Local.flag = 1;
                if dispvar(2)
                    commPlot.Local.CP = 1;
                else
                    commPlot.Local.CP = 0;
                end
                if dispvar(3)
                    commPlot.Local.EF = 1;
                else
                    commPlot.Local.EF = 0;
                end
                if dispvar(4)
                    commPlot.Local.DG = 1;
                else
                    commPlot.Local.DG = 0;
                end
                if dispvar(5)
                    commPlot.Local.EC = 1;
                else
                    commPlot.Local.EC = 0;
                end
                if dispvar(6)
                    commPlot.Local.MD = 1;
                else
                    commPlot.Local.MD = 0;
                end
                if dispvar(7)
                    commPlot.Local.BC = 1;
                else
                    commPlot.Local.BC = 0;
                end
                if dispvar(8)
                    commPlot.Local.CC = 1;
                else
                    commPlot.Local.CC = 0;
                end
            otherwise

        end
        commPlot.MAP = [];
        switch outputtype
            case 1
                GBPlot = commPlot;
                SurfPlot.flag = 0; 
                BoxPlot.flag = 0;
                GraphPlot = commPlot;
                if GraphPlot.Global == 1
                    GraphPlot.Global = 0;
                end
            case 2
                GBPlot.flag = 0;
                SurfPlot = commPlot; 
                BoxPlot.flag = 0;
                GraphPlot = commPlot;
                if GraphPlot.Global == 1
                    GraphPlot.Global = 0;
                end
            case 3
                GBPlot.flag = 0;
                SurfPlot.flag = 0; 
                GraphPlot.flag = 0;       
                if barplotoption(2)
                    GraphPlot.flag = 1;
                    GraphPlot.Global = 1;
                    GraphPlot.Local.flag = 0;
                end
                if barplotoption(1)
                    BoxPlot = commPlot;
                else
                    BoxPlot.flag = 0;
                end
            case 4
                GBPlot.flag = 0;
                SurfPlot.flag = 0;
                BoxPlot.flag = 0;
                GraphPlot.flag = 0;
                writeExcel(saveFolder, atlas, [saveFolder filesep 'ChaCo' num2str(atlas) '_MNI']);
            otherwise
        end
        figsave = ['_' num2str(atlas) '_AD'];
        if GBPlot.flag || SurfPlot.flag || commPlot.Local.flag
            fileIn.ChaCoResultsFile = [saveFolder filesep 'ChaCo' num2str(atlas) '_MNI'];
            fileIn.plotOutputOptions.SurfPlot = SurfPlot;
            fileIn.plotOutputOptions.GraphPlot = GraphPlot;
            fileIn.plotOutputOptions.GBPlot = GBPlot;
            fileIn.plotOutputOptions.plotlobecolor = 1;
            fileIn.plotOutputOptions.figstr = 'NeMo';
            brainography_nemo_lite(fileIn);
            SurfPlot.flag = 0;
            GBPlot.flag = 0;
            GraphPlot.Local.flag = 0;
        end
        %%%%%%%%%%%
        if BoxPlot.flag || GraphPlot.Global
            PlotChaCoResults([saveFolder filesep 'ChaCo' num2str(atlas) '_MNI'], GBPlot,SurfPlot,BoxPlot,GraphPlot, figsave, 1);
        end

        %PlotChaCoResults([saveFolder filesep 'ChaCo' num2str(atlas) '_MNI'], GBPlot,SurfPlot,BoxPlot,GraphPlot, figsave, 1);
    else
        set(handles.pushbuttonV, 'Enable', 'off');
    end

% --- Executes on button press in checkboxCC.
function checkboxCC_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxCC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxCC
        global dispvar
        dispvar(8) = get(handles.checkboxCC, 'Value');
