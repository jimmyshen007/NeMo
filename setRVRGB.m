function setRVRGB(handles,rawMap)
minVal = get(handles.slider1,'Value');
maxVal = get(handles.slider2,'Value');
rvTable = get(handles.uitable1,'Data');
entriez = cell2mat(rvTable(:,2));
rgbMapping = getRGBTriple(rawMap,minVal,maxVal,entriez);
rvTable(:,3:5) = num2cell(rgbMapping);
set(handles.uitable1,'Data',rvTable);