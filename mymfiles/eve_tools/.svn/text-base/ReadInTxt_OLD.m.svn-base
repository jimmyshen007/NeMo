function [scanned_text] = ReadInTxt(file_source)
%ReadInTxt: Reads in an EOL-delimited file (ie, a list of filenames) and returns a char array
%Allows for commenting out of rows beginning with the '%' symbol

MSL=500;
scanned_text=[];
fid=fopen(file_source,'r');
tline=fgets(fid);
while ischar(tline)
    C=textscan(tline,['%c %*[^\n]']);
    if ~strcmp(C{1},'%')
        scanned_text=[scanned_text; tline blanks(5000-length(tline))];
    end
    tline=fgets(fid);
end
fclose(fid);
return;
end