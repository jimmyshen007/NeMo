function [scanned_text] = ReadInTxt_special(file_source)
%ReadInTxt: Reads in an EOL-delimited file (ie, a list of filenames) and returns a char array

MSL=500;
scanned_text=[];
fid=fopen(file_source,'r');
tline=fgets(fid);
while ischar(tline)
    C=textscan(tline,['%*2c %' int2str(length(tline)-3) 'c']);
    scanned_text=[scanned_text; ['/home/eve/UCSF_data/Epilepsy_Subjects/Ashish_R01_Data' filesep C{1} blanks(MSL-length(tline)+1)]];
    disp(C{1});
    tline=fgets(fid);
end
fclose(fid);
return;
end

