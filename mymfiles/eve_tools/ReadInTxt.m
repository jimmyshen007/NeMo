function [scanned_text] = ReadInTxt(file_source,dir_tree,MaxStrLen)
%ReadInTxt: Reads in an EOL-delimited file (ie, a list of filenames) 
%which requires prefixing to create full path (usually from 'find' 
%results in Linux) and returns a char array

if nargin<3
    if nargin<2
        dir_tree='';
    end 
    MSL=500;
else
    MSL=MaxStrLen;
end

if exist(file_source,'file')~=2
    ftmp=dir_tree;
    dir_tree=file_source;
    file_source=ftmp;
end

scanned_text=[];
fid=fopen(file_source,'r');
tline=fgets(fid);
% while ischar(tline)
%     T=textscan(tline,'%c %*[^\n]');
%     if ~strcmp(T{1},'%')
%         if nargin==2
%             C=textscan(tline,['%*2c %' int2str(length(tline)-3) 'c']);
%             outline=[dir_tree filesep C{1}];
%             scanned_text=[scanned_text; [outline blanks(MSL-length(outline))]];
%             %disp(C{1});
%         else
%             scanned_text=[scanned_text; [tline blanks(MSL-length(tline))]];
%             %disp(tline);
%         end
%     end
%     tline=fgets(fid);
% end

while ischar(tline)
    T=textscan(tline,'%c %*[^\n]');
    if ~strcmp(T{1},'%')
        if nargin==2
            C=textscan(tline,['%*2c %' int2str(length(tline)-3) 'c']);
            outline=[dir_tree filesep C{1}];
            scanned_text=arrcat(scanned_text, outline, MSL);
            %disp(C{1});
        else
            scanned_text=arrcat(scanned_text, tline, MSL);
            %disp(tline);
        end
    end
    tline=fgets(fid);
end
    
fclose(fid);
return;
end

