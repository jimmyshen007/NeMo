function [scanned_text] = ReadInTxtPrefix(file_source,dir_tree)
%ReadInTxt: Reads in an EOL-delimited file (ie, a list of filenames) 
%which requires prefixing to create full path (usually from 'find' 
%results in Linux) and returns a char array

scanned_text=ReadInTxt(file_source,dir_tree);

return;
end

