function DirToFile(read_dir, modif, fileout)
%DirToFile does a dir search and strips out the filenames from the
%resulting structs, outputs as EOL-delimited text

[B, sz]=DirToArr(read_dir, modif);

ArrToFile(B,fileout);

disp(['File out: ' fileout '.']);
open(fileout);
end

