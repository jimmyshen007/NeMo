function ArrToFile(inputArray, destination)
%ArrToFile.m takes care of outputting a list to a text file

sz=size(inputArray);
fid=fopen(destination,'a');

for i=1:sz
    StringOut(inputArray(i,:),fid);
end

fclose(fid);

end

