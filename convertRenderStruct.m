function [TF, updatedStruct] = convertRenderStruct(oldStruct)

disp('Converting legacy format struct');

TF = 1;

updatedStruct = newRenderStruct;

oldFields = fields(oldStruct);
newFields = fields(updatedStruct);

for i=1:length(oldFields)
    if ismember(oldFields{i}, newFields)
        for j=1:length(oldStruct)
            eval(['updatedStruct(j).' oldFields{i} ' = oldStruct(j).' oldFields{i} ';']);
        end
    else
        TF = 0;
        disp(['Unable to convert struct: unknown field ' oldFields{i}]);
        break;
    end
end