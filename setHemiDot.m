function [hemidot, menuidx] = setHemiDot(atlassize,hemiSph)
% left116 = 1:2:107;
% right116 = 2:2:108;
% subcort116 = find(setLobeSchema(116, 2) == 5);
% left86 = [35:68 78:86];
% right86  = [1:34 68:77];
% subcort86 = find(setLobeSchema(86, 2) == 5);

if atlassize == 86
    if strcmp(hemiSph,'right')
        nanidx = [35:68 78:86]; % left86'
        menuidx = 3;
    elseif strcmp(hemiSph, 'left')
        nanidx = [1:34 69:77];
        menuidx = 2;
    elseif strcmp(hemiSph, 'subcort')
        nanidx =  find(setLobeSchema(86, 2) ~= 5);
        menuidx = 4;
    elseif strcmp(hemiSph, 'both')
        nanidx = ones(atlassize,1);
        menuidx = 1;
    end
elseif atlassize ==116
    if strcmp(hemiSph, 'right')
        nanidx = [1:2:107 109:116];
        menuidx = 3;
    elseif strcmp(hemiSph, 'left')
        nanidx = [2:2:108 109:116];
        menuidx = 2;
    elseif strcmp(hemiSph, 'subocrt')
        nanidx = find(setLobeSchema(86, 2) ~= 5);
        menuidx = 4;
    elseif strcmp(hemiSph, 'both')
        nanidx = ones(atlassize, 1);
        menuidx = 1;
    end
end
hemidot = ones(atlassize, 1);
hemidot(nanidx) = nan;