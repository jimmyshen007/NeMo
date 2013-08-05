function [nodeSchema, colorSchema] = setLobeSchema(atlassize,schemaOpt,nodeVals)

% schemaOpt 1==lobe colors, 2==7atlas

if nargin < 3
    nodeVals = ones(atlassize,1);
end

negPos = 3;
lobeColors = 1;
sevenAt = 2;

if schemaOpt == negPos
    nodeSchema = zeros(atlassize, 1);
    nodeSchema(find(nodeVals < 0)) = 1;
    nodeSchema(find(nodeVals > 0)) = 2; 
    colorSchema = [0 0 1; 1 0 0];
elseif atlassize == 116 && schemaOpt ~= negPos
    nodeSchema = zeros(116,1);
    nodeSchema(1:28) = 1;
    nodeSchema([29:40, 71:78]) = 5;
    nodeSchema(43:54) = 3;
    nodeSchema([55:56, 79:90]) = 4;
    nodeSchema(57:70) = 2;
    nodeSchema(91:116) = 6;
    % new mods in limbic system
    nodeSchema([29, 30]) = 1;  % insula
    nodeSchema([31, 32]) = 1;  % ant cing
    nodeSchema([33, 34]) = 1; % mid cing
    nodeSchema([35, 36]) = 2; % post cing
    nodeSchema([37, 38]) = 5; % hippo
    nodeSchema([39, 40]) = 4; % parahippo 
    nodeSchema([41, 42]) = 4; % parahippo
    colorSchema = [0 0 1; %blue: Frontal           
                1 0 1; %magenta: Parietal
                0 1 0; %green: Temporal
                1 0 0; %red: Occipital
                0 1 1; %cyan: Cingulate 
                0 0 0];%black: everyewhere else
elseif atlassize == 86 && schemaOpt == sevenAt
    load(['resource' filesep 'Convert_86to7atlas.mat'], 'Functional86_roi');
    nodeSchema = Functional86_roi;
    nodeSchema(find(nodeSchema==0)) = 8;
    colorSchema =   [0 0 1; %blue: Visual           
                1 0 1; %magenta: Somatomotor
                0 1 0; %green: Dorsal Attn
                1 0 0; %red: Ventral Attn
                0 1 1; %cyan: Limbic  
                1 1 0; %yellow: Fronto-pareital 
                0 0 0; %black: Default Mode   
                1 0.5 0]; %orange Cerebellar/SubCort
elseif atlassize ==86 && schemaOpt == lobeColors
% lobeColors86
lobeColors86 = zeros(1,86);
    % Frontal (blue)
    a = [3 11 13 16 17 18 19 23 26 27 31];
    lobeColors86(a) = 1;
    lobeColors86(a + 34) = 1;

    % Parietal (magenta)
    a =[7 21 24 28 30]; 
    lobeColors86(a) = 2;
    lobeColors86(a + 34) = 2;

    % Temporal (green)
    a = [1 5 6 8 14 15 29 32 33];
    lobeColors86(a) = 3;
    lobeColors86(a + 34) = 3;
% handles
    % Occipital (red)
    a = [4 10 12 20];
    lobeColors86(a) = 4;
    lobeColors86(a + 34) = 4;

    % Cingulate (cyan)
    a = [2 9 22 25];
    lobeColors86(a) = 5;
    lobeColors86(a + 34) = 5;

    % Everywhere else (black)
    ilobe = (lobeColors86 == 0);
    lobeColors86(ilobe) = 6;

    nodeSchema = lobeColors86';
    colorSchema = [0 0 1; %blue: Frontal           
                1 0 1; %magenta: Parietal
                0 1 0; %green: Temporal
                1 0 0; %red: Occipital
                0 1 1; %cyan: Cingulate 
                0 0 0];%black: everyewhere else
end
