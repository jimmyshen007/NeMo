function newStruct = importChaCoStruct(ChaCoResults, plotOutputOptions)

newStruct = [newRenderStruct newRenderStruct];

atlassize = length(ChaCoResults(1).Regions);
newStruct(1).numberROI = atlassize;
atFileName = getAtlas(atlassize);
newStruct(1).brain_at = flipdim(spm_read_vols(spm_vol(atFileName)),1);

newStruct(1).regionvalues = num2cell(zeros(atlassize,5));
newStruct(1).regionvalues(:,1) = num2cell([1:atlassize]');
newStruct(1).regionvalues(:,2) = num2cell(ChaCoResults(1).Regions');
newStruct(1).regionvalues(:,3:5) = num2cell(getRGBTriple(hot(200), min(ChaCoResults(1).Regions), max(ChaCoResults(1).Regions), ChaCoResults(1).Regions'));
newStruct(1).brain_colormaprange = [min(ChaCoResults(1).Regions) max(ChaCoResults(1).Regions)];
% colorMapOpts = {'autumn'; 'bone'; 'cool'; 'hot'; 'jet'; 'spring'; 'summer'; 'winter'};
newStruct(1).brain_colormapidx = 4;


% if GraphPlot.Local.flag == 1
        %figure;
        %imagesc(LoCoResults(1).nConMat)
%         if GraphPlot.Local.DG==1
%             DEG = sum(LoCoResults(1).nConMat);
%             brain_network_ploting(pth,at, zeros(size(LoCoResults(1).nConMat)), DEG', [], lobes, 'GB_LocDeg',[],0,strsave)
%         end
%         if GraphPlot.Local.EF == 1
            Eff = efficiency(ChaCoResults(1).nConMat,1);
            sEff = Eff-0.9*min(Eff(find(Eff)));
            sEff = sEff./max(sEff);
            sEff(find(sEff<0)) = 0;
%             brain_network_ploting(pth,at, zeros(size(LoCoResults(1).nConMat)), sEff, [], lobes, 'GB_LocEff',[],0,strsave)
%         end
%         if GraphPlot.Local.BC == 1
            BCs = betweenness_wei(1./ChaCoResults(1).nConMat);
%             brain_network_ploting(pth,at, zeros(size(LoCoResults(1).nConMat)), BCs, [], lobes, 'GB_BetwCent',[],0,strsave)
%            
%         end
%         if GraphPlot.Local.CP == 1
            D=distance_wei(1./ChaCoResults(1).nConMat);
            CPs = sum(D)';
            sCPs = CPs-0.9*min(CPs(find(CPs)));
            sCPs = sCPs./max(sCPs);
            sCPs(find(sCPs<0)) = 0;
%             brain_network_ploting(pth,at, zeros(size(LoCoResults(1).nConMat)), sCPs, [], lobes, 'GB_AvShortPath',[],0,strsave)
%         end
%         if GraphPlot.Local.MD == 1
            MD = modularity_und(ChaCoResults(1).nConMat);
%             brain_network_ploting(pth,at, zeros(size(LoCoResults(1).nConMat)), 0.5*ones(length(LoCoResults(1).nConMat),1), [], MD, 'GB_Modularity',[],0,strsave)
%         end
%         if GraphPlot.Local.EC == 1
            [~,~,Ecc,~,~] = charpath(D);
            sEcc = Ecc-0.9*min(Ecc(find(Ecc)));
            sEcc = sEcc./max(sEcc);
            sEcc(find(sEcc<0)) = 0;
%             brain_network_ploting(pth,at, zeros(size(LoCoResults(1).nConMat)), sEcc, [], lobes, 'GB_Ecc',[],0,strsave)
%         end
%         if GraphPlot.Local.CC == 1 %clustering coefficient
            CC=clustering_coef_wu(ChaCoResults(1).nConMat);
%             brain_network_ploting(pth,at, zeros(size(LoCoResults(1).nConMat)), CC, [], lobes, 'GB_ClustCoef',[],0,strsave)
%         end
%     end
% newStruct(1).pipeColorHyperCube = ChaCoResults(1).nLocalMets;
newStruct(1).pipeColorHyperCube = [sEff'; BCs'; sCPs'; MD'; sEcc'; CC'];

% % Set image properties
if plotOutputOptions.SurfPlot.flag
    newStruct(1).opacity = 1;
    newStruct(1).singleColorFlag = 0;
    
    % Functionality = hemisphere display choice
    [newStruct(1).custom_colormap, newStruct(2).custom_colormap] = setHemiDot(atlassize,plotOutputOptions.SurfPlot.PlotHemi);


elseif plotOutputOptions.GraphPlot.Local.flag
    newStruct(1).opacity = 0.08;
    newStruct(1).singleColorFlag = 1;
    newStruct(1).nodes = 1;
    
    
    localPlot = [plotOutputOptions.GraphPlot.Local.EF plotOutputOptions.GraphPlot.Local.BC ...
        plotOutputOptions.GraphPlot.Local.CP plotOutputOptions.GraphPlot.Local.MD ...
        plotOutputOptions.GraphPlot.Local.EC plotOutputOptions.GraphPlot.Local.CC];
    
%     
%     for i = 1:length(localPlot)
%         res = localPlot(i);
%         if res
%             res = i;    
%             break;
%         end
%     end

    res = find(localPlot,1);
    newStruct(1).nodeStyle = res;
    newStruct(1).nodeProps = num2cell(zeros(atlassize,3));
    
    newStruct(1).nodeProps(:,1) = num2cell([1:atlassize]');
    newStruct(1).nodeProps(:,2) = num2cell(newStruct(1).pipeColorHyperCube(res,:)');
    
    [nodeSchema, colorSchema] = setLobeSchema(atlassize,plotOutputOptions.plotlobecolor,newStruct(1).pipeColorHyperCube(res,:)');
    newStruct(1).nodeProps(:,3) = num2cell(nodeSchema);    
    newStruct(1).nodeSchema = colorSchema;
    
    
    [newStruct(1).custom_colormap, newStruct(2).custom_colormap] = setHemiDot(atlassize,'both');

end

newStruct(1).brain_colormaprange = [min(ChaCoResults.Regions) max(ChaCoResults.Regions)];

newStruct(2).nodeSchema = plotOutputOptions.plotlobecolor;
newStruct(2).renderRes = 1;
newStruct(2).saveImages = 1;
newStruct(2).figstr = ['ChaCo_' num2str(atlassize) '_'];

% plotOutputOptions.SurfPlot.flag = 0;
% plotOutputOptions.SurfPlot.PlotHemi = 'left';
% plotOutputOptions.SurfPlot.MAP = [];
% plotOutputOptions.SurfPlot.movie = 0;
% plotOutputOptions.GBPlot.flag = 1;
% plotOutputOptions.GBPlot.movie = 0;
% plotOutputOptions.BoxPlot.flag =0;
% plotOutputOptions.GraphPlot.flag = 1;
% plotOutputOptions.GraphPlot.Global = 0;
% plotOutputOptions.GraphPlot.Local.flag = 1;
% plotOutputOptions.GraphPlot.Local.EF = 0;
% plotOutputOptions.GraphPlot.Local.BC = 0;
% plotOutputOptions.GraphPlot.Local.CP = 0;
% plotOutputOptions.GraphPlot.Local.MD = 1;
% plotOutputOptions.GraphPlot.Local.EC = 0;
% plotOutputOptions.GraphPlot.Local.CC = 0;
% plotOutputOptions.figsave = '_116_AD';
% plotOutputOptions.plotlobecolor = 1 %,2 or 3

% 
% newStruct=struct('volString','','brain_at',[],'dim',[],'mat',[],'opacity',1.0, ...
%     'regionvalues',[],'singleColorFlag',0,'singleColor',[0.4792    0.5625    0.5625],'brain_colormap','bone','custom_colormap',[], ...
%     'brain_colormapidx',1,'brain_colormaprange',[],'nodes',0,'pipes',0, 'connectivityMatrix',[], ...
%     'nodeScale',2.5,'nodeSchema',[],'nodeProps',[],'nodeStyle',1,'pipeScale',1.5, ...
%     'pipeScheme',1,'pipeColorHyperCube',[],'pipeCouplet',[rand(1,3);rand(1,3)], ...
%     'pipeColorMap','jet','pipeCoupletThreshold',50,'pipeStyle',1,'pipeUniform',0,'renderRes',2, ...
%     'currentVol',1,'saveImages',1,'saveMovie',0,'figstr','NeMo','mainHandle',[],'numberROI',[]);

