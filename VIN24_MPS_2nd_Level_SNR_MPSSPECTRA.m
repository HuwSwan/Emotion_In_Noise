%% Output from VIN24_Flux_MPS (onl MPS) to look at effects of MPS SNR by signal at the points of detection.

ccc;
pTop = 'X:\PhD\03-Original_OIM'
pDat = fullfile(pTop,'02-analysis/2024-MPS_Flux');
pCells = 'X:\PhD\03-Original_OIM\02-analysis\2024-MPS_2nd_Level';
pOutTop = fullfile(pTop,'02-analysis/2024-MPS_2nd_Level_CorrectedSNR');

pIms = fullfile(pTop,'02-analysis/graphics');

pStim = fullfile(pTop,'04-stimuli');
types = {'normal','scram'};%



for iT = 1:length(types)

    disp(types{iT})
    prDat = fullfile(pDat,types{iT});

    pOutDat = fullfile(pOutTop,types{iT},'models');
    pOutPlot = fullfile(pOutTop,types{iT},'plots');

    mkdir(pOutDat);
    mkdir(pOutPlot);

    % Load data - loads emos and mps.
    % clear mps emos
    load(fullfile(prDat,'mps_data.mat'))

    %% plotting details

    B = max(mps{1}.MS.val(:))-(2*(std(mps{1}.MS.val(:))));
    C = min(mps{1}.MS.val(:))+(2*(std(mps{1}.MS.val(:))));

    xAx = mps{1}.MS.x;
    yAx = mps{1}.MS.y;

    %colum 1 signal. Colum 2, noise. Colum 3, combined.
    % Roguhness is 8 mamp mod by 15 freq mod.

    %1. MPS SNR
    %2. Combined MPS > noise MPS.


    emotions = {'neu','ang','hap'};

    %% MPS Spectra
    %Extract values for comparisons
    clear cellVals
    neuVals = cellfun(@(x) x.MS.val, mps(emos==1,:), 'UniformOutput', false);
    neuVals = cellfun(@(x) smoothdata2(x, 'loess', 10), neuVals, 'UniformOutput', false);

    cellVals{1}{1} = cat(3,neuVals{:,1}); % Signal
    cellVals{1}{2} = cat(3,neuVals{:,2}); % Noise
    cellVals{1}{3} = cat(3,neuVals{:,3}); % COmbined

    angVals = cellfun(@(x) x.MS.val, mps(emos==2,:), 'UniformOutput', false);
    angVals = cellfun(@(x) smoothdata2(x, 'loess', 10), angVals, 'UniformOutput', false);

    cellVals{2}{1} = cat(3,angVals{:,1});
    cellVals{2}{2} = cat(3,angVals{:,2});
    cellVals{2}{3} = cat(3,angVals{:,3});

    hapVals = cellfun(@(x) x.MS.val, mps(emos==3,:), 'UniformOutput', false);
    hapVals = cellfun(@(x) smoothdata2(x, 'loess', 10), hapVals, 'UniformOutput', false);

    cellVals{3}{1} = cat(3,hapVals{:,1});
    cellVals{3}{2} = cat(3,hapVals{:,2});
    cellVals{3}{3} = cat(3,hapVals{:,3});



    % % Save the smooth MPS spectra
    % mpsStr = fullfile(pCells,'individualMPS.mat');
    % save(mpsStr,'cellVals');


    mpsStr = fullfile(pCells,types{iT},'models','individualMPS.mat');
    load(mpsStr)

    %metohds/
    mStr = {'diff','ratio'};

    % PERMUTATION SETTING
    nP = 10^4;


    %% --------------------------- MODELS-----------------------------
    %% Differences and ratios -
    % Ratio is really most informative. Do I wanna do comp? hard to say.
    for iM = 2;%1:2
        return
        clear meanRoughRats ratCollatte angClusters hapClusters angHapClusters hapAngClusters
        for iR = 1:3 %Emotions.
            %SNR for each map
            clear matValues
            if iM == 1;
                % Just raw difference
                matValues = (cellVals{iR}{3})-(cellVals{iR}{2});

            else
                % ratio SNR db of Signal to Noise.
                % matValues = mag2db((abs(cellVals{iR}{1}))./(abs(cellVals{iR}{2})));

                holdSig = cellVals{iR}{1};
                holdNoise = cellVals{iR}{2};

                clipN = size(holdSig);

                %Normalise the data so it can be negative
                % Also adds then increases values by their own relative
                % minimum
                for iN = 1:clipN(3)

                    tempSig = holdSig(:,:,iN);
        
                    tempNoise = holdNoise(:,:,iN);
       

                    minVal = min([tempSig(:); tempNoise(:)]);
                    maxVal = max([tempSig(:); tempNoise(:)]);




                    A_normalized = (tempSig - minVal) / (maxVal - minVal);
                    B_normalized = (tempNoise - minVal) / (maxVal - minVal);


                    % RESCALE BY ADDING MIN. Rest doesn't work
                    scaledSig(:,:,iN) = rescale(tempSig)+mnTem;


                    scaledNoise(:,:,iN) = rescale(tempNoise)+mnTem;
                end
                %Scaled so we can do the ratios.
                matValues = mag2db(scaledSig./scaledNoise);                                                            

            end
            %For display image
            meanRoughRats = mean(matValues,3);
            %Combined stim MPS - noise sample MPS - For statistics.
            ratCollatte{iR} = matValues;
            results.meanOriginal{iR} = meanRoughRats;
        end

        dims = size(ratCollatte{1}(:,:,1));

        modVec = [2,1; 3,1; 2,3];
        modNames = {'AngNeu','HapNeu','AngHap'};

        % Through the models.
        for iE = 1:3

            cells1 = ratCollatte{modVec(iE,1)};
            cells2 = ratCollatte{modVec(iE,2)};

            [clusters, p, t, pDist ] = ...
                permutest(cells1, cells2,false,0.05,nP,false);

            clustMap = zeros(dims);
            clustMap(cat(1,clusters{:}))=1;           
            results.(modNames{iE}).clusters = clusters;
            results.(modNames{iE}).p = p;
            results.(modNames{iE}).t = t;
            results.(modNames{iE}).pDist = pDist;
            results.(modNames{iE}).meanFirst = mean(cells1,3);
            results.(modNames{iE}).meanSecond = mean(cells2,3);
            results.(modNames{iE}).clustMap = clustMap;
        end

        close all
        tiledlayout(3,3);

        %Mean values
        for iE = 1:3
            nexttile
            dat =  results.meanOriginal{iE};
            imagesc(dat); axis xy
        end

        nexttile
        for iE = 1:2
            nexttile
            dat = results.(modNames{iE}).clustMap;
            imagesc(dat); axis xy
        end

        %Maskedvalues > neu values
        nexttile
        for iE = 1:2
            nexttile
            check = results.(modNames{iE}).meanFirst>results.(modNames{iE}).meanSecond;
            map = results.(modNames{iE}).clustMap;


            temp = results.(modNames{iE}).meanFirst;
            temp(~check&~map) = nan;
            % Where correct is higher and significant.
            dat = temp;

            imagesc(dat); axis xy
        end


        % Save the results this time huw.
        resStr = fullfile(pOutDat,['results_',mStr{iM},'.mat']);
        save(resStr,'results');


        figStr = fullfile(pOutPlot,['mps_results_MPS_',mStr{iM},'.fig']);
        savefig(gcf,figStr)

        figStr = fullfile(pOutPlot,['mps_results_MPS_',mStr{iM},'.tiff']);
        saveas(gcf,figStr)

        clear results
        close all
    end


%% MPS Differences between the average points of detection.

%------------------- FOR THE STIMULI

clear meanRoughRats ratCollatte angClusters hapClusters angHapClusters hapAngClusters


%Significance models.
modVec = [2,1; 3,1; 2,3];
modNames = {'AngNeu','HapNeu','AngHap'};

% Through the models.
for iE = 1:3

    %Just the 1s clips of stimuli
    cells1 = cellVals{modVec(iE,1)}{1};
    cells2 = cellVals{modVec(iE,2)}{1};

    % TWO TAILED!
    [clusters, p, t, pDist ] = ...
        permutest(cells1, cells2,false,0.05,nP,true);

    clustMap = zeros(dims);
    clustMap(cat(1,clusters{:}))=1;
    resultsStim.(modNames{iE}).clusters = clusters;
    resultsStim.(modNames{iE}).p = p;
    resultsStim.(modNames{iE}).t = t;
    resultsStim.(modNames{iE}).pDist = pDist;
    resultsStim.(modNames{iE}).meanFirst = mean(cells1,3);
    resultsStim.(modNames{iE}).meanSecond = mean(cells2,3);
    resultsStim.(modNames{iE}).clustMap = clustMap;
end

close all
tiledlayout(3,3);

%Mean values
nexttile
dat =  resultsStim.AngNeu.meanSecond;
imagesc(dat); axis xy
set(gca,'CLim',[-5,1.9])
nexttile
dat =  resultsStim.AngNeu.meanFirst;
imagesc(dat); axis xy
set(gca,'CLim',[-5,1.9])
nexttile
dat =  resultsStim.HapNeu.meanFirst;
imagesc(dat); axis xy
set(gca,'CLim',[-5,1.9])



nexttile
for iE = 1:2
    nexttile
    dat = resultsStim.(modNames{iE}).clustMap;
    imagesc(dat); axis xy
end

%Maskedvalues > neu values
nexttile
for iE = 1:2
    nexttile
    check = resultsStim.(modNames{iE}).meanFirst>resultsStim.(modNames{iE}).meanSecond
    map = resultsStim.(modNames{iE}).clustMap;


    temp = resultsStim.(modNames{iE}).meanFirst;
    temp(~check&~map) = nan;
    % Where correct is higher and significant.
    dat = temp;

    imagesc(dat); axis xy
end


resStr = fullfile(pOutDat,['results_stimuli.mat']);
save(resStr,'resultsStim');


figStr = fullfile(pOutPlot,['mps_stimuli_',mStr{iM},'.fig']);
savefig(gcf,figStr)

figStr = fullfile(pOutPlot,['mps_stimuli_',mStr{iM},'.tiff']);
saveas(gcf,figStr)


end










return
%% Raw comparison of MPS - PROBABLY NO POINTS
%
% modVec = [2,1; 3,1; 2,3];
% modNames = {'AngNeu','HapNeu','AngHap'};
%
% % Through the models.
% for iE = 1:3
%
%     cells1 = cellVals{modVec(iE,1)}{3};
%     cells2 = cellVals{modVec(iE,2)}{3};
%     [clusters, p, t, pDist ] = ...
%         permutest(cells1, cells2,false,0.05,nP,true);
%
%     dims = size(cells2(:,:,1));
%
%
%     clustMap = zeros(dims);
%     clustMap(cat(1,clusters{:}))=1;
%     results.(modNames{iE}).clusters = clusters;
%     results.(modNames{iE}).p = p;
%     results.(modNames{iE}).t = t;
%     results.(modNames{iE}).pDist = pDist;
%     results.(modNames{iE}).meanFirstEmo = mean(cells1,3);
%     results.(modNames{iE}).meanSecondEmo = mean(cells2,3);
%     results.(modNames{iE}).clustMap = clustMap;
% end
%
%
% close all
% tiledlayout(3,3);
%
% %Mean values
% for iE = 1:3
%     nexttile
%     dat =   mean(cellVals{iE}{3},3);
%     imagesc(dat); axis xy
% end
%
% nexttile
% for iE = 1:2
%     nexttile
%     dat = results.(modNames{iE}).clustMap;
%     imagesc(dat); axis xy
% end
%
% %Maskedvalues > neu values
% nexttile
% for iE = 1:2
%     nexttile
%     check = results.(modNames{iE}).meanFirstEmo > results.(modNames{iE}).meanSecondEmo
%     map = results.(modNames{iE}).clustMap;
%
%
%     temp = results.(modNames{iE}).meanFirstEmo;
%     temp(~check&~map) = nan;
%     % Where correct is higher and significant.
%     dat = temp;
%
%     imagesc(dat); axis xy
% end
%
% figStr = fullfile(pOut,['mps_results_MPSCOMPS','.fig']);
% savefig(gcf,figStr)
%
% figStr = fullfile(pOut,['mps_results_MPSCOMPS','.tiff']);
% saveas(gcf,figStr)
%
% close all
%
%
%
%
%
%
%
%
%
%
%
%
%
%
%
%







%
%
%
%     % permutation testing
%     dims = size(ratCollatte{1}(:,:,1));
%
%     % ang to neutral.
%     [angClusters, p_ang, t_sums_ang, permutation_distribution_ang ] = ...
%         permutest(ratCollatte{2}, ratCollatte{1},false,0.05,10^4,false);
%
%     angClustMap = zeros(dims);
%     angClustMap(cat(1,angClusters{:}))=1;
%
%     % hap to neutral.
%     [hapClusters, p_hap, t_sums_hap, permutation_distribution_hap ] = ...
%         permutest(ratCollatte{3}, ratCollatte{1},false,0.05,10^4,false);
%
%     hapClustMap = zeros(dims);
%     hapClustMap(cat(1,hapClusters{:}))=1;
%
%     %Ang > Hap
%     [angHapClusters, p_AH, t_sums_hap, permutation_distribution_hap ] = ...
%         permutest(ratCollatte{2}, ratCollatte{3},false,0.05,10^4,false);
%
%     AHClustMap = zeros(dims);
%     AHClustMap(cat(1,angHapClusters{:}))=1;
%
%     %Hap > Ang
%     [hapAngClusters, p_AH, t_sums_hap, permutation_distribution_hap ] = ...
%         permutest(ratCollatte{3}, ratCollatte{2},false,0.05,10^4,false);
%
%     HAClustMap = zeros(dims);
%     HAClustMap(cat(1,hapAngClusters{:}))=1;
%
%
%     %Ang High than neu
%     angDif = (meanRoughRats{2}-meanRoughRats{1});
%     angDif(angDif<0 | angClustMap == 0)=nan;
%
%     %Hap higher than nue.
%     hapDif = (meanRoughRats{3}-meanRoughRats{2});
%     hapDif(hapDif<0 | hapClustMap == 0 )=nan;
%
%      %Ang High than HAp
%     AHDif = (meanRoughRats{2}-meanRoughRats{3});
%     AHDif(AHDif<0 | AHClustMap == 0)=nan;
%
%     %Hap higher than Ang.
%     HADif = (meanRoughRats{3}-meanRoughRats{2});
%     HADif(HADif<0 | HAClustMap == 0 )=nan;
%
%
%
%     %% Visualisation
%     close all
%
%     tiledlayout(5,3)
%     % MPS SNR
%     for iI = 1:3
%         nexttile
%         imagesc(meanRoughRats{iI}); axis xy
%     end
%
%     nexttile
%     nexttile
%     % Clusters
%     imagesc(angClustMap);axis xy
%     nexttile
%     imagesc(hapClustMap);axis xy
%
%     nexttile
%     nexttile
%     %Differences
%     imagesc(angDif);axis xy
%     nexttile
%     imagesc(hapDif);axis xy
%
%     nexttile
%     nexttile
%     % Clusters
%     imagesc(AHClustMap);axis xy
%     nexttile
%     imagesc(HAClustMap);axis xy
%
%     nexttile
%     nexttile
%     % Differences
%     imagesc(AHDif);axis xy
%     nexttile
%     imagesc(HADif);axis xy
%
%     figStr = fullfile(pIms,'mps',[types{iT},'_SNR_Signficance.tiff']);
%     saveas(gcf,figStr)
%
%     close all
%
%
%
%
% end
%
%
%     % Noise - COmbined Differences.
%
%
%     clear meanRoughRats ratCollatte
%     for iR = 1:3
%         %SNR for each map
%         matValues = ((abs(cellRough{iR}{3}).^2)-(abs(cellRough{iR}{2}).^2));
%         %For display image
%         meanRoughRats{iR} = mean(matValues,3);
%
%         %Combined stim MPS - noise sample MPS - For statistics.
%         ratCollatte{iR} = matValues;
%     end
%
%
%     % permutation testing
%     dims = size(ratCollatte{1}(:,:,1));
%
%     % ang to neutral.
%     [angClusters, p_ang, t_sums_ang, permutation_distribution_ang ] = ...
%         permutest(ratCollatte{2}, ratCollatte{1},false,0.05,10^4,false);
%
%     angClustMap = zeros(dims);
%     angClustMap(cat(1,angClusters{:}))=1;
%
%     % hap to neutral.
%     [hapClusters, p_hap, t_sums_hap, permutation_distribution_hap ] = ...
%         permutest(ratCollatte{3}, ratCollatte{1},false,0.05,10^4,false);
%
%     hapClustMap = zeros(dims);
%     hapClustMap(cat(1,hapClusters{:}))=1;
%
%     %Ang > Hap
%     [angHapClusters, p_AH, t_sums_hap, permutation_distribution_hap ] = ...
%         permutest(ratCollatte{2}, ratCollatte{3},false,0.05,10^4,false);
%
%     AHClustMap = zeros(dims);
%     AHClustMap(cat(1,angHapClusters{:}))=1;
%
%     %Hap > Ang
%     [hapAngClusters, p_AH, t_sums_hap, permutation_distribution_hap ] = ...
%         permutest(ratCollatte{3}, ratCollatte{2},false,0.05,10^4,false);
%
%     HAClustMap = zeros(dims);
%     HAClustMap(cat(1,hapAngClusters{:}))=1;
%
%
%     %Ang High than neu
%     angDif = (meanRoughRats{2}-meanRoughRats{1});
%     angDif(angDif<0 | angClustMap == 0)=nan;
%
%     %Hap higher than nue.
%     hapDif = (meanRoughRats{3}-meanRoughRats{2});
%     hapDif(hapDif<0 | hapClustMap == 0 )=nan;
%
%      %Ang High than HAp
%     AHDif = (meanRoughRats{2}-meanRoughRats{3});
%     AHDif(AHDif<0 | AHClustMap == 0)=nan;
%
%     %Hap higher than Ang.
%     HADif = (meanRoughRats{3}-meanRoughRats{2});
%     HADif(HADif<0 | HAClustMap == 0 )=nan;
%
%
%     figStr = fullfile(pIms,'mps',[types{iT},'_Difference_Signficance.tiff']);
%     saveas(gcf,figStr)
%     close all
%
%
%
%
% end
%
%
%
%





