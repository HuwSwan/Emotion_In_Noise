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
    clear mps emos
    % load(fullfile(prDat,'mps_data.mat'))
    % 
    % %% plotting details
    % 
    % B = max(mps{1}.MS.val(:))-(2*(std(mps{1}.MS.val(:))));
    % C = min(mps{1}.MS.val(:))+(2*(std(mps{1}.MS.val(:))));
    % 
    % xAx = mps{1}.MS.x;
    % yAx = mps{1}.MS.y;
    % 
    % save(fullfile(pOutDat,'axes_values.mat'),'xAx','yAx')

    %colum 1 signal. Colum 2, noise. Colum 3, combined.
    % Roguhness is 8 mamp mod by 15 freq mod.

    %1. MPS SNR
    %2. Combined MPS > noise MPS.


    emotions = {'neu','ang','hap'};

    %% Correct for the Log Transform done in the MPS toolbox.

    %parrelise
    % 
    % clear cellVals
    % parfor iP = 1:3
    %     neuVals = cellfun(@(x) exp(x.MS.val), mps(emos==iP,:), 'UniformOutput', false);
    %     neuValSmooth = cellfun(@(x) smoothdata2(x, 'loess', 10), neuVals, 'UniformOutput', false);
    %     neuValAdj = cellfun(@(x,y) max(x,min(y(:))), neuValSmooth, neuVals, 'UniformOutput', false);
    % 
    %     cellVals{iP}{1} = cat(3,neuValAdj{:,1}); % Signal
    %     cellVals{iP}{2} = cat(3,neuValAdj{:,2}); % Noise
    %     cellVals{iP}{3} = cat(3,neuValAdj{:,3}); % COmbined
    % end
    % 


    %% MPS Spectra
    %Extract values for comparisons
    % clear cellVals
    % neuVals = cellfun(@(x) exp(x.MS.val), mps(emos==1,:), 'UniformOutput', false);
    % neuValSmooth = cellfun(@(x) smoothdata2(x, 'loess', 10), neuVals, 'UniformOutput', false);
    % neuValAdj = cellfun(@(x,y) max(x,min(y(:))), neuValSmooth, neuVals, 'UniformOutput', false);
    % 
    % cellVals{1}{1} = cat(3,neuValAdj{:,1}); % Signal
    % cellVals{1}{2} = cat(3,neuValAdj{:,2}); % Noise
    % cellVals{1}{3} = cat(3,neuValAdj{:,3}); % COmbined
    % 
    % angVals = cellfun(@(x) exp(x.MS.val), mps(emos==2,:), 'UniformOutput', false);
    % angValSmooth = cellfun(@(x) smoothdata2(x, 'loess', 10), angVals, 'UniformOutput', false);
    % angValAdj = cellfun(@(x,y) max(x,min(y(:))), angValSmooth, angVals, 'UniformOutput', false);
    % 
    % cellVals{2}{1} = cat(3,angValAdj{:,1});
    % cellVals{2}{2} = cat(3,angValAdj{:,2});
    % cellVals{2}{3} = cat(3,angValAdj{:,3});
    % 
    % hapVals = cellfun(@(x) exp(x.MS.val), mps(emos==3,:), 'UniformOutput', false);
    % hapValSmooth = cellfun(@(x) smoothdata2(x, 'loess', 10), hapVals, 'UniformOutput', false);
    % hapValsAdju = cellfun(@(x,y) max(x,min(y(:))), hapValSmooth, hapVals, 'UniformOutput', false);
    % 
    % cellVals{3}{1} = cat(3,hapVals{:,1});
    % cellVals{3}{2} = cat(3,hapVals{:,2});
    % cellVals{3}{3} = cat(3,hapVals{:,3});

    % 
    % 
    % % % Save the smooth MPS spectra
    mpsStr = fullfile(pOutDat,'individualMPS.mat');
    % save(mpsStr,'cellVals');

    load(mpsStr)

    %metohds -
    % Difference - Difference in raw analytic power (combined - noise) then
    % log transform trhat.
    %SNR - log transform ratio.
    mStr = {'ratio','diff',''};

    % PERMUTATION SETTING
    nP = 10^4;


    %% --------------------------- MODELS-----------------------------
    %% Differences and ratios -
    % Ratio is really most informative. Do I wanna do comp? hard to say.
    for iM = 1:2

        clear meanRoughRats ratCollatte angClusters hapClusters angHapClusters hapAngClusters
        for iR = 1:3 %Emotions.
            %SNR for each map
            clear matValues
            if iM == 2;
                % Just raw difference
                % This exclusively looks where combined is higher than
                % noise alone.
                matValues = (cellVals{iR}{3})-(cellVals{iR}{2});
                matValues(matValues<0)=nan;
                matValues = log(matValues);     
                matValues(isnan(matValues))=0;

            else
                % ratio SNR db of Signal to Noise.
                % matValues = mag2db((abs(cellVals{iR}{1}))./(abs(cellVals{iR}{2})));

                holdSig = cellVals{iR}{1};
                holdNoise = cellVals{iR}{2};

                %Log transformed.
                matValues = mag2db(holdSig./holdNoise);       
            end
            %For display image
            meanRoughRats = nanmean(matValues,3);
            %Combined stim MPS - noise sample MPS - For statistics.
            ratCollatte{iR} = matValues;
            results.meanOriginal{iR} = meanRoughRats;
        end

        dims = size(ratCollatte{1}(:,:,1));

        modVec = [2,1; 3,1; 2,3];
        modNames = {'AngNeu','HapNeu','AngHap'};

        %predefine for parralelisation
        % Through the models.
        parfor iE = 1:3

            cells1 = ratCollatte{modVec(iE,1)};
            cells2 = ratCollatte{modVec(iE,2)};

            [clusters, p, t, pDist ] = ...
                permutest(cells1, cells2,false,0.05,nP,false);

            clustMap = zeros(dims);
            clustMap(cat(1,clusters{:}))=1;            

            clustersHold{iE} = clusters;
            pHold{iE} = p;
            tHold{iE} = t;
            pDistHold{iE} = pDist;
            meanFirstHold{iE} = nanmean(cells1,3);
            meanSecondHold{iE} = nanmean(cells2,3);
            clustMapHold{iE} = clustMap;
        end

        for iE = 1:3
            results.(modNames{iE}).clusters = clustersHold{iE};
            results.(modNames{iE}).p = pHold{iE};
            results.(modNames{iE}).t = tHold{iE};
            results.(modNames{iE}).pDist = pDistHold{iE};
            results.(modNames{iE}).meanFirst = meanFirstHold{iE};
            results.(modNames{iE}).meanSecond = meanSecondHold{iE};
            results.(modNames{iE}).clustMap = clustMapHold{iE};
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
% 
% for iE = 1:3
%     holdMean = cellVals{iE}{1};
%     rawStimuli{iE} = nanmean(log(holdMean),3);
% end
% 
% stStri = fullfile(pOutDat,['stimuli_mps_values.mat']);
% save(stStri,'rawStimuli');


    %Significance models.
    modVec = [2,1; 3,1; 2,3];
    modNames = {'AngNeu','HapNeu','AngHap'};

    % Through the models.
    parfor iE = 1:3

        %Just the 1s clips of stimuli
        cells1 = cellVals{modVec(iE,1)}{1};
        cells2 = cellVals{modVec(iE,2)}{1};

        % TWO TAILED! - P value adjusted accordingly as it is not in the script
        [clusters, p, t, pDist ] = ...
            permutest(cells1, cells2,false,0.025,nP,true);


        clustMap = zeros(dims);
        clustMap(cat(1,clusters{:}))=1;

        clustersHold{iE} = clusters;
        pHold{iE} = p;
        tHold{iE} = t;
        pDistHold{iE} = pDist;
        meanFirstHold{iE} = nanmean(cells1,3);
        meanSecondHold{iE} = nanmean(cells2,3);
        clustMapHold{iE} = clustMap;



    end


    for iE = 1:3
        resultsStim.(modNames{iE}).clusters = clustersHold{iE};
        resresultsStimults.(modNames{iE}).p = pHold{iE};
        resultsStim.(modNames{iE}).t = tHold{iE};
        resultsStim.(modNames{iE}).pDist = pDistHold{iE};
        resultsStim.(modNames{iE}).meanFirst = meanFirstHold{iE};
        resultsStim.(modNames{iE}).meanSecond = meanSecondHold{iE};
        resultsStim.(modNames{iE}).clustMap = clustMapHold{iE};
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









