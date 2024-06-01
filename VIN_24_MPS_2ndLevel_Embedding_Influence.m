%% The goal here is to look at how different the embedded VIN is from just the noise, and how that differs from scrambled to normal.
ccc
pTop = 'X:\PhD\03-Original_OIM'
pDat = fullfile(pTop,'02-analysis/2024-MPS_2nd_Level_CorrectedSNR');
pOut = fullfile(pTop,'02-analysis/2024-MPS_2nd_Level_CorrectedSNR/norm_scram');

pIms = fullfile(pTop,'02-analysis/graphics');

%load data

norm = load(fullfile(pDat,'normal','models','individualMPS.mat'));
scram = load(fullfile(pDat,'scram','models','individualMPS.mat'));


%This feeds from the cell values in _undoLogTransform so it's accurate
%% Sig testing

% THese are all two tailed so p=0.025 to account.

fnames ={'norm','scram'};
mStr = {'embed_difference'};

nP= 10^4;

for iM = 1%1:2
    clear meanRoughRats ratCollatte angClusters hapClusters angHapClusters hapAngClusters


    pOutDat = fullfile(pOut,'models');
    pOutPlot = fullfile(pOut,'plots');


    modNames = {'Neu','Ang','Hap'};

    for iT = 1:2
        % Through the models for norm vs scram differences.
        parfor iE = 1:3

            if iT ==1
                cells1 = norm.cellVals{iE}{3};
                cells2 = norm.cellVals{iE}{2};

            else
                cells1 = scram.cellVals{iE}{3};
                cells2 = scram.cellVals{iE}{2};
            end

            [clusters, p, t, pDist ] = ...
                permutest(cells1, cells2,false,0.025,nP,true);

            clustersHold{iE} = clusters;
            pHold{iE} = p;
            tHold{iE} = t;
            pDistHold{iE} = pDist;
            meanFirstHold{iE} = nanmean(cells1,3);
            meanSecondHold{iE} = nanmean(cells2,3);
        end

        for iE = 1:3
            results.(fnames{iT}).(modNames{iE}).clusters = clustersHold{iE};
            results.(fnames{iT}).(modNames{iE}).p = pHold{iE};
            results.(fnames{iT}).(modNames{iE}).t = tHold{iE};
            results.(fnames{iT}).(modNames{iE}).pDist = pDistHold{iE};
            results.(fnames{iT}).(modNames{iE}).meanFirst = meanFirstHold{iE};
            results.(fnames{iT}).(modNames{iE}).meanSecond = meanSecondHold{iE};
        end

    end
    mkdir(pOutDat)
    % Save the results this time huw.
    resStr = fullfile(pOutDat,['results_',mStr{iM},'.mat']);
    save(resStr,'results');
end







