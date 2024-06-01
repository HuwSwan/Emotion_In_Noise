%% The goal here is to look at how different the embedded VIN is from just the noise, and how that differs from scrambled to normal.
ccc
pTop = 'X:\PhD\03-Original_OIM'
pDat = fullfile(pTop,'02-analysis/2024-MPS_2nd_Level_CorrectedSNR');
pOut = fullfile(pTop,'02-analysis/2024-MPS_2nd_Level_CorrectedSNR/VIN_EmotionContrast');

pIms = fullfile(pTop,'02-analysis/graphics');

%load data

norm = load(fullfile(pDat,'normal','models','individualMPS.mat'));
scram = load(fullfile(pDat,'scram','models','individualMPS.mat'));

% for rearranging stats and plots etc.
axesValues = load('X:\PhD\03-Original_OIM\02-analysis\2024-MPS_2nd_Level_CorrectedSNR\normal\models\axes_values.mat');
[absX, xIdx] = sort(abs(axesValues.xAx));


%This feeds from the cell values in _undoLogTransform so it's accurate
%% Sig testing

% THese are all two tailed so p=0.025 to account.

%% TESTING NORMAL VS SCRAMBLED! INDIVIDUAL AND GROUPING EMOTIONAL SPEECH

fnames ={'normal','scramble'};
mStr = {'embed_difference'};

nP= 10^4;

for iF = 1:2
    clear meanRoughRats ratCollatte angClusters hapClusters angHapClusters hapAngClusters


    pOutDat = fullfile(pOut,'models');
    pOutPlot = fullfile(pOut,'plots');


    modNames = {'AngNeu','HapNeu'};

    mIdx = [2,1;3,1];


    % Through the models for norm vs scram differences.


    parfor iE = 1:2

        if iF ==1 %normal
            cells1 = norm.cellVals{iE+1}{3}; %affects
            cells2 = norm.cellVals{1}{3}; %neutral
        else
            cells1 = scram.cellVals{iE+1}{3}; %normin noise
            cells2 = scram.cellVals{1}{3}; %scram in noise
        end

        


            [clusters, p, t, pDist ] = ...
                permutest(cells1, cells2,false,0.025,nP,true);

            clustersHold{iE} = clusters;
            pHold{iE} = p;
            tHold{iE} = t;
            pDistHold{iE} = pDist;
            meanFirstHold{iE} = nanmean(cells1,3);
            meanSecondHold{iE} = nanmean(cells2,3);

        for iE = 1:3
            results.(fnames{iF}).(modNames{iE}).clusters = clustersHold{iE};
            results.(fnames{iF}).(modNames{iE}).p = pHold{iE};
            results.(fnames{iF}).(modNames{iE}).t = tHold{iE};
            results.(fnames{iF}).(modNames{iE}).pDist = pDistHold{iE};
            results.(fnames{iF}).(modNames{iE}).meanFirst = meanFirstHold{iE};
            results.(fnames{iF}).(modNames{iE}).meanSecond = meanSecondHold{iE};
        end

    end
    mkdir(pOutDat)
    % Save the results this time huw.
    resStr = fullfile(pOutDat,['results_oneSideMPS',fnames{iF},mStr{iM},'.mat']);
    save(resStr,'results');
end







