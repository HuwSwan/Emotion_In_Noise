%% The goal here is to look at how different the embedded VIN is from just the noise, and how that differs from scrambled to normal.
ccc
pTop = 'X:\PhD\03-Original_OIM'
pDat = fullfile(pTop,'02-analysis/2024-MPS_2nd_Level_CorrectedSNR');
pOut = fullfile(pTop,'02-analysis/2024-MPS_2nd_Level_CorrectedSNR/norm_scram_one_sided');

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

fnames ={'individual','affectNorm'};
mStr = {'embed_difference'};

nP= 10^4;

for iM = 1%1:2
    clear meanRoughRats ratCollatte angClusters hapClusters angHapClusters hapAngClusters


    pOutDat = fullfile(pOut,'models');
    pOutPlot = fullfile(pOut,'plots');


    modNames = {'Neu','Ang','Hap'};


            eIts = [3,2]; %for indexing par processing for grouping all emo together.


    for iT = 1:2
        % Through the models for norm vs scram differences.


        parfor iE = 1:eIts(iT)

            if iT ==1 %individual
                cells1 = norm.cellVals{iE}{3}; %normin noise
                cells2 = scram.cellVals{iE}{3}; %scram in noise

                %Create 1 sided version.
                cells1 = cells1(:,xIdx,:);
                cells2 = cells2(:,xIdx,:);

            else %combining emotions.

                if iE == 1
                    cells1 = norm.cellVals{iE}{3}; %normin noise
                    cells2 = scram.cellVals{iE}{3}; %scram in noise
                else %combine emotion
                    %ang
                    cells1A = norm.cellVals{iE}{3}; %normin noise
                    cells2A = scram.cellVals{iE}{3}; %scram in noise
                    %hap
                    cells1H = norm.cellVals{iE+1}{3}; %normin noise
                    cells2H = scram.cellVals{iE+1}{3}; %scram in noise

                    cells1 = cat(3,cells1A,cells1H); %all norm emotions
                    cells2 = cat(3,cells2A,cells2H); %all scram emotions.
                end
                %Create 1 sided version.
                cells1 = cells1(:,xIdx,:);
                cells2 = cells2(:,xIdx,:);
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
    resStr = fullfile(pOutDat,['results_oneSideMPS',mStr{iM},'.mat']);
    save(resStr,'results');
end







