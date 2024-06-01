%% Rerunning analysis for simplicity.

%Extract 1 second prior to detection. (maybe do -2,-1, and +1 for time
%saving)
clear all; close all; clc;

pTop = 'X:\PhD\02-OIM_VIN'
pDataT = fullfile(pTop,'01-data');
pRT = fullfile(pDataT,'behavioural_RT');
pOutT = fullfile(pTop,'analysis/dwgp_extracts');
mkdir(pOutT);
pStim = fullfile(pTop,'04-stimuli');

types = {'normal','scram'};%
types = {'normal'};%

%exclude subj 8 (nv09) and 18 (nv110)

pLevel = 70;

%% OIM Extraction

for iT = 1:length(types)
    disp(types{iT})
    prDat = fullfile(pDataT,types{iT});
    prOut = fullfile(pOutT,types{iT});
    subj = dir(fullfile(prDat,'nv*')); subj = {subj.name};
    %excluding nv9 and 110
    subj([8 18]) = [];

    mkdir(prOut);
    rtFiles = dir(fullfile(pRT,'nv*'));
    rtFiles = {rtFiles.name};

    %for comps
    stims = dir(fullfile(pStim,types{iT},'ramp*')); stims={stims.name};

    for iS = 1:length(stims)
        [y fs] = audioread(fullfile(pStim,types{iT},stims{iS}));
        %replicating presentation script details - normalise.
        stimuli{iS} = (y./rms(y));
    end

    noise = audioread(fullfile(pStim,types{iT},'pinknoise.wav'));
    %replicates presentation script for normalising.
    noise = noise./rms(noise);  % normalization
    noiseInt = 0.05; %taken from OG presentation scripts.
    noise = noise.*noiseInt;

    %RT compensation
    tR = 10;
    snrR = 13.5;
    bias = -26;
    steps = snrR/tR;

    for iS = 1:length(subj)

        pSubj = fullfile(prDat,subj{iS});
        subj{iS}
        runs = dir(fullfile(pSubj,'NV_fMRI*'));   runs = {runs.name};
        rtDat = load(fullfile(pRT,rtFiles{iS}));
        meanRT = nanmean(rtDat.Data(:,4));
        snrRT = meanRT*steps; %equivelant adjustment for RT as an SNR change. Subtract from runmat3.

        wavTable.part = {};
        wavTable.rt = [];
        wavTable.cat = [];
        wavTable.wav = {};
        wavTable.rtTrialStart = [];
        wavTable.rt_slope = [];
        wavTable.noiseRT = [];
        wavTable.stim_clip = []; %raw stim
        wavTable.slope_clip = []; %stim with snr slope
        wavTable.noise_clip = []; %noise
        wavTable.combined = []; %stim and noise combined
        wavTable.SNR = []; % The SNR of the overall 3s clips
        wavTable.SlDelta = [];
        wavTable.rt_SNR = [];
        wavTable.objscore_global = [];


        %Populate detectrion times for all tracks
        for iR = 1:length(runs)
            partCell={};
            load(fullfile(pSubj,runs{iR}));

            runMat(:,[1 2 3]) = table2array(Data(:,[4 17 19]));% added column 19, SNR result

            partCell(1:length(runMat)) = subj(iS);

            wavTable.part = cat(1,wavTable.part,partCell');
            wavTable.cat = cat(1,wavTable.cat,runMat(:,1));%categoriy emo
            wavTable.SlDelta = cat(1,wavTable.SlDelta,Data.sloDeltrial);
            wavTable.rt_SNR = cat(1,wavTable.rt_SNR,runMat(:,3));
            %RT in relation to start of slope onset.

            %% Account for RT offsets due to system lag - precision required.
            % measures events, not the schedules they were on.
            t0 = header.t0;
            stimOns = Data.stimOns2-t0;
            noiseOns = Data.noiseOns2-t0;
            trueRT = Data.tKey-stimOns; % already accounts for the slope problems.

            rtAdj = trueRT-meanRT; %account for average RT of participants.
            rtslope = rtAdj-Data.sloDeltrial; %account for slope start.
            wavTable.rt = cat(1,wavTable.rt,rtAdj);
            wavTable.rt_slope = cat(1,wavTable.rt_slope,rtslope);
            wavTable.wav = cat(1, wavTable.wav, Data.Stimulus);

            %time from start of voice track, not slope onset
            wavTable.rtTrialStart = cat(1,wavTable.rtTrialStart,rtAdj);

            %create a value that is RT in relation to start of noise.
            % RT values +
            rtNoise = Data.tKey-noiseOns;
            wavTable.noiseRT = cat(1,wavTable.noiseRT,rtNoise);


            %loop through for before and after analyses? Maybe no.


            stim_clip =[];
            noise_clip = [];
            for iW = 1:length(Data.Stimulus)

                wavStr = Data.Stimulus{iW};
                %the slope for shaping
                slopeShape = Data.slopeShape{iW};
                %response check
                k = isnan(wavStr); k=unique(k);

                %Any reactions that are implausible quick are filtered out.
                %I.e. within 2.5  second of starting the trial.
                %                 kRT = rtTrialStart(iW)>3;

                kRT = Data.RT(iW)>2.5;

                %Check to see that stimuli is not responded too late. No
                %later than 11s or else analysis breaks. - Use rtAj. This
                %is adjusted to stimuli, not RT.
                kLate = rtAdj(iW)<11;

                if k == 0 & kRT == 1 & kLate==1;
                    idx = find(contains(stims,wavStr));

                    %idx for sample of the stimuli

                    % Find index numbers based on the ajdusted reaction

                    sampIdxTr = [round(fs*rtAdj(iW))-(fs),round(fs*rtAdj(iW))];
                    sampIdxNoi = [round(fs*rtNoise(iW))-(fs),round(fs*rtNoise(iW))];

                    % load sound for extract and adjust to the maximum of
                    % rms slope.
                    holdSample = stimuli{idx};
                    stimSample = holdSample.*(max(slopeShape));

                    %raw sample at max of DB slope.
                    stim_clip{1,iW} = stimSample(sampIdxTr(1):sampIdxTr(2));
                    %shape by Data.slopeShape
                    shapeSample = holdSample.*slopeShape';
                    slope_clip{1,iW} = shapeSample(sampIdxTr(1):sampIdxTr(2));
                    %Noise extraction
                    noise_clip{1,iW} = noise(sampIdxNoi(1):sampIdxNoi(2));

                    %combine slope and noise.
                    combined{1,iW} = slope_clip{1,iW}+noise_clip{1,iW};

                    %Conduct the BI-DWGP here - full 3s sample..
                    ref = slope_clip{1,iW}; %needs to be the sloped one.
                    nseIn = noise_clip{1,iW};
                    sig = combined{1,iW};
                    %This is the true SNR of 1s prior to detection. Not the equivelant
                    %global SNR.
                    stimRMS = rms(ref);
                    noiseRMS = rms(nseIn);
                    SNR{iW} = 10*log(stimRMS/noiseRMS);

                    objscore_global{1,iW} = window_DWGP(sig, nseIn, fs, ref, pLevel);  %retunds DWGP weighting and glimpse portion
return

                else
                    stim_clip{iW} = NaN;
                    noise_clip{iW} = NaN;
                    combined{iW} = NaN;
                    slope_clip{iW} = NaN;
                    SNR{iW} = NaN;
                    objscore_global{iW} = NaN;
                end
            end

            wavTable.stim_clip = cat(1,wavTable.stim_clip,stim_clip');
            wavTable.noise_clip = cat(1,wavTable.noise_clip,noise_clip');
            wavTable.combined = cat(1,wavTable.combined,combined');
            wavTable.slope_clip = cat(1,wavTable.slope_clip,slope_clip');
            wavTable.SNR =  cat(1,wavTable.SNR,SNR');
            wavTable.objscore_global = cat(1,wavTable.objscore_global,objscore_global');

        end

        save(fullfile(prOut,['1s_Pre_detection',subj{iS},'.mat']),'wavTable');

    end
end


%% Collate all details into a single table for the DWGP results
% 
% clear all; close all;
% pTop = 'Z:\SWANBOROUGH_VinDec\Concept_Trials\OIM\OIM2022'
% pDatT = fullfile(pTop,'analysis/dwgp_extracts');
% pOutT = fullfile(pTop,'analysis/dwgp_collation');
% mkdir(pOutT);
% 
% types = {'normal'}%,'scram'};%
% 
% for iT = 1:length(types)
%     
%     pData = fullfile(pDatT,types{iT});
%     files = dir(fullfile(pData,'*.mat')); files = {files.name};
%     
%     RT = []; %Reaction time is already adjusted for trial onset
%     RT_slope = []; %RT adjusted to slope.
%     SlDelt = [];
%     local_weights =[];
%     global_weights = [];
%     local_gps = [];
%     global_gps = [];
%     part = {};
%     emos = [];
%     snr = [];
%     snr_rt = [];
%     mask = {};
%     rtvMask = {};
%     wavs = {};
%     
%     for iF=1:length(files)
%         clear snrHold localWeightHold globalWeightHold localGpHold globalGpHold emoHold mk rvMk rtHold slHold wavHold rtslHold snr_rtHold
%         
%         load(fullfile(pData,files{iF})); %loads wavTable
%         T=struct2table(wavTable);
%         k=1; %for removing blank rows.
%         disp(['subj ',num2str(iF)]);
%         
%         for iR = 1:length(wavTable.part)
%             
%             %excludingn catch trials, null responses, or sections not big
%             %enough for 3s OIM extraction.
%             if  isnan(wavTable.cat(iR)) | isnan(wavTable.rt_slope(iR)) |isnan(wavTable.stim_clip{iR})
%                 continue
%             end
%             
%             %collate details
%             globalWeightHold(k) = wavTable.objscore_global{iR}.DWGP;
%             
%             globalGpHold(k) = wavTable.objscore_global{iR}.GP;
%             
%             mk{k} = wavTable.objscore_global{iR}.mask;
%             rvMk{k} = wavTable.objscore_global{iR}.rtv_mask;
%             
%             emoHold(k) = wavTable.cat(iR);
%             snrHold(k) = wavTable.SNR(iR);
%             snr_rtHold(k) = wavTable.rt_SNR(iR);
%             rtHold(k) = wavTable.rtTrialStart(iR);
%             rtslHold(k) = wavTable.rt_slope(iR);
%             slHold(k) = wavTable.SlDelta(iR);
%             wavHold(k) = wavTable.wav(iR);
%             
%             k=k+1;
%         end
%         %cat everything for table assignment.
%         RT = [RT;rtHold'];
%         RT_slope = [RT_slope;rtslHold'];
%         SlDelt = [SlDelt;slHold'];
%         
%         global_weights = [global_weights;globalWeightHold'];
%         
%         global_gps = [global_gps;globalGpHold'];
%         
%         
%         snr = [snr;snrHold'];
%         snr_rt = [snr_rt;snr_rtHold'];
%         pCell = cell(length(rtHold),1); pCell(:) = wavTable.part(1);
%         part = [part;pCell];
%         emos = [emos;emoHold'];
%         mask = [mask;mk'];
%         rtvMask = [rtvMask;rvMk'];
%         wavs = [wavs;wavHold'];
%     end
%     
%     %Assign all to a single table for each trial.
%     dwgp_table = table(part,emos,global_weights,global_gps,snr,snr_rt,mask,...
%         rtvMask,RT,RT_slope,SlDelt,wavs);
%     
%     %save it all.
%     mkdir(fullfile(pOutT,types{iT}))
%     save(fullfile(pOutT,types{iT},'DWGP_table.mat'),'dwgp_table');
%     
% end

%% Extract all bubble information and export it for LCA.
clear all; close all;

pTop = 'Z:\SWANBOROUGH_VinDec\Concept_Trials\OIM\OIM2022'
pDatT = fullfile(pTop,'analysis/dwgp_collation');


types = {'normal'}%,'scram'};%

for iT = 1:length(types)
    
    load(fullfile(pDatT,types{iT},'dwgp_table.mat'));
    
    % Creating window size and coherence analysis.
    regionSize= {};
    timeWindow = {};
    hzWindow = {};
    timeCentre = {};
    hzCentre = {};
    windowN = [];
    windowN3K = [];
    windowPixels = {};
    
    for iD = 1:height(dwgp_table)
        
        maskHold = dwgp_table.mask{iD};
        %find coherent bubbles > 3 pixels.
        regs = regionprops(maskHold==1,'Area','PixelIdxList');
        %regs = regs([regs.Area]>=3); % no thresh currently.
        %threshold notes - 1 row = 1 ECB. 1 column = 10ms
        %Process for finding out the common size of bubbles for identification.
        clear regK tK hzk tkMed hzkMed
        for iR = 1:length(regs)
            % Take the indeces
            [row,col] = ind2sub([34,101],regs(iR).PixelIdxList);
            regK{iR} = regs(iR).Area;
            tk{iR} = [min(col),max(col)]; % save time bins
            hzk{iR} = [min(row),max(row)]; % save erb bins
            tkMed{iR} = median(tk{iR}); %median points of the circles.
            hzkMed{iR} = median(hzk{iR});
            pixMasks{iR} = regs(iR).PixelIdxList; %index of the pixels that are exposed.
        end
        
        winN = length(regs);
        winNK = length(regs([regs.Area]>=3));
        
        
        regionSize = [regionSize,{regK}];
        timeWindow = [timeWindow,{tk}];
        hzWindow = [hzWindow,{hzk}];
        timeCentre = [timeCentre,{tkMed}];
        hzCentre = [hzCentre,{hzkMed}];
        windowN = [windowN;winN];
        windowN3K = [windowN3K;winNK];
        windowPixels = [windowPixels; {pixMasks}];
    end
    
    dwgp_table.regionSize = regionSize';
    dwgp_table.timeWindow = timeWindow';
    dwgp_table.hzWindow = hzWindow';
    dwgp_table.timeCentre = timeCentre';
    dwgp_table.hzCentre = hzCentre';
    dwgp_table.windowN = windowN;
    dwgp_table.windowN3K = windowN3K;
    dwgp_table.windowPixels = windowPixels;
    
    dwgp_table.snr = cell2mat(dwgp_table.snr) ;
    
    %% Preperation of table for analysis on a "glimpse" basis. I.e. each cluster of visible frames.
    % Calculate AREA of glimpses and the properties the areas contain K>2
    % Freqency bands-size-weight.
    
    hzHeights = [];
    hzCentres = [];
    tLengths = [];
    windowWeights = [];
    windowSize = [];
    emo = [];
    RT = [];
    stimRef = [];  % to check co-occurance.
    time = [];
    
    for iD = 1:height(dwgp_table)
        
        %Thresholding windows.
        vecSize = cell2mat(dwgp_table.regionSize{iD});
        sizeIdx = vecSize>=3;
        
        twins = dwgp_table.timeWindow{iD}(sizeIdx);
        hzwins = dwgp_table.hzWindow{iD}(sizeIdx);
        frameIdx = dwgp_table.windowPixels{iD}(sizeIdx);
        
        %For summing power of visible frames.
        weights = cell2mat(cellfun(@(x) mean(dwgp_table.rtvMask{iD}(x)), frameIdx, 'UniformOutput',0));
        winSze = cell2mat(dwgp_table.regionSize{iD}(sizeIdx));
        thold = cell2mat(dwgp_table.timeCentre{iD}(sizeIdx));
        
        %% Work out the intercept of window length in time, and it's spectral centre for each emotion.
        % The important factor is how the length of window interacts with
        % HZ for emotional voice detection.
        
        tLengths = [tLengths;cellfun(@(x) x(2)-x(1)+1,twins)'];
        hzHeights = [hzHeights;cellfun(@(x) x(2)-x(1)+1,hzwins)'];
        windowWeights = [windowWeights;weights'];
        windowSize = [windowSize;winSze'];
        stimRef = [stimRef;repmat(iD,sum(sizeIdx),1)];
        time = [time; thold'];
        
        hold = cellfun(@(x) x(1) ,dwgp_table.hzWindow{iD}(sizeIdx)); %this returns just the lowest of the Bands.
        hzCentres = [hzCentres;hold'];
        
        %duplicate important info from DWGP for latre analysis ease.
        emo = [emo;repmat(dwgp_table.emos(iD),length(twins),1)];
        RT = [RT;repmat(dwgp_table.RT_slope(iD),length(twins),1)];
    end
    window_tab = table(stimRef,emo,time,RT,hzCentres,hzHeights,tLengths,windowSize,windowWeights);
    summary(window_tab)
    
    mkdir(fullfile(pDatT,types{iT},'window_features_speech'));
    save(fullfile(pDatT,types{iT},'window_features_speech','summary_windows.mat'),'window_tab');
    writetable(window_tab,fullfile(pDatT,types{iT},'window_features_speech','summary_windows.csv'));
    
end

% Next step - Analysis with LCA in R, plus general statics of window
% properties by class etc.





































