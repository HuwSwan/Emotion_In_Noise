%% Final Analysis for 3rd Paper
% 1) MPS ratio of speech in noise
% 2) Spectral flux either as:
%      a) Revealed segments from OIM contain high flux moments of speech
%           i)measure flux in windows and in a random sleection of obscured
%           points....
%      b) Revealed segments cause large spectral flux that deviates from
%      the time average spectral flux.

% NO SPECTRAL FLUX! ITS GUFF

% clear all; close all; clc;

pTop = 'X:\PhD\03-Original_OIM';
pDat = fullfile(pTop,'02-analysis/2024-verification-dwgp_extracts');
pOut = fullfile(pTop,'06-analysis/2024-MPS_First');

types = {'normal','scram'};%

delete(gcp('nocreate'))
parpool(3);
for iT = 1:length(types)

    prDat = fullfile(pDat,types{iT});

    subj = dir(fullfile(prDat,'*detection*')); subj = {subj.name};
    prOut = fullfile(pOut,types{iT});
    mkdir(prOut)

    fs = 44100;
    Roughness = {};
    TF = {};
    MS = {};
    emos = [];

    k=1;

    mps = [];

    for iS = 1:length(subj)

        disp(iS)
        pSubjOut = fullfile(prOut,['subj_',num2str(iS)]);
        mkdir(pSubjOut)

        load(fullfile(prDat,subj{iS})); %loads wavTable.

        trIdx = ~isnan(cell2mat(wavTable.SNR));
        em = wavTable.cat(trIdx);

        clips{1} =  wavTable.slope_clip(trIdx);
        clips{2} = wavTable.noise_clip(trIdx);
        clips{3} = wavTable.combined(trIdx);
        emos = [emos;em];

        clear wavTable
        
        % Beautiful efficiency and vectorisation. Just... Beautiful.
        parfor p = 1:3
            subjMps(:,p) = cellfun(@MPS_analysis_contrastimproved,clips{p},repmat({fs},length(clips{p}),1),'UniformOutput',false);
        end


       %% Test        
        %Saves emo and MMPS outputs
        save(fullfile(pSubjOut,['subj_',num2str(iS),'_mps.mat']),'subjMps','em');      
        mps = [mps;subjMps];
        clear subjMps
    end

    % Cat output for whole study cell array.

    save(fullfile(prOut,'mps_data.mat'),'mps','emos');
    clear mps
end



%
%
%
%             for iW =1 :length(wavTable)
%
%
%                 for iW = 1:length(wavTable.combined)
%
%                 if trIdx(iW) == 0 %skip nans - controls and missed rials.
%                     continue
%                 end
%
%                 disp(iW)
%
%                 y = wavTable.slope_clip{iW};
%                 yN = wavTable.noise_clip{iW};
%                 yC = wavTable.combined{iW};
%
%                 [MS{k,1}, TF{k,1}, Roughness{k,1}] = MPS_analysis_contrastimproved(y,fs); %n
%                 [MS{k,2}, TF{k,2}, Roughness{k,2}] = MPS_analysis_contrastimproved(yN,fs);
%                 [MS{k,3}, TF{k,3}, Roughness{k,3}] = MPS_analysis_contrastimproved(yC,fs);
%
%                 k = k+1;
%             end
%         end
%     end
%     mps.MS = MS;
%     mps.TF = TF;
%     mps.Roughness = Roughness;
%     mps.emos = emos;
%
%
%     save(fullfile(prOut,'mps_data.mat'),'mps');
%     clear mps
%
% end

%
%
%
%
%
%             [yMS, yTF, Roughness(:,:,1)] = cellfun(@MPS_analysis_contrastimproved,wavTable.combined(trIdx),repmat({fs},sum(trIdx),1))
%
%             %THis is a disgustingly over vectorised bit of code.
%             [yMS, yTF, Roughness] = cellfun(@MPS_analysis_contrastimproved,wavTable.combined(trIdx),repmat({fs},sum(trIdx),1),'UniformOutput',false);
%
%
%
%                         MPS_analysis_contrastimproved(y,fs);
%
%
%
%
%             y = wavTable.slope_clip{iW};
%             yN = wavTable.noise_clip{iW};
%             yC = wavTable.combined{iW};
%
%
%             %Deos not take TFlog. That gives nan. Is this important? Who knwos.
%             clear Roughness
%             [yMS, yTF, Roughness(:,:,1)] = MPS_analysis_contrastimproved(y,fs);
%             [nMS, nTF, Roughness(:,:,2)] = MPS_analysis_contrastimproved(yN,fs);
%             [cMS, cTF, Roughness(:,:,3)] = MPS_analysis_contrastimproved(yC,fs);
%
%
%         end
%
%
%         % % This infromation can become mean MPS per Cycle per Modulation
%         % % Rate
%         % tiledlayout(1,3)
%         % nexttile
%         % imagesc(Roughness(:,:,1)); axis xy
%         % nexttile
%         % imagesc(Roughness(:,:,2)); axis xy
%         % nexttile
%         % imagesc(Roughness(:,:,3)); axis xy

%
%
%         difRough = Roughness(3,:)-Roughness(2,:)
%
%         % MPS Spectra
%
%         sig_noise = MS.val-nMS.val;
%         comb_sig = cMS.val-MS.val;
%         comb_noise = cMS.val-nMS.val;
%
%
%         B = max(cMS.val(:))-(2*(std(cMS.val(:))));
%         C = min(cMS.val(:))+(2*(std(cMS.val(:))));
%
%
%
%         %Differences
%         close all
%         tiledlayout(2,3)
%         nexttile
%         imagesc(MS.x,MS.y,sig_noise,[C B]); axis xy
%         nexttile
%         imagesc(MS.x,MS.y,comb_sig,[C B]); axis xy
%         nexttile
%         imagesc(MS.x,MS.y,comb_noise,[C B]); axis xy
%
%
%         %originals
%         % % close all
%         % figure()
%         % tiledlayout(1,3)
%         nexttile
%         imagesc(MS.x,MS.y,MS.val,[C B]); axis xy
%         nexttile
%         imagesc(MS.x,MS.y,nMS.val,[C B]); axis xy
%         nexttile
%         imagesc(MS.x,MS.y,cMS.val,[C B]); axis xy
%
%
%
%
%     end
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
% end
