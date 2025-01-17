function [outStruct] = MPS_analysis_contrastimproved(wavefile,fs)
%%

maxfq = 200; % Temporal Modulation max frequency default 200
% [signal,fs]=audioread(wavefile);

signal = wavefile;

% resample signal at 16000 Hz
fs2 = 16000;
if fs~=16000
    [p,q] = rat(fs/fs2);
    fs = fs2;
    signal = resample(signal,q,p);
else
end

% calculate cochleogram
TF = STM_CreateTF_v2(signal',fs,maxfq,'FIR');

% calculate MPS
MS0 = STM_Filter_Mod(TF);
Args = STM_Filter_Mod;
Args.MS_log = 0;

MS0 = STM_Filter_Mod(TF,[],[],Args);
% This got log transformd... for no reason. Christ. Undo that
MS2.orig_MS=log(MS0.orig_MS(fix(length(MS0.y_axis)/2+1):end,:));
MS2.x_axis=MS0.x_axis;
MS2.y_axis=MS0.y_axis(fix(length(MS0.y_axis)/2+1):end);

% resize output matrix to compare across sounds
% Removing the log transform because that is very wrong. Why was it ever
% included? A lesson in why you comment your code. M
% % MS.val = imresize(MS2.orig_MS, [64 400]);
MS.val = imresize(MS0.orig_MS, [64 400]);



[p,q] = rat(length(MS2.x_axis)/400);
MS.x = resample(MS2.x_axis,q,p);
[p,q] = rat(length(MS2.y_axis)/64);
MS.y = resample(MS2.y_axis,q,p);

% % extract values in the (30?150Hz) roughness range
% xs = [-150 -30 30 150];% roughness Freq limits in Hz (both <0 and >0 values are taken into account)
% for u=1:4; xz(u) = find(MS.x>xs(u),1,'first'); end
% roughness = squeeze(mean(mean(MS.val(:,[xz(1):xz(2),xz(3):xz(4)]),2),1));


%% Make system for finding roughness band values.
xsModMat = [-20 0 0 20;
    -40 -20 20 40;
    -60 -40 40 60;
    -80 -60 60 80;
    -100 -80 80 100;
    -120 -100 100 120;
    -140 -120 120 140;
    -150 -30 30 150];% roughness Freq limits in Hz (both <0 and >0 values are taken into account)

sp1 = min(MS.y):(48/64):max(MS.y);
sp2 = sp1(2:end);

xsSpecMat = [sp1(1:end-1)',sp2'];


k=1;
for iS = 1:length(xsSpecMat)

    xsS = xsSpecMat(iS,:);
    for u=1:2; xzS(u) = find(MS.y>xsS(u),1,'first'); end


    for iM = 1:length(xsModMat)
        xsM = xsModMat(iM,:);
        for u2=1:4; xzM(u2) = find(MS.x>xsM(u2),1,'first'); end
        roughness(iS,iM) = squeeze(mean(mean(MS.val([xzS(1):xzS(2)],[xzM(1):xzM(2),xzM(3):xzM(4)]),2),1));
    end
end

% Put all in one out struct so I can vectorise.
outStruct.MS = MS;
% outStruct.TF = TF;
outStruct.Roughness = roughness;

% 
% %% plot figure
% figure;
% subplot(2,2,1)
% plot(1/fs:1/fs:length(signal)/fs,signal)
% xlabel('time'); ylabel('Amplitude'); 
% 
% subplot(2,2,3)
% ylst = [0,1000,5000];
% ilst = []; for i = 1:length(ylst);  ilst(i) = find(TF.y_axis > ylst(i),1); end
% 
% imagesc(TF.x_axis,1:length(TF.y_axis),TF.TFlog); axis xy
% set(gca,'YTick',ilst,'YTickLabel',arrayfun(@(x)num2str(x/1000),ylst,'UniformOutput',false))
% xlabel('time'); ylabel('frequency (kHz)');
% 
% subplot(2,2,[2,4])
% B = max(MS.val(:))-(2*(std(MS.val(:))));
% C = min(MS.val(:))+(2*(std(MS.val(:))));
% imagesc(MS.x,MS.y,MS.val,[C B]); axis xy
% xlabel('Temporal Mod. (Hz)'); ylabel('Spectral Mod. (cycle./octave)');
% title('Modulation Power Spectrum')
% % 


