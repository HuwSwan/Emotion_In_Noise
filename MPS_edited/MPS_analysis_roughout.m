function [roughness] = MPS_analysis_roughout(wavefile)
%%

maxfq = 200; % Temporal Modulation max frequency 
[signal,fs]=audioread(wavefile);

% resample signal at 16000 Hz
fs2 = 16000;
if fs~=16000
    [p,q] = rat(fs/fs2);
    fs = fs2;
    signal = resample(signal,q,p);
else
end

% calculate cochleogram
TF = STM_CreateTF_v2(signal',fs,maxfq,'gauss');

% calculate MPS
MS0 = STM_Filter_Mod(TF);
Args = STM_Filter_Mod;
Args.MS_log = 0;

MS0 = STM_Filter_Mod(TF,[],[],Args);
MS2.orig_MS=log(MS0.orig_MS(fix(length(MS0.y_axis)/2+1):end,:));
        % MS2.orig_MS=MS0.orig_MS(fix(length(MS0.y_axis)/2+1):end,:);
        % MS2.orig_MS=MS0.orig_MS;
MS2.x_axis=MS0.x_axis;
MS2.y_axis=MS0.y_axis(fix(length(MS0.y_axis)/2+1):end);
        % MS2.y_axis=MS0.y_axis(fix(length(MS0.y_axis)/2+1):end);
        % MS2.y_axis=MS0.y_axis;

% resize output matrix to compare across sounds
MS.val = imresize(MS2.orig_MS, [64 400]);
[p,q] = rat(length(MS2.x_axis)/400);
MS.x = resample(MS2.x_axis,q,p);
[p,q] = rat(length(MS2.y_axis)/64);
MS.y = resample(MS2.y_axis,q,p);

% extract values in the (30?150Hz) roughness range
xs = [-150 -30 30 150];% roughness Freq limits in Hz (both <0 and >0 values are taken into account)
for u=1:4; xz(u) = find(MS.x>xs(u),1,'first'); end
roughness = squeeze(mean(mean(MS.val(:,[xz(1):xz(2),xz(3):xz(4)]),2),1));


% plot figure
figure;
subplot(2,2,1)
plot(1/fs:1/fs:length(signal)/fs,signal)
xlabel('time'); ylabel('Amplitude'); 

subplot(2,2,3)
ylst = [0,1000,5000];
ilst = []; for i = 1:length(ylst);  ilst(i) = find(TF.y_axis > ylst(i),1); end

imagesc(TF.x_axis,1:length(TF.y_axis),TF.TFlog); axis xy
set(gca,'YTick',ilst,'YTickLabel',arrayfun(@(x)num2str(x/1000),ylst,'UniformOutput',false))
xlabel('time'); ylabel('frequency (kHz)');

subplot(2,2,[2,4])
B = max(MS.val(:))-(2*(std(MS.val(:))));
C = min(MS.val(:))+(2*(std(MS.val(:))));
imagesc(MS.x,MS.y,MS.val,[C B]); axis xy
xlabel('Temporal Mod. (Hz)'); ylabel('Spectral Mod. (cycle./octave)');
title('Modulation Power Spectrum')



