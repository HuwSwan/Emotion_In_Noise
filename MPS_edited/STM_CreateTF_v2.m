function [TFout] = STM_CreateTF_v2(signal,Fs,maxfq,CB_filter,Args)
%
%  Caluclate the time frequency representation (spectrogram)
%
%   TFout = STM_CreateTF(signal,Fs <,CB_filter,Args>)
%           signal          The input auditory waveform
%           Fs                 Sampling rate
%           CB_filter      Type of filtering ('gauss','FIR', 'Shamma')
%           Args             Arguments, defaults given by Args = STM_CreateTF()
%
%           TFout          output TF structure
%               .TF           Downsampled Time Frequency matrix
%               .TFfull      Time frequency matrix in Fs sampling
%               .TFlog      Time frequency matrix in log units (DB)
%               .x_axis     x axis labels for plotting
%               .y_axis     y axis labels for plotting
%               .Args        Arguments used to create the TF structure
%
%   Examples:
%       TF = STM_CreateTF(signal,16000); % calculate the TF matrix
%       imagesc(TF.TF); axis xy;
%   
%       Args = STM_CreateTF; Args.do_plot = 1;
%       TF = STM_CreateTF(signal,16000,'gauss',Args); % calculate the TF matrix and auto-plot
%       
%
%       Adeen Flinker, Jan 2013 (adeen.f@gmail.com)
%
if ~exist('Args') % default Args
        Args.Nchans = 128;
        Args.CB_Filter = 'gauss';                      % FIR design in the tempotal domain or Full Width Half Max guassian in the frequency domain
        Args.CB_FIR_order = 256;                    % Filter order if using FIR filtering
        CF = 440 * 2 .^ ((-31:97)/24);           % Use Shamma's CF 24 channels/octave
        Args.CB_CenterFrs = round(CF(1:128)+diff(CF)/2);
                                                                              % Center Frequencies, computed here to support manual entry outside the function
        BW = 24.7*(Args.CB_CenterFrs*4.37/1000+1);
        Args.CB_FrHigh      = round(Args.CB_CenterFrs+BW/2);   
                                                                              % High Cutoff Frequencies, computed here to support manual entry outside the function
        Args.CB_FrLow       = round(Args.CB_CenterFrs-BW/2);
                                                                              % Low Cutoff Frequencies, computed here to support manual entry outside the function
%         Args.TF_ReFs          = 125;                    % Resampled rate for the spectrogram time axis
        Args.TF_ReFs          = 2*maxfq;                    % Resampled rate for the spectrogram time axis

        %%HUW EDIT - REMOVE NOISE FLOOR COS I WILL MEASURE NOISE - Or
        %%perhaps.... can this be the analysis?
        % Args.DBNoise          = 63;                      % Log noise threshold

        Args.DBNoise          = 0;                      % Log noise threshold

        Args.do_plot           = 0;                         % plotting flag
        
        %% HUW Edit for Stim creation.
%         gtFilt = gammatoneFilterBank([50,maxfq],128);
%         Args.CB_CenterFrs = round(getCenterFrequencies(gtFilt));                                                                           % Center Frequencies, computed here to support manual entry outside the function
%         BW = 24.7*(Args.CB_CenterFrs*4.37/1000+1);
%         Args.CB_FrHigh      = round(Args.CB_CenterFrs+BW/2);   
%                                                                               % High Cutoff Frequencies, computed here to support manual entry outside the function
%         Args.CB_FrLow       = round(Args.CB_CenterFrs-BW/2);
%         
        
        
end
if exist('CB_filter')  & ~isempty('CB_filter')
    Args.CB_Filter = CB_filter;
end
if nargin<1, TFout = Args; return; end
if size(signal,1)>1, signal = signal'; end        % transpose if column
if size(signal,1)>1, signal = signal(1,:); end % mono if stereo
% if Fs ~= 16000 & Fs ~= 8000
%     warning('Automatic Center Frequencies of 24 channels/octave are supported for 16 and 8 Khz Fs only, please enter them manually via Args,CB_CenterFrs');
% end
        
SignalLen = length(signal);
[p,q] = rat(Args.TF_ReFs /Fs);

for i = 1:Args.Nchans
        switch Args.CB_Filter
            case 'gauss'
                [SigFiltered(i,:) ,Pl(i)]= band_pass(signal,Fs,Args.CB_FrLow(i),Args.CB_FrHigh(i),1,'HMFWgauss');
                SigFilteredAmp(i,:) = (hilbert(SigFiltered(i,:)));
                
                %HUW EDIT: remove downsampling.
                
                
%                 TF(i,:) = abs(SigFilteredAmp(i,:));
                TF(i,:) = resample(abs(SigFilteredAmp(i,:)),p,q);

                
                
                
                
            case 'FIR'
                SigFiltered(i,:) = CB_FIR_filter(signal, Fs, SignalLen, Args.CB_FIR_order,Args.CB_FrLow(i),Args.CB_FrHigh(i));
                SigFilteredAmp(i,:)  =(hilbert(SigFiltered(i,:)));
                TF(i,:) = resample(abs(SigFilteredAmp(i,:)),p,q);
            case 'Shamma'
                if ~isfield(Args,'CB_paras');
                    if Fs == 16000
                    Args.CB_paras = [8 8 -2 0];
                    elseif Fs== 8000
                        Args.CB_paras = [8 8 -2 -1]
                    end
                end
                if Fs ~= 8000 & Fs ~=16000, error('unsupported Fs for nsltools'); end
                try 
                    loadload
                    TF = sqrt(wav2aud(signal,Args.CB_paras))';
                catch ME
                    error('nsltools functions not found');
                end
                Args.TF_ReFs = 1000/Args.CB_paras(1);
                CF = 440 * 2 .^ ((-31:97)/24); 
                Args.CB_CenterFrs = round(CF(1:128)+diff(CF)/2);
                SigFilteredAmp = TF;
                break
            otherwise, fprintf('Unrecognized CB_filter argument\n');
        end
end
TFout.TF = TF+abs(min(min(TF)));
TFout.TFfull =  SigFilteredAmp;
% TFout.TFlog = 20*log10(TFout.TF)+Args.DBNoise;
TFout.TFlog = 20*log(TFout.TF)+Args.DBNoise; %This is the closest one so far

% TFout.TFlog(TFout.TFlog<0) = 0; %why is this here?

TFout.x_axis = 1/Args.TF_ReFs:1/Args.TF_ReFs:size(TF,2)/Args.TF_ReFs;
TFout.y_axis = Args.CB_CenterFrs;
TFout.Args = Args;
TFout.Args.Fs = Fs; % save Fs for upsampling in Spectrum Inversion

if Args.do_plot 
    Chans_disp = 6;
    jm =  fix(Args.Nchans/Chans_disp);
    Ytick = jm:jm:Args.Nchans;
    YtickStr = {}; for l = 1:length(Ytick), YtickStr = {YtickStr{:} [num2str(TFout.Args.CB_CenterFrs(Ytick(l)))]}; end
    set(gca,'Ytick',[20:20:120],'YtickLabel',YtickStr);
    subplot(121);
    pcolor(TFout.x_axis,TFout.y_axis,TFout.TFlog); shading interp; title('linear scale with interpolation');
    subplot(122);
    imagesc(TFout.x_axis,1:Args.Nchans,TFout.TFlog);set(gca,'Ytick',[20:20:120],'YtickLabel',YtickStr); title('log scale'); axis xy;
end
end

function output = CB_FIR_filter(input, Fs, len, order, FrLow, FrHigh)
    s0=zeros(1,len+order/2); 
    s0(1:len)=input(1:len); 
    FIRfilter = fir1(order,[FrLow FrHigh]./Fs.*2);
    s1=filter(FIRfilter,1,s0);
    output(1:len) = s1(order/2+1:order/2+len);
end
    