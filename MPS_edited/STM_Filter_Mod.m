function [MSout] = STM_Filter_Mod(TFstruct,TemModFilter, SpecModFilter,Args)
%
%  Caluclate modulation spectrum and filter
% 
%   MSout = STM_Filter_Mod(TFstruct, <TemModFilter, SpecModFilter,Args>)
%           TFstruct                Time frequency representation from STM_CreateTF
%           TemModFilter      Temporal Modulation filter cutoff
%           SpecModFilter      Spectral Modulation filter cutoff
%           Args                       Arguments, defaults given by Args = STM_Filter_Mod
%
%           MSout                       output MS structure
%               .new_TF                New TF after MS filtering
%               .new_MS               New MS amplitude after filtering
%               .new_MS_phase   New MS phase after filtering
%               .orig_TF                Original TF 
%               .orig_MS               Original MS amplitude (no filtering)
%               .orig_MS_phase   Original MS phase (no filtering)
%               .x_axis                   MS x axis labels for plotting
%               .y_axis                   MS y axis labels for plotting
%               .Args                     Arguments used to create the MS structure
%               .Args_TF               Arguments used to create the TF  structure (passed from STM_CreateTF)
%
%   Examples:
%       TF = STM_CreateTF(signal,16000); % calculate the TF matrix
%       MS = STM_Filter_Mod(TF);
%       imagesc(MS.x_axis,MS.y_axis,log(MS.orig_MS)); axis xy;
%       imagesc(MS.x_axis,MS.y_axis(fix(length(MS.y_axis)/2+1):end),log(MS.orig_MS(fix(length(MS.y_axis)/2+1):end,:))); axis xy;
%     
%       MS = STM_Filter_Mod(TF,10,15);
%       imagesc(MS.new_TF); axis xy;
%       
%
%       Adeen Flinker, Jan 2013 (adeen.f@gmail.com)
%

if ~exist('Args')
    Args.MS_log = 0;                 % work on log transformed spectrogram in DB (1) or analytic amplitude (0)
    Args.MS_filter = 'lowpass'; % lowpass, highpass or bandpass 
    Args.MS_fstep       = 1/24;        % Default 24 chans per ocatve (default in TFstruct.Args.CF_CenterFrs)
end

if nargin<1, MSout = Args; return; end
% 
% if Args.MS_log
    % TF = TFstruct.TFlog;
% else
    TF = TFstruct.TF;
% end

if ~exist('SpecModFilter') 
    SpecModFilter = Inf;
elseif isempty(SpecModFilter)
    SpecModFilter = Inf;
end

if ~exist('TemModFilter') 
    TemModFilter = Inf;
elseif isempty(TemModFilter)
    TemModFilter = Inf;
end

%THIS IS THE MODULATION SPECTRUM CALCULATION> WHY ISNT IT COMMENTED
orig_ms = fft2(TF);

[Nfr,Ntm] = size(TF);
fstep = Args.MS_fstep;
ampsrate = TFstruct.Args.TF_ReFs;
for i=1:ceil((Nfr+1)/2)
    dwf(i)= (i-1)*(1/(fstep*Nfr));    % positive spectral modulation frequencies
    if (i > 1)
        dwf(Nfr-i+2)=-dwf(i);          % negative spectral modulation frequencies
    end
end
for i=1:ceil((Ntm+1)/2)
    dwt(i) = (i-1)*(ampsrate/Ntm);
    if (i > 1 )
        dwt(Ntm-i+2) = -dwt(i);
    end
end

switch  Args.MS_filter
    case 'lowpass'
        wf_high = SpecModFilter;
        wt_high = TemModFilter;
        % code below modified from Eliott et al. modfilter.m
        % This defines the ramp of the gain from 0 to 1
        
        dfi=0.0001;   % 0.0001 cycle per octave ramp in frequency
        dti=1;        % One Hz ramp in time
        
        newamp = abs(orig_ms);
        newphase = angle(orig_ms);
        %the use of gain is misleading, it is not gain.
        % But this is where we manipulate the signal in "roughness
        % dimension".
        % Multiply all this stuff.
        
        % For some frequency bands, include gainmap for middle and high manipualtion
        % Make sure to remove the filtering and phase erasure.
        
        gainmap=ones(size(orig_ms)); % Define a gain by which to multiply the mod spectrum
        for f=1:Nfr
         for t=1:Ntm
            % Define the regions to set to zero gain - first along the
            % wf axis
            if ((abs(dwf(f)))>wf_high+dfi)
                gainmap(f, t) = 0.0;
                newphase(f,t) = (rand(1)-0.5)*2*pi;     % Randomize the phase
                newphase(1,1) = 0.0;  % The phase of the DC value has to be zero
            end

            if (wf_high~=0)
                if ((abs(dwf(f))>=wf_high) && (abs(dwf(f))<=(wf_high+dfi)))
                    gainmap(f,t)=gainmap(f,t)*cos((((abs(dwf(f))-wf_high)./dfi))*(pi/2))^2;
                end
            end

            % Define the regions to set to zero - along the wt axis
            if ((abs(dwt(t)))>wt_high+dti)
                gainmap(f, t) = 0.0;
                newphase(f,t) = (rand(1)-0.5)*2*pi;     % Randomize the phase
            end

            if (wt_high~=0)
                if ((abs(dwt(t))>=(wt_high) && (abs(dwt(t))<=(wt_high+dti))))
                    gainmap(f,t)=gainmap(f,t)*cos((((abs(dwt(t))-wt_high)./dti))*(pi/2))^2;
                end
            end
        end
    end
    newphase(1,1) = 0.0;  % The phase of the DC value has to be zero
end

newamp = newamp.*gainmap;

new_ms = newamp.*exp(complex(0,newphase));

MSout.new_TF = real(ifft2(new_ms));
MSout.new_MS = fftshift(abs(new_ms));
MSout.new_MS_phase =fftshift(angle(new_ms));
MSout.orig_TF = TF;
MSout.orig_MS = fftshift(abs(orig_ms)); % This shows AM power at Y=centre frequencies and X = AM rate.


MSout.orig_MS_phase = fftshift(angle(orig_ms));
MSout.x_axis = fftshift(dwt);
MSout.y_axis = fftshift(dwf);
MSout.Args = Args;
MSout.TF_Args = TFstruct.Args;

end

