
% Script to find roughness value of multiple audio files sorted by category
% Written by Caitlyn Trevor 11/15/19

% Step 1: Initial set-up
% Step 2: Pull audiofile names 
% Step 2: Run MPS function on each category of stimuli with roughness
% filling out the results matrix


% VERSION LOG:

% V.1 = Initial code generated - entire code working

%% STEP 1: INITIAL SET-UP

    % Set WD
    cd C:\Users\Caitlyn\Drive\Research\Screaming_Strings_Silver_Screen_Paper\Scripts
    
    % Adding paths to stimuli
    addpath(genpath('\stim_for_analyses')); % adds folder plus subfolders (genpath) to path
            
    % Clear the workspace and the screen
    clear all; % clears all workspace variables
    close all; % closes all open figures    
    sca % executes "Screen(CloseAll)", also unhides cursor
%     clc % clears command window
    
    % Create necessary variables
    nfiles = 200;
    mean_roughs = zeros(nfiles, 1);
    category = zeros(nfiles, 1);
        % 1 = Screams
        % 2 = Non-screams
        % 3 = Scream-like music
        % 4 = Non-scream-like music
    
    
%% STEP 2: PULL AUDIO FILE NAMES
    
    % Pull audio filenames
    audio_filenames = extract_filenames('stim_for_analyses', 200);
    

%% STEP 3: FIND ALL ROUGHNESS VALUES FOR ALL AUDIO FILES IN DIFFERENT CATEGORIES

for f = 1:nfiles
    
    mean_roughs(f) = MPS_analysis_roughout(audio_filenames{f});
    
    if strcmp(audio_filenames{f}(9:15),'scr_fea')
        mean_roughs(f) = MPS_analysis_roughout(audio_filenames{f});
        category(f) = 1;
    elseif strcmp(audio_filenames{f}(9:15),'scr_neu')
        mean_roughs(f) = MPS_analysis_roughout(audio_filenames{f});
        category(f) = 2;
    elseif strcmp(audio_filenames{f}(5:15),'nonscreamus')
        mean_roughs(f) = MPS_analysis_roughout(audio_filenames{f});
        category(f) = 4;
    elseif strcmp(audio_filenames{f}(5:10),'scrmus')
        mean_roughs(f) = MPS_analysis_roughout(audio_filenames{f});
        category(f) = 3;
    end
end


%% STEP 4: Save Data 

 % Change working directory to Data folder
        cd Data
        
    % Save workspace as MATLAB file
        save(sprintf('mean_roughness_data_%s', date));
        
    % Save datafile matrix as excel spreadsheet
        T = table(category, mean_roughs);
        filename = sprintf('mean_roughness_data_%s', date);
        writetable(T,filename,'FileType','spreadsheet','WriteVariableNames', true);
    
    % Leave Data folder for original directory
        cd ..
    
%% Analyze data quickly

    % Find means by category
    varfun(@mean,T,'InputVariables','mean_roughs','GroupingVariables', 'category')
    