
% FUNCTION called by Matched_Stimuli_Experiment_Interface
%   extracts filenames for audio stimuli
%   Written by Caitlyn Trevor
%   17 September 2019


function [filenames] = extract_filenames(extract_folder_name, nTrials)

            initial_directory = dir(extract_folder_name);
            initial_list_of_file_names = {initial_directory.name};
            filenames = initial_list_of_file_names(1,3:(nTrials+2));