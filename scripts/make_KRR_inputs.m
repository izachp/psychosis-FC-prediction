function make_KRR_inputs(repo, scale, timepoint, fc_type, preproc, parc, outcome_type, group)
%% Concatenate FC matrices into .mat and save subject list and change scores (clinical outcomes) in .csv files

    % Set directory according to model
    krr = strcat(repo, '/results/krr/', scale, '/', timepoint, '/', fc_type, '/');

    % Read in subjects list according to FC type
    if strcmp(fc_type, 'fc_baseline')
        fid = fopen(fullfile(repo, 'data', 'subjects_ses-1.txt'), 'r');
        subj_list = textscan(fid, '%s');
        fclose(fid);
        subj_list = subj_list{1};

    elseif strcmp(fc_type, 'fc_change_bl3m')
        fid = fopen(fullfile(repo, 'data', 'subjects_ses-2.txt'), 'r');
        subj_list = textscan(fid, '%s');
        fclose(fid);
        subj_list = subj_list{1};

    else
        error('Invalid FC type. Please enter "fc_baseline" or "fc_change_bl3m".');
    end

    % Load clinical data for scale and timepoint of interest
    if strcmp(group, 'all')
        group = '';
    end
    
    outcomes = readtable(strcat(repo, '/data/', scale, '_', outcome_type, group, '.txt'), ...
        'FileType', 'text', 'Delimiter', '\t');

    if strcmp(outcome_type, 'change_scores')
        % Convert 'timepoint' column to numeric if needed
        if iscell(outcomes.timepoint)
            outcomes.timepoint = str2double(outcomes.timepoint);
        end
        outcomes = outcomes(outcomes.timepoint == str2double(timepoint), :);
        outcomes.timepoint = [];  % Drop the timepoint column
    end

    % Keep subjects with both fMRI & clinical data
    % Convert to cell array of character vectors for matching
    bids_ids = string(outcomes.BIDS_ID);
    subj_list = subj_list(ismember(subj_list, bids_ids));

    % Filter the outcomes to only those subjects
    outcomes = outcomes(ismember(bids_ids, subj_list), :);

    % Set the row names to 'BIDS_ID' and remove the 'BIDS_ID' column
    outcomes.Properties.RowNames = outcomes.BIDS_ID;
    outcomes.BIDS_ID = [];

    % Save observed clinical outcomes
    writetable(outcomes, strcat(krr, 'outcomes_observed', group, '.csv'), ...
        'WriteRowNames', true);
    
    % Save subject list
    writecell(subj_list, strcat(krr, strcat('subjects', group, '.txt')), 'Delimiter', 'tab');
    
    % Set FC matrix directory
    if strcmp(fc_type, 'fc_baseline')
        folderPath = strcat(repo, '/data/conmats/', preproc, '/', parc, '/ses-1');
    
    elseif strcmp(fc_type, 'fc_change_bl3m')
        folderPath = strcat(repo, '/data/conmats/', preproc, '/', parc, '/change');
        
    end
       
    txtFiles = dir(folderPath);
    
    % Initialize an empty cell array to store the matrices
    n_regions = str2double(extractBefore(parc, '_'));
    
    matrix3D = zeros(n_regions, n_regions, length(subj_list));
    
    % Initialize a counter for valid subjects
    validFileCount = 0;
    
    % Loop through each .txt file
    for i = 1:length(txtFiles)
        % Get the full path and file name of the current file
        fullFileName = fullfile(folderPath, txtFiles(i).name);
        
        % Check if the file name contains a subject from the list
        isValidSubject = false;
        for j = 1:length(subj_list)
            if contains(fullFileName, subj_list{j})
                isValidSubject = true;
                break;
            end
        end
        
        % If the file corresponds to a valid subject, process it
        if isValidSubject
            % Increment valid file counter
            validFileCount = validFileCount + 1;
    
            matrixData = load(fullfile(folderPath, txtFiles(i).name));
    
            % Assign to 3D array
            matrix3D(:,:,validFileCount) = matrixData;
    
        end
    end
    
    % Save the 3D array to a .mat file
    save(strcat(krr, 'FC_', preproc, '_', parc, group), 'matrix3D');

end
