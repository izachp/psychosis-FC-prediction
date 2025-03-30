function make_KRR_inputs(repo, scale, timepoint, fc_type)
%% Concatenate FC matrices into .mat and save subject list and change scores (clinical outcomes) in .csv files

    % Set directory according to model
    krr = strcat(repo, '/krr/', scale, '/', timepoint, '/', fc_type, '/');
    
    % Read in subject list via observed change scores file in CPM results
    dataTable = readtable(strcat(repo, '/cpm/', scale, '/', timepoint, '/', fc_type, '/change_observed.xlsx'));
    subjectIDs = dataTable{:,1};
    
    % Save subject list & clinical outcomes
    writecell(subjectIDs, strcat(krr, 'subjects.txt'), 'Delimiter', 'tab');
    writetable(dataTable(:,1:2), strcat(krr, 'change-scores.csv'), 'Delimiter', 'tab');
    
    % Set FC matrix directory
    if strcmp(fc_type, 'fc_baseline')
        folderPath = strcat(repo, '/conmats/dt_AROMA_8Phys-4GMR_bpf/ses-1');
        filePattern = fullfile(folderPath, '*ses-1.txt');
    
    elseif strcmp(fc_type, 'fc_change_bl3m')
        folderPath = strcat(repo, '/conmats/dt_AROMA_8Phys-4GMR_bpf/change');
        filePattern = fullfile(folderPath, '*.txt');
    
    end
       
    txtFiles = dir(filePattern);
    
    % Initialize an empty cell array to store the matrices
    matrix3D = zeros(328, 328, length(subjectIDs));
    
    % Initialize a counter for valid subjects
    validFileCount = 0;
    
    % Loop through each .txt file
    for i = 1:length(txtFiles)
        % Get the full path and file name of the current file
        fullFileName = fullfile(folderPath, txtFiles(i).name);
        
        % Check if the file name contains a subject from the list
        isValidSubject = false;
        for j = 1:length(subjectIDs)
            if contains(fullFileName, subjectIDs{j})
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
    save(strcat(krr, 'FC'), 'matrix3D');

end
