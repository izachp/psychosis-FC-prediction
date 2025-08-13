function KRR_null(indir, cbig, perm, preproc, parc, group)
%% Permute clinical outcomes and run KRR over 100 different splits/seeds
% This script has been adapted from a general KRR implementation (Li et al,
% 2019) available at https://github.com/ThomasYeoLab/CBIG/tree/master/utilities/matlab/predictive_models/KernelRidgeRegression

    warning('off', 'all');
    
    cd(indir)
    addpath(genpath(cbig))
    
    if strcmp(group, 'all')
        group = '';
    end
    
    num_test_folds=4; % number of outer folds
    num_inner_folds=4; % number of inner folds (hyper parameter tuning)
    seeds=1:50;
    
    %% Set up other parameters
    data_csv=readtable(strcat(indir, '/outcomes_observed', group, '.csv'));
    setup_param.covariates = 'NONE';
    
    setup_param.feature_mat = load(strcat(indir, '/FC_', preproc, '_', parc, group, '.mat')); 
    fns=fieldnames(setup_param.feature_mat);
    setup_param.feature_mat = setup_param.feature_mat.(fns{1});
    
    setup_param.num_inner_folds = num_inner_folds;
    setup_param.with_bias = 1;
    setup_param.ker_param.type = 'corr';
    setup_param.ker_param.scale = NaN;
    setup_param.lambda_set = [ ...
    0, 1.00000000000000e-05, 0.000100000000000000, 0.00100000000000000, ...
    0.00400000000000000, 0.00700000000000000, 0.0100000000000000, ...
    0.0400000000000000, 0.0700000000000000, 0.100000000000000, ...
    0.400000000000000, 0.700000000000000, 1, 1.50000000000000, ...
    2, 2.50000000000000, 3, 3.50000000000000, 4, 5, 10, 15, 20, ...
    30, 40, 50, 60, 70, 80, 100, 150, 200, 300, 500, 700, ...
    1000, 10000, 100000, 1000000];
    setup_param.threshold_set = [];
    setup_param.metric='corr';
    setup_param.cov_X=[];
    setup_param.save_kernel = 0;
    
    permdir = strcat(indir, '/', preproc, '_', parc, group, '/nulls/', perm); % intermediate directory, is deleted later
    setup_param.outdir = permdir;
    
    % Permute and read in behavioural measures
    rng(str2double(perm))
    data_csv.outcome = data_csv.outcome(randperm(height(data_csv)));
    setup_param.y = data_csv;
    setup_param.y = table2array(setup_param.y(:,2:end));
    
    %% CV split   

    for rep=1:length(seeds)
        seed=seeds(rep);
        setup_param.outstem=['rep' +num2str(rep)];
        CBIG_cross_validation_data_split(strcat(indir, '/subjects', group, '.txt'), 'NONE', 'NONE', ...
        'NONE', num_test_folds, seed, fullfile(permdir, ['randseed_' num2str(seed)]), ',' );
        sub_fold_file = fullfile(permdir, ['randseed_' num2str(seed)], ...
            ['no_relative_' num2str(num_test_folds) '_fold_sub_list.mat']);
        sub_fold = load(sub_fold_file);
        setup_param.sub_fold = sub_fold.sub_fold;
    
        %% Run KRR
        CBIG_KRR_workflow_LITE(setup_param);
    end
    
    % Collate results
    cd(permdir)
	results=zeros(length(seeds), 1);
    for s=1:length(seeds)
        load(sprintf('%s/final_result_rep%d.mat', perm, s))
        results(s)=mean(optimal_acc);
    end
    
    % Save predicted-observed correlations for all KRR repeats
    cd(strcat(indir, '/', preproc, '_', parc, group, '/nulls'))
    
    filename=strcat('r_vals_perm', perm, '.txt');
    writematrix(results, filename, 'Delimiter', 'tab')

    % Delete intermediate folder to save space
    rmdir(permdir, 's');

end
