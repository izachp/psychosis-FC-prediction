function run_KRR_repeats(outdir, cbig)
%% Run KRR over 100 different splits/seeds
% This script has been adapted from a general KRR implementation (Li et al,
% 2019) available at https://github.com/ThomasYeoLab/CBIG/tree/master/utilities/matlab/predictive_models/KernelRidgeRegression

    cd(outdir)
    addpath(genpath(cbig))
    
    num_test_folds=4; % number of outer folds
    num_inner_folds=4; % number of inner folds (hyper parameter tuning)
    seeds=1:100;
    
    %% Read predictor variables
    
    setup_param.y = readtable(strcat(outdir, '/change-scores.csv')); % csv/txt of predictors with headings (should be in the same subject order as the FC)
    setup_param.y = table2array(setup_param.y(:,2:end));
    
    %% Set up other parameters
    
    setup_param.covariates = 'NONE';
    
    setup_param.feature_mat = load(strcat(outdir, '/FC.mat')); % FC is saved in the .mat as roi x roi x subjects
    fns=fieldnames(setup_param.feature_mat);
    setup_param.feature_mat = setup_param.feature_mat.(fns{1});
    
    setup_param.num_inner_folds = num_inner_folds;
    setup_param.with_bias = 1;
    setup_param.ker_param.type = 'corr';
    setup_param.ker_param.scale = NaN;
    setup_param.lambda_set=load('stable_projects/preprocessing/Li2019_GSR/KernelRidgeRegression/lambda_set.mat');
    setup_param.lambda_set = setup_param.lambda_set.lambda_set;
    setup_param.threshold_set = [];
    setup_param.metric='corr';
    setup_param.cov_X=[];
    
    outdir = strcat(outdir, '/results');
    setup_param.outdir = outdir;
    
    %% CV split
    
    for rep=1:length(seeds)
        seed=seeds(rep);
        setup_param.outstem=['rep' num2str(rep)];
        CBIG_cross_validation_data_split(strcat(outdir, '/subjects.txt'), 'NONE', 'NONE', ...
        'NONE', num_test_folds, seed, fullfile(outdir, ['randseed_' num2str(seed)]), ',' );
        sub_fold_file = fullfile(outdir, ['randseed_' num2str(seed)], ...
            ['no_relative_' num2str(num_test_folds) '_fold_sub_list.mat']);
        sub_fold = load(sub_fold_file);
        setup_param.sub_fold = sub_fold.sub_fold;
    
    %% Run KRR
        CBIG_KRR_workflow_LITE(setup_param)
    end

    % Collate results
    for s=1:length(seeds)
        load(sprintf('results/final_result_rep%d.mat', s));
        
        if s==1
            results=zeros(length(seeds), 1);
            predictions=zeros(length(y_predict_concat), length(seeds));
        end
        
        results(s)=mean(optimal_acc);
        predictions(:,s)=y_predict_concat;
    end
    
    % Save predicted-observed correlations for all KRR repeats
    filename='r_vals.xlsx';
    writematrix(results, filename);
    
    % Save individual predictions for all KRR repeats
    filename='change_predictions.xlsx';
    writematrix(predictions, filename);

end
