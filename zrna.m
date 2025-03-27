% RNA Structure Quality Assessment using Z-score Analysis
% This script calculates standardized scores (Z-scores) for various RNA structure quality metrics
% and combines them into a comprehensive ZRNA score using weighted averaging.

% Read the input CSV file containing RNA structure quality metrics
% The file should contain columns for RMSD, INF_all, TM, GDT, clash, and lDDT
data = readtable('data.csv', 'ReadVariableNames', true);

% Extract the relevant metrics for Z-score calculation
% Note: GDT appears twice in the original data, we'll use it once
metrics = data{:, {'RMSD', 'INF_all', 'TM', 'GDT', 'clash', 'lDDT'}};

% Calculate Z-scores for each metric using MATLAB's built-in zscore function
% Z-score = (X - μ) / σ, where μ is the mean and σ is the standard deviation
% This standardizes all metrics to have mean=0 and standard deviation=1
standardized_metrics = zscore(metrics);

% Create new column names for the Z-scores by adding 'Zscore_' prefix
% This helps distinguish original metrics from their standardized versions
new_column_names = strcat('Zscore_', metrics.Properties.VariableNames);

% Add the standardized metrics to the original data table
% Each metric gets its own Z-score column
for i = 1:length(new_column_names)
    data.(new_column_names{i}) = standardized_metrics(:, i);  
end

% Define weights for each metric in the final ZRNA score
% These weights reflect the relative importance of each metric:
% - TM (Template Modeling): 1/3 (33.3%) - High weight for overall structure similarity
% - GDT (Global Distance Test): 1/3 (33.3%) - High weight for global structure accuracy
% - lDDT (local Distance Difference Test): 1/8 (12.5%) - Medium weight for local accuracy
% - INF_all (Information Score): 1/8 (12.5%) - Medium weight for information content
% - clash: 1/12 (8.3%) - Lower weight for steric clashes
weights = struct('TM', 1/3, 'GDT', 1/3, 'lDDT', 1/8, 'INF_all', 1/8, 'clash', 1/12);

% Extract Z-scores for each metric
TM_Z = data.Zscore_TM;        % Template Modeling Z-score
GDT_Z = data.Zscore_GDT;      % Global Distance Test Z-score
lDDT_Z = data.Zscore_lDDT;    % local Distance Difference Test Z-score
INF_all_Z = data.Zscore_INF_all;  % Information Score Z-score
clash_Z = data.Zscore_clash;  % Clash Z-score

% Calculate the combined ZRNA score using weighted averaging
% Note: clash_Z is subtracted because lower clash values are better
% The weights sum to 1 (1/3 + 1/3 + 1/8 + 1/8 + 1/12 = 1)
combined_Z = (weights.TM * TM_Z + weights.GDT * GDT_Z + ...
               weights.lDDT * lDDT_Z + weights.INF_all * INF_all_Z - ...
               weights.clash * clash_Z);

% Add the combined ZRNA score to the data table
data.Combined_Zscore = combined_Z;

% The resulting Combined_Zscore represents the overall quality of the RNA structure
% Higher values indicate better overall structure quality
% The score is standardized (mean=0, std=1) due to the Z-score normalization 
