% Copyright (c) 2024, George Crowley (gcrowley1@sheffield.ac.uk)
% All rights reserved.

% This source code is licensed under the BSD-style license found in the
% LICENSE file in the root directory of this source tree.

% -----------------------------------------------------------------

function [Table_results, error_matrix_GLM_algo, error_matrix_GRNNET_algo] = Run_estimation(k_list_vector, optimal_sensor_selection_table, training_data, validation_data, k, n, MI_in_table, seq, node_totalsum_flow, Rule_selection)

warning('off', 'all')

% --------- Data which is to be used in each function ----- %

training_un_observed = training_data;
training_observed = training_data;

validation_un_observed = validation_data;
validation_observed = validation_data;

% Note these will be modified dependent on the sensor selection^^.

%  ---------------------------------------------------

max_chosen_selection = optimal_sensor_selection_table(k_list_vector(k), 5:5+k_list_vector(k)-1);
MI_algo = optimal_sensor_selection_table(k_list_vector(k), MI_in_table);

choice_number_sensors = k_list_vector(k);

% Max chosen selection is the k best selection.

[prediction_matrix_GLM_algo, validation_un_observed_algo, error_matrix_GLM_algo, NMSE_GLM_algo] = GLM_estimation1(n, max_chosen_selection, seq, training_un_observed, training_observed, validation_un_observed, validation_observed);
[prediction_matrix_GRNNET_algo, ~, error_matrix_GRNNET_algo, NMSE_GRNNET_algo] = GRNNET_estimation1(n, max_chosen_selection, seq, training_un_observed, training_observed, validation_un_observed, validation_observed);

% --------------------------------------------------------------%


% Biggest K sensor selection for total flow in training data %
chosen_totalsum_flow = sort(node_totalsum_flow(1:choice_number_sensors));

[prediction_matrix_GLM_totalsum_flow, validation_un_observed_totalsum_flow, error_matrix_GLM_totalsum_flow, NMSE_GLM_totalsum_flow] = GLM_estimation1(n, chosen_totalsum_flow, seq, training_un_observed, training_observed, validation_un_observed, validation_observed);
[prediction_matrix_GRNNET_totalsum_flow, ~, error_matrix_GRNNET_totalsum_flow, NMSE_GRNNET_totalsum_flow] = GRNNET_estimation1(n, chosen_totalsum_flow, seq, training_un_observed, training_observed, validation_un_observed, validation_observed);


% --------------------------------------------------------------%

% Rule based selection %

chosen_rule = sort(Rule_selection(1:choice_number_sensors));

[prediction_matrix_GLM_rule, validation_un_observed_rule, error_matrix_GLM_rule, NMSE_GLM_rule] = GLM_estimation1(n, chosen_rule, seq, training_un_observed, training_observed, validation_un_observed, validation_observed);
[prediction_matrix_GRNNET_rule, ~, error_matrix_GRNNET_rule, NMSE_GRNNET_rule] = GRNNET_estimation1(n, chosen_rule, seq, training_un_observed, training_observed, validation_un_observed, validation_observed);

% ------------------------------------------------------------------- %
% --------------------- Random selections --------------------------- %
% Load random matrix from Data folder

file_name_rand = 'Data/Random_selections/random_sensor_selections_k_%d.mat';
load(sprintf(file_name_rand, k_list_vector(k)), 'random_selection_matrix');
clear file_name_rand;

% Calculate GLM for each random selection

[prediction_matrix_GLM_rand1, validation_un_observed_rand1, error_matrix_GLM_rand1, NMSE_GLM_rand1] = GLM_estimation1(n, random_selection_matrix(1, :), seq, training_un_observed, training_observed, validation_un_observed, validation_observed);
[prediction_matrix_GLM_rand2, validation_un_observed_rand2, error_matrix_GLM_rand2, NMSE_GLM_rand2] = GLM_estimation1(n, random_selection_matrix(2, :), seq, training_un_observed, training_observed, validation_un_observed, validation_observed);
[prediction_matrix_GLM_rand3, validation_un_observed_rand3, error_matrix_GLM_rand3, NMSE_GLM_rand3] = GLM_estimation1(n, random_selection_matrix(3, :), seq, training_un_observed, training_observed, validation_un_observed, validation_observed);
[prediction_matrix_GLM_rand4, validation_un_observed_rand4, error_matrix_GLM_rand4, NMSE_GLM_rand4] = GLM_estimation1(n, random_selection_matrix(4, :), seq, training_un_observed, training_observed, validation_un_observed, validation_observed);
[prediction_matrix_GLM_rand5, validation_un_observed_rand5, error_matrix_GLM_rand5, NMSE_GLM_rand5] = GLM_estimation1(n, random_selection_matrix(5, :), seq, training_un_observed, training_observed, validation_un_observed, validation_observed);

% Calculate GRNNET for each random selection

[prediction_matrix_GRNNET_rand1, ~, error_matrix_GRNNET_rand1, NMSE_GRNNET_rand1] = GRNNET_estimation1(n, random_selection_matrix(1, :), seq, training_un_observed, training_observed, validation_un_observed, validation_observed);
[prediction_matrix_GRNNET_rand2, ~, error_matrix_GRNNET_rand2, NMSE_GRNNET_rand2] = GRNNET_estimation1(n, random_selection_matrix(2, :), seq, training_un_observed, training_observed, validation_un_observed, validation_observed);
[prediction_matrix_GRNNET_rand3, ~, error_matrix_GRNNET_rand3, NMSE_GRNNET_rand3] = GRNNET_estimation1(n, random_selection_matrix(3, :), seq, training_un_observed, training_observed, validation_un_observed, validation_observed);
[prediction_matrix_GRNNET_rand4, ~, error_matrix_GRNNET_rand4, NMSE_GRNNET_rand4] = GRNNET_estimation1(n, random_selection_matrix(4, :), seq, training_un_observed, training_observed, validation_un_observed, validation_observed);
[prediction_matrix_GRNNET_rand5, ~, error_matrix_GRNNET_rand5, NMSE_GRNNET_rand5] = GRNNET_estimation1(n, random_selection_matrix(5, :), seq, training_un_observed, training_observed, validation_un_observed, validation_observed);

% Merge all results in a table, and save table.

list_choice = [k_list_vector(k), "Algorithm 1", "Rule based selection", "Total sum of flows selection", "Random selection 1", "Random selection 2", "Random selection 3", "Random selection 4", "Random selection 5"]';
GLM_norm_MSE = ["NMSE - GLM", NMSE_GLM_algo, NMSE_GLM_rule, NMSE_GLM_totalsum_flow, NMSE_GLM_rand1, NMSE_GLM_rand2, NMSE_GLM_rand3, NMSE_GLM_rand4, NMSE_GLM_rand5]';
GRNNET_norm_MSE = ["NMSE - GRNN", NMSE_GRNNET_algo, NMSE_GRNNET_rule, NMSE_GRNNET_totalsum_flow, NMSE_GRNNET_rand1, NMSE_GRNNET_rand2, NMSE_GRNNET_rand3, NMSE_GRNNET_rand4, NMSE_GRNNET_rand5]';

Table_results = table(list_choice, GLM_norm_MSE, GRNNET_norm_MSE);

parsave_table(sprintf('Results_folder/Table_results_estimation/Table_estimation_k_%d.mat', k_list_vector(k)), Table_results);

    function parsave_table(fname, Table_results)
        save(fname, 'Table_results')
    end


end
