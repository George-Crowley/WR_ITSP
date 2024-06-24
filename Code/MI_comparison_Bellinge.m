% Copyright (c) 2024, George Crowley (gcrowley1@sheffield.ac.uk)
% All rights reserved.

% This source code is licensed under the BSD-style license found in the
% LICENSE file in the root directory of this source tree.

% -----------------------------------------------------------------

function [MI_totalsum_flow, MI_mean, MI_min, MI_max, MI_rule] = MI_comparison_Bellinge(max_number_sensors, number_rand_sims, node_totalsum_flow, covmatrix_training, sigmas, n, Rule_selection)

% Pre define vectors of mutual information for each heuristic

MI_mean = zeros(max_number_sensors, 1);
MI_min = zeros(max_number_sensors, 1);
MI_max = zeros(max_number_sensors, 1);
MI_totalsum_flow = zeros(max_number_sensors, 1);
MI_rule = zeros(max_number_sensors, 1);

for i = 1:n
    covmatrix_training(i, i) = covmatrix_training(i, i) + sigmas;
end

% -------------- Matrix normalisation for determinant calculations --------- %

% Divide each element of the matrix by the mean of values > 0.

minval = mean(covmatrix_training(covmatrix_training > 0));
covmatrix_training = covmatrix_training / minval;

% -------------------------------------------------------------------------- %

% Run first simulation for the Rule based and total sum mutual Information.

for i = 1:max_number_sensors

    det_matrix_totalsum_flow = det(covmatrix_training(node_totalsum_flow(1:i), node_totalsum_flow(1:i)));
    det_matrix_rule = det(covmatrix_training(Rule_selection(1:i), Rule_selection(1:i)));

    MI_totalsum_flow(i) = -i * 0.5 * log(sigmas) + 0.5 * (log(det_matrix_totalsum_flow) - i * log(1/minval));
    MI_rule(i) = -i * 0.5 * log(sigmas) + 0.5 * (log(det_matrix_rule) - i * log(1/minval));

end

% Now run $number_rand_sims$ random simulations for sensor selection and evaluate Mutual information

for j = 1:max_number_sensors

    number_sensors = j;
    perm_matrix_rand = zeros(number_rand_sims, number_sensors+1);

    for i = 1:number_rand_sims
        perm_matrix_rand(i, 1:number_sensors) = sort(randperm(n, number_sensors));
        det_rand_matrix = det(covmatrix_training(perm_matrix_rand(i, 1:number_sensors), perm_matrix_rand(i, 1:number_sensors))); % Calculate determinant of chosen selection
        perm_matrix_rand(i, number_sensors+1) = -number_sensors * 0.5 * (log(sigmas)) + 0.5 * (log(det_rand_matrix) - number_sensors * log(1/minval)); % Calculate MI of selection chosen randomly
    end

    MI_mean(j) = mean(perm_matrix_rand(:, number_sensors+1));
    MI_min(j) = min(perm_matrix_rand(:, number_sensors+1));
    MI_max(j) = max(perm_matrix_rand(:, number_sensors+1));

end

end