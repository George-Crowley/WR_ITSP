% Copyright (c) 2024, George Crowley (gcrowley1@sheffield.ac.uk)
% All rights reserved.

% This source code is licensed under the BSD-style license found in the
% LICENSE file in the root directory of this source tree.

% -----------------------------------------------------------------

function [prediction_matrix, validation_un_observed, error_matrix, NMSE] = GLM_estimation1(n, max_chosen_selection, seq, training_un_observed, training_observed, validation_un_observed, validation_observed)

% ---------------- GLM / Polynomial regession ------------------%

training_un_observed(:, max_chosen_selection) = []; %<- Set chosen sensors columns to be deleted.
training_observed(:, setdiff(seq, max_chosen_selection)) = []; %<- Set Un_selected sensors columns to be deleted.

validation_un_observed(:, max_chosen_selection) = []; %<- Set chosen sensors columns to be deleted.
validation_observed(:, setdiff(seq, max_chosen_selection)) = []; %<- Set Un_selected sensors columns to be deleted.

% -- Data matrix -- %

% -- Input data -- %

X_data_matrix = (zeros(height(training_observed), width(training_observed)+1));
% Insert beta_0;
X_data_matrix(:, 1) = 1;
% Insert observed training data
X_data_matrix(:, 2:width(training_observed)+1) = training_observed;

% ----------------------- %

X_data_matrix_validation = (zeros(height(validation_observed), width(validation_observed)+1));
% Insert beta_0;
X_data_matrix_validation(:, 1) = 1;
% Insert observed training data
X_data_matrix_validation(:, 2:width(validation_observed)+1) = validation_observed;


% ------------------ Calculate beta coefficients from training data

for i = 1:width(validation_un_observed)

    sensor_i_beta_coefficients = regress(training_un_observed(:, i), X_data_matrix);
    prediction_matrix(:, i) = X_data_matrix_validation * sensor_i_beta_coefficients;
    clear sensor_i_beta_coefficients;

end

% -- Call error function -- %

[error_matrix, NMSE] = error_calculations(prediction_matrix, validation_un_observed);


end