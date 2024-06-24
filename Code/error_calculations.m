% Copyright (c) 2024, George Crowley (gcrowley1@sheffield.ac.uk)
% All rights reserved.

% This source code is licensed under the BSD-style license found in the
% LICENSE file in the root directory of this source tree.

% -----------------------------------------------------------------

function [error_matrix, NMSE] = error_calculations(estimated_matrix_validation_un_observed, validation_un_observed)

% -- Calculate error matrix -- %

error_matrix = estimated_matrix_validation_un_observed - validation_un_observed;

% -- Calculate square error for each entry -- %

error_matrix_sq = error_matrix.^2;

% -- Normalised error metrics -- %

% Calculate energy of true un_observed data from the validation set

energy_matrix = validation_un_observed.^2;
energy = sum(energy_matrix, 'all');

% -- Calculate normalised total square error of all estimated nodes and readings -- %

NMSE = sum(error_matrix_sq, "all") / energy;

end