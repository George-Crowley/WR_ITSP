% Copyright (c) 2024, George Crowley (gcrowley1@sheffield.ac.uk)
% All rights reserved.

% This source code is licensed under the BSD-style license found in the
% LICENSE file in the root directory of this source tree.

% -----------------------------------------------------------------

function [prediction_matrix, validation_un_observed, error_matrix, NMSE] = GRNNET_estimation1(n, max_chosen_selection, seq, training_un_observed, training_observed, validation_un_observed, validation_observed)

training_un_observed(:, max_chosen_selection) = []; %<- Set chosen sensors columns to be deleted.
training_observed(:, setdiff(seq, max_chosen_selection)) = []; %<- Set Un_selected sensors columns to be deleted.

validation_un_observed(:, max_chosen_selection) = []; %<- Set chosen sensors columns to be deleted.
validation_observed(:, setdiff(seq, max_chosen_selection)) = []; %<- Set Un_selected sensors columns to be deleted.

% ------------------------------------------------------------ %

training_observed_normalized = training_observed;
validation_observed_normalized = validation_observed;

for i = 1:width(training_observed_normalized)

    stdev = std(training_observed(:, i));

    if stdev == 0 % We check because some nodes can have 0's for all time realizations in the simulated SWMM data.

        stdev = 1;

    end

    training_observed_normalized(:, i) = (training_observed_normalized(:, i) - mean(training_observed(:, i))) / stdev;
    validation_observed_normalized(:, i) = (validation_observed_normalized(:, i) - mean(training_observed(:, i))) / stdev;

    clear stdev;
end

spread = .5;

p = training_observed_normalized';
t = training_un_observed';

if iscell(p), p = cell2mat(p); end
if iscell(t), t = cell2mat(t); end

% Dimensions
[R, Q] = size(p);
[S, Q] = size(t);
net = network(1, 2, [1; 0], [1; 0], [0, 0; 1, 0], [0, 1]);

% Simulation
net.inputs{1}.size = R;
net.layers{1}.size = Q;
net.inputWeights{1, 1}.weightFcn = 'dist';
net.layers{1}.netInputFcn = 'netprod';
net.layers{1}.transferFcn = 'radbasn';
net.layers{2}.size = S;
net.layerWeights{2, 1}.weightFcn = 'dotprod';
net.b{1} = zeros(Q, 1) + sqrt(1/2) / spread;
net.iw{1, 1} = p';
net.lw{2, 1} = t;

prediction_matrix = net(validation_observed_normalized')';

% ---------------------------------------- %

% -- Call error function -- %

[error_matrix, NMSE] = error_calculations(prediction_matrix, validation_un_observed);


end
