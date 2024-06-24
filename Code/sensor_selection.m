% Copyright (c) 2024, George Crowley (gcrowley1@sheffield.ac.uk)
% All rights reserved.

% This source code is licensed under the BSD-style license found in the
% LICENSE file in the root directory of this source tree.

% -----------------------------------------------------------------

function [optimal_sensor_selection_table] = sensor_selection(n, covmatrix, sigmas, choice_number_sensors)

% ----- Including additive white Gaussian noise into covariance matrix ----- %

for i = 1:n
    covmatrix(i, i) = covmatrix(i, i) + sigmas;
    vectorsequence(i) = i;
end

% ----------- Simplifying notation and defining relevant variables ----------%

k = choice_number_sensors;
sensorselection = zeros(k);
dataset_sensorselection_mutualinformation = zeros(k^2, k+1);
placeholder = 0; % Placeholder for table positioning

% -------------- Matrix normalisation for determinant calculations --------- %

% Divide each element of the matrix by the smallest.

chosenval = mean(covmatrix(covmatrix > 0));
covmatrix_adjusted = covmatrix / chosenval;

% ---------------------------- Wait bar ---------------------------------- %

wait = waitbar(0, '1', 'Name', 'Sensor Placement Algorithm', ...
    'CreateCancelBtn', 'setappdata(gcbf,''canceling'',1)');

setappdata(wait, 'canceling', 0);

% -------------------------------------------------------------------------- %

for node = 1:n

    initialnode = node; % Initial node(s)
    nextiterationplacement = initialnode; % Setting up repeated loops.
    counter = 2; % Counter represents number of sensors the search starts with
    % for MI calculations in each loop

    if getappdata(wait, 'canceling')
        break
    end
    waitbar(node/n, wait, sprintf('%12.9f', node/n))
    % -------------------------------------------------------------------------- %

    while length(nextiterationplacement) < k

        choicenodes = setdiff(vectorsequence, nextiterationplacement);

        for choice = 1:length(choicenodes)

            seqchoicenodes = [nextiterationplacement, choicenodes(choice)]; % Possible combination (sequential on each loop)
            seqchoicenodes = sort(seqchoicenodes); % Ordered Possible combination (so doesnt affect det)
            collection_nodes(choice, 1:length(seqchoicenodes)) = seqchoicenodes; % Record keeping for sequence selection

            detmatrix = zeros(length(seqchoicenodes));


            for i = 1:length(seqchoicenodes)

                for j = 1:length(seqchoicenodes)

                    detmatrix(i, j) = covmatrix_adjusted(seqchoicenodes(i), seqchoicenodes(j));

                end

            end

            % Take determinants and calculate MI %

            detval = det(detmatrix);
            mutualinformationvalues(choice) = -counter * 0.5 * log(sigmas) + 0.5 * (log(detval) - counter * log(1/chosenval));

            % counter*log(1/minval)) -- is the adjustment term for normalizing in the determinant
        end

        % ----- Finding max determinant selection ---- %

        [valmax, idxmax] = max(mutualinformationvalues);
        nextiterationplacement = collection_nodes(idxmax, :);
        remainingnodes = setdiff(vectorsequence, nextiterationplacement);

        sensorselection(counter, 1:length(nextiterationplacement)) = nextiterationplacement;
        sensorselection_mutual_information(counter) = valmax;

        counter = counter + 1; % Update number of sensors for MI calculation for next loop.

        % ----- Clear variables for next loop ---- %

        clear collection_nodes;
        clear mutualinformationvalues;


    end

    % Puts first node in sensor selection matrix

    sensorselection(1, 1:length(initialnode)) = initialnode;

    % Save Placements and MI scores asscociated in general matrix from
    % sensorselection and "ss"_mutual_information.

    dataset_sensorselection_mutualinformation(placeholder+1:placeholder+length(sensorselection), 1:length(sensorselection)) = sensorselection;
    dataset_sensorselection_mutualinformation(placeholder+1:placeholder+k, k+1) = transpose(sensorselection_mutual_information);

    % Updating table placement

    placeholder = placeholder + k;

    clear sensorselection;

    mutualinformationtable(2:length(sensorselection_mutual_information)+1, node+1) = transpose(sensorselection_mutual_information);


end

delete(wait) % Delete wait bar

% ------------------------------------------------------------------------------------------- %

% Manually inputting values for first node MI, since we dont include first value in the loop%

for i = 1:n

    mutualinformationtable(2, i+1) = 0.5 * log((1 / (sigmas^(1)))) + 0.5 * (log(det(covmatrix_adjusted(i, i))) - log(1/chosenval));

end


% ------- Creating Tables -------- %

% Create table of information values to find max.

for i = 2:n + 1
    mutualinformationtable(1, i) = i - 1;
end

for i = 2:k + 1
    mutualinformationtable(i, 1) = i - 1;
end

% Find best starting node for each k

for i = 2:k + 1

    [valmax_MI_table, idxmax_MI_table] = max(mutualinformationtable(i, 2:n+1));
    mutualinformationtable(i, n+3) = valmax_MI_table;
    mutualinformationtable(i, n+4) = idxmax_MI_table;
    clear valmax_MI_table;
    clear idxmax_MI_table;

end

% Table for optimal sensor selection (max case) %

for i = 1:k

    optimal_sensor_selection_table(i, 1) = i; %Label number of sensor

end

% Optimal starting sensor - place in table

optimal_sensor_selection_table(:, 3) = mutualinformationtable(2:k+1, n+4);

for i = 1:k

    optimal_sensor_selection_table(i, 5:k+4) = dataset_sensorselection_mutualinformation((optimal_sensor_selection_table(i, 3) - 1)*k+i, 1:k);

end

optimal_sensor_selection_table(:, k+5) = mutualinformationtable(2:k+1, n+3);


end
