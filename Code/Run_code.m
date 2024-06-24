% Copyright (c) 2024, George Crowley (gcrowley1@sheffield.ac.uk)
% All rights reserved.

% This source code is licensed under the BSD-style license found in the
% LICENSE file in the root directory of this source tree.

% -----------------------------------------------------------------

clear all
close all
rng(12345);

% -- Make the folder with files the path - create some folders to save results in

mkdir Results_folder
mkdir Results_folder Sensor_selection
mkdir Results_folder Table_results_estimation
mkdir Results_folder Error_matrix_k_max

% Flow data - training and validation

load('sim_1_sim_2_merged_flow.mat');
load('sim_3_sim_4_merged_flow.mat');

data_matrix_Y1 = time_series_sim_1_and_2_master; % Training
data_matrix_Y2 = time_series_sim_3_and_4_master; % Validation

clear time_series_sim_1_and_2_master;
clear time_series_sim_3_and_4_master;

% Remove times -----------------------

training_data = table2array(data_matrix_Y1(:, 2:width(data_matrix_Y1)));
validation_data = table2array(data_matrix_Y2(:, 2:width(data_matrix_Y2)));

% -----------------------------

if isstring(training_data) == 1

    training_data = str2double(training_data);

end

if isstring(validation_data) == 1

    validation_data = str2double(validation_data);

end

% -----------------------------

covmatrix_training = cov(training_data);
n = width(covmatrix_training);

% Variance sigma^2 introduced from noise due to sensor
sigmas = 10^(-4);

% Input max number of sensor placements the algorithm stops at (note the rule based placement only goes to 250) - keep this as a multiple of 25.
max_number_sensors = 250;

if max_number_sensors > 250 | mod(max_number_sensors, 25)

    promptMessage_est = sprintf('To proceed with the sensor placement algorithm, make sure you have set max_number_sensors less than 250 and a multiple of 25 so that the rest of the script runs, press cancel to amend');
    button = questdlg(promptMessage_est, 'Continue', 'Continue', 'Cancel', 'Continue');
    if strcmpi(button, 'Cancel')
        return; % Or break or continue
    end

end

% Check to see if the file has already been run - note, if you want to run
% the sensor placement algorithm with more sensors, you'll need to amend
% the loop or change the name of the old file.

if isfile("Results_folder/Sensor_selection/sensor_selection.mat") == 0

    % Run sensor placement (Algorithm 1)
    [optimal_sensor_selection_table] = sensor_selection(n, covmatrix_training, sigmas, max_number_sensors);

    % Column one of optimal_sensor_selection_table is k (the number of sensors). Column 3 is the initial node that obtained the best MI from the modified one-step greedy algo.

    save("Results_folder/Sensor_selection/sensor_selection.mat", "optimal_sensor_selection_table") % save table so the script doesn't need to be run again
    disp("Sensor placement complete")

else

    load("Results_folder/Sensor_selection/sensor_selection.mat");

end
% --------- Generate script for comparing mutual information between:

% Algorithm 1 sensor placement
% Random sensor placement (min, mean and max of i.e. 10,000 simulations) - see number_rand_sims
% Heuristics sensor placement (total sum)

% First rank the heuristics ---- (i+1) in the array to account for date/time below
max_sum_flow_vector = zeros(1, width(training_data));

for i = 1:n

    max_sum_flow_vector(i) = sum((training_data(:, i)));

end

[order_totalsum_flow, node_totalsum_flow] = sort(max_sum_flow_vector, "descend");

% Node_rank_max_flow and node_totalsum_flow contain the vector of ranked
% sensors according to the calculation.

% Before comparing Mutual information between each heuristic, define how
% many random simulations will be compared for min,mean and max of each
% simulation for each k \leq max_chosen_sensor sensor placement.

number_rand_sims = 1000;

%load Rule based SS - Note it does not exceed 250 placements.

load('Data/Rule_sensor_selection.mat');

% Run MI for the heuristics

[MI_totalsum_flow, MI_mean, MI_min, MI_max, MI_rule] = MI_comparison_Bellinge(max_number_sensors, number_rand_sims, node_totalsum_flow, covmatrix_training, sigmas, n, Rule_selection);

% -- This is FIGURE 5 - Bellinge plot for mutual information comparisons -- %

t = tiledlayout('flow', 'TileSpacing', 'compact');

nexttile
plot((1:max_number_sensors), optimal_sensor_selection_table(:, max_number_sensors+5), 'DisplayName', 'Algorithm 1 placement', 'color', 'red')
hold on
plot((1:max_number_sensors), MI_max, 'displayname', 'Max of random placements', 'color', 'blue')
hold on
plot((1:max_number_sensors), MI_min, 'DisplayName', "Min of random placements", 'color', 'black')
hold on
plot((1:max_number_sensors), MI_mean, 'DisplayName', "Mean of random placements", 'color', 'green')
hold on
plot((1:max_number_sensors), MI_rule, 'DisplayName', "Rule based placement", 'color', 'magenta')
hold on
plot((1:max_number_sensors), MI_totalsum_flow, 'DisplayName', "Max total sum placement ", 'color', 'cyan')
xlabel('Number of sensors', "Fontsize", 13)
ylabel('Mutual information value', "Fontsize", 13)
xlim([0, max_number_sensors])
legend('Orientation', 'Horizontal', 'NumColumns', 3, 'Location', 'southoutside', 'Fontsize', 12)

nexttile
plot((optimal_sensor_selection_table(:, max_number_sensors+5) - MI_max)./optimal_sensor_selection_table(:, max_number_sensors+5), 'DisplayName', 'Max of random placements', 'color', 'blue')
hold on
plot((optimal_sensor_selection_table(:, max_number_sensors+5) - MI_min)./optimal_sensor_selection_table(:, max_number_sensors+5), 'DisplayName', 'Min of random placements', 'color', 'black')
hold on
plot((optimal_sensor_selection_table(:, max_number_sensors+5) - MI_mean)./optimal_sensor_selection_table(:, max_number_sensors+5), 'DisplayName', 'Mean of random placements', 'color', 'green')
hold on
plot((optimal_sensor_selection_table(:, max_number_sensors+5) - MI_rule)./optimal_sensor_selection_table(:, max_number_sensors+5), 'DisplayName', 'Rule based placement', 'color', 'magenta')
hold on
plot((optimal_sensor_selection_table(:, max_number_sensors+5) - MI_totalsum_flow)./optimal_sensor_selection_table(:, max_number_sensors+5), 'DisplayName', 'Max total sum placement', 'color', 'cyan')
xlabel('Number of sensors', "Fontsize", 13)
ylabel('MI gain', "Fontsize", 13)
xlim([0, max_number_sensors])
legend('Orientation', 'Horizontal', 'NumColumns', 3, 'Location', 'southoutside', 'Fontsize', 12)
set(gcf, 'Position', [200, 200, 600, 600])
saveas(gcf, "Results_folder/Figure5.fig");

%%

% ----------------------------------------------------------------- %

% Prepare for running estimation scripts

seq = (1:n); % For comparing locations unselected for sensor allocation in the below loop.

% Decide which 'k' sensor placements you are interested in for estimation: i.e. [5,10,20,150,....].

% Note, pick only values in (25:25:250) to pick the random matrices used in
% the paper. If you want to make your own, you will need to do so.

k_list_vector = (25:25:max_number_sensors);

MI_in_table = 5 + max(k_list_vector); % Table positioning of mutual information value

full_k_list_table_results = zeros(8*length(k_list_vector), 3);

k_max = 50; % -> For plot of FIGURES 7 and 8. If you change this value, the plots will not be the same.

if isempty(intersect(k_max, k_list_vector)) == 1

    disp("Make sure k_max is in the vector k_list_vector, set k_max = max chosen sensors")
    k_max = max(k_list_vector);

end

% Random placements are already saved and are loaded in the run_estimation script.

% ----------------------------------------

promptMessage_est = sprintf('To proceed with the estimation scripts, the first 4000 realizations are taken for training and validation sets. If you want to change these values, press cancel and change them, otherwise continue?');
button = questdlg(promptMessage_est, 'Continue', 'Continue', 'Cancel', 'Continue');
if strcmpi(button, 'Cancel')
    return; % Or break or continue
end

% NOTE: If you try to run the GRNN estimation with the full training data set - you will likely error due to insufficient RAM. Either run the simulation on a HPC or remove data from the data sets. I.e.
% Take the first 4000 data points. The GLM should run without the GRNN, so you can remove this from all the relevant places if you so wish.

training_data = training_data(1:4000, :);
validation_data = validation_data(1:4000, :);

% parfor works here if you have GPUs - just define your pool. Relevant data is saved in the scripts

for k = 1:length(k_list_vector)

    % In this case, the function Run_estimation will pull out relevant matrices to histogram

    if k_list_vector(k) == k_max

        [~, error_matrix_GLM_algo, error_matrix_GRRNET_algo] = Run_estimation(k_list_vector, optimal_sensor_selection_table, training_data, validation_data, k, n, MI_in_table, seq, node_totalsum_flow, Rule_selection);

        % Save error tables for consistency of estimator plots
        parsave_error_matrix(sprintf('Results_folder/Error_matrix_k_max/error_matrix_GLM_algo_k_%d.mat', k_max), error_matrix_GLM_algo);
        parsave_error_matrix(sprintf('Results_folder/Error_matrix_k_max/error_matrix_GRNNET_algo_k_%d.mat', k_max), error_matrix_GRRNET_algo);

    else

        [~, ~, ~] = Run_estimation(k_list_vector, optimal_sensor_selection_table, training_data, validation_data, k, n, MI_in_table, seq, node_totalsum_flow, Rule_selection);


    end

    sprintf("Sensor placement estimation with k = "+k_list_vector(k)+" has finished")

end

clear k;

% --------- Load results ------------------

file_name = 'Results_folder/Table_results_estimation/Table_estimation_k_%d.mat';

for k = 1:length(k_list_vector)

    load(sprintf(file_name, k_list_vector(k)))
    Results_table_Bellinge((k - 1)*9+1:k*9, 1:3) = Table_results;
    clear Table_results;

end

clear filename;
save("Results_folder/Full_table_results.mat", "Results_table_Bellinge");

% ------------Consistency of estimators plot ----------------- %

% Extract data from loading error_matrix_tables... saved

file_name_GLM = 'Results_folder/Error_matrix_k_max/error_matrix_GLM_algo_k_%d.mat';
file_name_GRNNET = 'Results_folder/Error_matrix_k_max/error_matrix_GRNNET_algo_k_%d.mat';

load(sprintf(file_name_GLM, k_max))
error_matrix_GLM_k_max = error_matrix.^2;
clear error_matrix; clear file_name_GLM;

load(sprintf(file_name_GRNNET, k_max))
error_matrix_GRNNET_k_max = error_matrix.^2;
clear error_matrix; clear file_name_GRNNET;

bins = 100;

% -- This is FIGURE 7 (note, that the plot will look different if you have not run the full training and validation data sets)

t = tiledlayout('flow', 'TileSpacing', 'compact');

nexttile
histogram(log10(error_matrix_GLM_k_max(:)), 'Normalization', 'pdf', 'EdgeColor', 'none')
xlabel('Log_{10} square error for each GLM data point estimate for unmonitored nodes', 'FontSize', 15)
ylabel('PDF value', 'FontSize', 15)
xlim([-18, -2])
nexttile
histogram(log10(error_matrix_GRNNET_k_max(:)), 'Normalization', 'pdf', 'EdgeColor', 'none')
xlabel('Log_{10} square error for each GRNN data point estimate for unmonitored nodes', 'FontSize', 15)
ylabel('PDF value', 'FontSize', 15)
xlim([-18, -2])
saveas(gcf, "Results_folder/Figure7.fig");
% ------------------------------------------------------------ %

% Load Bellinge shapefiles

load('Data/Bellinge_sf.mat')

junctions_table = struct2table(junctions);
outfalls_table = struct2table(outfalls);
storages_table = struct2table(storages);

junctions_table(:, 5:end) = [];
outfalls_table(:, 5:end) = [];
storages_table(:, 5:end) = [];

% Merge tables -> NOTE THE BELOW TABLE IS IN THE SAME ORDER AS THE
% TIMESERIES DATA IN Sim_1_Sim_2_flow.mat ETC. This is important otherwise the sensor seletion shown will not be correct.

merged_nodes_table = [junctions_table; outfalls_table; storages_table];
merged_nodes_table = table2struct(merged_nodes_table);

% Plot k_max sensor placement on Bellinge

k_seq = (1:n)';
nodes_to_delete = zeros(1, n-k_max);
nodes_to_delete(1, 1:width(nodes_to_delete)) = intersect(k_seq, setdiff(k_seq, optimal_sensor_selection_table(k_max, 5:5+k_max-1)));

nodes_k_max = merged_nodes_table;
nodes_k_max(nodes_to_delete(1, 1:n-k_max)) = [];

% -- This is FIGURE 8 -----

figure
b(1) = mapshow(links, 'Color', 'blue');
hold on
b(2) = mapshow(junctions, 'MarkerEdgeColor', 'black', 'Marker', '.', 'MarkerSize', 9);
hold on
b(3) = mapshow(orifices, 'MarkerEdgeColor', 'green', 'Marker', '*', 'MarkerSize', 9);
hold on
b(4) = mapshow(outfalls, 'MarkerEdgeColor', '[0.8500 0.3250 0.0980]', 'Marker', 'o', 'MarkerSize', 9);
hold on
b(5) = mapshow(pumps, 'MarkerEdgeColor', 'red', 'Marker', '.', 'MarkerSize', 10);
hold on
b(6) = mapshow(raingages, 'MarkerEdgeColor', '[0.4660 0.6740 0.1880]', 'Marker', '^', 'MarkerSize', 8);
hold on
b(7) = mapshow(storages, 'MarkerEdgeColor', 'black', 'Marker', 'square', 'MarkerSize', 10);
hold on
b(8) = mapshow(weirs, 'MarkerEdgeColor', 'magenta', 'Marker', 'x', 'MarkerSize', 10);
hold on
b(9) = mapshow(nodes_k_max, 'MarkerEdgeColor', "black", 'DisplayName', 'Algorithm 1 sensor placement', 'MarkerSize', 7, 'Marker', 'sq', 'MarkerFaceColor', 'yellow');
axis off
legend(b, {'Sewer pipes', 'Manholes', 'Orifices', 'Outfalls', 'Pumps', 'Raingages', 'Storages', 'Weirs', 'Algorithm 1 placement \newlinefor k = ' + string(k_max)}, 'Location', 'NorthWest', 'FontSize', 20)
set(gcf, 'Position', get(0, 'Screensize'));
saveas(gcf, "Results_folder/Figure8.fig");

% ------- Function for saving files in the par/for loop

function parsave_table(fname, Table_results)
save(fname, 'Table_results')
end

function parsave_error_matrix(fname, error_matrix)
save(fname, 'error_matrix')
end

function parsave_rand(fname, random_selection_matrix)
save(fname, 'random_selection_matrix')
end
