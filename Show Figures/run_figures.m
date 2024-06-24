% Copyright (c) 2024, George Crowley (gcrowley1@sheffield.ac.uk)
% All rights reserved.

% This source code is licensed under the BSD-style license found in the
% LICENSE file in the root directory of this source tree.

% ---------------------------------------------------------------------

% This code is intended to simulate all Figures in Case Study 2 without running the
% simulations. Each figure has its own section, so please go to the section
% you want to look at and use "Run Section".

% ---------------------------------------------------------------------

% Figue 4 - map of Bellinge

clear all
close all

load("Data/Bellinge_sf.mat"); % load shape-files already saved in structured arrays

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
axis off
legend(b, {'Sewer pipes', 'Manholes', 'Orifices', 'Outfalls', 'Pumps', 'Raingages', 'Storages', 'Weirs'}, 'Location', 'NorthWest', 'FontSize', 20)

% ---------------------------------------------------------------------

%%
% Figure 5 - Mutual information plots

clear all
close all

% Make sure the flow data is downloaded and added to the MATLAB path.

if isfile("sim_1_sim_2_merged_flow.mat") == 0 || isfile("sim_3_sim_4_merged_flow.mat") == 0

    sprintf("Make sure you have downloaded the relevant datasets and added them to the MATLAB path")
    return;

else

    load('sim_1_sim_2_merged_flow.mat');
    load('sim_3_sim_4_merged_flow.mat');

end

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

% Input max number of sensor selections the algorithm stops at (note the rule based selection only goes to 250)
max_number_sensors = 250;

% Load the sensor selection results, with the mutual information for the
% other sensor selections.

load("Data/Fig_5_MI_values.mat");
load("Data/flow_sensor_placement_250.mat")

t = tiledlayout('flow', 'TileSpacing', 'compact')

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
plot((1:max_number_sensors), MI_totalsum_flow, 'DisplayName', "Max total sum placement", 'color', 'cyan')
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
% --------------------------------------------------------

%%
% Figure 6 - Sensor placement plot - flow

clear all
close all

load("Data/Bellinge_sf.mat")
load("Data/flow_sensor_placement_250.mat")

junctions_table = struct2table(junctions);
outfalls_table = struct2table(outfalls);
storages_table = struct2table(storages);

junctions_table(:, 5:end) = [];
outfalls_table(:, 5:end) = [];
storages_table(:, 5:end) = [];

% Merge tables -> NOTE THE BELOW TABLE IS IN THE SAME ORDER AS THE
% TIMESERIES DATA IN TS_SIM_1_MASTER ETC.

merged_nodes_table = [junctions_table; outfalls_table; storages_table];
merged_nodes_table = table2struct(merged_nodes_table);

n = 1020;
list = optimal_sensor_selection_table(:, 5:end-1);

for i = 1:n

    percentage_vec(i) = length(find(list == i)) / height(optimal_sensor_selection_table);

end

percentage_vec = percentage_vec';
percent = percentage_vec;
percentage_vec = num2cell(percentage_vec);
[merged_nodes_table(:).PERCENTAGE] = percentage_vec{:};

figure
cmap = hot(250);
colorRange = makesymbolspec('Point', ...
    {'PERCENTAGE', [0.001, 1], 'MarkerEdgeColor', cmap, 'MarkerFaceColor', cmap, 'Marker', 'sq', 'MarkerEdgeColor', 'black', 'MarkerSize', 9}, ...
    {'PERCENTAGE', [0, 0.001], 'MarkerEdgeColor', cmap});
mapshow(merged_nodes_table, "SymbolSpec", colorRange)
colormap(cmap)
clim([0, 1])
hold on
mapshow(links, 'Color', 'blue')
axis off
box off
set(gcf, 'Position', [200, 200, 2000, 1000])
ylabel(colorbar(), '')

% --------------------------------------------------------

%%
% Figure 7 - Square error histogram

clear all
close all

% To run this script, you will need to download the file
% "flow_estimates_k_50.mat" from
% https://zenodo.org/doi/10.5281/zenodo.12517148. It is approximately
% 4.5gb.

if isfile("flow_estimates_k_50.mat") == 0

    promptMessage_est = sprintf('To proceed with plotting Figure 7, please download the dataset from https://zenodo.org/records/12517149');
    button = questdlg(promptMessage_est, 'Continue', 'Continue', 'Cancel', 'Continue');
    if strcmpi(button, 'Cancel')
        return; % Or break or continue
    end

    return;
else

    load("flow_estimates_k_50.mat");

end

error_matrix_GLM = validation_un_observed - prediction_matrix_GLM;
error_matrix_GRNNET = validation_un_observed - prediction_matrix_GRNNET;

error_sq_GLM = error_matrix_GLM.^2;
clear error_matrix_GLM;
error_sq_GRNNET = error_matrix_GRNNET.^2;
clear error_matrix_GRNNET;

t = tiledlayout('flow', 'TileSpacing', 'compact');

nexttile
histogram(log10(error_sq_GLM(:)), 'Normalization', 'pdf', 'EdgeColor', 'none')
xlabel('Log_{10} square error for each GLM data point estimate for unmonitored nodes (m^3s^{-1})^2', 'FontSize', 15)
ylabel('PDF value', 'FontSize', 15)
xlim([-17, -3])
nexttile
histogram(log10(error_sq_GRNNET(:)), 'Normalization', 'pdf', 'EdgeColor', 'none')
xlabel('Log_{10} square error for each GRNN data point estimate for unmonitored nodes (m^3s^{-1})^2', 'FontSize', 15)
ylabel('PDF value', 'FontSize', 15)
xlim([-17, -3])
set(gcf, 'Position', [200, 200, 715, 500])

% ----------------------------------------------------

%%
% Figure 8 - Flow sensor placement with k = 50.

clear all
close all

load('Data/Bellinge_sf.mat')
load("Data/flow_sensor_placement_250.mat")

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

n = 1020;
k_max = 50;

k_seq = (1:n)';
nodes_to_delete = zeros(1, n-k_max);
nodes_to_delete(1, 1:width(nodes_to_delete)) = intersect(k_seq, setdiff(k_seq, optimal_sensor_selection_table(k_max, 5:5+k_max-1)));

nodes_k_max = merged_nodes_table;
nodes_k_max(nodes_to_delete(1, 1:n-k_max)) = [];

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

% ----------------------------------------------------

%%
% Figure 9 - Estimate vs simulated for wet weather.

clear all
close all

% To run this script, you will need to download the file
% "flow_estimates_k_50.mat" from
% https://zenodo.org/doi/10.5281/zenodo.12517148. It is approximately
% 4.5gb.

if isfile("flow_estimates_k_50.mat") == 0

    promptMessage_est = sprintf('To proceed with plotting Figure 7, please download the dataset from https://zenodo.org/records/12517149 and add it to the MATLAB path');
    button = questdlg(promptMessage_est, 'Continue', 'Continue', 'Cancel', 'Continue');
    if strcmpi(button, 'Cancel')
        return; % Or break or continue
    end

    return;
else

    load("flow_estimates_k_50.mat");

end

load('sim_3_sim_4_merged_flow.mat');

storm_1_t1 = datetime("03-03-2019", 'Inputformat', 'dd-MM-yyyy');
storm_1_t2 = datetime("17-03-2019", 'Inputformat', 'dd-MM-yyyy');

storm_2_t1 = datetime("08-10-2019", 'Inputformat', 'dd-MM-yyyy');
storm_2_t2 = datetime("19-10-2019", 'Inputformat', 'dd-MM-yyyy');

storms = time_series_sim_3_and_4_master;
storms = renamevars(storms, "Date/time", "Date");

actual = validation_un_observed;
glm = prediction_matrix_GLM;
grnn = prediction_matrix_GRNNET;

% ------------- Section where we delete all unselected dates ------------- %

Dates_Keep_1 = isbetween(storms.Date, storm_1_t1, storm_1_t2);
Dates_Keep_2 = isbetween(storms.Date, storm_2_t1, storm_2_t2);

merge_dates_keep = Dates_Keep_1 + Dates_Keep_2;
merge_dates_delete = abs(1-merge_dates_keep);
merge_dates_delete = logical(merge_dates_delete);

% storms(merge_dates_delete,:) = [];
actual(merge_dates_delete, :) = [];
glm(merge_dates_delete, :) = [];
grnn(merge_dates_delete, :) = [];

% Reshape into a column vector for plotting
actual_reshape = reshape(actual, [], 1);
glm_reshape = reshape(glm, [], 1);
grnn_reshape = reshape(grnn, [], 1);

% Plot against y = x

yfunc = @(x) x;
xaxis = [0:0.01:1.5];
y_vals = yfunc(xaxis);

min_y_axis = min(min(glm_reshape), min(grnn_reshape));
max_y_axis = max(max(glm_reshape), max(grnn_reshape));

t = tiledlayout('flow', 'TileSpacing', 'compact')
nexttile
scatter(actual_reshape, glm_reshape)
hold on
plot(xaxis, y_vals, "Color", "black");
xlabel("Simulated flow values for unmonitored nodes (m^3s^{-1})", "Fontsize", 15)
ylabel("Estimated flow values \newline      GLM (m^3s^{-1})", "Fontsize", 15, "HorizontalAlignment", "center")
ylim([-0.5, 1.5]);
nexttile
scatter(actual_reshape, grnn_reshape)
hold on
plot(xaxis, y_vals, "Color", "black");
xlabel("Simulated flow values for unmonitored nodes  (m^3s^{-1})", "Fontsize", 15)
ylabel("Estimated flow values \newline      GRNN (m^3s^{-1})", "Fontsize", 15)
set(gcf, 'Position', [200, 200, 650, 500])
ylim([-0.5, 1.5]);

GLM_error = glm_reshape - actual_reshape;
GRNN_error = grnn_reshape - actual_reshape;
GLM_square_error_vec = GLM_error.^2;
GRNN_square_error_vec = GRNN_error.^2;
GLM_square_error = sum(GLM_square_error_vec, "all");
GRNN_square_error = sum(GRNN_square_error_vec, "all");
[glm_elements, ~] = size(glm_reshape);
[grnn_elements, ~] = size(grnn_reshape);

GLM_mean_square_error = GLM_square_error / glm_elements
GRNN_mean_square_error = GRNN_square_error / grnn_elements

% ----------------------------------------------------

%%
% Figure 10 - Water level sensor selection plot

clear all
close all

% Variance sigma^2 introduced from noise due to sensor
sigmas = 2.5 * 10^(-2);

% Input max number of sensor selections the algorithm stops at (note the rule based selection only goes to 250)
max_number_sensors = 250;

load("Data/Height_sensor_placement_250.mat");
load("Data/Bellinge_sf.mat")

junctions_table = struct2table(junctions);
outfalls_table = struct2table(outfalls);
storages_table = struct2table(storages);

junctions_table(:, 5:end) = [];
outfalls_table(:, 5:end) = [];
storages_table(:, 5:end) = [];

% Merge tables -> NOTE THE BELOW TABLE IS IN THE SAME ORDER AS THE
% TIMESERIES DATA IN TS_SIM_1_MASTER ETC.

merged_nodes_table = [junctions_table; outfalls_table; storages_table];
merged_nodes_table = table2struct(merged_nodes_table);

n = 1020;
list = optimal_sensor_selection_table_height(:, 5:end-1);

for i = 1:n

    percentage_vec_h(i) = length(find(list == i)) / height(optimal_sensor_selection_table_height);

end

percentage_vec_h = percentage_vec_h';
percent_h = percentage_vec_h;
percentage_vec_h = num2cell(percentage_vec_h);
[merged_nodes_table(:).PERCENTAGE_height] = percentage_vec_h{:};

figure
cmap = hot(250);
colorRange = makesymbolspec('Point', ...
    {'PERCENTAGE_height', [0.001, 1], 'MarkerEdgeColor', cmap, 'MarkerFaceColor', cmap, 'Marker', 'sq', 'MarkerEdgeColor', 'black', 'MarkerSize', 9}, ...
    {'PERCENTAGE_height', [0, 0.001], 'MarkerEdgeColor', cmap});
mapshow(merged_nodes_table, "SymbolSpec", colorRange)
colormap(cmap)
clim([0, 1])
hold on
mapshow(links, 'Color', 'blue')
axis off
box off
set(gcf, 'Position', [200, 200, 2000, 1000])
ylabel(colorbar(), '')

% ----------------------------------------------------

%%
% Figure 11 - Water level vs Flow (absolute difference)

clear all
close all

load("Data/flow_sensor_placement_250.mat")
load("Data/Height_sensor_placement_250.mat");
load("Data/Bellinge_sf.mat")

n = 1020;

junctions_table = struct2table(junctions);
outfalls_table = struct2table(outfalls);
storages_table = struct2table(storages);

junctions_table(:, 5:end) = [];
outfalls_table(:, 5:end) = [];
storages_table(:, 5:end) = [];

merged_nodes_table = [junctions_table; outfalls_table; storages_table];
merged_nodes_table = table2struct(merged_nodes_table);

list = optimal_sensor_selection_table(:, 5:end-1);

for i = 1:n

    percentage_vec(i) = length(find(list == i)) / height(optimal_sensor_selection_table);

end

percentage_vec = percentage_vec';
percent = percentage_vec;
percentage_vec = num2cell(percentage_vec);
[merged_nodes_table(:).PERCENTAGE] = percentage_vec{:};

clear list

list = optimal_sensor_selection_table_height(:, 5:end-1);

for i = 1:n

    percentage_vec_h(i) = length(find(list == i)) / height(optimal_sensor_selection_table_height);

end

percentage_vec_h = percentage_vec_h';
percent_h = percentage_vec_h;
percentage_vec_h = num2cell(percentage_vec_h);
[merged_nodes_table(:).PERCENTAGE_height] = percentage_vec_h{:};


absolute_percentage_error = abs(percent_h-percent);
absolute_percentage_error = num2cell(absolute_percentage_error);
[merged_nodes_table(:).PERCENTAGE_absolute_error] = absolute_percentage_error{:};

figure
cmap = hot(250);
colorRange = makesymbolspec('Point', ...
    {'PERCENTAGE_absolute_error', [0.001, 1], 'MarkerEdgeColor', cmap, 'MarkerFaceColor', cmap, 'Marker', 'sq', 'MarkerEdgeColor', 'black', 'MarkerSize', 9}, ...
    {'PERCENTAGE_absolute_error', [0, 0.001], 'MarkerEdgeColor', cmap});
mapshow(merged_nodes_table, "SymbolSpec", colorRange)
colormap(cmap)
clim([0, 1])
hold on
mapshow(links, 'Color', 'blue')
axis off
box off
set(gcf, 'Position', [200, 200, 2000, 1000])
ylabel(colorbar(), '')
