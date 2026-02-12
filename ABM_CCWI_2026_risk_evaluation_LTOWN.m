%% ABM_CCWI_2026 Risk evaluation
%% The user has the option to run the risk evaluation for both ABM and RAW scenarios by commenting/uncommenting the relevant sections

clear; clc; close all;
fprintf('=================================================================\n');
fprintf('   QMRA - ENTEROVIRUS - L-TOWN\n');
fprintf('=================================================================\n\n');
% %% Load data
%% BASELINE DETECTION | Contamination at 08:00

% load MSX_results_LTOWN_ABM_BASELINE_MODERATE_SOURCE.mat 
% load MSX_results_LTOWN_RAW_BASELINE_MODERATE_SOURCE.mat 
% load MSX_results_LTOWN_ABM_BASELINE_MODERATE_DETECTION.mat 
% load MSX_results_LTOWN_RAW_BASELINE_MODERATE_DETECTION.mat 

%% EARLY DETECTION | Contamination at 08:00

% load MSX_results_LTOWN_ABM_EARLY_MODERATE_SOURCE.mat 
% load MSX_results_LTOWN_RAW_EARLY_MODERATE_SOURCE.mat 
% load MSX_results_LTOWN_ABM_EARLY_MODERATE_DETECTION.mat 
% load MSX_results_LTOWN_RAW_EARLY_MODERATE_DETECTION.mat

%% BASELINE DETECTION | Contamination at 20:00

load MSX_results_LTOWN_ABM_BASELINE_MODERATE_SOURCE_INVERSED.mat 
% load MSX_results_LTOWN_RAW_BASELINE_MODERATE_SOURCE_INVERSED.mat 
% load MSX_results_LTOWN_ABM_BASELINE_MODERATE_DETECTION_INVERSED.mat 
% load MSX_results_LTOWN_RAW_BASELINE_MODERATE_DETECTION_INVERSED.mat 

%% EARLY DETECTION | Contamination at 20:00

% load MSX_results_LTOWN_ABM_EARLY_MODERATE_SOURCE_INVERSED.mat 
% load MSX_results_LTOWN_RAW_EARLY_MODERATE_SOURCE_INVERSED.mat 
% load MSX_results_LTOWN_ABM_EARLY_MODERATE_DETECTION_INVERSED.mat 
% load MSX_results_LTOWN_RAW_EARLY_MODERATE_DETECTION_INVERSED.mat



%%%% ================Load Stream demands for contamination at 08:00 ===================%%%%

% load Stream_demands_ABM_ltown_SOURCE.mat 
% load Stream_demands_ABM_ltown_DETECTION.mat
% load Stream_demands_ABM_ltown_EARLY_SOURCE.mat 
% load Stream_demands_ABM_ltown_EARLY_DETECTION.mat 

%%%% ================Load Stream demands for contamination at 20:00 ===================%%%%

load Stream_demands_ABM_ltown_INVERSED_SOURCE.mat 
% load Stream_demands_ABM_ltown_INVERSED_DETECTION.mat
% load Stream_demands_ABM_ltown_EARLY_INVERSED_SOURCE.mat 
% load Stream_demands_ABM_ltown_EARLY_INVERSED_DETECTION.mat 

P_conc = P_junctions;
N_people = round(People_per_node);

fprintf('Loaded:\n');
fprintf('  Timesteps: %d\n', size(P_conc, 1));
fprintf('  Nodes: %d\n', size(P_conc, 2));
fprintf('  Population: %d people\n', sum(N_people));

%% Parameters
DAILY_CONSUMPTION = 1.0;
MAX_SIP = 0.25;
TIMESTEP_HOURS = 5/60;
NUM_DAYS = 6;
r_entero = 0.014472;

total_timesteps = size(P_conc, 1);
TIMESTEPS_PER_DAY = floor(total_timesteps / NUM_DAYS);
has_extra_timestep = (total_timesteps > TIMESTEPS_PER_DAY * NUM_DAYS);

if has_extra_timestep
    fprintf('\n  NOTE: Data has %d timesteps (%d days × %d + %d extra)\n', ...
        total_timesteps, NUM_DAYS, TIMESTEPS_PER_DAY, total_timesteps - TIMESTEPS_PER_DAY * NUM_DAYS);
end

fprintf('\nQMRA Setup:\n');
fprintf('  Pathogen: Enterovirus\n');
fprintf('  Parameter r: %.6f\n', r_entero);
fprintf('  Daily consumption: %.2f L/person/day\n', DAILY_CONSUMPTION);
fprintf('  Simulation: %d days (%d timesteps)\n', NUM_DAYS, total_timesteps);

%% Unit Conversion
% fprintf('\n=== STEP 1: Unit Conversion ===\n');
% ABM
if max(Stream_faucet(:)) < 10
    Stream_faucet = Stream_faucet * 1000;
    fprintf('  ✓ Converted CMH to L/h: max = %.2f L/h\n', max(Stream_faucet(:)));
else
    fprintf('  Already in L/h: max = %.2f L/h\n', max(Stream_faucet(:)));
end 

% RAW
% if max(Stream_faucet_raw(:)) < 10
%     Stream_faucet_raw = Stream_faucet_raw * 1000;
%     fprintf('  ✓ Converted CMH to L/h: max = %.2f L/h\n', max(Stream_faucet_raw(:)));
% else
%     fprintf('  Already in L/h: max = %.2f L/h\n', max(Stream_faucet_raw(:)));
% end 

Stream_faucet_per_timestep = Stream_faucet * TIMESTEP_HOURS; % ABM
% Stream_faucet_raw_per_timestep = Stream_faucet_raw * TIMESTEP_HOURS; % RAW


%% Initialize
num_junctions = length(junction_names);
Volume = cell(1, num_junctions);
for m = 1:num_junctions
    Volume{m} = zeros(total_timesteps, N_people(m));
end

%% Process each day
fprintf('\n=== STEP 2-3: Daily Processing ===\n');
tic;

for day = 1:NUM_DAYS
    fprintf('\n--- Day %d / %d ---\n', day, NUM_DAYS);
    
    start_idx = (day - 1) * TIMESTEPS_PER_DAY + 1;
    end_idx = day * TIMESTEPS_PER_DAY;
    % % % 
    % ABM
    Stream_tap_day = Stream_faucet_per_timestep(start_idx:end_idx, :);
    people_tot = N_people;
    stream_tot = sum(Stream_tap_day, 1);
    fraction = people_tot ./ stream_tot;
    fraction(stream_tot < 0.001) = 0; 

    % RAW
    % Stream_tap_day_raw = Stream_faucet_raw_per_timestep(start_idx:end_idx, :);
    % people_tot = N_people;
    % stream_tot_raw = sum(Stream_tap_day_raw, 1);
    % fraction = people_tot ./ stream_tot_raw;
    % fraction(stream_tot_raw < 0.001) = 0; 

    Consumed_Stream = Stream_tap_day .* fraction; % ABM
    % Consumed_Stream_raw = Stream_tap_day_raw .* fraction; % RAW

    for m = 1:num_junctions
        if N_people(m) == 0
            continue;
        end
        
        Volume_day = zeros(TIMESTEPS_PER_DAY, N_people(m));
        consumed_daily = zeros(1, N_people(m));
        person_idx = 1;
        
        for t = 1:TIMESTEPS_PER_DAY
            remaining_water = Consumed_Stream(t, m); % ABM
            % remaining_water = Consumed_Stream_raw(t, m); % RAW
            
            while remaining_water > 0.00001
                person_remaining = DAILY_CONSUMPTION - consumed_daily(person_idx);
                
                if person_remaining > 0.001
                    drink = min([remaining_water, MAX_SIP, person_remaining]);
                    Volume_day(t, person_idx) = Volume_day(t, person_idx) + drink;
                    consumed_daily(person_idx) = consumed_daily(person_idx) + drink;
                    remaining_water = remaining_water - drink;
                end
                
                person_idx = person_idx + 1;
                if person_idx > N_people(m)
                    person_idx = 1;
                end
                
                if all(consumed_daily >= DAILY_CONSUMPTION - 0.001)
                    break;
                end
            end
        end
        
        Volume{m}(start_idx:end_idx, :) = Volume_day;
    end
    
    fprintf('  ✓ Day %d complete\n', day);
end

fprintf('\n✓ Volume allocation complete (%.1f s)\n', toc);

%% FORCE IDENTICAL VOLUME BEFORE DETECTION
% detection_ts = t_adv1;% DETECTION
detection_ts = t_adv2;% SOURCE

%FOR RAW RUN: Save Volume as reference
% save('Volume_RAW_reference_BASELINE_MODERATE_SOURCE_INVERSED.mat', 'Volume', 'detection_ts');
% fprintf('✓ Saved Volume as RAW reference\n');

% %FOR ABM RUN: Load and force-copy
if exist('Volume_RAW_reference_BASELINE_MODERATE_SOURCE_INVERSED.mat', 'file')
    RAW_data = load('Volume_RAW_reference_BASELINE_MODERATE_SOURCE_INVERSED.mat');
    Volume_RAW_ref = RAW_data.Volume;
    fprintf('Forcing Volume = RAW reference before detection...\n');
    for m = 1:num_junctions
        Volume{m}(1:detection_ts,:) = Volume_RAW_ref{m}(1:detection_ts,:);
    end
    fprintf('✓ Forced Volume identical before detection\n');
end
%% Dose Calculation
fprintf('\n=== STEP 4: Dose Calculation ===\n');
Dose = cell(1, num_junctions);
for m = 1:num_junctions
    Dose{m} = Volume{m} .* P_conc(:, m);
end
fprintf('✓ Dose calculated\n');

%% Infection Risk
fprintf('\n=== STEP 5: Infection Risk ===\n');
Risk_per_node = cell(1, num_junctions);
for m = 1:num_junctions
    Risk_per_node{m} = 1 - exp(-r_entero * Dose{m});
end
fprintf('✓ Risk calculated\n');

%% Per-Person and Time-Series Metrics
fprintf('\n=== STEP 6: Per-Person and Time-Series Metrics ===\n');

Total_risk_per_person = cell(1, num_junctions);
Total_risk_per_person_for_each_ts = cell(1, num_junctions);
Total_infections_per_timestep = zeros(total_timesteps, num_junctions);

for m = 1:num_junctions
    if N_people(m) == 0 || isempty(Risk_per_node{m})
        continue;
    end
    
    risk_matrix = Risk_per_node{m}';  % [people × timesteps]
    
    % Cumulative risk per person at each timestep
    Total_risk_per_person_for_each_ts{m} = zeros(size(risk_matrix));
    Total_risk_per_person_for_each_ts{m}(:, 1) = risk_matrix(:, 1);
    
    for t = 2:size(risk_matrix, 2)
        for p = 1:size(risk_matrix, 1)
            Total_risk_per_person_for_each_ts{m}(p, t) = ...
                1 - (1 - Total_risk_per_person_for_each_ts{m}(p, t-1)) * (1 - risk_matrix(p, t));
        end
    end
    
    % Final risk per person (at end of simulation)
    Total_risk_per_person{m} = Total_risk_per_person_for_each_ts{m}(:, end)';
    
    % Total infections at this node at each timestep
    for t = 1:total_timesteps
        Total_infections_per_timestep(t,m) = sum(Total_risk_per_person_for_each_ts{m}(:, t));
    end
end

% Aggregate metrics
Total_Infections_day = sum(cell2mat(cellfun(@(x) x(:)', Total_risk_per_person, 'UniformOutput', false)));
Total_risk_of_infection = (Total_Infections_day / sum(N_people)) * 100;

affected_nodes = find(max(P_conc, [], 1) > 1);
Downstream_pop = sum(N_people(affected_nodes));
Affected_Risk = (Total_Infections_day / Downstream_pop) * 100;

Total_infections_per_timestep_aggregated = sum(Total_infections_per_timestep, 2);

fprintf('✓ Metrics calculated\n');

%% Results
fprintf('\n=== RESULTS ===\n');
fprintf('  Total_Infections_day: %.1f people\n', Total_Infections_day);
fprintf('  Total_risk_of_infection: %.2f%%\n', Total_risk_of_infection);
fprintf('  Affected_Risk: %.2f%%\n', Affected_Risk);
fprintf('  Downstream_pop: %d people\n', Downstream_pop);

fprintf('\nTop 5 nodes by final infections:\n');
node_infections = cellfun(@sum, Total_risk_per_person);
[sorted_inf, sorted_idx] = sort(node_infections, 'descend');
for i = 1:min(5, num_junctions)
    idx = sorted_idx(i);
    if sorted_inf(i) > 0
        fprintf('  Node %s: %.1f / %d (%.1f%%)\n', ...
            junction_names{idx}, sorted_inf(i), N_people(idx), ...
            100*sorted_inf(i)/N_people(idx));
    end
end

fprintf('\nInfection progression by day:\n');
for day = 1:NUM_DAYS
    day_idx = day * TIMESTEPS_PER_DAY;
    if day_idx <= total_timesteps
        inf = Total_infections_per_timestep_aggregated(day_idx);
        fprintf('  Day %d: %.1f infections (%.2f%% attack rate)\n', ...
            day, inf, 100*inf/sum(N_people));
    end
end

% Save
save('QMRA_results_LTOWN_Enterovirus_FINAL_ABM_BASELINE_MODERATE_SOURCE_INVERSED.mat', ...
    'Volume', 'Dose', 'Risk_per_node', ...
    'Total_risk_per_person', 'Total_risk_per_person_for_each_ts', ...
    'Total_infections_per_timestep', 'Total_infections_per_timestep_aggregated', ...
    'Total_Infections_day', 'Total_risk_of_infection', 'Affected_Risk', 'Downstream_pop', ...
    'Time_hours', 'Time_days', 'junction_names', 'N_people', 'r_entero', '-v7.3'); % ABM % BASELINE SCENARIO
% % % 
% save('QMRA_results_LTOWN_Enterovirus_FINAL_RAW_BASELINE_MODERATE_SOURCE_INVERSED.mat', ...
%     'Volume', 'Dose', 'Risk_per_node', ...
%     'Total_risk_per_person', 'Total_risk_per_person_for_each_ts', ...
%     'Total_infections_per_timestep', 'Total_infections_per_timestep_aggregated', ...
%     'Total_Infections_day', 'Total_risk_of_infection', 'Affected_Risk', 'Downstream_pop', ...
%     'Time_hours', 'Time_days', 'junction_names', 'N_people', 'r_entero', '-v7.3'); % RAW % BASELINE SCENARIO


fprintf('\n✓ Saved: QMRA_results_LTOWN_Enterovirus_FINAL.mat\n');
fprintf('\n=================================================================\n');
fprintf('   COMPLETE                                                      \n');
fprintf('=================================================================\n');
