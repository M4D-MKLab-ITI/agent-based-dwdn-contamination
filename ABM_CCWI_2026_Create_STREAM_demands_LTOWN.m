% ABM modelling - Complete Working Version
%%%%%%%% Step 1: Create new Demands %%%%%%%%%
try
    d.unload
catch ERR
end
fclose all;clear class;clear all;clc;close all;

% ==========================================================================
%  SCENARIO CONFIGURATION - TUNABLE PARAMETERS FOR WHAT-IF ANALYSIS
% ==========================================================================

% ==========================================================================
% 1. TIMELINE PARAMETERS
% ==========================================================================

% Simulation duration
t_d = 6;                    % Total simulation days (6 days)

% Advisory trigger mode
% 'DETECTION' -> advisory at t_adv1
% 'SOURCE'    -> advisory at t_adv2
ADVISORY_MODE = 'SOURCE';   

% Contamination and response timeline
t_contam = 288 + 240;        % Contamination injection Day 2, 8:00= 96 - and Day 2, 20:00= 240
t_adv1 = t_contam + 288;     % Early Detection = 24, Baseline Detection = 288
t_adv2 = t_adv1 + 120;        % Source identification

switch ADVISORY_MODE
    case 'DETECTION'
        t_advisory = t_adv1;
    case 'SOURCE'
        t_advisory = t_adv2;
    otherwise
        error('Unknown ADVISORY_MODE');
end

t_all_clear = t_advisory + 288;   

% ==========================================================================
% 2. AWARENESS DYNAMICS PARAMETERS
% ==========================================================================

% Baseline awareness growth and social diffusion
alpha = 0.010;              % Baseline awareness growth rate (autonomous adoption)
beta = 0.12;               % Social diffusion strength (neighbor influence) 

% Initial awareness seeding (at detection)
seed_lvl_suspected = 0.20;  % Initial awareness in suspected contamination zone (20%)

% Awareness boost when source is found
bump_lvl_source = 0.30;     % Maximum awareness boost at source location (30%)
bump_decay_k = 3.0;         % Spatial decay rate of source bump (exponential falloff)

% Awareness decay after all-clear
ENABLE_AWARENESS_DECAY = true;  % Toggle: true = decay enabled, false = constant awareness
awareness_decay_rate = 0.02;   % Decay rate per timestep 

% ==========================================================================
% 3. ADOPTION CEILING PARAMETERS (Spatial Heterogeneity)
% ==========================================================================

% Maximum adoption rates vary with distance from contamination source
Amin = 0.70;                % Minimum adoption ceiling (far from source) 
Amax = 0.95;                % Maximum adoption ceiling (near source) - 
pA = 1.0;                   % Spatial gradient curvature (1.0 = linear, >1 = steeper)

% ==========================================================================
% 4. BEHAVIORAL ADOPTION PARAMETERS
% ==========================================================================

% Adherence probabilities (based on awareness threshold crossing)
PROB_FULL = 0.70;           % Probability of full adherence 
PROB_PARTIAL = 0.20;        % Probability of partial adherence 

% Water reduction factors
REDUCTION_FULL = 0.80;      % Full adherence
REDUCTION_PARTIAL = 0.45;   % Partial adherence
% ==========================================================================
% 5. VULNERABILITY PARAMETERS
% ==========================================================================

% Vulnerable populations 
vuln_min = 0.05;            % Minimum vulnerable population fraction (5%)
vuln_max = 0.40;            % Maximum vulnerable population fraction (40%)
vulnerability_boost = 0.40; % Extra adherence probability for vulnerable households (30%)

% Vulnerability zone centers (node IDs where vulnerable populations concentrate)
vulnerability_center_ids = {'n144', 'n359'};

% ==========================================================================
% 6. NETWORK TOPOLOGY PARAMETERS
% ==========================================================================

% Contamination source
true_src_ids = {'n129'};    % True contamination source node(s)

% Suspected contamination zone (if empty, will be auto-generated with spatial error)
suspected_ids = {'n69', 'n129', 'n498', 'n520', 'n519', 'n140', 'n144', 'n435', 'n498', 'n117', ...
                 'n415', 'n496', 'n497', 'n130', 'n521', 'n141', 'n133', 'n528', 'n529', 'n70', ...
                 'n67', 'n66', 'n413', 'n412', 'n56', 'n52', 'n407', 'n58', 'n406', 'n51', ...
                 'n128', 'n517', 'n518', 'n522', 'n138'};         % Leave empty for automatic clustered zone generation
error_radius = 0.05;        % Spatial uncertainty: fraction of network diameter (20%)
n_suspected = 35;           % Number of nodes in suspected zone

% ==========================================================================
% 7. POPULATION AND CONSUMPTION PARAMETERS
% ==========================================================================

% Household characteristics
Vpd = 150;                  % Volume per person per day (liters)
Population_unc = 0.1;       % Population uncertainty (±10%)

% Random seed control
rng_seed_main = 42;         % Main RNG seed for zone generation
rng_seed_vuln = 43;         % RNG seed for vulnerability zones
rng_seed_population = 12345;% RNG seed for household generation

% ==========================================================================
%  END OF CONFIGURATION SECTION
% ==========================================================================

LOG_HH_ABM = true;
KEEP_TS = true;
addpath(genpath(pwd));
disp('EPANET-MATLAB Toolkit Paths Loaded.');
inpname = 'L-TOWN.inp';
d = epanet(inpname);
t_d = 6;
duration_sec = t_d * 24 * 60 * 60;
d.setTimeSimulationDuration(duration_sec);
d.setTimeReportingStep(300)
d.setTimeHydraulicStep(300)
d.setTimePatternStep(300)
hydraulics = d.getComputedTimeSeries;
Demand = hydraulics.Demand(:, d.getNodeJunctionIndex);
Demand = Demand .* 1000;
Node = d.getNodeNameID;
junctionIndex = d.getNodeJunctionIndex;
Node = Node(junctionIndex);
Node_fields = Node;
K= 288;
Dt= double(d.getTimeHydraulicStep)/3600;
V_sc=cumsum(Demand(1:K,:),1)*Dt;
People_per_node= V_sc(end,:)/Vpd;
dt_min = 5;
steps_per_day = (24*60)/dt_min;
T = t_d*steps_per_day + 1;
L = T;
I = numel(junctionIndex);
Node = Node(:);
if numel(Node) ~= I
    Node = Node(1:I);
end
linkNodes = d.getLinkNodesIndex;
A_adj = sparse(I, I);
for e = 1:size(linkNodes,1)
    a_all = linkNodes(e,1);
    b_all = linkNodes(e,2);
    ia = find(junctionIndex==a_all, 1);
    ib = find(junctionIndex==b_all, 1);
    if ~isempty(ia) && ~isempty(ib)
        A_adj(ia, ib) = 1;
        A_adj(ib, ia) = 1;
    end
end
deg = sum(A_adj,2); deg(deg==0) = 1;
XY = d.getNodeCoordinates();
if iscell(XY)
    if numel(XY) < 2
        error('getNodeCoordinates returned a cell with <2 elements.');
    end
    X_all = XY{1};
    Y_all = XY{2};
elseif isnumeric(XY)
    if size(XY,2) < 2
        error('getNodeCoordinates returned numeric but not N x 2.');
    end
    X_all = XY(:,1);
    Y_all = XY(:,2);
elseif isstruct(XY)
    if ~isfield(XY,'X') || ~isfield(XY,'Y')
        error('getNodeCoordinates returned struct without X/Y fields.');
    end
    X_all = XY.X;
    Y_all = XY.Y;
else
    error('Unknown format from getNodeCoordinates().');
end
X_all = X_all(:);
Y_all = Y_all(:);
if any(junctionIndex > numel(X_all))
    error('junctionIndex exceeds node count reported by getNodeCoordinates().');
end
Xj = X_all(junctionIndex);
Yj = Y_all(junctionIndex);
t_day = @(d_) min(max(round((d_*24*60)/dt_min) + 1, 1), T);
true_src_idx = [];
for s = 1:numel(true_src_ids)
    k = find(strcmp(Node, true_src_ids{s}), 1);
    if ~isempty(k), true_src_idx(end+1) = k; end
end
if isempty(true_src_idx)
    warning('True source IDs not found among junctions. Using random fallback.');
    rng(rng_seed_main);
    true_src_idx = randi(I);
end
suspected_idx = [];
for s = 1:numel(suspected_ids)
    k = find(strcmp(Node, suspected_ids{s}), 1);
    if ~isempty(k), suspected_idx(end+1) = k; end
end
if isempty(suspected_idx)
    rng(42);
    src_x = Xj(true_src_idx(1));
    src_y = Yj(true_src_idx(1));
    network_diameter = sqrt((max(Xj)-min(Xj))^2 + (max(Yj)-min(Yj))^2);
    error_radius = error_radius * network_diameter;
    suspected_center_x = src_x + error_radius * (2*rand() - 1);
    suspected_center_y = src_y + error_radius * (2*rand() - 1);
    D_suspected = hypot(Xj - suspected_center_x, Yj - suspected_center_y);
    [~, sorted_idx] = sort(D_suspected);
    n_suspected = min(n_suspected, I);
    suspected_idx = sorted_idx(1:n_suspected);
    fprintf('Generated clustered suspected zone:\n');
    fprintf('  Center: (%.1f, %.1f)\n', suspected_center_x, suspected_center_y);
    fprintf('  True source: (%.1f, %.1f)\n', src_x, src_y);
    fprintf('  Spatial error: %.1f units (%.1f%%%% of network)\n', ...
            hypot(suspected_center_x - src_x, suspected_center_y - src_y), ...
            100 * hypot(suspected_center_x - src_x, suspected_center_y - src_y) / network_diameter);
    fprintf('  Zone contains %d nodes\n', n_suspected);
    fprintf('  Max radius: %.1f units\n', max(D_suspected(suspected_idx)));
end
vulnerability_idx = [];
for s = 1:numel(vulnerability_center_ids)
    k = find(strcmp(Node, vulnerability_center_ids{s}), 1);
    if ~isempty(k), vulnerability_idx(end+1) = k; end
end
if isempty(vulnerability_idx)
    rng(rng_seed_vuln);
    n_vuln_zones = 2;
    vulnerability_idx = randperm(I, min(n_vuln_zones, I));
    fprintf('Generated %d vulnerability zone centers randomly\n', n_vuln_zones);
else
    fprintf('Using manually specified vulnerability centers: %s, %s\n', vulnerability_center_ids{:});
end
D_vuln = inf(I, 1);
for v = 1:numel(vulnerability_idx)
    xv = Xj(vulnerability_idx(v));
    yv = Yj(vulnerability_idx(v));
    D_vuln = min(D_vuln, hypot(Xj - xv, Yj - yv));
end
D_vuln_max = max(D_vuln);
if D_vuln_max <= 0, D_vuln_max = 1; end
D_vuln_norm = D_vuln / D_vuln_max;
vuln_prob_node = vuln_max - (vuln_max - vuln_min) * D_vuln_norm;
fprintf('Vulnerability zones defined:\n');
fprintf('  Centers: %d nodes\n', numel(vulnerability_idx));
fprintf('  Vulnerable population ranges: %.1f%%%% to %.1f%%%%\n', ...
        100*min(vuln_prob_node), 100*max(vuln_prob_node));
D = inf(I,1);
for s = 1:numel(true_src_idx)
    xs = Xj(true_src_idx(s)); ys = Yj(true_src_idx(s));
    D = min(D, hypot(Xj - xs, Yj - ys));
end
Dmax = max(D); if Dmax <= 0, Dmax = 1; end
Dnorm = D / Dmax;
Amax_nodes = Amin + (Amax - Amin) * (1 - Dnorm).^pA;
% === AWARENESS DECAY PARAMETERS ===
A_now = zeros(I,1);
A_mat = zeros(T, I);
fprintf('\n=== Awareness Dynamics ===\n');
fprintf('Detection (t_advisory): %d (Day %.1f, 8:00 AM)\n', t_advisory, t_advisory/288);
fprintf('Source found (t_advisory): %d (Day %.1f, 12:30 PM)\n', t_advisory, t_advisory/288);
if ENABLE_AWARENESS_DECAY
    fprintf('Awareness decay: ENABLED (%.1f%% per timestep, starts at t=%d)\n', ...
        100*awareness_decay_rate, t_all_clear);
else
    fprintf('Awareness decay: DISABLED (constant awareness after adoption)\n');
end

for t = 1:T
    % Initial awareness seeding at detection
    if t == t_advisory
        A_now(suspected_idx) = max(A_now(suspected_idx), seed_lvl_suspected);
    end
    
    % Social diffusion (awareness spreads)
    if t >= t_advisory && t < t_all_clear
        neigh_mean = (A_adj * A_now) ./ deg;
        dA = (alpha + beta * neigh_mean) .* (1 - A_now);
        A_now = min(Amax_nodes, A_now + dA);
    end
    
    % === NEW: Awareness decay after all-clear ===
    if ENABLE_AWARENESS_DECAY && t >= t_all_clear
        % Exponential decay - awareness fades over time
        A_now = A_now * (1 - awareness_decay_rate);
    end
    
    A_mat(t, :) = A_now.';
end
fprintf('Max final awareness: %.3f\n', max(A_mat(end,:)));
load database.mat
People_per_node_rnd= round(People_per_node);
num_max_households = ceil(sum(People_per_node_rnd) / 1);
rng(rng_seed_population);
all_household_seeds = randi(1e9, num_max_households, 1);
all_household_sizes = randi(5, num_max_households, 1);
all_vuln_draws = rand(num_max_households, 1);
all_threshold_draws = rand(num_max_households, 1);
all_adherence_draws = rand(num_max_households, 1);
global_seed_counter = 0;
all_population_rand = rand(length(People_per_node_rnd), 1);
output = struct;
adoption = struct;
output_raw.output = struct();

for i=1:length(People_per_node_rnd)
    disp(['Simulating Node ',num2str(i),' of ',num2str(length(People_per_node_rnd))])
    Population=People_per_node_rnd(i);
    Population_l=Population-Population_unc*Population;
    Population_u=Population+Population_unc*Population;
    Population=Population_l+all_population_rand(i).*(Population_u-Population_l);
    Population = round(Population);
    Ai = A_mat(:, i);
    hh_thresholds = [];
    hh_t_adopt5 = [];
    hh_adherence_category = {};
    hh_reduction_factor = [];
    currHH = 0;
    hh_is_vulnerable = [];
    StToilet_raw = zeros(L,1);
    StShower_raw = zeros(L,1);
    StFaucet_raw = zeros(L,1);
    StClothesWasher_raw = zeros(L,1);
    StDishwasher_raw = zeros(L,1);
    TOTAL_raw = zeros(L,1);
    if LOG_HH_ABM
        hh_abm_log = [];
    end
    home=0;
    while Population>0
        home=home+1;
        param.HHsize = all_household_sizes(global_seed_counter + 1);
        Population=Population-param.HHsize;
        vuln_draw = all_vuln_draws(global_seed_counter + 1);
        if vuln_draw < vuln_prob_node(i)
            is_vulnerable = true;
        else
            is_vulnerable = false;
        end
        currHH = currHH + 1;
        hh_is_vulnerable(currHH,1) = is_vulnerable;
        u = all_threshold_draws(global_seed_counter + 1);
        t_hit = find(Ai >= u, 1, 'first'); if isempty(t_hit), t_hit = inf; end
        hh_thresholds(currHH,1) = u;
        hh_t_adopt5(currHH,1) = t_hit;
        if isinf(t_hit)
            hh_adherence_category{currHH,1} = 'none';
            hh_reduction_factor(currHH,1) = 0.0;
        else
            adherence_draw = all_adherence_draws(global_seed_counter + 1);
            spatial_bias = (Amax_nodes(i) - Amin) / (Amax - Amin);
            base_tier1 = PROB_FULL;
            base_tier2 = PROB_FULL + PROB_PARTIAL;
            if is_vulnerable
                vulnerability_boost = 0.30;
            else
                vulnerability_boost = 0.0;
            end
            bias_strength = 0.15;
            tier1_cutoff = base_tier1 + bias_strength * spatial_bias + vulnerability_boost;
            tier1_cutoff = min(max(tier1_cutoff, 0), 1);
            tier2_cutoff = base_tier2 + bias_strength * spatial_bias + vulnerability_boost * 0.5;
            tier2_cutoff = min(max(tier2_cutoff, 0), 1);
            if adherence_draw < tier1_cutoff
                hh_adherence_category{currHH,1} = 'full';
                hh_reduction_factor(currHH,1) = REDUCTION_FULL;
            elseif adherence_draw < tier2_cutoff
                hh_adherence_category{currHH,1} = 'partial';
                hh_reduction_factor(currHH,1) = REDUCTION_PARTIAL;
            else
                hh_adherence_category{currHH,1} = 'none';
                hh_reduction_factor(currHH,1) = 0.0;
            end
        end
        param.currHH = currHH;
        param.reduction_factor = hh_reduction_factor(currHH);
        param.t_adopt5_vec = hh_t_adopt5;
        param.appliances.StToilet = 1; param.appliances.HEToilet = 0;
        param.appliances.StShower = 1; param.appliances.HEShower = 0;
        param.appliances.StFaucet = 1; param.appliances.HEFaucet = 0;
        param.appliances.StClothesWasher = 1; param.appliances.HEClothesWasher = 0;
        param.appliances.StDishwasher = 1; param.appliances.HEDishwasher = 0;
        param.appliances.StBathtub = 1; param.appliances.HEBathtub = 0;
        param.H = 6;
        param.ts = 30;
        temp=checkInput(param);
        global_seed_counter = global_seed_counter + 1;
        seed_hh = all_household_seeds(global_seed_counter);
        rng(seed_hh);
        outputTrajectory = initializeTrajectories(param);
        appNames = fieldnames(outputTrajectory);
        for app=appNames'
            outputTrajectory.(char(app))=zeros(1,length(outputTrajectory.TOTAL)+30);
        end
        outputTrajectory = generateConsumptionEvents(outputTrajectory, param, database);
        outputTrajectory = sumToTotal(outputTrajectory);
        outputTrajectory = aggregateSamplingResolution(outputTrajectory, param);
        StToilet_raw = outputTrajectory.StToilet(:) + StToilet_raw;
        StShower_raw = outputTrajectory.StShower(:) + StShower_raw;
        StFaucet_raw = outputTrajectory.StFaucet(:) + StFaucet_raw;
        StClothesWasher_raw = outputTrajectory.StClothesWasher(:) + StClothesWasher_raw;
        StDishwasher_raw = outputTrajectory.StDishwasher(:) + StDishwasher_raw;
        TOTAL_raw = outputTrajectory.TOTAL(:) + TOTAL_raw;
    end
    node_field = Node_fields{i};
    output_raw.output.(node_field).StToilet = StToilet_raw;
    output_raw.output.(node_field).StShower = StShower_raw;
    output_raw.output.(node_field).StFaucet = StFaucet_raw;
    output_raw.output.(node_field).StClothesWasher = StClothesWasher_raw;
    output_raw.output.(node_field).StDishwasher = StDishwasher_raw;
    output_raw.output.(node_field).TOTAL = TOTAL_raw;
    adoption.hh_t_adopt5_by_node.(node_field) = hh_t_adopt5;
    adoption.hh_adherence_category_by_node.(node_field) = hh_adherence_category;
    adoption.hh_reduction_factor_by_node.(node_field) = hh_reduction_factor;
    adoption.hh_is_vulnerable_by_node.(node_field) = hh_is_vulnerable;
    adoption.hh_thresholds_by_node.(node_field) = hh_thresholds;  
end

fprintf('\nCreating ABM from RAW data (WITH DYNAMIC DECAY)...\n');
for i = 1:length(People_per_node_rnd)
    node_field = Node_fields{i};
    
    % Copy RAW data
    output.output.(node_field).StToilet = output_raw.output.(node_field).StToilet(:);
    output.output.(node_field).StShower = output_raw.output.(node_field).StShower(:);
    output.output.(node_field).StFaucet = output_raw.output.(node_field).StFaucet(:);
    output.output.(node_field).StClothesWasher = output_raw.output.(node_field).StClothesWasher(:);
    output.output.(node_field).StDishwasher = output_raw.output.(node_field).StDishwasher(:);
    output.output.(node_field).TOTAL = output_raw.output.(node_field).TOTAL(:);
    
    % Get awareness time series for this node
    Ai = A_mat(:, i);
    num_hh = length(adoption.hh_t_adopt5_by_node.(node_field));
    
    % Create time-varying reduction based on current awareness
    reduction_vector = zeros(L, 1);
    
    for t = 1:L
        active_reduction = 0;
        active_hh = 0;
        
        for hh = 1:num_hh
            threshold = adoption.hh_thresholds_by_node.(node_field)(hh);
            rf = adoption.hh_reduction_factor_by_node.(node_field)(hh);
            
            % Check if awareness is CURRENTLY above threshold
            if Ai(t) >= threshold && rf > 0
                active_reduction = active_reduction + rf;
                active_hh = active_hh + 1;
            end
        end
        
        if active_hh > 0
            reduction_vector(t) = active_reduction / num_hh;
        end
    end
    
    % Apply time-varying reduction
    output.output.(node_field).StFaucet = ...
        output.output.(node_field).StFaucet .* (1 - reduction_vector);
    output.output.(node_field).StShower = ...
        output.output.(node_field).StShower .* (1 - reduction_vector);
    
    % Recalculate TOTAL
    output.output.(node_field).TOTAL = ...
        output.output.(node_field).StToilet + ...
        output.output.(node_field).StShower + ...
        output.output.(node_field).StFaucet + ...
        output.output.(node_field).StClothesWasher + ...
        output.output.(node_field).StDishwasher;
end
expected_length = L;
J = numel(Node_fields);
Stream_demand_tot_raw = zeros(expected_length, J);
Stream_demand_Faucet_raw = zeros(expected_length, J);
Stream_demand_tot = zeros(expected_length, J);
Stream_demand_Faucet = zeros(expected_length, J);
for j = 1:J
    node_field = Node_fields{j};
    total_series_raw = output_raw.output.(node_field).TOTAL(:);
    faucet_series_raw = output_raw.output.(node_field).StFaucet(:);
    if isempty(total_series_raw), total_series_raw = zeros(expected_length,1); end
    if isempty(faucet_series_raw), faucet_series_raw = zeros(expected_length,1); end
    Stream_demand_tot_raw(:, j) = total_series_raw;
    Stream_demand_Faucet_raw(:, j) = faucet_series_raw;
    total_series = output.output.(node_field).TOTAL(:);
    faucet_series = output.output.(node_field).StFaucet(:);
    if isempty(total_series), total_series = zeros(expected_length,1); end
    if isempty(faucet_series), faucet_series = zeros(expected_length,1); end
    Stream_demand_tot(:, j) = total_series;
    Stream_demand_Faucet(:, j) = faucet_series;
end
Stream_raw = (Stream_demand_tot_raw .* 12) / 1000;
Stream_faucet_raw = (Stream_demand_Faucet_raw .* 12) / 1000;
Stream = (Stream_demand_tot .* 12) / 1000;
Stream_faucet = (Stream_demand_Faucet .* 12) / 1000;

% FORCE ABM = RAW before detection time
fprintf('\nForcing ABM = RAW for timesteps 1:t_advisory...\n');
Stream(1:t_advisory, :) = Stream_raw(1:t_advisory, :); 
Stream_faucet(1:t_advisory, :) = Stream_faucet_raw(1:t_advisory, :);

% Verify
fprintf('Verification: sum(Stream(1:t_adv1,:)) = %.2f, sum(Stream_raw(1:t_adv1,:)) = %.2f\n', ...
    sum(sum(Stream(1:t_advisory,:))), sum(sum(Stream_raw(1:t_advisory,:))));
fprintf('Difference: %.10f (should be exactly 0)\n', ...
    abs(sum(sum(Stream(1:t_advisory,:))) - sum(sum(Stream_raw(1:t_advisory,:)))));

for j = 1:J
    node_field = Node_fields{j};
    output.output.(node_field).TOTAL(1:t_advisory) = output_raw.output.(node_field).TOTAL(1:t_advisory);
    output.output.(node_field).StFaucet(1:t_advisory) = output_raw.output.(node_field).StFaucet(1:t_advisory);
    output.output.(node_field).StShower(1:t_advisory) = output_raw.output.(node_field).StShower(1:t_advisory);
end
fprintf('Also forced output struct to match before t=t_advisory.\n');

outfile = sprintf('Stream_demands_ABM_ltown_INVERSED_%s.mat', ADVISORY_MODE); % For contamination at 20:00 the name should be: Stream_demands_ABM_ltown_INVERSED_%s.mat | For contamination at 08:00 the name should be: Stream_demands_ABM_ltown_%s.mat

save(outfile, ...
    'output','output_raw','adoption','A_mat','Amax_nodes','t_adv1','t_adv2','t_contam', ...
    'Stream','Stream_faucet','Stream_raw','Stream_faucet_raw', ...
    'People_per_node','Node','junctionIndex','dt_min','T','I', ...
    'vuln_prob_node','vulnerability_idx','D_vuln_norm', ...
    'Xj', 'Yj', 'A_adj', 'deg', 'D', 'Dnorm','scenario_table');  

fprintf('✓ Saved: %s\n', outfile);


