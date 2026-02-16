%% Paper results: Section 3.1- Fig 2
clear; clc;

f_det = 'Stream_demands_ABM_ltown_DETECTION.mat';
f_src = 'Stream_demands_ABM_ltown_SOURCE.mat';

Sdet = load(f_det);
Ssrc = load(f_src);

Stream_det = Sdet.Stream;        % ABM at detection
Stream_src = Ssrc.Stream;        % ABM at source ID
Stream_raw = Sdet.Stream_raw;    % RAW baseline

dt_sec = 300;
nSteps = size(Stream_raw,1);
t_days = ((0:nSteps-1)*dt_sec)/86400;

Q_raw = sum(Stream_raw, 2); 
Q_det = sum(Stream_det, 2);
Q_src = sum(Stream_src, 2);

% --- Event timings (steps/day = 288) ---
steps_per_day = 86400/dt_sec;     % 288
t_contam = 288 + 96;              % Day 2 @ 08:00
t_adv1  = t_contam + 288;         % moderate detection: +24h
t_adv2  = t_adv1 + 120;           % source identified: +10h
detect_day = t_adv1 / steps_per_day;
source_day = t_adv2 / steps_per_day;

% styling
lw_raw = 2.8; lw_abm = 2.2; lw_vline = 2.0;

fig = figure('Color','w'); fig.Units = 'pixels'; fig.Position = [120 120 1200 520];
hold on; box on; grid on;

plot(t_days, Q_raw, 'k-', 'LineWidth', lw_raw, 'DisplayName','RAW (baseline)');
plot(t_days, Q_det, 'LineWidth', lw_abm, 'DisplayName','ABM: Advisory at detection');
plot(t_days, Q_src, 'LineWidth', lw_abm, 'DisplayName','ABM: Advisory at source ID');

xlabel('Time (days)','FontSize', 18);
ylabel('Total network demand (m^3/h)','FontSize', 18);
set(gca, 'FontSize', 18); 
xline(detect_day, 'k:', 'LineWidth', lw_vline, 'HandleVisibility','off');
xline(source_day, 'k:', 'LineWidth', lw_vline, 'HandleVisibility','off');
ax = gca;
yl = ylim(ax);
y_text_detect = yl(1) + 0.11 * range(yl);   
y_text_source = yl(1) + 0.20 * range(yl);

text(detect_day, y_text_detect, 'Detection', ...
    'Rotation', 90, 'VerticalAlignment','bottom', 'HorizontalAlignment','center', ...
    'FontSize', 18, 'BackgroundColor','w', 'Margin', 1);

text(source_day, y_text_source, 'Source identified', ...
    'Rotation', 90, 'VerticalAlignment','bottom', 'HorizontalAlignment','center', ...
    'FontSize', 18, 'BackgroundColor','w', 'Margin', 1);

% --- Legend ---
lgd = legend('show');
lgd.Box = 'off';
lgd.FontSize = 18;
lgd.Units = 'normalized';
lgd.Position = [0.68 0.125 0.25 0.12];

%% Section 3.3- Fig 3
%% Temporal Evolution of Infections: Linear Split (08:00 vs 20:00)
clear; clc; close all;

dt_min = 5;  
candI = {'Total_infections_per_timestep','total_infections_per_timestep', ...
         'Infections_per_timestep','expected_infections_ts'};

%% --- DATA DEFINITION ---
early_cases = {
    'Adv@Detection', 'QMRA_results_LTOWN_Enterovirus_FINAL_ABM_EARLY_MODERATE_DETECTION.mat',          'Stream_demands_ABM_ltown_EARLY_DETECTION.mat',          1, [0 0.447 0.741], '-', 'Adv@Detection', 'Morning';
    'Adv@SourceID',  'QMRA_results_LTOWN_Enterovirus_FINAL_ABM_EARLY_MODERATE_SOURCE.mat',             'Stream_demands_ABM_ltown_EARLY_SOURCE.mat',             2, [0.85 0.325 0.098], '-', 'Adv@SourceID', 'Morning';
    'Adv@Detection', 'QMRA_results_LTOWN_Enterovirus_FINAL_ABM_EARLY_MODERATE_DETECTION_INVERSED.mat', 'Stream_demands_ABM_ltown_EARLY_INVERSED_DETECTION.mat', 1, [0 0.447 0.741], ':', 'Adv@Detection', 'Inversed';
    'Adv@SourceID',  'QMRA_results_LTOWN_Enterovirus_FINAL_ABM_EARLY_MODERATE_SOURCE_INVERSED.mat',    'Stream_demands_ABM_ltown_EARLY_INVERSED_SOURCE.mat',    2, [0.85 0.325 0.098], ':', 'Adv@SourceID', 'Inversed';
};

early_raw = {
    'RAW 08:00', 'QMRA_results_LTOWN_Enterovirus_FINAL_RAW_EARLY_MODERATE_DETECTION.mat', 'Morning';
    'RAW 20:00', 'QMRA_results_LTOWN_Enterovirus_FINAL_RAW_EARLY_MODERATE_DETECTION_INVERSED.mat', 'Inversed';
};

mod_cases = {
    'Adv@Detection', 'QMRA_results_LTOWN_Enterovirus_FINAL_ABM_BASELINE_MODERATE_DETECTION.mat',          'Stream_demands_ABM_ltown_DETECTION.mat',          1, [0 0.447 0.741], '-', 'Adv@Detection', 'Morning';
    'Adv@SourceID',  'QMRA_results_LTOWN_Enterovirus_FINAL_ABM_BASELINE_MODERATE_SOURCE.mat',             'Stream_demands_ABM_ltown_SOURCE.mat',             2, [0.85 0.325 0.098], '-', 'Adv@SourceID', 'Morning';
    'Adv@Detection', 'QMRA_results_LTOWN_Enterovirus_FINAL_ABM_BASELINE_MODERATE_DETECTION_INVERSED.mat', 'Stream_demands_ABM_ltown_INVERSED_DETECTION.mat', 1, [0 0.447 0.741], ':', 'Adv@Detection', 'Inversed';
    'Adv@SourceID',  'QMRA_results_LTOWN_Enterovirus_FINAL_ABM_BASELINE_MODERATE_SOURCE_INVERSED.mat',    'Stream_demands_ABM_ltown_INVERSED_SOURCE.mat',    2, [0.85 0.325 0.098], ':', 'Adv@SourceID', 'Inversed';
};

mod_raw = {
    'RAW 08:00', 'QMRA_results_LTOWN_Enterovirus_FINAL_RAW_BASELINE_MODERATE_DETECTION.mat', 'Morning';
    'RAW 20:00', 'QMRA_results_LTOWN_Enterovirus_FINAL_RAW_BASELINE_MODERATE_DETECTION_INVERSED.mat', 'Inversed';
};

%% --- EXECUTE PLOTTING ---
plot_split_linear(early_raw, early_cases, dt_min, candI, 'Early Scenario', 'Split_Early_Linear.png');
plot_split_linear(mod_raw, mod_cases, dt_min, candI, 'Moderate Scenario', 'Split_Moderate_Linear.png');

%% ========================================================================
function plot_split_linear(rawFiles, cases, dt_min, candI, mainTitle, saveName)
    fig = figure('Color','w','Position', [100 100 1600 650]); 
    
    % Change to 1 row, 2 columns
    t = tiledlayout(1, 2, 'TileSpacing', 'Compact', 'Padding', 'Compact');
    title(t, mainTitle, 'FontSize', 24, 'FontWeight', 'bold');
    
    types = {'Morning', 'Inversed'};
    subTitles = {'Contamination at 08:00', 'Contamination at 20:00'};
    
    for s = 1:2
        nexttile; hold on; box on; grid on;
        currentType = types{s};
        
        % 1. Plot RAW Baseline
        idxRaw = find(strcmp(rawFiles(:,3), currentType));
        [I_raw, t_days] = load_and_sum(rawFiles{idxRaw,2}, candI, dt_min);
        plot(t_days, I_raw, 'k', 'LineWidth', 2.5, 'DisplayName', 'No Advisory (RAW)');
        
        % 2. Plot ABM Curves
        idxCases = find(strcmp(cases(:,8), currentType))';
        for i = idxCases
            [I_ts, t_days] = load_and_sum(cases{i,2}, candI, dt_min);
            plot(t_days, I_ts, 'LineWidth', 2.2, 'Color', cases{i,5}, ...
                'LineStyle', cases{i,6}, 'DisplayName', cases{i,1});
        end
        
        % Adjust Y-limit
        yl = ylim;
        ylim([0, yl(2)*1.25]); 
        yl = ylim; 
        
        % 3. Vertical Lines and Labels
        for i = idxCases
            Sd = load(cases{i,3});
            fld = ['t_adv', num2str(cases{i,4})];
            t_val = Sd.(fld);
            x_adv = (double(t_val(1)) - 1) * (dt_min/1440);
            
            line([x_adv x_adv], yl, 'Color', cases{i,5}, 'LineWidth', 2.0, ...
                'LineStyle', cases{i,6}, 'HandleVisibility', 'off');
            
            y_pos = yl(1) + 0.02 * diff(yl); 
            
            text(x_adv, y_pos, ['  ' cases{i,7}], ...
                'Rotation', 90, 'FontWeight', 'bold', 'FontSize', 18, ...
                'Color', cases{i,5}, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
        end
        
        ylabel('Infections');
        xlabel('Time (days)'); 
        title(subTitles{s}, 'FontSize', 20);
        legend('Location', 'northeast', 'FontSize', 18);
        set(gca, 'FontSize', 18);
    end
    
    saveas(fig, saveName);
end


function [I_out, t_out] = load_and_sum(filename, candI, dt_min)
    S = load(filename);
    data = [];
    for k = 1:numel(candI)
        if isfield(S, candI{k}), data = S.(candI{k}); break; end
    end
    if ~isvector(data)
        [r, c] = size(data);
        if r > c, I_out = sum(data, 2, 'omitnan'); else, I_out = sum(data, 1, 'omitnan'); end
    else
        I_out = data;
    end
    I_out = I_out(:);
    t_out = (0:length(I_out)-1)' * (dt_min/1440);
end

%% Additional info
% Total water demand reduction per water end-use

clear; clc;

% --- Load scenario file (moderate detection, contam at 08:00) ---
S = load('Stream_demands_ABM_ltown_DETECTION.mat');  % change if needed

output_raw = S.output_raw;
output     = S.output;

% --- Node list ---
if isfield(S,'Node_fields')
    nodes = S.Node_fields;
else
    nodes = fieldnames(output_raw.output);
end

% --- Pick one node to infer time resolution robustly ---
nf0 = nodes{1};
x0  = output_raw.output.(nf0).TOTAL(:);
L   = numel(x0);

% --- Infer timestep length (hours) from known simulation duration (6 days) ---
t_d = 6;                           % days
duration_hours = t_d * 24;         % 144 hours

% Use (L-1) when it fits, otherwise fall back to L.
dt_h_1 = duration_hours / max(L-1,1);
dt_h_2 = duration_hours / max(L,1);

common = [5/60, 10/60, 15/60, 30/60];
[~,i1] = min(abs(dt_h_1 - common));
[~,i2] = min(abs(dt_h_2 - common));

if abs(dt_h_1 - common(i1)) <= abs(dt_h_2 - common(i2))
    dt_h = dt_h_1;
else
    dt_h = dt_h_2;
end

fprintf('Detected series length L=%d -> inferred dt ≈ %.4f h (%.1f min)\n', L, dt_h, dt_h*60);

% --- Helper: sum end-use across nodes and convert to m3/h ---
sum_m3ph = @(Sstruct, fieldname) local_sum_m3ph(Sstruct, nodes, fieldname, dt_h);

% RAW (network totals, m3/h)
raw_faucet = sum_m3ph(output_raw.output, 'StFaucet');
raw_shower = sum_m3ph(output_raw.output, 'StShower');
raw_toilet = sum_m3ph(output_raw.output, 'StToilet');
raw_cw     = sum_m3ph(output_raw.output, 'StClothesWasher');
raw_dw     = sum_m3ph(output_raw.output, 'StDishwasher');
raw_total  = sum_m3ph(output_raw.output, 'TOTAL');

% ABM (network totals, m3/h)
abm_faucet = sum_m3ph(output.output, 'StFaucet');
abm_shower = sum_m3ph(output.output, 'StShower');
abm_toilet = sum_m3ph(output.output, 'StToilet');
abm_cw     = sum_m3ph(output.output, 'StClothesWasher');
abm_dw     = sum_m3ph(output.output, 'StDishwasher');
abm_total  = sum_m3ph(output.output, 'TOTAL');

% --- Time axis (days) ---
t_days = (0:L-1) * (dt_h/24);

% --- Plot (6 subplots) ---
figure('Color','w','Position',[120 80 1200 820]);

subplot(3,2,1);
plot(t_days, raw_faucet, 'k-', 'LineWidth', 2.0); hold on;
plot(t_days, abm_faucet,        'LineWidth', 1.8);
title('Faucet'); xlabel('Time (days)'); ylabel('m^3/h'); grid on; box on;
legend('RAW','ABM','Location','northeast');

subplot(3,2,2);
plot(t_days, raw_shower, 'k-', 'LineWidth', 2.0); hold on;
plot(t_days, abm_shower,        'LineWidth', 1.8);
title('Shower'); xlabel('Time (days)'); ylabel('m^3/h'); grid on; box on;
legend('RAW','ABM','Location','northeast');

subplot(3,2,3);
plot(t_days, raw_toilet, 'k-', 'LineWidth', 2.0); hold on;
plot(t_days, abm_toilet,        'LineWidth', 1.8);
title('Toilet'); xlabel('Time (days)'); ylabel('m^3/h'); grid on; box on;
legend('RAW','ABM','Location','northeast');

subplot(3,2,4);
plot(t_days, raw_cw, 'k-', 'LineWidth', 2.0); hold on;
plot(t_days, abm_cw,        'LineWidth', 1.8);
title('Clothes washer'); xlabel('Time (days)'); ylabel('m^3/h'); grid on; box on;
legend('RAW','ABM','Location','northeast');

subplot(3,2,5);
plot(t_days, raw_dw, 'k-', 'LineWidth', 2.0); hold on;
plot(t_days, abm_dw,        'LineWidth', 1.8);
title('Dishwasher'); xlabel('Time (days)'); ylabel('m^3/h'); grid on; box on;
legend('RAW','ABM','Location','northeast');

subplot(3,2,6);
plot(t_days, raw_total, 'k-', 'LineWidth', 2.2); hold on;
plot(t_days, abm_total,        'LineWidth', 2.0);
title('TOTAL'); xlabel('Time (days)'); ylabel('m^3/h'); grid on; box on;
legend('RAW','ABM','Location','northeast');

sgtitle('Moderate detection (contam 08:00): end-use demand, RAW vs ABM (m^3/h)');

% -------- local function --------
function y_m3ph = local_sum_m3ph(S, nodes, fieldname, dt_h)
    y_L_per_step = [];
    for k = 1:numel(nodes)
        nf = nodes{k};
        if ~isfield(S.(nf), fieldname), continue; end
        x = S.(nf).(fieldname)(:);   % liters per step
        if isempty(y_L_per_step)
            y_L_per_step = zeros(size(x));
        end
        m = min(numel(y_L_per_step), numel(x));
        y_L_per_step(1:m) = y_L_per_step(1:m) + x(1:m);
    end
    % Convert L/step -> m3/h: (L/step)/(1000)/dt_h
    y_m3ph = (y_L_per_step / 1000) / dt_h;
end
