 %% ABM_CCWI_2026 setup and run MSX

clear; clc;
% Start EPANET MATLAB TOOLKIT
addpath(genpath(pwd))

%%%==========Set up inpname==========%%%%%

%% BASELINE DETECTION | Contamination at 08:00

% inpname = 'networks\L-TOWN_stream_ABM_BASELINE_MODERATE_DETECTION.inp'; % ΑΒΜ BEST-CASE BASELINE  
% inpname = 'networks\L-TOWN_stream_RAW_BASELINE_MODERATE_DETECTION.inp'; % ΑΒΜ BEST-CASE BASELINE  
% inpname = 'networks\L-TOWN_stream_ABM_BASELINE_MODERATE_SOURCE.inp'; % ΑΒΜ BEST-CASE BASELINE  
% inpname = 'networks\L-TOWN_stream_RAW_BASELINE_MODERATE_SOURCE.inp'; % ΑΒΜ BEST-CASE BASELINE  

%% EARLY DETECTION| Contamination at 08:00

% inpname = 'networks\L-TOWN_stream_ABM_EARLY_MODERATE_DETECTION.inp'; % ΑΒΜ BEST-CASE BASELINE  
% inpname = 'networks\L-TOWN_stream_RAW_EARLY_MODERATE_DETECTION.inp'; % ΑΒΜ BEST-CASE BASELINE  
% inpname = 'networks\L-TOWN_stream_ABM_EARLY_MODERATE_SOURCE.inp'; % ΑΒΜ BEST-CASE BASELINE  
% inpname = 'networks\L-TOWN_stream_RAW_EARLY_MODERATE_SOURCE.inp'; % ΑΒΜ BEST-CASE BASELINE 

%% BASELINE DETECTION | Contamination at 20:00

% inpname = 'networks\L-TOWN_stream_ABM_BASELINE_MODERATE_DETECTION_INVERSED.inp'; % ΑΒΜ BEST-CASE BASELINE  
% inpname = 'networks\L-TOWN_stream_RAW_BASELINE_MODERATE_DETECTION_INVERSED.inp'; % ΑΒΜ BEST-CASE BASELINE  
inpname = 'networks\L-TOWN_stream_ABM_BASELINE_MODERATE_SOURCE_INVERSED.inp'; % ΑΒΜ BEST-CASE BASELINE  
% inpname = 'networks\L-TOWN_stream_RAW_BASELINE_MODERATE_SOURCE_INVERSED.inp'; % ΑΒΜ BEST-CASE BASELINE  

%% EARLY DETECTION| Contamination at 20:00

% inpname = 'networks\L-TOWN_stream_ABM_EARLY_MODERATE_DETECTION_INVERSED.inp'; % ΑΒΜ BEST-CASE BASELINE  
% inpname = 'networks\L-TOWN_stream_RAW_EARLY_MODERATE_DETECTION_INVERSED.inp'; % ΑΒΜ BEST-CASE BASELINE  
% inpname = 'networks\L-TOWN_stream_ABM_EARLY_MODERATE_SOURCE_INVERSED.inp'; % ΑΒΜ BEST-CASE BASELINE  
% inpname = 'networks\L-TOWN_stream_RAW_EARLY_MODERATE_SOURCE_INVERSED.inp'; % ΑΒΜ BEST-CASE BASELINE 



%%%==========Set up msxname==========%%%%%


%% BASELINE DETECTION | Contamination at 08:00

% msxname = 'networks\L-TOWN_Stream_ABM_contamination_BASELINE_MODERATE_DETECTION.msx'; % ABM % BEST-CASE BASELINE 
% msxname = 'networks\L-TOWN_Stream_RAW_contamination_BASELINE_MODERATE_DETECTION.msx'; % ABM % BEST-CASE BASELINE 
% msxname = 'networks\L-TOWN_Stream_ABM_contamination_BASELINE_MODERATE_SOURCE.msx'; % ABM % BEST-CASE BASELINE 
% msxname = 'networks\L-TOWN_Stream_RAW_contamination_BASELINE_MODERATE_SOURCE.msx'; % ABM % BEST-CASE BASELINE 

%% EARLY DETECTION| Contamination at 08:00

% msxname = 'networks\L-TOWN_Stream_ABM_contamination_EARLY_MODERATE_DETECTION.msx'; % ABM % BEST-CASE BASELINE 
% msxname = 'networks\L-TOWN_Stream_RAW_contamination_EARLY_MODERATE_DETECTION.msx'; % ABM % BEST-CASE BASELINE 
% msxname = 'networks\L-TOWN_Stream_ABM_contamination_EARLY_MODERATE_SOURCE.msx'; % ABM % BEST-CASE BASELINE 
% msxname = 'networks\L-TOWN_Stream_RAW_contamination_EARLY_MODERATE_SOURCE.msx'; % ABM % BEST-CASE BASELINE 

%% BASELINE DETECTION | Contamination at 20:00

% msxname = 'networks\L-TOWN_Stream_ABM_contamination_BASELINE_MODERATE_DETECTION_INVERSED.msx'; % ABM % BEST-CASE BASELINE 
% msxname = 'networks\L-TOWN_Stream_RAW_contamination_BASELINE_MODERATE_DETECTION_INVERSED.msx'; % ABM % BEST-CASE BASELINE 
msxname = 'networks\L-TOWN_Stream_ABM_contamination_BASELINE_MODERATE_SOURCE_INVERSED.msx'; % ABM % BEST-CASE BASELINE 
% msxname = 'networks\L-TOWN_Stream_RAW_contamination_BASELINE_MODERATE_SOURCE_INVERSED.msx'; % ABM % BEST-CASE BASELINE 

%% EARLY DETECTION| Contamination at 20:00

% msxname = 'networks\L-TOWN_Stream_ABM_contamination_EARLY_MODERATE_DETECTION_INVERSED.msx'; % ABM % BEST-CASE BASELINE 
% msxname = 'networks\L-TOWN_Stream_RAW_contamination_EARLY_MODERATE_DETECTION_INVERSED.msx'; % ABM % BEST-CASE BASELINE 
% msxname = 'networks\L-TOWN_Stream_ABM_contamination_EARLY_MODERATE_SOURCE_INVERSED.msx'; % ABM % BEST-CASE BASELINE 
% msxname = 'networks\L-TOWN_Stream_RAW_contamination_EARLY_MODERATE_SOURCE_INVERSED.msx'; % ABM % BEST-CASE BASELINE 


fprintf('=================================================================\n');
fprintf('   MSX SETUP - L-TOWN          \n');
fprintf('=================================================================\n\n');

%% Load network
d = epanet(inpname);

fprintf('Network: L-TOWN\n');
fprintf('  Junctions: %d\n', d.getNodeJunctionCount);
fprintf('  Links: %d\n', d.getLinkCount);

%% Timeline setup
t_d = 6;
Tts = 1729;
d.setTimeHydraulicStep(300);
d.setTimeReportingStep(300);
duration_hrs = t_d*24;
duration_sec = duration_hrs*60*60;
d.setTimeSimulationDuration(duration_sec) 

%%%==== Load Stream demands for contamination at 08:00=====%%%%

% load Stream_demands_ABM_ltown_DETECTION.mat
% load Stream_demands_ABM_ltown_SOURCE.mat
% load Stream_demands_ABM_ltown_EARLY_DETECTION.mat
% load Stream_demands_ABM_ltown_EARLY_SOURCE.mat


%%%==== Load Stream demands for contamination at 08:00=====%%%%

% load Stream_demands_ABM_ltown_INVERSED_DETECTION.mat
load Stream_demands_ABM_ltown_INVERSED_SOURCE.mat
% load Stream_demands_ABM_ltown_EARLY_INVERSED_DETECTION.mat
% load Stream_demands_ABM_ltown_EARLY_INVERSED_SOURCE.mat

dt_min = 5;
t_day = @(d_) min(max(round((d_*24*60)/dt_min) + 1, 1), Tts);

contamination_start_timestep = t_contam;          
detection_timestep = t_adv1;                   
contamination_end_timestep = t_adv2;         

fprintf('\nTimeline:\n');
fprintf('  Start: Day 2, 8 AM (ts %d)\n', contamination_start_timestep);
fprintf('  End: Day 3, 12:30 PM (ts %d)\n', contamination_end_timestep);
fprintf('  Detection: Day 3, 8 AM (ts %d)\n', detection_timestep);

%% CONTAMINATION PARAMETERS
Pathogen_concentration = 1e6;  % PFU/L (Enterovirus)
Injection_rate = 100;  % L/h

% MASS = concentration * flow rate (PFU/h)
Pathogen_mass = Pathogen_concentration * Injection_rate;  % PFU/h

% Organic reducing agents
TOC_concentration = 140;  % mg/L
C_FRA_fraction = 0.4;
C_SRA_fraction = 0.6;

C_FRA_mass = C_FRA_fraction * TOC_concentration * Injection_rate;  % mg/h
C_SRA_mass = C_SRA_fraction * TOC_concentration * Injection_rate;  % mg/h

fprintf('\nContamination (MASS INJECTION):\n');
fprintf('  Pathogen: %.2e PFU/L at %.0f L/h\n', Pathogen_concentration, Injection_rate);
fprintf('  Pathogen MASS: %.2e PFU/h \n', Pathogen_mass);
fprintf('  C_FRA MASS: %.2f mg/h\n', C_FRA_mass);
fprintf('  C_SRA MASS: %.2f mg/h\n', C_SRA_mass);

%% Create patterns
PathPAT = zeros(1, Tts);
C_FRAPAT = zeros(1, Tts);
C_SRAPAT = zeros(1, Tts);

PathPAT(contamination_start_timestep:contamination_end_timestep) = 1;
C_FRAPAT(contamination_start_timestep:contamination_end_timestep) = 1;
C_SRAPAT(contamination_start_timestep:contamination_end_timestep) = 1;

%% Define MSX structure
msx = {};
msx.FILENAME = msxname;
msx.TITLE = {'L-TOWN MSX - Mass injection with dilution'};
msx.AREA_UNITS = 'FT2';
msx.RATE_UNITS = 'HR';
msx.SOLVER = 'RK5';
msx.TIMESTEP = 300;
msx.RTOL = 0.001;
msx.ATOL = 0.001;

msx.SPECIES(1) = {'BULK CL2 MG 0.01 0.001'};
msx.SPECIES(2) = {'BULK P CFU 0.01 0.001'};
msx.SPECIES(3) = {'BULK C_FRA MG 0.01 0.001'};
msx.SPECIES(4) = {'BULK C_SRA MG 0.01 0.001'};


msx.COEFFICIENTS(1) = {'CONSTANT T 12'};
msx.COEFFICIENTS(2) = {'CONSTANT Kp 92.3'};  % Enterovirus inactivation rate
msx.COEFFICIENTS(3) = {'CONSTANT A 1'};
msx.COEFFICIENTS(4) = {'CONSTANT B 14'};

msx.TERMS(1) = {'Km (1.5826e-04 * RE^0.88) / D'};
msx.TERMS(2) = {'KWAL A*EXP(-B*CL2)'};
msx.TERMS(3) = {'K_FAST 0.2808'};
msx.TERMS(4) = {'K_SLOW 0.0071'};

msx.PIPES(1) = {'RATE CL2 -K_FAST*C_FRA*CL2-K_SLOW*C_SRA*CL2-((4/D)*(KWAL/(1+KWAL/Km)))*CL2'};
msx.PIPES(2) = {'RATE P -Kp*P*CL2'};
msx.PIPES(3) = {'RATE C_FRA -K_FAST*C_FRA*CL2'};
msx.PIPES(4) = {'RATE C_SRA -K_SLOW*C_SRA*CL2'};

msx.TANKS(1) = {'RATE CL2 -K_FAST*C_FRA*CL2-K_SLOW*C_SRA*CL2'};
msx.TANKS(2) = {'RATE P -Kp*P*CL2'};
msx.TANKS(3) = {'RATE C_FRA -K_FAST*C_FRA*CL2'};
msx.TANKS(4) = {'RATE C_SRA -K_SLOW*C_SRA*CL2'};

msx.SOURCES(1) = {'SETPOINT R1 CL2 0 CL2PAT'};
msx.SOURCES(2) = {'SETPOINT R1 C_SRA 1.85 C_SRA_REPAT'};
msx.SOURCES(3) = {'SETPOINT R2 CL2 0 CL2PAT'};
msx.SOURCES(4) = {'SETPOINT R2 C_SRA 1.85 C_SRA_REPAT'};

msx.PATTERNS(1) = {'PathPAT 1'};
msx.PATTERNS(2) = {'CL2PAT 1'};
msx.PATTERNS(3) = {'C_FRAPAT 1'};
msx.PATTERNS(4) = {'C_SRAPAT 1'};
msx.PATTERNS(5) = {'C_SRA_REPAT 1'};

%% Pipe parameters
nl = d.getLinkCount;
fprintf('\nGenerating pipe parameters for %d links...\n', nl);

for m = 1:nl
    pipe_name = d.getLinkNameID{m};

    if d.getLinkRoughnessCoeff(m) >= 140
        A = 0.01;
    else
        A = rand(1);
    end

    msx.PARAMETERS(m) = {['PIPE ', pipe_name, ' A ', num2str(A)]};
    msx.PARAMETERS(nl+m) = {['PIPE ', pipe_name, ' B 14']};

    if mod(m, 100) == 0
        fprintf('  %d/%d\n', m, nl);
    end
end

fprintf('✓ Pipe parameters complete\n');

%% Write and load MSX
fprintf('\nWriting MSX file...\n');
d.writeMSXFile(msx);
d.loadMSXFile(msxname);

%% Set contamination source (MASS INJECTION!)
source_node = 'n129';

fprintf('\nSetting contamination source at %s...\n', source_node);
fprintf('  USING MASS INJECTION (PFU/h)\n');

% MASS injection (PFU/h, mg/h) - gets diluted by nodal flow
d.setMSXSources(source_node, 'P', 'MASS', Pathogen_mass, 'PathPAT');
d.setMSXSources(source_node, 'C_FRA', 'MASS', C_FRA_mass, 'C_FRAPAT');
d.setMSXSources(source_node, 'C_SRA', 'MASS', C_SRA_mass, 'C_SRAPAT');

d.setMSXPattern('PathPAT', PathPAT);
d.setMSXPattern('C_FRAPAT', C_FRAPAT);
d.setMSXPattern('C_SRAPAT', C_SRAPAT);

fprintf('  P: %.2e CFU/h\n', Pathogen_mass);
fprintf('  C_FRA: %.2f mg/h\n', C_FRA_mass);
fprintf('  C_SRA: %.2f mg/h\n', C_SRA_mass);


%% Run MSX
fprintf('\nRunning MSX...\n');
tic;
CL2_all = d.getMSXComputedQualitySpecie('CL2');
P_all = d.getMSXComputedQualitySpecie('P');
C_FRA = d.getMSXComputedQualitySpecie('C_FRA');
C_SRA = d.getMSXComputedQualitySpecie('C_SRA');
fprintf('Complete (%.1f min)\n', toc/60);
%% Extract
CL2_matrix = CL2_all.NodeQuality;
P_matrix = P_all.NodeQuality;
C_FRA_matrix = C_FRA.NodeQuality;
C_SRA_matrix = C_SRA.NodeQuality;
junction_indices = d.getNodeJunctionIndex;
junction_names = d.getNodeNameID(junction_indices);
CL2_junctions = CL2_matrix(:, junction_indices);
P_junctions = P_matrix(:, junction_indices);
Time_seconds = (0:size(P_matrix,1)-1)' * 300;
Time_hours = Time_seconds / 3600;
Time_days = Time_hours / 24;


%% Save metadata
contamination_metadata = struct();
contamination_metadata.source_node_id = source_node;
contamination_metadata.contamination_start_timestep = contamination_start_timestep;
contamination_metadata.detection_timestep = detection_timestep;
contamination_metadata.contamination_end_timestep = contamination_end_timestep;
contamination_metadata.pathogen_concentration = Pathogen_concentration;
contamination_metadata.injection_rate_Lh = Injection_rate;
contamination_metadata.pathogen_mass_CFUh = Pathogen_mass;
contamination_metadata.dt_min = dt_min;

%% Save

% Save new input file:
    % The user has the following options for the newInpname:
    % 1. THey can switch between ABM and RAW, 
    % 2. They can switch between EARLY and BASELINE
    % 3. They can switch between SOURCE and DETECTION
    % For contamination at 20:00 the save should be:'MSX_results_LTOWN_ABM_BASELINE_MODERATE_SOURCE_INVERSED.mat'
    % For contamination at 08:00 the save should be:'MSX_results_LTOWN_ABM_BASELINE_MODERATE_SOURCE.mat'

save('MSX_results_LTOWN_ABM_BASELINE_MODERATE_SOURCE_INVERSED.mat','contamination_metadata', 'msx','msxname',...
    'CL2_matrix', 'P_matrix', 'CL2_junctions', 'P_junctions','C_FRA_matrix','C_SRA_matrix', ...
    'Time_seconds', 'Time_hours', 'Time_days', ...
    'junction_indices', 'junction_names', '-v7.3'); 
fprintf('\n✓ Saved: MSX_results_LTOWN_ABM_SHIFTED.mat\n'); 