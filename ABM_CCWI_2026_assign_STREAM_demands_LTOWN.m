%% ABM_CCWI_2026_assign_STREAM_demands_LTOWN.m
try
    d.unload
catch ERR
end
fclose all;clear class;
close all;clear all;clc;
% Start EPANET MATLAB TOOLKIT
addpath(genpath(pwd))
%%% Load Network:
inpname = 'L-TOWN.inp';
dispname = 'L-TOWN';
d = epanet(inpname);
nj = d.getNodeJunctionCount;
nn = d.getNodeCount;
fprintf('Loaded network: %d junctions, %d total nodes\n', nj, nn);
% Assign demands to network:
%%% Get existing patterns:
demInd1 = double(d.getNodeDemandPatternIndex{1}(:));
demInd1G = demInd1;
demInd2 = double(d.getNodeDemandPatternIndex{2}(:));
demInd3 = double(d.getNodeDemandPatternIndex{3}(:));
%%% Zero base demands
for i=1:nn
    disp(['Zero base demand ',num2str(i)])
    d.setNodeBaseDemands(i, 1, 0)
    d.deleteNodeJunctionDemand(i,2)
    d.deleteNodeJunctionDemand(i,2)
end
%%% Delete three(3) old patterns:
d.deletePattern(1)
d.deletePattern(1) % pattern number 2 becomes 1 after deletion of 1
d.deletePattern(1)
%%% Save new input file:
emptyInpName = ['networks\',dispname,'_empty','.inp'];
d.saveInputFile(emptyInpName);
disp('EMPTY NETWORK READY!')
%%% Load new actual demands:

%%%% ================Contamination at 08:00 ===================%%%%

% load Stream_demands_ABM_ltown_SOURCE.mat  
% load Stream_demands_ABM_ltown_DETECTION.mat
% load Stream_demands_ABM_ltown_EARLY_DETECTION.mat 
% load Stream_demands_ABM_ltown_EARLY_SOURCE.mat

%%%% ================Contamination at 20:00 ===================%%%%

load Stream_demands_ABM_ltown_INVERSED_SOURCE.mat
% load Stream_demands_ABM_ltown_INVERSED_DETECTION.mat
% load Stream_demands_ABM_ltown_EARLY_INVERSED_SOURCE.mat  
% load Stream_demands_ABM_ltown_EARLY_INVERSED_DETECTION.mat

%%% Calculate and assign new base demands:

%%% When calculating ABM, please uncomment the base_Stream_demand = mean(Stream); | When calculating RAW, please uncomment the base_Stream_demand = mean(Stream_raw);
 for m=1 
    d = epanet(emptyInpName);
    demInd1 = demInd1G; 
    % base_Stream_demand = mean(Stream); % ABM
    base_Stream_demand = mean(Stream_raw); % RAW
    for i=1:nj
        disp(['Assign base demand ',num2str(i)])
        d.setNodeBaseDemands(i, 1, base_Stream_demand(i))
    end 
    %%% Calculate new patterns:
    
    %%% When calculating ABM, please uncomment the pattern_Stream_demand = Stream ./ base_Stream_demand; | When calculating RAW, please uncomment the pattern_Stream_demand = Stream_raw ./ base_Stream_demand;

    % pattern_Stream_demand = Stream ./ base_Stream_demand; % ABM
    pattern_Stream_demand = Stream_raw ./ base_Stream_demand; % RAW
    pattern_Stream_demand(isnan(pattern_Stream_demand)) = 0;
    %%% Add new patterns:
    for i = 1:nj
        demands = pattern_Stream_demand(:,i);
        resDemPatInd(i) = d.addPattern(['P-Res-',num2str(i)], demands);
        disp(['Creating pattern Residential ',num2str(i)])
    end
    for i=1:nj
        disp(['Indexing pattern ',num2str(i),' category 1'])
        if demInd1(i)==0
            continue
        elseif demInd1(i)==1 % Residential
            demInd1(i)=i;
        else
            error('unknown demand pattern')
        end
    end
    %%% Assign new patterns:
    for i=1:nn
        disp(['Assigning pattern ',num2str(i),' out of ',num2str(nn)])
        d.setNodeDemandPatternIndex(i, 1, demInd1(i))
    end 
    % Correct times:
    d.setTimeReportingStep(300)
    d.setTimeHydraulicStep(300)
    d.setTimePatternStep(300)
    % Save new input file:
    % The user has the following options for the newInpname:
    % 1. THey can switch between ABM and RAW, 
    % 2. They can switch between EARLY and BASELINE
    % 3. They can switch between SOURCE and DETECTION
    % For contamination at 20:00 the newInpname should be:'networks\',dispname,'_stream_RAW_BASELINE_MODERATE_SOURCE_INVERSED','.inp'
    % For contamination at 08:00 the newInpname should be:'networks\',dispname,'_stream_RAW_BASELINE_MODERATE_SOURCE','.inp'

    newInpname = ['networks\',dispname,'_stream_RAW_BASELINE_MODERATE_SOURCE_INVERSED','.inp']; 
    d.saveInputFile(newInpname);
    disp('NETWORK READY!')
end
