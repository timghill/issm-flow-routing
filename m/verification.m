addpath(genpath([getenv('ISSM_DIR'), '/bin']))
addpath(genpath([getenv('ISSM_DIR'), '/lib']))

% Load ISSM model and set surface and melt fields
% loadmodel('ASE_2300_ks_1e3_GlaDS_Steady_State_ks.mat')
phi = md.geometry.surface;
% phi(md.mask.ocean_levelset<0) = nan;

%% load surface melt data that has been interpolated to the ISSM mesh

load('melt_2300.mat')
load('moulins_2300.mat')
%% Calculate strain


md = mechanicalproperties(md, md.initialization.vx, md.initialization.vy);
eff_strain = md.results.strainrate.effectivevalue;

%% CHECK 1: strain where we have placed moulins
disp('Min strain where moulins:')
disp(min(eff_strain(moulins)))

% Allow moulins where strain threshold and discharge threshold are met
condition = eff_strain> 0.01;
discharge_threshold = 5; % m3/s

% Recompute flow accumulation for final moulins
% compute_flow_accumulation needs element areas and
% edge-node connections
step = 1;
meshArea = GetAreas(md.mesh.elements, md.mesh.x, md.mesh.y);

connectEdge = reconstruct_edges(md);
connectEdge = connectEdge(:, 1:2);
[flowAcc, sinkDischarge, sinkCatchments] = compute_flow_accumulation(md, ...
    meshArea, connectEdge, phi, melt, 'step', step, 'sinks', moulins); %sinkDischarge is the discharge going into each moulin, sinkCatchment is a cell array of the elements that feed into the moulins


%% CHECK that we are not double counting elements towards multiple moulins
drained_elements = [];
for index=moulins
    drained_elements = [drained_elements, sinkCatchments{index}];
end
% Total length and unique length should be the same,
% otherwise we have duplicates
disp(length(drained_elements))
disp(length(unique(drained_elements)))

%% CHECK discharge computed two ways
discharge_error = zeros(length(moulins), 1);
elmelt = mean(melt(md.mesh.elements), 2);
% Check that sinkDischarge is the same as that computed otherwise
% from summing melt within catchments
for index=1:length(moulins)
    d1 = sinkDischarge(index);
    d2 = sum(elmelt(sinkCatchments{moulins(index)}).*meshArea(sinkCatchments{moulins(index)}));
    discharge_error(index) = d2-d1;
end
disp('Max discharge discrepancy:')
disp(max(abs(discharge_error)))



% Compute total melt
total_melt = sum(elmelt.*meshArea);
disp('Total melt:')
disp(total_melt)

disp('Sum of moulin discharge')
disp(sum(sinkDischarge))

disp('Masked out melt')
melt(md.mask.ocean_levelset<0) = 0;
elmelt = mean(melt(md.mesh.elements), 2);
disp(sum(elmelt.*meshArea))
