addpath(genpath([getenv('ISSM_DIR'), '/bin']))
addpath(genpath([getenv('ISSM_DIR'), '/lib']))

% Load ISSM model and set surface and melt fields
loadmodel('ASE_2300_ks_1e3_GlaDS_Steady_State_ks.mat')
phi = md.geometry.surface;

% Simple elevation-dependent melt rate
% 2 m/year at sea level
% 0 m/year above 1500 m
ela = 1500;
melt = 2/365/86400*(ela-phi)/ela; % UNITS: m/s
melt(melt<0) = 0;

% Allow moulins everywhere
condition = phi>0;
discharge_threshold = 5;
step = 1;
maxiter = 250;
[moulins, discharge, trace] = place_moulins(md, phi, melt, condition, ...
    discharge_threshold, 'step', step, 'maxiter', maxiter);


% Recompute flow accumulation for final moulins
% compute_flow_accumulation needs element areas and
% edge-node connections
meshArea = GetAreas(md.mesh.elements, md.mesh.x, md.mesh.y);
connectEdge = reconstruct_edges(md);
connectEdge = connectEdge(:, 1:2);
[flowAcc, sinkDischarge, sinkCatchments] = compute_flow_accumulation(md, ...
    meshArea, connectEdge, phi, melt, 'step', step, 'sinks', moulins);

% Save moulin indices and catchments (i.e., the elements)
% that feed into each moulins
moulin_map = struct;
moulin_map.moulins = moulins;
moulin_map.catchments = sinkCatchments;
save('moulins.mat', '-struct', 'moulin_map');

% Scatter plot flow accumulation
% Set marker size based on flow accumulation field
s = flowAcc;
s(s>5) = 5;
s(s<=0) = 0.1;
flowAcc(md.mask.ocean_levelset<1) = nan;

figure;
hold on
scatter(md.mesh.x, md.mesh.y, s, flowAcc, 'filled')
scatter(md.mesh.x(moulins), md.mesh.y(moulins), 5, 'color', 'red')
axis image
colorbar;
print('matlab_flowacc.png', '-dpng', '-r400')

% Trace figure for moulin discharge
trace_moulins = trace{1};
trace_discharge = trace{2};
discharge_arr = zeros(length(moulins), maxiter).*nan
for i=1:maxiter
    ni = length(trace_moulins{i});
    discharge_arr(1:ni, i) = trace_discharge{i};
end

figure;
hold on;
pcolor(discharge_arr)
shading flat
colorbar
print('matlab_trace.png', '-dpng', '-r400')
