addpath(genpath([getenv('ISSM_DIR'), '/bin']))
addpath(genpath([getenv('ISSM_DIR'), '/lib']))

% Load ISSM model and set surface and melt fields
loadmodel('ASE_2300_ks_1e3_GlaDS_Steady_State_ks.mat')
phi = md.geometry.surface;


%% load surface melt data that has been interpolated to the ISSM mesh

load('melt_2300.mat')
%% Calculate strain

md = mechanicalproperties(md, md.initialization.vx, md.initialization.vy);
eff_strain = md.results.strainrate.effectivevalue;

%%
% Allow moulins where strain threshold and discharge threshold are met
elx = mean(md.mesh.x(md.mesh.elements), 2);
ely = mean(md.mesh.y(md.mesh.elements), 2);
eff_strain_vertex = griddata(elx, ely, eff_strain, md.mesh.x, md.mesh.y);
condition = eff_strain_vertex> 0.01;
discharge_threshold = 5; % m3/s
step = 1;
maxiter = 10; % increase this!
[moulins, discharge, trace] = place_moulins_highelev(md, phi, melt, condition, ...
    discharge_threshold, 'step', step, 'maxiter', maxiter, ...
    'decay_factor', 1);
%%

% Recompute flow accumulation for final moulins
% compute_flow_accumulation needs element areas and
% edge-node connections
meshArea = GetAreas(md.mesh.elements, md.mesh.x, md.mesh.y);
connectEdge = reconstruct_edges(md);
connectEdge = connectEdge(:, 1:2);
[flowAcc, sinkDischarge, sinkCatchments] = compute_flow_accumulation(md, ...
    meshArea, connectEdge, phi, melt, 'step', step, 'sinks', moulins); %sinkDischarge is the discharge going into each moulin, sinkCatchment is a cell array of the elements that feed into the moulins

% save('moulin_discharge_2300.mat','sinkDischarge');
% Save moulin indices and catchments (i.e., the elements)
% that feed into each moulins
moulin_map = struct;
moulin_map.moulins = moulins;
moulin_map.catchments = sinkCatchments;
moulin_map.flowAcc = flowAcc;
save('moulins_highelev_2300.mat', '-struct', 'moulin_map');

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
print('matlab_flowacc_highelev.png', '-dpng', '-r400')

% Trace figure for moulin discharge
nbin = size(trace{1}, 2);
figure('Position', [100, 100, 600, 1200])
hold on
for binnum=1:nbin
    subplot(nbin, 1, binnum);
    hold on;
    trace_moulins = trace{1}{binnum};
    trace_discharge = trace{2}{binnum};
    niter = size(trace_moulins, 2);
    discharge_arr = zeros(length(moulins), maxiter).*nan;
    for i=1:niter
        ni = length(trace_moulins{i});
        discharge_arr(1:ni, i) = trace_discharge{i};
    end

    pcolor(discharge_arr)
    shading flat
    colorbar
end
print('matlab_trace_highelev.png', '-dpng', '-r400')
