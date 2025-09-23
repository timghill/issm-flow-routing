addpath(genpath([getenv('ISSM_DIR'), '/bin']))
addpath(genpath([getenv('ISSM_DIR'), '/lib']))

loadmodel('ASE_2300_ks_1e3_GlaDS_Steady_State_ks.mat')
phi = md.geometry.surface;
melt = ones(md.mesh.numberofvertices,1)/365/86400;

% Compute element areas
meshArea = GetAreas(md.mesh.elements, md.mesh.x, md.mesh.y);
connectEdge = reconstruct_edges(md);
connectEdge = connectEdge(:, 1:2);

condition = phi>0;
discharge_threshold = 5;
step = 1;
maxiter = 100;
[moulins, discharge, trace] = place_moulins(md, phi, melt, condition, ...
    discharge_threshold, 'step', step, 'maxiter', maxiter);

[flowAcc, sinkDischarge, sinkCatchments] = compute_flow_accumulation(md, ...
    meshArea, connectEdge, phi, melt, 'step', step, 'sinks', moulins);


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

for i=1:length(moulins)
    m1 = moulins(i);
    load adjacent_nodes;
    xx = md.mesh.x(adjacent_nodes{m1});
    yy = md.mesh.y(adjacent_nodes{m1});
    scatter(xx, yy, 5, 'color', 'magenta')
    scatter(md.mesh.x(m1), md.mesh.y(m1), 8, 'color', 'black')
end
print('matlab_flowacc.png', '-dpng', '-r400')


% Plot the trace
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
