function [moulins, discharge, trace] = place_moulins(md, phi, melt, condition, discharge_threshold, varargin)
% defaults
step = 1;
p = 1.0;
maxiter = 250;
decay_factor = 0.5;
nvarg = length(varargin);
for k=1:2:nvarg
    if strcmp(varargin{k}, 'step')
        step = varargin{k+1};
    elseif strcmp(varargin{k}, 'p')
        p = varargin{k+1};
    elseif strcmp(varargin{k}, 'maxiter')
        maxiter = varargin{k+1};
    elseif strcmp(varargin{k}, 'decay_factor')
        decay_factor = varargin{k+1};
    end
end

moulins = [];

% Compute element areas
% disp('place_moulins::GetAreas')
area = GetAreas(md.mesh.elements, md.mesh.x, md.mesh.y);
connect_edge = reconstruct_edges(md);
connect_edge = connect_edge(:, 1:2);

% initial flow routing map
% disp('place_moulins::compute_flow_accumulation')
[flow, discharge, ~] = compute_flow_accumulation(md, area, connect_edge, ...
    phi, melt, 'step', step, 'sinks', moulins);

% disp('place_moulins::flow')
flow(flow<discharge_threshold) = nan;
flow(~condition) = nan;
iteration = 1;
failed_candidates = [];
trace_moulins = {};
trace_discharge = {};
while not(all(isnan(flow))) && iteration<=maxiter
    iteration = iteration;
    % disp('iteration');
    % disp(iteration);

    % randomly propose a candidate moulin from nodes where
    % the accumulated flow is at least discharge_threshold
    all_choices = find(flow>=discharge_threshold);
    new_moulin = all_choices(randi(length(all_choices)));

    % accumulate flow including the candidate moulin
    newmoulins = [moulins, new_moulin];
    [newflow, newdischarge] = compute_flow_accumulation(md, area, connect_edge, ...
        phi, melt, 'step', step, 'sinks', newmoulins);

    % Optional acceptance probability
    accept = rand()<=p;

    % If adding the candidate moulin captures too much discharge
    % from an existing moulin, reject the candidate
    if min(newdischarge)<(discharge_threshold*decay_factor)
        disp(sprintf('rejected node %d', new_moulin))
        failed_candidates = [failed_candidates, new_moulin];
    elseif accept
        disp(sprintf('accepted node %d', new_moulin))
        moulins = newmoulins;
        discharge = newdischarge;
        flow = newflow;
    else
        disp(sprintf('randomly rejected node %d', new_moulin))
    end

    trace_moulins{iteration} = moulins;
    trace_discharge{iteration} = discharge;

    % mask out areas we dont want to consider in the
    % next iteration
    flow(flow<discharge_threshold) = nan;
    flow(~condition) = nan;
    flow(failed_candidates) = nan;
    flow(moulins) = nan;

    iteration = iteration + 1;
end

% Final reassignment
[flow, discharge] = compute_flow_accumulation(md, area, connect_edge, ...
    phi, melt, 'step', step, 'sinks', moulins);

trace = {trace_moulins, trace_discharge};

end
