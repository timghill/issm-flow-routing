function [moulins, discharge, trace] = place_moulins(md, phi, melt, condition, discharge_threshold, varargin)
% defaults
step = 1;   % for faster debugging set step>1 (e.g., step=250)
p = 1.0;    % acceptance probability for each candidate
maxiter = 250; % maximum number of flow accumulation iterations
decay_factor = 0.5; % absolute minimum moulin discharge is decay_factor * discharge_threshold
dz = 250; % height of elevation bins for elevation preference
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


% Set up elevation bins
zmin = 0;
zmax = max(phi(condition)) + 50;
disp('zmin:')
disp(zmin)
disp('zmax:')
disp(zmax)
zlower = 0:dz:zmax;
zupper = zlower + dz;
zcenter = 0.5*(zlower + zupper);


% initial flow routing map
% disp('place_moulins::compute_flow_accumulation')
[flow, discharge, ~] = compute_flow_accumulation(md, area, connect_edge, ...
    phi, melt, 'step', step, 'sinks', moulins);

all_moulins = {};
all_discharge = {};
nbin = length(zcenter);
for binnum=nbin:-1:1
    trace_moulins = {};
    trace_discharge = {};
    flow(flow<discharge_threshold) = nan;
    flow(~condition) = nan;
    
    elev_condition = phi>=zlower(binnum) & phi<zupper(binnum);
    flow(~elev_condition) = nan;

    disp('Bin number:')
    disp(binnum)
    disp(size(phi(elev_condition)))
    disp(sum(elev_condition))
    
    iteration = 1;
    failed_candidates = [];
    while not(all(isnan(flow))) && iteration<=maxiter
        iteration = iteration;

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
        flow(~elev_condition) = nan;

        iteration = iteration + 1;
    end

    % Final reassignment
    [flow, discharge] = compute_flow_accumulation(md, area, connect_edge, ...
        phi, melt, 'step', step, 'sinks', moulins);
    
    % all_moulins = {all_moulins, trace_moulins}
    % all_discharge = {all_discharge, trace_discharge};
    all_moulins{binnum} = trace_moulins;
    all_discharge{binnum} = trace_discharge;
end

trace = {all_moulins, all_discharge}

end
