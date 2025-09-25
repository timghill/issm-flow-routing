function [flowAcc, sinkDischarge, sinkCatchments] = compute_flow_accumulation(md, area, connect_edge, phi, melt, varargin)

% disp('compute_flow_accumulation::optional args')
% handle optional arguments
sinks = [];
step = 250;
if length(varargin)>0
    k = 1;
    while k<length(varargin)
        % disp(k)
        if strcmp(varargin{k}, 'sinks')
            sinks = varargin{k+1};
            k = k+2;
        elseif strcmp(varargin{k}, 'step')
            step = varargin{k+1};
            k = k+2;
        end
    end
end


% compute melt on elements
% disp('place_moulins::conn')
conn = md.mesh.elements;
mdot = mean(melt(conn), 2);

% Assign initial melt volume from elements to nodes:
% route all water produced in an element through the
% lowest potential node
% if all neighbouring nodes are nan (i.e., on ice shelf),
% do nothing
% tic;
node2element = {};
for element=1:md.mesh.numberofvertices
    phineigh = phi(conn(element, :));
%    if element<20
%        size(phineigh)
%        phineigh
%    end
    if not(all(isnan(phineigh)))
        [~, ixmin] = min(phineigh);
%        if element<20
%            ixmin
%        end
        acc(conn(element, ixmin)) = area(element).*mdot(element);
        kk = conn(element, ixmin);
        if kk>length(node2element)
            node2element{kk} = [];
        end
        node2element{conn(element, ixmin)} = [node2element{conn(element, ixmin)}, element];
    end
end
% toc;

% tic;
if not(isfile('adjacent_nodes.mat'))
    % compute node adjacency, the connection between
    % nodes and nodes
    disp('Computing node adjacency list')
    % t0 = tic;
    nv = md.mesh.numberofvertices;
    adjacent_nodes = {};
    for i=1:nv
        edgenums = find(any(connect_edge==i, 2));
        neigh_nodenums = connect_edge(edgenums,:);
        neigh_nodenums = neigh_nodenums(neigh_nodenums~=i);
        adjacent_nodes{i} = neigh_nodenums;
    end
    % toc
    save('adjacent_nodes.mat', 'adjacent_nodes')
end

adjacent_nodes = load('adjacent_nodes.mat');
adjacent_nodes = adjacent_nodes.adjacent_nodes;

% initialize function outputs
flowAcc = 0.*phi;
sinkDischarge = 0.*sinks;
sinkCatchments = cell(md.mesh.numberofvertices, 1);
maxiter = 1000;

% tic;
groundedice = find(md.mask.ocean_levelset==1);
for jj=1:step:length(groundedice)
    start = groundedice(jj);
    if acc(start)>0
        phicopy = phi(:,:);
        nodenum = groundedice(jj);
        iters = 0;
        done = false;
        if ismember(start, sinks)
            sinkIndex = find(sinks==start);
            sinkDischarge(sinkIndex) = sinkDischarge(sinkIndex) + acc(start);
            done = true;
        end

        while not(done) && iters<maxiter
            flowAcc(nodenum) = flowAcc(nodenum) + acc(start);
            phicopy(nodenum) = nan;

            neigh_nodenums = adjacent_nodes{nodenum};
            phi_neigh = phicopy(neigh_nodenums);
            if all(isnan(phi_neigh)) || any(md.mask.ocean_levelset(neigh_nodenums)==-1)
                done = true;
                % disp('setting done = true')
            else
                [~,ixmin] = min(phi_neigh, [], 'omitnan');
                next_nodenum = neigh_nodenums(ixmin);
                nodenum = next_nodenum;
            end
            
            if ismember(nodenum, sinks)
                sinkIndex = find(sinks==nodenum);
                sinkDischarge(sinkIndex) = sinkDischarge(sinkIndex) + acc(start);
                done = true;
                sinkCatchments{nodenum} = [sinkCatchments{nodenum}, node2element{start}];
            end        
    
            iters = iters + 1;
            %if mod(iters,100)==0
            %    disp('iters')
            %    iters
            %end
        end
    end
end
% toc;
