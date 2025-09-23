import time
import pickle
import os

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.tri import Triangulation

import cmocean


def _matrix_flow_accumulation(mesh, phi, melt, levelset, 
    sinks=np.array([]), step=250):
    """Deterministic flow routing on ISSM mesh

    Parameters:
    ----------
    mesh : dict, ISSM mesh
    phi  : (numberofvertices,) array hydraulic potential to route flow
    melt : (numberofvertices,) array melt rate in m/second
    levelset : (numberofvertices,) array ISSM ice/ocean mask
    sinks : (n,) array node indices of sinks

    Returns:
    -------
    flowacc : (numberofvertices) discharge array (m3/s)
    sink_discharge : (n,) array discharge into sinks
    """
    acc = 0*phi
    conn = mesh['elements']-1
    mdot = np.mean(melt[:,None][conn,:], axis=1).squeeze()
    
    sink_discharge = np.zeros(sinks.shape)
    
    # Assign initial melt volume from elements to nodes:
    # route all water produced in an element through the
    # lowest potential node
    # If all neighbouring nodes are nan (i.e., on ice shelf),
    # do nothing
    node2element = {}
    for element in range(mesh['numberofelements']):
        phi_neigh = phi[conn[element,:]]
        if np.all(np.isnan(phi_neigh)):
            pass
        else:
            ixmin = np.nanargmin(phi_neigh)
            acc[conn[element, ixmin]] = mesh['area'][element]*mdot[element]
            if conn[element, ixmin] in node2element.keys():
                node2element[conn[element, ixmin]].append(element)
            else:
                node2element[conn[element, ixmin]] = [element]
    # return
    flowacc = 0*phi

    if not os.path.exists('adjacent_nodes.pkl'):
        # Compute node adjacency
        print('Constructing node adjacency list')
        t0 = time.perf_counter()
        nv = mesh['numberofvertices']
        adjacent_nodes = []
        for i in range(nv):
            if i%1000==0:
                print(i)
            edgenums = np.where(np.any(mesh['connect_edge']==i, axis=1))[0]
            neigh_nodenums = mesh['connect_edge'][edgenums]
            neigh_nodenums = neigh_nodenums[neigh_nodenums!=i]
            adjacent_nodes.append(neigh_nodenums)
        t1 = time.perf_counter()
        dt = t1-t0
        print(f'{dt:.3f} s to compute node adjacency list')

        with open('adjacent_nodes.pkl', 'wb') as fout:
            pickle.dump(adjacent_nodes, fout)
    
    adjacent_nodes = np.load('adjacent_nodes.pkl', allow_pickle=True)

    # if verbose:
        # print('Plotting tricontourf')
    # fig,ax = plt.subplots()
    # mtri = Triangulation(mesh['x'], mesh['y'], conn)
    # ax.tricontourf(mtri, levelset, levels=(-1.5, 0, 1.5), colors=('gray', 'lightgray'))
    # ax.set_aspect('equal')
    
    t0 = time.perf_counter()
    paths = {}
    sink_catchments = {}
    maxiter = 1000
    groundedice = np.where(levelset==1)[0]
    for start in groundedice[::step]:
        if acc[start]>0:
            phicopy = phi.copy()
            paths[start] = []
            nodenum = start.copy()
            iters = 0
            done = False
            if start in sinks:
                sinkIndex = np.where(sinks==start)[0]
                sink_discharge[sinkIndex] += acc[start]
                done = True
            while not done and iters<maxiter:
                flowacc[nodenum] += acc[start]
                phicopy[nodenum] = np.nan
                # edgenums = np.where(np.any(mesh['connect_edge']==nodenum, axis=1))[0]
                # neigh_nodenums = mesh['connect_edge'][edgenums]
                neigh_nodenums = adjacent_nodes[nodenum]

                phi_neigh = phicopy[neigh_nodenums]
                if np.all(np.isnan(phi_neigh)) or np.any(levelset[neigh_nodenums]==-1):
                    done = True
                else:
                    next_nodenum = neigh_nodenums[np.nanargmin(phi_neigh)]
                    paths[start].append(next_nodenum)
                    nodenum = next_nodenum
            
                # Keep track of moulin discharge
                if nodenum in sinks:
                    sinkIndex = np.where(sinks==nodenum)[0][0]
                    sink_discharge[sinkIndex] += acc[start]
                    done = True
                    if not nodenum in sink_catchments.keys():
                        sink_catchments[nodenum] = []
                    sink_catchments[nodenum].extend(node2element[start])
                iters+=1
            if iters>=maxiter:
                    print('Reached maxiters')

    t1 = time.perf_counter()
    dt = t1 - t0
    print(f'{dt:.3f} s to compute flow routing')
    return flowacc, sink_discharge, sink_catchments

def place_moulins(mesh, phi, melt, levelset, condition,
    discharge_threshold=5, step=250, p=0.5, maxiter=25,
    return_trace=True, decay_factor=0.5):
    """
    Simulate moulin placement according to a discharge threshold and
    optional external constraints

    Parameters:
    ----------
        mesh : dict, ISSM mesh

        melt : (numberofvertices,) array melt rate (m/s)

        levelset : (numberofvertices,) ISSM ice/ocean mask

        condition : (numberofvertices,) boolean array
                    only consider moulins where condition is True
        
        discharge_threshold : float, minimum moulin discharge

        step : int, skip factor for rapid prototyping only

        p : float [0, 1], acceptance probability of candidate moulin

        return_trace : bool, True to return trace information

        decay_factor : float <=1, factor by which moulin discharge is
                        allowed to drop below discharge_threshold in
                        later iterations. The absoloute minimum moulin
                        discharge will be
                            discharge_threshold*decay_factor
    """
    moulins = []
    rng = np.random.default_rng()
    # initial flow-routing map
    flow,discharge,_ = _matrix_flow_accumulation(mesh, phi, melt, levelset,
        sinks=np.array([]), step=step)
    
    flow[flow<discharge_threshold] = np.nan
    flow[~condition] = np.nan
    iteration = 1
    failed_candidates = []
    trace_moulins = []
    trace_discharge = []
    while not(np.all(np.isnan(flow))) and iteration<maxiter:
        print('Iteration', iteration)

        # Randomly propose a candidate moulin from nodes where
        # the accumulated flow is at least discharge_threshold
        new_moulin = np.random.choice(np.where(flow>=discharge_threshold)[0], size=1)[0]

        # Accumulate flow including the new candidate moulin
        _moulins = np.array(moulins + [new_moulin])
        flow,_discharge,_ = _matrix_flow_accumulation(mesh, phi, melt,
            levelset, sinks=_moulins, step=step)

        # Optional acceptance probability
        accept = rng.random()<=p
        
        # If adding the candidate moulin captures too much discharge
        # from an existing moulin, reject the candidate
        if np.min(_discharge)<(discharge_threshold*decay_factor):
            failed_candidates.append(new_moulin)
            print('rejected node', new_moulin)
        elif accept:
            print('accepted node', new_moulin)
            moulins.append(new_moulin)
            discharge = _discharge
        # Random rejection
        else:
            print('randomly rejected node', new_moulin)

        # Keep track of the number of moulins and discharge
        # for each iteration
        trace_moulins.append(moulins)
        trace_discharge.append(discharge)

        # Mask out areas we don't want to consider in the
        # next iteration
        flow[flow<discharge_threshold] = np.nan
        flow[~condition] = np.nan
        flow[failed_candidates] = np.nan
        flow[moulins] = np.nan

        iteration+=1

    _, discharge,_ = _matrix_flow_accumulation(mesh, phi, melt,
            levelset, sinks=np.array(moulins), step=step)
    moulins = np.array(moulins)
    if return_trace:
        trace = (trace_moulins, trace_discharge)
        outputs = (moulins, discharge, trace)
    else:
        outputs = (moulins, discharge)
    return outputs

def plot_trace(trace):
    """
    Plot the traces produced by place_moulins
    """
    fig,ax = plt.subplots(ncols=1, figsize=(5, 7))

    moulins,discharge = trace

    nmoulin = len(moulins[-1])
    nsteps = len(moulins)
    discharge_array = np.nan*np.zeros((nmoulin, nsteps))
    for i in range(nsteps):
        for j,d in enumerate(discharge[i]):
            discharge_array[j,i] = d
    
    steps = np.arange(1, nsteps+2)
    moulins_phony = np.arange(nmoulin+1)
    pc = ax.pcolormesh(steps, moulins_phony, discharge_array, cmap=cmocean.cm.thermal,
        vmin=0)
    ax.set_yticks(moulins_phony[:-1]+0.5, moulins[-1])

    ax.set_ylabel('Moulin node number')
    ax.set_xlabel('Iteration')

    fig.colorbar(pc, label='Moulin discharge (m$^3$ s$^{-1}$)')

    return fig,ax
    

if __name__=='__main__':
    mesh = np.load('mesh.npy', allow_pickle=True)
    phi = np.load('surface.npy')
    levelset = np.load('ocean_levelset.npy')
    melt = np.ones(phi.shape)/365/86400 # 1 m/a in units of m/s


    flow,_, catchments = _matrix_flow_accumulation(mesh, phi, melt, levelset, step=1)
    print('Max discharge:', np.max(flow))

    ela = 1500
    melt = 2*np.ones(phi.shape)/365/86400*(ela-phi)/ela
    melt[melt<0] = 0

    print(melt.max()*365*86400)
    print(melt.mean()*365*86400)

    step = 1

    t1 = time.perf_counter()
    moulins, discharge, trace = place_moulins(mesh, phi, melt, levelset,
        phi>0, step=step, discharge_threshold=5, return_trace=True, 
        maxiter=500, p=1.)
    t2 = time.perf_counter()
    print('Final configuration in {:.2f} seconds:'.format(t2-t1))
    print(moulins, discharge)

    flow,_, catchments = _matrix_flow_accumulation(mesh, phi, melt, levelset, step=step,
        sinks=moulins)
    flow[levelset<1] = np.nan

    res = dict(moulins=moulins, discharge=discharge, 
        trace=trace, catchments=catchments)
    with open('moulins.pkl', 'wb') as fout:
        pickle.dump(res, fout)

    fig,ax = plot_trace(trace)
    fig.savefig('flowacc_moulins.png', dpi=400)

    print('catchments:', catchments[moulins[0]])

    fig, ax1 = plt.subplots(ncols=1, figsize=(5, 4))
    mtri = Triangulation(mesh['x'], mesh['y'], mesh['elements']-1)
    pc = pc = ax1.tripcolor(mtri, flow, vmin=0, vmax=np.nanmax(flow))
    ax1.set_aspect('equal')
    ax1.plot(mesh['x'][moulins], mesh['y'][moulins], marker='o', 
        markeredgecolor='red', markerfacecolor='none', linestyle='none')

    # To mark individual catchments for testing
    # conn = mesh['elements']-1
    # elx = np.mean(mesh['x'][conn], axis=1)
    # ely = np.mean(mesh['y'][conn], axis=1)
    # ax1.plot(elx[catchments[moulins[0]]], ely[catchments[moulins[0]]], 'b*')
    # ax1.plot(mesh['x'][moulins[0]], mesh['y'][moulins[0]], 'bo')

    fig.colorbar(pc, label='Discharge (m$^3$ s$^{-1}$)')
    
    fig.savefig('flowacc.png', dpi=400)

