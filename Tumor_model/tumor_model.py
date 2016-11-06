import numpy as np
import random
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from mpl_toolkits.mplot3d import Axes3D
import time
from datetime import timedelta
import csv
import sys
np.set_printoptions(suppress=True)

#TODO: Plot Population of tumor cells over times
#TODO: Run multiple trials, find relation between B(K) and Population
#TODO: cleanup code - DONE
#TODO: Function that saves state (includes positions, boundary) so that we can use those to run again
#TODO: # Look up literature values (mutation rate of cancer (and how that relates to number of mutations))
       # Number of oncogene/tumor suppressor genes
       # Growth rate of tumors
       # Death rate of tumors

##### FLAGS #####
run_with_mutation_flag = True
save_state_flag = False
# 0: end positions, 1: all time positions, 2: subpop growth curves, 3: borders
output_flag = 2

##### PARAMETERS #####
n = 20                                                    #number of subpopulations
d = [0.02]*n                                              #death rate
mb = [0.005 for i in range(n)]                            #backward mutation rate
t_steps = 100                                             #timesteps
g = [0.005 * i**1.01 for i in range(n)] if\
    run_with_mutation_flag else [0.005 for i in range(n)] #growth rate
mf = [.1/(1+2.718**(.2*(-i+(n/4.)))) for i in range(n)] if\
    run_with_mutation_flag else [0.01 for i in range(n)]  #forward mutation rate

##### INITIALIZE #####
# index 0: Dead positions, index 1 to n: subpopulation n
positions = [[]] + [[(0, 0, 0)]] + [[] for x in range(n-1)]
positions_history = []
boundary = [(0,0,0)]

# make matplotlib plots prettier - adapted code from Adrian Veres
def prettify_ax(ax):
    for spine in ax.spines.itervalues():
        spine.set_visible(False)
    ax.set_frameon=True
    ax.patch.set_facecolor('#eeeeef')
    ax.grid('on', color='w', linestyle='-', linewidth=1)
    ax.tick_params(direction='out')
    ax.set_axisbelow(True)
def simple_ax(figsize=(6,4), **kwargs):
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111, **kwargs)
    prettify_ax(ax)
    return fig, ax

# each live boundary cell has growth probability according to subpopulation type
def growth(sub_id):
    global g, positions, boundary
    sub_pop = positions[sub_id][:]
    if sub_pop != [()]:
        outer_cells = list(set(sub_pop).intersection(set(boundary)))
        for position in outer_cells:
            neighbors = get_neighbors(position)
            for neighbor in neighbors:
                if random.random() < g[sub_id]:
                    positions[sub_id].append(neighbor)
                    if len(get_neighbors(neighbor)) > 20:
                        boundary.append(neighbor)

# each live cell can forward mutate or back mutate
def mutate(sub_id):
    global positions, mf, mb, n
    sub_pop = positions[sub_id][:]
    if sub_pop != [()]:
        if 1 < sub_id < n:
            for position in sub_pop:
                random1 = random.random()
                random2 = random.random()
                if random1 < mf[sub_id] and random2 > mb[sub_id]:
                    positions[sub_id].remove(position)
                    positions[sub_id + 1].append(position)
                elif random1 > mf[sub_id] and random2 < mb[sub_id]:
                    positions[sub_id].remove(position)
                    positions[sub_id - 1].append(position)
        elif sub_id == 1:
            for position in sub_pop:
                if random.random() < mf[sub_id]:
                    positions[sub_id].remove(position)
                    positions[sub_id + 1].append(position)
        elif sub_id == n:
            for position in sub_pop:
                if random.random() < mb[sub_id]:
                    positions[sub_id].remove(position)
                    positions[sub_id - 1].append(position)

# return all open positions to grow into surrounding a cell
def get_neighbors(position):
    global positions
    possibles = []
    for delta_x in [-1, 0, 1]:
        for delta_y in [-1, 0, 1]:
            for delta_z in [-1, 0, 1]:
                possibles.append((position[0] + delta_x, position[1] + delta_y, position[2] + delta_z))
    neighbors = []
    all_positions = [item for sublist in positions for item in sublist]
    for x in possibles:
        if x not in all_positions:
            neighbors.append(x)
    return neighbors

# cells in the interior of the tumor have probability of dying from hypoxia
def death(sub_id):
    global d, positions, boundary
    sub_pop = positions[sub_id][:]
    inner = set(sub_pop) - set(sub_pop).intersection(set(boundary))
    inner_cells = [x for x in inner]
    for cell in inner_cells:
        if random.random() < d[sub_id]:
            positions[sub_id].remove(cell)
            if cell in boundary:
                boundary.remove(cell)
            positions[0].append(cell)

# display results according to output_flag, following the end of simulation
def display_results():
    global positions_history
    print
    if save_state_flag:
        pass   #TODO
    if output_flag == 0:
        print positions_history
    elif output_flag == 1:
        pass   #already printed!
    elif output_flag == 2:
        toPlot = [[x[i] for x in positions_history] for i in range(n+1)]
        fig, ax = simple_ax(figsize=(11,8))
        ax.set_title('Subpopulations Over Time')
        ax.set_xlabel('Time')
        ax.set_ylabel('Subpopulation Size')
        for j in range(n+1):
            label = "Dead cells" if j==0 else "Subpopulation " + str(j)
            ax.plot(toPlot[j], label=label)
        ax.legend(loc=2, prop={'size':10})
        fig_filename = 'mutSubpopCurves' if run_with_mutation_flag \
                    else 'sphereSubpopCurves'
        plt.savefig('Figures/' + fig_filename)
    elif output_flag == 3:
        pass #TODO

if __name__ == '__main__':
    end = time.clock()
    for t in range(t_steps):
        try:
            for i in range(1, n):
                growth(i)
                mutate(i)
                death(i)
            for k in range(1, n+2):
                np.savetxt("Tumor_mutation/tumor_model"+str(t+1)+"_"+str(k)+".csv",\
                positions[k-1], delimiter=",", fmt= '%i')
            positions_history.append([len(x) for x in positions])
            if output_flag == 1:
                print
                print "Timestep: " + str(t)
                print [len(x) for x in positions]
                print "Population: " + str(len([item for sublist in positions[1:] for item in sublist]))
                print "Time taken: " + str(str(timedelta(seconds=round(time.clock()-end))))
                print "Boundary length: " + str(len(boundary))
            else:
                sys.stdout.write("Timestep %i done. " % t)
                sys.stdout.flush()
        except KeyboardInterrupt:
            print "\nStopping simulation and exporting current results."
            display_results()
            sys.exit(0)
    display_results()
