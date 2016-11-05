import numpy as np
import random
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from mpl_toolkits.mplot3d import Axes3D
import pdb
import time
from datetime import timedelta
import csv
np.set_printoptions(suppress=True)

#

# Parameters
n = 20              #number of subpopulations
# g = [0.01 * i for i in range(n)]          #growth rate
# d = [0.1]*n                               #death rate
# mf = [0.015 * i for i in range(n)]          #forward mutation rate
# mb = [0.005 * i for i in range(n)]        #backward mutation rate

g = [0.005 for i in range(n)]          #growth rate
d = [0.1]*n         #death rate
mf = [0.01 for i in range(n)]          #forward mutation rate
mb = [0.005 for i in range(n)]        #backward mutation rate

# Initial condition: Single tumor cell at the origin
x, y, z = [0], [0], [0]

# Index 0: Dead positions, index 1 to n: subpopulation n
positions = [[]] + [[(0, 0, 0)]] + [[] for x in range(n-1)]
boundary = [(0,0,0)]

def make_axes(positions):
    axes = [[], [], []]
    for sub_pop in range(len(positions)):
        for point in range(len(positions[sub_pop])):
            for coord in [0, 1, 2]:
                axes[coord].append(positions[sub_pop][point][coord])
    return axes

def make_plot(positions):
    return [[x[coordinate] for x in positions] for coordinate in [0, 1, 2]]


def growth(sub_id):
    global g, positions, boundary
    sub_pop = positions[sub_id][:]
    if sub_pop != [()]:
        outer = set(sub_pop).intersection(set(boundary))
        outer_cells = [x for x in outer]
        for position in outer_cells:
            neighbors = get_neighbors(position)
            for neighbor in neighbors:
                if random.random() < g[sub_id]:
                    positions[sub_id].append(neighbor)
                    if len(get_neighbors(neighbor)) > 4:
                        boundary.append(neighbor)

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
#
# def get_boundary():
#     global positions
#     boundary = []
#     all_positions = [item for sublist in positions[1:] for item in sublist]
#     for position in all_positions:
#         if len(get_neighbors(position)) > 0:
#             boundary.append(position)
#     return boundary

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
            boundary_neighbors = get_neighbors(cell)
            for neighbor in boundary_neighbors:
                if len(get_neighbors(neighbor)) <= 4 and neighbor in boundary:
                    boundary.remove(neighbor)
            positions[0].append(cell)

if __name__ == "__main__":
    g = [0.005 for i in range(n)]          #growth rate
    d = [0.1]*n         #death rate
    mf = [0.01 for i in range(n)]          #forward mutation rate
    mb = [0.005 for i in range(n)]        #backward mutation rate
    boundary = [(0,0,0)]
    positions_over_time = []
    for t in range(300):
        for i in range(1, n):
            growth(i)
            mutate(i)
            death(i)
        for k in range(1, n+2):
            np.savetxt("Tumor_sphere/tumor_sphere" + str(t+1) + "_" + str(k) + ".csv", positions[k-1], delimiter=",", fmt= '%i')
        #print "TIME: " + str(t) + ", POPULATION: " + str(len([item for sublist in positions[1:] for item in sublist]))
        print "TIME: " + str(t)
        print [len(x) for x in positions]
        print "POPULATION: " + str(len([item for sublist in positions[1:] for item in sublist]))
        print str(timedelta(seconds=round(time.clock())))
        print len(boundary)
    # end = time.clock()
    # g = [0.005 * i**1.01 for i in range(n)]          #growth rate
    # d = [0.1]*n                               #death rate
    # mf = [0.0001 * 10 ** i for i in range(n)]          #forward mutation rate
    # mb = [0.005 * i for i in range(n)]        #backward mutation rate
    # positions_over_time = []
    # # Initial condition: Single tumor cell at the origin
    # x, y, z = [0], [0], [0]
    #
    # # Index 0: Dead positions, index 1 to n: subpopulation n
    # positions = [[]] + [[(0, 0, 0)]] + [[] for x in range(n-1)]
    # boundary = [(0,0,0)]
    #
    # for t in range(300):
    #     for i in range(1, n):
    #         growth(i)
    #         #print "Grew"
    #         mutate(i)
    #         #print "Mutated"
    #         death(i)
    #         #print "Died"
    #     for k in range(1, n+2):
    #         np.savetxt("Tumor_mutation/tumor_model" + str(t+1) + "_" + str(k) + ".csv", positions[k-1], delimiter=",", fmt= '%i')
    #     #print "TIME: " + str(t) + ", POPULATION: " + str(len([item for sublist in positions[1:] for item in sublist]))
    #     print "TIME: " + str(t)
    #     print [len(x) for x in positions]
    #     print "POPULATION: " + str(len([item for sublist in positions[1:] for item in sublist]))
    #     print str(timedelta(seconds=round(time.clock()-end)))
