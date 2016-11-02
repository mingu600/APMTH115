import numpy
import random
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

#

# Parameters
n = 20              #number of subpopulations
g = [.1]*n          #growth rate
d = [.01]*n         #death rate
mf = [1]*n        #forward mutation rate
mb = [0]*n        #backward mutation rate

# Initial condition: Single tumor cell at the origin
x, y, z = [0], [0], [0]

# Index 0: Dead positions, index 1 to n: subpopulation n
positions = [[()]] + [[(0, 0, 0)]] + [[()]] * (n-1)


def make_axes(positions):
    return [[x[coordinate] for x in positions] for coordinate in [0, 1, 2]]

def growth(sub_id):
    global g
    global positions
    sub_pop = positions[sub_id][:]
    if sub_pop != [()]:
        for position in sub_pop:
            neighbors = get_neighbors(position)
            for neighbor in neighbors:
                if random.random() < g[sub_id]:
                    positions[sub_id].append(neighbor)

def mutate(sub_id):
    global positions, mf, mb, n
    sub_pop = positions[sub_id][:]
    if sub_pop != [()]:
        if 1 < sub_id < n:
            for position in sub_pop:
                print position
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
    for x in possibles:
        if x not in positions:
            neighbors.append(x)
    return neighbors

for i in range(1, len(positions)):
    growth(i)
    mutate(i)
#print positions


# Plot the tumor as a 3D scatter plot
# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
# sep_coords = make_axes(positions)
# ax.scatter(sep_coords[0], sep_coords[1], sep_coords[2], zdir='z', c= 'red')
# plt.savefig("tumor.png")
