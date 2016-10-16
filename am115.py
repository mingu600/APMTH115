#
# ALEXANDER MUNOZ | MINGU KIM | ZHAODONG CHEN
# Kidney Dialysis: Stochastic Simulation
# Applied Math 115
#
# Models kidney dialysis as a function of flow rate and flow direction of
# dialysis fluid / patient blood. Urea is cleared from blood when dialysis
# fluid is hypotonic to patient blood.
#

#TODO: Figure of urea concentration in body over time
#TODO: Optimal velocity of dialysis fluid flow rate
#TODO: Individual blood velocity
#TODO: Robustness of results (change number of boxes- do results converge?)
#TODO: Parameter optimization


from __future__ import division
import numpy as np
import random
import matplotlib.pyplot as plt
import pdb
import math
np.set_printoptions(suppress=False)

# SPECIFY DIRECTION : TRUE FOR REVERSE FLOW | FALSE FOR FORWARD FLOW
reverse_direction = True

# SPECIFY INITIAL CONDITIONS FOR BLOOD TUBING : TRUE FOR FULL | FALSE FOR EMPTY
init_full_blood_flag = False

# SPECIFY OUTPUT : 0 = PRINT STATS | 1 = HEATMAP | 2 = PATIENT UREA SCATTER
figure_to_show = 1

#PARAMETERS
time_steps = 3 * 10**3 #seconds, assuming ~3hours dialysis
length_segments = 10**3 #assuming foot-long tube, each compartment is .3048mm
pt_blood_urea_init = 5*10**7 #initial urea molecules in patient's blood
frac_pt_blood_in_dial = .1 #10% of patient's blood volume fits in dialysizer
blood_velocity = 13 #assuming, 400ml/min blood pump rate
dialysis_velocity = 20 #assuming, 600ml/min blood pump rate
diffusion_constant = 0.3 #diffusion proportional to difference in concentration

def init_parameters():
    global time_steps, length_segments, pt_blood_urea_init, frac_pt_blood_in_dial, blood_velocity, dialysis_velocity, diffusion_constant
    time_steps = 3 * 10**3 #seconds, assuming ~3hours dialysis
    length_segments = 10**3 #assuming foot-long tube, each compartment is .3048mm
    pt_blood_urea_init = 5*10**7 #initial urea molecules in patient's blood
    frac_pt_blood_in_dial = .1 #10% of patient's blood volume fits in dialysizer
    blood_velocity = 13 #assuming, 400ml/min blood pump rate
    dialysis_velocity = 20 #assuming, 600ml/min blood pump rate
    diffusion_constant = 0.3 #diffusion proportional to difference in concentration

def prettify_ax(ax):
    """
    Nifty function we can use to make our axes more pleasant to look at
    """
    for spine in ax.spines.itervalues():
        spine.set_visible(False)
    ax.set_frameon=True
    ax.patch.set_facecolor('#eeeeef')
    ax.grid('on', color='w', linestyle='-', linewidth=1)
    ax.tick_params(direction='out')
    ax.set_axisbelow(True)

def simple_ax(figsize=(6,4), **kwargs):
    """
    Shortcut to make and 'prettify' a simple figure with 1 axis
    """
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111, **kwargs)
    prettify_ax(ax)
    return fig, ax


def initialize():
    global blood, dialysis, diffused_count, time_steps_count, urea_cleared_flag, blood_list, dialysis_list
    blood = np.zeros(length_segments)
    if init_full_blood_flag:
        blood = np.full(length_segments, pt_blood_urea_init * \
        frac_pt_blood_in_dial / length_segments)
    dialysis = np.zeros(length_segments)
    diffused_count = 0
    time_steps_count = 0
    urea_cleared_flag = False
    blood_list = np.zeros((time_steps, length_segments))
    dialysis_list = np.zeros((time_steps, length_segments))






if figure_to_show ==2:
    pt_urea_list = np.zeros((2,time_steps))
else:
    pt_urea_list = np.zeros(time_steps)

for i in range(min(2, math.factorial(figure_to_show))):
    initialize()
    init_parameters()
    for t in range(time_steps):
        time_steps_count += 1
        curr_pt_urea = pt_blood_urea_init - diffused_count
        if figure_to_show ==2:
            pt_urea_list[i, t] = curr_pt_urea
        else:
            pt_urea_list[t] = curr_pt_urea
        urea_step = curr_pt_urea * frac_pt_blood_in_dial / length_segments
        if urea_step < 1: #all urea cleared from patient
            urea_cleared_flag = True
            break
        blood = np.roll(blood, blood_velocity)
        blood[:blood_velocity] = urea_step
        if reverse_direction:
            dialysis = np.roll(dialysis, -dialysis_velocity)
            dialysis[(length_segments-dialysis_velocity):] = 0
        else:
            dialysis = np.roll(dialysis, dialysis_velocity)
            dialysis[:dialysis_velocity] = 0
        for j in range(length_segments):
            #diffusion_rate = random.random() / 3
            step_diffusion = diffusion_constant * (blood[j] - dialysis[j])
            blood[j] -= step_diffusion
            dialysis[j] += step_diffusion
            diffused_count += step_diffusion
        blood_list[t] = blood
        dialysis_list[t] = dialysis
    if figure_to_show == 2:
        reverse_direction = not reverse_direction

# print results
if figure_to_show == 0:
    print "Dialysis Fluid:"
    print dialysis
    print "\nBlood Fluid:"
    print blood
    print "\nTime steps passed: %i" % time_steps_count
    print "Amount of urea diffused: %f mmol" % (diffused_count/10**6)
    if urea_cleared_flag:
        print "ALL UREA CLEARED FROM PATIENT"
    else:
        urea_remaining = pt_blood_urea_init - diffused_count
        print "Urea remaining in patient: %f mmol" % (urea_remaining/10**6)

# plot heatmap
elif figure_to_show == 1:
    f, axarr = plt.subplots(1, sharex=True)
    # axarr[0].plot(list(range(length_segments)), dialysis, label='Dialysis Tube')
    # axarr[0].plot(list(range(length_segments)), blood, label='Blood')
    # axarr[0].legend()
    # axarr[0].set_title('Urea Concentration in Blood and Dialysis Tube')
    # axarr[0].set_xlabel('Position')
    # axarr[0].set_ylabel('Urea Concentration')
    # axarr[0].set_ylabel('Urea Concentration')
    heatmap = axarr.imshow(blood_list, extent=[0,length_segments,time_steps, 0])
    heatmap.set_cmap('Greys_r')
    axarr.set_ylim([0,time_steps])
    axarr.set_xlabel('Position')
    axarr.set_ylabel('Time Step')
    plt.colorbar(heatmap)
    # plt.figure(figsize = (10, 10))
    # heatmap = plt.imshow(blood_list, cmap='Greys_r', extent=[0,length_segments, 0, time_steps])
    # plt.colorbar(heatmap)
    plt.show()

# plot urea graph
elif figure_to_show == 2:
    fig, ax = simple_ax(figsize=(6,6))
    ax.scatter(range(time_steps), pt_urea_list[0]/10**6, c='blue', label = 'Reverse Flow', edgecolor='b')
    ax.scatter(range(time_steps), pt_urea_list[1]/10**6, c='red', label='Forward Flow', edgecolor='r')
    ax.legend()
    ax.set_title('Comparing Urea Levels in Patient for Forward and Reverse Flows')
    ax.set_xlabel('Time (in seconds)')
    ax.set_ylabel('Urea Concentration in Blood (mmol/L)')
    plt.show()
    # plt.scatter(range(time_steps), pt_urea_list[0]/10**6)
    # plt.scatter(range(time_steps), pt_urea_list[1]/10**6)
    # plt.show()
