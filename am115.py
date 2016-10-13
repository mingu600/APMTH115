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
np.set_printoptions(suppress=True)

# SPECIFY DIRECTION : TRUE FOR REVERSE FLOW | FALSE FOR FORWARD FLOW
reverse_direction = True

# PARAMETERS
time_steps = 10**2 #each time step advances blood/dialysis in tubings
length_segments = 10**2 #divide blood/dialysis tubing into compartments
pt_blood_urea_init = 10**4 #initial urea molecules in patient's blood
frac_pt_blood_in_dial = .1 #i.e 10% of patient's blood volume fits in dialysizer
blood_velocity = 1 #positive int
dialysis_velocity = 1 #positive int
diffusion_constant = 0.3 #diffusion proportional to difference in concentration

# initialize
blood = np.zeros(length_segments)
dialysis = np.zeros(length_segments)
diffused_count = 0
time_steps_count = 0
urea_cleared_flag = False
blood_list = np.zeros((time_steps, length_segments))
dialysis_list = np.zeros((time_steps, length_segments))

for t in range(time_steps):
    time_steps_count += 1
    curr_pt_urea = pt_blood_urea_init - diffused_count
    urea_step = curr_pt_urea * frac_pt_blood_in_dial / length_segments
    if urea_step < 1: #all urea cleared from patient
        urea_cleared_flag = True
        break
    blood = np.roll(blood, blood_velocity)
    blood[:blood_velocity] = urea_step
    if reverse_direction:
        dialysis = np.roll(dialysis, -dialysis_velocity)
        dialysis[dialysis_velocity:] = 0
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

# print results
print "Dialysis Fluid:"
print dialysis
print "\nBlood Fluid:"
print blood
print "\nTime steps passed: %i" % time_steps_count
print "Amount of urea diffused: %i" % diffused_count
if urea_cleared_flag:
    print "ALL UREA CLEARED FROM PATIENT"
else:
    urea_remaining = pt_blood_urea_init - diffused_count
    print "Urea remaining in patient: %i" % urea_remaining
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
