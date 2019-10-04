import numpy as np
#import functions
#import dynamics
from RK4 import RK4 as RK4

import matplotlib.pyplot as plt

m = 2 #mass of the quad in kg
g = 9.81 #gravitational acceleration

Ixx = 1
Iyy = 1
Izz = 1

#I = np.array([[Ixx, 0, 0], [0, Iyy, 0], [0, 0, Izz]])

#freq of position control loop is 100Hz
#freq of attitude control loop is 1000Hz
N = 5*100 #Simulation for 5 sec 

r_dot_dot_des = np.array([0, 0, 0])   #initial desired acceleration
r_dot = np.array([0, 0, 0])    #initializing the velocity
r = np.zeros([0, 0, 0])        #initializing the position

r_dot_dot_T = np.array([0, 0, 0])
r_dot_T = np.array([0, 0, 0])
r_T = np.array([0, 0, 0])

init_att = np.array([0, 0, 0])
init_att_rates = np.array([0, 0, 0]) #for initial value in RK4
init_pos = np.array([0, 0, 0])

##############################################################################


k_phi = np.array([1, 1])  # k_phi = [kp_phi, kd_phi]
k_theta = np.array([1, 1]) # k_theta = [kp_theta, kd_theta]
k_psi = np.array([1, 1])   # k_psi = [kp_psi, kd_psi]
k_r = np.array([1, 1, 1]) # k_r = [kp_r, kd_r, ki_r]

des_att = np.array([0, 0, 0])   #attitude in Euler angles in deg
des_att_rates = np.array([0, 0, 0])   #desired body angular rates [p, q, r]
des_pos = np.array([0, 0, 2])   #Position in m
omega_b = np.zeros(3)#np.array([0, 0, 0])   #[p, q, r]
control_inputs = np.array([0, 0, 0])  #[u2, u3, u4]
track_att = np.array([0, 0, 0])  # [phi_T, theta_T, psi_T]    atttude trajectory to be  tracked

omega_b_stored = np.zeros((3, 10*N))
omega_b_stored[:,0]=omega_b

time_att = np.linspace(0, 5, 5001)

def dynamics(t, y):
    a = control_inputs[0]/Ixx - omega_b[1]*omega_b[2]*(Izz-Iyy)
    b = control_inputs[1]/Iyy - omega_b[0]*omega_b[2]*(Ixx-Izz)
    c = control_inputs[2]/Izz
    return np.array([a, b, c])

for i in range(N):      #position control
    ## PID control to find the desired acceleration r_dot_dot_des ################
    r_dot_dot_des = r_dot_dot_T + k_r[0]*(r_T-r) + k_r[1]*(r_dot_T-r_dot) #+ k_r[2]*()       Need RK4 here
    u1_des = m*r_dot_dot_des[2]                       #only for position control
    for j in range(10):
        omega_b=RK4(dynamics, omega_b_stored[:, 10*i + j], time_att(i), time_att(i+1) ,0.0001/100 , 3)
        phi = omega_b[0]*0.001
        theta = omega_b[1]*0.001
        psi = omega_b[2]*0.001
        des_att[0]=(1/g)*(r_dot_dot_des[0]*np.sin(track_att[2]) - r_dot_dot_des[1]*np.cos(track_att[2]))
        des_att[1]=(1/g)*(r_dot_dot_des[0]*np.cos(track_att[2]) + r_dot_dot_des[1]*np.sin(track_att[2]))
        control_inputs[0] = k_phi[0]*(des_att[0]-phi) + k_phi[1](des_att_rates[0]-omega_b[0])         #for attitude as well as position control
        control_inputs[1] = k_theta[0]*(des_att[1]-theta) + k_theta[1](des_att_rates[1]-omega_b[1])         #for attitude as well as position control
        control_inputs[2] = k_psi[0]*(des_att[2]-psi) + k_psi[1](des_att_rates[2]-r)         #for attitude as well as position control
        omega_b_stored[:,10*i + j+1] = omega_b

plt.plot(time_att[::1], omega_b_stored[0,:] )