import numpy as np
import matplotlib.pyplot as plt
from functions import euler_to_rotm as rotm
from RK4 import RK4 as RK4

m = 2 #mass of the quad in kg
g = 9.81 #gravitational acceleration

Ixx = 1
Iyy = 1
Izz = 1

#I = np.array([[Ixx, 0, 0], [0, Iyy, 0], [0, 0, Izz]])

#freq of position control loop is 100Hz
#freq of attitude control loop is 1000Hz
N = 5*100 #Simulation for 5 sec 

r_dot_dot_T = np.zeros((3,N))           #acceleration to be traked at all times
r_dot_T     = np.zeros((3,N))           #velocity to be traked at all times
r_T         = np.zeros((3,N))           #position to be tracked at all times
for i in range(N):
    r_T[2,i] = 2                        #z coordinate is 2m; rest all 0 
    r_T[1,i] = 0
    r_T[0,i] = 0

r_dot_dot_des = np.zeros(3)             #desired acceleratio at that moment
r_dot_des = np.zeros(3)                 #desired velocity at that moment
r_des = np.zeros(3)                     #desired position at that moment

init_att        = np.zeros(3)           #initial values of [phi, theta, psi]
init_att_rates  = np.zeros(3)           #initial values of [p, q, r]

r               = np.zeros((3,N))       #position at all times        
r[:,0]          = np.array([0, 0, 0])   #initialising r at t=0

r_dot           = np.zeros((3,N))       #velocity at all times
r_dot[:,0]      = np.array([0, 0, 0])   #initialising r_dot at t=0

r_dot_dot       = np.zeros((3,N))       #acceleration at all times
r_dot_dot[:,0]  = np.array([0, 0, 0])   #initialising r_dot_dot at t=0

k_phi = np.array([0.5, 0.8])                # k_phi = [kp_phi, kd_phi]
k_theta = np.array([0.1, 0.7])              # k_theta = [kp_theta, kd_theta]
k_psi = np.array([0.5, 0.5])                # k_psi = [kp_psi, kd_psi]
k_r = np.array([0.5, 0.5, 0.5])               # k_r = [kp_r, kd_r, ki_r]

att_T = np.zeros((3, 10*N))             # [phi_T, theta_T, psi_T]    atttude trajectory to be  tracked
for i in range(10*N):
    att_T[0, i] = 0                     #phi_T = 0 at all times
    att_T[1, i] = 0                     #theta_T = 0 at all times
    att_T[2, i] = 0                     #psi_T = 0 at all times
    
att = np.zeros((3, 10*N))               #attitude [phi, theta, psi] at all times
att[:, 0] = np.array([0, 0, 0])         #initialising attitude at t=0
des_att = np.zeros(3)                   #desired attitude in that loop [phi, theta,  psi]
des_att_rates = np.zeros(3)             #desired body angular rates in that loop [p_des, q_des, r_des]
omega_b = np.zeros((3,10*N))               #[p, q, r] angular rates at each time
omega_b[:,0] = np.array([5*np.pi/180, 4*np.pi/180, 3*np.pi/180])
control_inputs = np.zeros((4, 10*N))    #[u1, u2, u3, u4]

def dynamics(t, y):
    p = control_inputs[0]/Ixx - omega_b[1]*omega_b[2]*(Izz-Iyy)
    q = control_inputs[1]/Iyy - omega_b[0]*omega_b[2]*(Ixx-Izz)
    r = control_inputs[2]/Izz
    return np.array([p, q, r])


time_att = np.linspace(0, 5, 5001)
time_att = time_att[0:5000]
time_pos = np.linspace(0, 5, 501)
time_pos = time_pos[0:500]

for i in range(N):  #position control loop
    r_dot_dot_des = r_dot_dot_T[:, i] + k_r[0]*(r_T[:, i]-r[:,i]) + k_r[1]*(r_dot_T[:, i]-r_dot[:, i]) + k_r[2]*((r_T[:, i]-r[:, i])*0.01)
    for j in range(10):
        des_att[0] = (r_dot_dot_des[0]*np.sin(att_T[2, 10*i+j]) - r_dot_dot_des[1]*np.cos(att_T[2, 10*i+j]))/g
        des_att[1] = (r_dot_dot_des[0]*np.cos(att_T[2, 10*i+j]) + r_dot_dot_des[1]*np.sin(att_T[2, 10*i+j]))/g
        control_inputs[0, 10*i+j] = m*r_dot_dot_des[2]
        control_inputs[1, 10*i+j] = k_phi[0]*(des_att[0]-att[0, 10*i+j]) + k_phi[1]*(des_att_rates[0] - omega_b[0, 10*i+j])
        control_inputs[2, 10*i+j] = k_theta[0]*(des_att[1]-att[1, 10*i+j]) + k_theta[1]*(des_att_rates[1] - omega_b[1, 10*i+j])
        control_inputs[3, 10*i+j] = k_psi[0]*(des_att[2]-att[2, 10*i+j]) + k_psi[1]*(des_att_rates[2] - omega_b[2, 10*i+j])
        if (i==N-1 and j==9):
            continue
#        else:
#            omega_b[:, 10*i+j+1] = RK4(dynamics, omega_b[:, 10*i + j], time_att(10*i+j), time_att(10*i+j+1) ,0.0001/100 , 3) #h=0.0001/100
        att[:, 10*i+j+1] = att[:, 10*i+j] + omega_b[:, 10*i+j]*0.0001
    R = rotm(att[0,10*i], att[1, 10*i], att[2, 10*i])
    r_dot_dot[:, i+1] = np.array([0, 0, g]) + (1/m)*np.matmul(R, np.array([0, 0, control_inputs[0]]))
    r_dot[:, i+1] = r_dot[:, i] + r_dot_dot[:, i]*.001
    r[:, i+1] = r[:, i] + r_dot[:, i]*0.001 + 0.5*r_dot_dot[:, i]*(0.001)**2


#plt.plot(time_pos, r_dot_dot[2,:])
print(R)
#plt.plot(time_att, att[2,:]*180/np.pi)
#plt.xlim(0, 0.001)
#plt.plot(time_pos, r[2,:])