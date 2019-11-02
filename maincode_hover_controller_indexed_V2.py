import numpy as np
import matplotlib.pyplot as plt
from functions import euler_to_rotm as rotm
from RK4 import RK4 as RK4
import time

m = 2.0 #mass of the quad in kg
g = 9.81 #gravitational acceleration

Ixx = 1.0
Iyy = 1.0
Izz = 1.0

#I = np.array([[Ixx, 0, 0], [0, Iyy, 0], [0, 0, Izz]])

#freq of position control loop is 100Hz
#freq of attitude control loop is 1000Hz
freq_pos = 100  #100Hz is the freq of position control loop
T = 1 #simulation time in seconds
N = T*freq_pos #Simulation for 5 sec 

loop_ratio = 10

r_dot_dot_T = np.zeros((3,N))           #acceleration to be traked at all times
r_dot_T     = np.zeros((3,N))           #velocity to be traked at all times
r_T         = np.zeros((3,N))           #position to be tracked at all times
r_T[2,:]    = 0                         #tracking z coorrdinate is always 0, so are x and y

r_dot_dot_des = np.zeros(3)             #desired acceleration at that moment
r_dot_des = np.zeros(3)                 #desired velocity at that moment
r_des = np.zeros(3)                     #desired position at that moment


r               = np.zeros((3,N+1))       #position at all times        
r[:,0]          = np.array([1, 0, 0])   #initialising r at t=0

r_dot           = np.zeros((3,N+1))       #velocity at all times
r_dot[:,0]      = np.array([0, 0, 0])   #initialising r_dot at t=0

r_dot_dot       = np.zeros((3,N+1))       #acceleration at all times
r_dot_dot[:,0]  = np.array([0, 0, 0])   #initialising r_dot_dot at t=0

k_phi           = 0.01*np.array([0.5, 0.8])                # k_phi = [kp_phi, kd_phi]
k_theta         = 0.01*np.array([0.1, 0.7])              # k_theta = [kp_theta, kd_theta]
k_psi           = 0.01*np.array([0.5, 0.5])                # k_psi = [kp_psi, kd_psi]
k_r             = 0.01*np.array([0.5, 0.5, 0.5])               # k_r = [kp_r, kd_r, ki_r]

att_T = np.zeros((3, 10*N))             # [phi_T, theta_T, psi_T]    atttude trajectory to be  tracked
# for i in range(10*N):
#     att_T[0, i] = 0                     #phi_T = 0 at all times
#     att_T[1, i] = 0                     #theta_T = 0 at all times
#     att_T[2, i] = 0                     #psi_T = 0 at all times
    
att             = np.zeros((3, 10*N))               #attitude [phi, theta, psi] at all times
att[:, 0]       = np.zeros(3)#np.array([0.1, 0.1, 0])         #initialising attitude at t=0
des_att         = np.zeros((3, 10*N+1))                   #desired attitude in that loop [phi, theta,  psi]
des_att_rates   = np.zeros((3, 10*N+1))             #desired body angular rates in that loop [p_des, q_des, r_des]
omega_b         = np.zeros((3,10*N+1))            #[p, q, r] angular rates at each time
omega_b[:,0]    = np.array([0,0,0])        #np.array([5*np.pi/180, 4*np.pi/180, 3*np.pi/180])    initialize omega_b
control_inputs  = np.zeros((4, 10*N))    #[u1, u2, u3, u4]
control_inputs[0,0] = m*g 

time_att = np.linspace(0, T, N*10+1)
time_att_plt = time_att[0:(N*10)]
time_pos = np.linspace(0, T, N+1)
time_pos_plt = time_pos[0:N]

start = time.time()
# print(r_dot_dot_T)

for i in range(0,N):  #position control loop
    
    #print(i)
    r_dot_dot_des = k_r[0]*(-r[:,i])  +   k_r[1]*(-r_dot[:, i])    +   k_r[2] *(-r[:, i])    #PID for des acceleratiom
    print(r)
    # print(r_dot[2,i])
    #print(r_dot_dot_des[2]+g)
    control_inputs[0, i] = m*(r_dot_dot_des[2]+g)

    for j in range(loop_ratio):   #attitude control loop
        #print("j=", j)
        if (10*i+j == 10*N-1):             #If statement to prevent the index error for the last index in propagation step
            break
        else:
            des_att[0, 10*i+j] = (r_dot_dot_des[0]*np.sin(att_T[2, 10*i+j]) - r_dot_dot_des[1]*np.cos(att_T[2, 10*i+j]))/g
            des_att[1, 10*i+j] = (r_dot_dot_des[0]*np.cos(att_T[2, 10*i+j]) + r_dot_dot_des[1]*np.sin(att_T[2, 10*i+j]))/g
            print(des_att[1, 10*i+j])
            control_inputs[1, 10*i+j] = k_phi[0]   * (des_att[0, 10*i+j] - att[0, 10*i+j])    +   k_phi[1]   * (des_att_rates[0, 10*i+j] - omega_b[0, 10*i+j])
            control_inputs[2, 10*i+j] = k_theta[0] * (des_att[1, 10*i+j] - att[1, 10*i+j])    +   k_theta[1] * (des_att_rates[1, 10*i+j] - omega_b[1, 10*i+j])
            control_inputs[3, 10*i+j] = k_psi[0]   * (des_att[2, 10*i+j] - att[2, 10*i+j])    +   k_psi[1]   * (des_att[2, 10*i+j] - att[2, 10*i+j])

            def dynamics(t,y):
                p_dot = control_inputs[1,10*i+j]/Ixx - y[1]*y[2]*(Izz-Iyy)/Ixx
                q_dot = control_inputs[2,10*i+j]/Iyy - y[0]*y[2]*(Ixx-Izz)/Iyy
                r_dot = control_inputs[3,10*i+j]/Izz
                return np.array([p_dot, q_dot, r_dot])

            omega_b[:, 10*i+j+1] = RK4(dynamics, omega_b[:, 10*i+j], time_att[10*i+j], time_att[10*i+j+1] , 0.0001/100, 3)[0]        # Propagation of omega_b
            att[:, 10*i+j+1]     = att[:, 10*i+j] + omega_b[:, 10*i+j]*0.0001     # is the time in sec for which att control loop runs

            if not(i ==0 and j==0):
                des_att_rates[:, 10*i+j+1] = 0.001*(des_att[:, 10*i+j]-des_att[:, 10*i+j-1])   # calculating phi_dot_des and so on using phi_des

    R = rotm(att[0,10*i], att[1, 10*i], att[2, 10*i])
    vec = np.array([0, 0, control_inputs[0, i]])
    r_dot_dot[:, i+1] = np.array([0, 0, -g]) + (1/m)*(R.dot(vec))
    r_dot[:, i+1] = r_dot[:, i] + r_dot_dot[:, i]*0.01#*0.001
    r[:, i+1] = r[:, i] + r_dot[:, i]*0.01 + 0.5*r_dot_dot[:, i]*(0.01**2)#*(0.001)**2
    # print(vec)


end=time.time()
print('time = ',end-start )

omega_b = omega_b[:, 0:10*N]
r = r[:, 0:N]
r_dot = r_dot[:, 0:N]
r_dot_dot = r_dot_dot[:, 0:N]

################### CSV ###############################
np.savetxt("omega_b.csv", omega_b, delimiter=",")
np.savetxt("r_dot_dot.csv", r_dot_dot, delimiter=",")
np.savetxt("r_dot.csv", r_dot, delimiter=",")
np.savetxt("r.csv", r, delimiter=",")
np.savetxt("control_inputs.csv", control_inputs, delimiter=",")
np.savetxt("r_dot_dot_des.csv", r_dot_dot, delimiter=",")
np.savetxt("r_dot_des.csv", r_dot, delimiter=",")
np.savetxt("r_des.csv", r, delimiter=",")


################## PLOTS ##############################
fig = plt.figure()
fig, ax = plt.subplots(2)
plt.xlabel('time')
plt.ylabel('omega_b')
ax[0].plot(time_att_plt, omega_b[0,:],'r')
ax[0].plot(time_att_plt, omega_b[1,:],'b')
ax[0].plot(time_att_plt, omega_b[2,:],'g')
plt.xlabel('time')
plt.ylabel('z(t)')
ax[1].plot(time_pos_plt, r[2,:])
plt.show()


#-----------------------------------------------------------------------------------------------------------------------------------


            # def dynamics(t,y):
            #     p = control_inputs[1,i]/Ixx - omega_b[1,i]*omega_b[2,10*i+j]*(Izz-Iyy)
            #     q = control_inputs[2,i]/Iyy - omega_b[0,i]*omega_b[2,10*i+j]*(Ixx-Izz)
            #     r = control_inputs[3,i]/Izz
            #     return np.array([p, q, r])
            # control_inputs[1:3, 10*i+j] = 0*control_inputs[1:3, 10*i+j]