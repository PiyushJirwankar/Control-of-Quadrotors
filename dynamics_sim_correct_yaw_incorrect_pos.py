import numpy as np
import matplotlib.pyplot as plt
from RK4 import RK4 as RK4
import functions

test = 1
g = 9.81
m = 2.0

T = 25	 #simulation time
N_pos = 100 * T    #number of sample points for outer loop
N_att = 10*N_pos   #number of sample points for inner loop
# print(N_att)
time_att = np.linspace(0, T, N_att, endpoint = False)
time_pos = np.linspace(0, T, N_pos, endpoint = False)

##### omega #############
omega = np.zeros((3, N_att+1))
omega[:, 0] = np.array([1.0, 1.0, 1.0])

##### euler angles ########
euler_angle = np.zeros((3, N_att))
euler_angle[:, 0] = (1/180.0)*(np.pi)*np.array([5, 6, 7])

psi_T = np.pi/(180.0)*60 
# print(psi_T) 

time_control = 5 #sec
psi_T_arr = np.linspace(euler_angle[0, 2], psi_T, N_att)


###### Moment of inertia #######
I = np.array([[1.0, 0.0, 0.0],
			  [0.0, 1.0, 0.0],
			  [0.0, 0.0, 2.0]])

Ixx = I[0,0]
Iyy = I[1,1]
Izz = I[2,2]


r_T = np.zeros((3, N_pos))
for i in range(N_pos):
	r_T[2] = 1

r = np.zeros((3, N_pos+1))
r[:, 0] = np.array([0.0, 0.0, 0.0])

r_dot = np.zeros((3, N_pos+1))
r_dot[:, 0] = np.array([0.0, 0.0, 0.0])

r_ddot = np.zeros((3, N_pos+1))
r_ddot[:, 0] = np.array([0.0, 0.0, 0.0])

r_des 		= np.zeros(3)
r_des_dot 	= np.zeros(3)
r_des_ddot 	= np.zeros(3)

control_inputs = np.zeros((4, N_att)) # [u1, u2, u3, u4]
des_euler = np.zeros((3, N_att))  # [phi, theta, psi]
des_omega = np.zeros((3, N_att))  # [p_des, q_des, r_des]

kp_pos = 10*np.array([1.0, 1.0, 1.0])
kd_pos = 1*np.array([1.0, 1.0, 10.0])
ki_pos = 1*np.array([1.0, 1.0, 1.0])

kp_ang = np.array([1.0, 1.0, 1.0])  # [kp_phi, kp_theta, kp_psi]
kd_ang = np.array([1.0, 1.0, 2.0])  # [kd_phi, kd_theta, kd_psi]

ei = np.zeros(3)
ei_ang = np.zeros(3)

# print(psi_T_arr)
for i in range(N_pos):
	print(i)
	ei = ei + (r_T[:, i] - r[:, i])*0.01#*(1/(1.0*N_pos*T))
	r_des_ddot[:] = -kd_pos*(r_dot[:, i]) + kp_pos*(r_T[:, i]-r[:, i]) + ki_pos*ei 
	print(r_dot[:,i])#r_des_ddot[:])

	for j in range(10):
		# print("j = ", j)

		if 10*i+j==N_att-1:
			break
		else:
			# print(control_inputs[1:, 10*i+j])
			des_euler[2, 10*i+j] = psi_T
			control_inputs[1:, 10*i+j] = kp_ang*(des_euler[:, 10*i+j] - euler_angle[:, 10*i+j]) + kd_ang*(des_omega[:, 10*i+j] - omega[:, 10*i+j])
			control_inputs[0, 10*i+j] = m*r_des_ddot[2]

			def rigid_body_dyn(t, w):  #returns w_dot	
				# w_dot = np.linalg.inv(I).dot(control_inputs[1:, 10*i+j] - np.cross(w, I.dot(w)))
				w_dot = np.zeros(3)
				w_dot[0] = control_inputs[1,10*i+j]/Ixx - w[1]*w[2]*(Izz-Iyy)/Ixx
				w_dot[1] = control_inputs[2,10*i+j]/Iyy - w[0]*w[2]*(Ixx-Izz)/Iyy
				w_dot[2] = control_inputs[3,10*i+j]/Izz
				return w_dot

			des_euler[0, 10*i+j] = (1/g)*(r_des_ddot[0]*np.sin(psi_T) - r_des_ddot[1]*np.cos(psi_T))
			des_euler[1, 10*i+j] = (1/g)*(r_des_ddot[0]*np.cos(psi_T) + r_des_ddot[1]*np.sin(psi_T))
			#(psi_T_arr[10*i+j] - euler_angle[2, 10*i+j]) #des_euler[2, 10*i+j-1] + omega[2, 10*i+j]*0.001
			# print(des_omega[2, 10*i+j])

			if not(i ==0 and j==0):
				des_omega[:, 10*i+j+1] = (des_euler[:, 10*i+j]-des_euler[:, 10*i+j-1])*0.001#/(1.0*N_att*T)
				# des_omega[2, 10*i+j+1]  = (psi_T_arr[10*i+j] - psi_T_arr[10*i+j-1])/0.001#des_euler[2, 10*i+j] + (omega[2, 10*i+j])*0.001   #calculating phi_dot_des

			omega[:, 10*i+j+1] = RK4(rigid_body_dyn, omega[:, 10*i+j], time_att[10*i+j], time_att[10*i+j+1] , 0.001/10, 3)[0]
			euler_angle[:, 10*i+j+1] = euler_angle[:, 10*i+j] + 0.001*omega[:, 10*i+j]

	R = functions.euler_to_rotm(euler_angle[0, 10*i], euler_angle[1, 10*i], euler_angle[2, 10*i])
	thrust = np.array([0.0, 0.0, control_inputs[0, 10*i]])
	r_ddot[:, i+1] = (1/m)*(np.dot(R, thrust)) - np.array([0.0, 0.0, g])
	r_dot [:, i+1] = r_dot[:, i] + r_ddot[:, i+1]*0.01
	r     [:, i+1] = r[:, i] + r_dot[:, i+1]*0.01
	# print(np.dot(R, thrust))
	# print(r[:, i])

omega 	= omega[:, :-1]
r 		= r[:, :-1]
r_dot 	= r_dot[:, :-1]
r_ddot 	= r_ddot[:, :-1]
################################################################################

# np.savetxt("quad_dynamics_kp1_kd1_ki1.csv", omega, delimiter=",")

fig, ax = plt.subplots(2,2)
ax[0,0].plot(time_att, omega[0,:],'r')
ax[0,0].plot(time_att, omega[1,:],'b')
ax[0,0].plot(time_att, omega[2,:],'g')
ax[0,0].set_xlabel('Time')
ax[0,0].set_ylabel('Angular rates')
ax[0,0].legend(['omega_x', 'omega_y', 'omega_z'], fontsize='xx-small')

# ax[0,0].plot(time_pos, r_ddot[0,:],'r')
# ax[0,0].plot(time_pos, r_ddot[1,:],'b')
# ax[0,0].plot(time_pos, r_ddot[2,:],'g')
# ax[0,0].set_xlabel('Time')
# ax[0,0].set_ylabel('Acceleration')

# ax[0,1].plot(time_pos, r_dot[0,:],'r')
# ax[0,1].plot(time_pos, r_dot[1,:],'b')
# ax[0,1].plot(time_pos, r_dot[2,:],'g')
# ax[0,1].set_xlabel('Time')
# ax[0,1].set_ylabel('velocity')


ax[0,1].plot(time_att, 180.0*(1/np.pi)*euler_angle[0,:], 'r')
ax[0,1].plot(time_att, 180.0*(1/np.pi)*euler_angle[1,:], 'b')
ax[0,1].plot(time_att, 180.0*(1/np.pi)*euler_angle[2,:], 'g')
ax[0,1].set_xlabel('Time')
ax[0,1].set_ylabel('Euler angles (deg)')
ax[0,1].legend(['Roll', 'Pitch', 'Yaw'], fontsize='xx-small')

ax[1,1].plot(time_pos, r[0,:], 'r')
ax[1,1].plot(time_pos, r[1,:], 'b')
ax[1,1].plot(time_pos, r[2,:], 'g')
ax[1,1].set_xlabel('Time (sec)')
ax[1,1].set_ylabel('Position (m)')
ax[1,1].legend(['X', 'Y', 'Z'], fontsize='xx-small')

ax[1,0].plot(time_att, control_inputs[0,:], 'r')
ax[1,0].plot(time_att, control_inputs[1,:], 'b')
ax[1,0].plot(time_att, control_inputs[2,:], 'g')
ax[1,0].plot(time_att, control_inputs[3,:], 'k')
ax[1,0].set_xlabel('Time')
ax[1,0].set_ylabel('Control inputs')
ax[1,0].legend(['Thrust (N)', 'Moment_x (N.m)', 'Moment_y (N.m)', 'Moment_z (N.m)'], fontsize='xx-small')

plt.show()