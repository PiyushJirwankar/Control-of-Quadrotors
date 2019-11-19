import numpy as np
import matplotlib.pyplot as plt
from RK4 import RK4 as RK4

test = 1
g = 9.81
m = 2.0

T = 2  #simulation time
N_pos = 100 * T    #number of sample points for outer loop
N_att = 10*N_pos   #number of sample points for inner loop
time_att = np.linspace(0, T, N_att, endpoint = False)
time_pos = np.linspace(0, T, N_pos, endpoint = False)

##### omega #############
omega = np.zeros((3, N_att))
omega[:, 0] = np.array([0.0, 0.0, 0.0])

##### euler angles ########
euler_angle = np.zeros((3, N_att))
euler_angle[:, 0] = np.array([0.0, 0.0, 0.0])

##### torques ########
torque = np.zeros((3, N_att))
for i in range(N_att): 
	torque[:, i] = np.array([1.0, 0.0, 0.0])

###### Moment of inertia #######
I = np.array([[1.0, 0.0, 0.0],
			  [0.0, 1.0, 0.0],
			  [0.0, 0.0, 1.0]])

r = np.zeros((3, N_pos))
r = np.array([0.0, 0.0, 0.0])

r_dot = np.zeros((3, N_pos))
r_dot = np.array([0.0, 0.0, 0.0])

r_ddot = np.zeros((3, N_pos))
r_ddot = np.array([0.0, 0.0, 0.0])


for i in range(N_pos):

	def rigid_body_dyn(t, w):  #returns w_dot
		w_dot = np.linalg.inv(I).dot(torque[:, i] - np.cross(w, I.dot(w)))
		return w_dot
	
	for j in range(10):
		print(i)
		if 10*i+j==N_att-1:
			break
		else:
			# omega[:, i+1] = RK4(rigid_body_dyn, omega[:, i], Time[i], Time[i+1], 0.1, 3)
			omega[:, 10*i+j+1] = RK4(rigid_body_dyn, omega[:, 10*i+j], time_att[10*i+j], time_att[10*i+j+1] , 0.001/10, 3)[0]
			euler_angle[:, 10*i+j+1] = euler_angle[:, 10*i+j] + 0.001*omega[:, 10*i+j]
	##########

np.savetxt("quad_dynamics%{test}.csv", omega, delimiter=",")

fig, ax = plt.subplots(2,2)
ax[0,0].plot(time_att, omega[0,:],'r')
ax[0,0].plot(time_att, omega[1,:],'b')
ax[0,0].plot(time_att, omega[2,:],'g')
ax[0,0].set_xlabel('Time')
ax[0,0].set_ylabel('Angular rates')
# plt.legend('abc')

ax[0,1].plot(time_att, euler_angle[0,:], 'r')
ax[0,1].plot(time_att, euler_angle[1,:], 'b')
ax[0,1].plot(time_att, euler_angle[2,:], 'g')
ax[0,1].set_xlabel('Time')
ax[0,1].set_ylabel('Euler angles')
# ax[0,1].set_title('Position vs Time')

plt.show()