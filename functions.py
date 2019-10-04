import numpy as np
import matplotlib.pyplot as plt

def euler_to_rotm(phi, theta, psi):
    R = np.identity(3)
    R[0,0]= np.cos(theta)*np.cos(psi)-np.sin(psi)*np.sin(theta)*np.sin(phi)
    R[0,1]= -np.cos(phi)*np.sin(psi)
    R[0,2]= np.cos(psi)*np.sin(theta)+np.cos(theta)*np.sin(phi)*np.sin(psi)
    R[1,0]= np.cos(theta)*np.sin(psi)+np.cos(psi)*np.sin(theta)*np.sin(phi)
    R[1,1]= np.cos(phi)*np.cos(psi)
    R[1,2]= np.sin(psi)*np.sin(theta)-np.cos(psi)*np.cos(theta)*np.sin(phi)
    R[2,0]= -np.cos(phi)*np.sin(theta)
    R[2,1]= np.sin(phi)
    R[2,2]= np.cos(phi)*np.cos(theta)
    R_new = np.array([[R[0,0], R[0,1], R[0,2]], [R[1,0], R[1,1], R[1,2]], [R[2,0], R[2,1], R[2,2]]])
    return R_new

def euler_rates_to_omega(phi, theta, psi, phi_dot, theta_dot, psi_dot):
    omega=np.zeros(3)
    M=np.identity(3)
    M[0,0]= np.cos(theta)
    M[0,1]= 0
    M[0,2]= -np.cos(phi)*np.sin(theta)
    M[1,0]= 0
    M[1,1]= 1
    M[1,2]= np.sin(phi)
    M[2,0]= np.sin(theta)
    M[2,1]= 0
    M[2,2]= np.cos(phi)*np.cos(theta)
    M_new = np.array([[M[0,0], M[0,1], M[0,2]], [M[1,0], M[1,1], M[1,2]], [M[2,0], M[2,1], M[2,2]]])
    euler_rates = np.array([phi_dot, theta_dot, psi_dot])
    omega = np.matmul(M_new, euler_rates)
    return omega

#v_e_p = {[(r_T-r).n]n + [(r_T-r).b]b}  v_e_p is the error vector of the position
#e_v = r_T_dot - r_dot
