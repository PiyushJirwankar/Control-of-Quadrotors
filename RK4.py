import numpy as np

#########################################################
### y_dot = f(t, y);    y(t0) = y0
#########################################################
#h is the time step between two instants
#n is the dimension of y 

# y is an (n,M) shaped matrix, with y_i at time i in the ith column
def RK4(f, y0, t0, tf, h, n):        #h is the time step between two instants
    N=int((tf-t0)/h)
    k1=np.array(np.zeros(n))
    k2=np.array(np.zeros(n))
    k3=np.array(np.zeros(n))
    k4=np.array(np.zeros(n))
    time=np.array(np.zeros(N))
    time[0]=t0
    y=np.array(np.zeros((n,N)))
    y[:, 0] = y0
    
    for j in range(N-1):
        for i in range(n):
            k1[i]=float(h*f(time[j],y[:,j])[i])
            k2[i]=float(h*f(time[j]+h/2,y[:,j]+k1/2)[i])
            k3[i]=float(h*f(time[j]+h/2,y[:,j]+k2/2)[i])
            k4[i]=float(h*f(time[j]+h,y[:,j]+k3)[i])
        y[:,j+1]=y[:,j]+(k1+2*(k2+k3)+k4)/6
        time[j+1]=time[j]+h
    return y, time

#time[0]=t0
#time[j+1]=time[j]+h
#       if j==N-1:
#           time[j+1]=time[j]