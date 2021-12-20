import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt


N = 15000# Population
# initial number of infected and recovered individuals  
# s->suceptible
# e->exposed
# i->infected
# q->isolated
# r->recovered

e0=2/N
s0=1-e0
i0=0
q0=0
r0=0
x0 = [s0,e0,i0,q0,r0]



tsie=0.00004     # susceptible->(infected+exposed)
tei= 0.05*1    # exposed->infected
for i in range(0,8):

    teq=0.8-0.1*i  # exposed->isolation
    tiq=0.08      # infected->isolated
    thetha=0.08     # isolation->recovery
    tndd=0.01*3     # (natuaral+diseaseed)_death rate
    tir=0.006*3
    ter=0.03


    def covid(x,t):
        s,e,i,q,r = x
        dx = np.zeros(5)
        dx[0] = tndd-tndd*s-tsie*N*s*(e+i)# dS/dt equation
        dx[1] = tsie*N*s*(e+i)-tei*e-(tndd+teq)*e-ter*e# dE/dt equation
        dx[2] = tei*e-tiq*i-tndd*i-tir*i # dI/dt equation
        dx[3] = teq*e+tiq*i-thetha*q-tndd*q# dQ/dt equation
        dx[4]=thetha*q-tndd*r +ter*e +tir*i  # dR/dt
        return dx

    t = np.linspace(0, 200, 1000)
    x = odeint(covid,x0,t)
    s = x[:,0]; e = x[:,1]; i = x[:,2]; r = x[:,4]

    # plot the data
    plt.figure(figsize=(8,5))

    plt.subplot(2,1,1)
    plt.title('Ro = '+str('''tsie*N*(tiq+tndd+tei)/((tei+tndd+teq)*(tiq+tndd)  ) ''')+ ', teq='+str(teq)  )
    plt.plot(t,s, color='blue', lw=3, label='Susceptible')
    plt.plot(t,r, color='red',  lw=3, label='Recovered')
    plt.ylabel('Fraction')
    plt.legend()

    plt.subplot(2,1,2)
    plt.plot(t,i, color='orange', lw=3, label='Infective')
    plt.plot(t,e, color='purple', lw=3, label='Exposed')
    plt.ylim(0, 0.8)
    plt.xlabel('Time (days)')
    plt.ylabel('Fraction')
    plt.legend()

    plt.show()