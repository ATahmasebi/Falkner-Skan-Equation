# -*- coding: utf-8 -*-
"""
Created on Wed Nov  8 22:45:26 2017

@author: Dan
"""

# Exercise 2.13 of Hastings-McLeod:
# Show that if beta > 0 then the boundary value problem
#   f"' + f f" + beta(1 - f'^2) =0
#   f(0) = f'(0) = 0, f'(infinity) = 1
# has a solution such that f' > 0 on (0,infty).
# Hint: Suppose that f"(0) = alpha.
# Vary alpha. Base your shooting sets on the behaviors of f' and f". You may
# wish to use a numerical ode solver to examine f' for different values of alpha
#
# Falkner-Skan equation

import numpy as np
import matplotlib.pyplot as plt
import os.path

def fs(inter,ic,n):
        [a,b] = inter; h = float(b-a)/n         # plot n points in total
        y = np.empty((n+1,3)); t = np.empty(n+1)
        y[0,:] = ic                                                     # enter initial conditions in y
        print (' ic=',ic)
        t[0] = a
        for i in range(0,n):
                t[i+1] = t[i]+h
                y[i+1,:] = rk4step(t[i],y[i,:],h)
        fig = plt.figure()
        for k in range(0,3):
                axk=fig.add_subplot(3,1,k+1)
                axk.plot(t,y[:,k])
                plt.grid(True);
#               ykmax=max(y[:,k]); ykmin=min(y[:,k]) # ok if not animating
                if(k == 0):  ykmin=0.;  ykmax=5.; varlabel="f"
                if(k == 1):  ykmin=-0.2;  ykmax=1.4; varlabel="f '"
                if(k == 2):  ykmin=-0.2;  ykmax=1.4; varlabel='f "'
                axk.axis([0,b,ykmin,ykmax])
                axk.set_ylabel(varlabel)
                if(k == 0):
                        plt.title("Falkner-Skan f\"' + f f\" + (1-(f')^2) = 0\n\
f(0) = f'(0)=0, f\"(0) = {0:.3f}".format(alpha))
        fig.savefig('falknerskan_{0:02d}'.format(fignum))
        return y

def rk4step(t,w,h):     # one step of Runge-Kutta order 4 method
        s1 = ydot(t,w)
        s2 = ydot(t+h/2,w+h*s1/2.)
        s3 = ydot(t+h/2,w+h*s2/2.)
        s4 = ydot(t+h,  w+h*s3)
        return w + h*(s1+2*s2+2*s3+s4)/6.

def ydot(t,w): # [w[0],w[1],w[2]]=[f,f',f"]
        z = np.zeros(3)
        z[0] = w[1]
        z[1] = w[2]
#   f"' = -( f f" + beta(1 - f'^2) )
        z[2] = -w[0]*w[2] - beta*(1-w[1]**2) # = f"'
        return z
os.system("rm falknerskan.gif")
os.system("rm falknerskan_??.png")
for fignum in range(0,51): # will write falknerskan_00.png, .., falknerskan_50.png
        alpha=1.2+0.001*fignum; beta=1
        y = fs([0,5],[0.,0.,alpha],200)
os.system("convert falknerskan_??.png falknerskan.gif") # create multigif animation