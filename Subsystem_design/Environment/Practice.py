import matplotlib.pyplot as plt
import numpy as np
import scipy as sp

def f(x):
  A=.5
  x0=0.7
  if -x0<x<x0:
    return(np.cos(x)+x)
  else:
    return(0)

def g(x):
  A=0.5
  x0=0.7

  if x<x0:
    return(A)
  else:
    return(0)

def convo(x):
  fog=0
  xt=-3
  dx=0.01
  while xt<3:
    fog=fog+f(xt)*g(x-xt)*dx
    xt=xt+dx
  return(fog)


def log(x):
  fx = np.ones(len(x))
  fx[x<0.5]= 0
  fx[x>1.5] = 0

  return fx
def step(x):
  gx = np.ones(len(x))
  gx[x<0.5] = 0
  gx[x>1.5] = 0

  return gx
t = np.arange(-1,6,0.01)

x = np.arange(-1,6,0.01)
x[(t>0)&(t<2)]= np.linspace(0,5,len(x[(t>0)&(t<2)]))
print('x',x)
output = np.convolve(step(x),log(x),'full')/100

# print(output[:step(x).size])
# x2 = np.linspace(-1,6,len(output))
# x2[(x<0)&(x>2)]= 5
# print(x2)
plt.plot(x,log(x),label='Exp of x')
plt.plot(x,step(x),label='step function')
plt.plot(x2,output,label='convolution')
plt.legend()
plt.show()
# #g*f = int_0^t f(x)g(t-x)dx
# dt=0.005
# t= np.arange(-2.,2.+dt,dt)
#
# f_t = []
# g_t = []
# conv = []
# for time in t:
#   f_t.append(f(time))
#   g_t.append(g(time))
#   conv.append(convo(time))
#
#
# plt.plot(t,f_t,label='f(t)')
# plt.plot(t,g_t,label='g(t)')
# plt.plot(t,conv,label='convolution(t)')
# plt.legend()
# plt.show()