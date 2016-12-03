#  seismic_prototype.py
#  Omkar H. Ramachandran
#  omkar.ramachandran@colorado.edu
#
#  Simulation of a 2D seismic problem. This is a python prototype that will be 
#  rewritten in C and parallelized with MPI and CUDA
#
#  This project is licensed under  For more information about the 
#  license go online to <http://www.gnu.org/licenses>. Any reproduction of code
#  in this file MUST contain this header in its entirety.
#
from matplotlib.pylab import *
import time
pi = 3.14159265
# Only using square domains for now. Will change to include both Nx and Nz
N = 200
Tmax = 40.0
Xmax = 10000
vs = 500.0
vp = 1000.0
rho = 1000.0
# Use vs to determine Courant condition
dx = Xmax/N
dt = 1.0/sqrt(2)*dx/vs
# Define time vector
t = linspace(0,Tmax,Tmax/dt)
# Define X and Y vector for defining forcing function
x = linspace(0,Xmax,N)
# define displacement vectors
Nstep = 1
Nsave = int((Tmax/dt)/Nstep)+1
vnew = zeros([N,N])
v = zeros([N,N])
vold = zeros([N,N])
wnew = zeros([N,N])
w = zeros([N,N])
wold = zeros([N,N])
Uxvals = zeros([Nsave,N,N])
Uzvals = zeros([Nsave,N,N])
# define source
nsave = 0
T0 = 3
omega = 1
Xc = 15
Yc = 15
Xl = 1000
Yl = 1000
# Define Ricker time function
Jsrc = sin(0.5*t)
# Define gaussian forcing function
X,Y = meshgrid(x,x)
Jshape = 10*exp(-(X-x[N/2])**2/Xl**2-(Y-x[N/2])**2/Yl**2)
#fig2 = figure(2)
#contourf(x,x,Jshape)
#show()
# Forcing function f = Jsrc[t]*Jshape
Nt=int(Tmax/dt)
def grad_div(v,w):
	grad_div = zeros([N,N,2])
	grad_div[1:-1,1:-1,0] = ((v[2:,1:-1]-v[1:-1,1:-1])/dx-(v[1:-1,1:-1]-v[:-2,:-2])/dx)/dx \
+((w[2:,2:]-w[2:,:-2])/(2*dx)+(w[:-2,2:]-w[:-2,:-2])/(2*dx))/(2*dx)
	grad_div[1:-1,1:-1,1] = ((v[1:-1,2:]-v[1:-1,1:-1])/dx-(v[1:-1,1:-1]-v[1:-1,:-2])/dx)/dx \
+((w[2:,2:]-w[:-2,2:])/(2*dx)+(w[2:,:-2]-w[:-2,:-2])/(2*dx))/(2*dx)
	return grad_div
def delsq(v,w):
	delsq = zeros([N,N,2])
	delsq[1:-1,1:-1,0] = (v[2:,1:-1]-v[:-2,1:-1])/dx**2 + (v[1:-1,2:]-v[1:-1,:-2])/dx**2
	delsq[1:-1,1:-1,1] = (w[2:,1:-1]-w[:-2,1:-1])/dx**2 + (w[1:-1,2:]-w[1:-1,:-2])/dx**2
	return delsq
def ucalc(v,vnew,vold,u,unew,uold):
	unew[1:-1,1:-1] = dt**2/rho*(grad_div(v,w)[1:-1,1:-1,0]*(vs**2+vp**2)-vs**2*delsq(v,w)[1:-1,1:-1,0])+2*u[1:-1,1:-1]-uold[1:-1,1:-1]
	vnew[1:-1,1:-1] = dt**2/rho*(grad_div(v,w)[1:-1,1:-1,1]*(vs**2+vp**2)-vs**2*delsq(v,w)[1:-1,1:-1,1])+2*v[1:-1,1:-1]-vold[1:-1,1:-1]
#	print(unew.max())	
for n in range(Nt):
	ucalc(v,vnew,vold,w,wnew,wold)
	v = vnew
	vold = v
	w = wnew + Jsrc[n]*Jshape
	wold = w
	if (n%Nstep==0):
#		fig1 = figure(1)
#		contourf(x,x,log10(w[:,:]))
#		show()
		Uxvals[nsave,:,:] = v
		Uzvals[nsave,:,:] = w
		nsave += 1
		#print(Uxvals[nsave,:,:].max())
vals = linspace(-4,6,21)
#fig1 = figure(1)
#print t[1*Nstep]
#contourf(x,x,Uzvals[2,:,:])
#colorbar()
#show()
y = Uzvals[:,120,100]
xx = linspace(0,Tmax,Nsave)
plot(xx,y)
show()
#for i in range(Nsave):
#	figure(i)
#	contourf(x,x,log10(Uzvals[i,:,:]))
#	time.sleep(2)
#	colorbar()
#	show()
