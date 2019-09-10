# Numerical solution of the diffusion equation using
# different numerical schemes

import numpy as np   # External library for numerical calculations
import matplotlib.pyplot as plt   # Plotting library

# Function defining the initial and analytic solution
def initialBell(x):
    return np.where(x%1. < 0.5, np.power(np.sin(2*x*np.pi),2),0)

# Put everything inside a main function to avoid global variables
def main():
   # Setup parameters
   K = 1     #diffusion constant
   nt = 400000
   nx = 100      # Grid size
   dx = 1./nx   # horizontal steps
   dt = 10**(-3)/nt     # time step in seconds
   # Spatial variable going from zero to one inclusive
   x = np.linspace(0.0,1.0,nx+1)

   # Initial conditions (in meters or m/s)
   phi = initialBell(x)
   phiNew = phi.copy()
   phiOld = phi.copy()

   ### established solution

   # Initial conditions (in meters or m/s)
   phi = initialBell(x)
   phiNew = phi.copy()
   phiOld = phi.copy()

   # FTCS for the first (two) time-step(s), looping over space
   for j in range(1,nx):
      phi[j] = phiOld[j] + K*dt/(dx)**2*(phiOld[j+1] -2*phiOld[j] + phiOld[j-1])
   # apply periodic boundary conditions
   phi[0] = phiOld[0] + K*dt/(dx)**2*(phiOld[1] - 2*phiOld[0] + phiOld[nx-1])
   phi[nx] = phi[0]

   # Loop over remaining time-steps (nt) using CTCS
   for n in range(1,nt):
      # loop over space
      for j in range(1,nx):
         phiNew[j] = phiOld[j] + 2*K*dt/(dx)**2*(phi[j+1] -2*phi[j] + phi[j-1])
      # apply periodic boundary conditions
      phiNew[0] = phiOld[0] + 2*K*dt/(dx)**2*(phi[1] -2*phi[0] + phi[nx-1])
      phiNew[nx] = phiNew[0]
      # update phi for the next time-step
      phiOld = phi.copy()
      phi = phiNew.copy()

   ### Shaun's solution

   # Initial conditions (in meters or m/s)
   phi_shaun = initialBell(x)
   phiNew_shaun = phi_shaun.copy()
   phiOld_shaun = phi_shaun.copy()

   # FTCS for the first (two) time-step(s), looping over space
   for j in range(2,nx-1):
      phi_shaun[j] = phiOld_shaun[j] + K*dt/(2*dx)**2*(phiOld_shaun[j+2] -2*phiOld_shaun[j] + phiOld_shaun[j-2])
   # apply periodic boundary conditions
   phi_shaun[nx-1] = phiOld_shaun[nx-1] + K*dt/(2*dx)**2*(phiOld_shaun[1] - 2*phiOld_shaun[nx-1] + phiOld_shaun[nx-3])
   phi_shaun[0] = phiOld_shaun[0] + K*dt/(2*dx)**2*(phiOld_shaun[2] - 2*phiOld_shaun[0] + phiOld_shaun[nx-2])
   phi_shaun[1] = phiOld_shaun[1] + K*dt/(2*dx)**2*(phiOld_shaun[3] - 2*phiOld_shaun[1] + phiOld_shaun[nx-1])
   phi_shaun[nx] = phi_shaun[0]

   # Loop over remaining time-steps (nt) using CTCS
   for n in range(1,nt):
      # loop over space
      for j in range(2,nx-1):
         phiNew_shaun[j] = phiOld_shaun[j] + 2*K*dt/(2*dx)**2*(phi_shaun[j+2] -2*phi_shaun[j] + phi_shaun[j-2])
      # apply periodic boundary conditions
      phiNew_shaun[nx-1] = phiOld_shaun[nx-1] + 2*K*dt/(2*dx)**2*(phi_shaun[1] -2*phi_shaun[nx-1] + phi_shaun[nx-3])
      phiNew_shaun[0] = phiOld_shaun[0] + 2*K*dt/(2*dx)**2*(phi_shaun[2] -2*phi_shaun[0] + phi_shaun[nx-2])
      phiNew_shaun[1] = phiOld_shaun[1] + 2*K*dt/(2*dx)**2*(phi_shaun[3] -2*phi_shaun[1] + phi_shaun[nx-1])
      phiNew_shaun[nx] = phiNew_shaun[0]
      # update phi for the next time-step
      phiOld_shaun = phi_shaun.copy()
      phi_shaun = phiNew_shaun.copy()

#   ### my solution
#
#   # Initial conditions (in meters or m/s)
#   phi_my = initialBell(x)
#   phiNew_my = phi_my.copy()
#   phiOld_my = phi_my.copy()
#
#   # FTCS for the first (two) time-step(s), looping over space
#   for j in range(2,nx-1):
#      phi_my[j] = phiOld_my[j] + K*dt/(2*dx)**2*(phiOld_my[j+2] -phiOld_my[j] -phiOld_my[j-1] + phiOld_my[j-2])
#   # apply periodic boundary conditions
#   phi_my[0] = phiOld_my[0] + K*dt/(2*dx)**2*(phiOld_my[2] - phiOld_my[0] - phiOld_my[nx-1] + phiOld_my[nx-2])
#   phi_my[nx] = phi_my[0]
#
#   # Loop over remaining time-steps (nt) using CTCS
#   for n in range(1,nt):
#      # loop over space
#      for j in range(2,nx-1):
#         phiNew_my[j] = phiOld_my[j] + 2*K*dt/(2*dx)**2*(phi_my[j+2] -phi_my[j] - phi_my[j-1] + phi_my[j-2])
#      # apply periodic boundary conditions
#      phiNew_my[0] = phiOld_my[0] + 2*K*dt/(2*dx)**2*(phi_my[2] -phi_my[0] -phi_my[nx-1] + phi_my[nx-2])
#      phiNew_my[nx] = phiNew_my[0]
#      # update phi for the next time-step
#      phiOld_my = phi_my.copy()
#      phi_my = phiNew_my.copy()

   # derived quantities
   t = nt*dt

   # Plot the solution in comparison to the analytic solution
   plt.plot(x,initialBell(x),'k',label='initial')
   plt.plot(x,phi,'b',label='established solution')
   plt.plot(x,phi_shaun,'g',label='Shaun\'s solution')
#   plt.plot(x,phi_my,c='purple',label='My solution')
   plt.legend(loc='best')
   plt.ylabel('$\phi$')
   #plt.axhline(0,linestyle=':',color='black')
   plt.show()

#Execute the code
main()
