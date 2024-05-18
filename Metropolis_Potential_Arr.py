import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from drawnow import drawnow
import time

nmax = 20
N = 20
l = 10
L = np.array([2,4,6,8,10,12])
Ta = [10, 100, 500, 1000, 1500, 2000]
Tf = [10,60,100,250,400,600,900,1400,1800,2100,2500]
T = np.arange(10,1500,50)
L = np.arange(1,20,1)



def R(x, y, f, a):
    """
    Calculate the R fit parameter.

    Parameters:
        x (array-like): Independent variable data.
        y (array-like): Dependent variable data.
        f (function): Model function to fit.
        a (array-like): Parameters of the fited model function.

    Returns:
        float: The R fit parameter.
    """
    r = sum((y - f(x, *a))**2) / (len(y) - 2)
    v = np.std(y)**2
    R = 1 - (r / v)
    return R

def Pot_energy(T, N, l, nmax, plot, live):
    """
    Perform a simulation to compute the potential energy of a system of infinite potential wells.

    Parameters:
        T (float): Temperature of the system.
        N (int): Number of potential wells in the system.
        l (float): Length of the wells.
        nmax (int): Maximum level of occupancy for a well.
        plot (bool): Whether to generate a plot of the simulation results.
        live (bool): Whether to generate an animation of the distribution of energy states.

    Returns:
        tuple: A tuple containing the average energy and the heat capacit.
    """
    # Constants
    kb = 8.617e-5  # Boltzmann constant in eV*K
    h = 4.11e-15   # Planck constant in eV*s
    m = 0.510      # Electron mass in MeV*c**2
    c = 3e8        # Speed of light in m/s
    k = h**2/8/l**2/m*c**2/1e6  #Energy of the ground state

    # Variables for energy and convergence
    avg_energy = 0.
    E2 = 0.
    tc = 0
    nc = 0

    # Initialize the system with random initial states
    n = [np.random.randint(1, nmax) for _ in range(N)]
    n = np.array(n)

    # Energy of the initial state
    E = [sum(n**2)*k]

    # Convergence condition
    R = True

    # Epoch initialization
    t = 1

    # Array of epochs
    y = [t]

    # Array with information about occupancy levels
    nt = []

    # Array to store cumulative energy
    EC = [sum(n**2)*k]

    # Main simulation loop
    for j in range(5000):

        Ee = 0  # Initialize energy in epoch t

        if live and t<100: # The t limit can be change to observe the corresponding distribution at equilibrium, for a higher T this t must be greater
          # Generate live animation if enabled
            plt.clf()
            plt.hist(n, bins=15, range=(1,15), label='t = %s'%(t))
            plt.title('Histogram of Energy Levels Occupation at T = %s' %(T))
            plt.xlabel('Energy Level')
            plt.ylabel('Frequency')
            plt.ylim(0,50)
            plt.legend()
            plt.pause(0.5)
            #plt.savefig(f'histogram_animation_frame_{t}.png')

        # Perform random walk for each potential well
        for i in range(N):
            r = np.random.uniform()

            # Randomly choose to decrease or increase occupancy
            if r > 0.5 and n[i] > 1:
                nc = n[i] - 1
            else:
                nc = n[i] + 1

            # Calculate energy changes
            Ec = k*nc**2
            Ei = k*n[i]**2
            dE = Ec - Ei

            # Accept or reject the transition based on temperature
            if dE <= 0:
                n[i] = nc
            else:
                r = np.random.uniform()
                W = np.exp(-dE/T/kb)
                if r <= W:
                    n[i] = nc

            # Update energy in epoch t
            Ei = k*n[i]**2
            Ee += Ei

        # Calculate cumulative energy
        Ecum = sum(E)/t
        EC.append(Ecum)

        # Check for convergence
        if t > 1000:
            nt.append(n[i])
            Emean1 = sum(E[-1000:-500])/500   #Calculate the mean of the second last 500 iterations
            Emean2 = sum(E[-500:])/500        #Calculate the mean of the last 500 iterations
            e2 = np.array(E[-600:])           #Takes the last 600 values of energy to later calculate its cuadratic value
            R = abs(Emean2-Emean1) > 1e-3     #Sees if the means are close each other as a condition of convergence

        # Update variables if converged
        if not R:
            avg_energy = (Emean1+Emean2+avg_energy)/3   #Calculate the average energy
            E2 = (sum(e2**2)+E2)/(len(e2)+1)            #Calculate the cuadratic average energy
            if tc == 0:
                tc = t
                
        # Store energy for epoch t
        E.append(Ee)

        t += 1  # Increase epoch
        
        # Calculate cumulative energy
        Ecum = sum(E)/t
        EC.append(Ecum)
        
        # Append epoch to array
        y.append(t)
        
        # Append epoch to array
        y.append(t)

    # Convergence info
    convergence_info = f'Convergence Epoch: {tc}, Average Energy: {avg_energy:.4}'
    cv = (E2-avg_energy**2)/kb**2/T**2

    # Plot results if required
    if plot:
        fig = plt.figure(figsize=(10,4))
        fig.suptitle('Evolución de la energía de un sistema de 20 pozos infinitos de potencial \ny la distribución de los niveles de energía durante esta')

        # Plot energy evolution
        plt.subplot(1, 2, 1)
        plt.plot(y, E, label=f'Temperature: {T}')
        plt.xlabel('Epoch')
        plt.ylabel('Total energy (eV)')
        plt.title('Total energy evolution')
        plt.legend(loc='lower left')
        plt.text(0.07, 0.9, convergence_info, transform=plt.gca().transAxes, fontsize=10, verticalalignment='bottom', bbox=dict(facecolor='white', alpha=0.5))

        # Plot histogram of nt values
        plt.subplot(1, 2, 2)
        nt = np.array(nt)
        nt = nt - 1 #Redifine the ground level from 1 to 0 in order to make the adjustment
        counts, bins, _ = plt.hist(nt, bins=20, range=(0,20), density=True)
        bins = bins[:-1]

        plt.title('Histogram of Energy Levels')
        plt.xlabel('Energy Level')
        plt.ylabel('Frequency')
        plt.tight_layout()

        def boltz(E,a,b):
          return np.exp(-E*a/kb/T)*E*k/b

        #Make a curve fitting of the histogram to a Maxwell-Boltzman distribution
        resul, err = curve_fit(boltz, bins, counts, p0 = (k, 5 ) )
        b = np.arange(0,20,0.1)
        r=sum((counts-boltz(bins,*resul))**2)/(len(counts)-2)
        v=np.std(counts)**2
        R=1-(r/v)
        plt.plot(b, boltz(b, *resul ))
        plt.text(0.7,0.8, "R^2 = %.2f" %(R),  transform=plt.gca().transAxes, fontsize=10, verticalalignment='bottom', bbox=dict(facecolor='white', alpha=0.5))


        plt.show()

    return avg_energy, cv

def Potential_reservoir(Ta, N, l, nmax, plot=False, live=False):
    """
    Perform a simulation to compute the total energy of an array of N potential
    wells when the system reach the equilibrium with a thermal reservoir at a temperature T.

    Parameters:
        Ta (float, list, np.ndarray): Temperature or array of temperatures.
        N (int): Number of potential wells in the system.
        l (float, list, np.ndarray): Length of the wells or array of lengths in nm.
        nmax (int): Maximum level of occupancy for a well.
        plot (bool): Whether to generate a plot of the simulation results.
        live (bool): Whether to generate an animation of the distribution of energy states.

    Returns:
        tuple: A tuple containing arrays of average energy and heat capacity.
    """

    # Convert well width to meters
    l = l * 1e-9

    # Lists to store average energy and heat capacity
    E_avg = []
    Cv = []

    # Check if Ta is an array of temperatures
    if isinstance(Ta, (list, np.ndarray)):
        for i in Ta:
            E, C = Pot_energy(i, N, l, nmax, plot, live)
            E_avg.append(E)
            Cv.append(C)
        return E_avg

    # Check if l is an array of lengths
    elif isinstance(l, (list, np.ndarray)):
        for i in l:
            E, C = Pot_energy(Ta, N, i, nmax, plot, live)
            E_avg.append(E)
            Cv.append(C)
        return E_avg

    # Otherwise, perform simulation with single temperature and length
    else:
        E_avg, Cv = Pot_energy(Ta, N, l, nmax, plot, live)
        return E_avg

E_avg = Potential_reservoir(Tf, N, l, nmax, plot = True )

E_avgl = Potential_reservoir(400, N, L, nmax, plot = False )


#Definition of the functions to fit
def line(x,a,b):
  return a*x + b

def e(x,a,b,c):
  return np.exp(-x*a)*b+c

#Make a curve fit of the results

result, errt = curve_fit(line,Tf,E_avg)
Rt = R(Tf, E_avg, line, result)
print(result)
resull, errl = curve_fit(e,L,E_avgl)
Rl = R(L, E_avgl, e, resull)
print(resull)


fig = plt.figure(figsize=(10,4))
fig.suptitle('Average energy respect to the temperature and the well size')

plt.subplot(1, 2, 1)
plt.plot(Tf, E_avg, label='l= %s nm'%(l))
plt.plot(Tf, line(Tf,*result ), label='Fit')
plt.title('Total energy vs Temperature')
plt.text(0.7,0.6, "R^2 = %.3f" %(Rt),  transform=plt.gca().transAxes, fontsize=10, verticalalignment='bottom', bbox=dict(facecolor='white', alpha=0.5))

plt.ylabel('Total average energy (eV)')
plt.xlabel('Temperature')
plt.legend()

plt.subplot(1, 2, 2)
plt.plot(L, E_avgl, label='T = 400')
plt.plot(L,e(L,*resull), label='Fit')
plt.title('Total energy vs Well size')
plt.text(0.7,0.6, "R^2 = %.3f" %(Rl),  transform=plt.gca().transAxes, fontsize=10, verticalalignment='bottom', bbox=dict(facecolor='white', alpha=0.5))

plt.ylabel('Total average energy (eV)')
plt.xlabel('Potential well size (nm)')
plt.legend()


plt.show()



#Plot of the calculated heat capacities. It's hiden beacuse the result is not convincing
'''fig = plt.figure(figsize=(10,4))
fig.suptitle('Average energy respect to the potential well size at T = 400 ')'''

'''plt.subplot(1, 2, 2)
plt.plot(Tf, Cv, label='l= %s nm'%(l))
plt.xlabel('Temperature')
plt.ylabel('Cv')

plt.legend()
plt.tight_layout()'''

#Plot of the calculated heat capacity. It's hiden beacuse the result is not convincing
'''plt.subplot(1, 2, 2)
plt.plot(L, Cvl, label='T = 400')
plt.xlabel('Potential well size (nm)')
plt.ylabel('Cv')

plt.legend()
plt.tight_layout()
plt.show()'''
