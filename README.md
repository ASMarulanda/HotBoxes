# Monte Carlo Simulation of Particle Distribution in Potential Wells with Thermal Reservoir

In this work, we study the scenario of an array of infinite potential wells in contact with a thermal reservoir, each containing one electron. It is well known that in certain physical systems, the likelihood of being in a particular state is determined by the Boltzmann distribution, which is a function of energy and temperature.

$P(n_i) = \frac{e^{-E_i/k_b T}}{Z}$

This project aims to explore the distribution of particles within the energy levels of the system upon reaching equilibrium. Additionally, it aims to unveil how the system's behaviour varies concerning different parameters, such as the length of the wells. To achieve this, we will use the Monte Carlo approach with the Metropolis algorithm in which we make the system evolve according to certain probabilities. We aim to show how temperature affects the final distribution of particles on energy levels, the stability of the system, and the convergence of simulation results.

![Potenciales](https://github.com/samuelquitiang/HotBoxes/assets/53834570/66a53846-9845-46b1-90f9-aaeea9430fb4)

In this case, the initial state of the system can be expressed as:
$n_0 = (4, 3, 2, ..., 4)$

The code will evolve the system from a state $n_i$ to a state $n_j$ by doing a "random walk" in which the probability of doing a transition to a higher level of energy is given by:

$W=P\left(n_i \rightarrow n_j\right)=\frac{p\left(E_{n_j}\right)}{p\left(E_{n_i}\right)}=\frac{e^{-\beta E_{n_j}} / Z}{e^{-\beta E_{n_i} / Z}}=e^{\beta\left(E_{n_j}-E_{n_i}\right)}$

After taking into account the probability of moving to a lower energy state due to the system's tendency to reach its ground state. The transition conditions are computed by generating a random number and the following test:

*    $r: random.uniform()<0.5$: Downward transition
*    if $r>0.5$ $and$ $r': random.uniform() < W$ : Upward transition

The code has the function:

Potential_reservoir(Ta, N, l, nmax, plot=False, live=False)

    Perform a simulation to compute the total energy of an array of N potential wells when the system reach the equilibrium with a thermal reservoir at a temperature T. The functions make use of the function Pot_energy(T, N, l, nmax, plot, live).

    Parameters:
        Ta (float, list, np.ndarray): Temperature or array of temperatures.
        N (int): Number of potential wells in the system.
        l (float, list, np.ndarray): Length of the wells or array of lengths in nm.
        nmax (int): Maximum level of occupancy for a well.
        plot (bool): Whether to generate a plot of the simulation results. Total energy vs epoch and a histogram of the occupation distribution of the energy states. Default: False
        live (bool): Whether to generate an animation of the distribution of energy states. Default: False

    Returns:
        Tuple: Array or integer representing the average energy. The format of the tuple varies depending on the type of the Ta parameter: if Ta is an integer, both average energy and heat capacity are single values; if Ta is a list, each element of the tuple corresponds to an array of values, and similarly for the l parameter.
    

This is the main function for the user, from it the user can obtain how is the behaviour of the system energy respect to the different parameters.


# Heat capacity problem
The heat capacity of a system can be computed from its energy behaviour when it is in thermal equilibrium as:

$C_v=\frac{\left\langle E^2\right\rangle - \langle E\rangle^2}{k_b T^2}$.

Where $k_b$ is the Boltzmann constant. 

In this project, we tried to get the heat capacity from the evolution of the system but the results were not convincing and also varied a lot in each run due to the Markovian-stochasticity of the problem, this could be due to a bad calculation of the heat capacity or maybe by taking into account more potential wells this heat capacity stabilizes.

# Results - Graphics

![T100](https://github.com/samuelquitiang/HotBoxes/assets/53834570/f226dc94-0798-41c8-9edc-067b5c17841d)

This figure shows the behaviour of the system when the reservoir temperature is at $100 K$. At the beginning, the system has high energy, but after several iterations, the energy starts to decay, and eventually, it oscillates around a medium energy value. This oscillation shows minimal variation around the average value for this temperature. 

![T1800](https://github.com/samuelquitiang/HotBoxes/assets/53834570/350c0258-28c3-4e14-84fe-1ceddcd671d8)

The figure shows the system's behaviour when it is in contact with a thermal bath at a temperature of $2500K$. Similar to the previous case, the system starts with high energy, decreases gradually, and oscillates around a certain energy value. However, in this case, the oscillation is more significant as compared to the ones observed when the initial temperature was $100K$.

![Histogram](https://github.com/samuelquitiang/HotBoxes/assets/53834570/d1eef453-e8aa-4cf0-9dcd-fd57b03780ac)
Upon reaching the equilibrium, the system the occupation state distribution follows the Maxwell-Boltzmann distribution,
this distribution describes the probability of finding particles at various energy levels within the system, reflecting the thermal equilibrium achieved under the specified conditions. On this case, the higher energy levels are less likely to be occupied.

![TvsE](https://github.com/samuelquitiang/HotBoxes/assets/53834570/6e2b25c4-15bd-4497-af06-921ccb22b78a)
In this Figure is presented how the total energy of the system exhibits a linear increase with the temperature of the thermal reservoir, indicating a direct correlation between temperature and energy content. This observation follows the fundamental principle of thermodynamics, where higher temperatures lead to greater energy contributions within the system.

![EvsL](https://github.com/samuelquitiang/HotBoxes/assets/53834570/864fa088-643f-41a5-92d5-89099d15a628)

The relationship between the total energy and the length of the potential well, as the previous image, reveals an inverse exponential function. This trend suggests that narrower wells result in significantly higher energy compared to wider ones, implying a sensitivity of the system's energy state to variations in well dimensions.
