# SISSA Summer School

This repo includes all the codes I have done during the summer school at SISSA (International School for Advanced Study), Trieste, Italy in July of 2015. The summer school lasted for four weeks. 

## Preschool
The first week is about training for programmming. Although most of them had already learnt how to programming before attend the summer school, the preschool is still very interesting in that we were obeseed into the intense courses studied new tools. 

### Day 1&2
The first two days for the summer school was about basis of Fortran, the old yet superior language often used in science community. We learned basically all the important components and synatax in Fortran 90, including subroutine and dynamic allocation, which might be unique in Fortran. 

### Day 3 
The third day during the preschool was about AWK, a scripting language that is very useful when analysing data. The cool thing about AWK is that it reads the data line by line automatically, and it use '$n' to spcify the nth column. Therefore if we have a datasheet with several columns, it will be a good idea to use AWK to analyse them.

### Day 4
The last day of preschool was meant to study numerical intergration and numerical differention, and try to apply them in solving harmonic Schordinger Equation. This day was pretty hard since that some problem that might be trivial analytically turn out to be problematic when treated numerically. 


## Week 1 - Monte Carlo
After the preshool, the first week is used to study Monte Carlo, expecially Metropolis algorithm. 

### Day 1 
The first day started from basic probability theory, such as mean and variance of random variable. They it was disscussed that how to generating random numbers in a computer (a process often called 'sampling'). It should be bared in mind that uniform distribution is the ONLY one that can be sampled by the computer automatically. All the other distributions require some sort of efforts. Then we studied how to sample non-uniform distribution, with two methods, rejection method and
non-rejection method, both for discrete distribution and continuous ones. 

Some efforts were made to prove Central Limit Theorem (CLM), which gives the distribution of the average of N random variables. And as N grows very large, it essentially follows a Gaussian distribution. 

The last but not least extension was about importance sampling, which is often used to deal with distribution that has narrow and high peaks. The idea was that sample random variable that has higher possibility to fall into the peaks instead of uniform distributed. 

### Day 2

Tuesday was even more mathematical and physical, since we were exhaused to study Markov Chains rigorously. 

The simple example about using Markov Chain is finding the integral, say the size of a circle. Starting from a random walker and iterate pace by pace. Every time the walker fall into the circle the step is accepted, otherwise the step is rejected. Therefore the location of after a step depends and only depends on the local of last step. 

This idea can be put in more non-trivial and intriguing systems to find the equilibrium state of a many-body system, such as the Boltzman distribution of Lennard-Jones atoms. We proved Perron-Frobenus theorem that starting from every (or almost every) initial configuration we will end up in Boltzman distribution. 

Then we found the solution of transition probability that gives rise to such an evolution, which is given by the lovely Metropolis algorithem. 

### Day 3

Comparing to the compact and mathematical way we used on Tuesday, on wednesday we spent a large amout of time reviewing Markov chain in a more pedological and graphical way. Then we generalize to metropolis algorithem in to two cases: annealing algorithem and the so called "talks between different temperatures", in a way that to make the computation more efficient. 

Last we studied a physics idea that fluctuation at a given temperature actually tells information of the response of the system with respect to changing temperature. For example, specific heat is given by the fluctuation of internal energy, susceptibility is given by fluctuation of magnetic field. And from the fluctuation depending on temperature, we can actually found the density of state which DOES NOT depend on temperature. Therfore we can just use discrete temperature simulations
to find continuous information of the system. 

### Day 4 

On Thursday we focused on Stochastic Process, expecially Langvein Equation, where a random force and a force depending on the potential compete with each other. Use the idea of Markov Chain, we found the master equation of the probability of taking some position (or configuration) at discrete time steps (Fokker-Planck Equation). And then we proved that FPE satisfies the equilibrium state. 

Now a question is that can we reach the equilibrium state from any initialized configuration or position. The short answer is yes. If we compart with Shordinger equation, we will draw the conclusion that the system will evolve to graound state (GS). And since equilibrium distribution probability has NO NODES, it is essentially GS. 

Finally we move MC a little forward to study variational Monte Carlo, which is used to treate Quantum System where the wave function is not known in advance so we cannont sample the system accoring to exact distribuion. The idea is that we start from some trial wave function with some parameters, and use metropolis theorem to reach to "pseudo-equilibrium" state and fluctuate around. Then we change the parameters with the trial wfc to reach another "pseudo-equilibrium" state.
Eventually we will find the lowest possible energy as we adjust the parameters. And that should be our approximation for the GS energy. 














