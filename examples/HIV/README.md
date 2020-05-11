## Using Langil class

Langil is a C++ class that produces realisations of a Master Equation with the Gillespie algorithm or with Langevin dynamics. It uses the auxiliar class `reaction` defined in the header `langil.h`.

`HIVlangil.cpp` contains an example of usage of `langil`. The code is structured as follows:

###lines(1-65) headers, and definition of variables. 
Note that the dimension (number of species) is `DIM=4` in this case. And the species vectors will all of them consist of vectors (x1,x2,x3,x4).

###lines(78-155) creation of the stoichiometry and propensity of the reactions

For instance, the first reaction (Production of T) reads

```c++
// Production of T 
vector<int> Tpro_r(){ return {1,0,0,0};};
double Tpro(vector<species>& x){
return Omega*lambdaT;
};
```

The stoichiometry contains the values of change of each species when the reaction channel fires, in the example case, only the first of the 4 species is increased by 1. Note that there is already a set order in the species.

The second definition is the propensity of that reaction, that returns the probability per unit of time of the reaction to occur.  

###lines(157-168) command line call of HIVlangil

This deines how HIVlangil will be called from the termina

`HIVlangil file.in  file.out`

This is useful for the parallel calls from `pyHIVlangil.py`, which is a python script that creates different input files and runs HIVlangil with them in parallel. 

###lines(169-235) read of the input file

The inpunt file is created from `pyHIVlangil.py`. The file `CLEN0.in` is one example of input file, it can be modified by hand, but is highly recommended to automate this with `pyHIVlangil.py`.

###line(238) creation of instance `HIV` of class `langil` 

```c++
gillespie HIV(outputfilename,Omega,SEED); 
```

Note that the constructor only requires the name of the outputfile, the value of the system volume `Omega` (relating intensive with extensive magnitudes) and the `SEED` of the RNG (SEED=-1 uses a random seed).

### lines(249-252) adding species to `HIV`

Each species is added separately stating a string with its name, and its initial value as,

```c++
HIV.AddSpecies("T",T0*Omega); // Careful with the order
```

Note that the species have to be added in the same order that the stated in the stoichiometry vector. This is something that has to be changed in the future. Also note that the function requires absolute numbers of molecules. So if your input is concentration of molecules you have to rescale them with the parameters `Omega` (as in the example).

### lines (258-269) adding the reactions to `HIV`

Each reaction requires a string with the name of the reaction, the vector with the stoichiometries and the propensity function.

```c++
HIV.Add_Reaction("T Production",Tpro_r(),Tpro);
```

### lines (271-285) Declaring adiabadicity in reactions

Some reactions have a much faster dynamic than others. Such reactions can be considered being always in equilibrium with the rest. In Gillespie simulations they can be a problem since they fire much more often than others, without any change in the dynamics. For this reason, setting a reaction as adiabatic will remove the reaction from the propensity calculus and compute the changes that its mass action law approximation produces after a Gillespie step. A reaction can be set adiabatic after the reactions are declared (as happens in the code)

```c++
HIV.Set_AdiabaticReaction("T Production");
```

but could also be stated directly from the reaction declaration as

```c++
HIV.Add_Reaction("T Production",Tpro_r(),Tpro,ADIABATIC);
```

### line (288) Set RunType

RunType is the integration scheme chosen. So far, there are only 3 working options which are:

- *GILLESPIE=0* : Gillespie algorithm
- *CLE=1* : Chemical Langevin Equation
- *MACROSCOPIC=2* : Deterministic Mass Action law

(The code has some other Milstein integrators for the CLE, but they have not been fully tested). The RunType can be indicated by the interger or by their macro name as:

```c++
HIV.SetRunType(GILLESPIE,dt); 
```

The argument `dt` has no effect on Gillespie integration. For Macroscopic simulations is the step time of the integration of the mass action law. Finally, in the CLE case, `dt` if `dt>0` it is used as the time step for the CLE, while `dt<0`, a time step time is computed so the typical advance distance in a integration step do not involve a change in the concentration of the species greater than |`dt`| (use this last option with caution). 

### line (289) Set WriteState

WriteState states how the output data from the simulation should be stored,

```c++
HIV.SetWriteState(ffwrite,timelapse);
```

It takes two arguments. The options for `ffwrite` are:

 - `TIMESTEP = 0`, write the state of the system every `timelapse` amount of time
 - `ALL = 1`, write after each reaction (integration step) takes place
 - `NOWRITE = 2` no output
 - `TIMESTEPMULT`, stores steps following a geometrical time series, is very useful when the simulation can get very long and a set of events are expected to follow a exponential time distribution, such as scape times. The formula for the record times is: $T(i)=dt*(timelapse^i-1)$


### line (291) Add Temporal actions

Temporal actions are events that can change the behaviour of the system when some condition is met. This functions will be called after each integration. In the current case, there is a time action defined together with the reactions in lines (151-155)

```c++
// Change of aL in time
double aLtime(double t){
aL=aLs+aLC*(1-aLC2*cos(2*M_PI*t/aLT+M_PI));
return aL;
}
```

which defined the change of the parameter aL in time and is introduced in `HIV` in line (291) as,

```c++
HIV.Add_TimeAction(&aLtime);
```

### line (293) Run the integration

There are 3 ways to Run the integration, depending on the desired output

```c++
HIV.Run(totaltime)
```

Makes one single trajectory with a time length `totaltime`

```c++
HIV.RunTimes(N,totaltime)
```

Makes `N` trajectories and stores the last point of each trajectory. Useful to study decision over multiple stable states.

```c++
HIV.RunTransition(initialstate, finalstate, 
                  initialstatedist, finalstatedist, 
                  N,totaltime);
```

Thist last option is useful when you want the trajectory to end when it reaches a certain region in the phase space. It untegrates `N` trajectories starting at `initialstate` that stop when `finalstate` is reached or when time reaches `totaltime`. The output is the time required to reach `finalstate` (Mean First Passage Time) for each trajectory and a flag (1 or 0) stating is the integration terminated by reaching `finalstate` (0) or `totaltime` (1).

`initialstate`, `finalstate`, `initialstatedist`, `finalstatedist`, are vectors in the species space. The trajectory finishes when it enters the parallelepiped defined by `finalstate`+`finalstatedist` and `finalstate`-`finalstatedist`.

If the number of trajectories computed is `N`=1, then the whole trajectory is stored. The trajectory stored will span from the last time the trajectory left the initial parallelepided before reaching the final parallelepided. 

