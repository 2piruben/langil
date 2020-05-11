# Stochastic integration of trajectories of molecular content

'langil' is a C++ code that integrates stochastic trajectories of different interacting molecular species using the Gillespie algorithm or a continuous stochastic differential equation description.

## Description of files

* [langil.cpp](integration_bs.py) Definition of the `langil` class containing the main integration routines. Particular description of its usage can be found for the [HIV example](examples/HIV)
* [langil.h](langil.f) Headers for the `langil` class and definition of the `reaction` class defining the different types of reactions that can be used by `langil`
* [species.h](species.h) Definition of the `species` class that contains the different types of the moleculas species used by `langil`
* [cell_cycle.h](cell_cycle.h) Definition of the `cell_cycle` class used to incorporate cell cycle properties (cell cycle phases of different stochastic duration including material division at bipartition). Examples using it can be found in the folder `examples/CellCycleVariability`
* [examples](examples) Folder containing different projects using the `langil` class

## Compilation
There is a Makefile provided to compile the code, to compile the `langil` class run
```
make
```
To compile the code and the examples provided run
```
make examples
```
To clean all the object files run
```
make clean
```

## Examples

* [HIV](examples/HIV) Example for a simple model of HIV cell population infection. It contains a detailed description of the use of the class `langil` in the [HIV example](examples/HIV/HIV.html)
* [Cell Cycle Variability](examples/CellCycleVariability) Code used in the paper *R. Perez-Carrasco, C. Beentjes, R. Grima, Effects of cell cycle variability on lineage and population measurements of mRNA abundance. BioRxiv (2019)* to simulate mRNA abundance in cells with a stochastic cell cycle duration


### Prerequisites

`langil` uses the random number generation routines of the [GNU Scientific Library](https://www.gnu.org/software/gsl/).  

## Author

* **Ruben Perez-Carrasco** - [2piruben](https://github.com/2piruben)

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

    