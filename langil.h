# ifndef LANGIL_H
# define LANGIL_H

///////////////////////////////////////////////////////
//////   GILLESPIE CLASS
//////
//////						R. Perez-Carrasco
//////                               Created: 16 Nov '13
///////////////////////////////////////////////////////
/* Class for simulating master equation through langil algorythm. It defines three classes:
		- langil: In charge of the integration
		- species: contains a name a number and methods to change number
		- reaction: contains a name, a stoiciometry associated and a pointer to a propensity computing function		
*/

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <string>
#include <fstream>
#include <iostream>
#include <unordered_map>
#include <vector>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include <iomanip>
//#include <omp.h>

#include "species.h"
#include "cell_cycle.h"

#define VERBOSE 1
#define NON_VERBOSE 0

#define TIMESTEPMULT 3
#define TIMESTEP 0
#define ALL 1
#define NOWRITE 2
#define LASTPOINT 4

#define GILLESPIE 0
#define MACROSCOPIC 2
#define CLEEULER 1
#define CLEMILSTEIN 3		

#define ADIABATIC true
#define NONADIABATIC false

#define GEOMETRIC true
#define NONGEOMETRIC false

#define DETERMINISTIC_TIME true
#define STOCHASTIC_TIME false

#define EXPONENTIAL_REACTION 0
#define GEOMETRIC_REACTION 1
#define DETERMINISTIC_REACTION 2
#define ADIABATIC_REACTION 3


using namespace std;

//   struct reaction{ //Structure containing the info of each reaction. Each reaction added 
// through the method Add_Reaction(reaction) will be stored as an element in a vector
// by langil class.
//	string name; // name of the reaction
//	double (*propensity)(); // function of the propensity of the current reaction
//	vector<int> stoichiometry; // stoichiometry of each reaction
//	double prop; // actual value of the propensity (updated through propensity function)
//	reaction(string aname,double (*aprop)(),vector<int> asto, double apropx)
//		:name(aname), propensity(aprop),stoichiometry(asto),prop(apropx) {}
//};
// Since the initialization does not require the value prop, an initializator omiting it
// should be defined (it is not urgent though, any value in the initialization will give the
// same result). 


struct StructDivisionState{
	double time;
	double species0before;
	double species1before;
	double species0after;
	double species1after;
};


class reaction; // defined after langil
class species; // defined after reaction
class expressionzone; // defined in a separate file

class langil{
	
	friend class expressionzone;

//	friend class reaction; // reaction have access to the state of the system 

	public:

		langil(string outputfilename, double vol=100, int seed=-1);
		void SetState(string aname, float anumber); //Set Initial conditions
		double GetState(string a); // get the state of the selected species
		void SetTime(double time); //Set current time
		void AddSpecies(string aname, float anum); // Add a new molecular species
		void AddSpeciesTimeTracking(string aname, float anum, double lifetime); //  Add a molecular species with a timer
		void AddCellCycleSpecies(); // add a cell cycle to the system
		void addCellPhase(double duration, int type_phase); // add a cell phase to the cell cycle
		void PrintSummarySpecies();
		string StateString();
		void Add_Reaction(string name,vector<int> stoich,double (*prop_f)(v_species&),int typeofrec = EXPONENTIAL_REACTION); // Add a reaction to the system overloaded to create directly
		void UpdateTime(); // Update times for the deterministic events 
		void SetRunType(int rt,double dt=0);
		double GillespieStep(double timelimit = -1); // Advance one step in the Gillespie algorithm. Returns time of the step. Timelimit does not react if the reaction would take more than timelimit
		double LangevinEulerStep(double timelimit = -1); // Make an Euler-Maruyama integration step
		double LangevinMilsteinStep(double timelimit = -1); // Make an SDE integration step with the Milstein algorithm (not tested)
		double MacroStep(double timlimit = -1); // Make an Euler (deterministic) step
		double MacroAdiabaticStep(double dtt); //  Make an Euler step if the Adiabaitc species		
		double Run(double time, bool verbose=true); // React until reach time time
		double RunTimes(int N,double time, bool verbose=true); // Run Run() N times
		double RunTransition(v_species& initstate, v_species& endstate,
			v_species& initradius, v_species& endradius,
			 int loops, double maxT); // Run transition between two points, if loops=1 it also records the transition trajectory
		void WriteState(); // Write state in file 
		void WriteTempState(); // Write state in file 
		void WritePrevrecState(); // Write previous recorded state in file 
		void SetWriteState(int flag, double tl=0);// Switch betweem different states:
		double (langil::*MakeStep)(double); // pointer to the actual integration chosen
			// TIMESTEP: write the state every tl time lapse
			// ALL: write after each reaction takes place
			// NOWRITE: don't write anything in output file
		double SetLangevinTimeStep(double dx); // set dt to a characteristic time that will increase with the size of the system
		void Add_TimeAction(double (*actionfunc)(double));
		void RunTimeAction();
		void setPhaseDuration(int phase, double duration, int type_phase);
		void Set_AdiabaticReaction(string aname);
		void Set_GeometricReaction(string aname);
		void Set_Boundary_Behaviour(void (*b_b)(v_species&, v_species&)); // Boundary behaviour
		void Set_StoreDivisionTimes();// activate the storage of division times in DivisionTimes vector
		void (*Boundary_Behaviour)(v_species&, v_species&); // boundary behaviour function
		vector<StructDivisionState> Get_HistoryDivision();
		void Reset_HistoryDivision();
		void SetDivisionHistory(int b,int a);


	protected:

		v_species x; //value of each species. shared_ptr is used since the vector can contain species and also childrens of that class (i.e. different classes)
		vector<reaction> r;//vector with all the reactions
 		double time; // current time of simulation
		gsl_rng * rng; // allocator for the rng generator
		vector<double> rnd; // rnd numbers
		vector<double> xtemp; // temporal vector for integrations
		v_species xtemp0,xtemp1; // temporal vector for integrations
		double pro,pro0; // auxiliar variables for the program
		int sto;  // auxiliar value for stoichiometry
		double totalprop,cumprop, detprop; // sum of the propensities
		double nexttau; // time of the next step
		vector<reaction>::iterator nextreaction; // reaction selected to react
		vector<int>* nextstoichiometry;
		double Omega; // volume of the system to change between concentrations and absolute numbers
		double signal; // variable that can modulate an external signal its time course
						// may be changed through time actions with Add_TimeAction
		int celldivided; // flag to track possible celldivision
		int geostoichiometry; // dummy variable to set random geometric stoichiometries
	
		bool record; // Set record of the trajectory
		ofstream trajfile; // File for output traj
		string outputfilename; // File for output name
		int fwrite; // flag for write state when integrating
		double timelapse; // dt for recording if TIMESTEP
		double nextrectime; // t for recording if TIMESTEP
		double totaltime; // total integration time for the trajectory
		double dt; // timestep used in integrations for macroscopic and langevin (if negative, it is used by Langevin as a prediction for dx)
		
		int cellcyclespeciesidx;  // index of the species vector that tracks cell cycle phases
 		cell_cycle cell; // cell_cycle object to manage cell cycle events

 		bool fStoreDivisionTimes; // flag to store Division Times
 		vector<StructDivisionState> DivisionHistoryVector; // Division history vector. Each component is a vector 
 		StructDivisionState divisionstate;
		vector<double(*)(double)> actionfuncvec; // vector to time action functions

};




////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
class reaction{//Class containing the info of each reaction. Each reaction added 
// through the method Add_Reaction(reaction) will be stored as an element in a vector
// by langil class. The pointer x, points out to the state of the langil system and
// is managed automatically by langil when a reaction is added to a langil instance
// It has also the property "adiabatic" meaning that it is a fast variable and can always be simulated
// from its deterministic behaviour

	friend class langil;

	public:
		
		
	reaction(string aname,vector<int> asto,double (*aprop_f)(v_species&)=NULL,int reacttype= EXPONENTIAL_REACTION, 
		    vector<reaction>* x0=NULL){
			name=aname;
			stoichiometry=asto;
			prop_f=aprop_f;
			reactiontype=reacttype;
		}
		double GetPropensity(){
			return prop_f(*x);
		}
		double GetPropensity(v_species& x0){
			return prop_f(x0);
		}
		string GetName(){
			return name;
		}
		void SetPropensity(double (*p_f)(v_species&)){
			prop_f=p_f;
		}
		void SetState(v_species& x0){
			x=&x0;
		}
		vector<int>* GetStoichiometry(){
			return &stoichiometry;
		}
		bool IsAdiabatic(){
			return (reactiontype == ADIABATIC_REACTION);
		}
		bool IsNotAdiabatic(){
			return (reactiontype != ADIABATIC_REACTION);
		}
		bool IsGeometric(){
			return (reactiontype == GEOMETRIC_REACTION);
		}
		bool IsNotGeometric(){
			return (reactiontype != GEOMETRIC_REACTION);
		}
		bool IsNotDetermTime(){
			return (reactiontype != (DETERMINISTIC_REACTION));
		}
		bool IsDetermTime(){
			return (reactiontype == (DETERMINISTIC_REACTION));
		}
		void SetAdiabatic(){
			reactiontype = ADIABATIC_REACTION;
		}
		void SetGeometric(){
			reactiontype = GEOMETRIC_REACTION;
		}
		void SetDeterministicTime(){
			reactiontype = DETERMINISTIC_REACTION;
		}	

	protected:

		string name; // name of the reactio
		vector<int> stoichiometry; // stoichiometry
		double (*prop_f)(v_species&); // propensity function
		v_species* x;  //value of each species. shared_ptr is used since the vector can contain species and also childrens of that class (i.e. different classes)
		int reactiontype; // can be any of the macros EXPONENTIAL_REACTION, DETERMINISITIC_REACTION, etc.

 // of the system


};



#endif
