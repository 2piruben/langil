# ifndef CELL_CYCLE_H
# define CELL_CYCLE_H

#include <vector>
#include <deque>
#include <string>
#include <memory>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <math.h>


using namespace std;

#define CONSTANT_PHASE 0
#define EXPONENTIAL_PHASE 1

#endif

//////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////	
class cell_cycle{//Class containing behaviour related to cell_cycle.
////////////
//////////

		private:
	
		int cell_state; // current state of the cycle
		double time_state; // current time in the state of the cycle
		double time_cycle; // total time in this cell_cycle
		double currentPhaseDuration; // time for the next phase
		vector<double> state_durations; // array with the cell states durations
		vector<int> state_types; // array with the cell phases types (exponential, uniform, ... )
		int divisions; // number of divisions so far
		int num_states; // total number of cell phases included
		bool divide_material; // flag to see if the species in the cell are dividied after cell division
		gsl_rng * rng; // allocator for the rng generator


		public:

		cell_cycle(){ // default constructor is a constructor with empty info
			time_state = 0;
			time_cycle = 0;
			divisions = 0;
			cell_state = 0;
			num_states = 0;
			rng= gsl_rng_alloc (gsl_rng_mt19937);
			gsl_rng_set (rng,::time(NULL)*getpid());
		}	

		void addCellPhase(double duration, int type_phase){
			cout<<"Adding Cell Phase of duration"<< duration<<" \n";
			state_durations.push_back(duration);
			state_types.push_back(type_phase);
			num_states++;
			if (num_states == 1){// if this is the first state set, also use it for the current cell cycle time
				getPhaseDuration();
			}
			cout<<"Cell Phase added\n";
		}

		void ResetCellCycle(){
			time_state = 0;
			time_cycle = 0;
			cell_state = 0;
			getPhaseDuration();
			divisions = 0;
			divide_material = false;
		}

		// cell_cycle(int nstates,float duration, int type){
		// 	// The default constructor generates n cycle phases with same duration of certain type
		// 	time_state = 0;
		// 	time_cycle = 0;
		// 	divisions = 0;
		// 	cell_state = 0;
		// 	num_states = nstates;
		// 	state_durations = new double[num_states];
		// 	state_types = new int[num_states];
		// 	for(int i=0;i<nstates;i++){
		// 		state_durations[i] = duration;
		// 		state_types[i] = type;
		// 	} 
		// }

		void setPhaseDuration(int phase, double duration, int type_phase){
			state_durations[phase] = duration;
			state_types[phase] = type_phase;
		}

		double getPhaseDuration(){
		// compute the value of the time for the next phase and store it in NextPhaseTime to avoid recomputing
			if (state_types[cell_state] == CONSTANT_PHASE){
				currentPhaseDuration = state_durations[cell_state]; 
			}
			else if (state_types[cell_state] == EXPONENTIAL_PHASE){
				currentPhaseDuration = state_durations[cell_state] * log(1.0/gsl_rng_uniform_pos(rng)); 
			}

			// cout<<"Next Phase duration: "<<currentPhaseDuration<<'\n';
			return currentPhaseDuration;
		}

		double getTimeForNextCellPhase(){
			return currentPhaseDuration - time_state;
		}

		int AdvanceNextPhase(){
			cell_state++;
			if (cell_state==num_states){// if the cell reaches the final state
				cell_state = 0;
				time_cycle = 0;
			}
			time_state = 0;
			getPhaseDuration();
			if (cell_state == 0){ // if we are dividing
				divisions += 1;
				divide_material = true;
			}
			else{
				divide_material = false;
			}
			return cell_state;
		}

		bool MaterialIsDivided(){
			return divide_material; 
		}

		void updateCelltime(double time){
			time_cycle += time;
			time_state += time;
		}

		void setCellCyclePositionTime(int phase, double time){
			cell_state = phase;
			time_state = time;
			time_cycle = time;
		}

		void printCellProperties(){
			cout<<"total states "<<num_states<<'\n';
			cout<<"of duration ";
			for (auto d: state_durations){ cout << d << ' ';};	
			cout<<"\n of type";
			for (auto d: state_types){ cout << d << ' ';};		
			cout<<'\n';
		}
	
		int getCellPhase(){
			return cell_state;
		}

		double getPhaseTime(){
			return time_state;
		}

		void PrintState(){
			cout<< "cell state: "<<cell_state<< " time_state: "<<time_state<<'\n';
			cout<< "currentPhaseDuration: "<<currentPhaseDuration<<'\n';
			printCellProperties();
			cout<< "time cell has divided: "<<divisions<<'\n';
			cout<< "cell is ready for division: "<< divide_material<<'\n';
		}

};