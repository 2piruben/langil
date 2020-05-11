# ifndef SPECIES_H
# define SPECIES_H

#include <vector>
#include <deque>
#include <string>
#include <memory>

using namespace std;

#endif

//////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////	
class species{//Class containing the info of each species.
////////////
//////////
		// I know that n and name should be private for safety,
		// but is the only way of passing "n" in a short readable way
		// actually I could use friendship to n

		public:
	
		string name; // name of the species
		double n; // amount of species (it is chosen to be double to allow macroscopic and CLE integration)

		bool ftracktime = false; // flag to track time of species that start with a lifetime, lifetime
		deque<double> determtimeque; // queue of deterministic times for deterministc time reaction species (e.g. constant production time for mRNA) 
		double lifetime; // lifetime of the species


		species(string aname,float anum){
			name=aname;
			n=anum;
		}

		species(const species &obj) { // copy constructor
			name = obj.name;
			n = obj.n;
			ftracktime = obj.ftracktime;
			determtimeque = obj.determtimeque;
		}

		species* clone(){ //function that returns a pointer to a new copied instance
			return new species(*this); // alocates memory and return pointer
		// is it necessary to repeat declaration in child class to avoid splicing??? Please check
		}

		void MakeSpeciesTimeTracking(double lt){
			ftracktime = true;
			lifetime = lt;
		}

		void SetNum(float anum){
			n=anum;
		}

		double GetNum(){
			return n;
		}

		string GetName(){
			return name;
		}

		void React(double ast){
			n+=ast;
			if (ftracktime){ //  if we keep track of time, the times are updated
				if (ast > 0){ // reaction thta adds elements to the queue
					for (int i=0; i<ast; i++){
						determtimeque.push_back(lifetime);
					}
				}
				else if (ast < 0){ // reaction that removes elements form the queue
					for (int i=0; i>ast; i--){
						if (!(determtimeque.empty())){ // only remove if possible
							determtimeque.pop_front();
						}
					}			
				}
			}
		}

		double GetNextTime(){
			if (ftracktime){
				if (n>0){
					return determtimeque.front();
				}
				else{
					return -1; // Negative -1 time means that there is not elements in the queue 
				}
			}
			else{
				return -2; // Negative -2 time means that the reaction is not supposed to be asked for times
			}

		}

		void Updatetime(double tau){
			if (ftracktime){
				for(auto &time: determtimeque){
					time -= tau;
				}	
			}
			
		}
};


typedef std::unique_ptr<species> p_species;
typedef std::vector<std::unique_ptr<species> > v_species; // shortname for vectorsv_species
