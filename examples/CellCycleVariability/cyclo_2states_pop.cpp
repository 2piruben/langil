
/////////////////////////////////////////////////////////////////////////
//
// GILLESPIE SIMULATION CLASS OF CYCLOSTATIONARY PARTITION
//
//							R Perez-Carrasco
//							Created: 4 May '19
//////////////////////////////////////////////////////////////////////
//
// The system contains 8 species that states the possible 4 states of each of the two genes Olig2 and Irx3
// Each of the species are the following:
//
// promO: promoter for Olig2 gene, can be active or inactive i.e (can take 0 or 1). If we consider diploidity
//        it can take larger values
// nascO: number of nascent Olig2 genes, they are treated separatedly to mature mRNA since they cannot be transcribed, 
//        can be related to FISH measurements, and can have deterministic lifetimes 
// mrnaO: mature Olig2 mRNAs
// protO: number of Olig2 proteins
//
// similarly we can define promI, nacI, mrnaI, protI
//
// They will be ordered as (0: promO, 1:nascO, 2:mrnaO, 3:protO ,4: promI, 5:nascI, 6:mrnaI, 7:protI

#include "langil.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

#define DIM 3 //// Dimension, number of reactants (species)

using namespace std;

	const int dim=DIM;
	// Main control parameters
	double Omega;// Volume of the system, relates concentrations of mRNA and protein to actual numbers

	double k0,k1; // rates of pproduction/death mRNA
	double k2,k3; // rates of production/death protein
	double T; // cell cycle time
	double r1,r2; // rates to transit between


	int stochasticcycle; // are the cell pahses stochastic?
	int phasenum; // number of phases per cell cycle
	int presynnum;  // number of cell cycle phases before gene doubling
	vector<double> phaserates; // rates of the different cell cycle phases

	double m0,p0; // Initial conditions

	int mRNAgeo; // flag to determine if mRNA is produced geometrically (-1 is constitutive, >0 is the burst average size)


	// Model solving parameters

	int ffwrite; // kind of recording protocol (defined in .in or .h)
	double timelapse; // timestep for recording
	double fathertime; // 
	double totaltime;
	double runtimes; 
	double dt;
	int runtype;
	int SEED; // RNG SEED

	////////////////////
	// DEFINITION OF PROPENSITIES AND STOICHIOMETRIES OF THE DIFFERENT REACTIONS
	//////////////////////////////////////////////////	
	// This methods are used in the constructor of the class to define the reaction structures
	// both of them are declared as functions dependending on global parameters tha	t can be changed at runtime

	// For speed reasons x[0], x[1] and x[2] are assigned numerally, it would be interesting to access them by
	// their names, nevertheless it may cost computationally, it is worth to check it in the future creating a function
	// void assign (A,B, ) or a check function that checks that x[0].name==N etc.

	// The variables in the propensity functions are given in number of molecules (X=x*Omega)

	// x[0] will be mRNA and p[0] protein

	// x[0] mRNA
	// x[1] protein
	// x[2] cell state

	// mRNA production
	vector<int> m_pro_r(){ return {mRNAgeo,0,0};};
	double m_pro_0(v_species& x){
		if (x[2]->n<presynnum){
  			return k0;
  		} 
  		else{
  			return 2*k0; 
  		}
	};

	// mRNA degradation
	vector<int> m_deg_r(){ return {-1,0,0};};
	double m_deg(v_species& x){
  		return (x[0]->n)*k1;
	};

	// Deactivation of Olig2 promoter
	vector<int> p_pro_r(){ return {0,1,0};};
	double p_pro(v_species& x){
  		return x[0]->n*k2;
	};

	// Deactivation of Irx3 promoter
	vector<int> p_deg_r(){ return {0,-1,0};};
	double p_deg(v_species& x){
  		return x[1]->n*k3;
	};

	int main(int argc, char* argv[]){

		string inputfilename,outputfilename; 

		if(argc==3){
			inputfilename=argv[1];
			outputfilename=argv[2];
		}
		else {
			cout<<"Error: Two files (input,output) are mandatory!\n";	
			exit(0);
		}
		/////////////////////////
		// READING THE INPUT FILE
		////////////////////////////////////////

		// Opening input file
		ifstream inputfile;
		inputfile.open(inputfilename);

		// Parameters reading from input file
		string stringline;
		string par_name;
		string par_value;

		// initi random variable
		gsl_rng * rngz; // allocator for the rng generator
		rngz= gsl_rng_alloc (gsl_rng_mt19937);
		gsl_rng_set (rngz,::time(NULL)*getpid());

		while(getline(inputfile,stringline)){
			if (stringline[0]!='#' && stringline!=""){// If not a comment or blank line
				istringstream streamline(stringline);
				streamline>>par_name>>par_value;
				if(par_name=="Omega"){Omega=atof(par_value.c_str());}
				else if(par_name=="k0"){k0=atof(par_value.c_str());}
				else if(par_name=="k1"){k1=atof(par_value.c_str());}
				else if(par_name=="k2"){k2=atof(par_value.c_str());}
				else if(par_name=="k3"){k3=atof(par_value.c_str());}
				else if(par_name=="T"){T=atof(par_value.c_str());}
				else if(par_name=="m0"){m0=atof(par_value.c_str());}
				else if(par_name=="p0"){p0=atof(par_value.c_str());}
				else if(par_name=="ffwrite"){ffwrite=atoi(par_value.c_str());}
				else if(par_name=="timelapse"){timelapse=atof(par_value.c_str());}
				else if(par_name=="totaltime"){totaltime=atof(par_value.c_str());}
				else if(par_name=="dt"){dt=atof(par_value.c_str());}
				else if(par_name=="runtype"){runtype=atoi(par_value.c_str());}
				else if(par_name=="SEED"){SEED=atoi(par_value.c_str());}
				else if(par_name=="runtimes"){runtimes=atoi(par_value.c_str());}
				else if(par_name=="stocycle"){stochasticcycle=atoi(par_value.c_str());}
				else if(par_name=="phasenumber"){phasenum=atoi(par_value.c_str());}
				else if(par_name=="mRNAgeo"){mRNAgeo=atoi(par_value.c_str());}
				else if(par_name=="presynthesis"){presynnum=atoi(par_value.c_str());}
				else if(par_name=="cellphaserates"){
				// cellphaserates can contain more than one rate	
					phaserates.push_back(atof(par_value.c_str()));
						cout<<"reading phase rate:"<<par_value.c_str()<<'\n';
					while (streamline>>par_value){
						phaserates.push_back(atof(par_value.c_str()));
						cout<<"reading phase rate:"<<par_value.c_str()<<'\n';
					}
				}
				else {cout<<"Error: Reading Unknown Parameter "<<par_name<<" !\n"; exit(1);}
			}// End of line reading
		}// End of while for reading lines
		
		inputfile.close();

		////// CREATING GILLESPIE INSTANCE
		langil cyclo(outputfilename,Omega,SEED); // No argumenTS for random seed

	    ////// SETTING OF INTIIAL CONDITIONS

		// Order of species as rest of the program (T,L,TS,V)
		cyclo.AddSpecies("mRNA",m0); // Careful with the order
		cyclo.AddSpecies("prot",p0); // Careful with the order
		cyclo.AddCellCycleSpecies();
		for (auto k: phaserates){
			cyclo.addCellPhase(1.0/k,stochasticcycle); // CONSTANT_PHASE OR EXPONENTIAL_PHASE
		}
		cyclo.PrintSummarySpecies();

		/////////////////////////
		// DEFINING REACTIONS 
		////////////////////////////////////////////////
		// Here the different reactions are created, the propensity functions
		// are dfined as global functions in the top of the file

		// add last argument ADIABATIC to make the reaction resilient to noise

		if (mRNAgeo < 0){ // Normal constitutive case
			mRNAgeo = 1;
			cyclo.Add_Reaction("mRNA production phase",m_pro_r(),m_pro_0);
			}
		else{
			cyclo.Add_Reaction("mRNA production phase",m_pro_r(),m_pro_0);
			cyclo.Set_GeometricReaction("mRNA production phase");
		}
		cyclo.Add_Reaction("mRNA degradation",m_deg_r(),m_deg);

		cyclo.Add_Reaction("protein production",p_pro_r(),p_pro);
		cyclo.Add_Reaction("protein degradation",p_deg_r(),p_deg);


		cyclo.SetRunType(runtype,dt); 
		cyclo.SetWriteState(ffwrite,timelapse);


		cyclo.Set_StoreDivisionTimes();

	
//		cyclo.Add_TimeAction(&aLtime);
		vector<StructDivisionState> last_divisions; // divisions of the last Run
		vector<StructDivisionState> queue_divisions; // queue of new cells to simulate
		// In queue divisions the time is the global time of the population whilst in last_divisions is local to the parent time
		// StructDivisionState current_cell = {time : 0, species0before: 0, species1before:0, species0after : m0,species1after :0}; // cell evaluated at each loop
		// StructDivisionState sister_cell = {time : 0, species0before : 0, species1before:0, species0after: 0, species1after: 0}; // first cell 
		StructDivisionState current_cell = {.time = 0, .species0before = 0, .species1before = 0, .species0after = m0, .species1after = 0}; // cell evaluated at each loop
		StructDivisionState sister_cell = {.time = 0, .species0before = 0, .species1before = 0, .species0after = 0, .species1after = 0}; // first cell 
		queue_divisions.push_back(current_cell);
		vector<double> mRNAvec; // vector with the population of mRNAs at end of run
		vector<double> phasesvec; // vector with the population of phases at end of run
		int counter_divisions = 0;

		gsl_rng * rng; // allocator for the rng generator
		rng= gsl_rng_alloc (gsl_rng_mt19937);
		gsl_rng_set (rng,::time(NULL)*getpid());



		while( (queue_divisions.empty() == false)){
			counter_divisions ++;

			//cout<<'\n'<<"Simulating cell "<<counter_divisions<<'\n';
			// Selecting the next cell to simulate
			current_cell = queue_divisions.back();
			fathertime = current_cell.time;
			cyclo.SetState("mRNA", current_cell.species0after); //Set Initial conditions
			//cyclo.SetDivisionHistory(current_cell.species0before,current_cell.species0after);
			// cyclo.SetState("prot", current_cell.species0after); //Set Initial conditions
			//Setting Initial conditions
			cyclo.SetState("cell_cycle", 0); //Set Initial conditions
			cyclo.SetTime(0);
			cyclo.Reset_HistoryDivision();
			// cout<<"rest of integration: "<<totaltime-current_cell.time<<'\n';
			// Run simulation and store the outputs of division and final state
			cyclo.Run(totaltime-current_cell.time, false);
			last_divisions = cyclo.Get_HistoryDivision();
			queue_divisions.pop_back();
			mRNAvec.push_back(cyclo.GetState("mRNA"));
			phasesvec.push_back(cyclo.GetState("cell_phase"));
			// cout<<"Initial state "<<current_cell.species0after<<'\n';
			// cout<<"Final state "<<cyclo.GetState("mRNA")<<'\n';
			// cout<<"Results in the following divisions: "<<'\n';

			for (std::vector<StructDivisionState>::const_iterator i = last_divisions.begin(); i != last_divisions.end(); ++i){
				// cout<<"time_local:  "<<(*i).time<<" time_global "<<((*i).time + current_cell.time)<<" mRNAdaughter1 "<<((*i).species0after)<<" mRNAdaughter2 "<<((*i).species0before - (*i).species0after)<<'\n';
				sister_cell.time = (*i).time + fathertime;
				sister_cell.species0before = (*i).species0before; // not really necessary
				//sister_cell.species0after = gsl_ran_binomial(rng, 0.5, (*i).species0before); // not really necessary
				sister_cell.species0after = (*i).species0before-(*i).species0after; // remaining cellular content is passed to sister
				//sister_cell.species0after = 34; // remaining cellular content is passed to sister
				// Update the queue of divisions
				//cout<<"Resulting in divisions with mRNAS: "<<(*i).species0before-(*i).species0after <<' ';
				queue_divisions.push_back(sister_cell);
				//cout<<(*i).species0before<<'\n';
				//cout<<'\n';
				//cout<<"Temporary size of queue: "<<queue_divisions.size()<<"\n\n\n";
			}
			//cout<<"Current size of queue: "<<queue_divisions.size()<<"\n\n\n";
		}

		// cout<<"Finished population with mRNA population:";
		// for (std::vector<double>::const_iterator i = mRNAvec.begin(); i != mRNAvec.end(); ++i){
				// cout<<*i<<' ';
		// }
		// cout<<'\n';
		// cout<<"Finished population with phases:";
		// for (std::vector<double>::const_iterator i = phasesvec.begin(); i != phasesvec.end(); ++i){
				// cout<<*i<<' ';
		// }
		// cout<<'\n';				

	}
	


