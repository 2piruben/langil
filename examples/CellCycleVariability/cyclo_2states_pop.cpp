
/////////////////////////////////////////////////////////////////////////
//
// GILLESPIE SIMULATION CLASS OF CYCLOSTATIONARY PARTITION
//
//							R Perez-Carrasco
//							Created: 4 May '19
//////////////////////////////////////////////////////////////////////
//

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
	double T; // average cell cycle time duration

	int stochasticcycle; // flag for cell cycle stochasticity
	int phasenum; // number of phases per cell cycle
	int presynnum;  // number of cell cycle phases before gene doubling
	vector<double> phaserates; // rates of the different cell cycle phases

	double m0,p0; // Initial conditions

	int mRNAgeo; // flag to determine if mRNA is produced geometrically (-1 is constitutive, >0 is the burst average size)


	// Model solving parameters

	int ffwrite; // kind of output
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
	// These methods are used in the constructor of the class to define the reactions
	
	// x[0] mRNA
	// x[1] protein
	// x[2] cell state
	// cell state ractions are managed directly with the class cell_cycle.h

	// mRNA production
	vector<int> m_pro_r(){ return {mRNAgeo,0,0};};
	double m_pro_0(v_species& x){
		if (x[2]->n<presynnum){
  			return k0;
  		} 
  		else{
  			return 2*k0;  // post-replication
  		}
	};

	// mRNA degradation
	vector<int> m_deg_r(){ return {-1,0,0};};
	double m_deg(v_species& x){
  		return (x[0]->n)*k1;
	};

	// protein production
	vector<int> p_pro_r(){ return {0,1,0};};
	double p_pro(v_species& x){
  		return x[0]->n*k2;
	};

	// protein degradation
	vector<int> p_deg_r(){ return {0,-1,0};};
	double p_deg(v_species& x){
  		return x[1]->n*k3;
	};

//////////////////////////// MAIN /////////////////////////////////////////////////

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

		////// CREATING LANGIL INSTANCE
		langil cyclo(outputfilename,Omega,SEED); // No argumenTS for random seed

	    ////// ADDING SSPECIES 

		cyclo.AddSpecies("mRNA",m0); // Careful with the order
		cyclo.AddSpecies("prot",p0); // Careful with the order
		cyclo.AddCellCycleSpecies();
		for (auto k: phaserates){
			cyclo.addCellPhase(1.0/k,stochasticcycle); // // stochasticcycle can be CONSTANT_PHASE OR EXPONENTIAL_PHASE
		}
		cyclo.PrintSummarySpecies();

		///// ADDING REACTIONS

		if (mRNAgeo < 0){ // constitutive case
			mRNAgeo = 1; // number of mRNAs produced per reaction
			cyclo.Add_Reaction("mRNA production phase",m_pro_r(),m_pro_0);
			}
		else{ // bursty
			cyclo.Add_Reaction("mRNA production phase",m_pro_r(),m_pro_0);
			cyclo.Set_GeometricReaction("mRNA production phase");
		}
		cyclo.Add_Reaction("mRNA degradation",m_deg_r(),m_deg);

		cyclo.Add_Reaction("protein production",p_pro_r(),p_pro);
		cyclo.Add_Reaction("protein degradation",p_deg_r(),p_deg);

		//// DEFINING RUN TYPE AND OUTPUT
		cyclo.SetRunType(runtype,dt); 
		cyclo.SetWriteState(ffwrite,timelapse);


		/// CREATING STORAGE TO TRACK DIVIDING POPULATION
		cyclo.Set_StoreDivisionTimes();
		vector<StructDivisionState> last_divisions; // divisions of the last trajectory
		vector<StructDivisionState> queue_divisions; // queue of new cells to simulate from previous trajectory
		StructDivisionState current_cell = {.time = 0, .species0before = 0, .species1before = 0, .species0after = m0, .species1after = 0}; // cell evaluated at each loop
		StructDivisionState sister_cell = {.time = 0, .species0before = 0, .species1before = 0, .species0after = 0, .species1after = 0}; // first cell 
		queue_divisions.push_back(current_cell);
		vector<double> mRNAvec; // vector with the population of mRNAs at end of run
		vector<double> phasesvec; // vector with the population of phases at end of run
		int counter_divisions = 0;

		while( (queue_divisions.empty() == false)){ // keep dividing as long as there are cells younger than the final time
			counter_divisions ++;
			// Selecting the next cell to simulate
			current_cell = queue_divisions.back();
			fathertime = current_cell.time;
			cyclo.SetState("mRNA", current_cell.species0after); //Set Initial mRNA from mother cell
			cyclo.SetState("cell_cycle", 0); // Initiate cell cycle
			cyclo.SetTime(0);
			cyclo.Reset_HistoryDivision();
			cyclo.Run(totaltime-current_cell.time, false); // Run a trajectory until final time and store untracked daughters
			last_divisions = cyclo.Get_HistoryDivision();
			queue_divisions.pop_back();
			mRNAvec.push_back(cyclo.GetState("mRNA"));
			phasesvec.push_back(cyclo.GetState("cell_phase"));

			for (std::vector<StructDivisionState>::const_iterator i = last_divisions.begin(); i != last_divisions.end(); ++i){
				// loop to add each untracked sister to teh queue of trajectories to simulate
				sister_cell.time = (*i).time + fathertime;
				sister_cell.species0before = (*i).species0before; // not really necessary
				sister_cell.species0after = (*i).species0before-(*i).species0after; // remaining cellular content is passed to sister
				// Update the queue of divisions
				queue_divisions.push_back(sister_cell);
			}
		}
	}
	


