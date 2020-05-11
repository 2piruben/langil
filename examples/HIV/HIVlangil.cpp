
/////////////////////////////////////////////////////////////////////////
//
// GILLESPIE SIMULATION CLASS OF BISTABLE SWITCH SYSTEM
//
//							R Perez-Carrasco
//							Created: 12 Aug '15
//////////////////////////////////////////////////////////////////////
//
// Class of the HIV-1 infection model from Sanchez-Taltavull switch. The species considered are
//   T : Healthy CD4+ T-cells, susceptible to be infected
//   L : Latently infected cells
//   TS : Productively inected cells
//   V: Viral load
//
// This class is derived from the Gillespie class, defining its own parameters and reactions
// reading them from a input file.

// The dependence with volume is introduced because the analyis is followed backwards from 
// the deterministic concentration macroscopic equation. Volumetric nature of each parameter should be stated and
// introduce it directly in the reaction equations to see the final Master Equation and 
// Mass action Equation

#include "langil.h"
#include <iostream>
#include <fstream>
#include <sstream>


#define DIM 4 //// Dimension, number of reactants (species)

using namespace std;

	const int HIVdim=DIM;
	// Main control parameters
	double Omega;// Volume of the system

	double lambdaT, deltaT; // production and degradation of healthy cells
	double k,eps,eta; // infection rate, Drug efficacy, fraction of resulting in L
	double deltaL,r,Lmax; // death rate of L, L homeostatis speed, carrying capacity of homeostatic prolif.

	double aL; // rate of activation from latency (will depend on parameters described bellow)
	double aLs, aLC, aLC2, aLT; // rate without drug, parameters controlling change in time

	double deltaTS; // death of active infected cells
	double lambdaV,deltaV; // virus production and clearance rate 

	int adiabaticmethod;
	// Initial and Final species number

	double T0,L0,TS0,V0; // Initial conditions
	double Tf,Lf,TSf,Vf; // Final conditions (if required)

	double DT0,DL0,DTS0,DV0; // Ball on Initial conditions
	double DTf,DLf,DTSf,DVf; // Ball on Final conditions (if required)

	// Model solving parameters

	int ffwrite; // kind of recording protocol (defined in .in or .h)
	double timelapse; // timestep for recording
	double totaltime;
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
	// void assign (A,B, ) or a check function that checks that x[0]->name==N etc.

	// The variables in the propensity functions are given in number of molecules (X=x*Omega)

	// Production of T 
	vector<int> Tpro_r(){ return {1,0,0,0};};
	double Tpro(v_species& x){
  		return Omega*lambdaT;
	};

	// Death of T
	vector<int> Tdeg_r(){ return {-1,0,0,0};};
	double Tdeg(v_species& x){
		return deltaT*x[0]->n; 
	};

	// Infection -> latency
	vector<int> Linf_r(){ return {-1,1,0,-1};};
	double Linf(v_species& x){
  		return eta*(1-eps)*k*x[0]->n*x[3]->n/Omega;
  	};

	// Infection -> active cell
	vector<int> TSinf_r(){ return {-1,0,1,-1};};
	double TSinf(v_species& x){
		return (1-eta)*(1-eps)*k*x[0]->n*x[3]->n/Omega; 
	};

	// Proliferation of L
	vector<int> Lpro_r(){ return {0,1,0,0};};
	double Lpro(v_species& x){
		return r*x[1]->n; 
	};

	// Death of L
	vector<int> Ldeg_r(){ return {0,-1,0,0};};
	double Ldeg(v_species& x){
		return deltaL*x[1]->n; 
	};

	// Homeostatic death of L
	vector<int> Ldeghom_r(){ return {0,-2,0,0};};
	double Ldeghom(v_species& x){
		return r/Lmax*x[1]->n*(x[1]->n-1)/2.0/Omega; 
	};

	// Activation L->TS
	vector<int> TSact_r(){ return {0,-1,1,0};};
	double TSact(v_species& x){
		return aL*x[1]->n; 
	};

	// TS death
	vector<int> TSdeg_r(){ return {0,0,-1,0};};
	double TSdeg(v_species& x){
		return deltaTS*x[2]->n; 
	};

	// V production from TS cells
	vector<int> Vpro_r(){ return {0,0,0,1};};
	double Vpro(v_species& x){
		return lambdaV*x[2]->n; 
	};

	// V death
	vector<int> Vdeg_r(){ return {0,0,0,-1};};
	double Vdeg(v_species& x){
		return deltaV*x[3]->n; 
	};

	// V infection failure
	vector<int> Vfail_r(){ return {0,0,0,-1};};
	double Vfail(v_species& x){
		return eps*k*x[0]->n*x[3]->n/Omega; 
	};


	// Change of aL in time
	double aLtime(double t){
		aL=aLs+aLC*(1-aLC2*cos(2*M_PI*t/aLT+M_PI));
		return aL;
	}

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
				else if(par_name=="lambdaT"){lambdaT=atof(par_value.c_str());}
				else if(par_name=="deltaT"){deltaT=atof(par_value.c_str());}
				else if(par_name=="k"){k=atof(par_value.c_str());}
				else if(par_name=="eps"){eps=atof(par_value.c_str());}
				else if(par_name=="eta"){eta=atof(par_value.c_str());}
				else if(par_name=="deltaL"){deltaL=atof(par_value.c_str());}
				else if(par_name=="r"){r=atof(par_value.c_str());}
				else if(par_name=="Lmax"){Lmax=atof(par_value.c_str());}
				else if(par_name=="aL"){aL=atof(par_value.c_str());}
				else if(par_name=="aLs"){aLs=atof(par_value.c_str());}
				else if(par_name=="aLC"){aLC=atof(par_value.c_str());}
				else if(par_name=="aLC2"){aLC2=atof(par_value.c_str());}
				else if(par_name=="aLT"){aLT=atof(par_value.c_str());}				
				else if(par_name=="deltaTS"){deltaTS=atof(par_value.c_str());}
				else if(par_name=="lambdaV"){lambdaV=atof(par_value.c_str());}
				else if(par_name=="deltaV"){deltaV=atof(par_value.c_str());}
				else if(par_name=="adiabaticmethod"){adiabaticmethod=atoi(par_value.c_str());}				
				else if(par_name=="T0"){T0=atof(par_value.c_str());}
				else if(par_name=="L0"){L0=atof(par_value.c_str());}
				else if(par_name=="TS0"){TS0=atof(par_value.c_str());}
				else if(par_name=="V0"){V0=atof(par_value.c_str());}
				else if(par_name=="Tf"){Tf=atof(par_value.c_str());}
				else if(par_name=="Lf"){Lf=atof(par_value.c_str());}
				else if(par_name=="TSf"){TSf=atof(par_value.c_str());}
				else if(par_name=="Vf"){Vf=atof(par_value.c_str());}
				else if(par_name=="DT0"){DT0=atof(par_value.c_str());}
				else if(par_name=="DL0"){DL0=atof(par_value.c_str());}
				else if(par_name=="DTS0"){DTS0=atof(par_value.c_str());}
				else if(par_name=="DV0"){DV0=atof(par_value.c_str());}
				else if(par_name=="DTf"){DTf=atof(par_value.c_str());}
				else if(par_name=="DLf"){DLf=atof(par_value.c_str());}
				else if(par_name=="DTSf"){DTSf=atof(par_value.c_str());}
				else if(par_name=="DVf"){DVf=atof(par_value.c_str());}
				else if(par_name=="ffwrite"){ffwrite=atoi(par_value.c_str());}
				else if(par_name=="timelapse"){timelapse=atof(par_value.c_str());}
				else if(par_name=="totaltime"){totaltime=atof(par_value.c_str());}
				else if(par_name=="dt"){dt=atof(par_value.c_str());}
				else if(par_name=="runtype"){runtype=atoi(par_value.c_str());}
				else if(par_name=="SEED"){SEED=atoi(par_value.c_str());}
				else {cout<<"Error: Reading Unknown Parameter "<<par_name<<" !\n"; exit(1);}
			}// End of line reading
		}// End of while for reading lines
		
		inputfile.close();

		////// CREATING GILLESPIE INSTANCE
		langil HIV(outputfilename,Omega,SEED); // No argumenTS for random seed

	    ////// SETTING OF INTIIAL CONDITIONS

		// Order of species as rest of the program (T,L,TS,V)
		HIV.AddSpecies("T",T0*Omega); // Careful with the order
		HIV.AddSpecies("L",L0*Omega);
		HIV.AddSpecies("TS",TS0*Omega);
		HIV.AddSpecies("V",V0*Omega);

		//cout<<"Starting Integration point   ";NOP.PrinTState();

		/////////////////////////
		// DEFINING REACTIONS FOR NOIP SYSTEM
		////////////////////////////////////////////////
		// Here the different reactions are created, the propensity functions
		// are dfined as global functions in the top of the file

		// add last argument ADIABATIC to make the reaction resilient to noise

		HIV.Add_Reaction("T Production",Tpro_r(),Tpro);
		HIV.Add_Reaction("T Death",Tdeg_r(),Tdeg);
		HIV.Add_Reaction("T infection into L",Linf_r(),Linf);
		HIV.Add_Reaction("T infection into TS",TSinf_r(),TSinf);
		HIV.Add_Reaction("L Production",Lpro_r(),Lpro);
		HIV.Add_Reaction("L death",Ldeg_r(),Ldeg);
		HIV.Add_Reaction("L homeostatic death",Ldeghom_r(),Ldeghom);
		HIV.Add_Reaction("L activation into TS",TSact_r(),TSact);
		HIV.Add_Reaction("TS death",TSdeg_r(),TSdeg);
		HIV.Add_Reaction("V production from TS",Vpro_r(),Vpro);
		HIV.Add_Reaction("V death",Vdeg_r(),Vdeg);
		HIV.Add_Reaction("V infection failure",Vfail_r(),Vfail);

		// adiabaticmethod 4 is the default and includes all the reactions
		if (adiabaticmethod>4 || adiabaticmethod<1){
			cout<<" Wrong adiabatic method: "<<adiabaticmethod;
		}
		if (adiabaticmethod<4){ // T is considered constant
			HIV.Set_AdiabaticReaction("T Production");
			HIV.Set_AdiabaticReaction("T Death");
			if (adiabaticmethod<3){ // V is considered adiabatic
				HIV.Set_AdiabaticReaction("V production from TS");
				HIV.Set_AdiabaticReaction("V death");
				HIV.Set_AdiabaticReaction("V infection failure");
				if (adiabaticmethod<2){ // L is considered adiabatic
					HIV.Set_AdiabaticReaction("TS death");
				}
			}
		}

		HIV.SetRunType(runtype,dt); 
		HIV.SetWriteState(ffwrite,timelapse);
	
//		HIV.Add_TimeAction(&aLtime);

		HIV.Run(totaltime);
/*
		species Tinit("T",T0*Omega); species Linit("L",L0*Omega);
		species TSinit("TS",TS0*Omega); species Vinit("V",V0*Omega); 

		v_species initialstate;
		initialstate.push_back(Tinit); initialstate.push_back(Linit);
		initialstate.push_back(TSinit); initialstate.push_back(Vinit);

		species Tend("T",Tf*Omega); species Lend("L",Lf*Omega);
		species TSend("TS",TSf*Omega); species Vend("V",Vf*Omega); 

		v_species finalstate;
		finalstate.push_back(Tend); finalstate.push_back(Lend);
		finalstate.push_back(TSend); finalstate.push_back(Vend);

		species Tinitdist("T",DT0*Omega); species Linitdist("L",DL0*Omega);
		species TSinitdist("TS",DTS0*Omega); species Vinitdist("V",DV0*Omega); 

		v_species initialstatedist;
		initialstatedist.push_back(Tinitdist); initialstatedist.push_back(Linitdist);
		initialstatedist.push_back(TSinitdist); initialstatedist.push_back(Vinitdist);

		species Tenddist("T",DTf*Omega); species Lenddist("L",DLf*Omega);
		species TSenddist("TS",DTSf*Omega); species Venddist("V",DVf*Omega); 

		v_species finalstatedist;
		finalstatedist.push_back(Tenddist); finalstatedist.push_back(Lenddist);
		finalstatedist.push_back(TSenddist); finalstatedist.push_back(Venddist);

		cout<<"Time step "<<dt<<'\n';
		cout.flush();
		HIV.RunTransition(initialstate, finalstate, initialstatedist, finalstatedist, 1,totaltime);
*/
		}
	


