///////////////////////////////////////////////////////
//////   LANGIL CLASS
//////
//////						R. Perez-Carrasco
//////                               Created: 18 Nov '13
///////////////////////////////////////////////////////
/* Class for simulating master equation through Gillespie algorythm or Langevin equation.
 It uses the reaction class that defines each reaction taking place in the simulation 
*/

#include "langil.h"

	// Constructor of the langil class	
	langil::langil(string file,double vol, int seed){
		rng= gsl_rng_alloc (gsl_rng_mt19937);
		if(seed==-1){ // if seed==-1 (default) take a random number
			gsl_rng_set (rng,::time(NULL)*getpid());
			cout<<"RNG Seed used: "<<::time(NULL)*getpid()<<'\n';
		}
		else{
			gsl_rng_set (rng,seed);
		}
		fwrite=ALL; // by default write reaction after each step
		trajfile.open(file); // output trajectory
		nextrectime=0;
		time = 0;
		signal = 0.0;
		fStoreDivisionTimes = false;
		cellcyclespeciesidx = -1; // negative means that there will be no cell phases
		MakeStep= &langil::GillespieStep; // by deafult gillespie algorythm is used
		rnd.resize(2); // by default only two elements are necessary in Gillespie Simulation
		Omega=vol; //Setting Volume of the system (relation between concentration and numbre of molecules)
	}
		
	void langil::AddSpecies(string aname, float anum){
			x.push_back(p_species(new species(aname,anum)));
			xtemp0.push_back(p_species(new species(aname,anum)));
			xtemp1.push_back(p_species(new species(aname,anum)));
	}


	void langil::AddSpeciesTimeTracking(string aname, float anum, double lifetime){
			langil::AddSpecies(aname, anum); // same as previous function
			x.back()->MakeSpeciesTimeTracking(lifetime); // but know we also set the time and flag for time tracking
			xtemp0.back()->MakeSpeciesTimeTracking(lifetime);
			xtemp1.back()->MakeSpeciesTimeTracking(lifetime);
	}

	void langil::AddCellCycleSpecies(){
			x.push_back(p_species(new species("cell cycle",0)));
			xtemp0.push_back(p_species(new species("cell cycle",0)));
			xtemp1.push_back(p_species(new species("cell cycle",0)));
			cellcyclespeciesidx = x.size()-1;
	}

	void langil::addCellPhase(double duration, int type_phase){
		cell.addCellPhase(duration,type_phase);
	}

	void langil::PrintSummarySpecies(){
		cout<<"Chemical Species\n";
		for(auto& s: x){
			cout<< s->GetName() <<  ':' << s->n << " with next time: " << s->GetNextTime() << '\n';
		}
		cout<<"Cell state\n";
		cell.printCellProperties();
	}

	void langil::SetState(string aname, float anum){
		for(auto& s: x){
			if((s->GetName()).compare(aname) == 0)  { // strcmp  compares strings and returns 0 if equal 
				s->SetNum(anum);
				return;
			}
		}
		if(aname.compare("cell_cycle") == 0){
			// cout<<"Resetting cell cycle"<<'\n';
			cell.ResetCellCycle();
			x.back()->SetNum(0);
			return;
		}
		cout<<"Error: Non-existing species to be accessed by SetState() \n";
		exit(0);
	}

	double langil::GetState(string aname){
		for(auto& s: x){
			if((s->GetName()).compare(aname) == 0)  { // strcmp  compares strings and returns 0 if equal 
				return s->GetNum();
			}
		}
		if(aname.compare("cell_phase") == 0){
			return cell.getCellPhase();
		}
		cout<<"Error: Non-existing species to be accessed by SetState() \n";
		exit(0);
	}

	void langil::SetTime(double t){
		time = t;
	}

	string langil::StateString(){
		string statestring;
		for (auto& s : x){
			 statestring+=s->GetName()+": "+to_string(s->GetNum())+"  ";
		}
		return statestring;
	}


	void langil::Add_Reaction(string name,vector<int> stoich,double (*prop_f)(v_species&),int reacttype){
		reaction rec(name,stoich,prop_f,reacttype);
		rec.SetState(x);
		r.push_back(rec);
		cout<<"Added Reaction: "<<rec.GetName()<<" "<<rec.IsAdiabatic()<<"\n";
	}

	void langil::Set_AdiabaticReaction(string aname){
		for(reaction& s: r){
			if(s.GetName()==aname){ 
				s.SetAdiabatic();
				cout<<"Setting adiabatic reaction: "<<s.GetName()<<'\n';
				break;
			}
//	else{
//			cout<<"Error: Non-existing reaction to be accessed by SetState(): "<<aname<<"\n";
//			exit(0);
//			}
		}
	}

	void langil::Set_GeometricReaction(string aname){
		for(reaction& s: r){
			if(s.GetName()==aname){ 
				s.SetGeometric();
				cout<<"Setting geometric stoichiometry reaction: "<<s.GetName()<<'\n';
				break;
			}
//	else{
//			cout<<"Error: Non-existing reaction to be accessed by SetState(): "<<aname<<"\n";
//			exit(0);
//			}
		}
	}


	void langil::SetRunType(int rt, double adt){
		switch (rt){
			case (GILLESPIE):
				cout<<"Preparing Gillespie algorithm...\n";
				MakeStep=&langil::GillespieStep;
				rnd.resize(2);
				xtemp.resize(x.size()); // Euler integration requires a temporal vector
				dt = adt;
				break;
			case (MACROSCOPIC):
				MakeStep=&langil::MacroStep;
				rnd.resize(0);
				xtemp.resize(x.size()); // Euler integration requires a temporal vector
				dt=adt;
				break;
			case (CLEEULER):
				MakeStep=&langil::LangevinEulerStep;
				rnd.resize(r.size()); // CLE requires one random per reaction
				xtemp.resize(x.size()); // Euler integration requires a temporal vector
				cout<<"Integration CLE timestep: "<< SetLangevinTimeStep(adt)<<'\n'; // Compute optimum time with characteritic advance adx
//				dt=adt;
				break;
			case (CLEMILSTEIN):
				MakeStep=&langil::LangevinMilsteinStep;
				rnd.resize(r.size()); // CLE requires one random per reaction
				cout<<"Integration CLE timestep: "<< SetLangevinTimeStep(adt)<<'\n'; // Compute optimum time with characteritic advance adx
//				dt=adt;
				break;
			default:
				cout<<"Error: Integration Method Non-existing";
				exit(0);
		}
	}

	double langil::GillespieStep(double timelimit){

		// cout<<"Making Gillespie Step"<<'\n';
		rnd[0]=gsl_rng_uniform_pos (rng);
		rnd[1]=gsl_rng_uniform_pos (rng);

		totalprop=0;
		detprop=-1;
		celldivided = 0;

		for(unsigned int i=0;i<x.size();i++){ // update deterministic species times
			xtemp0[i]->n = x[i]->n; // memorize the last step
		}

		// propensity sum calculation
		for(vector<reaction>::iterator it = r.begin(); it != r.end(); ++it) {
			 // cout<<(it->GetName())<<' '<<(it->IsAdiabatic())<<"\n";
			 // cout.flush();
			if (it->IsNotAdiabatic() && it->IsNotDetermTime() ){
			  // cout<<"Propensity of "<<it->GetName()<<':'<<it->GetPropensity()<<'\n';
			  // cout.flush();
				totalprop+= it->GetPropensity();
			}
			else if (it->IsDetermTime() && it->GetPropensity()>0){// if Propensity is negative in a deterministic reaction, means that it has no elements
				if (detprop<0) {
					detprop = it->GetPropensity();
					// cout<<"Det. Propensity of "<<it->GetName()<<':'<<it->GetPropensity()<<'\n';
					nextreaction = it; // by default, next reaction pointer poitns to the first possible deterministic time reaction it founds
				}
				else if (it->GetPropensity()>detprop){
						detprop = it->GetPropensity();
						nextreaction = it; //  if a faster deterministic event is coming, point to it
				}
			}
		}

		nexttau=(1./totalprop)*log(1./rnd[0]); // next stochastic reaction time lapse
		// cout<<"Comparing detprop with stoprop:  "<<detprop<<"   "<<1.0/nexttau<<'\n';
		// cout<<"Next cell cycle propensity: "<<1.0/cell.getTimeForNextCellPhase()<<'\n';
		// cout<<"Next cell phase event: "<<cell.getTimeForNextCellPhase()<<'\n';

		//vector<double> vectordivision;
		if (cellcyclespeciesidx >= 0 && (1.0/cell.getTimeForNextCellPhase()>1.0/nexttau) && (1.0/cell.getTimeForNextCellPhase() > 1.0/detprop) ){
		// if the next event is a cell cycle phase event
			nexttau = cell.getTimeForNextCellPhase();
			if(nexttau>timelimit){
				return nexttau; // if time is greater than last allowed time, return before updating species
			}
			x[cellcyclespeciesidx] -> SetNum(cell.AdvanceNextPhase());
			// cout<<"Cell cycle advanced to "<<cell.getCellPhase()<<'\n';
			if (cell.MaterialIsDivided()) {
				celldivided = 1;
				if (fStoreDivisionTimes){
					divisionstate.time = time+nexttau;
					divisionstate.species0before = x[0]->GetNum();
					divisionstate.species1before = x[1]->GetNum();
					// cout<<"mRNA before division: "<<x[0]->GetNum()<<'\n';
				}
			// advance to next state, and if it is the start raise the division flag
			// the celldivided flag is used in the output file to mark cell division

			// cout<<"Dividing cell at time "<<time+nexttau<<"\n";
				for(unsigned int i=0;i<x.size();i++){ // update new variables
						if (i != cellcyclespeciesidx){
							x[i]->SetNum(gsl_ran_binomial(rng, 0.5, x[i]->GetNum())); 
							//x[i]->SetNum(x[i]->GetNum()/2.0);
							//x[i]->SetNum(0);
						}
				}
				if (fStoreDivisionTimes){
					divisionstate.species0after = x[0]->GetNum();
					divisionstate.species1after = x[1]->GetNum();
					DivisionHistoryVector.push_back(divisionstate);
					// cout<<"mRNA after division: "<<x[0]->GetNum()<<'\n';
				}
			}
			for(unsigned int i=0;i<x.size();i++){ // update deterministic species times
				x[i]->Updatetime(nexttau);
			}
		}

		else if(detprop>1.0/nexttau){
			// If the next reaction is a determinisitc one
				nextstoichiometry=(nextreaction->GetStoichiometry());
				// cout<<"Next reaction is "<< nextreaction->GetName()<<'\n';
				nexttau = 1.0/detprop;
				if(nexttau>timelimit){
					return nexttau;
				}
				// cout<<"Deterministic reaction with nexttau "<<nexttau<<'\n';
				for(unsigned int i=0;i<x.size();i++){ // update deterministic times
					x[i]->Updatetime(nexttau);
				}
				for(unsigned int i=0;i<x.size();i++){ // update new variables
					x[i]->React((*nextstoichiometry)[i]);
				}
				if(cellcyclespeciesidx >=0){
				cell.updateCelltime(nexttau);
				}		
		}

		else{ // if the next time is a stochastic time, then follow on with the 
		// selection of reaction
			if(nexttau>timelimit){
				cell.updateCelltime(timelimit);		
				return nexttau;
			}
			for(unsigned int i=0;i<x.size();i++){ // update deterministic species times
				x[i]->Updatetime(nexttau);
			}

			cell.updateCelltime(nexttau);		
			cumprop=0;
			for(vector<reaction>::iterator it = r.begin(); it != r.end(); ++it) {
				if (it->IsNotAdiabatic() && it->IsNotDetermTime()){
	//				cout<<"cumprop"<<"\n";
	//				cout.flush();
					cumprop+=it->GetPropensity();
					if(cumprop >= totalprop*rnd[1]){
						nextreaction=it;
						break;
					}
				}
			}
			// cout<<"Next reaction is "<< nextreaction->GetName()<<"\n\n";
			// Actualization of the state of the system
			nextstoichiometry=(nextreaction->GetStoichiometry());	
			if (nextreaction->IsNotGeometric()){
				for(unsigned int i=0;i<x.size();i++){
					x[i]->React((*nextstoichiometry)[i]);
				}
			}
			else { // if it is geometric
				// cout<<"A geometric reaction is happening..."<<'\n';
				for(unsigned int i=0;i<x.size();i++){
					// The probability of fail for a burst size b is 1/(1+b)
					if ((*nextstoichiometry)[i] > 0){
						geostoichiometry = gsl_ran_geometric(rng,1.0/(1.0+(*nextstoichiometry)[i]))-1;
					// the -1 is for the convention that k=0 successses is possible
						// cout<<"Species "<< i<<" with sto "<<(*nextstoichiometry)[i]<<" is increased "<< geostoichiometry<<'\n';
						x[i]->React(geostoichiometry);
					}
				}
			}
		}
		MacroAdiabaticStep(nexttau);
	//		time+=nexttau;
		// cout<<"Time to reaction: "<<nexttau<<"   (Time: "<< time+nexttau <<")\n";

		// update time fror deterministic reactions


		// cout<<"State after reaction: "<<'\n';
		// for(unsigned int i=0;i<x.size();i++){
			// cout<< x[i]->GetNum() << ' ';
			// }
		// cout<<'\n';
		

		return nexttau;

	}

	double langil::LangevinEulerStep(double timelimit){
//		cout<<time<<' '<<dt<<' '<<StateString()<<'\n';
		for(uint i=0;i<rnd.size();i++) rnd[i]=gsl_ran_gaussian_ziggurat (rng,1.0);
		for(uint i=0;i<x.size();i++){ // for each species
			cout<<x[i]->GetName()<<" = "<< x[i]->n <<"\n";
			xtemp[i]=x[i]->n;
			for(uint j=0;j<r.size();j++){ // for each reaction
				sto = (*(r[j].GetStoichiometry()))[i]; 
				if (abs(sto)>0){
					cout<<r[j].GetName()<<":"<<r[j].GetPropensity()<<"\n";
				}
				if(sto!=0){
					xtemp[i]+=sto*r[j].GetPropensity()*dt; // deterministic part
					xtemp[i]+=abs(sto)*sqrt(r[j].GetPropensity()*dt)*rnd[j]; //stochastic
				}
			}
		}
		cout<<"Proposed step:\n";
		for(uint i=0;i<x.size();i++){ // for each species
			xtemp0[i]->n=xtemp[i]; //casting from array to vector of species
			cout<<x[i]->GetName()<<' '<<xtemp0[i]->n<<'\n';
		}
		//Boundary_Behaviour(x,xtemp0); // Apply defined rules when system is out of boundary
		// boundary behaviour is still not fully implemented
		cout<<"Corrected step:\n";
		for(uint i=0;i<x.size();i++){
			x[i]->n=xtemp0[i]->n; // Update
			cout<<x[i]->GetName()<<' '<<xtemp0[i]->n<<'\n';
		}
		cout<<"\n\n\n";
		return dt;
	}


	double langil::LangevinMilsteinStep(double timelimit){
//		cout<<time<<' '<<dt<<' '<<StateString()<<'\n';
		// First predictor step is deterministic (this substitutes an explicit derivative calculation)
		for(uint i=0;i<x.size();i++){ // for each species
			xtemp0[i]->n=x[i]->n;
			for(uint j=0;j<r.size();j++){
//					cout<<r[j].GetName()<<":"<<r[j].stoichiometry[i]<<"     ";
				sto = (*(r[j].GetStoichiometry()))[i];
				if(sto!=0){
					pro = r[j].GetPropensity();  
					xtemp0[i]->n += sto*pro*dt; // deterministic part
					xtemp0[i]->n += abs(sto)*sqrt(pro*dt); //stochastic
				}
			}
		}
		// 
		for(uint i=0;i<rnd.size();i++) rnd[i]=gsl_ran_gaussian_ziggurat (rng,1.0);
		for(uint i=0;i<x.size();i++){ // for each species
			xtemp1[i]->n = x[i]->n;
			for(uint j=0;j<r.size();j++){
//					cout<<r[j].GetName()<<":"<<r[j].stoichiometry[i]<<"     ";
				sto = (*(r[j].GetStoichiometry()))[i]; 
				if(sto!=0){
					pro = r[j].GetPropensity();
					pro0 = r[j].GetPropensity(xtemp0);
					xtemp1[i]->n += sto*pro*dt; // deterministic part
					xtemp1[i]->n += abs(sto) * sqrt(pro*dt)*rnd[j];  // first order noise
					xtemp1[i]->n += abs(sto) * 0.5 * (sqrt(pro0)-sqrt(pro))*sqrt(dt)*(rnd[j]*rnd[j]-1.0);
					if (xtemp1[i]->n<0) xtemp1[i]->n =0;
					cout<<"pro "<<pro<<" pro0 "<<pro0<<'\n';
				}
			}
		}
//		cout<<'\n';
		for(uint i=0;i<x.size();i++) x[i]->n=xtemp1[i]->n; // Update
		return dt;
	}



	double langil::MacroStep(double timelimit){
		for(uint i=0;i<x.size();i++){ // for each species
			xtemp[i]=x[i]->n;
			for(uint j=0;j<r.size();j++){
				sto = (*(r[j].GetStoichiometry()))[i]; 
				if(sto!=0){
					xtemp[i]+=sto*r[j].GetPropensity()*dt; // deterministic part
				}
			}
		}
		for(uint i=0;i<x.size();i++) x[i]->n=xtemp[i]; // Update
		return dt;
	}

	double langil::MacroAdiabaticStep(double dtime){
		for(uint i=0;i<x.size();i++){ // for each species
			xtemp[i]=x[i]->n;
			for(uint j=0;j<r.size();j++){ // for each adiabatic reaction
				sto = (*(r[j].GetStoichiometry()))[i];
				if (r[j].IsAdiabatic()){
					if(sto!=0){
						xtemp[i]+=sto*r[j].GetPropensity()*dtime; // deterministic part
					}
				}
			}
		}
		for(uint i=0;i<x.size();i++) x[i]->n=xtemp[i]; // Update
		return dtime;
	}


	double langil::Run(double T,bool verbose){
		
		if(verbose){
			cout<<"Running Integration algorithm...";
			cout<<"Total time:"<<T<<"  Initial time:"<<time<<'\n';
			cout<<"Initial State:"<<x[0]->n<<' '<<x[1]->n<<' '<<'\n';
			cout<<"ffwrite : "<<fwrite<<'\n';
			if (cellcyclespeciesidx>=0){cell.PrintState();};
		}

		while(time<T){
			
//		cout<<t<<' '<<T<<'\n';
			RunTimeAction();
			time+=(this->*MakeStep)(T-time); // T-time is the time limit

			if(fwrite==ALL) WriteState();
			else if(fwrite==TIMESTEP){
				while(time>nextrectime){
					nextrectime+=timelapse;
					WritePrevrecState();
				}
			}
			else if(fwrite==TIMESTEPMULT){
				while(time>nextrectime){
					nextrectime=timelapse*(nextrectime+dt)-dt;
					// This returns recording times following formula (T(i)=dt*(timelapse^i-1))
					// This has the properties T(0)=0, T(1)=dt and T(i>>1)=timelapse*T(i-1)
					// It is useful for exponential variations in time events
					WritePrevrecState();
				}
			}
		}
		if(fwrite == LASTPOINT) WriteState(); 
//		cout<<t<<' '<<T<<'\n';
		if(verbose){
			cout<<"Done!"<<'\n';
			cout<<"Final state: (";
			cout<<StateString();
			cout<<")\n";
		}
		return time;
	}


	double langil::RunTimes(int N, double T,bool verbose){
// Runs N integrations and records the last point of the integration

		v_species xinit;
		for (auto &s: x){
			xinit.push_back(p_species(s->clone())); // xinit = x
		}
		
		if(verbose){
			cout<<"Running Integration algorithm...";
		}

		for(int Ncount=0;Ncount<N;Ncount++){
			nextrectime=0;
			time=-1; // negative times are used to encode actions that only occur at the start of a simulation
			// such as chosen parameters in random way
			// this is done this way to avoid comparison with time == 0
			RunTimeAction();

			time=0;
			while(time<T){		
				RunTimeAction();
				time+=(this->*MakeStep)(T-time);
			}
		
			WriteTempState();
			nextrectime=0;
			time=0;

			for(uint i=0;i<x.size();i++){//  x = xinit
				*x[i] = *xinit[i];
				*xtemp0[i] = *xinit[i];
				x[i]->determtimeque.clear(); // clear any deterministic time in the queue

			}
		}

//		cout<<t<<' '<<T<<'\n';
		if(verbose){
			cout<<"Done!"<<'\n';
		}
		return time;
	}

	double langil::RunTransition(bool (*inside_init)(v_species&), bool (*inside_end)(v_species&), 
							     int loopcount,double maxT){ // loopcount=1 also records the transition trajectory

		clock_t start,end;

		bool maxTreached=false;	
		bool recording=false;// Is now x being recorded in xtran?
		bool nearfromend=false; // Is x close to the initstate?
		bool nearfrombeginning=true; // Is x close to the endstate?
		trajfile.precision(10);
		double nextrecordingtime = 0; 
	

		v_species xinit; // saving initial position to restart the system after each transition	
		for (auto &s: x){
			xinit.push_back(p_species(s->clone())); // xtran.push_back(x)
		}
		// The following two could be joined in a compact way, but I will not complicate the code
		vector< v_species > xtran; // transition trajectory is a vector of species vectors
		vector<double> ttran; // transition trajecotry times		

		trajfile<<"# Simulation time     Computation time (seconds)\n";
		trajfile.flush();



		for(int i=0;i<loopcount;i++){
			start=::clock();
			time=0;
			for(uint i=0;i<x.size();i++){//  x = xinitstate
				*x[i] = *xinit[i];
			}	
			nearfromend=false;
			nearfrombeginning=true;
			maxTreached=false;	
			xtran.clear();
			ttran.clear();
//			int count=0;

			while(!nearfromend){ // The integration will end when nearfromend=true;
				if (time>maxT){ // Or the maximum time is reached
					maxTreached=true;
					break;
				} 
				time+=(this->*MakeStep)(maxT-time);
	 			nearfrombeginning=true;
				for(uint d=0;d<x.size();d++){
					if( nearfrombeginning && inside_init(x) ) {
						nearfrombeginning=true;
					}
					else{ nearfrombeginning=false;}
		 		}
				if(loopcount==1 && !nearfrombeginning){ 
				// if the trajectory is far enough from the beginning and you want to record
					if(!recording){ 
						recording=true;} // start recording
	//					for(int ii=0;ii<3;ii++){
		//					cout<<ii<<' '<<x[ii].n<<' '<<endradius[ii].n<<' '<<initstate[ii].n;
			//				cout<<' '<<fabs(x[ii].n-initstate[ii].n)<<' '<< nearfrombeginning<<'\n';
			//			}
					if (time > nextrecordingtime){
						xtran.push_back(v_species()); // 
						for (auto &s: x){
							(xtran.back()).push_back(p_species(s->clone())); // xtran.push_back(x)
						}
						ttran.push_back(time);
						nextrecordingtime = nextrecordingtime + timelapse;
					}
				}
				else if(recording){// if it is close the init point and recording, stop recording
								// and erase all the trajectory recorded
					recording=false;
					xtran.clear();
					ttran.clear();
				}

				nearfromend=true; // let's assume is close to the end and check if its true
				for(uint d=0;d<x.size();d++){
					if(nearfromend && inside_end(x)){
						nearfromend=true;
					}
					else nearfromend=false;
				} 
			}// end of nearfromend while
			if(loopcount==1){trajfile<<"#  ";} // if the trajectory is recorded,
											//the times are a comment
			end=::clock();
			trajfile<<time<<' '<<(float)(end-start)/CLOCKS_PER_SEC<<' '<<maxTreached<<'\n';	
			trajfile.flush();
		}// end of recording loop (for loopcount==1)
		
		if(loopcount==1){// Print trajectory	
			for (uint i=0;i<xtran.size();i++){//				
					trajfile<<ttran[i]<<' ';
				for (auto &s : xtran[i]){
					trajfile<< s->GetNum()/Omega<<' ';
				}
				trajfile<<'\n';
			}
		}
			return time;			
	}
		
	void langil::WriteState(){
//		if (string="cout"){
//			cout<<'(';
			trajfile<<std::setprecision(12)<<time<<' ';
			for (auto &s : x){
				trajfile<< s->GetNum()/Omega<<' ';
			}
			if (fStoreDivisionTimes){
				trajfile<<cell.getPhaseTime()<<' ';
				trajfile<< divisionstate.species0before<<' ';
				trajfile<< divisionstate.species0after<<' ';
			}
			trajfile<<celldivided;
			trajfile<<'\n';
			trajfile.flush();
//			cout<<')';
//			}
	}

	void langil::WriteTempState(){ // write temporary state (usually previous time step)
//		if (string="cout"){
//			cout<<'(';
			trajfile<<std::setprecision(12)<<time-nexttau<<' ';
			for (auto &s : xtemp0){
				trajfile<< s->GetNum()/Omega<<' ';
			}
			trajfile<<celldivided;
			trajfile<<'\n';
			trajfile.flush();
//			cout<<')';
//			}
	}

	void langil::WritePrevrecState(){ // write temporary state (usually previous time step)
//		if (string="cout"){
//			cout<<'(';
			trajfile<<std::setprecision(12)<<nextrectime<<' ';
			for (auto &s : xtemp0){
				trajfile<< s->GetNum()/Omega<<' ';
			}
			trajfile<<celldivided;
			trajfile<<'\n';
			trajfile.flush();
//			cout<<')';
//			}
	}


	void langil::SetWriteState(int flag, double tl){
			fwrite=flag;
			timelapse=tl;
		}

	double langil::SetLangevinTimeStep(double xthreshold){
		// This function tries to stimate the timestep necessary to run the Langevin simulation
		// through evaluating the propensities against a certain variation in concentration xthreshold
		// It is not accurate but the measure of the optimum timestep is not a fixed quantity either
		int maxsto=0; 
		double maxprop=0;
		int stotemp=0;
		int proptemp=0;
		double timestep=0;

		if (xthreshold>0){
			dt=xthreshold;
			return xthreshold;
		}

		else{ //A negative xthreshold computres dt to get a maximum "distance" to cover per step.

		for(uint j=0;j<r.size();j++){
			for(uint i=0;i<x.size();i++){
				if ((*(r[j].GetStoichiometry()))[i] > maxsto) maxsto=(*(r[j].GetStoichiometry()))[i];
			}				
		}

		if(maxsto*maxsto/Omega/xthreshold >1){ // Case in which stochastic forces are greater
			for(uint j=0;j<r.size();j++){
				proptemp=r[j].GetPropensity();
				stotemp=0;
				for(uint i=0;i<x.size();i++){
					if ( (*(r[j].GetStoichiometry()))[i]>stotemp) stotemp=(*(r[j].GetStoichiometry()))[i] ;
				}
				if (stotemp*stotemp*proptemp > maxprop){ // look for minimum time
					maxprop=stotemp*stotemp*proptemp;
				}
			}
			timestep=Omega*Omega*xthreshold*xthreshold/maxprop;
		}

		else{   // Case in which deterministic terms are greater			
			for(vector<reaction>::iterator it = r.begin(); it != r.end(); ++it) {
				proptemp=it->GetPropensity();
				if (proptemp > maxprop){
					maxprop=proptemp;
				}
			}
			timestep=Omega*xthreshold/maxprop;
		}
		dt=timestep;
		return timestep;
		}}

void langil::Add_TimeAction(double (*actionfunc)(double)){
		actionfuncvec.push_back(actionfunc);
	}

// This function has to be rewritten carefully taking into account which is the expected 
// phase population distribution 
// void langil::setCellatRndTime(int phase){
// // The random time is set to obtain cyclostationary 	
// 	double newtime;
// 	c.cell_state = phase;
// 	c.getPhaseDuration();
// 	newtime =  

// }

	void langil::setPhaseDuration(int phase, double duration, int type_phase){
		cell.setPhaseDuration(phase,duration,type_phase);
	}

void langil::RunTimeAction(){
// if in a future several signals are desired, this should be changed by a 
// signal arrays
		// cout<<"Entering RunTimeAction at time: "<< time<<'\n';
		for (auto f: actionfuncvec){
			// cout<<"Running function inside\n";
			signal=(*f)(time);
		}
	}
void langil::Set_Boundary_Behaviour(void (*b_b)(v_species &, v_species &)){
	cout<<"Setting Boundary Behaviour\n";
	cout.flush();
	Boundary_Behaviour = b_b;
	cout.flush();
	cout<<"Boundary Behaviour set\n";	
} // Boundary behaviour
	
void langil::Set_StoreDivisionTimes(){
	fStoreDivisionTimes = true;
}

vector<StructDivisionState> langil::Get_HistoryDivision(){
	return DivisionHistoryVector;
}

void langil::Reset_HistoryDivision(){
	DivisionHistoryVector.clear();
}

void langil::SetDivisionHistory(int bef,int aft){
	divisionstate.species0before=bef;
	divisionstate.species0after=aft;
}
