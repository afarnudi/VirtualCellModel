#include "Membrane.h"
#include "OpenMM_structs.h"
#include "OpenMM_funcs.hpp"
#include "General_functions.hpp"

void Monte_Carlo_Reinitialize(MyOpenMMData* omm, Bonds* bonds, Dihedrals* dihedrals, Membrane &mem, MyAtomInfo atoms[], int &MC_total_tries, int&Accepted_Try_Counter,  double &MC_Acceptance_Rate){
    bool preservestate=1;
     double time, initenergy,finalenergy,initpenergy, finalpenergy;
     double localDeltaE=0;
     int pyramid_counter=0;
     for(int i=0; i<GenConst::MC_step; i++){

        //myGetOpenMMState(omm, time, initenergy,initpenergy, atoms);
        mem.monte_carlo_flip(omm, bonds, dihedrals, atoms,localDeltaE, Accepted_Try_Counter, pyramid_counter,MC_total_tries, MC_Acceptance_Rate);
        /*
        omm->context->reinitialize(preservestate);
        
        myGetOpenMMState(omm, time, finalenergy,finalpenergy, atoms);
        double globalDeltaE=finalpenergy-initpenergy;
        if(globalDeltaE!=0){
            if(abs(localDeltaE-globalDeltaE)>0.0000001){
                cout<<"warning! local and global DeltaE are different  "<<"global Delta potential energy  "<<globalDeltaE<<"   "<<"locall  "<<localDeltaE<<endl;
                cout<<"Test  "<<globalDeltaE-localDeltaE<<endl;
            }

        }*/
    
     }
     omm->context->reinitialize(preservestate);

     //cout<<"num_of_accepted tries  "<<Accepted_Try_Counter<<"  out of  "<<GenConst::MC_step<<"  pyramid_counter  "<<pyramid_counter<<endl;
     
}