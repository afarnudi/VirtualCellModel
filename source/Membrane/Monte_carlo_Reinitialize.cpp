#include "Membrane.h"
#include "OpenMM_structs.h"
#include "OpenMM_funcs.hpp"
#include "General_functions.hpp"

void Monte_Carlo_Reinitialize(MyOpenMMData* omm, Bonds* bonds, Dihedrals* dihedrals, Membrane &mem, MyAtomInfo atoms[], int &MC_total_tries, int&Accepted_Try_Counter,  double &MC_Acceptance_Rate){
    bool preservestate=1;
//     double time, initenergy, finalenergy, initpenergy, finalpenergy;
     double localDeltaE=0;
     int pyramid_counter=0;
     for(int i=0; i<GenConst::MC_step; i++){


       
        mem.monte_carlo_flip(omm, bonds, dihedrals, atoms,localDeltaE, Accepted_Try_Counter, pyramid_counter,MC_total_tries, MC_Acceptance_Rate);
       

     }
     omm->context->reinitialize(preservestate);

     //cout<<"num_of_accepted tries  "<<Accepted_Try_Counter<<"  out of  "<<GenConst::MC_step<<"  pyramid_counter  "<<pyramid_counter<<endl;
     
}

