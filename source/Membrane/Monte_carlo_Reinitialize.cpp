#include "Membrane.h"
#include "OpenMM_structs.h"
#include "OpenMM_funcs.hpp"
#include "General_functions.hpp"

void Monte_Carlo_Reinitialize(MyOpenMMData* omm, Bonds* bonds, Dihedrals* dihedrals, Membrane &mem,   MyAtomInfo atoms[]){
    bool preservestate=1;
    double time, initenergy,initpenergy;//,finalenergy, finalpenergy;
    double localDeltaE=0;
    int Accepted_Try_Counter=0;
    int pyramid_counter=0;
    for(int i=0; i<GenConst::MC_step; i++){
        
        myGetOpenMMState(omm, time, initenergy,initpenergy, atoms);
        if ( mem.monte_carlo_flip(omm, bonds, dihedrals, atoms,localDeltaE, Accepted_Try_Counter, pyramid_counter) ){
            omm->context->reinitialize(preservestate);
        }
        
        
        //        myGetOpenMMState(omm, time, finalenergy,finalpenergy, atoms);
        
    }
    cout<<"Accepted tries : "<<Accepted_Try_Counter<<" [out of "<<GenConst::MC_step<<" ]"<<endl;
}
