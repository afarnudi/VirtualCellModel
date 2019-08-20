#include <sstream>
#include "Membrane.h"

void Membrane::set_state(MyAtomInfo all_atoms[], int atom_count){
    for (int i=0; i<Num_of_Nodes; i++) {
        for (int j=0; j<3; j++) {
            Node_Position[i][j] = all_atoms[i+atom_count].posInAng[j];
            Node_Velocity[i][j] = all_atoms[i+atom_count].velocityInAngperPs[j];
        }
    }
}
