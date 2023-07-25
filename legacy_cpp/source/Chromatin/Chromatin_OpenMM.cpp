#include <sstream>
#include "Chromatin.h"

void Chromatin::set_state(MyAtomInfo all_atoms[], int atom_count){
    for (int i=0; i<Num_of_Nodes; i++) {
        for (int j=0; j<3; j++) {
            Node_Position[i][j] = all_atoms[i+atom_count].posInNm[j];
            Node_Velocity[i][j] = all_atoms[i+atom_count].velocityInNmperPs[j];
        }
    }
}
