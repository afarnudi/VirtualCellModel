#include "Membrane.h"

void Membrane::Rescale(double rescale_factor){
    for(int j=0 ; j<Num_of_Nodes ; j++){
        for(int i=0; i<3; i++){
            Node_Position[j][i]/=rescale_factor;}
            }
}