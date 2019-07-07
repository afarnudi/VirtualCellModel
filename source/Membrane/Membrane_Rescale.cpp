#include "Membrane.h"

//void Membrane::Rescale(double rescale_factor){
//    for(int j=0 ; j<Num_of_Nodes ; j++){
//        for(int i=0; i<3; i++){
//            Node_Position[j][i]*=rescale_factor;}
//            }
//}

void Membrane::shift_position (double x, double y, double z) {
      for (int i=0; i<Num_of_Nodes; i++) {
            Node_Position[i][0]+=x;
            Node_Position[i][1]+=y;
            Node_Position[i][2]+=z;
        }
        X_in+=x;
        Y_in+=y;
        Z_in+=z;
}
