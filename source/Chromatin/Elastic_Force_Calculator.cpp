#include "Chromatin.h"

void Chromatin::Force_Calculator(void)
{
    Total_Potential_Energy=0.0;
	if (spring_model==0) {potential_1();}
    if (spring_model==1) {FENE();}
}

// This is the old code and it is designed to work with a chromatin chain with only 2 node types.
void Chromatin::Force_Calculator_2(void)
{
    Total_Potential_Energy=0.0;
    if (spring_model==0) {
        potential_1();
    } else if (spring_model==1) {
        FENE();
    }
    
    double deltax,deltay,deltaz,Node_distance,temp_force;
    double temp_potential_energy = 0.0;//, epsilon;
    int Node_A, Node_B;
//    double le1=2.3*Node_radius;
//    double lmin=2.*Node_radius;
    
//    double Threshold_force=10;
    
    for (int k=0 ; k< Num_of_Nodes-2 ; k++)
    {
        Node_B=k;
        for (int j=k+2; j<Num_of_Nodes; j++) {
            Node_A=j;
            
            deltax=Node_Position[Node_A][0]-Node_Position[Node_B][0];
            deltay=Node_Position[Node_A][1]-Node_Position[Node_B][1];
            deltaz=Node_Position[Node_A][2]-Node_Position[Node_B][2];
            
            Node_distance=sqrt(deltax*deltax+deltay*deltay+deltaz*deltaz);
            Node_distance /= Node_radius;
            temp_force=0.0;
            
//            if (Node_distance < 4.0 && Node_distance > 2.0) {
//                if (Node_distance<2.5) {
//                    Contact_Matrix[k][j]++;
//                }
//                if (ABC_index[k]*ABC_index[j] != 0) {
//                    double sigma = 0.35;
//                    double r_1 = (Node_distance -2.0)/sigma;
//                    double r_3 = r_1*r_1*r_1;
//                    double r_5 = r_3*r_1*r_1;
//                    epsilon = 4 * GenConst::K * GenConst::MD_T * 0.01;
//
//                    temp_force = 2*epsilon*( 2.0/r_5 - 1.0/r_3 )/sigma;
//                    temp_potential_energy = epsilon * ( r_1*1.0/r_5 - r_1*1.0/r_3  );
////
////
//                }
//                    else if (Node_distance < 2.3){
//                        double exp_le1=exp(1.0/(Node_distance-le1));
//                        temp_force = ( (Spring_coefficient*exp_le1)/(Node_distance-lmin) )*( 1/(Node_distance-lmin)+1/( (Node_distance-le1)*(Node_distance-le1) ) );
//                        temp_potential_energy = Spring_coefficient*exp_le1/(Node_distance-lmin);
//                }
//
//            }
//            else if(temp_force<-Threshold_force   ||  Node_distance < 2 )
//            {
////                if (Node_distance<1) {
////                    cout<<"warning!\n";
////                }
//
//                temp_force = Threshold_force*Node_distance/2.0 - 2*Threshold_force;
//                temp_potential_energy=   Threshold_force*Node_distance*Node_distance/4.0 -2*Threshold_force*Node_distance;
//
////                double c=1000;
////                temp_force = -2.0*c*Node_distance/(Node_radius) + 3*c;
////                temp_potential_energy=   c*Node_distance*Node_distance/(Node_radius) - 3*c*Node_distance;
//            }
            
            Total_Potential_Energy += temp_potential_energy;
            
            temp_force=temp_force/Node_distance;
            
            // implimentation of forces:
            Node_Force[Node_A][0] +=  temp_force*deltax;
            Node_Force[Node_A][1] +=  temp_force*deltay;
            Node_Force[Node_A][2] +=  temp_force*deltaz;
            
            Node_Force[Node_B][0] += -temp_force*deltax;
            Node_Force[Node_B][1] += -temp_force*deltay;
            Node_Force[Node_B][2] += -temp_force*deltaz;
        }
    }
}
