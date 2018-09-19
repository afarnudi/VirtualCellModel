#include "Membrane.h"

void Membrane::Elastic_Force_Calculator(double theta_0)
{
	//modifications on this function:
	// 1- changing Membrane_Node_pair_list to Membrane_Edges
	// 2- changing Membrane_Num_of_Node_Pairs to Membrane_num_of_Triangle_Pairs (because these 2 numbers are equal)
	
//    double le0,le1,lmax,lmin;
//    double deltax,deltay,deltaz,temp_Node_distance,temp_force;
//    int pos1,pos2,pos3,pos4;  // to making calculation of surface force easier
////    double temp_potential_energy = 0.0;
//    double temp_p1[3], temp_p2[3], temp_p3[3], temp_p4[3];
//    double  N1[3], N2[3], N3[3], p3p1[3], p3p2[3], p4p2[3], p4p1[3], Ep2p1[3], sinus, F0, F1[3], F2[3], F3[3], F4[3];// for exmple p3p1 is p3-p1 and ....

    /// calculate network force:
//    int temp_Node_A, temp_Node_B;
//    le0=1.15000*Node_radius;
//    lmax=1.33000*Node_radius;
//    le1=0.85000*Node_radius;
//    lmin=0.67000*Node_radius;
    Total_Potential_Energy=0.0;
	if (spring_model=1) {potential_1();}
    if (spring_model=2) {potential_2();}

//    Bending_potetial();
    Bending_potetial_2(0);
    //    exit (EXIT_FAILURE);
}
