#include "Membrane.h"

void Membrane::Elastic_Force_Calculator()
{
	//modifications on this function:
	// 1- changing Membrane_Node_pair_list to Membrane_Edges
	// 2- changing Membrane_Num_of_Node_Pairs to Membrane_num_of_Triangle_Pairs (because these 2 numbers are equal)
	
    double le0,le1,lmax,lmin;
    double deltax,deltay,deltaz,temp_Node_distance,temp_force;
    int pos1,pos2,pos3,pos4;  // to making calculation of surface force easier
    double temp_potential_energy = 0.0;
    double temp_p1[3], temp_p2[3], temp_p3[3], temp_p4[3];
    double  N1[3], N2[3], N3[3], p3p1[3], p3p2[3], p4p2[3], p4p1[3], Ep2p1[3], sinus, F0, F1[3], F2[3], F3[3], F4[3];// for exmple p3p1 is p3-p1 and ....

    /// calculate network force:
    int temp_Node_A, temp_Node_B;
    le0=1.15000*Node_radius;
    lmax=1.33000*Node_radius;
    le1=0.85000*Node_radius;
    lmin=0.67000*Node_radius;
    Total_Potential_Energy=0.0;

    for (int k=0 ; k< Membrane_num_of_Triangle_Pairs ; k++)
    {
        temp_Node_B=Membrane_Edges[k][0];
        temp_Node_A=Membrane_Edges[k][1];

        deltax=Membrane_Node_Position[temp_Node_A][0]-Membrane_Node_Position[temp_Node_B][0];
        deltay=Membrane_Node_Position[temp_Node_A][1]-Membrane_Node_Position[temp_Node_B][1];
        deltaz=Membrane_Node_Position[temp_Node_A][2]-Membrane_Node_Position[temp_Node_B][2];


        temp_Node_distance=sqrt(deltax*deltax+deltay*deltay+deltaz*deltaz);
        temp_force=0.0;
        double temp_exp_le0=exp(1.0/(le0-temp_Node_distance));
        double temp_exp_le1=exp(1.0/(temp_Node_distance-le1));
        //*******************************************************************************************************
        /*BUG
         |---\   |    |  /---\
         |    |  |    |  |
         |---<   |    |  |  -\
         |    |  |    |  |   |
         |---/   \----/  \---/
         */
        //*******************************************************************************************************
        //***************** Potential BUG: F=-dU/dr but in many cases I cannot determin wheather ****************
        //***************** the '-' has been implemented or not. Since the potential energy is   ****************
        //***************** never used in the code it does not a threat. ****************************************
        //*******************************************************************************************************

        if(temp_Node_distance >le1  & temp_Node_distance < le0 )  //free zone
        {
            temp_potential_energy=0 ; // free zone
        }

        if(temp_Node_distance > le0  & temp_Node_distance <lmax )  //bondforce
        {
            temp_force = (Membrane_spring_coefficient*temp_exp_le0/(lmax-temp_Node_distance))*( 1.0/(lmax-temp_Node_distance) +  1.0/((le0-temp_Node_distance)*(le0-temp_Node_distance)));
            temp_potential_energy= Membrane_spring_coefficient*temp_exp_le0/(lmax-temp_Node_distance);

        }

        if(temp_Node_distance < le1   &  temp_Node_distance > lmin  )  // repulsive force
        {
            temp_force= -(Membrane_spring_coefficient*temp_exp_le1/(temp_Node_distance-lmin))*( 1.0/(temp_Node_distance-lmin) + 1.0/((temp_Node_distance-le1)*(temp_Node_distance-le1)));                 // force on i th from j
            temp_potential_energy= Membrane_spring_coefficient*temp_exp_le1/(temp_Node_distance-lmin);
        }
        /// my cutoff for force amplitute and for avoiding leting particle scape from force trap
        if(temp_force>965.31  || temp_Node_distance>lmax )
        {
            temp_force = 965.31+Membrane_spring_force_cutt_off* ( temp_Node_distance - 1.3280*Node_radius );
            temp_potential_energy=   1.81599  + 965.31 * ( temp_Node_distance - 1.3280*Node_radius )+0.5*Membrane_spring_force_cutt_off * ( temp_Node_distance - 1.3280*Node_radius ) * ( temp_Node_distance - 1.3280*Node_radius );
        }


        if(temp_force<-1000.05   ||  temp_Node_distance<lmin )
        {
            temp_force =-1000.05-Membrane_spring_force_cutt_off* ( 0.671965*Node_radius - temp_Node_distance );
            temp_potential_energy = 1.85038 + 1005.05 * ( 0.671965*Node_radius - temp_Node_distance )+0.5*Membrane_spring_force_cutt_off*( 0.671965*Node_radius - temp_Node_distance )*( 0.671965*Node_radius - temp_Node_distance );
        }

        Total_Potential_Energy += temp_potential_energy;

        // implimentation of forces:
        Membrane_Node_Force[temp_Node_A][0] += temp_force*deltax/temp_Node_distance+Membrane_damping_coefficient*(Membrane_Node_Velocity[temp_Node_A][0]-Membrane_Node_Velocity[temp_Node_B][0]);
        Membrane_Node_Force[temp_Node_A][1] += temp_force*deltay/temp_Node_distance+Membrane_damping_coefficient*(Membrane_Node_Velocity[temp_Node_A][1]-Membrane_Node_Velocity[temp_Node_B][1]);
        Membrane_Node_Force[temp_Node_A][2] += temp_force*deltaz/temp_Node_distance+Membrane_damping_coefficient*(Membrane_Node_Velocity[temp_Node_A][2]-Membrane_Node_Velocity[temp_Node_B][2]);

        Membrane_Node_Force[temp_Node_B][0] += -temp_force*deltax/temp_Node_distance-Membrane_damping_coefficient*(Membrane_Node_Velocity[temp_Node_A][0]-Membrane_Node_Velocity[temp_Node_B][0]); //from j  to i
        Membrane_Node_Force[temp_Node_B][1] += -temp_force*deltay/temp_Node_distance-Membrane_damping_coefficient*(Membrane_Node_Velocity[temp_Node_A][1]-Membrane_Node_Velocity[temp_Node_B][1]);
        Membrane_Node_Force[temp_Node_B][2] += -temp_force*deltaz/temp_Node_distance-Membrane_damping_coefficient*(Membrane_Node_Velocity[temp_Node_A][2]-Membrane_Node_Velocity[temp_Node_B][2]);
    }
    // End of Membrane Node Pair forces

    // Beginning of the  triangle-triangle (bending) force calculations
    for(int i=0 ;i<Membrane_num_of_Triangle_Pairs;i++)  // who are neighbors?
    {
        pos1=Membrane_Triangle_Pair_Nodes[i][0];
        pos2=Membrane_Triangle_Pair_Nodes[i][1];
        pos3=Membrane_Triangle_Pair_Nodes[i][2];
        pos4=Membrane_Triangle_Pair_Nodes[i][3];
        for (int index=0; index<3; index++) {
            temp_p1[index]=Membrane_Node_Position[pos1][index];
            temp_p2[index]=Membrane_Node_Position[pos2][index];
            temp_p3[index]=Membrane_Node_Position[pos3][index];
            temp_p4[index]=Membrane_Node_Position[pos4][index];
            p3p1[index]=temp_p3[index]-temp_p1[index];
            p3p2[index]=temp_p3[index]-temp_p2[index];
            p4p2[index]=temp_p4[index]-temp_p2[index];
            p4p1[index]=temp_p4[index]-temp_p1[index];
            Ep2p1[index]=temp_p2[index]-temp_p1[index];
        }
        crossvector(N1, p3p1, p3p2);
        crossvector(N2, p4p2, p4p1);
        crossvector(N3, N2, N1);
        sinus=vectorlength(N3)/(vectorlength(N2)*vectorlength(N1));
        crossvector(N3, N1, N2);
        double temp_Ep2p1_length=vectorlength(Ep2p1);
        //        if ((sinus- temp_sinus)>0.0001) {
        //            cout<< "oops!"<<endl;
        //            cout<<"temp_sinus="<<temp_sinus<<"\nsinus="<<sinus<<endl;
        //            exit (EXIT_FAILURE);
        //        }
        F0 = -(1-2*signbit(innerproduct(N3, Ep2p1)))*Membrane_bending_coefficient*sinus;
        //        cout<<"\nF0="<<F0<<endl;
        //        if( parallelORantiparallel(xpos ) == +1 )
        //        {
        //            F0=-F0;
        //        }
        double temp_N1_length_squared=vectorlength(N1)*vectorlength(N1);
        double temp_N2_length_squared=vectorlength(N2)*vectorlength(N2);
        //**************************************************** force calculation
        for (int l=0; l<3; l++) {
            F3[l]= F0 * temp_Ep2p1_length* N1[l]/ temp_N1_length_squared;

            F4[l]= F0 * temp_Ep2p1_length* N2[l]/ temp_N2_length_squared;

            F1[l]= (F0/temp_Ep2p1_length)*( innerproduct(p3p2,Ep2p1)*N1[l]/temp_N1_length_squared + innerproduct(p4p2,Ep2p1)*N2[l]/temp_N2_length_squared );

            F2[l]= (-F0/temp_Ep2p1_length)*( innerproduct(p3p1,Ep2p1)*N1[l]/temp_N1_length_squared + innerproduct(p4p1,Ep2p1)*N2[l]/temp_N2_length_squared );
            Membrane_Node_Force[pos1][l] += F1[l];
            Membrane_Node_Force[pos2][l] += F2[l];
            Membrane_Node_Force[pos3][l] += F3[l];
            Membrane_Node_Force[pos4][l] += F4[l];

        }

        //*******************************************************************************************************
        /*BUG
         |---\   |    |  /---\
         |    |  |    |  |
         |---<   |    |  |  -\
         |    |  |    |  |   |
         |---/   \----/  \---/
         */
        //*******************************************************************************************************
        //***************** Potential BUG: Not very sure about the calculations *****************************
        //***************** for 'Total_Potential_Energy. Need to double check *********************************
        //*******************************************************************************************************
        Total_Potential_Energy += Membrane_bending_coefficient*(1.0 -  ( innerproduct(N1,N2)/(vectorlength(N1)*vectorlength(N2)  )   ) );
    }  // end of 'for-i'
    //    exit (EXIT_FAILURE);
}
