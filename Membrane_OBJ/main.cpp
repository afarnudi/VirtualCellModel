#include <stdio.h>
#include "Membrane.h"
#include <ctime>
#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <limits>
#include <cstdlib>
#include <random>
#include <string>
#include <math.h>
#include "General_functions.hpp"

/*
 //this should be written in membrane class
 void Membrane_Normal_direction_Identifier( Membrane &m)
 {
 double AC[3], AB[3], ABxAC[3], xyz[3];
 int Point_A, Point_B, Point_C;
 
 for(  int i=0;i<m.Membrane_new_triangle_list.size();i++)
 {
 Point_A=m.Membrane_new_triangle_list[i][0];
 Point_B=m.Membrane_new_triangle_list[i][1];
 Point_C=m.Membrane_new_triangle_list[i][2];
 
 AB[0]=m.Membrane_Node_Position[Point_B][0]-m.Membrane_Node_Position[Point_A][0];
 AB[1]=m.Membrane_Node_Position[Point_B][1]-m.Membrane_Node_Position[Point_A][1];
 AB[2]=m.Membrane_Node_Position[Point_B][2]-m.Membrane_Node_Position[Point_A][2];
 
 AC[0]=m.Membrane_Node_Position[Point_C][0]-m.Membrane_Node_Position[Point_A][0];
 AC[1]=m.Membrane_Node_Position[Point_C][1]-m.Membrane_Node_Position[Point_A][1];
 AC[2]=m.Membrane_Node_Position[Point_C][2]-m.Membrane_Node_Position[Point_A][2];
 
 xyz[0]=m.Membrane_Node_Position[m.Membrane_new_triangle_list[i][0]][0]/2.0;
 xyz[1]=m.Membrane_Node_Position[m.Membrane_new_triangle_list[i][0]][1]/2.0;
 xyz[2]=m.Membrane_Node_Position[m.Membrane_new_triangle_list[i][0]][2];
 
 crossvector(ABxAC, AB, AC);
 //        Throughout the code the ABC vertexes of the membrane triangles are taken as A=Membrane_triangle_list[][0], B=Membrane_triangle_list[][1], C=Membrane_triangle_list[][2]. Also we often use the ABxAC cross product and we want the triangles on the membrane to point out of the cell. At the beginning of the cell construction, the centre of the cell is the same as the origin. So for the ABxAC product to point outwards, the inner product of the position of the triangle and the ABxAC should be positive. We put +/- 1 in the 'Membrane_Normal_direction[][1]' list for each triangle and define the normal direction of each triangle as Membrane_Normal_direction[i][1]*ABxAC that will always be positive, hence pointing out of the cell.
 
 if(innerproduct(ABxAC, xyz)<0 )
 {
 m.Membrane_new_triangle_list[i][1]=Point_C;
 m.Membrane_new_triangle_list[i][2]=Point_B;
 //cout<<"min"<<endl;
 }
 } // END OF: for(  int i=0;i<Membrane_num_of_Triangles;i++  )
 }// END OF: Membrane_Normal_direction_Identifier function
 
 //end of normal_direction_identifire***
 
 // i guess this should be written in membrane class
 int Membrane_triangle_pair_counter(Membrane m)
 {
 //In this function we count the total number of triangles that have a common edge (we count them twice, hence report half the number at the end).
 int temp_triangle_node_A, temp_triangle_node_B, temp_triangle_node_C;
 int triangle_pairs=0;  // This counts the number of triangle pairs that have an edge in common.
 for(int i=0 ;i<m.Membrane_new_triangle_list.size();i++)  // who are neighbors??
 {
 temp_triangle_node_A=m.Membrane_new_triangle_list[i][0];  // read the tree lable number of nodes  of every triangle
 temp_triangle_node_B=m.Membrane_new_triangle_list[i][1];
 temp_triangle_node_C=m.Membrane_new_triangle_list[i][2];
 
 for(int j=0;j<m.Membrane_new_triangle_list.size();j++)
 {
 if      ( m.Membrane_new_triangle_list[j][0]==temp_triangle_node_A  &  m.Membrane_new_triangle_list[j][1]==temp_triangle_node_B  & m.Membrane_new_triangle_list[j][2]!=temp_triangle_node_C ){
 triangle_pairs++;
 }
 if     ( m.Membrane_new_triangle_list[j][0]==temp_triangle_node_B  &  m.Membrane_new_triangle_list[j][1]==temp_triangle_node_A  & m.Membrane_new_triangle_list[j][2]!=temp_triangle_node_C ){
 triangle_pairs++;
 }
 if      ( m.Membrane_new_triangle_list[j][0]==temp_triangle_node_A  &  m.Membrane_new_triangle_list[j][1]!=temp_triangle_node_C  & m.Membrane_new_triangle_list[j][2]==temp_triangle_node_B ){
 triangle_pairs++;
 }
 if      ( m.Membrane_new_triangle_list[j][0]==temp_triangle_node_B  &  m.Membrane_new_triangle_list[j][1]!=temp_triangle_node_C  & m.Membrane_new_triangle_list[j][2]==temp_triangle_node_A ){
 triangle_pairs++;
 }
 if      ( m.Membrane_new_triangle_list[j][0]!=temp_triangle_node_C  &  m.Membrane_new_triangle_list[j][1]==temp_triangle_node_A  & m.Membrane_new_triangle_list[j][2]==temp_triangle_node_B ){
 triangle_pairs++;
 }
 if      ( m.Membrane_new_triangle_list[j][0]!=temp_triangle_node_C  &  m.Membrane_new_triangle_list[j][1]==temp_triangle_node_B  & m.Membrane_new_triangle_list[j][2]==temp_triangle_node_A ){
 triangle_pairs++;
 }
 // neibors of temp_triangle_node_B-temp_triangle_node_C :
 if      ( m.Membrane_new_triangle_list[j][0]==temp_triangle_node_B  &  m.Membrane_new_triangle_list[j][1]==temp_triangle_node_C  & m.Membrane_new_triangle_list[j][2]!=temp_triangle_node_A ){
 triangle_pairs++;
 }
 if     ( m.Membrane_new_triangle_list[j][0]==temp_triangle_node_C  &  m.Membrane_new_triangle_list[j][1]==temp_triangle_node_B  & m.Membrane_new_triangle_list[j][2]!=temp_triangle_node_A ){
 triangle_pairs++;
 }
 if      ( m.Membrane_new_triangle_list[j][0]==temp_triangle_node_B  &  m.Membrane_new_triangle_list[j][1]!=temp_triangle_node_A  & m.Membrane_new_triangle_list[j][2]==temp_triangle_node_C ){
 triangle_pairs++;
 }
 if      ( m.Membrane_new_triangle_list[j][0]==temp_triangle_node_C  &  m.Membrane_new_triangle_list[j][1]!=temp_triangle_node_A  & m.Membrane_new_triangle_list[j][2]==temp_triangle_node_B ){
 triangle_pairs++;
 }
 if      ( m.Membrane_new_triangle_list[j][0]!=temp_triangle_node_A  &  m.Membrane_new_triangle_list[j][1]==temp_triangle_node_B  & m.Membrane_new_triangle_list[j][2]==temp_triangle_node_C ){
 triangle_pairs++;
 }
 if      ( m.Membrane_new_triangle_list[j][0]!=temp_triangle_node_A  &  m.Membrane_new_triangle_list[j][1]==temp_triangle_node_C  & m.Membrane_new_triangle_list[j][2]==temp_triangle_node_B ){
 triangle_pairs++;
 }
 // neibors of temp_triangle_node_C-temp_triangle_node_A :
 if      ( m.Membrane_new_triangle_list[j][0]==temp_triangle_node_C  &  m.Membrane_new_triangle_list[j][1]==temp_triangle_node_A  & m.Membrane_new_triangle_list[j][2]!=temp_triangle_node_B ){
 triangle_pairs++;
 }
 if     ( m.Membrane_new_triangle_list[j][0]==temp_triangle_node_A  &  m.Membrane_new_triangle_list[j][1]==temp_triangle_node_C  & m.Membrane_new_triangle_list[j][2]!=temp_triangle_node_B ){
 triangle_pairs++;
 }
 if      ( m.Membrane_new_triangle_list[j][0]==temp_triangle_node_C  &  m.Membrane_new_triangle_list[j][1]!=temp_triangle_node_B  & m.Membrane_new_triangle_list[j][2]==temp_triangle_node_A ){
 triangle_pairs++;
 }
 if      ( m.Membrane_new_triangle_list[j][0]==temp_triangle_node_A  &  m.Membrane_new_triangle_list[j][1]!=temp_triangle_node_B  & m.Membrane_new_triangle_list[j][2]==temp_triangle_node_C ){
 triangle_pairs++;
 }
 if      ( m.Membrane_new_triangle_list[j][0]!=temp_triangle_node_B  &  m.Membrane_new_triangle_list[j][1]==temp_triangle_node_C  & m.Membrane_new_triangle_list[j][2]==temp_triangle_node_A ){
 triangle_pairs++;
 }
 if      ( m.Membrane_new_triangle_list[j][0]!=temp_triangle_node_B  &  m.Membrane_new_triangle_list[j][1]==temp_triangle_node_A  & m.Membrane_new_triangle_list[j][2]==temp_triangle_node_C ){
 triangle_pairs++;
 }
 }
 }
 return triangle_pairs/2;
 }
 
 //end of triangle_pair_counter****
 
 // i guess this should define in membrane class
 /*
 void Membrane_Force_Calculator (Membrane m,int Membrane_Node_Pair_list[][2],int Membrane_Triangle_Pair_Nodes[][4],double &Total_Potential_Energy, int Membrane_num_of_Triangle_Pairs, int Membrane_num_of_Node_Pairs, int Membrane_num_of_Nodes)
 {
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
 
 for (int k=0 ; k< Membrane_num_of_Node_Pairs ; k++)
 {
 temp_Node_B=Membrane_Node_Pair_list[k][0];
 temp_Node_A=Membrane_Node_Pair_list[k][1];
 
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
/* below should not be comment (Hoda)
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
 Membrane_Node_Force[temp_Node_A][0] += temp_force*deltax/temp_Node_distance+membrane_damping_coefficient*(Membrane_Node_Velocity[temp_Node_A][0]-Membrane_Node_Velocity[temp_Node_B][0]);
 Membrane_Node_Force[temp_Node_A][1] += temp_force*deltay/temp_Node_distance+membrane_damping_coefficient*(Membrane_Node_Velocity[temp_Node_A][1]-Membrane_Node_Velocity[temp_Node_B][1]);
 Membrane_Node_Force[temp_Node_A][2] += temp_force*deltaz/temp_Node_distance+membrane_damping_coefficient*(Membrane_Node_Velocity[temp_Node_A][2]-Membrane_Node_Velocity[temp_Node_B][2]);
 
 Membrane_Node_Force[temp_Node_B][0] += -temp_force*deltax/temp_Node_distance-membrane_damping_coefficient*(Membrane_Node_Velocity[temp_Node_A][0]-Membrane_Node_Velocity[temp_Node_B][0]); //from j  to i
 Membrane_Node_Force[temp_Node_B][1] += -temp_force*deltay/temp_Node_distance-membrane_damping_coefficient*(Membrane_Node_Velocity[temp_Node_A][1]-Membrane_Node_Velocity[temp_Node_B][1]);
 Membrane_Node_Force[temp_Node_B][2] += -temp_force*deltaz/temp_Node_distance-membrane_damping_coefficient*(Membrane_Node_Velocity[temp_Node_A][2]-Membrane_Node_Velocity[temp_Node_B][2]);
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
/* below should not be comment (Hoda)
 
 Total_Potential_Energy += Membrane_bending_coefficient*(1.0 -  ( innerproduct(N1,N2)/(vectorlength(N1)*vectorlength(N2)  )   ) );
 }  // end of 'for-i'
 //    exit (EXIT_FAILURE);
 }///end of function
 
 
 void ConstantSurfaceForceLocalTriangles(double Membrane_Node_Position[][3],double Membrane_Node_Force[][3], vector<vector<int> > Membrane_new_triangle_list)
 {
 
 int temp_node_A, temp_node_B, temp_node_C;
 double AB[3], AC[3], temp_AB_length_squared, temp_AC_length_squared, f0, temp_ABxAC[3];
 
 // cout <<surfacearea(x,tri )<<"   s0="<< s0<<endl;
 double s0_i=0.41*Node_radius*Node_radius; //1.732=3^0.5
 double s_i=0;
 
 
 for(  int i=0;i<Membrane_new_triangle_list.size();i++)
 {
 temp_node_A=Membrane_new_triangle_list[i][0];
 temp_node_B=Membrane_new_triangle_list[i][1];
 temp_node_C=Membrane_new_triangle_list[i][2];
 
 AB[0]=Membrane_Node_Position[temp_node_B][0]-Membrane_Node_Position[temp_node_A][0];
 AB[1]=Membrane_Node_Position[temp_node_B][1]-Membrane_Node_Position[temp_node_A][1];
 AB[2]=Membrane_Node_Position[temp_node_B][2]-Membrane_Node_Position[temp_node_A][2];
 
 AC[0]=Membrane_Node_Position[temp_node_C][0]-Membrane_Node_Position[temp_node_A][0];
 AC[1]=Membrane_Node_Position[temp_node_C][1]-Membrane_Node_Position[temp_node_A][1];
 AC[2]=Membrane_Node_Position[temp_node_C][2]-Membrane_Node_Position[temp_node_A][2];
 
 crossvector(temp_ABxAC, AB, AC);
 s_i = vectorlength(temp_ABxAC)/2.0;
 f0= K_surfaceConstant_local*(s_i -  s0_i )/2.0*vectorlength(temp_ABxAC);
 
 for (int j=0; j<3; j++) {
 temp_node_A=Membrane_new_triangle_list[i][j%3];
 temp_node_B=Membrane_new_triangle_list[i][(j+1)%3];
 temp_node_C=Membrane_new_triangle_list[i][(j+2)%3];
 
 AB[0]=Membrane_Node_Position[temp_node_B][0]-Membrane_Node_Position[temp_node_A][0];
 AB[1]=Membrane_Node_Position[temp_node_B][1]-Membrane_Node_Position[temp_node_A][1];
 AB[2]=Membrane_Node_Position[temp_node_B][2]-Membrane_Node_Position[temp_node_A][2];
 
 AC[0]=Membrane_Node_Position[temp_node_C][0]-Membrane_Node_Position[temp_node_A][0];
 AC[1]=Membrane_Node_Position[temp_node_C][1]-Membrane_Node_Position[temp_node_A][1];
 AC[2]=Membrane_Node_Position[temp_node_C][2]-Membrane_Node_Position[temp_node_A][2];
 
 temp_AB_length_squared=vectorlength(AB)*vectorlength(AB);
 temp_AC_length_squared=vectorlength(AC)*vectorlength(AC);
 
 Membrane_Node_Force[temp_node_A][0] +=  f0 *2.0* ( -temp_AC_length_squared *  AB[0]  -temp_AB_length_squared *  AC[0] + innerproduct(AB,AC)*  ( AB[0]+AC[0] )   ) ;
 Membrane_Node_Force[temp_node_A][1] +=  f0 *2.0* ( -temp_AC_length_squared *  AB[1]  -temp_AB_length_squared *  AC[1] + innerproduct(AB,AC)*  ( AB[1]+AC[1] )   ) ;
 Membrane_Node_Force[temp_node_A][2] +=  f0 *2.0* ( -temp_AC_length_squared *  AB[2]  -temp_AB_length_squared *  AC[2] + innerproduct(AB,AC)*  ( AB[2]+AC[2] )   ) ;
 }
 }
 
 }
 
 */

void trajectory (Membrane membrane, string label)
{
    ofstream Trajectory;
    Trajectory.open(membrane.output_file_neme.c_str());
    Trajectory << std:: fixed;
    Trajectory << membrane.Membrane_num_of_Nodes<<endl;
    Trajectory << " nodes  "<<endl;
    for(int j=0; j< membrane.Membrane_num_of_Nodes;j++) // saving trajectory
    {
        Trajectory << label <<setprecision(5)<< setw(20)<<membrane.Membrane_Node_Position[j][0]<< setw(20)<<membrane.Membrane_Node_Position[j][1]<< setw(20)<<membrane.Membrane_Node_Position[j][2]<<endl;
    }
}


int main(int argc, char **argv)
{
    /*clock_t tStart = clock();//Time the programme
     time_t t = time(0);   // get time now
     struct tm * now = localtime( & t );
     char buffer [80];
     strftime (buffer,80,"%Y_%m_%d_%H:%M",now);
     //outputfiles:
     /*string traj_file_name;
     
     traj_file_name="results/membrane_";
     traj_file_name +=buffer;
     traj_file_name +=".xyz";*/
    
    //just adding a comment
    Membrane  membrane("membrane_2D_mesh");
    //    Membrane cage("Cage");
//    for (int i=0;i<200;i++)
//    {
        trajectory(membrane, "membrane");
        //    trajectory(cage, "cage");
//    }
    
    
    
    return 0;
}

