#include "Membrane.h"
#include "OpenMM_structs.h"
#include "OpenMM_funcs.hpp"
#include "General_functions.hpp"

double Membrane::calculating_the_bend_energy(int uncommon1, int common2, int common3, int uncommon4, bool initial_or_final,MyAtomInfo  atoms[], int number_of_privious_mem_nodes){
    double p1[3],p2[3],p3[3],p4[3];
    double p3p4[3],p2p4[3],p3p1[3],p2p1[3];
    double outward_vector[3];
    double N1[3], N2[3];
    double N1_length,N2_length;
    double cosine=1;
    double bending_energy;
    double origin_point[3]={0,0,0}; //in order to use the center of mass it is not working if we use membrane's COM, because node positions in memberane class dont update
    for (int index=0; index<3; index++){
        p1[index]=atoms[uncommon1+number_of_privious_mem_nodes].posInAng[index];
        p2[index]=atoms[common2+number_of_privious_mem_nodes].posInAng[index];
        p3[index]=atoms[common3+number_of_privious_mem_nodes].posInAng[index];
        p4[index]=atoms[uncommon4+number_of_privious_mem_nodes].posInAng[index];
        
        p3p4[index]=p4[index]-p3[index];
        p2p4[index]=p4[index]-p2[index];
        p3p1[index]=p1[index]-p3[index];
        p2p1[index]=p1[index]-p2[index];
        
        
        outward_vector[index]=p1[index]- origin_point[index];
    }
    crossvector(N1, p3p1, p2p1);
    crossvector(N2, p3p4, p2p4);
    if (innerproduct(N1,outward_vector)<0){
        N1[0]=-N1[0];
        N1[1]=-N1[1];
        N1[2]=-N1[2];
    }
    if (innerproduct(N2,outward_vector)<0){
        N2[0]=-N2[0];
        N2[1]=-N2[1];
        N2[2]=-N2[2];
    }
    N1_length=vector_length(N1);
    N2_length=vector_length(N2);
    if(N1_length !=0 and N2_length !=0){
        cosine= innerproduct(N1, N2)/(vector_length(N1)*vector_length(N2));}
    else if(initial_or_final==0){
        cout<<"warning! The trianglur mesh messed up. its impossible to calculate the bending energy." <<endl;
        EXIT_FAILURE;
    }
    else{
        bending_energy=-1;
        cout<<"monte_carlo flip rejected because of bad configuration"<<endl;
    }
    bending_energy= Bending_coefficient //* OpenMM::KJPerKcal
                                        *(1.00001-cosine);
    return(bending_energy);
}

double Membrane::calculating_the_bend_energy_2(int uncommon1, int common2, int common3, int uncommon4,MyAtomInfo atoms[],int number_of_privious_mem_nodes){
    double bending_energy=0;
    int node_A, node_B, node_C, node_D;
    
    double points[3][3];
    
    double A, B, C, E, F, G;
    node_C=uncommon1; 
    node_A=common2;
    node_B=common3;
    node_D=uncommon4;
    
    //Update the triangle node positions from OpenMM
    for (int k=0; k<3; k++) {
            points[0][k]=atoms[node_A+number_of_privious_mem_nodes].posInAng[k];
            points[1][k]=atoms[node_B+number_of_privious_mem_nodes].posInAng[k];
            points[2][k]=atoms[node_C+number_of_privious_mem_nodes].posInAng[k];
           
    }
    calc_surface_coefficeints_2(points, A, B, C);
            
    for (int k=0; k<3; k++) {
            points[2][k]=atoms[node_D+number_of_privious_mem_nodes].posInAng[k];
            
    }
    calc_surface_coefficeints_2(points, E, F, G);
   
    double cosine=-1;
    double denominator = sqrt(A*A + B*B + C*C) * sqrt(E*E + F*F + G*G);
    if (denominator > 0.001) {
            
            cosine = ( A*E + B*F + C*G )/denominator;
           
    }
            
    bending_energy= Bending_coefficient //* OpenMM::KJPerKcal
                                        *(1+ cosine);
   
    return(bending_energy);
    
}    
