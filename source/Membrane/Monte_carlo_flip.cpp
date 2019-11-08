#include "Membrane.h"
#include "OpenMM_structs.h"
#include "OpenMM_funcs.hpp"
#include "General_functions.hpp"

void Membrane::monte_carlo_flip(MyOpenMMData* omm, Bonds* bonds, Dihedrals* dihedrals, MyAtomInfo atoms[], double& localDeltaE, int& Accepted_Try_Counter, int& pyramid_counter){
   
    bool accept=0; 
    
    //Randomly choosing a Dihedral and finding the neighbours.

    int triangle_A, triangle_B,temp_triangle, u1;
    int initial_pair;
    vector<int> A_neighbours_dihedral_index;
    vector<int> B_neighbours_dihedral_index;
    double initial_bend_energy=0;
    double initial_bond_energy=0;
    double final_bend_energy=0;
    double final_bond_energy=0;
    int number_of_privious_mem_nodes=0;
    
    int new_neighbour_dihedrals[4][6];
    
    //** the new pair must have 4 new neighbours.
    //**each array contains the information of new neighbours. the first 4 elements are the index of nodes which constructs the dihedral. the fifths element indicates the index of the dihedral in both Triangle_Pair_Nodes and Triangle_pair_list. the last element is a boolian which 0 indicates that we dont have to update the indexes of triangles which are connected to each other and one indicates that we have to update the lable of new connected triangles.
    
    initial_pair=rand()%(Num_of_Triangle_Pairs-1);
    
    //initial_pair=2;
    triangle_A=Triangle_pair_list[initial_pair][0];
    triangle_B=Triangle_pair_list[initial_pair][1];
    //the next part is for making sure that the first uncommon node in dihedral belongs to triangle A;
    u1=Triangle_Pair_Nodes[initial_pair][0];
    
    if(Triangle_list[triangle_A][0]!=u1 and Triangle_list[triangle_A][1]!=u1 and Triangle_list[triangle_A][2]!=u1){
        temp_triangle=triangle_A;
        triangle_A=triangle_B;
        triangle_B=temp_triangle;
        
    }//end of if(Triangle_list[triangle_A][0]!=u1 and Triangle_list[triangle_A][1]!=u1 and Triangle_list[triangle_A][2]!=u1)
 
    for (int i=0; i<Num_of_Triangle_Pairs ; i++){
        
        if ((Triangle_pair_list[i][0]==triangle_A and Triangle_pair_list[i][1]!=triangle_B)  || (Triangle_pair_list[i][1]==triangle_A and Triangle_pair_list[i][0]!=triangle_B)){

            A_neighbours_dihedral_index.push_back(i);}
        if ((Triangle_pair_list[i][0]==triangle_B and Triangle_pair_list[i][1]!=triangle_A)  || (Triangle_pair_list[i][1]==triangle_B and Triangle_pair_list[i][0]!=triangle_A)){
        
            B_neighbours_dihedral_index.push_back(i);}
        
    }//end of for (int i=0; i<Num_of_Triangle_Pairs ; i++)
    
    
    if (A_neighbours_dihedral_index.size()==2 and B_neighbours_dihedral_index.size()==2){
        if(check_Pyramid_2(initial_pair,A_neighbours_dihedral_index, B_neighbours_dihedral_index)==false and check_Pyramid(A_neighbours_dihedral_index, B_neighbours_dihedral_index)==false){
            accept=1;
        }
        else{
            pyramid_counter+=1;
        }
    }//end of  if (A_neighbours_dihedral_index.size()==2 and B_neighbours_dihedral_index.size()==2)
        
 
   else{
       cout<<"no enough neighbours,    num of naighbours  "<<A_neighbours_dihedral_index.size()<<"  and  "<<B_neighbours_dihedral_index.size()<<endl;
        cout<<"index  "<<initial_pair<<endl;
    }//end of else 
    
    
    
    if (accept){ //(checking the conditions)
        //cout<<"starting monte_carlo try"<<endl;
        
        int uncommon1,common2,common3, uncommon4;
        //calculating the initial state energy
        ////calculating the bond energy
        bool initial_or_final=0;
        initial_bond_energy=calculating_the_bond_energy(initial_pair, initial_or_final, atoms, number_of_privious_mem_nodes);
        //cout<<"initial_bond_energy  "<<initial_bond_energy<<endl;
   
        ////calculating the bend energy of initial pairs
        //cout<<"calculating the bend energy of initial pairs"<<endl;
        uncommon1=Triangle_Pair_Nodes[initial_pair][0];
        common2=Triangle_Pair_Nodes[initial_pair][1];
        common3=Triangle_Pair_Nodes[initial_pair][2];
        uncommon4=Triangle_Pair_Nodes[initial_pair][3];
        initial_bend_energy+=calculating_the_bend_energy_2(uncommon1, common2, common3, uncommon4, atoms, number_of_privious_mem_nodes);
    
        ////calculating the bend energy of initial neighbours.
        for(int i=0; i<2; i++){
            //triangle_A neighbours
            uncommon1=Triangle_Pair_Nodes[A_neighbours_dihedral_index[i]][0];
            common2=Triangle_Pair_Nodes[A_neighbours_dihedral_index[i]][1];
            common3=Triangle_Pair_Nodes[A_neighbours_dihedral_index[i]][2];
            uncommon4=Triangle_Pair_Nodes[A_neighbours_dihedral_index[i]][3];
          //  cout<<"previous neigbours   "<<uncommon1<<"  "<< common2<<" "<<common3<<"  "<<uncommon4<<"  "<<A_neighbours_dihedral_index[i]<<endl;
            initial_bend_energy+=calculating_the_bend_energy_2(uncommon1, common2, common3, uncommon4,atoms, number_of_privious_mem_nodes);
            //triangle_B neighbours
            uncommon1=Triangle_Pair_Nodes[B_neighbours_dihedral_index[i]][0];
            common2=Triangle_Pair_Nodes[B_neighbours_dihedral_index[i]][1];
            common3=Triangle_Pair_Nodes[B_neighbours_dihedral_index[i]][2];
            uncommon4=Triangle_Pair_Nodes[B_neighbours_dihedral_index[i]][3];
           // cout<<"previous neigbours   "<<uncommon1<<"  "<< common2<<" "<<common3<<"  "<<uncommon4<<"  "<<B_neighbours_dihedral_index[i]<<endl;
            initial_bend_energy+=calculating_the_bend_energy_2(uncommon1, common2, common3, uncommon4, atoms, number_of_privious_mem_nodes);
        } //end of for(int i=0; i<2; i++) (calculating the bending energy of initial neighbors)
        
        //cout<<"initial_bend_energy  "<<initial_bend_energy<<endl;
        //calculating the final state energy
        ////calculating the bond energy
        initial_or_final=1;
        final_bond_energy=calculating_the_bond_energy(initial_pair, initial_or_final, atoms, number_of_privious_mem_nodes);
        //cout<<"final_bond_energy  "<<final_bond_energy<<endl;
        ////calculating the bend energy of  newpairs
       // cout<<"calculating the bend energy of  newpairs"<<endl;
        uncommon1=Triangle_Pair_Nodes[initial_pair][1];
        common2=Triangle_Pair_Nodes[initial_pair][0];
        common3=Triangle_Pair_Nodes[initial_pair][3];
        uncommon4=Triangle_Pair_Nodes[initial_pair][2];
        final_bend_energy+=calculating_the_bend_energy_2(uncommon1, common2, common3, uncommon4, atoms, number_of_privious_mem_nodes);
        //find the neighbours of new pairs
        
        int n=0;
        bool A_or_B;
        int previous_dihedral_index;
        for(int i=0; i<2; i++){
            A_or_B=0;
            previous_dihedral_index=A_neighbours_dihedral_index[i];
            find_the_new_neighbour(new_neighbour_dihedrals[n],  previous_dihedral_index ,  initial_pair, A_or_B);
             //cout<<"new_neighbour_dihedrals"<< new_neighbour_dihedrals[n][0]<<"  "<< new_neighbour_dihedrals[n][1]<<"  "<< new_neighbour_dihedrals[n][2]<<"  "<< new_neighbour_dihedrals[n][3]<<"  "<< new_neighbour_dihedrals[n][4]<<"  "<< new_neighbour_dihedrals[n][5]<<endl;
            n++;
            A_or_B=1;
            previous_dihedral_index=B_neighbours_dihedral_index[i];
            find_the_new_neighbour(new_neighbour_dihedrals[n],  previous_dihedral_index ,  initial_pair, A_or_B);
            //cout<<"new_neighbour_dihedrals"<< new_neighbour_dihedrals[n][0]<<"  "<< new_neighbour_dihedrals[n][1]<<"  "<< new_neighbour_dihedrals[n][2]<<"  "<< new_neighbour_dihedrals[n][3]<<"  "<< new_neighbour_dihedrals[n][4]<<"  "<< new_neighbour_dihedrals[n][5]<<endl;
            n++;
        }//end of for(int i=0; i<2; i++) (finding the neighbors of new configuration)

        
        
        /*
         * ////calculating the bend energy of final neighbours with the function 1
         * for(int i=0; i<4; i++){
                if(accept){
                    double energy=calculating_the_bend_energy_2(new_neighbour_dihedrals[i][0], new_neighbour_dihedrals[i][1], new_neighbour_dihedrals[i][2], new_neighbour_dihedrals[i][3], atoms, number_of_privious_mem_nodes);
                    if(energy>=0){
                    final_bend_energy+=energy;}
                    else{
                    accept=0;}
                }//end of if(accept) related to the function calculation of bending energy (1)
            }*/
            
            
        for (int i=0; i<4; i++){
            final_bend_energy+=calculating_the_bend_energy_2(new_neighbour_dihedrals[i][0], new_neighbour_dihedrals[i][1], new_neighbour_dihedrals[i][2], new_neighbour_dihedrals[i][3], atoms, number_of_privious_mem_nodes);
        }
       // cout<<"final_bend_energy  "<<final_bend_energy<<endl;
        
        double deltaE=(final_bend_energy+final_bond_energy)-(initial_bend_energy+initial_bond_energy);
        localDeltaE=deltaE;
        //cout<<"local delta E  "<<deltaE<<endl;
        double fluidity=0.3;
        
        if(deltaE>0){
            double random_number=(double) rand() / (RAND_MAX);
            double metropolice_condition =exp((-deltaE)/(0.001987*GenConst::temperature));
            //the coefficient 0.001987 is R (gas constant) and its unit is kcal*K^(-1)mole^(-1)
            //our energy calculatin unit is kcal*mol(-1)
            //cout<<"random number  "<<random_number<<" fluidity  "<<metropolice_condition<<"  temprature  "<< GenConst::temperature <<endl;
            if(random_number>metropolice_condition){
                accept=0;}
                
        }//end of if(deltaE>0)
    }// end of if(accept) (checking the conditions)
    
    //if(accept){
      //  check_before_update(triangle_A, triangle_B,new_neighbour_dihedrals,pyramid_counter, accept);
        //}
    
    if(accept){ // (updating the mem and parameters in openmm)
            Accepted_Try_Counter+=1;
            
            update_Membrane_class_and_openmm(initial_pair,triangle_A,triangle_B, new_neighbour_dihedrals, omm, bonds,  dihedrals);
            
  
    } //end of if(accept){ // (updating the mem and parameters in openmm)
}

double Membrane::calculating_the_bond_energy(int index, bool initial_or_final, MyAtomInfo  atoms[],int number_of_privious_mem_nodes){
    vector<double> bond;
    double bond_length;
    double bond_energy=0;
    int A,B;
    bond.resize(3);
    //In Triangle_Pair_Nodes nodes of the common bond are sorted in a way that locate between the uncommon nodes of triangles (indices 1 and 2) 
    if (initial_or_final==0){// initial, common bond 
        A=Triangle_Pair_Nodes[index][1];
        B=Triangle_Pair_Nodes[index][2];
        }
    else{
        A=Triangle_Pair_Nodes[index][0];
        B=Triangle_Pair_Nodes[index][3];
    }
    for(int i=0; i<3; i++){    
    bond[i]= (atoms[A+number_of_privious_mem_nodes].posInAng[i] -atoms[B+number_of_privious_mem_nodes].posInAng[i]) *OpenMM::NmPerAngstrom;
    //cout<<"Atom A  "<<atoms[A+number_of_privious_mem_nodes].posInAng[i]<<endl;
    //cout<<"Atom B  "<<atoms[B+number_of_privious_mem_nodes].posInAng[i]<<endl;
    }
    
    bond_length=sqrt(bond[0]*bond[0]+bond[1]*bond[1]+bond[2]*bond[2]);
    
    //cout<<"bond between nodes  "<<A<<"  and   "<<B<<"  lenght   "<<bond_length<<endl;
    switch (spring_model){
        case 1: //Fene
        {
            cout<<"FEne is under construction"<<endl;
            
        }
        case 2: //harmonic
        {
            bond_energy= 0.5*Spring_coefficient
                                      //  * OpenMM::KJPerKcal
                                        * OpenMM::AngstromsPerNm
                                        * OpenMM::AngstromsPerNm
                                        *(bond_length-(Average_node_pair_length* OpenMM::NmPerAngstrom))
                                        *(bond_length-(Average_node_pair_length* OpenMM::NmPerAngstrom));

        }
        case 5: //realharmonic :)) (x4harmonic is the name but its potential is really 1/2*k*x^2)
        {
            bond_energy= 0.5*Spring_coefficient
                                      //  * OpenMM::KJPerKcal
                                        * OpenMM::AngstromsPerNm
                                        * OpenMM::AngstromsPerNm
                                        *(bond_length-(Average_node_pair_length* OpenMM::NmPerAngstrom))
                                        *(bond_length-(Average_node_pair_length* OpenMM::NmPerAngstrom));

        }
    }
    return(bond_energy);
}



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

void Membrane::find_the_new_neighbour(int neighbour_id[6], int previous_dihedral_index , int initial_pair, bool A_or_B){
    //cout<<"find_the_new_neighbours"<<endl;
    int Uncommon_node;
    int uncommon1=Triangle_Pair_Nodes[initial_pair][0];
    int common2=Triangle_Pair_Nodes[initial_pair][1];
    int common3=Triangle_Pair_Nodes[initial_pair][2];
    int uncommon4=Triangle_Pair_Nodes[initial_pair][3];
    neighbour_id[4]=previous_dihedral_index;
    if(A_or_B==0){
        if(Triangle_Pair_Nodes[previous_dihedral_index][0]==common2){
          neighbour_id[0]=uncommon4;
          for (int i=1; i<4; i++){
              neighbour_id[i]= Triangle_Pair_Nodes[previous_dihedral_index][i];}
          neighbour_id[5]=1;}
        
        if(Triangle_Pair_Nodes[previous_dihedral_index][3]==common2){
          neighbour_id[3]=uncommon4;
          for (int i=0; i<3; i++){
              neighbour_id[i]= Triangle_Pair_Nodes[previous_dihedral_index][i];}
          neighbour_id[5]=1;}
        
        if(Triangle_Pair_Nodes[previous_dihedral_index][0]==common3){
          neighbour_id[0]=uncommon4;
          for (int i=1; i<4; i++){
              neighbour_id[i]= Triangle_Pair_Nodes[previous_dihedral_index][i];}
          neighbour_id[5]=0;}
          
        if(Triangle_Pair_Nodes[previous_dihedral_index][3]==common3){
          neighbour_id[3]=uncommon4;
          for (int i=0; i<3; i++){
              neighbour_id[i]= Triangle_Pair_Nodes[previous_dihedral_index][i];}
          neighbour_id[5]=0;}
        
    }
    else{
        if(Triangle_Pair_Nodes[previous_dihedral_index][0]==common2){
          neighbour_id[0]=uncommon1;
          for (int i=1; i<4; i++){
              neighbour_id[i]= Triangle_Pair_Nodes[previous_dihedral_index][i];}
          neighbour_id[5]=0;}
        
        if(Triangle_Pair_Nodes[previous_dihedral_index][3]==common2){
          neighbour_id[3]=uncommon1;
          for (int i=0; i<3; i++){
              neighbour_id[i]= Triangle_Pair_Nodes[previous_dihedral_index][i];}
          neighbour_id[5]=0;}
        
        if(Triangle_Pair_Nodes[previous_dihedral_index][0]==common3){
          neighbour_id[0]=uncommon1;
          for (int i=1; i<4; i++){
              neighbour_id[i]= Triangle_Pair_Nodes[previous_dihedral_index][i];}
          neighbour_id[5]=1;}
          
        if(Triangle_Pair_Nodes[previous_dihedral_index][3]==common3){
          neighbour_id[3]=uncommon1;
          for (int i=0; i<3; i++){
              neighbour_id[i]= Triangle_Pair_Nodes[previous_dihedral_index][i];}
          neighbour_id[5]=1;}
        
    }
    
   
        
}

void Membrane::update_Membrane_class_and_openmm(int initial_pair,int triangle_A,int triangle_B, int new_neighbour_dihedrals[4][6], MyOpenMMData* omm, Bonds* bonds, Dihedrals* dihedrals){
   // check_monte_carlo=1;
    
    int u1=Triangle_Pair_Nodes[initial_pair][0];
    int c2=Triangle_Pair_Nodes[initial_pair][1];
    int c3= Triangle_Pair_Nodes[initial_pair][2];
    int u4=Triangle_Pair_Nodes[initial_pair][3];
    //cout<<"u1 "<<u1<<"  c2  "<<c2<<"  c3  "<<c3<<"  u4  "<<u4<<endl;
    
        
    //updating Node_bond_list
    int bond_index;
    Update_Membrane(initial_pair,triangle_A,triangle_B, new_neighbour_dihedrals, bond_index);

    
    //updating the openmm
    ////updating the struct bonds and setbondparameter of spring
    vector<double> parameters;
    bonds[bond_index].atoms[0]=u1;
    bonds[bond_index].atoms[1]=u4;
    switch (spring_model){
        case 1: //Fene
            { 
                    cout<<"FEne is under construction"<<endl; //openmm setparameter
            }
        case 2: //harmonic
            {
                parameters.resize(2);
                parameters[0]=bonds[bond_index].nominalLengthInAngstroms* OpenMM::NmPerAngstrom; //I must check the bond index in case of having other objects
                parameters[1]=bonds[bond_index].stiffnessInKcalPerAngstrom2* OpenMM::KJPerKcal* OpenMM::AngstromsPerNm * OpenMM::AngstromsPerNm;
                omm->harmonic->setBondParameters(bond_index,u1,u4, parameters[0], parameters[1] );
                    
            }
    }
    
    ////updating dihedral structs and setBondParameters of Dihedral Forces.
    vector<double> bendingparameter;
    bendingparameter.resize(1);
    bendingparameter[0]=Bending_coefficient* OpenMM::KJPerKcal;
    
    dihedrals[initial_pair].atoms[0]=c2;
    dihedrals[initial_pair].atoms[1]=u1;
    dihedrals[initial_pair].atoms[2]=u4;
    dihedrals[initial_pair].atoms[3]=c3;
    omm->Dihedral[0]->setBondParameters(initial_pair,Triangle_Pair_Nodes[initial_pair], bendingparameter);
     
    for(int i=0; i<4; i++){
        int index= new_neighbour_dihedrals[i][4];
        omm->Dihedral[0]->setBondParameters(index,Triangle_Pair_Nodes[index], bendingparameter);
    }
  

 /* for(int k=0; k<Num_of_Triangles; k++){
        cout<<"newww Triangle list  "<<k<<"   "<<  Triangle_list[k][0]<<"  "<<Triangle_list[k][1]<<"   "<<Triangle_list[k][2]<<endl;}

    for(int p=0; p<Num_of_Triangle_Pairs; p++){
      cout<<"new neighbours list   "<<Triangle_pair_list[p][0]<<"  "<<Triangle_pair_list[p][1]<<endl; 
      cout<<"new neighbour nodes   "<<Triangle_Pair_Nodes[p][0]<<"  "<<Triangle_Pair_Nodes[p][1]<<"  "<<Triangle_Pair_Nodes[p][2]<<"  "<<Triangle_Pair_Nodes[p][3]<<endl;
      
    }*/
     
}

void Membrane::Update_Membrane(int initial_pair,int triangle_A,int triangle_B, int new_neighbour_dihedrals[4][6], int& bond_index){
    int u1=Triangle_Pair_Nodes[initial_pair][0];
    int c2=Triangle_Pair_Nodes[initial_pair][1];
    int c3= Triangle_Pair_Nodes[initial_pair][2];
    int u4=Triangle_Pair_Nodes[initial_pair][3];
    //cout<<"u1 "<<u1<<"  c2  "<<c2<<"  c3  "<<c3<<"  u4  "<<u4<<endl;
    
        
    //updating Node_bond_list
    
    for(int i; i<Num_of_Node_Pairs; i++){
        if((Node_Bond_list[i][0]==c2 and Node_Bond_list[i][1]==c3) ||(Node_Bond_list[i][1]==c2 and Node_Bond_list[i][0]==c3)){
            bond_index=i;
            Node_Bond_list[i][0]=u1;
            Node_Bond_list[i][1]=u4;
        }
    }
    
    
    //updating triangle list
    Triangle_list[triangle_A][0]=u1;
    Triangle_list[triangle_A][1]=c2;
    Triangle_list[triangle_A][2]=u4;
    
    
    Triangle_list[triangle_B][0]=u1;
    Triangle_list[triangle_B][1]=c3;
    Triangle_list[triangle_B][2]=u4;
    
    
    //updating Triangle_Pair_Nodes
    Triangle_Pair_Nodes[initial_pair][0]=c2;
    Triangle_Pair_Nodes[initial_pair][1]=u1;
    Triangle_Pair_Nodes[initial_pair][2]=u4;
    Triangle_Pair_Nodes[initial_pair][3]=c3;
     
    for(int i=0; i<4; i++){
        int index= new_neighbour_dihedrals[i][4];
        for (int j=0; j<4; j++){
            Triangle_Pair_Nodes[index][j]=new_neighbour_dihedrals[i][j];}
            
    //updating Triangle_Pair_list
        if(new_neighbour_dihedrals[i][5]==1){
            if(Triangle_pair_list[index][0]==triangle_A){
                Triangle_pair_list[index][0]=triangle_B;
            }else if(Triangle_pair_list[index][1]==triangle_A){
                Triangle_pair_list[index][1]=triangle_B;
            }else if(Triangle_pair_list[index][0]==triangle_B){
                Triangle_pair_list[index][0]=triangle_A;
            }else if(Triangle_pair_list[index][1]==triangle_B){
                Triangle_pair_list[index][1]=triangle_A;
            }
        }
    }
}


bool Membrane::check_Pyramid(vector <int> A_neighbours, vector <int> B_neighbours ){
    

        bool check_the_butterfly_effect=0;
        for(int i=0; i<2; i++){
            for(int j=0; j<2; j++){
            if      (Triangle_pair_list[A_neighbours[i]][0]==Triangle_pair_list[B_neighbours[j]][0]
                    or Triangle_pair_list[A_neighbours[i]][0]==Triangle_pair_list[B_neighbours[j]][1]
                    or Triangle_pair_list[A_neighbours[i]][1]==Triangle_pair_list[B_neighbours[j]][1]
                    or Triangle_pair_list[A_neighbours[i]][1]==Triangle_pair_list[B_neighbours[j]][0]){
                    check_the_butterfly_effect=1;
                }
            }
        }
        
        if(check_the_butterfly_effect){
            return (true);
        }
        else{
            return (false);
        }
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


double Membrane::calculating_the_bond_energy_check(int p1, int p2, MyAtomInfo atoms[]){
    
    vector<double> bond;
    double Avg=get_avg_node_dist();
    double bond_length;
    double bond_energy=0;
    bond.resize(3);    
    for(int i=0; i<3; i++){    
    bond[i]= (atoms[p1].posInAng[i] -atoms[p2].posInAng[i]) *OpenMM::NmPerAngstrom;
   }
    
    bond_length=sqrt(bond[0]*bond[0]+bond[1]*bond[1]+bond[2]*bond[2]);
   // cout<<"bond_length  "<<bond_length<<endl;
    //cout<<"bond_length-Avg  "<< bond_length-(Avg* OpenMM::NmPerAngstrom)<<endl;
    
    switch (spring_model){
//        case 1: //Fene
//        {
//            cout<<"FEne is under construction"<<endl;
//            
//        }
        case 2: //harmonic
        {
            bond_energy= 0.5*Spring_coefficient
                                      //  * OpenMM::KJPerKcal
                                        * OpenMM::AngstromsPerNm
                                        * OpenMM::AngstromsPerNm
                                        *(bond_length-(Avg* OpenMM::NmPerAngstrom))
                                        *(bond_length-(Avg* OpenMM::NmPerAngstrom));
            

        }
//    case 3: //realharmonic :)) (x4harmonic is the name but its potential is really 1/2*k*x^2)
//        {
//            bond_energy= 0.5*Spring_coefficient
//                                      //  * OpenMM::KJPerKcal
//                                        * OpenMM::AngstromsPerNm
//                                        * OpenMM::AngstromsPerNm
//                                        *(bond_length-(Avg* OpenMM::NmPerAngstrom))
//                                        *(bond_length-(Avg* OpenMM::NmPerAngstrom));
//
//        }
    }
    return(bond_energy);
    
}

double Membrane::calculating_the_bond_length_check(int p1, int p2, MyAtomInfo atoms[]){
    
    vector<double> bond;
    double Avg=get_avg_node_dist();
    double bond_length;
    bond.resize(3);    
    for(int i=0; i<3; i++){    
    bond[i]= (atoms[p1].posInAng[i] -atoms[p2].posInAng[i]) *OpenMM::NmPerAngstrom;
   }
    
    bond_length=sqrt(bond[0]*bond[0]+bond[1]*bond[1]+bond[2]*bond[2]);

    return(bond_length);
    
}

void Monte_Carlo_Reinitialize(MyOpenMMData* omm, Bonds* bonds, Dihedrals* dihedrals, Membrane &mem, MyAtomInfo atoms[]){
    bool preservestate=1;
     double time, initenergy,finalenergy,initpenergy, finalpenergy;
     double localDeltaE=0;
     int Accepted_Try_Counter=0;
     int pyramid_counter=0;
     for(int i=0; i<GenConst::MC_step; i++){

        myGetOpenMMState(omm, time, initenergy,initpenergy, atoms);
        mem.monte_carlo_flip(omm, bonds, dihedrals, atoms,localDeltaE, Accepted_Try_Counter, pyramid_counter);

        omm->context->reinitialize(preservestate);
        myGetOpenMMState(omm, time, finalenergy,finalpenergy, atoms);
        double globalDeltaE=finalpenergy-initpenergy;
        if(globalDeltaE!=0){
            if(abs(localDeltaE-globalDeltaE)>0.0000001){
                cout<<"warning! local and global DeltaE are different  "<<"global Delta potential energy  "<<globalDeltaE<<"   "<<"locall  "<<localDeltaE<<endl;
                cout<<"Test  "<<globalDeltaE-localDeltaE<<endl;
            }

        }
    
     }

     cout<<"num_of_accepted tries  "<<Accepted_Try_Counter<<"  out of  "<<GenConst::MC_step<<"  pyramid_counter  "<<pyramid_counter<<endl;
     
}

void Membrane::check_before_update(int triangle_A,int triangle_B, int new_neighbour_dihedrals[4][6],int& pyramid_counter, bool& accept){
    vector<int> A_neighbors;
    vector<int> B_neighbors;
    for(int i=0; i<4; i++){
        int index= new_neighbour_dihedrals[i][4];
        if(new_neighbour_dihedrals[i][5]==0){
            if(Triangle_pair_list[index][0]==triangle_A or Triangle_pair_list[index][1]==triangle_A){
                A_neighbors.push_back(index);
                
            
            }if(Triangle_pair_list[index][0]==triangle_B or Triangle_pair_list[index][1]==triangle_B){
                B_neighbors.push_back(index);
            }
        }
        else{
            //updating Triangle_Pair_list
            if(Triangle_pair_list[index][0]==triangle_A){
                Triangle_pair_list[index][0]=triangle_B;
            }else if(Triangle_pair_list[index][1]==triangle_A){
                Triangle_pair_list[index][1]=triangle_B;
            }else if(Triangle_pair_list[index][0]==triangle_B){
                Triangle_pair_list[index][0]=triangle_A;
            }else if(Triangle_pair_list[index][1]==triangle_B){
                Triangle_pair_list[index][1]=triangle_A;
            }
            
            
            if(Triangle_pair_list[index][0]==triangle_A or Triangle_pair_list[index][1]==triangle_A){
                B_neighbors.push_back(index);
            }if(Triangle_pair_list[index][0]==triangle_B or Triangle_pair_list[index][1]==triangle_B){
                A_neighbors.push_back(index);
            }
        }
     
       
    }

    
    
    if(A_neighbors.size()!=2 or B_neighbors.size()!=2){
        cout<<"warning- check before update, monte carlo flip is rejected"<<endl;
        accept=0;
    }else{
        if(check_Pyramid(A_neighbors, B_neighbors)){
            pyramid_counter+=1;
        cout<<"warning- check before update,pyramid formation, monte carlo flip is rejected"<<endl;
        accept=0;
        
        }
    
    }
    if(accept==0){
        for(int i=0; i<4; i++){
            //undo Triangle_Pair_list
            int index=new_neighbour_dihedrals[i][4];
        if(new_neighbour_dihedrals[i][5]==1){
            if(Triangle_pair_list[index][0]==triangle_A){
                Triangle_pair_list[index][0]=triangle_B;
            }else if(Triangle_pair_list[index][1]==triangle_A){
                Triangle_pair_list[index][1]=triangle_B;
            }else if(Triangle_pair_list[index][0]==triangle_B){
                Triangle_pair_list[index][0]=triangle_A;
            }else if(Triangle_pair_list[index][1]==triangle_B){
                Triangle_pair_list[index][1]=triangle_A;
            }
        }
    }

    }
}

void Membrane::check_the_flip(MyOpenMMData* omm, Bonds* bonds, Dihedrals* dihedrals){
    bool preservestate=1;
    int triangle_A, triangle_B,temp_triangle, u1;
    int initial_pair;
    vector<int> A_neighbours_dihedral_index;
    vector<int> B_neighbours_dihedral_index;
    double initial_bend_energy=0;
    double initial_bond_energy=0;
    double final_bend_energy=0;
    double final_bond_energy=0;
    int number_of_privious_mem_nodes=0;
    
    int new_neighbour_dihedrals[4][6];
    
    //** the new pair must have 4 new neighbours.
    //**each array contains the information of new neighbours. the first 4 elements are the index of nodes which constructs the dihedral. the fifths element indicates the index of the dihedral in both Triangle_Pair_Nodes and Triangle_pair_list. the last element is a boolian which 0 indicates that we dont have to update the indexes of triangles which are connected to each other and one indicates that we have to update the lable of new connected triangles.
    
    //initial_pair=rand()%(Num_of_Triangle_Pairs-1);
    //cout<<"initial pair  "<<initial_pair<<endl;
    initial_pair=2;
    triangle_A=Triangle_pair_list[initial_pair][0];
    triangle_B=Triangle_pair_list[initial_pair][1];
   


 
    /*
     * cout<<"Tiangle A  "<<triangle_A<<endl;
    
    //cout<<"Tiangle B  "<<triangle_B<<endl;
    for(int k=0; k<Num_of_Triangles; k++){
        cout<<"initial Triangle list  "<<k<<"   "<<  Triangle_list[k][0]<<"  "<<Triangle_list[k][1]<<"   "<<Triangle_list[k][2]<<endl;}

    for(int p=0; p<Num_of_Triangle_Pairs; p++){
      cout<<"initial neighbours list   "<<Triangle_pair_list[p][0]<<"  "<<Triangle_pair_list[p][1]<<endl; 
      cout<<"initial neighbour nodes   "<<Triangle_Pair_Nodes[p][0]<<"  "<<Triangle_Pair_Nodes[p][1]<<"  "<<Triangle_Pair_Nodes[p][2]<<"  "<<Triangle_Pair_Nodes[p][3]<<endl;
      
    }
    */




    //the next part is for making sure that the first uncommon node in dihedral belongs to triangle A;
    u1=Triangle_Pair_Nodes[initial_pair][0];
    if(Triangle_list[triangle_A][0]!=u1 and Triangle_list[triangle_A][1]!=u1 and Triangle_list[triangle_A][2]!=u1){
        temp_triangle=triangle_A;
        triangle_A=triangle_B;
        triangle_B=temp_triangle;
        
    }
   
   for (int i=0; i<Triangle_pair_list.size() ; i++){
        
        if ((Triangle_pair_list[i][0]==triangle_A and Triangle_pair_list[i][1]!=triangle_B)  || (Triangle_pair_list[i][1]==triangle_A and Triangle_pair_list[i][0]!=triangle_B)){

            A_neighbours_dihedral_index.push_back(i);}
        if ((Triangle_pair_list[i][0]==triangle_B and Triangle_pair_list[i][1]!=triangle_A)  || (Triangle_pair_list[i][1]==triangle_B and Triangle_pair_list[i][0]!=triangle_A)){
        
            B_neighbours_dihedral_index.push_back(i);}
        
    }
    
    ///////check the fliping 
    int n=0;
        bool A_or_B;
        int previous_dihedral_index;
        for(int i=0; i<2; i++){
            A_or_B=0;
            previous_dihedral_index=A_neighbours_dihedral_index[i];
            find_the_new_neighbour(new_neighbour_dihedrals[n],  previous_dihedral_index ,  initial_pair, A_or_B);
             //cout<<"new_neighbour_dihedrals"<< new_neighbour_dihedrals[n][0]<<"  "<< new_neighbour_dihedrals[n][1]<<"  "<< new_neighbour_dihedrals[n][2]<<"  "<< new_neighbour_dihedrals[n][3]<<"  "<< new_neighbour_dihedrals[n][4]<<"  "<< new_neighbour_dihedrals[n][5]<<endl;
            n++;
            A_or_B=1;
            previous_dihedral_index=B_neighbours_dihedral_index[i];
            find_the_new_neighbour(new_neighbour_dihedrals[n],  previous_dihedral_index ,  initial_pair, A_or_B);
            //cout<<"new_neighbour_dihedrals"<< new_neighbour_dihedrals[n][0]<<"  "<< new_neighbour_dihedrals[n][1]<<"  "<< new_neighbour_dihedrals[n][2]<<"  "<< new_neighbour_dihedrals[n][3]<<"  "<< new_neighbour_dihedrals[n][4]<<"  "<< new_neighbour_dihedrals[n][5]<<endl;
            n++;
        }
 
     
    update_Membrane_class_and_openmm(initial_pair,triangle_A,triangle_B, new_neighbour_dihedrals,  omm,  bonds,  dihedrals);
    omm->context->reinitialize(preservestate);
    cout<<"fliped!"<<endl;
}


bool Membrane::check_Pyramid_2(int initial_pair, vector<int> A_neighbours_dihedral_index, vector<int> B_neighbours_dihedral_index){
    bool ButerflyEffect=0;
    int u1= Triangle_Pair_Nodes[initial_pair][0];
    int c2= Triangle_Pair_Nodes[initial_pair][1];
    int c3= Triangle_Pair_Nodes[initial_pair][2];
    int u4= Triangle_Pair_Nodes[initial_pair][3];
    vector<int> uncommen_nodes;
    for(int i=0; i<2; i++){
        if (Triangle_Pair_Nodes[A_neighbours_dihedral_index[i]][0]!=c2 and Triangle_Pair_Nodes[A_neighbours_dihedral_index[i]][0]!=c3){
            uncommen_nodes.push_back(Triangle_Pair_Nodes[A_neighbours_dihedral_index[i]][0]);
        }
        if (Triangle_Pair_Nodes[A_neighbours_dihedral_index[i]][3]!=c2 and Triangle_Pair_Nodes[A_neighbours_dihedral_index[i]][3]!=c3){
            uncommen_nodes.push_back(Triangle_Pair_Nodes[A_neighbours_dihedral_index[i]][3]);
        }
        if (Triangle_Pair_Nodes[B_neighbours_dihedral_index[i]][0]!=c2 and Triangle_Pair_Nodes[B_neighbours_dihedral_index[i]][0]!=c3){
            uncommen_nodes.push_back(Triangle_Pair_Nodes[B_neighbours_dihedral_index[i]][0]);
        }
        if (Triangle_Pair_Nodes[B_neighbours_dihedral_index[i]][3]!=c2 and Triangle_Pair_Nodes[B_neighbours_dihedral_index[i]][3]!=c3){
            uncommen_nodes.push_back(Triangle_Pair_Nodes[B_neighbours_dihedral_index[i]][3]);
        }
        
    }
    for (int j=0; j<4; j++){
        for(int k=j+1; k<4; k++){
            if(uncommen_nodes[j]==uncommen_nodes[k]){
                ButerflyEffect=1;
            }
        }
    }
    if (ButerflyEffect){
        //cout<<"ButerflyEffect2"<<endl;
        return(true);
    }
    else{
        return(false);
    }
}

