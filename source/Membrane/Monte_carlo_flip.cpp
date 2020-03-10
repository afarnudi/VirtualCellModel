#include "Membrane.h"
#include "OpenMM_structs.h"
#include "OpenMM_funcs.hpp"
#include "General_functions.hpp"

bool Membrane::monte_carlo_flip(MyOpenMMData* omm, Bonds* bonds, Dihedrals* dihedrals, MyAtomInfo atoms[], double& localDeltaE, int& Accepted_Try_Counter, int& pyramid_counter){
   
    bool accept = false, update = false;
    
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
    u1 = Triangle_Pair_Nodes[initial_pair][0];
    
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
       cout<<"not enough neighbours,    num of naighbours  "<<A_neighbours_dihedral_index.size()<<"  and  "<<B_neighbours_dihedral_index.size()<<endl;
        cout<<"index  "<<initial_pair<<endl;
    }//end of else 
    
    
    
    if (accept){ //(checking the conditions)
        
        int uncommon1,common2,common3, uncommon4;
        //calculating the initial state energy
        ////calculating the bond energy
        bool initial_or_final=0;
        initial_bond_energy=calculating_the_bond_energy(initial_pair, initial_or_final, atoms, number_of_privious_mem_nodes);
        //cout<<"initial_bond_energy  "<<initial_bond_energy<<endl;
   
        ////calculating the bend energy of initial pairs
        //cout<<"calculating the bend energy of initial pairs"<<endl;
        uncommon1 = Triangle_Pair_Nodes[initial_pair][0];
        common2   = Triangle_Pair_Nodes[initial_pair][1];
        common3   = Triangle_Pair_Nodes[initial_pair][2];
        uncommon4 = Triangle_Pair_Nodes[initial_pair][3];
        initial_bend_energy+=calculating_the_bend_energy_2(uncommon1, common2, common3, uncommon4, atoms, number_of_privious_mem_nodes);
    
        ////calculating the bend energy of initial neighbours.
        for(int i=0; i<2; i++){
            //triangle_A neighbours
            uncommon1 = Triangle_Pair_Nodes[A_neighbours_dihedral_index[i]][0];
            common2   = Triangle_Pair_Nodes[A_neighbours_dihedral_index[i]][1];
            common3   = Triangle_Pair_Nodes[A_neighbours_dihedral_index[i]][2];
            uncommon4 = Triangle_Pair_Nodes[A_neighbours_dihedral_index[i]][3];
          //  cout<<"previous neigbours   "<<uncommon1<<"  "<< common2<<" "<<common3<<"  "<<uncommon4<<"  "<<A_neighbours_dihedral_index[i]<<endl;
            initial_bend_energy+=calculating_the_bend_energy_2(uncommon1, common2, common3, uncommon4,atoms, number_of_privious_mem_nodes);
            //triangle_B neighbours
            uncommon1 = Triangle_Pair_Nodes[B_neighbours_dihedral_index[i]][0];
            common2   = Triangle_Pair_Nodes[B_neighbours_dihedral_index[i]][1];
            common3   = Triangle_Pair_Nodes[B_neighbours_dihedral_index[i]][2];
            uncommon4 = Triangle_Pair_Nodes[B_neighbours_dihedral_index[i]][3];
           // cout<<"previous neigbours   "<<uncommon1<<"  "<< common2<<" "<<common3<<"  "<<uncommon4<<"  "<<B_neighbours_dihedral_index[i]<<endl;
            initial_bend_energy+=calculating_the_bend_energy_2(uncommon1, common2, common3, uncommon4, atoms, number_of_privious_mem_nodes);
        } //end of for(int i=0; i<2; i++) (calculating the bending energy of initial neighbors)
        
        //calculating the final state energy
        ////calculating the bond energy
        initial_or_final=1;
        final_bond_energy = calculating_the_bond_energy(initial_pair, initial_or_final, atoms, number_of_privious_mem_nodes);
        
        ////calculating the bend energy of  newpairs
        uncommon1 = Triangle_Pair_Nodes[initial_pair][1];
        common2   = Triangle_Pair_Nodes[initial_pair][0];
        common3   = Triangle_Pair_Nodes[initial_pair][3];
        uncommon4 = Triangle_Pair_Nodes[initial_pair][2];
        
        final_bend_energy+=calculating_the_bend_energy_2(uncommon1, common2, common3, uncommon4, atoms, number_of_privious_mem_nodes);
        //find the neighbours of new pairs
        
        int n=0;
        bool A_or_B;
        int previous_dihedral_index;
        for(int i=0; i<2; i++){
            A_or_B=0;
            previous_dihedral_index = A_neighbours_dihedral_index[i];
            find_the_new_neighbour(new_neighbour_dihedrals[n],  previous_dihedral_index ,  initial_pair, A_or_B);
             //cout<<"new_neighbour_dihedrals"<< new_neighbour_dihedrals[n][0]<<"  "<< new_neighbour_dihedrals[n][1]<<"  "<< new_neighbour_dihedrals[n][2]<<"  "<< new_neighbour_dihedrals[n][3]<<"  "<< new_neighbour_dihedrals[n][4]<<"  "<< new_neighbour_dihedrals[n][5]<<endl;
            n++;
            A_or_B=1;
            previous_dihedral_index=B_neighbours_dihedral_index[i];
            find_the_new_neighbour(new_neighbour_dihedrals[n],  previous_dihedral_index ,  initial_pair, A_or_B);
            //cout<<"new_neighbour_dihedrals"<< new_neighbour_dihedrals[n][0]<<"  "<< new_neighbour_dihedrals[n][1]<<"  "<< new_neighbour_dihedrals[n][2]<<"  "<< new_neighbour_dihedrals[n][3]<<"  "<< new_neighbour_dihedrals[n][4]<<"  "<< new_neighbour_dihedrals[n][5]<<endl;
            n++;
        }//end of for(int i=0; i<2; i++) (finding the neighbors of new configuration)

        
        for (int i=0; i<4; i++){
            final_bend_energy+=calculating_the_bend_energy_2(new_neighbour_dihedrals[i][0], new_neighbour_dihedrals[i][1], new_neighbour_dihedrals[i][2], new_neighbour_dihedrals[i][3], atoms, number_of_privious_mem_nodes);
        }
       
        
        double deltaE=(final_bend_energy+final_bond_energy)-(initial_bend_energy+initial_bond_energy);
        localDeltaE=deltaE;
       
        
        if(deltaE>0){
            //the coefficient 0.001987 is R (gas constant) and its unit is kcal*K^(-1)mole^(-1)
            //our energy calculatin unit is kcal*mol(-1)
            
            if( double(rand()/RAND_MAX) < exp((-deltaE)/(0.0083*GenConst::temperature))){
                Accepted_Try_Counter+=1;
                update_Membrane_class_and_openmm(initial_pair,triangle_A,triangle_B, new_neighbour_dihedrals, omm, bonds,  dihedrals);
                update =  true;
            }
                
        }//end of if(deltaE>0)
    }// end of if(accept) (checking the conditions)
    return update;
}



