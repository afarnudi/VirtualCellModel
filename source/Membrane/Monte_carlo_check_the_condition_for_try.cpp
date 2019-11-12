#include "Membrane.h"
#include "OpenMM_structs.h"
#include "OpenMM_funcs.hpp"
#include "General_functions.hpp"

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

void Membrane::check_the_flip(MyOpenMMData* omm, Bonds* bonds, Dihedrals* dihedrals){
    bool preservestate=1;
    int triangle_A, triangle_B,temp_triangle, u1;
    int initial_pair;
    vector<int> A_neighbours_dihedral_index;
    vector<int> B_neighbours_dihedral_index;
//    double initial_bend_energy=0;
//    double initial_bond_energy=0;
//    double final_bend_energy=0;
//    double final_bond_energy=0;
//    int number_of_privious_mem_nodes=0;
    
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
//    int u1= Triangle_Pair_Nodes[initial_pair][0];
    int c2= Triangle_Pair_Nodes[initial_pair][1];
    int c3= Triangle_Pair_Nodes[initial_pair][2];
//    int u4= Triangle_Pair_Nodes[initial_pair][3];
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