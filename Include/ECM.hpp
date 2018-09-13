//
//  ECM.hpp
//  Mem
//
//  Created by Ali Farnudi on 18/08/2018.
//  Copyright © 2018 Ali Farnudi. All rights reserved.
//

#ifndef ECM_hpp
#define ECM_hpp

#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include "General_functions.hpp"

using namespace std;

class ECM {
private:
    
    int Num_of_Nodes=0;
    int Num_of_Triangle_Pairs=0;
    int Num_of_Node_Pairs=0;
    int Num_of_Triangles=0;
    
    
    double interaction_range=1.0;
    double epsilon=0.6;
    double sigma=15.0;
    
    
    vector<vector<double> > Node_Velocity;
    vector<vector<int> > Node_Pair_list;
    
    void read_input(string input_file);
    void read_gmesh_file (string gmesh_file);
    void normal_direction_Identifier(double x, double y, double z);
    void node_pair_identifier(void);
    
public:
    vector<vector<double> >Node_Force;
    vector<vector<double> >ECM_Node_Position;
    vector<vector<int> > ECM_triangle_list;
    
    
    string output_file_neme;
    
    
    ECM(string ECM_mesh_file_name, double x, double y, double z)
    {
//        read_ECM_input(input_file_name);
        read_gmesh_file(ECM_mesh_file_name);
        output_file_neme=ECM_mesh_file_name ;// it is for generating trajectory file. it can be modifyed to have date and time in it.this modification can be done in main.
        cout<<"\n\nECM class initiated"<<endl;
        normal_direction_Identifier(x, y, z);
        node_pair_identifier();
        
//        ECM_Triangle_Pair_and_Edges_Identifier();
        
        
        
    }
    int return_num_of_nodes(void){
        return Num_of_Nodes;
    }
    int return_num_of_triangles(void){
        return Num_of_Triangles;
    }
    double return_interaction_range(void){
        return interaction_range;
    }
    void set_interaction_range(double range){
        interaction_range=range;
    }
    double return_epsilon(void){
        return epsilon;
    }
    double return_sigma(void){
        return sigma;
    }
};

#endif /* ECM_hpp */
