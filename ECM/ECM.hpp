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
    
    int ECM_num_of_Nodes=0;
    int ECM_num_of_Triangle_Pairs=0;
    int ECM_num_of_Node_Pairs=0;
    int ECM_num_of_Triangles=0;
    
    
    
    
    
    
    vector<vector<double> >ECM_Node_Velocity;
    vector<vector<double> >ECM_Node_Force;
    vector<vector<int> > ECM_triangle_list;
    vector<vector<int> > ECM_Node_Pair_list;
    
    void read_ECM_input(string input_file);
    void read_gmesh_file (string gmesh_file);
    void ECM_Normal_direction_Identifier(double x, double y, double z);
    void ECM_Node_Pair_Identifier(void);
    
public:
    string output_file_neme;
    vector<vector<double> >ECM_Node_Position;
    
    
    ECM(string ECM_mesh_file_name, double x, double y, double z)
    {
//        read_ECM_input(input_file_name);
        read_gmesh_file(ECM_mesh_file_name);
        output_file_neme=ECM_mesh_file_name ;// it is for generating trajectory file. it can be modifyed to have date and time in it.this modification can be done in main.
        cout<<"ECM class initiated"<<endl;
        ECM_Normal_direction_Identifier(x, y, z);
        ECM_Node_Pair_Identifier();
        
//        ECM_Triangle_Pair_and_Edges_Identifier();
        
        
        
    }
    int return_num_of_nodes(void){
        return ECM_num_of_Nodes;
    }
    
};

#endif /* ECM_hpp */
