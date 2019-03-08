//
//  ECM.hpp
//  Mem
//
//  Created by Ali Farnudi on 18/08/2018.
//  Copyright © 2018 Ali Farnudi. All rights reserved.
//

#ifndef ECM_h
#define ECM_h

#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <map>
#include "General_functions.hpp"

using namespace std;

class ECM {
private:
    
    string file_time;
    int index;
    
    map<string, double> param_map;
    
    double Node_Mass=1;
    double Node_radius=1;
    int spring_model=1;
    double Spring_coefficient=10;
    
    
    double Shift_in_X_direction=0;
    double Shift_in_Y_direction=0;
    double Shift_in_Z_direction=0;
    
    int Num_of_Nodes=0;
    int Num_of_Triangle_Pairs=0;
    int Num_of_Node_Pairs=0;
    int Num_of_Triangles=0;
    
    double Downward_speed=0;
    double Kelvin_Damping_Coefficient=100;
    double Dashpot_Viscosity=0.02;
    
    vector<vector<double> > Node_Force;
    vector<vector<double> > Node_Position;
    vector<vector<int> > Triangle_List;
    
    double interaction_range=1.0;
    double epsilon=0.6;
    double sigma=15.0;
    
    
    vector<vector<double> > Node_Velocity;
    vector<vector<int> > Node_Pair_list;
    vector<vector<int> > Node_neighbour_list;
    void Node_Bond_identifier(void);
    
    void initialise(string Mesh_file_name, int dimension);
    void set_map_parameter(string param_name, double param_value);
    void Node_neighbour_list_constructor(void);
    
//    void read_input(string input_file);
    void read_gmesh_file_2D (string gmesh_file);
    void normal_direction_Identifier(double x, double y, double z);
    
    
    
public:
    
    void import_config(string config_file_name);
    void MD_Evolution_beginning (double MD_Time_Step);
    void MD_Evolution_end (double MD_Time_Step);
    
//    void triangle_normal_calculator(int triangle_index, double ABxAC[3]);
    
    string output_file_neme;
    
    
//    ECM(string ECM_mesh_file_name, double x, double y, double z)
//    {
////        read_ECM_input(input_file_name);
//        read_gmesh_file(ECM_mesh_file_name);
//        output_file_neme=ECM_mesh_file_name ;// it is for generating trajectory file. it can be modifyed to have date and time in it.this modification can be done in main.
//        cout<<"\n\nECM class initiated"<<endl;
//        normal_direction_Identifier(x, y, z);
//        node_pair_identifier();
//        
////        ECM_Triangle_Pair_and_Edges_Identifier();
//        
//        
//        
//    }
    void set_file_time(char* buffer){
        file_time=buffer;
    }
    void set_index(int ind){
        index=ind;
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
