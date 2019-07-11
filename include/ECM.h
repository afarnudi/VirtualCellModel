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
#include <iomanip>
#include <iterator>
#include "General_functions.hpp"

using std::vector;
using std::cout;
using std::endl;

class ECM {
private:
    
    std::string label;
    
    std::string file_time;
    int index;
    
    std::map<std::string, double> param_map;
    
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
    double COM_velocity[3]={0};
    double COM_position[3]={0};
    
    vector<vector<double> > Node_Velocity;
    vector<vector<int> > Node_Pair_list;
    vector<vector<int> > Node_neighbour_list;
    void Node_Bond_identifier(void);
    
    void initialise(std::string Mesh_file_name, int dimension);
    void set_map_parameter(std::string param_name, double param_value);
    void Node_neighbour_list_constructor(void);
    
//    void read_input(string input_file);
    void read_gmesh_file_2D (std::string gmesh_file);
    void read_gmesh_file_3D (std::string gmesh_file);
    void normal_direction_Identifier (double x, double y, double z);
    void normal_direction_Identifier (void);
    
    
public:
    
    void import_config(std::string config_file_name);
    void MD_Evolution_beginning (double MD_Time_Step);
    void MD_Evolution_end (double MD_Time_Step);
    void write_traj (std::string traj_name, std::string label);
    void generate_report(void);
    
    
    /** Assigns the label(pdb) used to write to the trajectory files. */
    void set_label(std::string lab){
        label=lab;
    }
    /** return the label(pdb) used to write to the trajectory files. */
    std::string get_label(void){
        return label;
    }
    
    std::string output_file_neme;
   
    void set_file_time(char* buffer){
        file_time=buffer;
    }
    void set_index(int ind){
        index=ind;
    }
    
    /** \brief public access to total number of ECM nodes.
     * \return integer number of nodes in the ECM.
     */
    int get_num_of_nodes(void){
        return Num_of_Nodes;
    }
    /**public access to the total number of triangles in the ECM.*/
    int get_num_of_triangles(void){
        return Num_of_Triangles;
    }
    /**Return the node mass. At the current stage of this code all membrane nodes have the same mass. */
    double get_node_mass(void){
        return Node_Mass;
    }
    double get_node_radius(void){
        return Node_radius;
    }
    double get_interaction_range(void){
        return interaction_range;
    }
    void set_interaction_range(double range){
        interaction_range=range;
    }
    double get_epsilon(void){
        return epsilon;
    }
    double get_sigma(void){
        return sigma;
    }
    
    double get_node_position(int node_number, int node_coordinate){
        return Node_Position[node_number][node_coordinate];
    }
    void add_to_force(double force,int index, int coor){
        Node_Force[index][coor]+=force;
    }
    void shift_node_positions(void){
        for (int i=0; i<Num_of_Nodes; i++) {
            Node_Position[i][0]+=Shift_in_X_direction;
            Node_Position[i][1]+=Shift_in_Y_direction;
            Node_Position[i][2]+=Shift_in_Z_direction;
        }
    }
    void update_COM_velocity(void){
        COM_velocity[0]=0;
        COM_velocity[1]=0;
        COM_velocity[2]=0;
        for (int i=0; i<Num_of_Nodes; i++) {
            COM_velocity[0]+=Node_Velocity[i][0];
            COM_velocity[1]+=Node_Velocity[i][1];
            COM_velocity[2]+=Node_Velocity[i][2];
        }
        COM_velocity[0]/=Num_of_Nodes;
        COM_velocity[2]/=Num_of_Nodes;
        COM_velocity[1]/=Num_of_Nodes;
    }
    void update_COM_position(void){
        COM_position[0]=0;
        COM_position[1]=0;
        COM_position[2]=0;
        for (int i=0; i<Num_of_Nodes; i++) {
            COM_position[0]+=Node_Position[i][0];
            COM_position[1]+=Node_Position[i][1];
            COM_position[2]+=Node_Position[i][2];
        }
        COM_position[0]/=Num_of_Nodes;
        COM_position[2]/=Num_of_Nodes;
        COM_position[1]/=Num_of_Nodes;
    }
};

#endif /* ECM_hpp */
