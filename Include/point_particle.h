#ifndef POINT_PARTICLE_H
#define POINT_PARTICLE_H
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include "General_constants.h"
#include "General_functions.hpp"
#include <math.h>
#include <map>
#include <iomanip>
#include <iterator>

using namespace std;

class point_particle
{
  
private:
    string file_time;
    map<string, double> param_map;
    
    int index;
    double Node_Mass=1;
    double X_0 , Y_0 , Z_0 ;
    vector<double> Node_Position;
    vector<double>  Node_Velocity;
    vector<double>  Node_Force;
    
    void set_map_parameter(string param_name, double param_value);
    
    double P_P_sigma=2;
    double P_P_epsilon=100;
    double P_P_cut_off=5;
    
    double P_Membrane_sigma=0.8;
    double P_Membrane_epsilon=100;
    double P_Membrane_cut_off=3;
    
    
    
    
public:
    vector<double>  Neighbour;
    int Neighbournumber;
    double Neighbour_dist;
    double cut_off;
    bool on_or_off_MD_evolution=1;
    
    
    void initialize_point_particle(double x, double y, double z);
    void write_traj (string traj_name);
    void MD_Evolution_beginning(double MD_Time_Step);
    void MD_Evolution_end(double MD_Time_Step);
    void generate_report();
    void import_config(string config_file_name);
    
    void set_index(int ind){
    index=ind;
    }
    void set_file_time(char* buffer){
        file_time=buffer;
    }
    double return_P_P_sigma(void){
        return P_P_sigma ;
    }
    
    double return_P_P_epsilon(void){
        return P_P_epsilon ;
    }
    
    
    double return_P_P_cut_off(void){
        return P_P_cut_off ;
    }
    
    
    double return_P_Membrane_sigma(void){
        return P_Membrane_sigma ;
    }
    
    
    double return_P_Membrane_epsilon(void){
        return P_Membrane_epsilon ;
    }
    
    double return_P_Membrane_cut_off(void){
        return P_Membrane_cut_off ;
    }
    double return_position(int node_coordinate){
        return Node_Position[node_coordinate];
    }
    void add_to_force(double force, int coor){
        Node_Force[coor]+=force;
    }
    

};

#endif // POINT_PARTICLE_H
