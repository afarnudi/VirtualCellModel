#ifndef ACTIN_H
#define ACTIN_H
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <map>
#include <iomanip>
#include <iterator>

#include "General_functions.hpp"

using namespace std;

class Actin
{
private:
    string file_time;
    int index;
    
    map<string, double> param_map;
    
    double Node_Mass=1;
    double Node_radius=1;
    int spring_model=2;
    double Spring_coefficient=100;
    double Damping_coefficient=0.5;
    
    double Shift_in_X_direction=0, Shift_in_Y_direction=0, Shift_in_Z_direction=0;
    double Downward_speed=0;
    
    int Num_of_Nodes=0;
    vector<vector<double> >Node_Position;
    vector<vector<double> >Node_Velocity;
    vector<vector<double> >Node_Force;
    
    vector<vector<int> > Pyramid_Nodes;
    
    
    //Private members:
    void initialise(string Mesh_file_name);
    void read_gmesh_file (string gmesh_file);
    
public:
    //Member headers
    void import_config(string config_file_name);
    void set_map_parameter(string param_name, double param_value);
    
    
    //General members:
    void set_file_time(char* buffer){
        file_time=buffer;
    }
    void set_index(int ind){
        index=ind;
    }
};

#endif // MEMBRANE_H
