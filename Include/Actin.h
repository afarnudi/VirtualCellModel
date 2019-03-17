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
//    double Damping_coefficient=0.5;
    
    double Shift_in_X_direction=0, Shift_in_Y_direction=0, Shift_in_Z_direction=0;
    double Downward_speed=0;
    
    int Num_of_Nodes=0;
    vector<vector<double> >Node_Position;
    vector<vector<double> >Node_Velocity;
    vector<vector<double> >Node_Force;
    vector<vector<int> > Pyramid_Nodes;
    
    string output_file_neme;
    
    vector<vector<int> > Node_Bond_list;
    int Num_of_Node_Pairs=0;
    vector<double> Node_Bond_relaxed_length;
    
    double Total_Potential_Energy=0;

    double Kelvin_Damping_Coefficient=100;
    double Total_Kinetic_Energy=0;
    double COM_velocity[3]={0};
    double COM_position[3]={0};
    double Dashpot_Viscosity;
    double tau_Maxwell_relax=0;
    double exp_tau=0;
    
    
    //Private members:
    void initialise(string Mesh_file_name);
    void read_gmesh_file (string gmesh_file);
    void Node_Bond_identifier(void);
    void Node_Bond_relaxed_length_initialiser(void);
    double Hookian(double distance, double initial_distance);
    double Kelvin(double distance, int bond_index);
    void initialise_node_bond_relaxed_length(void);
    double Maxwell(double distance, int bond_index);
    

public:
    //Shared List
    vector<vector<int> > Actin_Membrane_shared_Node_list;
    int Num_of_Actin_Membrane_shared_Nodes=0;
    
    //Member headers
    void import_config(string config_file_name);
    void set_map_parameter(string param_name, double param_value);
    void Elastic_Force_Calculator(void);
    void MD_Evolution_beginning (double MD_Time_Step);
    void MD_Evolution_end (double MD_Time_Step);
    void Thermostat_Bussi(double MD_T);
    void write_traj (string traj_name, string label);
    void generate_report();
    
    //General members:
    void set_file_time(char* buffer){
        file_time=buffer;
    }
    void set_index(int ind){
        index=ind;
    }
    int return_num_of_nodes(void){
        return Num_of_Nodes;
    }
    double return_node_position(int node_number, int node_coordinate){
        return Node_Position[node_number][node_coordinate];
    }
    
    
    void shift_position (double x, double y, double z){
        for (int i=0; i<Num_of_Nodes; i++) {
            Node_Position[i][0]+=x;
            Node_Position[i][1]+=y;
            Node_Position[i][2]+=z;
        }
    }
    void shift_velocity (double vx, double vy, double vz){
        for (int i=0; i<Num_of_Nodes; i++) {
            Node_Velocity[i][0]+=vx;
            Node_Velocity[i][1]+=vy;
            Node_Velocity[i][2]+=vz;
        }
    }
    int return_num_of_actin_membrane_shared_nodes(void){
        return Num_of_Actin_Membrane_shared_Nodes;
    }
    void add_to_force(double force,int index, int coor){
        Node_Force[index][coor]+=force;
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
        COM_velocity[0]/=Num_of_Nodes;
        COM_velocity[2]/=Num_of_Nodes;
        COM_velocity[1]/=Num_of_Nodes;
    }
};

#endif // MEMBRANE_H
