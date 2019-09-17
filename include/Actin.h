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
#include "OpenMM_structs.h"

using std::vector;
using std::string;

class Actin
{
private:
    std::string file_time;
    int index;
    
    
    
    std::string label;
    std::string Mesh_file_name;
    std::map<std::string, double> param_map;
    
    double Node_Mass=1;
    double Node_radius=1;
    int spring_model=2;
    double Spring_coefficient=100;
    
    double rescale_factor=1;
    double Damping_coefficient=0.5;
    
    double Shift_in_X_direction=0, Shift_in_Y_direction=0, Shift_in_Z_direction=0;
    double Downward_speed=0;
    
    int Num_of_Nodes=0;
    vector<vector<double> >Node_Position;
    vector<vector<double> >Node_Velocity;
    vector<vector<double> >Node_Force;
    vector<vector<int> > Pyramid_Nodes;
    
    std::string output_file_neme;
    
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
    double Min_node_pair_length=1000;
    double Max_node_pair_length=0;
    double Average_node_pair_length=10;
    
    
    //Private members:
    void initialise();
    void read_gmesh_file (std::string gmesh_file);
    void Node_Bond_identifier(void);
    void Node_Bond_relaxed_length_initialiser(void);
    double Hookian(double distance, double initial_distance);
    double Kelvin(double distance, int bond_index);
    void initialise_node_bond_relaxed_length(void);
    double Maxwell(double distance, int bond_index);
    

public:
    //Shared List
    vector<vector<vector<int> > > Actin_Membrane_shared_Node_list;
    vector<int> Num_of_Actin_Membrane_shared_Nodes;
    
    //Member headers
    void import_config(std::string config_file_name);
    void set_map_parameter(std::string param_name, double param_value);
    void Elastic_Force_Calculator(void);
    void MD_Evolution_beginning (double MD_Time_Step);
    void MD_Evolution_end (double MD_Time_Step);
    void Thermostat_Bussi(double MD_T);
    void write_traj (std::string traj_name, std::string label);
    void generate_report();
    void export_for_resume(int MD_step);
    void export_for_resume(int MD_step, MyAtomInfo atoms[], int atom_count);
    /**Set the current state (OpenMM) of the class.*/
    void set_state(MyAtomInfo all_atoms[], int atom_count);
    /** Assigns the label(pdb) used to write to the trajectory files. It is also used to identify the class object throught the programme */
    void set_label(std::string lab){
        label=lab;
    }
    /** \brief public access to total number of ECM nodes.
     * \return integer number of nodes in the ECM.
     */
    int get_num_of_nodes(void){
        return Num_of_Nodes;
    }
    /** return the label(pdb) used to write to the trajectory files. */
    std::string get_label(void){
        return label;
    }
    double get_node_position(int node_number, int node_coordinate){
        return Node_Position[node_number][node_coordinate];
    }
    /**Return the node mass. At the current stage of this code all membrane nodes have the same mass. */
    double get_node_mass(void){
        return Node_Mass;
    }
    double get_node_radius(void){
        return Node_radius;
    }
    void shift_node_positions(void){
        for (int i=0; i<Num_of_Nodes; i++) {
            Node_Position[i][0]+=Shift_in_X_direction;
            Node_Position[i][1]+=Shift_in_Y_direction;
            Node_Position[i][2]+=Shift_in_Z_direction;
        }
    }
    /**Returns the x (0), y (1), and z (2) velocities of the node index (number).*/
    double get_node_velocity(int node_number, int node_coordinate){
        return Node_Velocity[node_number][node_coordinate];
    }
    /**Return the number of bonds between membrane nodes.*/
    int get_num_of_node_pairs(void){
        return Num_of_Node_Pairs;
    }
    /**Return input spring model, used to setup the openmm system for the bonds.*/
    int get_spring_model(void){
        return spring_model;
    }
    /**Return the id(int) of the nodes in the bonds. The id of each node in the bond list is stored in the first (0) and the second (1) id slot.*/
    int get_node_pair(int bond_num, int node_id){
        return Node_Bond_list[bond_num][node_id];
    }
    /**Return the average distance of the nodes as calculated from the mesh. */
    double get_avg_node_dist(void){
        return Average_node_pair_length;
    }
    /**Return spring stiffness coefficient. */
    double get_spring_stiffness_coefficient(void){
        return Spring_coefficient;
    }
    /**Return kelvin damping coefficient. */
    double get_kelvin_damping_coefficient(void){
        return Kelvin_Damping_Coefficient;
    }
    void check(void);
    
    //General members:
    void set_file_time(char* buffer){
        file_time=buffer;
    }
    void set_index(int ind){
        index=ind;
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
    int return_num_of_actin_membrane_shared_nodes(int j){
        return Num_of_Actin_Membrane_shared_Nodes[j];
    }
    
    int return_ActMem_shared_act_atom(int mem, int bond)
    {
        return Actin_Membrane_shared_Node_list[mem][bond][0];
    }
    int return_ActMem_shared_mem_atom(int mem, int bond)
    {
        return Actin_Membrane_shared_Node_list[mem][bond][1];
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
        COM_position[0]/=Num_of_Nodes;
        COM_position[2]/=Num_of_Nodes;
        COM_position[1]/=Num_of_Nodes;
    }
};

#endif // MEMBRANE_H
