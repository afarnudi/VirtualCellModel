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
    double Contractile_force =10;
    int contractile_model = 1;
    double contractile_hill_co = 0;
    double Contractile_k1=50;
    double Contractile_k2=25;
    double Contractile_rmin=0;
    double Contractile_rmax=1000000;
    double act_r0factor = 1;
    
    
    double abp_force =0;
    double abp_k1=0;
    double abp_k2=0;
    double abp_rmin=0;
    double abp_rmax=1000000;
    double abp_r0factor = 1;
    double abp_hill_co = 100 ;
    double abp_model = 1;
    double abp_spring_model =2;
    double abp_Spring_coefficient =0;
    
    double MT_force =0;
    double MT_k1=0;
    double MT_k2=0;
    double MT_rmin=0;
    double MT_rmax=1000000;
    double MT_r0factor = 1;
    double MT_hill_co = 100 ;
    double MT_model = 1;
    double MT_spring_model =2;
    double MT_Spring_coefficient =0;
    
    double rescale_factor=1;
    double Damping_coefficient=0.5;
    
    double Shift_in_X_direction=0, Shift_in_Y_direction=0, Shift_in_Z_direction=0;
    double x_speed=0.0; //???
    double y_speed=0.0;
    double z_speed=0.0;
    
    int ext_force_model=0;
    double kx=10;
    double ky=10;
    double kz=10;
    
    int Num_of_Nodes=0;
    vector<vector<double> >Node_Position;
    vector<vector<double> >Node_Velocity;
    vector<vector<double> >Node_Force;
    vector<vector<int> > Pyramid_Nodes;
    vector<vector<int>> filaments;
    vector<vector<int>> abps;
    vector<vector<int>> MTs;
    
    std::string output_file_neme;
    
    vector<vector<int> > Node_Bond_list;
    int Num_of_Node_Pairs=0;
    vector<vector<int> > abp_Bond_list;
    int Num_of_abp_Pairs=0;
    vector<vector<int> > MT_Bond_list;
    int Num_of_MT_Pairs=0;
    vector<double> Node_Bond_relaxed_length;
    vector<double> abp_Bond_relaxed_length;
    vector<double> MT_Bond_relaxed_length;
    
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
    double Average_node_pair_length=0;
    double Min_abp_pair_length=1000;
    double Max_abp_pair_length=0;
    double Average_abp_pair_length=0;
    double Min_MT_pair_length=1000;
    double Max_MT_pair_length=0;
    double Average_MT_pair_length=0;
    
    
    //Private members:
    void initialise(int type);
    void read_gmesh_file (std::string gmesh_file);
    void read_gmesh_file_2 (std::string gmesh_file);
    void Node_Bond_identifier(void);
    void Node_Bond_identifier_2(void);
    void Node_Bond_relaxed_length_initialiser(void);
    void abp_Bond_relaxed_length_initialiser(void);
    void MT_Bond_relaxed_length_initialiser(void);
    double Hookian(double distance, double initial_distance);
    double Kelvin(double distance, int bond_index);
    void initialise_node_bond_relaxed_length(void);
    void initialise_abp_bond_relaxed_length(void);
    void initialise_MT_bond_relaxed_length(void);
    double Maxwell(double distance, int bond_index);
    

public:
    //Shared List
    vector<vector<vector<int> > > Actin_Membrane_shared_Node_list;
    vector<int> Num_of_Actin_Membrane_shared_Nodes;
    
    //Member headers
    void import_config(std::string config_file_name);
    void import_config(vector<string> configlines);
    void set_map_parameter(std::string param_name, double param_value);
    void Elastic_Force_Calculator(void);
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
    
    double get_act_relaxlength(int bond_number)
    {
        return Node_Bond_relaxed_length[bond_number];
    }
    
    double get_abp_relaxlength(int bond_number)
      {
          return abp_Bond_relaxed_length[bond_number];
      }
    
    double get_MT_relaxlength(int bond_number)
      {
          return MT_Bond_relaxed_length[bond_number];
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
    /**Return the number of bonds between membrane nodes.*/
    int get_num_of_abp_pairs(void){
        return Num_of_abp_Pairs;
    }
    
    int get_num_of_MT_pairs(void){
        return Num_of_MT_Pairs;
    }
    
    /**Return input spring model, used to setup the openmm system for the bonds.*/
    int get_spring_model(void){
        return spring_model;
    }
    /**Return input abp spring model, used to setup the openmm system for the bonds.*/
    int get_abp_spring_model(void){
        return abp_spring_model;
    }
    
    int get_MT_spring_model(void){
           return MT_spring_model;
       }
    
    /**Return the id(int) of the nodes in the bonds. The id of each node in the bond list is stored in the first (0) and the second (1) id slot.*/
    int get_node_pair(int bond_num, int node_id){
        return Node_Bond_list[bond_num][node_id];
    }
    /**Return the id(int) of the nodes in the bonds. The id of each node in the bond list is stored in the first (0) and the second (1) id slot.*/
    int get_abp_pair(int bond_num, int node_id){
        return abp_Bond_list[bond_num][node_id];
    }
    
    int get_MT_pair(int bond_num, int node_id){
        return MT_Bond_list[bond_num][node_id];
    }
    
    /**Return the average distance of the nodes as calculated from the mesh. */
    double get_avg_node_dist(void){
        return Average_node_pair_length;
    }
    /**Return the average distance of the nodes as calculated from the mesh. */
    double get_avg_abp_dist(void){
        return Average_abp_pair_length;
    }
    
    double get_avg_MT_dist(void){
        return Average_MT_pair_length;
    }
    
    /**Return the abp r0 factor for rest length of bond. */
       double get_abp_r0factor(void){
           return abp_r0factor;
       }
    
    double get_MT_r0factor(void){
        return MT_r0factor;
    }
    
    /**Return the actin r0 factor for rest length of bond. */
       double get_act_r0factor(void){
           return act_r0factor;
       }
    
    /**Return spring stiffness coefficient. */
    double get_spring_stiffness_coefficient(void){
        return Spring_coefficient;
    }
    /**Return abp spring stiffness coefficient. */
    double get_abp_spring_stiffness_coefficient(void){
        return abp_Spring_coefficient;
    }
    
    double get_MT_spring_stiffness_coefficient(void){
        return MT_Spring_coefficient;
    }
    
    /**Return contractile constant force. */
    double get_contractile_force(void){
        return Contractile_force;
    }
    /**Return contractile stiffness_1. */
    double get_contractile_k1(void){
        return Contractile_k1;
    }
    /**Return contractile stiffness_2. */
    double get_contractile_k2(void){
        return Contractile_k2;
    }
    /**Return contractile r_min. */
    double get_contractile_rmin(void){
        return Contractile_rmin;
    }
    /**Return contractile r_max. */
    double get_contractile_rmax(void){
        return Contractile_rmax;
    }
    
    /**Return hill coefficient. */
    double get_abp_hill_co(void){
        return abp_hill_co;
    }
    
    double get_MT_hill_co(void){
        return MT_hill_co;
    }
    
    /**Return hill coefficient. */
    double get_contractile_hill_co(void){
        return contractile_hill_co;
    }
    
    /**Return contractile type. */
    int get_contractile_model(void){
        return contractile_model;
    }
    
    /**Return abp type. */
       int get_abp_model(void){
           return abp_model;
       }
    
    int get_MT_model(void){
        return MT_model;
    }
    
    
    /**Return contractile constant force. */
    double get_abp_force(void){
        return abp_force;
    }
    double get_MT_force(void){
        return MT_force;
    }
    /**Return contractile stiffness_1. */
    double get_abp_k1(void){
        return abp_k1;
    }
    double get_MT_k1(void){
        return MT_k1;
    }
    /**Return contractile stiffness_2. */
    double get_abp_k2(void){
        return abp_k2;
    }
    double get_MT_k2(void){
        return MT_k2;
    }
    /**Return contractile r_min. */
    double get_abp_rmin(void){
        return abp_rmin;
    }
    double get_MT_rmin(void){
        return MT_rmin;
    }
    /**Return contractile r_max. */
    double get_abp_rmax(void){
        return abp_rmax;
    }
    
    double get_MT_rmax(void){
        return MT_rmax;
    }
    
    
    /**Return kelvin damping coefficient. */
    double get_kelvin_damping_coefficient(void){
        return Kelvin_Damping_Coefficient;
    }
    /**Return external force model. */
    int get_ext_force_model(void){
        return ext_force_model;
    }
    double get_kx(void){
        return kx;
    }
    double get_ky(void){
        return ky;
    }
    double get_kz(void){
        return kz;
    }
    void check(void);
    void check_2(void);
    void check_3(void);

    
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
    
    std::map<string, vector<string> > Params;
    vector<string> insertOrder;
    vector<string> values;
    Actin(){
        values.resize(2);
        
        values[0] ="value 0";
        values[1] ="#This is a parameter example for Actin with default value 'value 0'.";
        Params["ActinSampleParam0"] = values;
        insertOrder.push_back("ActinSampleParam0");
        
        values[0] ="value 1";
        values[1] ="#This is a parameter example for Actin with default value 'value 1'.";
        Params["ActinSampleParam1"] = values;
        insertOrder.push_back("ActinSampleParam1");
        
    }
    
    std::map<string, vector<string> > get_map(){
        return Params;
    }
    vector<string > get_insertOrder(){
        return insertOrder;
    }
    void assign_key_value(string key, string value){
        Params[key][0]=value;
    }
};

#endif // MEMBRANE_H
