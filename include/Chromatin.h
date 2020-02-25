#ifndef CHROMATIN_H
#define CHROMATIN_H

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

using std::string;
using std::vector;

class Chromatin
{
private: //(if we define these constants as private members of the class, we can't put them in the final report)
    
    double Node_Mass=1.0;//  also use in MD loop and should not be private unless we write some functions to get it outside the class
    double Total_Potential_Energy;
    double Node_radius=1;
    
    int index;
    int spring_model=0;
    
    double Total_Kinetic_Energy;
    double Total_potential_Energy=0.0;
    double Spring_coefficient=10.0; // stretching constant
    double Damping_coefficient=0.0; // Viscosity of the Membrane. It is applied in Force calculation for the Membrane Node pairs. I have commented out these parts in the 'Membrane_Force_Calculator' because I think the current code does not need it (some energy consuming array calculations were involved).
    double Spring_force_cut_off=1000.0;
    double Shift_in_X_direction=0.0; //???
    double Shift_in_Z_direction=0.0; //???
    double Shift_in_Y_direction=0.0; //???
    double x_speed=0.0; //???
    double y_speed=0.0;
    double z_speed=0.0;
    vector<vector<int > > virtual_bond_pairs;
    
    double bond_length=0;
    double bond_radius=0;
    bool   optimise_bond_radius=false;
    int    num_virtual_sites_per_bond=0;
    int    num_of_extra_VS_bonds=0;
    int    num_of_total_bonds=0;
    
    double rescale_factor=1;
    
    double com[3]; //center of mass
    double Min_node_pair_length, Max_node_pair_length, Average_node_pair_length;
    
    string output_file_name;
    string file_time;
    string label;
    
    vector<vector<int> > Vsite_and_bindings;
    vector<vector<double> > Vsite_binding_weights;
    
    vector<vector<double> > Node_Position;
    vector<vector<double> > Node_Velocity;// also update in MD loop and should not be private unless we write some functions to get it outside the class
    vector<vector<double> > Node_Force;// also update in MD loop and should not be private unless we write some functions to get it outside the class
    /*variables*/
    
    vector<int> ABC_index;
    
    /**Different node types may reside on a chain with customised interactions. Default 1.*/
    int num_of_node_types=1;
    /**The Lenard Jones 12 6 epsilon for different node types*/
    vector<double> epsilon_LJ;
    /**The Lenard Jones 12 6 sigma for different node types*/
    vector<double> sigma_LJ;
    
    int Num_of_Nodes=0;
    /*constants*/
    //This is the number of nodes on the membrane (Both the outer membrane and the Nucleus). This is the first number that appears in the 'membrane' file (once opened with a text editor)
    std::map<string, double> param_map;
    
//    double Average_Node_Distance();
    void read_membrane_input(string input_file);

    void potential_1 (void);
    void potential_2 (void);
//    void FENE (void);
    void check(void);
    void calculate_mesh_properties(void);
    void node_distance_correction(void);
    void initialise(void);
    void initialise(double min_radius);
    void Pack(double min_radius);
    double chromatin_prepack(void);
    void packing_potential(double Sphere_Radius);
    void packing_traj (void);
    void reset_com_velocity(void);
    void rescale_velocities(double scale);
    int generate_virtual_sites(void);
    
public: //these are using in Monte Carlo flip function. for defining them as private variables, we have tow ways: defining monte_carlo_flip as a member of this class or writing some functions to make them accessible out of membrane class.
    
    void pdb_label_check(void);
    
    double COM_velocity[3];
    double COM_position[3];

    vector<vector<int> > Node_neighbour_list;
    vector<vector<int> > Membrane_neighbour_node;
//    vector<vector<double> > Contact_Matrix;
    
    
    
    void MD_Evolution_beginning (double MD_Time_Step);
    void MD_Evolution_end (double MD_Time_Step);
    void Node_neighbour_list_constructor();
    void export_for_resume(int MD_step);
    void export_for_resume(int MD_step, MyAtomInfo atoms[], int atom_count);
    void write_traj (string traj_name, string label);
    
    void import_resume(string import_file_name);
    void import_coordinates(string import_file_name);
    void import_config(string config_file_name);
    void import_config(string config_file_name, double min_radius);
    void set_map_parameter(string param_name, double param_value);
    void generate_report(void);
    void Thermostat_2(double MD_KT);
    void Thermostat_N6(double MD_KT);
    void Thermostat_Bussi(double MD_T);
    void Results (string label);
    void build_random_chain(void);
    void Force_Calculator();
    void Force_Calculator_2();
    void FENE(void);
    void hard_sphere (void);
    void Strong_spring(void);
//    void write_parameters(int MD_Step);
    void export_pack(int MD_step);
    /**Set the current state (OpenMM) of the class.*/
    void set_state(MyAtomInfo all_atoms[], int atom_count);
    
    void calculate_extra_virtual_bonds(void);
    //=========================================================================================================
    //Function definitions:
    //=========================================================================================================
    
    /**return the chromatin node-node distance.*/
    double get_bond_length(void){
        return bond_length;
    }
    
    /**return a 2 dimentional list where the first column of each row is the index  of the virtual site on the chromatin and the second and third column correspond to the real sites that are needed for coordiante calculations. the weights can be obtained through the "get_Vsite_binding_weight_list" member.  */
    vector<vector<int> > get_Vsite_and_bindings_list(void){
        return Vsite_and_bindings;
    }
    
    /**return the weights used to calculate the coordinates of the virtual sites. the coresponding real node indecies can be obtained through the "get_Vsite_and_bindings_list" member.*/
    vector<vector<double> > get_Vsite_binding_weight_list(void){
        return Vsite_binding_weights;
    }
    
    
    /**return chromatin number of bonds. If there are no virtual site this number is  Num_of_Nodes-1. In the presence of Virtual Sites this number is reduced to [(Num_of_Nodes + num_virtual_sites_per_bond)/(num_virtual_sites_per_bond+1)]-1*/
    double get_num_of_bonds(void){
        if (!GenConst::ChromatinVirtualSites) {
            num_of_total_bonds = Num_of_Nodes-1;
            return num_of_total_bonds;
        } else {
            if (num_of_extra_VS_bonds==0){
                calculate_extra_virtual_bonds();
            }
            return num_of_total_bonds;
        }
    }
    
    /**return the number of  chromatin's real (non virtual) nodes .*/
    int get_num_of_real_site(void){
        return (Num_of_Nodes+num_virtual_sites_per_bond)/(num_virtual_sites_per_bond+1);
    }
    
    /**return the number of  chromatin's real (non virtual) nodes .*/
    int get_num_of_virtual_sites_per_bond(void){
        return num_virtual_sites_per_bond;
    }
    
    /**return the list of all virtual site pairs that should be excluded from non bonded interactions.*/
    vector< vector < int> > get_virtual_bonds (void){
        return virtual_bond_pairs;
    }
    
    /**return the assigned sigma to the provided node type.*/
    double get_sigma_LJ_12_6(int node_type){
        return sigma_LJ[node_type];
    }
    /**return the assigned sigma to the provided node type.*/
    double get_epsilon_LJ_12_6(int node_type){
        return epsilon_LJ[node_type];
    }
    
    
    /** Assigns the label(pdb) used to write to the trajectory files.
     */
    void set_label(std::string lab){
        label=lab;
    }
    
    /** return the label(pdb) used to write to the trajectory files.
     */
    std::string get_label(void){
        return label;
    }
    /**Return the node mass. At the current stage of this code all membrane nodes have the same mass.
     */
    double get_node_mass(void){
        return Node_Mass;
    }
    /**Return input spring model, used to setup the openmm system for the bonds.
     */
    int get_spring_model(void){
        return spring_model;
    }
    /**Return the number of node types in a chromatin chain.
     */
    int get_num_of_node_types(){
        return num_of_node_types;
    }
    /**Return node type. node types: A (0), B (1)
     */
    int get_node_type(int index){
        return ABC_index[index];
    }
    /**Return spring stiffness coefficient.
     */
    double get_spring_stiffness_coefficient(void){
        return Spring_coefficient;
    }
    void shift_node_positions(void){
        for (int i=0; i<Num_of_Nodes; i++) {
            Node_Position[i][0]+=Shift_in_X_direction;
            Node_Position[i][1]+=Shift_in_Y_direction;
            Node_Position[i][2]+=Shift_in_Z_direction;
        }
    }
    void shift_node_velocities(void){
        for (int i=0; i<Num_of_Nodes; i++) {
            Node_Velocity[i][0]+=x_speed;
            Node_Velocity[i][1]+=y_speed;
            Node_Velocity[i][2]+=z_speed;
        }
    }
    int get_num_of_nodes(void){
        return Num_of_Nodes;
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
    /**Returns the x (0), y (1), and z (2) velocities of the node index (number).*/
    double get_node_velocity(int node_number, int node_coordinate){
        return Node_Velocity[node_number][node_coordinate];
    }
    double get_node_position(int node_number, int node_coordinate){
        return Node_Position[node_number][node_coordinate];
    }
    
    void calculate_average_force(void){
        double average_force_x=0, average_force_y=0, average_force_z=0;
        for(int j=0 ; j<Num_of_Nodes ; j++){
            average_force_x+=Node_Force[j][0];
            average_force_y+=Node_Force[j][1];
            average_force_z+=Node_Force[j][2];
            
        }
    }

    
    void set_file_time(char* buffer){
        file_time=buffer;
    }
    void set_index(int ind){
        index=ind;
    }
    
    double get_node_radius(void){
        return Node_radius;
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

#endif // CHROMATIN_H

