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
    std::string MeshType;
    std::map<std::string, double> param_map;
    
    double Node_Mass=1;
    string Node_radius_stat;
    vector<double> Node_radius;
    int spring_model=2;
    double Spring_coefficient=100;
    double Contractile_force =0;
    int contractile_model = 1;
    double contractile_hill_co = 0;
    double Contractile_k1=0;
    double Contractile_k2=0;
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
    
    vector<double> Shift_position_xyzVector;
    vector<double> Shift_velocities_xyzVector;
//    double Shift_in_X_direction=0, Shift_in_Y_direction=0, Shift_in_Z_direction=0;
//    double x_speed=0.0; //???
//    double y_speed=0.0;
//    double z_speed=0.0;
    
    int ext_force_model=0;
    vector<double> ext_force_rigidity ;
    //double kx=10;
    //double ky=10;
    //double kz=10;
    
    void read_actin_file(void);
    
    /**Construct a vector that lists all the neighbours of  each node and their respective index in the 'Node_Bond_list'.
     */
    void Node_neighbour_list_constructor();
    /**This list is used to store all the neighbours of nodes:
     *Example:
     *If node i has 4 neighbours:  m, n, o, p
     *Then we find these neighbours  in the 'Node_neighbour_list' as follows:
     *Node_neighbour_list [ i ][0] = m
     *Node_neighbour_list [ i ][1] = n
     *Node_neighbour_list [ i ][2] = o
     *Node_neighbour_list [ i ][3] = p
     */
    vector<vector<int> > Node_neighbour_list;
    /**This list is used to cross reference the neighbouring nodes to their respective index in the 'Node_neighbour_list'.
     *Example:
     *Suppose  node i has 4 neighbours:  m, n, o, p
     *Then  we find  these  neighbours  in the 'Node_neighbour_list' as follows:
     *Node_neighbour_list [ i ][0] = m
     *Node_neighbour_list [ i ][1] = n
     *Node_neighbour_list [ i ][2] = o
     *Node_neighbour_list [ i ][3] = p
     *
     *Then the i'th index of the 'Node_neighbour_list_respective_bond_index' will also contain 4 elements:
     *Node_neighbour_list_respective_bond_index [ i ][0] = index_1
     *Node_neighbour_list_respective_bond_index [ i ][1] = index_2
     *Node_neighbour_list_respective_bond_index [ i ][2] = index_3
     *Node_neighbour_list_respective_bond_index [ i ][3] = index_4
     *
     *The index_#s will point to the index in the 'Node_Bond_list' that has stored the respective nodes as pairs:
     *Node_Bond_list [ index_1 ][ 0 ] and Node_Bond_list [ index_1 ][ 1 ]  = i, m
     *Node_Bond_list [ index_2 ][ 0 ] and Node_Bond_list [ index_2 ][ 1 ]  = i, n
     *Node_Bond_list [ index_3 ][ 0 ] and Node_Bond_list [ index_3 ][ 1 ]  = i, o
     *Node_Bond_list [ index_4 ][ 0 ] and Node_Bond_list [ index_4 ][ 1 ]  = i, p
     */
    vector<vector<int> > Node_neighbour_list_respective_bond_index;
    
    string mesh_format;
    
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
    string  Node_Bond_Nominal_Length_stat;
    double  Node_Bond_user_defined_Nominal_Length_in_Nm;
    vector<double> Node_Bond_Nominal_Length_in_Nm;
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
    void initialise(std::string Mesh_file_name);
    void read_gmesh_file (std::string gmesh_file);
    void read_gmesh_file_2 (std::string gmesh_file);
    void Node_Bond_identifier(void);
    void Node_Bond_identifier_2(void);
    void Node_Bond_identifier_3(void);
    
    double Hookian(double distance, double initial_distance);
    double Kelvin(double distance, int bond_index);
    void set_bond_nominal_length(void);
    void set_node_radius(void);
    void initialise_node_bond_relaxed_length(void);
    void initialise_abp_bond_relaxed_length(void);
    void initialise_MT_bond_relaxed_length(void);
    double Maxwell(double distance, int bond_index);
    

public:
    //Shared List
    vector<vector<vector<int> > > Actin_Membrane_shared_Node_list;
    vector<int> Num_of_Actin_Membrane_shared_Nodes;
    
    //Member headers

    
    
    void import_config(vector<string> configlines);
    void Elastic_Force_Calculator(void);
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
    
    double get_node_pair_Nominal_Length_in_Nm(int bond_number)
    {
        return Node_Bond_Nominal_Length_in_Nm[bond_number];
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
    double get_node_radius(int index){
        return Node_radius[index];
    }
    /** This function shifts the whole actin.*/
    void shift_position (double x, double y, double z) {
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
        return ext_force_rigidity[0];
    }
    double get_ky(void){
        return ext_force_rigidity[1];
    }
    double get_kz(void){
        return ext_force_rigidity[2];
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
        
        values[0] ="Path/to/my/meshfile.extension";
        values[1] ="#Path to the mesh file. Supported formats: Blender's ply, actin format .actin, and Gmsh 2. The Actin class cannot be initilised without a mesh file.";
        Params["MeshFile"] = values;
        insertOrder.push_back("MeshFile");
        
        values[0] ="normal";
        values[1] ="#Mesh type defined by Sajjad takes options 'normal' and 'sajjad'.  default value normal";
        Params["MeshType"] = values;
        insertOrder.push_back("MeshType");
        
        values[0] ="1";
        values[1] ="#Mass asssigned to each node. Default value 1";
        Params["NodeMass"] = values;
        insertOrder.push_back("NodeMass");
        
        values[0] ="Av";
        values[1] ="#Radius assigned to each node. Legal Values 1) Av, half of the average bond distance of nodes or 2) \"a value\". The node radius is used to calculate the cutt-off and minimum energy distance for the 'Excluded Volume' and the 'Lennard-Jones' potential.";
        Params["NodeRadius"] = values;
        insertOrder.push_back("NodeRadius");
        
//        values[0] ="Av";
//        values[1] ="#Set the rest length of the mesh springs using: 1) Au: The initial bond lengths of the mesh; 2) Av: The average node pair lengths; 3) \"value\": Where you type a specific  value.";
//        Params["NominalLengthInNm"] = values;
//        insertOrder.push_back("NominalLengthInNm");
        
        values[0] ="H";
        values[1] ="#Set the bond potential. 'H' for harmonic. 'FENE' for a finitely extensible nonlinear elastic model. 'kelvin' for the Kelvin-Voigt potential. 'N' for no potential. Default H";
        Params["SpringModel"] = values;
        insertOrder.push_back("SpringModel");
        
        values[0] ="0";
        values[1] ="#Set the bond potential rigidity coefficient. Default value 0.";
        Params["SpringCoeff"] = values;
        insertOrder.push_back("SpringCoeff");
        
        
        values[0] ="0";
        values[1] ="#Set the Contractile force for filaments in cytoskeleton. Default value 0.";
        Params["Contractile_force"] = values;
        insertOrder.push_back("Contractile_force");
        
        values[0] ="0";
        values[1] ="#Set the spring coefficient for spring parallel to contractile element.this spring exists in compression. Default value 0.";
        Params["Contractile_k1"] = values;
        insertOrder.push_back("Contractile_k1");
        
        values[0] ="0";
        values[1] ="#Set the spring coefficient for spring parallel to contractile element. this spring exists in extension. Default value 0.";
        Params["Contractile_k2"] = values;
        insertOrder.push_back("Contractile_k2");
        
        values[0] ="0";
        values[1] ="#Set the minimum length of contractile filaments. Default value 0.";
        Params["Contractile_rmin_factor"] = values;
        insertOrder.push_back("Contractile_rmin_factor");
        
        values[0] ="100";
        values[1] ="#Set the maximum length of contractile filaments. Default value 100.";
        Params["Contractile_rmax_factor"] = values;
        insertOrder.push_back("Contractile_rmax_factor");
        
        
        values[0] ="1";
        values[1] ="#Set the Contractile elements nominal length in sajjad type actin. Default value 1.";
        Params["actin_r0factor"] = values;
        insertOrder.push_back("actin_r0factor");
        
        values[0] ="0";
        values[1] ="#Set hill coefficient for hill contractiles. Default value 0.";
        Params["Contractile_hill_co"] = values;
        insertOrder.push_back("Contractile_hill_co");
        
        values[0] ="0";
        values[1] ="#Under development. Do not use this flag.";
        Params["ExtForceModel"] = values;
        insertOrder.push_back("ExtForceModel");
        
        values[0] ="0 0 0";
        values[1] ="#Under development. Do not use this flag.";
        Params["ExtForceRigidity"] = values;
        insertOrder.push_back("ExtForceRigidity");
        
        
        values[0] ="0";
        values[1] ="#Set the Kelvin Damping Coefficient. Default value 0.";
        Params["KelvinDampingCoeff"] = values;
        insertOrder.push_back("KelvinDampingCoeff");
        
        values[0] ="0";
        values[1] ="#Set the Dashpot Viscosity. Default value 0.";
        Params["DashpotViscosity"] = values;
        insertOrder.push_back("DashpotViscosity");
        
        values[0] ="Au";
        values[1] ="#Set the rest length of the mesh springs using: 1) Au: The initial bond lengths of the mesh; 2) Av: The average node pair lengths; 3) \"value\": Where you type a specific  value.";
        Params["NominalLengthInNm"] = values;
        insertOrder.push_back("NominalLengthInNm");
        
        values[0] ="1";
        values[1] ="#Used to scale the Actin coordinates. Default 1";
        Params["Scale"] = values;
        insertOrder.push_back("Scale");
        
        values[0] ="0 0 0";
        values[1] ="#X, Y, Z components of a vector used to translate all coordinates of the Mesh befor beginning the simluation.";
        Params["CoordinateTranslateVector"] = values;
        insertOrder.push_back("CoordinateTranslateVector");
        
        values[0] ="0 0 0";
        values[1] ="#Vx, Vy, Vz components of a vector used to add to all initial node velocities befor beginning the simluation.";
        Params["VelocityShiftVector"] = values;
        insertOrder.push_back("VelocityShiftVector");
        
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
    void assign_parameters(void);
    void consistancy_check(void);
    
    
    
    
    
};

#endif // MEMBRANE_H
