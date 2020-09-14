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
    /**Trnaslate the initial location of all nodes before the simulaiton begins.*/
    vector<double> Shift_position_xyzVector;
    /**Add velocity vector to the initial velocity of nodes before the simulaiton begins.*/
    vector<double> Shift_velocities_xyzVector;
//    double Shift_in_X_direction=0.0; //???
//    double Shift_in_Z_direction=0.0; //???
//    double Shift_in_Y_direction=0.0; //???
//    double x_speed=0.0; //???
//    double y_speed=0.0;
//    double z_speed=0.0;
    vector<vector<int > > virtual_bond_pairs;
    
    double bond_length=0;
    double bond_radius=0;
    bool   optimise_bond_radius=false;
    int    num_virtual_sites_per_bond=0;
    int    num_of_extra_VS_bonds=0;
    int    num_of_total_bonds=0;
    
    double rescale_factor=1;
    
    bool importcoordiantes = false;
    string import_file_name;
    bool ExportGeneratedCoordinates = false;
    
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
//    void Pack(double min_radius);
//    double chromatin_prepack(void);
//    void packing_potential(double Sphere_Radius);
//    void packing_traj (void);
//    void reset_com_velocity(void);
    int generate_virtual_sites(void);
    
public: //these are using in Monte Carlo flip function. for defining them as private variables, we have tow ways: defining monte_carlo_flip as a member of this class or writing some functions to make them accessible out of membrane class.
    
    void pdb_label_check(void);
    
    double COM_velocity[3];
    double COM_position[3];

    vector<vector<int> > Node_neighbour_list;
//    vector<vector<int> > Membrane_neighbour_node;
//    vector<vector<double> > Contact_Matrix;
    
    
    
//    void MD_Evolution_beginning (double MD_Time_Step);
//    void MD_Evolution_end (double MD_Time_Step);
    void Node_neighbour_list_constructor();
    void export_for_resume(int MD_step);
    void export_for_resume(int MD_step, MyAtomInfo atoms[], int atom_count);
    
    
    void import_resume(string import_file_name);
    void import_coordinates(string import_file_name);
    void import_config(string config_file_name);
    void import_config(vector<string> configlines);
    void set_map_parameter(string param_name, double param_value);
    void generate_report(void);
//    void Thermostat_2(double MD_KT);
//    void Thermostat_N6(double MD_KT);
//    void Thermostat_Bussi(double MD_T);
    void Results (string label);
    void build_random_chain(void);
    void random_walk_gen(double velocity_COM[3]);
    void Force_Calculator();
    void Force_Calculator_2();
    void FENE(void);
//    void hard_sphere (void);
//    void Strong_spring(void);
//    void write_parameters(int MD_Step);
//    void export_pack(int MD_step);
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
    
    /** Export position and velocities to txt file. This file can be later imported to initiate the chromatin.*/
    void export_coordinates(void);
    
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
        if ( (spring_model == 2) && (Spring_coefficient < 0.00001) ) {
            return 0;
        }
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
            Node_Position[i][0]+=Shift_position_xyzVector[0];
            Node_Position[i][1]+=Shift_position_xyzVector[1];
            Node_Position[i][2]+=Shift_position_xyzVector[2];
        }
//        std::cout<<"Chromatin positions shifted to:\n"<<Shift_in_X_direction<<" "<<Shift_in_Y_direction<<" "<<Shift_in_Z_direction<<std::endl;
    }
    void shift_node_velocities(void){
        for (int i=0; i<Num_of_Nodes; i++) {
            Node_Velocity[i][0]+=Shift_velocities_xyzVector[0];
            Node_Velocity[i][1]+=Shift_velocities_xyzVector[1];
            Node_Velocity[i][2]+=Shift_velocities_xyzVector[2];
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
    
    std::map<string, vector<string> > Params;
    vector<string> insertOrder;
    vector<string> values;
    Chromatin(){
        values.resize(2);
        
        values[0] ="path/to/my/coordinates.txt";
        values[1] ="#Path to a text file from which the X, Y, Z, Vx, Vy, and Vz of the chromatin nodes can be imported. Default path/to/my/coordinates.txt";
        Params["ImportCoordinates"] = values;
        insertOrder.push_back("ImportCoordinates");
        
        values[0] ="0";
        values[1] ="#Number of Chromatin nodes on a generated self avoiding random chain . If you are importing coordinates, this parameter will be ignored. Default 10";
        Params["GenerateRandomChain"] = values;
        insertOrder.push_back("GenerateRandomChain");
        
        values[0] ="false";
        values[1] ="#If \"true\" The programme will write the generated coordinates that can be imported later. Default false";
        Params["ExportGeneratedCoordinates"] = values;
        insertOrder.push_back("ExportGeneratedCoordinates");
        
        values[0] ="1";
        values[1] ="#Number of node types on the chain. For more than one node type, multiple LJsigma and LJepsillon can be defined to customise long range interactions. Default 10";
        Params["NodeTypes"] = values;
        insertOrder.push_back("NodeTypes");
        
        values[0] ="0";
        values[1] ="#Set the Lennard Jones 12-6 sigma. Default value 0. If the Chromatin is interacting with another class, the Sigma between them will be calculated as the average of the class's sigmas: Sigma {Chromatin & A} = 0.5(sigma{Chromatin}+Sigma{A}).  \n#If the chromatin has more than 1 node type, additional LJ sigmas should be provided, if not the default value (2.5*node radius) will be assigned to them. Example for 3 node types: LJsigma 1 3 2.5";
        Params["LJsigma"] = values;
        insertOrder.push_back("LJsigma");
        
        values[0] ="0";
        values[1] ="#Set the Lennard Jones 12-6 epsillon. Default value 0. If the Chromatin is interacting with another class, the Epsillon between them will be calculated as the geometrical average of the class's epsilons: Epsillon {Chromatin & A} = sqrt(epsillon{Chromatin}*+epsillon{A}). \n#If the chromatin has more than 1 node type, additional LJ epsiolnes should be provided, if not the default value (0) will be assigned to them. Example for 3 node types: LJepsilon 2 4 7";
        Params["LJepsilon"] = values;
        insertOrder.push_back("LJepsilon");
        
        values[0] ="1";
        values[1] ="#Mass asssigned to each node. Default value 1";
        Params["NodeMass"] = values;
        insertOrder.push_back("NodeMass");
        
        values[0] ="0";
        values[1] ="#Radius assigned to each node. If 0, half of the average bond distance of nodes will be used. The node radius is used to calculate the cutt-off and minimum energy distance for the 'Excluded Volume' and the 'Lennard-Jones' potential.";
        Params["NodeRadius"] = values;
        insertOrder.push_back("NodeRadius");
        
        values[0] ="H";
        values[1] ="#Set the bond potential. 'H' for harmonic. Default H.";
        Params["SpringModel"] = values;
        insertOrder.push_back("SpringModel");
        
        values[0] ="2000";
        values[1] ="#Set the bond potential rigidity coefficient. Default value 2000.";
        Params["SpringCoeff"] = values;
        insertOrder.push_back("SpringCoeff");
        
        values[0] ="0";
        values[1] ="#Set the damping coefficient for non harmonic potentials. Default value 0.";
        Params["DampingCoeff"] = values;
        insertOrder.push_back("DampingCoeff");
        
        values[0] ="0 0 0";
        values[1] ="#X, Y, Z components of a vector used to translate all coordinates befor beginning the simluation.";
        Params["CoordinateTranslateVector"] = values;
        insertOrder.push_back("CoordinateTranslateVector");
        
        values[0] ="0 0 0";
        values[1] ="#Vx, Vy, Vz components of a vector used to add to all initial node velocities befor beginning the simluation.";
        Params["VelocityShiftVector"] = values;
        insertOrder.push_back("VelocityShiftVector");
        
        values[0] ="1";
        values[1] ="#Used to scale the Chromatin coordinates. Default 1";
        Params["Scale"] = values;
        insertOrder.push_back("Scale");
        
        values[0] ="0";
        values[1] ="#If the chromatin node radius is smaller than the node distances along the chain, The VirtualBondLength (the chain bond length) will be used to put virtual chromatin nodes between the real nodes to mimic a cylinder. Default 0";
        Params["VirtualBondLength"] = values;
        insertOrder.push_back("VirtualBondLength");
        
        values[0] ="0";
        values[1] ="#If the VirtualBondLength is set, the VirtualBondRadius will be used to defing the radius of the virtual chromatin nodes. Default 0";
        Params["VirtualBondRadius"] = values;
        insertOrder.push_back("VirtualBondRadius");
        
        values[0] ="false";
        values[1] ="#If \"true\" The programme will find the best VirtualBondRadius that fits the parameters (NodeRadius and VirtualBondLength). Default 0";
        Params["OptimiseBondRadius"] = values;
        insertOrder.push_back("OptimiseBondRadius");
        
        
        
        
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

#endif // CHROMATIN_H

