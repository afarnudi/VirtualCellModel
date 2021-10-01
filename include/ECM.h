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
#include "OpenMM_structs.h"

using std::vector;
using std::cout;
using std::endl;
using std::string;

class ECM {
private:
    
    std::string label;
    std::string Mesh_file_name;
    std::string MeshType;
    
    std::string file_time;
    int index;
    
    std::map<std::string, double> param_map;
    
    double Node_Mass=1;
    string Node_radius_stat;
    //double Node_radius=1;
    vector<double> Node_radius;
    int spring_model=1;
    double Spring_coefficient=10;
    double stiffness_gradient_x=0;
    double stiffness_gradient_y=0;
    double stiffness_gradient_z=0;
    
    int receptor_type = 1;
    double receptor_density = 0.5;
    double receptor_center_x = 0;
    double receptor_center_y = 0;
    double receptor_center_z = 0;
    double receptor_gradient_x = 0;
    double receptor_gradient_y = 0;
    double receptor_gradient_z = 0;
    
    vector<double> Shift_position_xyzVector;
    vector<double> Shift_velocities_xyzVector;
//    double Shift_in_X_direction=0;
//    double Shift_in_Y_direction=0;
//    double Shift_in_Z_direction=0;
    
    string mesh_format;
    
    
    int Num_of_Nodes=0;
    int Num_of_Triangle_Pairs=0;
    int Num_of_Node_Pairs=0;
    int Num_of_Triangles=0;
    
   // double x_speed=0.0; //???
   // double y_speed=0.0;
   // double z_speed=0.0;
    
    int ext_force_model=0;
    vector<double> ext_force_rigidity ;
    //double kx=10;
    //double ky=10;
    //double kz=10;
    
    double Kelvin_Damping_Coefficient=100;
    double Dashpot_Viscosity=0.02;
    
    vector<vector<double> > Node_Force;
    vector<vector<double> > Node_Position;
    vector<vector<int> > Triangle_List;
    vector<vector<int> > Pyramid_Nodes;
    vector<vector<int> > square_Nodes;
    
    double interaction_range=1.0;
    double epsilon=0.6;
    double sigma=15.0;
    double COM_velocity[3]={0};
    double COM_position[3]={0};
    
    double Min_node_pair_length, Max_node_pair_length, Average_node_pair_length;
    double rescale_factor;
    
    double sigma_LJ_12_6=0;
    double epsilon_LJ_12_6=0;
    
    vector<vector<double> > Node_Velocity;
    vector<vector<int> > Node_Bond_list;
    vector<vector<int> > Node_neighbour_list;
    string  Node_Bond_Nominal_Length_stat;
    double  Node_Bond_user_defined_Nominal_Length_in_Nm;
    void Node_Bond_identifier_2D(void);
    void Node_Bond_identifier_3D(void);
    void Node_Bond_identifier_3D_square(void);
    
    void initialise(int dimension);
    void initialise(std::string Mesh_file_name);
    void Node_neighbour_list_constructor(void);
    
//    void read_input(string input_file);
    void read_gmesh_file_2D (std::string gmesh_file);
    void read_gmesh_file_3D (std::string gmesh_file);
    void read_gmesh_file_3D_square (std::string gmesh_file);
    void normal_direction_Identifier (double x, double y, double z);
    void normal_direction_Identifier (void);
    void check(void);
    
    
    
public:
    
    void import_config(vector<string> configlines);
    void export_for_resume(int MD_step);
    void export_for_resume(int MD_step, MyAtomInfo atoms[], int atom_count);
    
    /**Set the current state (OpenMM) of the class.*/
    void set_state(MyAtomInfo all_atoms[], int atom_count);
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
    /**Return spring stiffness coefficient. */
    double get_spring_stiffness_coefficient(void){
        return Spring_coefficient;
    }
    
    /**Return spring stiffness gradient_x */
    double get_stiffness_gradient_x(void){
        return stiffness_gradient_x;
    }
    
    /**Return spring stiffness gradient_y */
    double get_stiffness_gradient_y(void){
        return stiffness_gradient_y;
    }
    
    /**Return spring stiffness gradient_z */
    double get_stiffness_gradient_z(void){
        return stiffness_gradient_z;
    }
    
    /**Return receptor type. */
    int get_receptor_type(void){
        return receptor_type;
    }
    
    /**Return receptor density. */
    double get_receptor_density(void){
        return receptor_density;
    }
    
    /**Return receptor center_x. */
    double get_receptor_center_x(void){
        return receptor_center_x;
    }
    
    /**Return receptor center_y. */
    double get_receptor_center_y(void){
        return receptor_center_y;
    }
    
    /**Return receptor center_z. */
    double get_receptor_center_z(void){
        return receptor_center_z;
    }
    
    /**Return receptor density gradient_x */
    double get_receptor_gradient_x(void){
        return receptor_gradient_x;
    }
    
    /**Return receptor density gradient_y */
    double get_receptor_gradient_y(void){
        return receptor_gradient_y;
    }
    
    /**Return receptor density gradient_z */
    double get_receptor_gradient_z(void){
        return receptor_gradient_z;
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
    /**Return the Lenard Jones 12 6 sigma, used to setup the openmm system for the LJ interaction.*/
    double get_sigma_LJ_12_6(void){
        return sigma_LJ_12_6;
    }
    /**Return the Lenard Jones 12 6 sigma, used to setup the openmm system for the LJ interaction.*/
    double get_epsilon_LJ_12_6(void){
        return epsilon_LJ_12_6;
    }
    
//    double get_node_radius(void){
//        return Node_radius;
//    }
    
    double get_node_radius(int index){
           return Node_radius[index];
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
    
    double get_node_position(int node_number, int node_coordinate){
        return Node_Position[node_number][node_coordinate];
    }
    void add_to_force(double force,int index, int coor){
        Node_Force[index][coor]+=force;
    }
//    void shift_node_positions(void){
//        for (int i=0; i<Num_of_Nodes; i++) {
//            Node_Position[i][0]+=Shift_in_X_direction;
//            Node_Position[i][1]+=Shift_in_Y_direction;
//            Node_Position[i][2]+=Shift_in_Z_direction;
//        }
//    }
    
    /** This function shifts the whole ecm.*/
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
    /**Set FENE calculated parameters.*/
    void set_FENE_param(double &le0, double &le1, double &lmin, double &lmax){
        
        double width= Average_node_pair_length; // it defines a minimum limit for weil width.
        double delta_length=Max_node_pair_length-Min_node_pair_length;
        double width_scaling =0.2; // this variable adjusts the log_barrier potential width. if the edges have almost the same length, then it shlould be  redefined to have an appropriate weil width
        if ((1+ 2*width_scaling)*delta_length < width ){ //this condition shows the case when the mesh is almost ordered.
            width_scaling= (width-delta_length)/(2*delta_length); //in this case the width tuned in a way that the weil width becomes exactly 0.66*Average_node_pair_lenght
            //            cout<<"\n\nif == true\n\n";
        }
        else{
            width=  (1+ 2*width_scaling)*delta_length; //for disordered meshes it sets the weil witdh by scaling factor 0.2. in this case, the width may become larger than 0.66Average which was the minimum limit.
            //            cout<<"\n\nif == false\n\n";
        }
        lmin= Min_node_pair_length - width_scaling*delta_length;
        lmax= Max_node_pair_length +  width_scaling*delta_length;
        
        le0 = lmin + 3.0*(lmax-lmin)/4.0;
        le1 = lmin +   (lmax-lmin)/4.0;
        
        //        lmax=Max_node_pair_length*1.05;
        //        lmin=Min_node_pair_length*0.95;
        //        le0=lmin+3*(lmax-lmin)/4;
        //        le1=lmin+(lmax-lmin)/4;
    }
    
    std::map<string, vector<string> > Params;
    vector<string> insertOrder;
    vector<string> values;
    ECM(){
        values.resize(2);
        
      values[0] ="Path/to/my/meshfile.extension";
        values[1] ="#Path to the mesh file. Supported formats: Blender's ply, and Gmsh 2. The ECM class cannot be initilised without a mesh file.";
        Params["MeshFile"] = values;
        insertOrder.push_back("MeshFile");
        
        values[0] ="2d";
        values[1] ="#2d or 3d ecm.  default value 2d";
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
        
        values[0] ="Au";
        values[1] ="#Set the rest length of the mesh springs using: 1) Au: The initial bond lengths of the mesh; 2) Av: The average node pair lengths; 3) \"value\": Where you type a specific  value.";
        Params["NominalLengthInNm"] = values;
        insertOrder.push_back("NominalLengthInNm");
        
        values[0] ="H";
        values[1] ="#Set the bond potential. 'H' for harmonic. 'FENE' for a finitely extensible nonlinear elastic model. 'kelvin' for the Kelvin-Voigt potential. 'N' for no potential. Default H";
        Params["SpringModel"] = values;
        insertOrder.push_back("SpringModel");
        
        values[0] ="0";
        values[1] ="#Set the bond potential rigidity coefficient. Default value 0.";
        Params["SpringCoeff"] = values;
        insertOrder.push_back("SpringCoeff");
        
        values[0] ="0";
        values[1] ="#Set the damping coefficient. Default value 0.";
        Params["DampingCoeff"] = values;
        insertOrder.push_back("DampingCoeff");
        
        values[0] ="0";
        values[1] ="#Under development. Do not use this flag.";
        Params["ExtForceModel"] = values;
        insertOrder.push_back("ExtForceModel");
        
        values[0] ="0 0 0";
        values[1] ="#Under development. Do not use this flag.";
        Params["ExtForceRigidity"] = values;
        insertOrder.push_back("ExtForceRigidity");
        
        
        values[0] ="1";
        values[1] ="#Used to scale the ECM coordinates. Default 1";
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
        
        values[0] ="0";
        values[1] ="#Set the Lennard Jones 12-6 epsillon. Default value 0. If the ECM is interacting with another class, the Epsillon between them will be calculated as the geometrical average of the class's epsilons: Epsillon {ECM & A} = sqrt(epsillon{ECM}*+epsillon{A}).";
        Params["LJepsilon"] = values;
        insertOrder.push_back("LJepsilon");
        
//        values[0] ="10";
//        values[1] ="#Set the Lennard Jones 12-6 epsillon. Default value 0. If the ECM is interacting with another class, the Epsillon between them will be calculated as the geometrical average of the class's epsilons: Epsillon {ECM & A} = sqrt(epsillon{ECM}*+epsillon{A}).";
//        Params["LJepsilon"] = values;
//        insertOrder.push_back("LJsigma");
        
        
        values[0] ="0";
        values[1] ="#Set gradient in ecm rigidity in x direction. Default 0";
        Params["stiffness_gradient_x"] = values;
        insertOrder.push_back("stiffness_gradient_x");
        
        values[0] ="0";
        values[1] ="#Set gradient in ecm rigidity in y direction. Default 0";
        Params["stiffness_gradient_y"] = values;
        insertOrder.push_back("stiffness_gradient_y");
        
        values[0] ="0";
        values[1] ="#Set gradient in ecm rigidity in z direction. Default 0";
        Params["stiffness_gradient_z"] = values;
        insertOrder.push_back("stiffness_gradient_z");
        
        
        values[0] ="1";
        values[1] ="#Set receptor type on ecm. '1' for normal ecm with receptor density defined by 'receptor_density' parameter, 2 for stripe pattern receptors. Default 1";
        Params["receptor_type"] = values;
        insertOrder.push_back("receptor_type");
        
        
        values[0] ="1";
        values[1] ="#Set receptor surface density on ECM. Default 1";
        Params["receptor_density"] = values;
        insertOrder.push_back("receptor_density");
        
        
        values[0] ="0";
        values[1] ="#Set gradient in ecm receptor density in x direction. Default 0";
        Params["receptor_gradient_x"] = values;
        insertOrder.push_back("receptor_gradient_x");
        
        values[0] ="0";
        values[1] ="#Set gradient in ecm receptor density in y direction. Default 0";
        Params["receptor_gradient_y"] = values;
        insertOrder.push_back("receptor_gradient_y");
        
        values[0] ="0";
        values[1] ="#Set gradient in ecm receptor density in z direction. Default 0";
        Params["receptor_gradient_z"] = values;
        insertOrder.push_back("receptor_gradient_z");
        
        
        values[0] ="0";
        values[1] ="#will be added . Default 0";
        Params["receptor_center_x"] = values;
        insertOrder.push_back("receptor_center_x");
        
        values[0] ="0";
        values[1] ="#will be added . Default 0";
        Params["receptor_center_y"] = values;
        insertOrder.push_back("receptor_center_y");
        
        values[0] ="0";
        values[1] ="#will be added . Default 0";
        Params["receptor_center_z"] = values;
        insertOrder.push_back("receptor_center_z");
        
        
        
        
        
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

#endif /* ECM_hpp */
