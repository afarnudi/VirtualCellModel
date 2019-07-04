/**
 * \class Membrane Membrane.h Include/Membrane.h
 *
 *
 * \brief The class 'Membrane' manages the required enviroment (system, forces, etc) of coordinates on a triangular mesh.
 *
 * The Membrane class deciphers the mesh coordinates from the input mesh file provided in the configuration file. The bond information is also extracted from the mesh file and stored in the class.
 * The various interactions of membrane nodes and all the membrane node info is also stored in this class.
 *
 * \note Work in progress.
 *
 *
 *
 *
 * Contact: a.farnudi@gmail.com
 *
 */

#ifndef MEMBRANE_H
#define MEMBRANE_H
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <map>
#include <iomanip>
#include <iterator>


#include "OpenMM_structs.h"
#include "General_functions.hpp"

using std::vector;
using std::cout;
using std::endl;

class Membrane
{

private:
    /** Store the label(pdb) used to write to the trajectory file */

    std::string label;
    
//    /** Convert the node position data to OpenMM format*/
//    void convert_position_to_openmm(void);
//
//    /** Convert the bond info to OpenMM format*/
//    void convert_bond_info_to_openmm(void);
//
//    //OpenMM data structures.
//    MyAtomInfo* myatominfo;
//    Bonds* bonds;
//
//    /** This is our opaque "handle" class containing all the OpenMM objects that
//     * must persist from call to call during a simulation. The main programme gets
//     * a pointer to one of these but sees it as essentially a void* since it
//     * doesn't know the definition of this class.
//     */
//    struct MyOpenMMData {
//        MyOpenMMData() : system(0), context(0), integrator(0) {}
//        ~MyOpenMMData() {delete context; delete integrator; delete system;}
//        OpenMM::System*         system;
//        OpenMM::Integrator*     integrator;
//        OpenMM::Context*  context;
//    };
//    //OpenMM platform
//    std::string   platformName;
//    // Allocate space to hold OpenMM objects while we're using them.
//    MyOpenMMData* omm = new MyOpenMMData();
//
//    /** -----------------------------------------------------------------------------
//     *                      INITIALIZE OpenMM DATA STRUCTURES
//     * -----------------------------------------------------------------------------
//     * We take these actions here:
//     * (1) Load any available OpenMM plugins, e.g. Cuda and Brook.
//     * (2) Allocate a MyOpenMMData structure to hang on to OpenMM data structures
//     *     in a manner which is opaque to the caller.
//     * (3) Fill the OpenMM::System with the force field parameters we want to
//     *     use and the particular set of atoms to be simulated.
//     * (4) Create an Integrator and a Context associating the Integrator with
//     *     the System.
//     * (5) Select the OpenMM platform to be used.
//     * (6) Return the MyOpenMMData struct and the name of the Platform in use.
//     *
//     * Note that this function must understand the calling MD code's molecule and
//     * force field data structures so will need to be customized for each MD code.
//     */
//     MyOpenMMData* myInitializeOpenMM();
    
    
    int index;
    /*variables*/
    int Num_of_Nodes;
    /*constants*/
    ///This is the number of nodes on the membrane (Both the outer membrane and the Nucleus). This is the first number that appears in the 'membrane' file (once opend with a text editor)
    int Num_of_Triangles; ///This is the number of triangles on the membrane (Both the outer membrane and the Nucleus). This is the number that appears in the 'membrane' file after the node position list is finished and before Gmesh lists the nodes that make a triangle.
    int MD_num_of_Relaxation_steps=200000;
    int MD_correction_steps=2000;

    std::map<std::string, double> param_map;
    
    std::string Mesh_file_name="None";
    std::string resume_file_name="None";
    
    bool Relaxation=false;
    bool Shift= false;
    bool Relax_with_actin=false;
    int Relaxation_Process_Model=1; /// 1 represents the relaxation Processes without node correction and 2 includes node corrections.
    int correction_progress;
    double ECM_interaction_cut_off=0;
    double vesicle_interaction_cut_off=10;
    double Node_Mass=1.0;//  also use in MD loop and should not be private unless we write some functions to get it outside the class
    double Total_Potential_Energy;
    double Total_Kinetic_Energy;
    double Radius=0;
    double Node_radius=1;
    double COM_velocity[3]={0};
    double COM_position[3]={0};

    int Num_of_Node_Pairs; //??? (This variable should be defined and explained)
    int Num_of_Triangle_Pairs;
    double X_in=0;
    double Y_in=0;
    double Z_in=0;
    double X_scale=0;
    double Y_scale=0;
    double Z_scale=0;

    double ECM_interaction_strength=1;


    double vesicle_interaction_strength=30;
    double vesicle_interaction_sigma=3;
    double Average_Node_Distance();
    void read_gmesh_file (std::string gmesh_file);
    void read_ply_file (std::string ply_file);
    void read_membrabe_input(std::string input_file);
    void Triangle_pair_counter ();
    void Normal_direction_Identifier();
    void node_distance_correction(void);
    //    void Normal_direction_Identifier(double x, double y, double z);

    
    void FENE_log (void);
    void Hookian (void);
    void custom (void);
    void Relaxation_potential(void);
    void Node_Bonds_identifier(void);
    void Triangle_pair_identifier(void);
    void Bending_potetial(void);
    //    void Bending_potetial_2(void);
    void Bending_potetial_2(double theta_0);

    void export_relaxed(int MD_step);
    int membrane_counter;
    bool rescale=0;
    
    
public:
    ///Call all initilisation members and initilise openmm handles.
    void initilise_openmm(void);
    
    ///node-node hard sphere interaction.
    void excluded_volume(void);

    double min_radius_after_relaxation;

    std::string output_file_neme;
    std::string file_time;
    bool particle_type=0; // True means the object is a particle and false means it is a vesicle
    
    vector<vector<double> >Node_Position;
    vector<vector<int> > Triangle_list;
    //List that stores the IDs of triangles that are neighbours.
    vector<vector<int> > Triangle_pair_list;
    //vector<vector<int> > Membrane_Node_Pair_list;
    vector<vector<int> > Node_Bond_list;// this variable is  the same as Membrane_Node_pair_list. I think  the name "Membrane_Edges" is less confusing. and also we fill it in a different way.
    vector<vector<int> > Triangle_Pair_Nodes;
    vector<vector<double> > Node_Velocity;// also update in MD loop and should not be private unless we write some functions to get it outside the class
    vector<vector<double> > Node_Force;// also update in MD loop and should not be private unless we write some functions to get it outside the class
    vector<vector<int> > Node_neighbour_list;
    //    vector<double>DamperCheck;
    vector<double>SinusCheck;
    void Damper_check(int MD_step);
    void check(void);
    vector<vector<int> > ECM_Node_neighbour_list;
    vector<vector<int> > Vesicle_Node_neighbour_list;
    void Triangle_Pair_and_Node_Bonds_Identifier(); //I guess this will use in MD loop and thus it should define as a public membere of class.
    //int Membrane_num_of_Node_Pair_Counter();// Hoda: no need to this function after modifying Membrane_Triangle_Pair_and_Edges_Identifier
    //void Membrane_num_of_Node_Pair_Counter_2();//Hoda: no need to this function after modifying Membrane_Triangle_Pair_and_Edges_Identifier
    void Elastic_Force_Calculator(double theta_0);
    void MD_Evolution_beginning (double MD_Time_Step);
    void MD_Evolution_end (double MD_Time_Step);
    void ConstantSurfaceForceLocalTriangles ();
    void Node_neighbour_list_constructor();
    void export_for_resume(int MD_step);

    //    void initialise(string input_file_name , string Mesh_file_name);
    void initialise(std::string Mesh_file_name);
    //    void initialise(string Mesh_file_name, double x, double y, double z);
    void import(std::string import_file_name);
    void import_config(std::string config_file_name);
    void set_map_parameter(std::string param_name, double param_value);
    void generate_report(void);
    void export_velocities(int MD_step);
    void Thermostat_2(double MD_KT);
    void Thermostat_N6(double MD_KT);
    void Thermostat_Bussi(double MD_KT);
    void calculate_mesh_properties(void);
    void relaxation_traj (void);
    void write_traj (std::string traj_name);
    void write_traj (std::string traj_name, std::string label);
    void write_parameters(int MD_Step);
    void omega_calculator(void);
    void omega_calculator_2(void);
    void equilibrate (void);
    void write_pov_traj(std::string traj_name, std::string label, int currentstep);
    double Average_velocity_squared();
    double Omega[3]={0};



    int **Normal_direction; //??? (These 2 elements should be defined and explained)
    int spring_model=2;
    int mesh_format=1;// 1 represents gmsh generated mesh and 2 represents blender genereted mesh exported as a ply file.
    //    vector <double> T_Kinetic_Energy;
    double Total_potential_Energy=0.0;
    double Spring_coefficient=10.0*GenConst::MD_T*GenConst::K; // streching constant
    double Bending_coefficient=20.0*GenConst::MD_T*GenConst::K; // bending constant
    double Damping_coefficient=0.0; // Viscosity of the Mmmbrane. It is applied in Force calculation for the Membrane Node pairs. I have commented out these parts in the 'Membrane_Force_Calculator' because I think the current code does not need it (some energy consuming array calculations were invloved).
    double K_surfaceConstant_local=100.0;
    double Spring_force_cutt_off=1000.0;
    double Shift_in_X_direction=0.0; //???
    double Shift_in_Z_direction=0.0; //???
    double Shift_in_Y_direction=0.0; //???
    double Downward_speed=0.0; //???
    //bool =0;
    double com[3]; //center of mass
    double Min_node_pair_length, Max_node_pair_length, Average_node_pair_length;

    bool on_or_off_Spring_force_cutt_off=0; //??? I add it myself because virus should not have cut off


private:

    int mem_index;
    /*variables*/
    //    void Normal_direction_Identifier(double x, double y, double z);
    void Rescale(double rescale_factor);
    void potential_1 (void);
    void potential_2 (void);

    /** This function shifts the whole membrane.*/
    void shift_position (double x , double y, double z);  
public:
    double rescale_factor;
    /** Assigns the label(pdb) used to write to the trajectory files. */
    void set_label(std::string lab){
        label=lab;
    }
    /** return the label(pdb) used to write to the trajectory files. */
    std::string return_label(void){
        return label;
    }
    void Relax_1(void);
    void Relax_2(void);

    /** \brief public access to total number of membrane nodes.
     * \return integer number of nodes in the membrane
     */
    int get_num_of_nodes(void){
        return Num_of_Nodes;
    }
    /**Return the number of bonds between membrane nodes.*/
    int get_num_of_node_pairs(void){
        return Num_of_Node_Pairs;
    }
    /**Return the id(int) of the nodes in the bonds. The id of each node in the bond list is stored in the first (0) and the second (1) id slot.*/
    int get_node_pair(int bond_num, int node_id){
        return Node_Bond_list[bond_num][node_id];
    }
    /**Return the node mass. At the current stage of this code all membrane nodes have the same mass. */
    double get_node_mass(void){
        return Node_Mass;
    }
    /**Return the average distance of the nodes as calculated from the mesh. */
    double get_avg_node_dist(void){
        return Average_node_pair_length;
    }
    /**Return spring stiffness coefficient. */
    double get_spring_stiffness_coefficient(void){
        return Spring_coefficient;
    }
    /**Return bending stiffness coefficient. */
    double get_bending_stiffness_coefficient(void){
        return Bending_coefficient;
    }
    /**Returns the calculated number of triangles in the imported mesh file.*/
    int get_num_of_triangle_pairs(){
        return int(Triangle_Pair_Nodes.size());
    }
    /**Return the node IDs of the dihedral angles.*/
    vector<int> get_traingle_pair_nodes_list(int triangle_pair){
        return Triangle_Pair_Nodes[triangle_pair];
    }
    /**Return the node ID of the dihedral angles member.*/
    int get_traingle_pair_node(int triangle_pair, int index){
        return Triangle_Pair_Nodes[triangle_pair][index];
    }
    /**Return input spring model, used to setup the openmm system for the bonds.*/
    int get_spring_model(void){
        return spring_model;
    }
    
    
    
    /**Set FENE calculated parameters.*/
    void set_FENE_param(double &le0, double &le1, double &lmin, double &lmax){
        lmax=Max_node_pair_length*1.05;
        lmin=Min_node_pair_length*0.95;
        le0=lmin+3*(lmax-lmin)/4;
        le1=lmin+(lmax-lmin)/4;
    }
    
    int  get_Relaxation_Process_Model(void){
        return Relaxation_Process_Model;
    }
    bool  get_Relaxation_flag(void){
        return Relaxation;
    }



    void shift_velocity (double vx, double vy, double vz){
        for (int i=0; i<Num_of_Nodes; i++) {
            Node_Velocity[i][0]+=vx;
            Node_Velocity[i][1]+=vy;
            Node_Velocity[i][2]+=vz;
        }
    }
    /**Returns the calculated number of triangles in the imported mesh file.*/
    int get_num_of_triangle(){
        return Num_of_Triangles;
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
        cout<<"\n\naverage_force_x="<<average_force_x/Num_of_Nodes<<"\naverage_force_y="<<average_force_y/Num_of_Nodes<<"\naverage_force_z="<<average_force_z/Num_of_Nodes<<endl;
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
    double get_ECM_interaction_cut_off(void){
        return ECM_interaction_cut_off;
    }
    double get_ECM_interaction_strength(void){
        return ECM_interaction_strength;
    }
    double get_vesicle_interaction_cut_off(void){
        return vesicle_interaction_cut_off;
    } 
    double get_vesicle_interaction_sigma(void){
        return vesicle_interaction_sigma;
    }
    double get_vesicle_interaction_strength(void){
        return vesicle_interaction_strength;
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
    double get_min_radius_after_relaxation(void){
        return min_radius_after_relaxation;
    }
    bool get_relax_with_actin_flag(void){
        return Relax_with_actin;
    }
    int get_correction_progress(void){
        return correction_progress;}
    bool bending_coefficient_status(void){
        if (Bending_coefficient !=0) {
            return true;
        } else {
            return false;
        }
    }

    
};

#endif // MEMBRANE_H
