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
#include <complex>

#include "OpenMM_structs.h"
#include "General_functions.hpp"
#include "OpenMM_structs.h"
#include "Arg_pars.hpp"

using std::vector;
using std::cout;
using std::endl;

class Membrane
{
    
protected:
    /**Store the label(pdb) used to write to the trajectory file */
    std::string label;
    /**Store the mem index the instance of the class has in the vector of Membranes in the main programme.*/
    int index;
    /*variables*/
    /**Store the volume of the membrane*/
    double volume;
    /**Store the surface area of the membrane*/
    double surface_area;
    
    
    /*variables*/
    /**Number of nodes in the membrane.*/
    int Num_of_Nodes=0;
    /*constants*/
    ///This is the number of nodes on the membrane (Both the outer membrane and the Nucleus). This is the first number that appears in the 'membrane' file (once opend with a text editor)
    int Num_of_Triangles; ///This is the number of triangles on the membrane (Both the outer membrane and the Nucleus). This is the number that appears in the 'membrane' file after the node position list is finished and before Gmesh lists the nodes that make a triangle.
    
    std::map<std::string, double> param_map;
    
    std::string Mesh_file_name="None";
    std::string resume_file_name="None";
    
    double sigma_LJ_12_6=0;
    double epsilon_LJ_12_6=0;
    
    double Node_Mass=1.0;//  also use in MD loop and should not be private unless we write some functions to get it outside the class
    double Total_Potential_Energy=0.0;
    double Total_Kinetic_Energy;
    double Radius=0;
    double Node_radius=1;
    double New_node_radius=-1;
    double Begin_update_time_in_Ps=0;
    double Update_nominal_length=-1;
    double End_update_time_in_Ps=0;
    double COM_velocity[3]={0};
    double COM_position[3]={0};
    
    double FENE_min=0, FENE_max=0, FENE_epsilon=0, FENE_k=0;
    
    int Num_of_Node_Pairs; //??? (This variable should be defined and explained)
    int Num_of_Triangle_Pairs;
    double X_in=0;
    double Y_in=0;
    double Z_in=0;
    double X_scale=0;
    double Y_scale=0;
    double Z_scale=0;
    
    double Average_Node_Distance();
    void read_gmesh_file (std::string gmesh_file);
    void read_ply_file (std::string ply_file);
    
    void Triangle_pair_counter ();
    void Normal_direction_Identifier();
    void FENE_log (void);
    void custom (void);
    void Node_Bonds_identifier(void);
    void Triangle_pair_identifier(void);
    void Bending_potetial_2(double theta_0);
    
    double Total_Bending_Energy = 0.0;
    
public:
    //Analysis funcs/vars:
    bool initial_random_rotation_coordinates = false;
    void rotate_coordinates(double theta, double phi);
    void rotate_particle_to_axes(ArgStruct args);
    void update_spherical_positions();
    void analysis_init(std::string Mesh_path);
    void calculate_dOmega(void);
    vector<double> node_dOmega;
    
    void load_pdb_frame(int frame, ArgStruct args);
    void import_pdb_frames(ArgStruct args, int file_index);
    void generate_ulm_mode(int ell, int m, double ulm, double radius);
    
    void generate_ulm_mode_real(int Ell, int M, double uLM, double radius);
    
    void myWritePDBFrame(int frameNum, std::string traj_name);
    std::complex<double> calc_vectorlist_vectorlist_surface_integral(vector<std::complex<double> > vectorlist1, vector<std::complex<double> > vectorlist2);
    std::complex<double> calc_vectorlist_vectorlist_surface_integral(vector<std::complex<double> > vectorlist1, vector<double> vectorlist2);
    double calc_vectorlist_vectorlist_surface_integral(vector<double> vectorlist1, vector<double> vectorlist2);
    vector<std::complex<double> > get_ylm_vectorlist_for_mesh(int ell, int m, bool complex_conjugate);
    vector<double> get_real_ylm_vectorlist_for_mesh(int ell, int m);
    vector<double> get_ulmYlm_vectorlist_for_mesh();
    
    
    std::complex<double> calc_complex_ylm_surface_integral(int ell, int m, double radius);
    /**return the  complex spherical harmonic for the provided parameters: Y_l,m (theta, phi). Where l is a positiv integer, m is defined -l <= m <= l, theta is 0 <= theta <= pi, and phi is defined 0 <= phi <= 2pi. s*/
    std::complex<double> calc_complex_ylmthetaphi(int l,  int  m, double theta, double phi);
    double calc_real_ylmthetaphi(int l,  int  m, double theta, double phi);
    
    
    
    void calculate_ulm(ArgStruct args);
    void calculate_real_ulm(ArgStruct args);
    void calculate_ulm_radiustest(int ell_max, int analysis_averaging_option);
    void calculate_ulm_radiustest_real(int ell_max, int analysis_averaging_option);
    void calculate_ulm_sub_particles(int ell_max, int analysis_averaging_option);
    void write_ulm(ArgStruct args, int file_index);
    
    vector<vector<double> > ulm_avg;
    vector<vector<double> > ulm_std;
  
    
    
    vector<vector<vector<double> > > pdb_frames;
    vector<vector<double> > spherical_positions;
    
    void set_com_to_zero(){
        update_COM_position();
        for (int i=0; i<Num_of_Nodes; i++) {
            for (int j=0; j<3; j++) {
                Node_Position[i][j] -= COM_position[j];
            }
        }
    }
    /**Returns the total bending energy  of the membrane (indepentant of OpenMM calculations) */
    double calculate_bending_energy(void);
    
    /**Returns the node radius update value (set value in the configuration file). */
    double get_new_node_radius(void){
        return New_node_radius;
    }
    /**Returns the spring nominal length update value (set value in the configuration file). */
    double get_new_nominal_length(void){
        return Update_nominal_length;
    }
    
    std::string output_file_neme;
    std::string file_time;
    
    vector<vector<double> > Node_Position;
    vector<vector<int> >    Triangle_list;
    //List that stores the IDs of triangles that are neighbours.
    vector<vector<int> >    Triangle_pair_list;
    //vector<vector<int> > Membrane_Node_Pair_list;
    vector<vector<int> >    Node_Bond_list;// this variable is  the same as Membrane_Node_pair_list. I think  the name "Membrane_Edges" is less confusing. and also we fill it in a different way.
    vector<vector<int> >    Triangle_Pair_Nodes;
    vector<vector<double> > Node_Velocity;// also update in MD loop and should not be private unless we write some functions to get it outside the class
    vector<vector<double> > Node_Force;// also update in MD loop and should not be private unless we write some functions to get it outside the class
    
    
    
    
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
    
    vector<double> node_voronoi_area;
    double surface_area_voronoi=0;
    
    vector<int> Bond_triangle_neighbour_indices ;
    
    void check(void);
    void check_radius_update_values(void);
    void Triangle_Pair_and_Node_Bonds_Identifier(); //I guess this will use in MD loop and thus it should define as a public membere of class.
    void Elastic_Force_Calculator(double theta_0);
    /**Construct a vector that lists all the neighbours of  each node and their respective index in the 'Node_Bond_list'.
     */
    void Node_neighbour_list_constructor();
    /**Construct a vector that lists the triangle neighbours of each node  pair.
     */
    void Bond_triangle_neighbour_list_constructor();
    
    void export_for_resume(int MD_step);
    
    //monte carlo flip functions
    bool check_monte_carlo=0;
    void find_the_new_neighbour(int neighbour_id[6], int previous_dihedral_index , int initial_pair, bool A_or_B);


//    bool monte_carlo_flip(MyOpenMMData* omm, Bonds* bonds, Dihedrals* dihedrals, MyAtomInfo atoms[], double& localDeltaE, int& Accepted_Try_Counter,int& pyramid_counter);

    void monte_carlo_flip(MyOpenMMData* omm, Bonds* bonds, Dihedrals* dihedrals, MyAtomInfo atoms[], double& localDeltaE, int& Accepted_Try_Counter,int& pyramid_counter, int &MC_total_tries, double &MC_Acceptance_Rate);

    double calculating_the_bond_energy(int index, bool initial_or_final, MyAtomInfo  atoms[],int number_of_privious_mem_nodes);
    double calculating_the_bond_energy_check(int p1, int p2, MyAtomInfo atoms[]);
    double calculating_the_bend_energy(int uncommn1, int common2, int common3, int uncommon4, bool initial_or_final, MyAtomInfo  atoms[], int number_of_privious_mem_nodes);
    double calculating_the_bend_energy_2(int uncommon1, int common2, int common3, int uncommon4, MyAtomInfo  atoms[], int number_of_privious_mem_nodes);
    double calculating_the_bond_length_check(int p1, int p2, MyAtomInfo atoms[]);
    void check_before_update(int triangle_A,int triangle_B, int new_neighbour_dihedrals[4][6],int& pyramid_counter, bool& accept);
    void update_Membrane_class_and_openmm(int initial_pair,int triangle_A,int triangle_B, int new_neighbour_dihedrals[4][6], MyOpenMMData* omm, Bonds* bonds, Dihedrals* dihedrals);
    void Update_Membrane(int initial_pair,int triangle_A,int triangle_B, int new_neighbour_dihedrals[4][6], int& bond_index);
    bool check_Pyramid(vector <int> A_neighbors, vector <int> B_neighbors);
    bool check_Pyramid_2(int initial_pair, vector<int> A_neighbours_dihedral_index, vector<int> B_neighbours_dihedral_index);
    void check_the_flip(MyOpenMMData* omm, Bonds* bonds, Dihedrals* dihedrals);
    
    //end of monte carlo flip functions
    
    void export_for_resume(int MD_step, MyAtomInfo atoms[], int atom_count);
    
    
    //    void initialise(string input_file_name , string Mesh_file_name);
    void initialise(std::string Mesh_file_name);
    //    void initialise(string Mesh_file_name, double x, double y, double z);
    void import(std::string import_file_name);
    void import_config(std::string config_file_name);
    void set_map_parameter(std::string param_name, double param_value);
    void generate_report(void);

    void write_pov_traj(std::string traj_name, std::string label, int currentstep);
    
    int spring_model=2;
    int mesh_format=1;// 1 represents gmsh generated mesh and 2 represents blender genereted mesh exported as a ply file.
    //    vector <double> T_Kinetic_Energy;
    double Spring_coefficient=10.0*GenConst::MD_T*GenConst::K; // streching constant
    double Bending_coefficient=20.0*GenConst::MD_T*GenConst::K; // bending constant
    double Damping_coefficient=0.0; // Viscosity of the Mmmbrane. It is applied in Force calculation for the Membrane Node pairs. I have commented out these parts in the 'Membrane_Force_Calculator' because I think the current code does not need it (some energy consuming array calculations were invloved).
    
    double Shift_in_X_direction=0.0; //???
    double Shift_in_Z_direction=0.0; //???
    double Shift_in_Y_direction=0.0; //???
    double x_speed=0.0; //???
    double y_speed=0.0;
    double z_speed=0.0;
    
    int ext_force_model=0;
    double kx=10;
    double ky=10;
    double kz=10;
    
    double com[3]; //center of mass
    double Min_node_pair_length, Max_node_pair_length, Average_node_pair_length;
    
    /**Returns the last saved volume*/
    double return_volume(void){
        return volume;
    }
    /**Calculate the voronoi area of  each  node and return a list that indicated the voronoi area of each node.*/
    vector<double> return_voronoi_node_area(){
        calculate_surface_area_with_voronoi();
        return node_voronoi_area;
    }
    /**Calculate the volume of a closed membrane by summing triangular pyramids.*/
    void calculate_volume_and_surface_area(void);
    
    /**Calculate the volume of a closed membrane by summing triangular pyramids.*/
    void calculate_surface_area_with_voronoi(void);
    
    /**Returns the last saved surface area*/
    double get_surface_area(void){
        return surface_area;
    }
    /**Return the new node radius 'end update time' in Ps.*/
    double get_End_update_time_in_Ps(void){
        return End_update_time_in_Ps;
    }
    /**Return the new node radius 'begin update time' in Ps.*/
    double get_Begin_update_time_in_Ps(void){
        return Begin_update_time_in_Ps;
    }
    
    double rescale_factor=1;
    /** Assigns the label(pdb) used to write to the trajectory files. */
    void set_label(std::string lab){
        label=lab;
    }
    /** return the label(pdb) used to write to the trajectory files. */
    std::string get_label(void){
        return label;
    }
//    void Relax_1(void);
    
    
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
    /**Return damp coefficient. */
    double get_damping_coefficient(void){
        return Damping_coefficient;
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
    /**Return the Lenard Jones 12 6 sigma, used to setup the openmm system for the LJ interaction.*/
    double get_sigma_LJ_12_6(void){
        return sigma_LJ_12_6;
    }
    /**Return the Lenard Jones 12 6 sigma, used to setup the openmm system for the LJ interaction.*/
    double get_epsilon_LJ_12_6(void){
        return epsilon_LJ_12_6;
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
    void set_FENE_param_2(double &lmin, double &lmax, double &epsilon, double &k){
        lmin = FENE_min;
        lmax = FENE_max;
        epsilon = FENE_epsilon;
        k = FENE_k;
    }
    /**Calculate and return the cot of the  angle between nodes A(middle of the angle), B, and C.*/
    double calc_theta_angle_ABC(int node_A, int node_B, int node_C);
    
    /**Set FENE calculated parameters.*/
    void set_FENE_param(double &le0, double &le1, double &lmin, double &lmax){
        
        
        double width= 0.66*Average_node_pair_length; // it defines a minimum limit for weil width.
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
    /**Returns the r (0), theta (1), and phi (2) spherical coordinate of the node index (number). Note: update_spherical_positions(); should be called to upddate the values.*/
    double get_spherical_position(int node_number, int node_coordinate){
        return spherical_positions[node_number][node_coordinate];
    }
    /**Returns the x (0), y (1), and z (2) coordinate of the node index (number).*/
    double get_node_position(int node_number, int node_coordinate){
        return Node_Position[node_number][node_coordinate];
    }
    /**Returns the x (0), y (1), and z (2) velocities of the node index (number).*/
    double get_node_velocity(int node_number, int node_coordinate){
        return Node_Velocity[node_number][node_coordinate];
    }
    /**Set the current state (OpenMM) of the class.*/
    void set_state(MyAtomInfo all_atoms[], int atom_count);
    
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
    /** This function shifts the whole membrane.*/
    void shift_position (double x, double y, double z) {
        for (int i=0; i<Num_of_Nodes; i++) {
            Node_Position[i][0]+=x;
            Node_Position[i][1]+=y;
            Node_Position[i][2]+=z;
        }
        X_in+=x;
        Y_in+=y;
        Z_in+=z;
    }
    
};

#endif // MEMBRANE_H
