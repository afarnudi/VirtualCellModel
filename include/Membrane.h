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
using std::string;

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
    string Node_radius_stat;
    vector<double> Node_radius;
    double New_Radius=-1;
    double New_node_radius=-1;
    double Begin_update_time_in_Ps=0;
    double End_update_time_in_Ps=0;
    bool   update_radius_stat=false;
    
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
    vector<vector< int > > dihedral_atoms;
    
    /**write geometrical properties if it is true.*/
    bool WantGeometricProps=false;
    
public:
    //Analysis funcs/vars:
    bool initial_random_rotation_coordinates = false;
    void rotate_coordinates(double theta, double phi);
    void rotate_particle_to_axes(ArgStruct_Analysis args);
    void update_spherical_positions();
    void convert_spherical_positions_to_cartisian();
    void analysis_init(std::string Mesh_path);
    void calculate_dOmega(void);
    vector<double> node_dOmega;
    
    void load_pdb_frame(int frame, ArgStruct_Analysis args);
    void import_pdb_frames(ArgStruct_Analysis args, int file_index);
    void generate_ulm_mode(int ell, int m, double ulm, double radius);
    
    void generate_ulm_mode_real(int Ell, int M, double uLM, double radius);
    void add_ulm_mode_real(int Ell, int M, double uLM, double radius);
    
    void myWritePDBFrame(int frameNum, std::string traj_name);
    std::complex<double> calc_vectorlist_vectorlist_surface_integral(vector<std::complex<double> > vectorlist1, vector<std::complex<double> > vectorlist2);
    std::complex<double> calc_vectorlist_vectorlist_surface_integral(vector<std::complex<double> > vectorlist1, vector<double> vectorlist2);
    double calc_vectorlist_vectorlist_surface_integral(vector<double> vectorlist1, vector<double> vectorlist2);
    vector<std::complex<double> > get_ylm_vectorlist_for_mesh(int ell, int m, bool complex_conjugate);
    vector<double> get_real_ylm_vectorlist_for_mesh(int ell, int m);
    vector<double> get_ulmYlm_vectorlist_for_mesh();
    vector<double> get_ulmYlm_vectorlist_for_mesh(char Requiv);
    
    
    std::complex<double> calc_complex_ylm_surface_integral(int ell, int m, double radius);
    /**return the  complex spherical harmonic for the provided parameters: Y_l,m (theta, phi). Where l is a positiv integer, m is defined -l <= m <= l, theta is 0 <= theta <= pi, and phi is defined 0 <= phi <= 2pi. s*/
    std::complex<double> calc_complex_ylmthetaphi(int l,  int  m, double theta, double phi);
    double calc_real_ylmthetaphi(int l,  int  m, double theta, double phi);
    
    
    
    void calculate_ulm(ArgStruct_Analysis args);
    void calculate_real_ulm(ArgStruct_Analysis args);
    void calculate_real_ulm(ArgStruct_Analysis args, char Requiv, bool clear);
    void calculate_ulm_radiustest(int ell_max, int analysis_averaging_option);
    void calculate_ulm_radiustest_real(int ell_max, int analysis_averaging_option);
    void calculate_ulm_sub_particles(int ell_max, int analysis_averaging_option);
    void write_ulm(ArgStruct_Analysis args, int file_index);
    
    vector<vector<double> > ulm_avg;
    vector<vector<double> > ulm_std;
    vector<vector<double> > ulm_temp_for_analysis;
  
    
    
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
    /**Return the geometrical properties write status*/
    bool get_GeometricProps_flag(){
        return WantGeometricProps;
    }
    
    void write_geometrics();
    
    /**Returns the total bending energy  of the membrane (indepentant of OpenMM calculations) */
    double calculate_bending_energy(void);
    
    /**Returns the node radius update value (set value in the configuration file). */
    double get_new_radius(void){
        return New_Radius;
    }
    /**Returns whether the Membrane is setup for radius change (true) or not (false). */
    bool get_update_status(void){
        return update_radius_stat;
    }
    
    std::string output_file_neme;
    std::string file_time;
    
    vector<vector<double> > Node_Position;
    vector<vector<int> >    Triangle_list;
    //List that stores the IDs of triangles that are neighbours.
    vector<vector<int> >    Triangle_pair_list;
    //vector<vector<int> > Membrane_Node_Pair_list;
    vector<vector<int> >    Node_Bond_list;// this variable is  the same as Membrane_Node_pair_list. I think  the name "Membrane_Edges" is less confusing. and also we fill it in a different way.
    vector<double>          Node_Bond_distances_in_Nm;
    string                  Node_Bond_distances_stat;
    double                  Node_Bond_Nominal_Length_in_Nm=-1;
    
    vector<double>          Triangle_pair_angles_in_radians;
    string                  Triangle_pair_angle_stat;
    double                  Triangle_pair_Nominal_angle_in_degrees=-1;
    
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
    vector<vector<double> > node_voronoi_normal_vec;
    vector<double> node_voronoi_area;
    double surface_area_voronoi=0;
    
    vector<int> Bond_triangle_neighbour_indices ;
    
    void check(void);
    void set_bond_nominal_length(void);
    void set_node_radius(void);
    void set_bending_nominal_angle(void);
    void set_dihedral_atoms(void);
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
    void import_config(vector<string> configlines);

    void write_pov_traj(std::string traj_name, std::string label, int currentstep);
    
    int spring_model=2;
    int mesh_format=1;// 1 represents gmsh generated mesh and 2 represents blender genereted mesh exported as a ply file.
    //    vector <double> T_Kinetic_Energy;
    double Spring_coefficient=0.;
    double Bending_coefficient=0.;
    double SpontaneousTriangleBendingInDegrees=0.;
    double Damping_coefficient=0.0; // Viscosity of the Mmmbrane. It is applied in Force calculation for the Membrane Node pairs. I have commented out these parts in the 'Membrane_Force_Calculator' because I think the current code does not need it (some energy consuming array calculations were invloved).
    
    vector<double> Shift_position_xyzVector;
    vector<double> Shift_velocities_xyzVector;
//    double Shift_in_X_direction=0.0; //???
//    double Shift_in_Z_direction=0.0; //???
//    double Shift_in_Y_direction=0.0; //???
//    double x_speed=0.0; //???
//    double y_speed=0.0;
//    double z_speed=0.0;
    
    int ext_force_model=0;
    double kx=10;
    double ky=10;
    double kz=10;
    
    double Min_node_pair_length, Max_node_pair_length, Average_node_pair_length;
    
    /**Returns the last saved volume*/
    double get_volume(void){
        return volume;
    }
    /**Calculate the voronoi area of  each  node and return a list that indicated the voronoi area of each node.*/
    vector<double> get_voronoi_node_area(){
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
    /**return normal vector of each node. The normal vector of node N is the average of vectors normal to triangles that have node N as one of their vertices. It is garanteed that the vector will point away from the Membrane centre of mass.*/
    vector<double> get_node_normal_vec(int node){
        return node_voronoi_normal_vec[node];
    }
    
    /**return the voronoi area associated with each node.*/
    double get_node_voronoi_area(int node){
        return node_voronoi_area[node];
    }
    /**Return the sum of the voronoi area of each node.*/
    double get_surface_area_voronoi(){
        return surface_area_voronoi;
    }
    double rescale_factor=1;
    
    void rescale_membrane(double factor){
        for(auto &pos: Node_Position){
            for(auto &coord: pos){
                coord *=factor;
            }
        }
        update_average_Membrane_radius();
    }
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
    /**Return the distance between node pair index. This list is initiated at the begining of the simulation.*/
    double get_node_pair_distance_in_Nm(int bond_num){
        return Node_Bond_distances_in_Nm[bond_num];
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
    /**Return the spontaneous bending angle between triangle pairs in radians. */
    double get_spontaneous_angle_in_Rad(int nodePairIndex){
        return Triangle_pair_angles_in_radians[nodePairIndex];
    }
    /**Return the spontaneous bending angle between triangle pairs in radians. */
    
    
    
    /**Returns the calculated number of triangles in the imported mesh file.*/
    int get_num_of_triangle_pairs(){
        return int(Triangle_Pair_Nodes.size());
    }
    /**Return the node IDs of the dihedral angles.*/
    vector<int> get_traingle_pair_nodes_list(int triangle_pair){
        return Triangle_Pair_Nodes[triangle_pair];
    }
    
    /**Return the node IDs of the dihedral angles.*/
    vector<int> get_dihedral_atoms_list(int index){
        return dihedral_atoms[index];
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
    double get_node_radius(int index){
        return Node_radius[index];
    }
    
    double get_membrane_weighted_radius(void){
        
        vector< double> unit;
        unit.resize(Num_of_Nodes,1);
        vector<double> radii;
        radii.resize(Num_of_Nodes,0);
        for (int i=0; i<Num_of_Nodes; i++) {
            radii[i]= spherical_positions[i][0];
        }
        double weighted_radius_sum = calc_vectorlist_vectorlist_surface_integral(radii,unit);
        weighted_radius_sum/=4*M_PI;
        cout<<weighted_radius_sum<<endl;
//        exit(0);
        return 0;
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
    
    void update_average_Membrane_radius(void){
        Radius=0;
        update_COM_position();
        for (int i=0; i<Num_of_Nodes; i++) {
            Radius += sqrt( (Node_Position[i][0]-COM_position[0])*(Node_Position[i][0]-COM_position[0])
                           +(Node_Position[i][1]-COM_position[1])*(Node_Position[i][1]-COM_position[1])
                           +(Node_Position[i][2]-COM_position[2])*(Node_Position[i][2]-COM_position[2]) );
        }
        Radius/=Num_of_Nodes;
    }
    double get_average_Membrane_radius(void){
        return Radius;
    }
    double get_new_Membrane_radius(void){
        return New_Radius;
    }
    
    std::map<string, vector<string> > Params;
    vector<string> insertOrder;
    vector<string> values;
    Membrane(){
        values.resize(2);
        
        values[0] ="Path/to/my/meshfile.extension";
        values[1] ="#Path to the mesh file. Supported formats: Blender's ply and Gmsh 2. The Membrane class cannot be initilised without a mesh file.";
        Params["MeshFile"] = values;
        insertOrder.push_back("MeshFile");
        
        values[0] ="1";
        values[1] ="#Mass asssigned to each node. Default value 1.";
        Params["NodeMass"] = values;
        insertOrder.push_back("NodeMass");
        
        values[0] ="Av";
        values[1] ="#Radius assigned to each node. If Av, half of the average bond distance of nodes will be used. If Au, half of the average of bonds connected to each node in the initial mesh will be used. If a value is inputed, the value will be used. The node radius is used to calculate the cutt-off and minimum energy distance for the 'Excluded Volume' and the 'Lennard-Jones' potential.";
        Params["NodeRadius"] = values;
        insertOrder.push_back("NodeRadius");
        
        values[0] ="0 0 0";
        values[1] ="#X, Y, Z components of a vector used to translate all coordinates of the Mesh befor beginning the simluation.";
        Params["CoordinateTranslateVector"] = values;
        insertOrder.push_back("CoordinateTranslateVector");
        
        values[0] ="0 0 0";
        values[1] ="#Vx, Vy, Vz components of a vector used to add to all initial node velocities befor beginning the simluation.";
        Params["VelocityShiftVector"] = values;
        insertOrder.push_back("VelocityShiftVector");
        
        values[0] ="H";
        values[1] ="#Set the bond potential. 'H' for harmonic. Default H";
        Params["SpringModel"] = values;
        insertOrder.push_back("SpringModel");
        
        values[0] ="0";
        values[1] ="#Set the bond potential rigidity coefficient. Default value 0.";
        Params["SpringCoeff"] = values;
        insertOrder.push_back("SpringCoeff");
        
        values[0] ="0";
        values[1] ="#Set the damping coefficient for non harmonic potentials. Default value 0.";
        Params["DampingCoeff"] = values;
        insertOrder.push_back("DampingCoeff");
        
        values[0] ="0";
        values[1] ="#Set bending potential (harmonic dihedral) rigidity coefficient. Default 0";
        Params["BendingCoeff"] = values;
        insertOrder.push_back("BendingCoeff");
        
        values[0] ="Av";
        values[1] ="#Set the rest length of the mesh springs using: 1) Au: The initial bond lengths of the mesh; 2) Av: The average node pair lengths; 3) \"value\": Where you type a specific  value.";
        Params["NominalLengthInNm"] = values;
        insertOrder.push_back("NominalLengthInNm");
        
        values[0] ="180";
        values[1] ="#Set the nominal angle (in degrees) of the mesh triangle pair potential using: 1) Au: The initial angles  of the mesh triangles; 2) Av: The average triangle pair angles; 3) \"value\": Where you type a specific  value.";
        Params["SpontaneousTriangleBendingAngleInDegrees"] = values;
        insertOrder.push_back("SpontaneousTriangleBendingAngleInDegrees");
        
        values[0] ="1";
        values[1] ="#Used to scale the Membrane coordinates. Default 1";
        Params["Scale"] = values;
        insertOrder.push_back("Scale");
        
        values[0] ="1 1 1";
        values[1] ="#Used to scale the Mesh in the X, Y, znd Z direction.";
        Params["XYZscale"] = values;
        insertOrder.push_back("XYZscale");
        
        values[0] ="false";
        values[1] ="#begin the simulation with a random orientation. Default value false";
        Params["InitRandomRotation"] = values;
        insertOrder.push_back("InitRandomRotation");
        
        values[0] ="false";
        values[1] ="#Write geometrical characteristics of the Membrane. It includes The Membrane volume, surface area, voronoi area of each node, normal vector of each node at the savinig time step. Default value false";
        Params["WantGeometricProps"] = values;
        insertOrder.push_back("WantGeometricProps");
        
        values[0] ="0 0 0";
        values[1] ="#X, Y, Z coordinates of a point inside the Membrane (usually the geometrical centre) used to build normal vectors on the trianular mesh that point to the outside of the Membrane. This proccess is only performed once when the Membrane is imported from the mesh file.";
        Params["XYZinMembrane"] = values;
        insertOrder.push_back("XYZinMembrane");
        
        values[0] ="0";
        values[1] ="#Set the Lennard Jones 12-6 sigma. Default value 0. If the Memebrane is interacting with another class, the Sigma between them will be calculated as the average of the class's sigmas: Sigma {Membrane & A} = 0.5(sigma{Membrane}+Sigma{A})";
        Params["LJsigma"] = values;
        insertOrder.push_back("LJsigma");
        
        values[0] ="0";
        values[1] ="#Set the Lennard Jones 12-6 epsillon. Default value 0. If the Memebrane is interacting with another class, the Epsillon between them will be calculated as the geometrical average of the class's epsilons: Epsillon {Membrane & A} = sqrt(epsillon{Membrane}*+epsillon{A}).";
        Params["LJepsilon"] = values;
        insertOrder.push_back("LJepsilon");
        
        values[0] ="0";
        values[1] ="#Under development. Do not use this flag.";
        Params["ExtForceModel"] = values;
        insertOrder.push_back("ExtForceModel");
        
        values[0] ="0 0 0";
        values[1] ="#Under development. Do not use this flag.";
        Params["ExtForceRigidity"] = values;
        insertOrder.push_back("ExtForceRigidity");
        
        values[0] ="-1";
        values[1] ="#Set if you want to change (linear) the Radius of the Membrane to change to a new value during the simulation. Default value is negative, indicating no change.";
        Params["UpdateRadius"] = values;
        insertOrder.push_back("UpdateRadius");
        
        values[0] ="0";
        values[1] ="#The Membrane radius will begin changing to the user provided value at this time point during the simulation (measured in pico seconds).";
        Params["UpdateBeginTimeInPs"] = values;
        insertOrder.push_back("UpdateBeginTimeInPs");
        
        values[0] ="0";
        values[1] ="#The Membrane radius will finish updating at this time point during the simulation (measured in pico seconds).";
        Params["UpdateEndTimeInPs"] = values;
        insertOrder.push_back("UpdateEndTimeInPs");
        
        
        
        
    }
    
    std::map<string, vector<string> > get_map(){
        return Params;
    }
    vector<string > get_insertOrder(){
        return insertOrder;
    }
    void assign_parameters(void);
    void consistancy_check(void);
};

#endif // MEMBRANE_H
