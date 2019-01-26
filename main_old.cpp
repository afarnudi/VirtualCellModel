//-----------------------HEADERS
#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <limits>
#include <cstdlib>
#include <random>
#include <string>
#include <math.h>
#include "General_constants.h"
#include "General_functions.hpp"
#include "General_Membrane.h"
#include "Membrane_functions.hpp"
#include "General_Actin.h"
#include "Actin_functions.hpp"
#include "Actin_membrane_shared_functions.hpp"
#include "General_Actin_Membrane_shared.h"
#include "General_Chromatin.h"
#include "Chromatin_functions.hpp"
#include "General_ECM.h"
#include "ECM_functions.hpp"

using namespace std;
//-----------------------HEADERS
void Membrane_Force_Calculator (double Membrane_Node_Position[][3],double Membrane_Node_Velocity[][3],double Membrane_Node_Force [][3],int Membrane_Node_Pair_list[][2],int Membrane_Triangle_Pair_Nodes[][4],double &Total_Potential_Energy, int Membrane_num_of_Triangle_Pairs, int Membrane_num_of_Node_Pairs);// updates forces + relavant potential energy

int parallelORantiparallel( double xpos[3][5] );// + for parallel  -for anti parallel (used in force cacculator of membrane)

double potentialenergy(double Membrane_Node_Position[][3],int Membrane_Node_Pair_list[][2],int Membrane_Triangle_Pair_Nodes[][4] ,int Membrane_triangle_list[Membrane_num_of_Triangles][3], int Outer_Membrane_num_of_triangles, int Nucleus_Membrane_num_of_triangles, int Membrane_num_of_Triangle_Pairs, int Membrane_num_of_Node_Pairs);//  potential energy

void MonteCarlo(  int Membrane_triangle_list[Membrane_num_of_Triangles][3],int Membrane_Triangle_Pair_Nodes[][4],int Membrane_Node_Pair_list[][2],double &Total_Potential_Energy,double Membrane_Node_Position[][3],int Membrane_Normal_direction[Membrane_num_of_Triangles][2], int Outer_Membrane_num_of_triangles, int Nucleus_Membrane_num_of_triangles, int Membrane_num_of_Triangle_Pairs, int Membrane_num_of_Node_Pairs);
void trapeziumUpdateMonteCarlo(int Membrane_Triangle_Pair_Nodes[][4],int p1, int p2,int p3,int p4, int &linepos,int &lineA,int &lineB,int &lineC,int &lineD,int &Aa,int &Bb,int &Cc,int &Dd, int Membrane_num_of_Triangle_Pairs);
void updatetriangle(int p1, int p2,int p3,int p4 ,int Membrane_triangle_list[Membrane_num_of_Triangles][3],int &l1,int &l2, double  Membrane_Node_Position [][3], int Membrane_Normal_direction[Membrane_num_of_Triangles][2] );
double periodiccondition(double dx );
double Membrane_surface_area_calculator(double  Membrane_Node_Position [][3],  int Membrane_triangle_list[Membrane_num_of_Triangles][3], int Outer_Membrane_num_of_triangles);

double surfaceareaNucleus(double  Membrane_Node_Position [][3],  int Membrane_triangle_list[Membrane_num_of_Triangles][3], int Outer_Membrane_num_of_triangles);
void ConstantSurfaceForceLocalTriangles(double Membrane_Node_Position[][3],double Membrane_Node_Force[][3],int Membrane_triangle_list[Membrane_num_of_Triangles][3], int Outer_Membrane_num_of_triangles);

double VolumeNucleus(double Membrane_Node_Position[][3],int Membrane_triangle_list[Membrane_num_of_Triangles][3],int Membrane_Normal_direction[Membrane_num_of_Triangles][2], int Outer_Membrane_num_of_triangles);

double VolumeMemerane(double Membrane_Node_Position[][3],int tri[3][Membrane_num_of_Triangles],int Membrane_Normal_direction[Membrane_num_of_Triangles][2], int Outer_Membrane_num_of_triangles);
// ------------------------------membrane



//---------------------------------actin
void  Actin_Force_calculator( double  Actin_Node_Position [][3],double  Actin_Node_VelocityRungKuta [][3],double  Actin_Node_Force [][3],double Actin_Node_Pair_List[][3], int Actin_num_of_Bonds, double &Total_Potential_Energy);
//double  energyActin( double  Actin_Node_Position [][3],double  Actin_Node_Velocity [][3],double  Actin_Node_Force [][3],double Actin_Node_Pair_List[][3]   );
void Actin_Membrane_Barrier(double Actin_Node_Position[][3], double Actin_Node_Velocity[][3], double Membrane_Node_Position[][3], double Membrane_Node_Velocity[][3],  int Membrane_triangle_list[Membrane_num_of_Triangles][3], int Membrane_Normal_direction[Membrane_num_of_Triangles][2]);
void Actin_Membrane_Barrier_2(double Actin_Node_Position[][3], double Actin_Node_Velocity[][3], double Membrane_Node_Position[][3], double Membrane_Node_Velocity[][3],  int Membrane_triangle_list[Membrane_num_of_Triangles][3], int Membrane_Normal_direction[Membrane_num_of_Triangles][2]);

void Nucleus_Membrane_Barrier(int Nucleus_Membrane_list_of_Nodes[], double Membrane_Node_Position[][3], double Membrane_Node_Velocity[][3], int Membrane_triangle_list[Membrane_num_of_Triangles][3], int Membrane_Normal_direction[Membrane_num_of_Triangles][2], int Outer_Membrane_num_of_triangles, int  Outer_Membrane_num_of_Nodes, int Membrane_num_of_Nodes);
void Nucleus_Membrane_Barrier_2(int Nucleus_Membrane_list_of_Nodes[], double Membrane_Node_Position[][3], double Membrane_Node_Velocity[][3], int Membrane_triangle_list[Membrane_num_of_Triangles][3], int Membrane_Normal_direction[Membrane_num_of_Triangles][2], int Outer_Membrane_num_of_triangles, int  Outer_Membrane_num_of_Nodes, int Membrane_num_of_Nodes);

//-------------------------------chromaton--------------------------

void Chromatin_Force_calculator(double (&Chromatin_Bead_Position)[Chromatin_num_of_Beads][3],double (&Chromatin_Bead_Velocity)[Chromatin_num_of_Beads][3],double (&Chromatin_Bead_Force)[Chromatin_num_of_Beads][3], double &Total_Potential_Energy);  // calculate force of inside the chromatin beads on each other
void hardsphereforcechromatin(double Chromatin_Bead_Position[Chromatin_num_of_Beads][3], double Chromatin_Bead_Velocity[Chromatin_num_of_Beads][3], double Chromatin_Bead_Force[Chromatin_num_of_Beads][3], bool chromatin_contant_matrix_calculation_flag, float CmLaststeps[Chromatin_num_of_Beads][Chromatin_num_of_Beads], float Cm[Chromatin_num_of_Beads][Chromatin_num_of_Beads], float Cm1[Chromatin_num_of_Beads][Chromatin_num_of_Beads], float Cm3[Chromatin_num_of_Beads][Chromatin_num_of_Beads], double   Dissimilaritycut);  // calculate force of inside the chromatin beads on each other

void Chromatin_membrane_Barrier(double Chromatin_Bead_Position[Chromatin_num_of_Beads][3], double Chromatin_Bead_Velocity[Chromatin_num_of_Beads][3], double Chromatin_Bead_Force[Chromatin_num_of_Beads][3], double Membrane_Node_Position[][3], double Membrane_Node_Velocity[][3], int Membrane_triangle_list[Membrane_num_of_Triangles][3], int Membrane_Normal_direction[Membrane_num_of_Triangles][2], int Outer_Membrane_num_of_triangles);

void Chromatin_membrane_Barrier_2(double Chromatin_Bead_Position[Chromatin_num_of_Beads][3], double Chromatin_Bead_Velocity[Chromatin_num_of_Beads][3], double Chromatin_Bead_Force[Chromatin_num_of_Beads][3], double Membrane_Node_Position[][3], double Membrane_Node_Velocity[][3], int Membrane_triangle_list[Membrane_num_of_Triangles][3], int Membrane_Normal_direction[Membrane_num_of_Triangles][2], int Outer_Membrane_num_of_triangles);



//______loop
#define loopchain_on   1  //1=yes 0=no
#define number_of_loops_in_each_chain  5
#define force_loop_speing_treashould  30000
void loop_force_chromatin(double Chromatin_Bead_Position[Chromatin_num_of_Beads][3],double Chromatin_Bead_Force[Chromatin_num_of_Beads][3], double &Total_Potential_Energy);
//______loop



//contact matrix:

int  StepsCmLaststeps= 30000; //(step-30000) ; // calculate  contact matrix for For the last these Steps]

//dissimilarity distance
#define StepCm0  10000  // first steps to calculate avraged contact matrix for first StepCm0 steps
#define StepCm1  10000  // calculate contact every StepCm1 steps
#define ReadCm0rfomFile  0 // 1=yes 0=no
#define StepCm2  10000  // avaging steps for window data collecting

#define Calculate_spatialcontactmatrix  1 //0=no 1=yes
#define size_spatialcontactmatrix  100  //size(resuloution of matrix)
#define effective_thickness_of_spatial_dis  1.0  // calculate an efective disk aornd the com on nuclus
int spatial_contact_matrix_x [ size_spatialcontactmatrix] [size_spatialcontactmatrix]={0};
int spatial_contact_matrix_y [ size_spatialcontactmatrix] [size_spatialcontactmatrix]={0};
int spatial_contact_matrix_z [ size_spatialcontactmatrix] [size_spatialcontactmatrix]={0};
void update_spatial_contact_matrix(double contact_coordinates[3]);
double COM_of_nucleus_membrane[3];
void COM_of_nucleus_membrane_function(double Membrane_Node_Position[][3],int Nucleus_Membrane_list_of_Nodes[], int  Outer_Membrane_num_of_Nodes, int Membrane_num_of_Nodes);
//-------------------------------chromaton-------------------------


//--------------------------------ECM--------------------------------
#define AdhesiveIslandOffOrOn   0.0  // 1.0=on  0.0=off
#define NumberOfAdhesiveIslandNodes 4
#define ECM_Flexibility 1.0 // 1.0=flexible 0.0=rigid

#define sigmaECM 20.0// force cont of ECM-membrane  interactio
#define ECM_LJ_just_Repultion  0 // 0=no 1=yes

#define epsilonECM  1.0//0.600 // force range of ECM-membrane interaction
#define ECM_kelvin_damping_coefficient   1000.0
#define miuDampECM 4
#define fECMcuttoff 500
#define ECM_Min_Gradient_Coefficient 5.0
#define ECM_Max_Gradient_Coefficient 10*ECM_Min_Gradient_Coefficient
#define ECM_Gradient_length 70.0
#define ECM_Node_Mass 5.0//10.0//5.0



void  ECM_rigidity_constructor(double  ECM_Node_Position [][3],double ECM_Node_Pair_List[][3], double ECM_varying_stiffness_coefficient[], int ECM_num_of_Bonds);

void  ECM_Force( double  ECM_Node_Position [][3],double  ECM_Node_Velocity [][3],double  ECM_Node_Force [][3],double ECM_Node_Pair_List[][3]  ,double ECM_upper_surface_Node_Pairs[] ,double ECM_varying_stiffness_coefficient[], int ECM_num_of_Bonds, double &Total_Potential_Energy);
int interactingECM,interactingMemrane;

void cellshift( double  Membrane_Node_Position [][3],double Actin_Node_Position[Actin_num_of_Nodes][3],double Chromatin_Bead_Position[Chromatin_num_of_Beads][3],double  ECM_Node_Position [][3] ,double  Membrane_Node_Velocity [][3],double Actin_Node_Velocity[Actin_num_of_Nodes][3],double Chromatin_Bead_Velocity[Chromatin_num_of_Beads][3], int Membrane_num_of_Nodes);

#define membraneshiftinXdirection 0
#define membraneshiftinZdirection 0
#define cell_downward_speed 0.0
#define extengingtriangleECM   0.0

///__new interaction strategy :
#define strechingInteractionFrom  1.55//* epsilonECM
#define strechingInteractionTill   3.30   //  cutoff distance between element and ECM to actiave force * epsilonECM
#define Membrane_ECM_interaction_update_step   50 // should be 50 in contactescase but here is suspended cell

#define considerTimesHowFarElemets  7 // it should be in this scope (1,5)
#define maximumHowmanyECMCouldInteractWithmembrane  200//70 //
#define ECMsurfaceDeformationConstraint    0  // 0=no constraint  1=constraint
#define ECMstrechingforce   80.0
#define ECMrepultionforceToTrianglesOrnode  0 // 0=to triangles  1=to nodes

double asymmetryConstInZ=1.0;


void Membrane_ECM_interaction(int istep,double Membrane_Node_Position [][3],double Membrane_Node_Velocity [][3],double Membrane_Node_Force [][3], double  ECM_Node_Position [][3],double  ECM_Node_Velocity [][3],double  ECM_Node_Force [][3],int ECM_surface_triangle_list[ECM_Surface_num_of_Triangles][3],int Membrane_triangle_list[Membrane_num_of_Triangles][3],int Membrane_Normal_direction[Membrane_num_of_Triangles][2],  int ECMandmembrane[][ 1+maximumHowmanyECMCouldInteractWithmembrane] ,double Cellcom[3], int Outer_Membrane_num_of_triangles, int Membrane_num_of_Nodes);


void Membrane_ECM_interaction_4(int MD_Step, double Membrane_Node_Position[][3], double Membrane_Node_Velocity[][3], double Membrane_Node_Force[][3], double ECM_Node_Position[][3], double ECM_Node_Velocity[][3], double ECM_Node_Force[][3], int ECM_surface_triangle_list[ECM_Surface_num_of_Triangles][3],int Membrane_triangle_list[Membrane_num_of_Triangles][3], int Membrane_Normal_direction[Membrane_num_of_Triangles][2], vector<vector<int> > &ECM_Membrane_Trinagle_neighbours_2, double Cellcom[3], int Outer_Membrane_num_of_triangles, int  Outer_Membrane_num_of_Nodes, int Outer_Membrane_list_of_Nodes[], vector<int>  &Membrane_Edge_triangle_and_ECM_neighbours, double &Total_Potential_Energy);

void Membrane_ECM_interaction_5(int MD_Step, double Membrane_Node_Position[][3], double Membrane_Node_Velocity[][3], double Membrane_Node_Force[][3], double ECM_Node_Position[][3], double ECM_Node_Velocity[][3], double ECM_Node_Force[][3], int ECM_surface_triangle_list[ECM_Surface_num_of_Triangles][3],int Membrane_triangle_list[Membrane_num_of_Triangles][3], int Membrane_Normal_direction[Membrane_num_of_Triangles][2], vector<int> &ECM_Membrane_Trinagle_neighbours_3, double Cellcom[3], int Outer_Membrane_num_of_triangles, int  Outer_Membrane_num_of_Nodes, int Outer_Membrane_list_of_Nodes[], vector<int>  &Membrane_Edge_triangle_and_ECM_neighbours, double &Total_Potential_Energy);

void Membrane_ECM_interactionHELPER(double  ECM_Node_Position [][3],int indexA,int indexB,int indexC,double n[3],double D[3],double (&a1),double (&a2),double (&a3));
void Membrane_ECM_interactionHELPER_newstrategy(double ECM_Node_Position [][3], int temp_ECM_node_A, int temp_ECM_node_B, int temp_ECM_node_C, double ECM_triangle_ABxAC_unit_vector[3], double D[3], double (&a1), double (&a2), double (&a3));


//_________Migration

#define migratingcell  0  //0=no  1=yes
double cell_polarizatioz_direction_migration[3]={0.0,0.0,1.0};  // for more info go to the paper: Dynamics of Cell Ensembles on Adhesive Micropatterns:   Bridging the Gap between Single Cell Spreading and Collective Cell Migration
double etha_migration= 1.0/8.0;
double miu_migration=  5.0 ;
double Back_Contraction_Coefficient = 0.40;
#define  controll_migration_with_total_trakhing_force  1  //0=no   1=yes
#define total_tracking_foece    1000.0


//_______________Integrins: attachment cytoskeleton and ECM:

#define calculateAngleFromCOM  1// 0no 1yes
#define ECM_membrane_attachment  0  // 0n 1yes
#define kactintoECM       10.0
#define UpdateAttachement  100
#define Actin_Node_ECM_detachment_penalty     5000 // The time penalty Actin Nodes on the Membrane (not Nucleus) have to wait until they can reattach (if the conditions are right) to the ECM.
#define ActinECMattachrange  2.0  //distance of attaching
#define rangeactinECM    1.60 //1.6  // distance of deataching
//ACTIVE force disrtibution  >> to make  Fcom=0 , distribute force to what parts?   toButtomNodes+toOuuterNodes+nucleus=1
#define toButtomNodes    0.0
#define toOuuterNodes    0.50
#define toNucleusnodes   0.50
double total_tracking_foece_now;



void CellCOM( double com[3],double  Membrane_Node_Position [][3],double  Actin_Node_Position [][3] ,double Chromatin_Bead_Position[Chromatin_num_of_Beads][3], int Membrane_num_of_Nodes);

void Vector_transformation (double MV[3],double  M[3][3] ,double V[3]);



//_________MIgration

//--------------------------------EMC--------------------------------






//-------------------------------Sovlent-----------------------------
#define solventon  0 //0=off  1=on
#define showsolvent   1 //  1=show   0=dont show they are too many!
#define BoundaryType  6 //(1=slip and no force) (2=no slip and no force) (3=slip  and no force  and initial flow) (4=noslip and no force and  initial velocity)  (5=noslip and constant force) (6=shear flow)
#define Fsolvent   0.01 // if BoundaryType==5
#define Vinitialflow -2.0 // initial flow if  BoundaryType==3,4
#define Vshear   30.0 // if BoundaryType==6
#define collisiontimestep  30  /// collision/inreaction == must be integer
#define interactionstep    10  /// collision/inreaction == must be integer
#define  msolvent 0.5  //mass of solvent particle
#define  lbs  0.4 // thikness of membrane film in solvebt-triangle interaction
#define wallthickness 0.1 // thickness of wall in slip or no-slip boundary condition
int const Lx = 2.0;/// must be even number!!!         size of simulation box
int const Ly = 2.0;/// must be even number!!!         size of simulation box
int const Lz = 2.0;/// must be even number!!!         size of simulation box
#define VmaxSolvent 1 // range of randomly distributed velocities in initialization
int const solventdensityinxox=1;// number of solvent particles in eachb box in initialization
const int nb = Lx*Ly*Lz ; //number of boxes
const int nsolvent = Lx*Ly*Lz*solventdensityinxox; // number of solvent particles
void initialzesolvent();
void solventBoundaries();
void collision();
void updateHowmanySolventsineachBox();
void findcom(double vcom[3][Membrane_num_of_Triangles],double pcom[3][Membrane_num_of_Triangles],int Membrane_triangle_list[Membrane_num_of_Triangles][3], double Membrane_Node_Position[][3], double (&Membrane_Node_Velocity)[][3]);
void interaction(double Membrane_Node_Position[][3], double Membrane_Node_Velocity[][3], int trianglememebrane[3][Membrane_num_of_Triangles],int Membrane_Normal_direction[Membrane_num_of_Triangles][2]);

double xsolvent[3][nsolvent], vsolvent[3][nsolvent];
int eachSolventinwichBox [nb][2*solventdensityinxox]; //  wich particle is in each box!
int howmanySolventsineachBox[nb];                       // how many particle in each box
#define spliteachboxfordiagram   1
//-------------------------------Sovlent-----------------------------


//-----------------------------thermostat------------------
double kineticenergymembrane (double  Membrane_Node_Velocity [][3], int Membrane_num_of_Nodes);
double kineticenergyactin (double  Actin_Node_Velocity [][3] );
double kineticenergysolvent ();
double kineticenergychromatin (double Chromatin_Bead_Velocity[Chromatin_num_of_Beads][3]);
double kineticenergyECM (double  ECM_Node_Velocity [][3]);
void Thermostat(int istep, double Membrane_Node_Velocity[][3], double Actin_Node_Velocity[][3], double Chromatin_Bead_Velocity[Chromatin_num_of_Beads][3], double ECM_Node_Velocity[][3], int Membrane_num_of_Nodes);
//-----------------------------thermostat------------------


//---------------------------restart-----------------------

void restartsave(double Membrane_Node_Position[][3], double Membrane_Node_Velocity [][3], int Membrane_triangle_list[Membrane_num_of_Triangles][3], int Membrane_Normal_direction[Membrane_num_of_Triangles][2], int Membrane_Triangle_Pair_Nodes[][4], int bondslist[][2], double Actin_Node_Position[Actin_num_of_Nodes][3], double Actin_Node_Velocity[Actin_num_of_Nodes][3], double Actin_Node_Pair_List[][3], double Chromatin_Bead_Position[Chromatin_num_of_Beads][3], double Chromatin_Bead_Velocity[Chromatin_num_of_Beads][3], double ECM_Node_Position[][3], double ECM_Node_Velocity [][3], int Membrane_num_of_Triangle_Pairs, int Actin_num_of_Bonds, int Membrane_num_of_Nodes);


void restartread(double Membrane_Node_Position [][3], double Membrane_Node_Velocity [][3], int Membrane_triangle_list[Membrane_num_of_Triangles][3], int Membrane_Normal_direction[Membrane_num_of_Triangles][2], int Membrane_Triangle_Pair_Nodes[][4], int bondslist[][2], double Actin_Node_Position[Actin_num_of_Nodes][3], double Actin_Node_Velocity[Actin_num_of_Nodes][3],double Actin_Node_Pair_List[][3], double Chromatin_Bead_Position[Chromatin_num_of_Beads][3], double Chromatin_Bead_Velocity[Chromatin_num_of_Beads][3], double ECM_Node_Position [][3], double ECM_Node_Velocity [][3], int Membrane_num_of_Triangle_Pairs, int Actin_num_of_Bonds, int Membrane_num_of_Nodes);
//---------------------------restart-----------------------


//Povray------
# define export_povray_step_distance 10000  // clear!
void povray_output_creator(int currentStep, double Membrane_Node_Position[][3], int  Membrane_triangle_list[Membrane_num_of_Triangles][3], int Membrane_Normal_direction[Membrane_num_of_Triangles][2], int Membrane_Node_Pair_list[][2], double Actin_Node_Position[Actin_num_of_Nodes][3], double Actin_Node_Pair_List[][3], double Chromatin_Bead_Position[Chromatin_num_of_Beads][3], double  ECM_Node_Position[][3], double ECM_Node_Pair_List[][3], int ECM_surface_triangle_list[ECM_Surface_num_of_Triangles][3], int Outer_Membrane_num_of_triangles, int Membrane_num_of_Node_Pairs, int Outer_Membrane_num_of_Node_Pairs, int Actin_num_of_Bonds, int ECM_num_of_Bonds, int Membrane_num_of_Nodes);

//Ali
void Membrane_Actin_shared_Node_Force_calculator (double Membrane_Node_Position[][3],double  Actin_Node_Position [Actin_num_of_Nodes][3], double Membrane_Node_Force [][3], double  Actin_Node_Force [Actin_num_of_Nodes][3], int Membrane_Actin_shared_Node_list[Actin_Membrane_shared_num_of_Nodes][2], double Membrane_Node_velocity[][3], double Actin_Node_velocity[Actin_num_of_Nodes][3]);// updates forces + relavant potential energy
void generate_initial_condition_report (string initial_condition_file_name, int Membrane_num_of_Nodes);

#define Actin_membrane_stiff_spring_coefficient 400
#define Actin_membrane_damping_coefficient 0.0
#define Membrane_ECM_Max_cos_triangle_interaction_angle 0.2
#define ECM_Membrane_Radius_of_Hard_Sphere_Interaction 1.0//0.50
#define Membrane_ECM_triangle_interaction_update_step   50 // should be 50 in contactescase but here is suspended cell
#define Membrane_ECM_triangle_unbinding_update_step   50 // should be 50 in contactescase but here is suspended cell
#define Membrane_ECM_triangle_interaction_rate 10
#define Membrane_ECM_interaction_strength 20.0
#define Membrane_ECM_Binding_prob 0.5
#define Membrane_ECM_UnBinding_prob 0.2
#define spreading_flag 0.0
#define spreading_force_magnitude -100.0
#define spreading_force_min_range 1.0
#define spreading_force_max_range 4.0
#define spreading_force_cos_triangle_interaction_angle 0.4
#define energy_calculation_flag 0.0

int Membrane_Num_of_Nodes_reader (string membrane_mesh_file_name);
int Membrane_Num_of_Nodes_reader (string membrane_mesh_file_name){
    ifstream read; //This is the main ifstream that will read the Gmesh-Membrane generated file
    read.open(membrane_mesh_file_name.c_str() ); //It should be noted that the name of the file should not contain '-'. I don't know why but the memory managnet of the arrays (at the very least) in the programme will collapse when we use '-' in the file name.
    int temp_int; // This is just a temp intiger charachter that we use to read unnecessary Gmesh generated intigers. We never use these intigers in the actual programme.
    string temp_string;
    for (int i=0; i<6; i++) {
        read>> temp_string;
    }
    read>> temp_int;
    cout<<"Membrane number of Nodes is=\t"<<temp_int<<endl;
    return temp_int;
}



int main() //main**
{
    clock_t tStart = clock();//Time the programme
    time_t t = time(0);   // get time now
    struct tm * now = localtime( & t );
    int num_of_test_particles=13;
    int test_particles[num_of_test_particles][3];
    //direction identifier
    for (int i=0; i<num_of_test_particles; i++) {
        test_particles[i][0]=12+i;
        test_particles[i][1]=-11.5+i*1;
        test_particles[i][2]=0;
        
        //        test_particles[i][0]=12-i;
        //        test_particles[i][1]=-11.5+i*1;
        //        test_particles[i][2]=i*i;
        
        //        if (i<(num_of_test_particles/2)) {
        //            test_particles[i][0]=12;
        //        }   else {
        //            test_particles[i][0]=15;
        //        }
        //        test_particles[i][1]=-11.5+i*1;
        //        test_particles[i][2]=0;
    }
    //Edges
    test_particles[num_of_test_particles-3][0]=35;
    test_particles[num_of_test_particles-3][1]=10;
    test_particles[num_of_test_particles-2][0]=-35;
    test_particles[num_of_test_particles-2][1]=10;
    //Centre
    test_particles[num_of_test_particles-1][0]=0;
    test_particles[num_of_test_particles-1][1]=30;
    test_particles[num_of_test_particles-1][2]=0;
    
//    double test_cout=123456.78910111213;
    
    char buffer [80];
    strftime (buffer,80,"%Y-%m-%d-%H:%M",now);
    //outputfiles:
    string traj_file_name, initial_condition_file_name, energy_file_name;
    
    traj_file_name="results/RBC_";
    traj_file_name +=buffer;
    initial_condition_file_name = traj_file_name;
    energy_file_name=initial_condition_file_name;
    traj_file_name +=".xyz";
    initial_condition_file_name+=".txt";
    
    ofstream trajectory;
    trajectory.open(traj_file_name.c_str() );
    trajectory << std:: fixed;
    
    
    //Energy Calculation variables-Begin
    double Membrane_total_potential_Energy=0.0, Actin_total_potential_Energy=0.0, Chromatin_total_potential_Energy=0.0, ECM_total_potential_Energy=0.0, ECM_membrane_total_potential_Energy=0.0;
    energy_file_name+="energy.txt";
    ofstream write_energy;
    write_energy.open(energy_file_name.c_str());
    write_energy<<"Time\tTotal Energy\tMembrane\tActin\tChromatin\tECM\tInteraction_4"<<endl;
    //Energy Calculation variables-End
    
    
    bool resume=false;//If true the programme will call upon the 'restartread' function and read the positions and velocities of the nodes from an external file that is determined inside the function.
    bool chromatin_contant_matrix_calculation_flag=false;
    //-----------------------------------------------membrane:
    int counter = 0;//This counter counts the number of MD steps and wil be used to determin on which steps the preogramme is saved ('restartsave').
    string membrane_mesh_file_name="membrane";
    //*******************************************************************************************************
    /*BUG
     |---\   |    |  /---\
     |    |  |    |  |
     |---<   |    |  |  -\
     |    |  |    |  |   |
     |---/   \----/  \---/
     */
    //*******************************************************************************************************
    //***************** the membrane file is modefied *****************************
    //*******************************************************************************************************
    const int Membrane_num_of_Nodes=Membrane_Num_of_Nodes_reader(membrane_mesh_file_name);
    
    double Membrane_Node_Position[Membrane_num_of_Nodes][3], Membrane_Node_Velocity[Membrane_num_of_Nodes][3],Membrane_Node_Force[Membrane_num_of_Nodes][3];
    int Membrane_triangle_list[Membrane_num_of_Triangles][3];
    double Total_Kinetic_Energy,Total_Potential_Energy; // total kineti , potential ,and mechanical energy
    
    generate_initial_condition_report(initial_condition_file_name, Membrane_num_of_Nodes);
    
    cout << "Please double check the following parameters\nMembrane Radius= "<< Membrane_Radius <<"\tand Nucleus_Membrane_radius= " <<Nucleus_Membrane_radius<<endl;
    cout <<"Initialising the programme ..."<<endl;
    
    Membrane_constructor(Membrane_Node_Position, Membrane_Node_Velocity, Membrane_Node_Force, Membrane_triangle_list, membrane_mesh_file_name);
    
    int Membrane_Normal_direction[Membrane_num_of_Triangles][2];// n[i][0]= +1(-1) if the membrane triangle belongs to the outer membrane(Nucleus); And  n[i][1]= +1(-1) if the normal product of the triangle nodes (with the default order) is positive (negative).
    
    int  Outer_Membrane_num_of_triangles;       // will be updated in Membrane_Normal_direction_Identifier
    int  Nucleus_Membrane_num_of_triangles;     // will be updated in normalsorter
    Membrane_Normal_direction_Identifier(Membrane_Node_Position, Membrane_triangle_list, Membrane_Normal_direction, Outer_Membrane_num_of_triangles, Nucleus_Membrane_num_of_triangles);
    
    int  Outer_Membrane_num_of_Nodes;  // Is calculated in 'Membrane_Normal_direction_Identifier'
    Outer_Membrane_Identifier(Membrane_Normal_direction, Membrane_triangle_list, Outer_Membrane_num_of_triangles, Outer_Membrane_num_of_Nodes);
    
    int Outer_Membrane_list_of_Nodes[Outer_Membrane_num_of_Nodes];
    int Nucleus_Membrane_list_of_Nodes[Membrane_num_of_Nodes-Outer_Membrane_num_of_Nodes];
    
    Membrane_and_Nucleus_Node_list_builder( Membrane_Node_Position,Nucleus_Membrane_list_of_Nodes,Outer_Membrane_list_of_Nodes,Membrane_triangle_list, Outer_Membrane_num_of_triangles);
    
    int  Membrane_num_of_Triangle_Pairs;
    Membrane_num_of_Triangle_Pairs=Membrane_triangle_pair_counter( Membrane_triangle_list);
    int Membrane_Triangle_Pair_Nodes[Membrane_num_of_Triangle_Pairs][4]; // pos1 pos2 pos3 and po4 of all interactions are stored here
    Membrane_Triangle_Pair_Identifier(Membrane_triangle_list, Membrane_Triangle_Pair_Nodes, Membrane_num_of_Triangle_Pairs);
    
    
    int Outer_Membrane_num_of_Node_Pairs ; // usefull in POV ray- modulate in void sortingbonds(int bondslist[][2],int tri[3][Membrane_num_of_Triangles])
    int Membrane_num_of_Node_Pairs;
    Membrane_num_of_Node_Pairs=Membrane_num_of_Node_Pair_Counter(Membrane_triangle_list, Outer_Membrane_num_of_triangles, Outer_Membrane_num_of_Node_Pairs, Membrane_num_of_Nodes);
    
    int Membrane_Node_Pair_list[Membrane_num_of_Node_Pairs][2];
    //*******************************************************************************************************
    /*BUG
     |---\   |    |  /---\
     |    |  |    |  |
     |---<   |    |  |  -\
     |    |  |    |  |   |
     |---/   \----/  \---/
     */
    //*******************************************************************************************************
    //***************** Since we use the 'Membrane_num_of_Node_Pair_Counter' just to count the number of Node pairs, we have to write a similar function to build the actual 'Membrane_Node_Pair_list', which I have assigned 'Membrane_num_of_Node_Pair_Counter' to do so. Let t be noted that not calling this function will not cause any problems if the 'resume' is on!. *****************************
    //*******************************************************************************************************
    if (resume==false) {
        Membrane_num_of_Node_Pair_Counter_2(Membrane_Node_Pair_list, Membrane_triangle_list, Outer_Membrane_num_of_triangles, Membrane_num_of_Node_Pairs);
    }
    
    //-----------------------------------------------membrane:
    
    //*******************************************************************************************************
    /*BUG
     |---\   |    |  /---\
     |    |  |    |  |
     |---<   |    |  |  -\
     |    |  |    |  |   |
     |---/   \----/  \---/
     */
    //*******************************************************************************************************
    //***************** We spend so much energy and time to seperate the Outer  *****************************
    //***************** Membrane from the Nucleus. Why not just build 2? ************************************
    //*******************************************************************************************************
    
    //--------------------------------Migration---------------------------
    double centerOfmassoftheCell[3];
    
    ofstream movement_of_cell_COM;
    movement_of_cell_COM.open("results/movement_of_cell_COM.txt");
    movement_of_cell_COM << std:: fixed;
    
    
    ofstream total_tracking_foece_ofstream;
    total_tracking_foece_ofstream.open("results/total_tracking_foece.txt");
    total_tracking_foece_ofstream << std:: fixed;
    //--------------------------------Migration---------------------------
    
    //---------------------------------------------------actin:
    int Actin_num_of_Bonds=Actin_Node_Pair_Identifier(); // the value will be set automatically in actin parameter function
    
    double Actin_Node_Pair_List[Actin_num_of_Bonds][3]; //In this array 3 numbers are stored for each Actin Node pairs. The first two are of course the label (number) of the nodes that are a pair. The third element is used to store the distance between these pairs. The distance is used in the force calculations for the Maxwell spring initial length. The initial length is updated during each step so we have a diferent initial length, hence the Maxwell spring.
    
    double Actin_Node_Position[Actin_num_of_Nodes][3];
    double Actin_Node_Velocity[Actin_num_of_Nodes][3];
    double Actin_Node_Force[Actin_num_of_Nodes][3];
    double Actin_Node_VelocityRungKuta[Actin_num_of_Nodes][3];  //only used in  Rungekuta step
    Actin_constructor(Actin_Node_Position, Actin_Node_Velocity, Actin_Node_Force, Actin_Node_Pair_List, Actin_num_of_Bonds);
    
    //---------------------------------------------------actin
    
    
    //-----------------------------actin-membrane part: shared beads(suppose to have masses like mass of membrane)
    int Membrane_Actin_shared_Node_list[Actin_Membrane_shared_num_of_Nodes][2];// IMPORTANT: this list contains the indices of shared beads on the membrane and actin. [0] for membrane and [1] for actin.
    Membrane_Actin_shared_Node_Identifier(Membrane_Actin_shared_Node_list, Membrane_Node_Position, Actin_Node_Position, Membrane_num_of_Nodes);
    
    //-------------------------------chromaton--------------------------
    double Chromatin_Bead_Position[Chromatin_num_of_Beads][3];
    double Chromatin_Bead_Velocity[Chromatin_num_of_Beads][3];
    double Chromatin_Bead_Force[Chromatin_num_of_Beads][3];
    float CmLaststeps[Chromatin_num_of_Beads][Chromatin_num_of_Beads];
    double   Dissimilaritycut = sigmachromatin*Chromatin_Scaling_Factor*1.5;  //  within this distance chromatins are consedered in contact
    Chromatin_constructor(Chromatin_Bead_Position, Chromatin_Bead_Velocity, Chromatin_Bead_Force, chromatin_contant_matrix_calculation_flag, CmLaststeps);
    
    ofstream d2_10ofstream;
    d2_10ofstream.open("results/d2_10.txt"); // refrecne at the beginning
    d2_10ofstream << std:: fixed;
    
    ofstream d2_32ofstream;
    d2_32ofstream.open("results/d2_32.txt"); // refrecne at the beginning
    d2_32ofstream << std:: fixed;
    //--d2_32_REDUCED.txt
    float ContactProbablilityNumber; //stores the total prob number of contacts
    ofstream REDUCEDd2_32ofstream;
    REDUCEDd2_32ofstream.open("results/d2_32_REDUCED.txt");
    REDUCEDd2_32ofstream << std:: fixed;
    
    ofstream cm3ofstream;
    cm3ofstream.open("results/Cm3.txt");
    cm3ofstream << std:: fixed;
    float Cm[Chromatin_num_of_Beads][Chromatin_num_of_Beads];  // contact matrix at now
    float Cm0[Chromatin_num_of_Beads][Chromatin_num_of_Beads];  // avraged contact matrix for first StepCm0 steps
    float Cm1[Chromatin_num_of_Beads][Chromatin_num_of_Beads];  // avraged contact matrix for last StepCm1 steps
    float Cm2[Chromatin_num_of_Beads][Chromatin_num_of_Beads];  // avraged contact matrix for last 2*StepCm2 steps
    float Cm3[Chromatin_num_of_Beads][Chromatin_num_of_Beads];  // avraged contact matrix for last StepCm2 steps
    float d2_10;  // d^2 between cm0 and cm1
    float d2_32;  // d^2 between cm2 and cm3
    float ContactNumber; //stores the number of contacts
    // spatialcontact matrix
    
    int Flag=0; // useful in contact prob matrix
    
    if (chromatin_contant_matrix_calculation_flag==true) {
        
        
        for(int i1=0;i1<Chromatin_num_of_Beads;i1++)
        {
            for(int j1=0;j1<Chromatin_num_of_Beads;j1++)
            {
                Cm0[i1][j1]=0;
                Cm1[i1][j1]=0;
                Cm2[i1][j1]=0;
                Cm3[i1][j1]=0;
            }
        }
    }
    //-------------------------------chromaton--------------------------
    
    
    //--------------------------------ECM--------------------------------
    double ECM_Node_Position[ECM_num_of_Nodes][3];
    double ECM_Node_Velocity[ECM_num_of_Nodes][3];
    double ECM_Node_Force[ECM_num_of_Nodes][3];
    int ECM_surface_triangle_list[ECM_Surface_num_of_Triangles][3]; // first Membrane_num_of_triangles elements are for outer membrane and after that nucleus elementsk
    int ECM_num_of_Bonds;
    ECM_num_of_Bonds=ECM_Node_Pair_Identifier();
    double ECM_upper_surface_Node_Pairs[ECM_num_of_Bonds]; // 0.0 = Node pair is located on the upper surface of th ECM. 1= Node pair is not on the upper surface, hence in the volume or the edge or buttom surface.
    double ECM_varying_stiffness_coefficient[ECM_num_of_Bonds];
    double ECM_Node_Pair_List[ECM_num_of_Bonds][3]; // In the first two elements, the label (intigers) of the node pairs are stored, and in the last element their distance. The distance is used in the force calculations for the Maxwell spring initial length. The initial length is updated during each step so we have a diferent initial length, hence the Maxwell spring.
    
    ECM_constructor(ECM_Node_Position, ECM_Node_Velocity, ECM_Node_Force, ECM_surface_triangle_list, ECM_Node_Pair_List, ECM_upper_surface_Node_Pairs, ECM_num_of_Bonds);
    
    // resume**
    ECM_rigidity_constructor(ECM_Node_Position, ECM_Node_Pair_List, ECM_varying_stiffness_coefficient, ECM_num_of_Bonds);
    vector<vector<int> > ECM_Membrane_Trinagle_neighbours_2;//A 2D list of the trinagles (first index) and their ECM neighbours (seconde index)
    ECM_Membrane_Trinagle_neighbours_2.resize(Outer_Membrane_num_of_triangles);
    
    vector<int>  ECM_Membrane_Trinagle_neighbours_3;//A 1D list of the trinagles (first index) and its ECM neighbour
    ECM_Membrane_Trinagle_neighbours_3.resize(Outer_Membrane_num_of_triangles);
    
    vector<int>  Membrane_Edge_triangle_and_ECM_neighbours;// A 2d list of Membrane triangles that meke up the contact edge of the membrane (first element) and their corresponding ECM neighbours (second element). The membrane triangles are located right above the ECM triangles at a specific distance and at an angle (perpindicular) to the ECMs.
    Membrane_Edge_triangle_and_ECM_neighbours.resize(Outer_Membrane_num_of_triangles);
    
    //--------------------------------ECM--------------------------------
    
    
    
    
    
    if(solventon==1)
    {
        //-------------------------Solvent
        initialzesolvent();
        //-------------------------Solvent
    }
    
    
    //*******************************************************************************************************
    /*BUG
     |---\   |    |  /---\
     |    |  |    |  |
     |---<   |    |  |  -\
     |    |  |    |  |   |
     |---/   \----/  \---/
     */
    //*******************************************************************************************************
    //***************** Improvement: I should chek to see if we can avoid some of the 'constructors' If we are using the 'restartread' function.
    //*******************************************************************************************************
    
    
    
    
    cellshift(Membrane_Node_Position, Actin_Node_Position, Chromatin_Bead_Position, ECM_Node_Position, Membrane_Node_Velocity, Actin_Node_Velocity, Chromatin_Bead_Velocity, Membrane_num_of_Nodes);
    if( resume==true )
    {
        restartread(Membrane_Node_Position,Membrane_Node_Velocity, Membrane_triangle_list, Membrane_Normal_direction, Membrane_Triangle_Pair_Nodes, Membrane_Node_Pair_list, Actin_Node_Position, Actin_Node_Velocity, Actin_Node_Pair_List, Chromatin_Bead_Position, Chromatin_Bead_Velocity, ECM_Node_Position, ECM_Node_Velocity, Membrane_num_of_Triangle_Pairs, Actin_num_of_Bonds, Membrane_num_of_Nodes);
        cellshift(Membrane_Node_Position, Actin_Node_Position, Chromatin_Bead_Position, ECM_Node_Position, Membrane_Node_Velocity, Actin_Node_Velocity, Chromatin_Bead_Velocity, Membrane_num_of_Nodes);
    }
    //*******************************************************************************************************
    /*BUG
     |---\   |    |  /---\
     |    |  |    |  |
     |---<   |    |  |  -\
     |    |  |    |  |   |
     |---/   \----/  \---/
     */
    //*******************************************************************************************************
    //***************** It would seem that the ECM has been constructed twice, once here and once the ECM section.
    // We should probably find the reason why or get rid of one.
    //*******************************************************************************************************
    
    ECM_constructor(ECM_Node_Position, ECM_Node_Velocity, ECM_Node_Force, ECM_surface_triangle_list, ECM_Node_Pair_List, ECM_upper_surface_Node_Pairs, ECM_num_of_Bonds);
    
    //==========================================================================================================================
    //==========================================================================================================================
    //============================================== Beginning of the main MD loop ==============================================
    //==========================================================================================================================
    //==========================================================================================================================
    
    //============================================== adhesive bead ==============================================
    double adhesive_bead_position[3];
    int membrane_triangle_adhesive_bead_adhesion_point;
    vector < int > membrane_triangle_adhesive_bead_neighbour_list;
    double adhesion_distance_threshold;
    bool membrane_bead_adhesion=false;
    int bead_membrane_neighbour_update_step;
    
    //Check membrane for contact point
    if (membrane_bead_adhesion==false) {
        for (int membrane_tri_index=0; membrane_tri_index<Membrane_num_of_Triangles; membrane_tri_index++) {
            double membrane_triangle_com[3];
            double temp_membrane_bead_distance;
            membrane_triangle_com[0]=(Membrane_Node_Position[Membrane_triangle_list[membrane_tri_index][0]][0]+Membrane_Node_Position[Membrane_triangle_list[membrane_tri_index][1]][0]+Membrane_Node_Position[Membrane_triangle_list[membrane_tri_index][2]][0])/3.0;
            membrane_triangle_com[1]=(Membrane_Node_Position[Membrane_triangle_list[membrane_tri_index][0]][1]+Membrane_Node_Position[Membrane_triangle_list[membrane_tri_index][1]][1]+Membrane_Node_Position[Membrane_triangle_list[membrane_tri_index][2]][1])/3.0;
            membrane_triangle_com[2]=(Membrane_Node_Position[Membrane_triangle_list[membrane_tri_index][0]][2]+Membrane_Node_Position[Membrane_triangle_list[membrane_tri_index][1]][2]+Membrane_Node_Position[Membrane_triangle_list[membrane_tri_index][2]][2])/3.0;
            
            temp_membrane_bead_distance=(membrane_triangle_com[0]-adhesive_bead_position[0])*(membrane_triangle_com[0]-adhesive_bead_position[0])+(membrane_triangle_com[1]-adhesive_bead_position[1])*(membrane_triangle_com[1]-adhesive_bead_position[1])+ (membrane_triangle_com[2]-adhesive_bead_position[2])*(membrane_triangle_com[2]-adhesive_bead_position[2]);
            
            if (temp_membrane_bead_distance<adhesion_distance_threshold*adhesion_distance_threshold) {
                membrane_triangle_adhesive_bead_adhesion_point=membrane_tri_index;
                membrane_bead_adhesion=true;
                break;
            }
        }
    }
    
    if (membrane_bead_adhesion==true) {
        if (MD_step % bead_membrane_neighbour_update_step==0) {
            membrane_triangle_adhesive_bead_neighbour_list.clear();
            for (int membrane_node_index=0; membrane_node_index<Membrane_num_of_Nodes; membrane_node_index++) {
                double temp_bead_membrane_node_distance;
                temp_bead_membrane_node_distance=(Membrane_Node_Position[membrane_node_index][0]-adhesive_bead_position[0])*(Membrane_Node_Position[membrane_node_index][0]-adhesive_bead_position[0])+(Membrane_Node_Position[membrane_node_index][1]-adhesive_bead_position[1])*(Membrane_Node_Position[membrane_node_index][1]-adhesive_bead_position[1])+(Membrane_Node_Position[membrane_node_index][2]-adhesive_bead_position[2])*(Membrane_Node_Position[membrane_node_index][2]-adhesive_bead_position[2]);
                
                if (temp_bead_membrane_node_distance<adhesion_distance_threshold*adhesion_distance_threshold) {
                    membrane_triangle_adhesive_bead_neighbour_list.push_back(membrane_node_index);
                }
                
            }
        }
    }
    
    //bead membrane Force calculations:
    for (int neighbour_index=0; neighbour_index<membrane_triangle_adhesive_bead_neighbour_list.size(); neighbour_index++) {
        //apply force between nodes.
    }
    
    
    int counter2=0;
    cout<<"Beginning the MD loop\n";
    for(int MD_Step=0 ;MD_Step<=MD_num_of_steps ; MD_Step++)
    {
//        cout<<test_cout<<endl;
        //------------------------------------------------------------------------------------
        //------------------------------- Beginning of Membrane position update --------------
        //---------------------------------- Velocity Verlet ---------------------------------
        //------------------------------------------------------------------------------------
        
        for(int j=0 ; j<Membrane_num_of_Nodes ; j++)
        {
            Membrane_Node_Position[j][0] += Membrane_Node_Velocity[j][0]*MD_Time_Step - Membrane_Node_Force[j][0]*MD_Time_Step*MD_Time_Step/(Membrane_Node_Mass*2.0);
            Membrane_Node_Position[j][1] += Membrane_Node_Velocity[j][1]*MD_Time_Step - Membrane_Node_Force[j][1]*MD_Time_Step*MD_Time_Step/(Membrane_Node_Mass*2.0);
            Membrane_Node_Position[j][2] += Membrane_Node_Velocity[j][2]*MD_Time_Step - Membrane_Node_Force[j][2]*MD_Time_Step*MD_Time_Step/(Membrane_Node_Mass*2.0);
        }
        
        //------------------------------------------------------------------------------------
        //------------------------------- beginning of Actin position update -------------------------------
        //---------------------------------- Velocity Verlet ---------------------------------
        //------------------------------------------------------------------------------------
        
        for(int j=0 ; j<Actin_num_of_Nodes ; j++)  //actin
        {
            Actin_Node_Position[j][0] += Actin_Node_Velocity[j][0]*MD_Time_Step - Actin_Node_Force[j][0]*MD_Time_Step*MD_Time_Step/(Actin_Node_Mass*2.0);
            Actin_Node_Position[j][1] += Actin_Node_Velocity[j][1]*MD_Time_Step - Actin_Node_Force[j][1]*MD_Time_Step*MD_Time_Step/(Actin_Node_Mass*2.0);
            Actin_Node_Position[j][2] += Actin_Node_Velocity[j][2]*MD_Time_Step - Actin_Node_Force[j][2]*MD_Time_Step*MD_Time_Step/(Actin_Node_Mass*2.0);
        }
        //------------------------------------------------------------------------------------
        //------------------------------- End of Actin-Membrane Shared Nodes position update -----------------
        //------------------------------------------------------------------------------------
        
        //------------------------------------------------------------------------------------
        //------------------------------- End of Actin position update ---------------------------------------
        //------------------------------------------------------------------------------------
        
        //------------------------------------------------------------------------------------
        //------------------------------- End of Membrane position update ------------------------------------
        //------------------------------------------------------------------------------------
        
        //------------------------------------------------------------------------------------
        //------------------------------- Beginning of Chromatin position update -----------------------------
        //------------------------------------------------------------------------------------
        
        for(int j=0 ; j<Chromatin_num_of_Beads ; j++)
        {
            Chromatin_Bead_Position[j][0] += Chromatin_Bead_Velocity[j][0]*MD_Time_Step-Chromatin_Bead_Force[j][0]*MD_Time_Step*MD_Time_Step/(2.0*Chromatin_Bead_Mass);
            Chromatin_Bead_Position[j][1] += Chromatin_Bead_Velocity[j][1]*MD_Time_Step-Chromatin_Bead_Force[j][1]*MD_Time_Step*MD_Time_Step/(2.0*Chromatin_Bead_Mass);
            Chromatin_Bead_Position[j][2] += Chromatin_Bead_Velocity[j][2]*MD_Time_Step-Chromatin_Bead_Force[j][2]*MD_Time_Step*MD_Time_Step/(2.0*Chromatin_Bead_Mass);
        }
        
        //------------------------------------------------------------------------------------
        //------------------------------- End of Chromatin position update -----------------------------
        //------------------------------------------------------------------------------------
        
        //------------------------------------------------------------------------------------
        //------------------------------- Beginning of ECM position update -----------------------------
        //------------------------------------------------------------------------------------
        
        for(int j=0 ; j<ECM_num_of_Nodes ; j++)
        {
            // Ptential Bug. Because I have yet to understand the reason behind this condition.
            //*******************************************************************************************************
            //***************** From what I understand, this condition implies that the position of Nodes on the ECM
            // that pass a certain depth (In this code y determins depth/hight) will not be updated, hence, the Nodes
            // will freez after moving to a 'Membrane_Centre_distance_from_ECM + ECM_Thickness' depth.
            //*******************************************************************************************************
            
            if( ECM_Node_Position[j][1]!= -( Membrane_Centre_distance_from_ECM + ECM_Thickness) && ECM_Node_Position[j][1]>= -( Membrane_Centre_distance_from_ECM + 0.8*ECM_Thickness))
            {
                ECM_Node_Position[j][0] += ECM_Flexibility * (ECM_Node_Velocity[j][0]*MD_Time_Step - ECM_Node_Force[j][0]*MD_Time_Step*MD_Time_Step/(ECM_Node_Mass*2.0) );
                ECM_Node_Position[j][1] += ECM_Flexibility * (ECM_Node_Velocity[j][1]*MD_Time_Step -ECM_Node_Force[j][1]*MD_Time_Step*MD_Time_Step/(ECM_Node_Mass*2.0) );
                ECM_Node_Position[j][2] += ECM_Flexibility * (ECM_Node_Velocity[j][2]*MD_Time_Step -ECM_Node_Force[j][2]*MD_Time_Step*MD_Time_Step/(ECM_Node_Mass*2.0) );
            }
            
            
        }
        //------------------------------------------------------------------------------------
        //------------------------------- Beginning of the periodic condition -----------------------------
        //------------------------------------------------------------------------------------
        // Potential Bug.
        //*******************************************************************************************************
        //***************** From what I understand, this condition implies that the position of Nodes on the Membrane
        // and the Chromatin beads. I think that the Actin network hasto get a periodic condition as well.
        //*******************************************************************************************************
        
        if (Periodic_condtion_status == 1.0) {
            for(int ww=0 ; ww<2 ; ww++)
            {
                for(int j=0 ; j<Membrane_num_of_Nodes ; j++)  //
                {
                    if( Membrane_Node_Position[j][ww]> (Lbox+1.0)/2.0  )
                    {
                        
                        Membrane_Node_Position[j][ww]=Membrane_Node_Position[j][ww]-(Lbox+1.0 );
                    }
                    if( Membrane_Node_Position[j][ww]< (-Lbox-1.0)/2.0 )
                    {
                        
                        Membrane_Node_Position[j][ww]=Membrane_Node_Position[j][ww]+(Lbox+1.0 );
                    }
                }
            }
            
            for(int ww=0 ; ww<2 ; ww++)
            {
                for(int j=0 ; j<Chromatin_num_of_chains ; j++)  //
                {
                    if( Chromatin_Bead_Position[j][ww]> (Lbox+1.0)/2.0  )
                    {
                        
                        Chromatin_Bead_Position[j][ww] += -(Lbox+1.0 );
                    }
                    if( Chromatin_Bead_Position[j][ww]< (-Lbox-1.0)/2.0 )
                    {
                        Chromatin_Bead_Position[j][ww] += (Lbox+1.0 );
                    }
                }
            }
        }
        
        
        //------------------------------------------------------------------------------------
        //------------------------------- End of the periodic condition -----------------------------
        //------------------------------------------------------------------------------------
        
        //--------------------------------------------------------------------
        
        for(int j=0 ; j<Membrane_num_of_Nodes ; j++)  // loop to count every particle  and update its velocity
        {
            Membrane_Node_Velocity[j][0] += - Membrane_Node_Force[j][0]*MD_Time_Step/(Membrane_Node_Mass*2.0);
            Membrane_Node_Velocity[j][1] += - Membrane_Node_Force[j][1]*MD_Time_Step/(Membrane_Node_Mass*2.0);
            Membrane_Node_Velocity[j][2] += - Membrane_Node_Force[j][2]*MD_Time_Step/(Membrane_Node_Mass*2.0);
        }
        
        for(int j=0 ; j<Actin_num_of_Nodes ; j++)  //actin
        {
            
            Actin_Node_Velocity[j][0] += -Actin_Node_Force[j][0]*MD_Time_Step/(Actin_Node_Mass*2.0);
            Actin_Node_Velocity[j][1] += -Actin_Node_Force[j][1]*MD_Time_Step/(Actin_Node_Mass*2.0);
            Actin_Node_Velocity[j][2] += -Actin_Node_Force[j][2]*MD_Time_Step/(Actin_Node_Mass*2.0);
            
        }
        
        for(int j=0 ; j<Chromatin_num_of_Beads ; j++)  // loop to encount every particle  and update its velocity
        {
            Chromatin_Bead_Velocity[j][0] += -Chromatin_Bead_Force[j][0]*MD_Time_Step/(2.0*Chromatin_Bead_Mass);
            Chromatin_Bead_Velocity[j][1] += -Chromatin_Bead_Force[j][1]*MD_Time_Step/(2.0*Chromatin_Bead_Mass);
            Chromatin_Bead_Velocity[j][2] += -Chromatin_Bead_Force[j][2]*MD_Time_Step/(2.0*Chromatin_Bead_Mass);
        }
        //*******************************************************************************************************
        /*Improvment
         |---\   |    |  /---\
         |    |  |    |  |
         |---<   |    |  |  -\
         |    |  |    |  |   |
         |---/   \----/  \---/
         */
        //*******************************************************************************************************
        //***************** It is a good idea to set a condition for the 'Flexibility' of the *******************
        //***************** Membrane, so that the velocity and Force are not run altogether   *******************
        //*******************************************************************************************************
        for(int j=0 ; j<ECM_num_of_Nodes ; j++)  //ECM
        {
            ECM_Node_Velocity[j][0] += -ECM_Node_Force[j][0]*MD_Time_Step/(ECM_Node_Mass*2.0);
            ECM_Node_Velocity[j][1] += -ECM_Node_Force[j][1]*MD_Time_Step/(ECM_Node_Mass*2.0);
            ECM_Node_Velocity[j][2] += -ECM_Node_Force[j][2]*MD_Time_Step/(ECM_Node_Mass*2.0);
        }
        
        //------------------------------------------------------------------------------------
        //------------------------------- Beginning of the Force Calculations -----------------------------
        //------------------------------------------------------------------------------------
        
        Total_Potential_Energy=0.0;
        Total_Kinetic_Energy=0.0;
        //        cout<<"Total_Potential_Energy= "<<Total_Potential_Energy<<endl;
        // All Forces are set to zero. Throughout the programme, the relative force function will 'add' the required fore to the respective array elements.
        for(int j=0 ; j<Membrane_num_of_Nodes ; j++)
        {
            Membrane_Node_Force[j][0]=0.0;
            Membrane_Node_Force[j][1]=0.0;
            Membrane_Node_Force[j][2]=0.0;
        }
        for(int j=0 ; j<Actin_num_of_Nodes ; j++)
        {
            Actin_Node_Force[j][0]=0.0;
            Actin_Node_Force[j][1]=0.0;
            Actin_Node_Force[j][2]=0.0;
        }
        for(int j=0 ; j<Chromatin_num_of_Beads ; j++)
        {
            Chromatin_Bead_Force[j][0]=0.0;
            Chromatin_Bead_Force[j][1]=0.0;
            Chromatin_Bead_Force[j][2]=0.0;
        }
        for(int j=0 ; j<ECM_num_of_Nodes ; j++)
        {
            ECM_Node_Force[j][0]=0.0;
            ECM_Node_Force[j][1]=0.0;
            ECM_Node_Force[j][2]=0.0;
        }
        
        
        
        Membrane_Force_Calculator(Membrane_Node_Position, Membrane_Node_Velocity, Membrane_Node_Force, Membrane_Node_Pair_list, Membrane_Triangle_Pair_Nodes, Total_Potential_Energy,  Membrane_num_of_Triangle_Pairs, Membrane_num_of_Node_Pairs); // **********calling aceel update
        
        
        ConstantSurfaceForceLocalTriangles( Membrane_Node_Position, Membrane_Node_Force, Membrane_triangle_list, Outer_Membrane_num_of_triangles);
        if (MD_Step % mcstep == 0)// collsi
        {
            for(int s=0; s< fluidity*Membrane_num_of_Nodes ;s++)
            {
                MonteCarlo(Membrane_triangle_list, Membrane_Triangle_Pair_Nodes, Membrane_Node_Pair_list, Total_Potential_Energy, Membrane_Node_Position, Membrane_Normal_direction, Outer_Membrane_num_of_triangles, Nucleus_Membrane_num_of_triangles, Membrane_num_of_Triangle_Pairs, Membrane_num_of_Node_Pairs);
                counter2=0;
            }
            counter2++;
        }
        //*******************************************************************************************************
        /*Potetial bug
         |---\   |    |  /---\
         |    |  |    |  |
         |---<   |    |  |  -\
         |    |  |    |  |   |
         |---/   \----/  \---/
         */
        //*******************************************************************************************************
        //***************** Potetial bug: The Rung-Kuta method is used to further the velocity  *****************
        //***************** of the Actin nodes to be used in the ACtin_Force_calculator in the  *****************
        //***************** next step. First of all why do we need to used the Rung-Kuta for this ***************
        //***************** particular step (if it is good we should use it everywhere). And second, ************
        //***************** I should double check to see if we need the previous, hence no need for *************
        //*****************  the Rung-Kuta step, or the current step for the disipating force  ******************
        //*****************  calculations in the Actin_Force_calculator function.  ******************************
        //*******************************************************************************************************
        if (energy_calculation_flag==1.0) {
            Membrane_total_potential_Energy=Total_Potential_Energy;
            //        cout<<"Total_Potential_Energy= "<<Total_Potential_Energy<<endl;
        }
        
        
        for(int j=0 ; j<Actin_num_of_Nodes ; j++)  //Rung kuta update
        {
            Actin_Node_VelocityRungKuta[j][0] += -Actin_Node_Force[j][0]*MD_Time_Step/(Actin_Node_Mass);
            Actin_Node_VelocityRungKuta[j][1] += -Actin_Node_Force[j][1]*MD_Time_Step/(Actin_Node_Mass);
            Actin_Node_VelocityRungKuta[j][2] += -Actin_Node_Force[j][2]*MD_Time_Step/(Actin_Node_Mass);
        }
        
        Actin_Force_calculator(Actin_Node_Position, Actin_Node_VelocityRungKuta, Actin_Node_Force, Actin_Node_Pair_List, Actin_num_of_Bonds, Total_Potential_Energy); // updates with runge kuta
        if (energy_calculation_flag==1.0) {
            Actin_total_potential_Energy=Total_Potential_Energy-Membrane_total_potential_Energy;
            //        cout<<"Total_Potential_Energy= "<<Total_Potential_Energy<<endl;
        }
        
        
        Membrane_Actin_shared_Node_Force_calculator(Membrane_Node_Position, Actin_Node_Position, Membrane_Node_Force,Actin_Node_Force, Membrane_Actin_shared_Node_list, Membrane_Node_Velocity, Actin_Node_Velocity);
        
        
        Chromatin_Force_calculator(Chromatin_Bead_Position, Chromatin_Bead_Velocity, Chromatin_Bead_Force, Total_Potential_Energy);
        
        loop_force_chromatin(Chromatin_Bead_Position, Chromatin_Bead_Force, Total_Potential_Energy);
        
        if (energy_calculation_flag==1.0) {
            Chromatin_total_potential_Energy=Total_Potential_Energy-Membrane_total_potential_Energy-Actin_total_potential_Energy;
            //        cout<<"Total_Potential_Energy= "<<Total_Potential_Energy<<endl;
        }
        
        
        COM_of_nucleus_membrane_function(Membrane_Node_Position, Nucleus_Membrane_list_of_Nodes,  Outer_Membrane_num_of_Nodes, Membrane_num_of_Nodes);
        hardsphereforcechromatin(Chromatin_Bead_Position, Chromatin_Bead_Velocity, Chromatin_Bead_Force, chromatin_contant_matrix_calculation_flag, CmLaststeps, Cm, Cm1, Cm3, Dissimilaritycut);
        
        //        Chromatin_membrane_Barrier(Chromatin_Bead_Position, Chromatin_Bead_Velocity, Chromatin_Bead_Force, Membrane_Node_Position, Membrane_Node_Velocity, Membrane_triangle_list, Membrane_Normal_direction, Outer_Membrane_num_of_triangles);
        Chromatin_membrane_Barrier_2(Chromatin_Bead_Position, Chromatin_Bead_Velocity, Chromatin_Bead_Force, Membrane_Node_Position, Membrane_Node_Velocity, Membrane_triangle_list, Membrane_Normal_direction, Outer_Membrane_num_of_triangles);
        
        if(MD_Step%Membrane_barrier_calculation_rate==0)
        {
            //            Actin_Membrane_Barrier(Actin_Node_Position, Actin_Node_Velocity, Membrane_Node_Position, Membrane_Node_Velocity, Membrane_triangle_list, Membrane_Normal_direction);
            Actin_Membrane_Barrier_2(Actin_Node_Position, Actin_Node_Velocity, Membrane_Node_Position, Membrane_Node_Velocity, Membrane_triangle_list, Membrane_Normal_direction);
            //            Nucleus_Membrane_Barrier(Nucleus_Membrane_list_of_Nodes, Membrane_Node_Position, Membrane_Node_Velocity, Membrane_triangle_list, Membrane_Normal_direction, Outer_Membrane_num_of_triangles, Outer_Membrane_num_of_Nodes);
            Nucleus_Membrane_Barrier_2(Nucleus_Membrane_list_of_Nodes, Membrane_Node_Position, Membrane_Node_Velocity, Membrane_triangle_list, Membrane_Normal_direction, Outer_Membrane_num_of_triangles, Outer_Membrane_num_of_Nodes, Membrane_num_of_Nodes);
            
        }
        
        ECM_Force(ECM_Node_Position, ECM_Node_Velocity, ECM_Node_Force, ECM_Node_Pair_List, ECM_upper_surface_Node_Pairs, ECM_varying_stiffness_coefficient, ECM_num_of_Bonds, Total_Potential_Energy);
        
        if (energy_calculation_flag==1.0) {
            ECM_total_potential_Energy=Total_Potential_Energy-Membrane_total_potential_Energy-Actin_total_potential_Energy-Chromatin_total_potential_Energy;
            //        cout<<"Total_Potential_Energy= "<<Total_Potential_Energy<<endl;
        }
        
        
        if(calculateAngleFromCOM ==1)
        {
            CellCOM( centerOfmassoftheCell,Membrane_Node_Position,Actin_Node_Position ,Chromatin_Bead_Position, Membrane_num_of_Nodes);
        }
        
        
        //        Membrane_ECM_interaction_4(MD_Step, Membrane_Node_Position,Membrane_Node_Velocity, Membrane_Node_Force,ECM_Node_Position, ECM_Node_Velocity, ECM_Node_Force, ECM_surface_triangle_list,Membrane_triangle_list, Membrane_Normal_direction, ECM_Membrane_Trinagle_neighbours_2,centerOfmassoftheCell, Outer_Membrane_num_of_triangles, Outer_Membrane_num_of_Nodes, Outer_Membrane_list_of_Nodes, Membrane_Edge_triangle_and_ECM_neighbours, Total_Potential_Energy);
        Membrane_ECM_interaction_5(MD_Step, Membrane_Node_Position,Membrane_Node_Velocity, Membrane_Node_Force,ECM_Node_Position, ECM_Node_Velocity, ECM_Node_Force, ECM_surface_triangle_list,Membrane_triangle_list, Membrane_Normal_direction, ECM_Membrane_Trinagle_neighbours_3,centerOfmassoftheCell, Outer_Membrane_num_of_triangles, Outer_Membrane_num_of_Nodes, Outer_Membrane_list_of_Nodes, Membrane_Edge_triangle_and_ECM_neighbours, Total_Potential_Energy);
        
        if (energy_calculation_flag==1.0) {
            ECM_membrane_total_potential_Energy=Total_Potential_Energy-Membrane_total_potential_Energy-Actin_total_potential_Energy-Chromatin_total_potential_Energy-ECM_total_potential_Energy;
            //        cout<<"Total_Potential_Energy= "<<Total_Potential_Energy<<endl;
        }
        
        
        /*
         if (ECM_membrane_attachment ==1)
         {
         
         ActinECMattachment(MD_Step,Actin_Node_Position,Actin_Node_Velocity,Actin_Node_Force,ECM_Node_Position,ECM_Node_Velocity,ECM_Node_Force,ECM_surface_triangle_list,Num_of_Actin_Nodes_on_Membrane_list,Actin_Membrane_shared_Node_list,Actin_Node_ECM_attachment_info );
         
         update_cartesian_from_triangular_coordiante(  ECM_Node_Position , ECM_surface_triangle_list, Actin_Node_ECM_attachment_info,Num_of_Actin_Nodes_on_Membrane_list, Actin_Membrane_shared_Node_list);
         
         
         ActinECMattachmentForce( MD_Step,  Actin_Node_Position ,  Actin_Node_Velocity ,  Actin_Node_Force,  ECM_Node_Position ,ECM_Node_Velocity ,ECM_Node_Force , ECM_surface_triangle_list,Num_of_Actin_Nodes_on_Membrane_list,Actin_Membrane_shared_Node_list, Actin_Node_ECM_attachment_info);
         }
         */
        
        //------------------------------------------------------------------------------------
        //------------------------------- End of the Force Calculations -----------------------------
        //------------------------------------------------------------------------------------
        
        //------------------------------------------------------------------------------------
        //------------------------------- Beginning of the Velocity Calculations -----------------------------
        //------------------------------------------------------------------------------------
        
        for(int j=0 ; j<Membrane_num_of_Nodes ; j++)  // loop to count every particle  and update its velocity
        {
            Membrane_Node_Velocity[j][0] += - Membrane_Node_Force[j][0]*MD_Time_Step/(Membrane_Node_Mass*2.0);
            Membrane_Node_Velocity[j][1] += - Membrane_Node_Force[j][1]*MD_Time_Step/(Membrane_Node_Mass*2.0);
            Membrane_Node_Velocity[j][2] += - Membrane_Node_Force[j][2]*MD_Time_Step/(Membrane_Node_Mass*2.0);
        }
        
        for(int j=0 ; j<Actin_num_of_Nodes ; j++)  //actin
        {
            Actin_Node_Velocity[j][0] += -Actin_Node_Force[j][0]*MD_Time_Step/(Actin_Node_Mass*2.0);
            Actin_Node_Velocity[j][1] += -Actin_Node_Force[j][1]*MD_Time_Step/(Actin_Node_Mass*2.0);
            Actin_Node_Velocity[j][2] += -Actin_Node_Force[j][2]*MD_Time_Step/(Actin_Node_Mass*2.0);
        }
        
        
        for(int j=0 ; j<Chromatin_num_of_Beads ; j++)  // loop to encount every particle  and update its velocity
        {
            Chromatin_Bead_Velocity[j][0] += -Chromatin_Bead_Force[j][0]*MD_Time_Step/(2.0*Chromatin_Bead_Mass);
            Chromatin_Bead_Velocity[j][1] += -Chromatin_Bead_Force[j][1]*MD_Time_Step/(2.0*Chromatin_Bead_Mass);
            Chromatin_Bead_Velocity[j][2] += -Chromatin_Bead_Force[j][2]*MD_Time_Step/(2.0*Chromatin_Bead_Mass);
            //            if (isnan(Chromatin_Bead_Force[j][0])==true || isnan(Chromatin_Bead_Force[j][1])==true || isnan(Chromatin_Bead_Force[j][2])==true) {
            //                cout<<"here\n";
            //            }
        }
        
        for(int j=0 ; j<ECM_num_of_Nodes ; j++)  //ECM
        {
            
            ECM_Node_Velocity[j][0]= ECM_Node_Velocity[j][0]-ECM_Node_Force[j][0]*MD_Time_Step/(ECM_Node_Mass*2.0);
            ECM_Node_Velocity[j][1]= ECM_Node_Velocity[j][1]-ECM_Node_Force[j][1]*MD_Time_Step/(ECM_Node_Mass*2.0);
            ECM_Node_Velocity[j][2]= ECM_Node_Velocity[j][2]-ECM_Node_Force[j][2]*MD_Time_Step/(ECM_Node_Mass*2.0);
        }
        
        if(solventon==1)
        {
            //--------------------------------------------------------------------Solvent verlet
            
            if (MD_Step % interactionstep == 0)
            {
                // update position:
                if(BoundaryType==1 || BoundaryType==2 ||BoundaryType==3 ||BoundaryType==4 || BoundaryType==6 ) // no force
                {
                    for (int j = 0; j < nsolvent; j++) // update place of solvent particles + Randomshift
                    {
                        xsolvent[0][j] = xsolvent[0][j] + MD_Time_Step*interactionstep * vsolvent[0][j]; // Streaming
                        xsolvent[1][j] = xsolvent[1][j] + MD_Time_Step*interactionstep * vsolvent[1][j];
                        xsolvent[2][j] = xsolvent[2][j] + MD_Time_Step*interactionstep * vsolvent[2][j];
                    }
                    solventBoundaries();
                }
                
                else if(BoundaryType==5  ) // with force(Velocity verlet)
                {
                    for (int j = 0; j < nsolvent; j++) // update place of solvent particles + Randomshift
                    {
                        xsolvent[0][j] = xsolvent[0][j] + MD_Time_Step*interactionstep * vsolvent[0][j]+(MD_Time_Step*interactionstep*MD_Time_Step*interactionstep )*Fsolvent/msolvent; // Streaming
                        xsolvent[1][j] = xsolvent[1][j] + MD_Time_Step*interactionstep * vsolvent[1][j];
                        xsolvent[2][j] = xsolvent[2][j] + MD_Time_Step*interactionstep * vsolvent[2][j];
                    }
                    
                    for (int j = 0; j < nsolvent; j++) // update place of solvent particles + Randomshift
                    {
                        vsolvent[0][j] = vsolvent[0][j] + MD_Time_Step*interactionstep*Fsolvent/msolvent; // Streaming
                    }
                    solventBoundaries();
                    
                }
                
            }
            
            
            
            
            if (MD_Step % collisiontimestep == 0)// collsion step
            {
                double randomshftvector[3]; // stores the random shift vector in order to inverse it after collision
                // choose random vector
                randomshftvector[0] = (double) (rand() % 1000) / 1000 - 0.5; // each compenent between -0.5 and 0.5
                randomshftvector[1] = (double) (rand() % 1000) / 1000 - 0.5;
                randomshftvector[2] = (double) (rand() % 1000) / 1000 - 0.5;
                
                
                for (int j = 0; j < nsolvent; j++) // update place of solvent particles + Randomshift
                {
                    
                    
                    //_______________ Randomshift BLOCK :
                    xsolvent[0][j]=xsolvent[0][j] +  randomshftvector[0] ;
                    xsolvent[1][j]=xsolvent[1][j] +  randomshftvector[1] ;
                    xsolvent[2][j]=xsolvent[2][j] +  randomshftvector[2] ;
                    //keep solent in BOX
                    if (xsolvent[0][j] > Lx/2.0)
                        xsolvent[0][j] = xsolvent[0][j] - Lx ;
                    if (xsolvent[0][j] < -Lx/2.0)
                        xsolvent[0][j] = xsolvent[0][j] + Lx ;
                    
                    if (xsolvent[1][j] > Ly/2.0)
                        xsolvent[1][j] = xsolvent[1][j] - Ly ;
                    if (xsolvent[1][j] < -Ly/2.0)
                        xsolvent[1][j] = xsolvent[1][j] + Ly ;
                    
                    if (xsolvent[2][j] > Lz/2.0)
                        xsolvent[2][j] = xsolvent[2][j] - Lz ;
                    if (xsolvent[2][j] < -Lz/2.0)
                        xsolvent[2][j] = xsolvent[2][j] + Lz ;
                    //_______________ Randomshift BLOCK :
                }
                
                collision(); // update velocity in solvent
                
                
                for (int j = 0; j < nsolvent; j++) // Randomshift - INVERSE
                {
                    
                    //_______________ Randomshift BLOCK  : inversing
                    xsolvent[0][j]=xsolvent[0][j] -  randomshftvector[0] ;
                    xsolvent[1][j]=xsolvent[1][j] -  randomshftvector[1] ;
                    xsolvent[2][j]=xsolvent[2][j] -  randomshftvector[2] ;
                    //keep solent in BOX
                    if (xsolvent[0][j] > Lx/2.0)
                        xsolvent[0][j] = xsolvent[0][j] - Lx ;
                    if (xsolvent[0][j] < -Lx/2.0)
                        xsolvent[0][j] = xsolvent[0][j] + Lx ;
                    
                    if (xsolvent[1][j] > Ly/2.0)
                        xsolvent[1][j] = xsolvent[1][j] - Ly ;
                    if (xsolvent[1][j] < -Ly/2.0)
                        xsolvent[1][j] = xsolvent[1][j] + Ly ;
                    
                    if (xsolvent[2][j] > Lz/2.0)
                        xsolvent[2][j] = xsolvent[2][j] - Lz ;
                    if (xsolvent[2][j] < -Lz/2.0)
                        xsolvent[2][j] = xsolvent[2][j] + Lz ;
                    //_______________ Randomshift BLOCK :
                }
                
                updateHowmanySolventsineachBox();
                
            }
            
            if (MD_Step % interactionstep == 0) // interaction step
            {
                
                //    interaction(Membrane_Node_Position,Membrane_Node_Velocity,Membrane_triangle_list,Membrane_Normal_direction);
                
            }
            
            
        }
        
        //Energy calculations
        //        if (MD_Step%10000==0 & MD_Step!=0  ){
        //            //Kinetic Energies
        //            //Potential Energies
        //
        //        }
        // ______________________________________________________________SAVING AND CALCULATION SECTION
        if (MD_Step%10000==0 & MD_Step!=0  )
        {
            
            cout<<"Saving temp..."<<endl;
            restartsave(Membrane_Node_Position, Membrane_Node_Velocity, Membrane_triangle_list, Membrane_Normal_direction, Membrane_Triangle_Pair_Nodes, Membrane_Node_Pair_list,Actin_Node_Position, Actin_Node_Velocity, Actin_Node_Pair_List, Chromatin_Bead_Position,Chromatin_Bead_Velocity, ECM_Node_Position,ECM_Node_Velocity, Membrane_num_of_Triangle_Pairs, Actin_num_of_Bonds, Membrane_num_of_Nodes);
        }
        
        //--------------------------------Thermostate
        //*******************************************************************************************************
        /*BUG
         |---\   |    |  /---\
         |    |  |    |  |
         |---<   |    |  |  -\
         |    |  |    |  |   |
         |---/   \----/  \---/
         */
        //*******************************************************************************************************
        //***************** Potential BUG: from what i gather, the Thermostat is only applied *******************
        //***************** to the actin network and not the membrane *******************************************
        //*******************************************************************************************************
        if(counter%RunThermostatePerstep==0)
        {
            Thermostat(MD_Step, Membrane_Node_Velocity, Actin_Node_Velocity, Chromatin_Bead_Velocity, ECM_Node_Velocity, Membrane_num_of_Nodes);
        }
        
        
        
        
        
        //--------------------------------------------------------------------saving and friends!
        if (counter%savingstep==0  )  // SAVING AND CALCULATION SECTION
        {
            if (energy_calculation_flag==1.0) {
                write_energy<<counter<<"\t"<<Total_Potential_Energy<<"\t"<<Membrane_total_potential_Energy<<"\t"<<Actin_total_potential_Energy<<"\t"<<Chromatin_total_potential_Energy<<"\t"<<ECM_total_potential_Energy<<"\t"<<ECM_membrane_total_potential_Energy<<"\n";
            }
            
            //            cout << Membrane_surface_area_calculator(Membrane_Node_Position, Membrane_triangle_list, Outer_Membrane_num_of_triangles) << "  "<< s0membrane<<endl;
            //            cout << surfaceareaNucleus(Membrane_Node_Position, Membrane_triangle_list, Outer_Membrane_num_of_triangles ) << "  "<< s0nucleus<<endl;
            total_tracking_foece_ofstream << "total_tracking_foece_now= "<<total_tracking_foece_now<<endl;
            //*******************************************************************************************************
            /*BUG
             |---\   |    |  /---\
             |    |  |    |  |
             |---<   |    |  |  -\
             |    |  |    |  |   |
             |---/   \----/  \---/
             */
            //*******************************************************************************************************
            //***************** Potential BUG: why is the membrane damping force terminated after this step? ********
            //*******************************************************************************************************
            //
            //            if(i>1500)
            //            {
            //                membrane_damping_coefficient=0.0;
            //            }
            
            cout <<MD_Step<<endl;
            
            
            
            ///__________________________________________________Traj__________________________
            int lambda;
            
            trajectory << Membrane_num_of_Nodes+Actin_num_of_Nodes+nsolvent*showsolvent+Chromatin_num_of_Beads+ECM_num_of_Nodes+num_of_test_particles<<endl;//+Num_of_Actin_Nodes_on_Membrane_list <<endl;    // saving trajectories
            
            trajectory << " nodes  "<<endl;
            
            for(int j=0; j<Outer_Membrane_num_of_Nodes;j++) // saving trajectory
            {   lambda=Outer_Membrane_list_of_Nodes[j];
                trajectory <<"membrane"  <<setprecision(5)<< setw(20)<<Membrane_Node_Position[lambda][0]<< setw(20)<<Membrane_Node_Position[lambda][1]<< setw(20)<<Membrane_Node_Position[lambda][2]<<endl;
            }
            
            for(int j=0; j<Membrane_num_of_Nodes- Outer_Membrane_num_of_Nodes;j++) // saving trajectory
            {   lambda=Nucleus_Membrane_list_of_Nodes[j];
                trajectory <<"nucleus"  <<setprecision(5)<< setw(20)<<Membrane_Node_Position[lambda][0]<< setw(20)<<Membrane_Node_Position[lambda][1]<< setw(20)<<Membrane_Node_Position[lambda][2]<<endl;
            }
            
            
            for (int nchain=0;nchain<Chromatin_num_of_chains;nchain++)
            {
                for(int j=nchain*(Chromatin_num_of_Beads/Chromatin_num_of_chains)  ;j< (nchain+1)*(Chromatin_num_of_Beads/Chromatin_num_of_chains) ; j++ )  // all beads interaction whit the next one
                {
                    // trajectory << "c"<<nchain<<  setprecision(5)<< setw(20)<<Chromatin_Bead_Position[0][j]<< setw(20)<<Chromatin_Bead_Position[1][j]<< setw(20)<<Chromatin_Bead_Position[2][j]<<endl;
                    trajectory << "chromatin_"<< nchain <<" "<<setprecision(5)<< setw(20)<<Chromatin_Bead_Position[j][0]<< setw(20)<<Chromatin_Bead_Position[j][1]<< setw(20)<<Chromatin_Bead_Position[j][2]<<endl;
                }
            }
            
            
            for(int j=0; j< ECM_num_of_Nodes ;j++) // saving trajectory
            {
                trajectory << "ecm" <<  setprecision(5)<< setw(20)<<ECM_Node_Position[j][0]<< setw(20)<<ECM_Node_Position[j][1]<< setw(20)<<ECM_Node_Position[j][2]<<endl;
            }
            
            
            
            for(int j=0; j< Actin_num_of_Nodes ;j++) // saving trajectory
            {
                trajectory << "actin" <<  setprecision(5)<< setw(20)<<Actin_Node_Position[j][0]<< setw(20)<<Actin_Node_Position[j][1]<< setw(20)<<Actin_Node_Position[j][2]<<endl;
                
            }
            
            
            
            
            
            if(showsolvent==1)
            {
                for(int j=0; j< nsolvent;j++) // saving trajectory
                {
                    trajectory << "solvent"<<  setprecision(5)<< setw(20)<<xsolvent[0][j]<< setw(20)<<xsolvent[1][j]<< setw(20)<<xsolvent[2][j]<<endl;
                }
                
            }
            
            for(int j=0; j<num_of_test_particles;j++) // saving trajectory
            {
                trajectory <<"test"  <<setprecision(5)<< setw(20)<<test_particles[j][0]<< setw(20)<<test_particles[j][1]<< setw(20)<<test_particles[j][2]<<endl;
            }
            /*
             for(int j=0; j< Num_of_Actin_Nodes_on_Membrane_list;j++) // saving trajectory
             {
             lambda=Actin_Membrane_shared_Node_list[j];
             if( Actin_Node_ECM_attachment_info[lambda][0]==-1.0 ) //NOT ATTACHED : just hide them in chromatins!!!
             {
             trajectory << "g"<<  setprecision(5)<< setw(20)<<Chromatin_Bead_Position[0][0]<< setw(20)<<Chromatin_Bead_Position[0][1]<< setw(20)<<Chromatin_Bead_Position[0][2]<<endl;
             }
             else if( Actin_Node_ECM_attachment_info[lambda][0]==1.0 )
             {
             trajectory << "g" <<  setprecision(5)<< setw(20)<<Actin_Node_ECM_attachment_info[lambda][1]<< setw(20)<<Actin_Node_ECM_attachment_info[lambda][2]<< setw(20)<<Actin_Node_ECM_attachment_info[lambda][3]<<endl;
             }
             
             }
             */
            
            
            ///__________________________________________________Traj__________________________
            
            CellCOM( centerOfmassoftheCell,Membrane_Node_Position,Actin_Node_Position ,Chromatin_Bead_Position, Membrane_num_of_Nodes);
            
            movement_of_cell_COM << MD_Step << " "<< centerOfmassoftheCell[0]<<" "<<centerOfmassoftheCell[1]<<" "<<centerOfmassoftheCell[2]<<endl;
            
            
            
            counter=0;
        }
        counter=counter+1;
        // ______________________________________________________________SAVING AND CALCULATION SECTION
        
        
        
        
        
        //d2_10_______________________________________dissimilarity
        if (chromatin_contant_matrix_calculation_flag==true) {
            // *********build cm0
            if(MD_Step<=StepCm0 & ReadCm0rfomFile==0)
            {
                for(int i1=0;i1<Chromatin_num_of_Beads;i1++)
                {
                    for(int j1=0;j1<Chromatin_num_of_Beads;j1++)
                    {
                        Cm0[i1][j1]=Cm0[i1][j1]+Cm[i1][j1];
                    }
                }
                if(MD_Step==StepCm0)
                {
                    for(int i1=0;i1<Chromatin_num_of_Beads;i1++)
                    {
                        for(int j1=0;j1<Chromatin_num_of_Beads;j1++)
                        {
                            Cm0[i1][j1]=Cm0[i1][j1]/(StepCm0+1); //normalizing
                        }
                    }
                    ofstream Cm0save;  //saving
                    Cm0save.open("results/Cm0.txt");
                    for(int i1=0;i1<Chromatin_num_of_Beads;i1++)
                    {
                        for(int j1=0;j1<Chromatin_num_of_Beads;j1++)
                        {
                            Cm0save  <<  Cm0[i1][j1] << "  ";
                        }
                        Cm0save  <<endl;
                    }
                }
                
            }
            else if (MD_Step==0 & ReadCm0rfomFile==1)
            {
                
                ifstream Cm0read;  //reading
                Cm0read.open("Cm0.txt");
                for(int i1=0;i1<Chromatin_num_of_Beads;i1++)
                {
                    for(int j1=0;j1<Chromatin_num_of_Beads;j1++)
                    {
                        Cm0read >> Cm0[i1][j1];
                    }
                }
            }
            
            // *********build cm0
            
            
            
            if(MD_Step%StepCm1 ==0 & MD_Step!=0)
            {
                
                d2_10=0.0;
                for(int i1=0;i1<Chromatin_num_of_Beads;i1++)
                {
                    for(int j1=0;j1<Chromatin_num_of_Beads;j1++)
                    {
                        Cm1[i1][j1]=Cm1[i1][j1]/StepCm1; //normalize contact matrix
                        d2_10=d2_10+( (Cm1[i1][j1]-Cm0[i1][j1]) )*( (Cm1[i1][j1]-Cm0[i1][j1]) ); //calculate distance
                        Cm1[i1][j1]=0.0; //refresh cm1
                    }
                }
                
                if (ReadCm0rfomFile==0 & MD_Step>StepCm0)
                {
                    d2_10ofstream<<  MD_Step <<"         "<<d2_10 << endl;
                }
                else if (ReadCm0rfomFile==1)
                {
                    d2_10ofstream<<  MD_Step <<"         "<<d2_10 << endl;
                }
                
            }
            //d2_10_______________________________________dissimilarity
            
            
            
            //d2_32_______________________________________dissimilarity
            
            
            
            if(MD_Step%StepCm2 ==0 & MD_Step!=0)
            {
                d2_32=0.0;
                ContactNumber=0;
                ContactProbablilityNumber=0;
                
                for(int i1=0;i1<Chromatin_num_of_Beads;i1++)
                {
                    for(int j1=0;j1<Chromatin_num_of_Beads;j1++)
                    {
                        if(Cm2[i1][j1]!=0.0 & i1!=j1 ) //consider the contacts in offdiagonal
                        {
                            ContactNumber=ContactNumber+1.0;
                            ContactProbablilityNumber=ContactProbablilityNumber+Cm2[i1][j1];
                        }
                        
                        
                        Cm3[i1][j1]=Cm3[i1][j1]/StepCm2; //normalize contact matrix
                        d2_32=d2_32+( (Cm3[i1][j1]-Cm2[i1][j1]) )*( (Cm3[i1][j1]-Cm2[i1][j1]) ); //calculate distance
                        Cm2[i1][j1]=Cm3[i1][j1]; //refresh cm2
                        Cm3[i1][j1]=0.0; //refresh cm3
                    }
                }
                
                if(MD_Step>StepCm2) //let StepCm2 to become calculated
                {
                    d2_32ofstream<<  MD_Step <<"         "<< "         "<<d2_32 << endl;  //save
                    REDUCEDd2_32ofstream <<" "<<VolumeNucleus(Membrane_Node_Position, Membrane_triangle_list, Membrane_Normal_direction, Outer_Membrane_num_of_triangles)
                    << " "<< ContactNumber   <<  " "  <<d2_32/ContactProbablilityNumber <<endl;
                    
                    
                    for(int i1=0;i1<Chromatin_num_of_Beads;i1++)
                    {
                        for(int j1=0;j1<Chromatin_num_of_Beads;j1++)
                        {
                            cm3ofstream << Cm2[i1][j1] << "   ";  //note thar cm3 is reset and  stored in cm2
                        }
                        cm3ofstream <<"  " <<endl;
                    }
                    cm3ofstream <<"  " <<endl;
                    cm3ofstream <<"  " <<endl;
                    cm3ofstream <<"  " <<endl;
                }
            }
            
            //d2_32_______________________________________dissimilarity
            
            
            
            
            
            //______________________________________contact matrix part__________________________________________
            if (MD_Step==(MD_num_of_steps-StepsCmLaststeps)  &  MD_Step!=0  ) // for last steps cacculate contact matrix
            {
                for(int i1=0;i1<Chromatin_num_of_Beads;i1++)
                {
                    for(int j1=0;j1<Chromatin_num_of_Beads;j1++)
                    {
                        CmLaststeps[i1][j1]=0.0;
                        Flag=1;
                        //                    cout <<"StepsCmLaststeps flaged" <<endl;
                    }
                }
            }
            //______________________________________contact matrix part__________________________________________
            
        }
        
        if (MD_Step%(export_povray_step_distance)==0    ) // for last steps cacculate contact matrix
        {
            povray_output_creator(MD_Step, Membrane_Node_Position, Membrane_triangle_list, Membrane_Normal_direction,  Membrane_Node_Pair_list, Actin_Node_Position, Actin_Node_Pair_List, Chromatin_Bead_Position,   ECM_Node_Position, ECM_Node_Pair_List, ECM_surface_triangle_list, Outer_Membrane_num_of_triangles, Membrane_num_of_Node_Pairs, Outer_Membrane_num_of_Node_Pairs, Actin_num_of_Bonds, ECM_num_of_Bonds, Membrane_num_of_Nodes);
        }
    }
    
    //==========================================================================================================================
    //==========================================================================================================================
    //================================================ End of the main MD loop =================================================
    //==========================================================================================================================
    //==========================================================================================================================
    
    
    //pov
    povray_output_creator( MD_num_of_steps, Membrane_Node_Position, Membrane_triangle_list,
                          Membrane_Normal_direction, Membrane_Node_Pair_list, Actin_Node_Position,Actin_Node_Pair_List, Chromatin_Bead_Position, ECM_Node_Position, ECM_Node_Pair_List, ECM_surface_triangle_list, Outer_Membrane_num_of_triangles, Membrane_num_of_Node_Pairs, Outer_Membrane_num_of_Node_Pairs, Actin_num_of_Bonds, ECM_num_of_Bonds, Membrane_num_of_Nodes);
    //pov
    
    
    //_____________________________spatial contact matrix:
    
    // x
    ofstream spatial_contact_matrix_x_ofstram;
    spatial_contact_matrix_x_ofstram.open("results/spatial_contact_matrix_x.txt");
    for(int i=0;i<size_spatialcontactmatrix;i++)
    {
        
        for(int j=0;j<size_spatialcontactmatrix;j++)
        {
            spatial_contact_matrix_x_ofstram << ((double)spatial_contact_matrix_x[i][j]) /((double) MD_num_of_steps) <<" ";
        }
        spatial_contact_matrix_x_ofstram <<endl;
    }
    
    
    // y
    ofstream spatial_contact_matrix_y_ofstram;
    spatial_contact_matrix_y_ofstram.open("results/spatial_contact_matrix_y.txt");
    for(int i=0;i<size_spatialcontactmatrix;i++)
    {
        
        for(int j=0;j<size_spatialcontactmatrix;j++)
        {
            spatial_contact_matrix_y_ofstram << ((double)spatial_contact_matrix_y[i][j]) /((double) MD_num_of_steps)  <<" ";
        }
        spatial_contact_matrix_y_ofstram <<endl;
    }
    
    
    // z
    ofstream spatial_contact_matrix_z_ofstram;
    spatial_contact_matrix_z_ofstram.open("results/spatial_contact_matrix_z.txt");
    for(int i=0;i<size_spatialcontactmatrix;i++)
    {
        
        for(int j=0;j<size_spatialcontactmatrix;j++)
        {
            spatial_contact_matrix_z_ofstram << ((double)spatial_contact_matrix_z[i][j]) /((double) MD_num_of_steps) <<" " ;
        }
        spatial_contact_matrix_z_ofstram <<endl;
    }
    
    
    //_____________________________spatial contact matrix:
    
    
    //_____________________________ contact matrix:
    ofstream contact;
    contact.open("results/contactLaststeps.txt");
    
    for(int i=0;i<Chromatin_num_of_Beads;i++)
    {
        for(int j=0;j<Chromatin_num_of_Beads;j++)
        {
            if(Flag==0)
            {
                contact<< CmLaststeps[i][j]/MD_num_of_steps<<" ";
            }
            else if(Flag==1)
            {
                contact<< ((double)CmLaststeps[i][j])/((double)StepsCmLaststeps)<<" ";
            }
        }
        contact<<endl;
    }
    
    
    //_______________________________restart:
    
    cout<<"Saving..."<<endl;
    restartsave(Membrane_Node_Position, Membrane_Node_Velocity, Membrane_triangle_list, Membrane_Normal_direction, Membrane_Triangle_Pair_Nodes, Membrane_Node_Pair_list, Actin_Node_Position, Actin_Node_Velocity, Actin_Node_Pair_List, Chromatin_Bead_Position, Chromatin_Bead_Velocity, ECM_Node_Position, ECM_Node_Velocity, Membrane_num_of_Triangle_Pairs, Actin_num_of_Bonds, Membrane_num_of_Nodes);
    
    
    
    //_____________________________ topology
    ofstream topo1;
    topo1.open("results/topology.txt");
    
    topo1<< "package require topotools ;"<<endl;
    topo1 <<"topo clearbonds ; "<<"  "<<endl;
    int alpha = 0,beta = 0,ettha,gammmma;
    for(int i=0;i<Membrane_num_of_Node_Pairs;i++)
    {
        ettha=Membrane_Node_Pair_list[i][0];
        gammmma=Membrane_Node_Pair_list[i][1];
        
        for(int j=0; j<Outer_Membrane_num_of_Nodes;j++) // saving trajectory
        {
            if(Outer_Membrane_list_of_Nodes[j]==ettha)
            {
                alpha=j;
            }
            if(Outer_Membrane_list_of_Nodes[j]==gammmma)
            {
                beta=j;
            }
        }
        
        for(int j=0; j<Membrane_num_of_Nodes- Outer_Membrane_num_of_Nodes;j++) // saving trajectory
        {
            if(Nucleus_Membrane_list_of_Nodes[j]==ettha)
            {
                alpha=j+Outer_Membrane_num_of_Nodes;
            }
            if(Nucleus_Membrane_list_of_Nodes[j]==gammmma)
            {
                beta=j+Outer_Membrane_num_of_Nodes;
            }
        }
        
        topo1<< "topo addbond "<< alpha<< "  "<<beta<<" ;"; // add bond for membane
        
        
        
    }
    
    for (int nchain=0;nchain<Chromatin_num_of_chains;nchain++)
    {
        for(int i=nchain*(Chromatin_num_of_Beads/Chromatin_num_of_chains)  ;i< (nchain+1)*(Chromatin_num_of_Beads/Chromatin_num_of_chains) -1; i++ )  // all beads interaction whit the next one
        {
            topo1 <<  "topo addbond "  << i+Membrane_num_of_Nodes << "   " << i+1+Membrane_num_of_Nodes <<"    " << ";"  ;
        }
    }
    
    
    
    
    int w1,w2;
    for(int i=0;i<Actin_num_of_Bonds;i++)
    {
        
        w1=(int)  Actin_Node_Pair_List[i][0];
        w2=(int)  Actin_Node_Pair_List[i][1];
        
        if( (w1<Actin_Membrane_shared_num_of_Nodes & w2>Actin_Membrane_shared_num_of_Nodes) ||  (w2<Actin_Membrane_shared_num_of_Nodes & w1>Actin_Membrane_shared_num_of_Nodes)  )
        {
            // topo1<< "topo addbond "<< (int)  Actin_Node_Pair_List[i][0] +Membrane_num_of_Nodes << "  "<<( int) Actin_Node_Pair_List[i][1]+Membrane_num_of_Nodes<<" ;";
        }
    }
    
    
    
    
    
    
    for(int i=0;i<ECM_num_of_Bonds;i++)
    {
        w1=(int)  ECM_Node_Pair_List[i][0];
        w2=(int)  ECM_Node_Pair_List[i][1];
        
        topo1<< "topo addbond "<< w1 +Membrane_num_of_Nodes + Chromatin_num_of_Beads << "  "<<w2+Membrane_num_of_Nodes + Chromatin_num_of_Beads <<" ;";
    }
    
    topo1 << " draw delete all " <<endl;
    topo1 << " draw color red " <<endl;
    
    //    for(int i=0;i<ECM_Surface_num_of_Triangles;i++)
    //    {
    // topo1 << " draw triangle {"<< ECM_Node_Position[ ECM_surface_triangle_list[i][0] ] [0]  << "  "<< ECM_Node_Position[ ECM_surface_triangle_list[i][0] ] [1] << "  "<< ECM_Node_Position[ ECM_surface_triangle_list[i][0] ] [2]  << "}";
    // topo1 << " {"<< ECM_Node_Position[ ECM_surface_triangle_list[i][1] ] [0]  << "  "<< ECM_Node_Position[ ECM_surface_triangle_list[i][1] ] [1] << "  "<< ECM_Node_Position[ ECM_surface_triangle_list[i][1] ] [2]  << "}";
    // topo1 << " {"<< ECM_Node_Position[ ECM_surface_triangle_list[i][2] ] [0]  << "  "<< ECM_Node_Position[ ECM_surface_triangle_list[i][2] ] [1] << "  "<< ECM_Node_Position[ ECM_surface_triangle_list[i][2] ] [2]  << "}" <<endl;
    //    }
    cout<<"done"<<endl;
    printf("Time taken: %.2fminutes\n", (double)((clock() - tStart)/CLOCKS_PER_SEC)/60.0);
}

// ==============================================================================================================
// ============================================== END OF MAIN ===================================================
// ==============================================================================================================



///______________________basic functions:






void thermostat(double (&Chromatin_Bead_Velocity)[Chromatin_num_of_Beads][3] )//rescaling V to reach true temprature
{
    
    double alpha,K,change;
    K= kineticenergychromatin(Chromatin_Bead_Velocity);
    alpha=Chromatin_num_of_Beads*KT/K;
    
    change=sqrt(alpha);
    for(int i=0 ; i<Chromatin_num_of_Beads ;i++)
    {
        
        Chromatin_Bead_Velocity[i][0] *= change;
        Chromatin_Bead_Velocity[i][1] *= change;
        Chromatin_Bead_Velocity[i][2] *= change;
        
    }
    
}





//______________________membrane functions

//Also the potential Energy is calculated here
void Membrane_Force_Calculator (double Membrane_Node_Position[][3],double Membrane_Node_Velocity[][3],double Membrane_Node_Force [][3],int Membrane_Node_Pair_list[][2],int Membrane_Triangle_Pair_Nodes[][4],double &Total_Potential_Energy, int Membrane_num_of_Triangle_Pairs, int Membrane_num_of_Node_Pairs)
{
    double le0,le1,lmax,lmin;
    double deltax,deltay,deltaz,temp_Node_distance,temp_force;
    int pos1,pos2,pos3,pos4;  // to making calculation of surface force easier
    double temp_potential_energy = 0.0;
    double xpos[3][5],N1[3],N2[3], N3[3],p3p1[3],p3p2[3],p4p2[3],p4p1[3],Ep2p1[3],sinus,F0,F1[3],F2[3],F3[3],F4[3];// for exmple p3p1 is p3-p1 and ....
    
    /// calculate network force:
    int temp_Node_A, temp_Node_B;
    le0=1.15000*Node_radius;
    lmax=1.33000*Node_radius;
    le1=0.85000*Node_radius;
    lmin=0.67000*Node_radius;
    Total_Potential_Energy=0.0;
    
    for (int k=0 ; k< Membrane_num_of_Node_Pairs ; k++)
    {
        temp_Node_B=Membrane_Node_Pair_list[k][0];
        temp_Node_A=Membrane_Node_Pair_list[k][1];
        
        deltax=Membrane_Node_Position[temp_Node_A][0]-Membrane_Node_Position[temp_Node_B][0];
        deltay=Membrane_Node_Position[temp_Node_A][1]-Membrane_Node_Position[temp_Node_B][1];
        deltaz=Membrane_Node_Position[temp_Node_A][2]-Membrane_Node_Position[temp_Node_B][2];
        
        
        if (Periodic_condtion_status == 1.0) {
            deltax=periodiccondition(deltax);
            deltay=periodiccondition(deltay);
        }
        
        temp_Node_distance=sqrt(deltax*deltax+deltay*deltay+deltaz*deltaz);
        temp_force=0.0;
        double temp_exp_le0=exp(1.0/(le0-temp_Node_distance));
        double temp_exp_le1=exp(1.0/(temp_Node_distance-le1));
        //*******************************************************************************************************
        /*BUG
         |---\   |    |  /---\
         |    |  |    |  |
         |---<   |    |  |  -\
         |    |  |    |  |   |
         |---/   \----/  \---/
         */
        //*******************************************************************************************************
        //***************** Potential BUG: F=-dU/dr but in many cases I cannot determin wheather ****************
        //***************** the '-' has been implemented or not. Since the potential energy is   ****************
        //***************** never used in the code it does not a threat. ****************************************
        //*******************************************************************************************************
        
        if(temp_Node_distance >le1  & temp_Node_distance < le0 )  //free zone
        {
            temp_potential_energy=0 ; // free zone
        }
        
        if(temp_Node_distance > le0  & temp_Node_distance <lmax )  //bondforce
        {
            temp_force = (Membrane_spring_coefficient*temp_exp_le0/(lmax-temp_Node_distance))*( 1.0/(lmax-temp_Node_distance) +  1.0/((le0-temp_Node_distance)*(le0-temp_Node_distance)));
            temp_potential_energy= Membrane_spring_coefficient*temp_exp_le0/(lmax-temp_Node_distance);
            
        }
        
        if(temp_Node_distance < le1   &  temp_Node_distance > lmin  )  // repulsive force
        {
            temp_force= -(Membrane_spring_coefficient*temp_exp_le1/(temp_Node_distance-lmin))*( 1.0/(temp_Node_distance-lmin) + 1.0/((temp_Node_distance-le1)*(temp_Node_distance-le1)));                 // force on i th from j
            temp_potential_energy= Membrane_spring_coefficient*temp_exp_le1/(temp_Node_distance-lmin);
        }
        /// my cutoff for force amplitute and for avoiding leting particle scape from force trap
        if(temp_force>965.31  || temp_Node_distance>lmax )
        {
            temp_force = 965.31+Membrane_spring_force_cutt_off* ( temp_Node_distance - 1.3280*Node_radius );
            temp_potential_energy=   1.81599  + 965.31 * ( temp_Node_distance - 1.3280*Node_radius )+0.5*Membrane_spring_force_cutt_off * ( temp_Node_distance - 1.3280*Node_radius ) * ( temp_Node_distance - 1.3280*Node_radius );
        }
        
        
        if(temp_force<-1000.05   ||  temp_Node_distance<lmin )
        {
            temp_force =-1000.05-Membrane_spring_force_cutt_off* ( 0.671965*Node_radius - temp_Node_distance );
            temp_potential_energy = 1.85038 + 1005.05 * ( 0.671965*Node_radius - temp_Node_distance )+0.5*Membrane_spring_force_cutt_off*( 0.671965*Node_radius - temp_Node_distance )*( 0.671965*Node_radius - temp_Node_distance );
        }
        
        Total_Potential_Energy += temp_potential_energy;
        
        // implimentation of forces:
        Membrane_Node_Force[temp_Node_A][0] += temp_force*deltax/temp_Node_distance+membrane_damping_coefficient*(Membrane_Node_Velocity[temp_Node_A][0]-Membrane_Node_Velocity[temp_Node_B][0]);
        Membrane_Node_Force[temp_Node_A][1] += temp_force*deltay/temp_Node_distance+membrane_damping_coefficient*(Membrane_Node_Velocity[temp_Node_A][1]-Membrane_Node_Velocity[temp_Node_B][1]);
        Membrane_Node_Force[temp_Node_A][2] += temp_force*deltaz/temp_Node_distance+membrane_damping_coefficient*(Membrane_Node_Velocity[temp_Node_A][2]-Membrane_Node_Velocity[temp_Node_B][2]);
        
        Membrane_Node_Force[temp_Node_B][0] += -temp_force*deltax/temp_Node_distance-membrane_damping_coefficient*(Membrane_Node_Velocity[temp_Node_A][0]-Membrane_Node_Velocity[temp_Node_B][0]); //from j  to i
        Membrane_Node_Force[temp_Node_B][1] += -temp_force*deltay/temp_Node_distance-membrane_damping_coefficient*(Membrane_Node_Velocity[temp_Node_A][1]-Membrane_Node_Velocity[temp_Node_B][1]);
        Membrane_Node_Force[temp_Node_B][2] += -temp_force*deltaz/temp_Node_distance-membrane_damping_coefficient*(Membrane_Node_Velocity[temp_Node_A][2]-Membrane_Node_Velocity[temp_Node_B][2]);
    }
    // End of Membrane Node Pair forces
    // Beginning of the  triangle-triangle (bending) force calculations
    for(int i=0 ;i<Membrane_num_of_Triangle_Pairs;i++)  // who are neighbors?
    {
        
        pos1=Membrane_Triangle_Pair_Nodes[i][0];
        pos2=Membrane_Triangle_Pair_Nodes[i][1];
        pos3=Membrane_Triangle_Pair_Nodes[i][2];
        pos4=Membrane_Triangle_Pair_Nodes[i][3];
        
        for(int i2=0;i2<3;i2++)   // storing coordiates of poses in new array in order to imposing periodic boundary
        {
            
            xpos[i2][1]=Membrane_Node_Position[pos1][i2];
            xpos[i2][2]=Membrane_Node_Position[pos2][i2];
            xpos[i2][3]=Membrane_Node_Position[pos3][i2];
            xpos[i2][4]=Membrane_Node_Position[pos4][i2];
        }
        
        for(int i2=0;i2<2;i2++)   // puting pos1 in (0,0,?)  --- only x and y  ,z remains untouch
        {
            
            xpos[i2][2] += -xpos[i2][1];
            xpos[i2][3] += -xpos[i2][1];
            xpos[i2][4] += -xpos[i2][1];
            xpos[i2][1]=0.0;
        }
        if (Periodic_condtion_status == 1.0 ) {
            if(  xpos[0][2]> 6.0*Node_radius  )
                xpos[0][2]= xpos[0][2]-(Lbox+1.0);
            
            if(  xpos[1][2]> 6.0*Node_radius  )
                xpos[1][2]= xpos[1][2]-(Lbox+1.0);
            
            if(  xpos[0][2]< -6.0*Node_radius  )
                xpos[0][2]= xpos[0][2]+(Lbox+1.0);
            
            if(  xpos[1][2]< -6.0*Node_radius  )
                xpos[1][2]= xpos[1][2]+(Lbox+1.0);
            
            
            
            if(  xpos[0][3]> 6.0*Node_radius  )
                xpos[0][3]= xpos[0][3]-(Lbox+1.0);
            
            if(  xpos[1][3]> 6.0*Node_radius  )
                xpos[1][3]= xpos[1][3]-(Lbox+1.0);
            
            if(  xpos[0][3]< -6.0*Node_radius  )
                xpos[0][3]= +xpos[0][3]+(Lbox+1.0);
            
            if(  xpos[1][3]<- 6.0*Node_radius  )
                xpos[1][3]= +xpos[1][3]+(Lbox+1.0);
            
            
            if(  xpos[0][4]> 6.0*Node_radius  )
                xpos[0][4]= xpos[0][4]-(Lbox+1.0);
            
            if(  xpos[1][4]> 6.0*Node_radius  )
                xpos[1][4]= xpos[1][4]-(Lbox+1.0);
            
            if(  xpos[0][4]< -6.0*Node_radius  )
                xpos[0][4]= +xpos[0][4]+(Lbox+1.0);
            
            if(  xpos[1][4]<- 6.0*Node_radius  )
                xpos[1][4]= +xpos[1][4]+(Lbox+1.0);
            ///PERIODIC
        }
        
        
        
        ///  FORCE CACCULATION:
        
        //****************************** pre calculation:
        p3p1[0]=xpos[0][3]-xpos[0][1];
        p3p1[1]=xpos[1][3]-xpos[1][1];
        p3p1[2]=xpos[2][3]-xpos[2][1];
        
        p3p2[0]=xpos[0][3]-xpos[0][2];
        p3p2[1]=xpos[1][3]-xpos[1][2];
        p3p2[2]=xpos[2][3]-xpos[2][2];
        
        p4p2[0]=xpos[0][4]-xpos[0][2];
        p4p2[1]=xpos[1][4]-xpos[1][2];
        p4p2[2]=xpos[2][4]-xpos[2][2];
        
        p4p1[0]=xpos[0][4]-xpos[0][1];
        p4p1[1]=xpos[1][4]-xpos[1][1];
        p4p1[2]=xpos[2][4]-xpos[2][1];
        
        Ep2p1[0]=xpos[0][2]-xpos[0][1];
        Ep2p1[1]=xpos[1][2]-xpos[1][1];
        Ep2p1[2]=xpos[2][2]-xpos[2][1];
        
        
        crossvector(N1,p3p1,p3p2);
        crossvector(N2,p4p2,p4p1);
        crossvector(N3,N2,N1);
        sinus=vectorlength(N3)/(vectorlength(N2)*vectorlength(N1));
        F0 = Membrane_bending_coefficient*sinus;
        //        cout<<"\nF0="<<F0<<endl;
        if( parallelORantiparallel(xpos ) == +1 )
        {
            F0=-F0;
        }
        double temp_N1_length_squared=vectorlength(N1)*vectorlength(N1);
        double temp_N2_length_squared=vectorlength(N2)*vectorlength(N2);
        double temp_Ep2p1_length=vectorlength(Ep2p1);
        //**************************************************** force calculation
        for (int l=0; l<3; l++) {
            F3[l]= F0 * temp_Ep2p1_length* N1[l]/ temp_N1_length_squared;
            
            F4[l]= F0 * temp_Ep2p1_length* N2[l]/ temp_N2_length_squared;
            
            F1[l]= (F0/temp_Ep2p1_length)*( innerproduct(p3p2,Ep2p1)*N1[l]/temp_N1_length_squared + innerproduct(p4p2,Ep2p1)*N2[l]/temp_N2_length_squared );
            
            F2[l]= (-F0/temp_Ep2p1_length)*( innerproduct(p3p1,Ep2p1)*N1[l]/temp_N1_length_squared + innerproduct(p4p1,Ep2p1)*N2[l]/temp_N2_length_squared );
            
            
            Membrane_Node_Force[pos1][l] += F1[l];
            Membrane_Node_Force[pos2][l] += F2[l];
            Membrane_Node_Force[pos3][l] += F3[l];
            Membrane_Node_Force[pos4][l] += F4[l];
            
        }
        
        //*******************************************************************************************************
        /*BUG
         |---\   |    |  /---\
         |    |  |    |  |
         |---<   |    |  |  -\
         |    |  |    |  |   |
         |---/   \----/  \---/
         */
        //*******************************************************************************************************
        //***************** Potential BUG: Not very sure about the calculations *****************************
        //***************** for 'Total_Potential_Energy. Need to double check *********************************
        //*******************************************************************************************************
        
        Total_Potential_Energy += Membrane_bending_coefficient*(1.0 -  ( innerproduct(N1,N2)/(vectorlength(N1)*vectorlength(N2)  )   ) );
    }  // end of 'for-i'
    //    exit (EXIT_FAILURE);
}///end of function

void Membrane_Actin_shared_Node_Force_calculator (double Membrane_Node_Position[][3], double  Actin_Node_Position[][3], double Membrane_Node_Force[][3], double Actin_Node_Force[][3],  int Membrane_Actin_shared_Node_list[][2], double Membrane_Node_velocity[][3], double Actin_Node_velocity[Actin_num_of_Nodes][3])
{
    double delta_x,delta_y,delta_z,temp_Node_distance,temp_force;
    
    int temp_Node_mem,temp_Node_act;
    
    for (int act_mem_temp_node=0 ; act_mem_temp_node< Actin_Membrane_shared_num_of_Nodes ; act_mem_temp_node++)
    {
        temp_Node_mem = Membrane_Actin_shared_Node_list[act_mem_temp_node][0];
        temp_Node_act = Membrane_Actin_shared_Node_list[act_mem_temp_node][1];
        
        delta_x = Membrane_Node_Position[temp_Node_mem][0]-Actin_Node_Position[temp_Node_act][0];
        delta_y = Membrane_Node_Position[temp_Node_mem][1]-Actin_Node_Position[temp_Node_act][1];
        delta_z = Membrane_Node_Position[temp_Node_mem][2]-Actin_Node_Position[temp_Node_act][2];
        
        if (Periodic_condtion_status == 1.0) {
            delta_x=periodiccondition(delta_x);
            delta_y=periodiccondition(delta_y);
        }
        
        temp_Node_distance=sqrt(delta_x*delta_x + delta_y*delta_y + delta_z*delta_z);
        temp_force=0.0;
        
        if (temp_Node_distance > 0.0001 || temp_Node_distance < -0.0001) {
            temp_force = Actin_membrane_stiff_spring_coefficient*temp_Node_distance;
            
            Membrane_Node_Force[temp_Node_mem][0] += temp_force*delta_x/temp_Node_distance;
            Membrane_Node_Force[temp_Node_mem][1] += temp_force*delta_y/temp_Node_distance;
            Membrane_Node_Force[temp_Node_mem][2] += temp_force*delta_z/temp_Node_distance;
            
            Actin_Node_Force[temp_Node_act][0] += -temp_force*delta_x/temp_Node_distance;
            Actin_Node_Force[temp_Node_act][1] += -temp_force*delta_y/temp_Node_distance;
            Actin_Node_Force[temp_Node_act][2] += -temp_force*delta_z/temp_Node_distance;
        }
        Membrane_Node_Force[temp_Node_mem][0] += Actin_membrane_damping_coefficient*(Membrane_Node_velocity[temp_Node_mem][0]-Actin_Node_velocity[temp_Node_act][0]);
        Membrane_Node_Force[temp_Node_mem][1] += Actin_membrane_damping_coefficient*(Membrane_Node_velocity[temp_Node_mem][1]-Actin_Node_velocity[temp_Node_act][1]);
        Membrane_Node_Force[temp_Node_mem][2] += Actin_membrane_damping_coefficient*(Membrane_Node_velocity[temp_Node_mem][2]-Actin_Node_velocity[temp_Node_act][2]);
        Actin_Node_Force[temp_Node_act][0] += -Actin_membrane_damping_coefficient*(Membrane_Node_velocity[temp_Node_mem][0]-Actin_Node_velocity[temp_Node_act][0]);
        Actin_Node_Force[temp_Node_act][1] += -Actin_membrane_damping_coefficient*(Membrane_Node_velocity[temp_Node_mem][1]-Actin_Node_velocity[temp_Node_act][1]);
        Actin_Node_Force[temp_Node_act][2] += -Actin_membrane_damping_coefficient*(Membrane_Node_velocity[temp_Node_mem][2]-Actin_Node_velocity[temp_Node_act][2]);
        
    }
    //    exit (EXIT_FAILURE);
}

int parallelORantiparallel( double xpos[3][5] )// + for parallel  -for anti parallel (used in force cacculator of membrane)
{
    
    //************************************ Parallel or anti parallel calculator
    double deltax,deltay;
    double temp1[3],temp3[3],temp4[3],normmembrane[3]; //used to determine if the Membrane_num_of_Nodes and (p1-p2) are not parallel
    int sign;
    
    // calculating n1 and n2   :  normal to first and second triangle
    temp1[0]=(xpos[0][1]+xpos[0][2])/2.0;// Now i'm using temp1 to fond the middle POINT of the edge line of two triangle
    temp1[1]=(xpos[1][1]+xpos[1][2])/2.0;// Now i'm using temp1 to fond the middle POINT of the edge line of two triangle
    temp1[2]=(xpos[2][1]+xpos[2][2])/2.0;// Now i'm using temp1 to fond the middle POINT of the edge line of two triangl2
    
    temp3[0]=xpos[0][3]-temp1[0];// (m1 in doc1) this finds  vector between the midle of the edge and pos3
    temp3[1]=xpos[1][3]-temp1[1];// this finds  vector between the midle of the edge and pos3
    temp3[2]=xpos[2][3]-temp1[2];// this finds  vector between the midle of the edge and pos3
    
    temp4[0]=xpos[0][4]-temp1[0];// (m2 in doc1)this finds  vector between the midle of the edge and pos4
    temp4[1]=xpos[1][4]-temp1[1];// this finds  vector between the midle of the edge and pos4
    temp4[2]=xpos[2][4]-temp1[2];// this finds  vector between the midle of the edge and pos4
    
    temp1[0]=(xpos[0][1]-xpos[0][2]);// (P1-P2) Now i'm using temp1 to fond the  vector between pos1 and pos2  (this is useful in finding  Membrane_num_of_Nodes)
    temp1[1]=(xpos[1][1]-xpos[1][2]);// Now i'm using temp1 to fond the  vector between pos1 and pos2  (this is useful in finding  Membrane_num_of_Nodes)
    temp1[2]=(xpos[2][1]-xpos[2][2]);// Now i'm using temp1 to fond the  vector between pos1 and pos2  (this is useful in finding  Membrane_num_of_Nodes)
    
    deltax=innerproduct(temp1,temp3);// im usind deltax and deltay to store innerporduct of edge and m1 or edge nd m2: by this i can find M1 and M2
    deltay=innerproduct(temp1,temp4);// im usind deltax and deltay to store innerporduct of edge and m1 or edge nd m2: by this i can find M1 and M2
    
    temp3[0]=temp3[0]-deltax*temp1[0];// (M1 in doc1) this finds  vector between the midle of the edge and pos3
    temp3[1]=temp3[1]-deltax*temp1[1];// this finds  vector between the midle of the edge and pos3
    temp3[2]=temp3[2]-deltax*temp1[2];// this finds  vector between the midle of the edge and pos3
    
    temp4[0]=temp4[0]-deltay*temp1[0];// (M2 in doc1)this finds  vector between the midle of the edge and pos4
    temp4[1]=temp4[1]-deltay*temp1[1];// this finds  vector between the midle of the edge and pos4
    temp4[2]=temp4[2]-deltay*temp1[2];// this finds  vector between the midle of the edge and pos4
    
    
    // IF Membrane_num_of_Nodes||P1-P2  THEN n1 and n2 are  -(P1-p2)xM1 and  (P1-p2)xM2 .
    // But if Membrane_num_of_Nodes is anti paraller to P1-P2 THEN n1 and n2 are  (P1-p2)xM1 and  -(P1-p2)xM2
    crossvector(normmembrane,temp3,temp4);//(Membrane_num_of_Nodes in doc1) now i am using  for diffrent purpose  : storing the vector of edge of two triangles
    if(innerproduct(normmembrane,temp1)<0.0)
    {
        sign=-1;
    }
    else
    {
        sign=+1;
    }
    return sign;
    
    
    
    //************************************ Parallel or anti parallel calculator  calculator
}

void MonteCarlo(  int Membrane_triangle_list[Membrane_num_of_Triangles][3],int Membrane_Triangle_Pair_Nodes[][4] ,int Membrane_Node_Pair_list[][2],double &Total_Potential_Energy,double Membrane_Node_Position[][3], int Membrane_Normal_direction[Membrane_num_of_Triangles][2], int Outer_Membrane_num_of_triangles, int Nucleus_Membrane_num_of_triangles, int Membrane_num_of_Triangle_Pairs, int Membrane_num_of_Node_Pairs)
{
    int t1=-1,t2=-1,t3=-1,t4=-1,l1=-1,l2=-1,l3=-1; // just numbers
    
    int  linepos,lineA,lineB,lineC,lineD,Aa,Bb,Cc,Dd ; // IN ORDER TO REVERSE CHANGES IF STEP NOT ACCEPTED IN trapezium fuction
    
    //    Bug: I commented out this part because it just goes through the same calculations as we had in the 'Membrane_Force_Calculator' function.
    Total_Potential_Energy=potentialenergy(Membrane_Node_Position, Membrane_Node_Pair_list, Membrane_Triangle_Pair_Nodes, Membrane_triangle_list, Outer_Membrane_num_of_triangles, Nucleus_Membrane_num_of_triangles, Membrane_num_of_Triangle_Pairs, Membrane_num_of_Node_Pairs);
    
    double temp_total_potential_energy=0.0, Metropolis_condition=0.0;
    //C++11
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<> rd1(0,Membrane_num_of_Triangle_Pairs-1);
    uniform_int_distribution<> rd2(1,3);
    
    
    /// selecting a bond randomly  (bond is supposed to be between t1 and t2 )  and flippping bonds
    l1=rd1(gen);
    t1=Membrane_Triangle_Pair_Nodes[l1][0];
    t2=Membrane_Triangle_Pair_Nodes[l1][1];
    t3=Membrane_Triangle_Pair_Nodes[l1][2];
    t4=Membrane_Triangle_Pair_Nodes[l1][3];
    
    
    
    
    
    //      -----------------------------------------update trapezium
    trapeziumUpdateMonteCarlo(Membrane_Triangle_Pair_Nodes,t1,t2,t3,t4,linepos,lineA,lineB,lineC,lineD,Aa,Bb,Cc,Dd, Membrane_num_of_Triangle_Pairs);
    if( Aa!=Bb & Aa!=Cc & Aa!=Dd & Bb!=Cc & Bb!=Dd & Cc!=Dd  ) //updaing condition
    {
        for(int i=0;i<Membrane_num_of_Node_Pairs;i++)
        {
            if ((Membrane_Node_Pair_list[i][0]==t1 & Membrane_Node_Pair_list[i][1]==t2)  || (Membrane_Node_Pair_list[i][0]==t2 & Membrane_Node_Pair_list[i][1]==t1))
            {
                Membrane_Node_Pair_list[i][0]=t3;
                Membrane_Node_Pair_list[i][1]=t4;
                l3=i;
                temp_total_potential_energy=potentialenergy(Membrane_Node_Position, Membrane_Node_Pair_list, Membrane_Triangle_Pair_Nodes, Membrane_triangle_list, Outer_Membrane_num_of_triangles, Nucleus_Membrane_num_of_triangles, Membrane_num_of_Triangle_Pairs, Membrane_num_of_Node_Pairs);
                // cout<<"change  "<<endl;
                
            }
        }
    }
    else
    {
        temp_total_potential_energy=10000000000;
    }
    ///    END   -----------------------------------------update trapezium
    
    //C++11
    uniform_real_distribution<> rd3(0,1);
    
    Metropolis_condition=( exp(   (-temp_total_potential_energy+Total_Potential_Energy)/KT ));
    if (    rd3(gen) > Metropolis_condition   )  // if the changes are not accepted  ==== reverse the changes
    {
        //        I have to check to see if I do this the next montecarlo step will have bugs or not.
        //        Total_Potential_Energy=temp_total_potential_energy;
        if( Aa!=Bb & Aa!=Cc & Aa!=Dd & Bb!=Cc & Bb!=Dd & Cc!=Dd  )
        {
            
            // if(linepos=!-1)
            {
                
                Membrane_Node_Pair_list[l3][0]=t1;   //revesing Bonds
                Membrane_Node_Pair_list[l3][1]=t2;
                Membrane_Triangle_Pair_Nodes[linepos][0]=t1;  //revesing trapezium
                Membrane_Triangle_Pair_Nodes[linepos][1]=t2;
                Membrane_Triangle_Pair_Nodes[linepos][2]=t3;
                Membrane_Triangle_Pair_Nodes[linepos][3]=t4;
                
            }
            // if( lineA=!-1)
            {
                Membrane_Triangle_Pair_Nodes[lineA][3]=t2;
                Membrane_Triangle_Pair_Nodes[lineA][2]=Aa;
            }
            
            
            //if( lineB=!-1)
            {
                Membrane_Triangle_Pair_Nodes[lineB][3]=t2;
                Membrane_Triangle_Pair_Nodes[lineB][2]=Bb;
            }
            
            
            
            
            // if( lineC=!-1)
            {
                Membrane_Triangle_Pair_Nodes[lineC][3]=t1;
                Membrane_Triangle_Pair_Nodes[lineC][2]=Cc;
            }
            
            
            //  if( lineD=!-1)
            {
                
                Membrane_Triangle_Pair_Nodes[lineD][3]=t1;
                Membrane_Triangle_Pair_Nodes[lineD][2]=Dd;
            }
            
        }
        
    }
    else  // if changes are acepted we update trianglesa and normalvector inside this function:
    {
        updatetriangle(t1,t2,t3,t4,Membrane_triangle_list,l1,l2,Membrane_Node_Position,Membrane_Normal_direction);
        //   cout << "accept  "<< -U2+U <<endl;
    }
    
}

void trapeziumUpdateMonteCarlo(int Membrane_Triangle_Pair_Nodes[][4],int p1, int p2,int p3,int p4, int &linepos,int &lineA,int &lineB,int &lineC,int &lineD,int &Aa,int &Bb,int &Cc,int &Dd, int Membrane_num_of_Triangle_Pairs)
{
    
    linepos=-1;
    lineA=-1;
    lineB=-1;
    lineC=-1;
    lineD=-1;
    Aa=-1;
    Bb=-1;
    Cc=-1;
    Dd=-1;  // IN ORDER TO REVERSE CHANGES IF STEP NOT
    int findpos=0,finda=0,findb=0,findc=0,findd=0;
    
    
    int counter;
    
    counter=0;
    // find  relevant trapezium
    for(int i=0;i<Membrane_num_of_Triangle_Pairs ;i++)
    {
        if(  (Membrane_Triangle_Pair_Nodes[i][0]==p1 &  Membrane_Triangle_Pair_Nodes[i][1]==p2 ) || (Membrane_Triangle_Pair_Nodes[i][0]==p2 &  Membrane_Triangle_Pair_Nodes[i][1]==p1 )    )
        {
            if(  (Membrane_Triangle_Pair_Nodes[i][2]==p3 &  Membrane_Triangle_Pair_Nodes[i][3]==p4 ) || (Membrane_Triangle_Pair_Nodes[i][2]==p4 &  Membrane_Triangle_Pair_Nodes[i][3]==p3 )    )
            {
                linepos=i;
                findpos=+1;/////
                counter=counter+1;
            }
        }
    }
    
    //flipping the bond
    
    for(int i=0;i<Membrane_num_of_Triangle_Pairs ;i++)
    {
        
        
        // find A f
        if(  (Membrane_Triangle_Pair_Nodes[i][0]==p1 &  Membrane_Triangle_Pair_Nodes[i][1]==p4 ) || (Membrane_Triangle_Pair_Nodes[i][0]==p4 &  Membrane_Triangle_Pair_Nodes[i][1]==p1 )    )
        {
            if(  Membrane_Triangle_Pair_Nodes[i][2]==p2  )
            {
                Aa= Membrane_Triangle_Pair_Nodes[i][3];//     Membrane_Triangle_Pair_Nodes[i][2]=p3;
                lineA=i;
                counter=counter+1;
                finda=+1;
            }
            else  if(  Membrane_Triangle_Pair_Nodes[i][3]==p2 )
            {
                
                Aa= Membrane_Triangle_Pair_Nodes[i][2];//    Membrane_Triangle_Pair_Nodes[i][3]=p3;
                lineA=i;
                counter=counter+1;
                finda=+1;
            }
            
            // else
            //  cout << " A ";
        }
        
        
        // find    B
        if(  (Membrane_Triangle_Pair_Nodes[i][0]==p1 &  Membrane_Triangle_Pair_Nodes[i][1]==p3 ) || (Membrane_Triangle_Pair_Nodes[i][0]==p3 &  Membrane_Triangle_Pair_Nodes[i][1]==p1 )    )
        {
            if(  Membrane_Triangle_Pair_Nodes[i][2]==p2  )
            {
                Bb= Membrane_Triangle_Pair_Nodes[i][3];//       Membrane_Triangle_Pair_Nodes[i][2]=p4;
                lineB=i;
                counter=counter+1;
                findb=+1;
            }
            else if(  Membrane_Triangle_Pair_Nodes[i][3]==p2 )
            {
                
                Bb= Membrane_Triangle_Pair_Nodes[i][2];//       Membrane_Triangle_Pair_Nodes[i][3]=p4;
                lineB=i;
                counter=counter+1;
                findb=+1;
            }
            
            //         else
            //    cout << " B ";
        }
        
        
        // find    C
        if(  (Membrane_Triangle_Pair_Nodes[i][0]==p3 &  Membrane_Triangle_Pair_Nodes[i][1]==p2 ) || (Membrane_Triangle_Pair_Nodes[i][0]==p2 &  Membrane_Triangle_Pair_Nodes[i][1]==p3 )    )
        {
            if(  Membrane_Triangle_Pair_Nodes[i][2]==p1  )
            {
                
                Cc= Membrane_Triangle_Pair_Nodes[i][3];//                 Membrane_Triangle_Pair_Nodes[i][2]=p4;
                lineC=i;
                counter=counter+1;
                findc=+1;
            }
            else  if(  Membrane_Triangle_Pair_Nodes[i][3]==p1 )
            {
                Cc= Membrane_Triangle_Pair_Nodes[i][2];//                 Membrane_Triangle_Pair_Nodes[i][3]=p4;
                lineC=i;
                counter=counter+1;
                findc=+1;
            }
            //             else
            //   cout << " C ";
        }
        
        
        
        // find    D
        if(  (Membrane_Triangle_Pair_Nodes[i][0]==p2 &  Membrane_Triangle_Pair_Nodes[i][1]==p4 ) || (Membrane_Triangle_Pair_Nodes[i][0]==p4 &  Membrane_Triangle_Pair_Nodes[i][1]==p2 )    )
        {
            if(  Membrane_Triangle_Pair_Nodes[i][2]==p1  )
            {
                
                Dd= Membrane_Triangle_Pair_Nodes[i][3];//                 Membrane_Triangle_Pair_Nodes[i][2]=p3;
                lineD=i;
                counter=counter+1;
                findd=+1;
            }
            else  if(  Membrane_Triangle_Pair_Nodes[i][3]==p1 )
            {
                
                Dd= Membrane_Triangle_Pair_Nodes[i][2]; //                 Membrane_Triangle_Pair_Nodes[i][3]=p3;
                lineD=i;
                counter=counter+1;
                findd=+1;
            }
            //          else
            //   cout << " D ";
            
        }
        
    }// END OF FOR
    
    
    
    
    if(Aa!=Bb & Aa!=Cc & Aa!=Dd & Bb!=Cc & Bb!=Dd & Cc!=Dd  )
    {
        
        
        
        if(findpos==+1)
        {
            Membrane_Triangle_Pair_Nodes[linepos][0]=p3;
            Membrane_Triangle_Pair_Nodes[linepos][1]=p4;
            Membrane_Triangle_Pair_Nodes[linepos][2]=p1;
            Membrane_Triangle_Pair_Nodes[linepos][3]=p2;
            
            
            if(finda==+1)
            {
                Membrane_Triangle_Pair_Nodes[lineA][3]=Aa;
                Membrane_Triangle_Pair_Nodes[lineA][2]=p3;
            }
            if(findb==+1)
            {
                Membrane_Triangle_Pair_Nodes[lineB][3]=Bb;
                Membrane_Triangle_Pair_Nodes[lineB][2]=p4;
            }
            
            if(findc==+1)
            {
                Membrane_Triangle_Pair_Nodes[lineC][3]=Cc;
                Membrane_Triangle_Pair_Nodes[lineC][2]=p4;
            }
            
            
            if(finda==+1)
            {
                Membrane_Triangle_Pair_Nodes[lineD][3]=Dd;
                Membrane_Triangle_Pair_Nodes[lineD][2]=p3;
            }
            
        }
    }//endif
    
    
    
}

void updatetriangle(int p1, int p2,int p3,int p4 ,int Membrane_triangle_list[Membrane_num_of_Triangles][3],int &l1,int &l2 , double  Membrane_Node_Position [][3], int Membrane_Normal_direction[Membrane_num_of_Triangles][2])
{
    
    double a1a3[3],a1a2[3],nl1[3],nl2[3],ml1[3],ml2[3],Dl1,Dl2,v31[3],v32[3],v42[3],v41[3];
    
    
    int counting=0;
    for(  int i=0;i<Membrane_num_of_Triangles;i++  )  /// finding two triangles in trapezium
    {
        if(   (Membrane_triangle_list[i][0]==p1 & Membrane_triangle_list[i][1]==p2 & Membrane_triangle_list[i][2]==p3 ) ||
           (Membrane_triangle_list[i][0]==p1 & Membrane_triangle_list[i][1]==p3 & Membrane_triangle_list[i][2]==p2 ) ||
           (Membrane_triangle_list[i][0]==p2 & Membrane_triangle_list[i][1]==p1 & Membrane_triangle_list[i][2]==p3 ) ||
           (Membrane_triangle_list[i][0]==p2 & Membrane_triangle_list[i][1]==p3 & Membrane_triangle_list[i][2]==p1 ) ||
           (Membrane_triangle_list[i][0]==p3 & Membrane_triangle_list[i][1]==p1 & Membrane_triangle_list[i][2]==p2 ) ||
           (Membrane_triangle_list[i][0]==p3 & Membrane_triangle_list[i][1]==p2 & Membrane_triangle_list[i][2]==p1 )  )
        {
            l1=i;
            counting=counting+1;
        }
        if(   (Membrane_triangle_list[i][0]==p1 & Membrane_triangle_list[i][1]==p2 & Membrane_triangle_list[i][2]==p4 ) ||
           (Membrane_triangle_list[i][0]==p1 & Membrane_triangle_list[i][1]==p4 & Membrane_triangle_list[i][2]==p2 ) ||
           (Membrane_triangle_list[i][0]==p2 & Membrane_triangle_list[i][1]==p1 & Membrane_triangle_list[i][2]==p4 ) ||
           (Membrane_triangle_list[i][0]==p2 & Membrane_triangle_list[i][1]==p4 & Membrane_triangle_list[i][2]==p1 ) ||
           (Membrane_triangle_list[i][0]==p4 & Membrane_triangle_list[i][1]==p1 & Membrane_triangle_list[i][2]==p2 ) ||
           (Membrane_triangle_list[i][0]==p4 & Membrane_triangle_list[i][1]==p2 & Membrane_triangle_list[i][2]==p1 )  )
        {
            l2=i;
            counting=counting+1;
        }
    }
    
    //cout<< counting<<endl;
    //-------------------------------------------------------
    a1a2[0]=Membrane_Node_Position[ Membrane_triangle_list[l1][1]][0]-Membrane_Node_Position[ Membrane_triangle_list[l1][0]][0];
    a1a2[1]=Membrane_Node_Position[ Membrane_triangle_list[l1][1]][1]-Membrane_Node_Position[ Membrane_triangle_list[l1][0]][1];
    a1a2[2]=Membrane_Node_Position[ Membrane_triangle_list[l1][1]][2]-Membrane_Node_Position[ Membrane_triangle_list[l1][0]][2];
    a1a3[0]=Membrane_Node_Position[ Membrane_triangle_list[l1][2]][0]-Membrane_Node_Position[ Membrane_triangle_list[l1][0]][0];
    a1a3[1]=Membrane_Node_Position[ Membrane_triangle_list[l1][2]][1]-Membrane_Node_Position[ Membrane_triangle_list[l1][0]][1];
    a1a3[2]=Membrane_Node_Position[ Membrane_triangle_list[l1][2]][2]-Membrane_Node_Position[ Membrane_triangle_list[l1][0]][2];
    crossvector(nl1,a1a2,a1a3);
    nl1[0]=nl1[0];
    nl1[1]=nl1[1];
    nl1[2]=nl1[2];
    
    a1a2[0]=Membrane_Node_Position[ Membrane_triangle_list[l2][1]][0]-Membrane_Node_Position[ Membrane_triangle_list[l2][0]][0];
    a1a2[1]=Membrane_Node_Position[ Membrane_triangle_list[l2][1]][1]-Membrane_Node_Position[ Membrane_triangle_list[l2][0]][1];
    a1a2[2]=Membrane_Node_Position[ Membrane_triangle_list[l2][1]][2]-Membrane_Node_Position[ Membrane_triangle_list[l2][0]][2];
    a1a3[0]=Membrane_Node_Position[ Membrane_triangle_list[l2][2]][0]-Membrane_Node_Position[ Membrane_triangle_list[l2][0]][0];
    a1a3[1]=Membrane_Node_Position[ Membrane_triangle_list[l2][2]][1]-Membrane_Node_Position[ Membrane_triangle_list[l2][0]][1];
    a1a3[2]=Membrane_Node_Position[ Membrane_triangle_list[l2][2]][2]-Membrane_Node_Position[ Membrane_triangle_list[l2][0]][2];
    crossvector(nl2,a1a2,a1a3);
    nl2[0]=nl2[0];
    nl2[1]=nl2[1];
    nl2[2]=nl2[2];
    
    
    
    v31[0]=Membrane_Node_Position[ p1][0]-Membrane_Node_Position[ p3][0];
    v31[1]=Membrane_Node_Position[ p1][1]-Membrane_Node_Position[ p3][1];
    v31[2]=Membrane_Node_Position[ p1][2]-Membrane_Node_Position[ p3][2];
    v32[0]=Membrane_Node_Position[ p2][0]-Membrane_Node_Position[ p3][0];
    v32[1]=Membrane_Node_Position[ p2][1]-Membrane_Node_Position[ p3][1];
    v32[2]=Membrane_Node_Position[ p2][2]-Membrane_Node_Position[ p3][2];
    crossvector(ml1,v31,v32);
    
    v42[0]=Membrane_Node_Position[ p2][0]-Membrane_Node_Position[ p4][0];
    v42[1]=Membrane_Node_Position[ p2][1]-Membrane_Node_Position[ p4][1];
    v42[2]=Membrane_Node_Position[ p2][2]-Membrane_Node_Position[ p4][2];
    v41[0]=Membrane_Node_Position[ p1][0]-Membrane_Node_Position[ p4][0];
    v41[1]=Membrane_Node_Position[ p1][1]-Membrane_Node_Position[ p4][1];
    v41[2]=Membrane_Node_Position[ p1][2]-Membrane_Node_Position[ p4][2];
    crossvector(ml2,v42,v41);
    
    if( ml1[0]!=0.0  & nl1[0]!=0.0)
    {
        Dl1=nl1[0]/ml1[0];
    }
    else if( ml1[1]!=0.0  & nl1[1]!=0.0 )
    {
        Dl1=nl1[1]/ml1[1];
    }
    else if( ml1[2]!=0.0  & nl1[2]!=0.0)
    {
        Dl1=nl1[2]/ml1[2];
    }
    
    
    
    if( ml2[0]!=0.0    & nl2[0]!=0.0 )
    {
        Dl2=nl2[0]/ml2[0];
    }
    else if( ml2[1]!=0.0  & nl2[1]!=0.0 )
    {
        Dl2=nl2[1]/ml2[1];
    }
    else if( ml1[2]!=0.0   & nl2[2]!=0.0 )
    {
        Dl2=nl2[2]/ml2[2];
    }
    
    
    
    //-------------------------------------------------------
    if(Dl1>0)
    {
        Membrane_triangle_list[l1][0]=p3;
        Membrane_triangle_list[l1][1]=p1;
        Membrane_triangle_list[l1][2]=p4;
    }
    else if(Dl1<0)
    {
        Membrane_triangle_list[l1][0]=p3;
        Membrane_triangle_list[l1][1]=p4;
        Membrane_triangle_list[l1][2]=p1;
    }
    
    
    if(Dl2>0)
    {
        Membrane_triangle_list[l2][0]=p4;
        Membrane_triangle_list[l2][1]=p2;
        Membrane_triangle_list[l2][2]=p3;
    }
    else if(Dl2<0)
    {
        Membrane_triangle_list[l2][0]=p4;
        Membrane_triangle_list[l2][1]=p3;
        Membrane_triangle_list[l2][2]=p2;
    }
    
    
    
}
double potentialenergy(double Membrane_Node_Position[][3],int Membrane_Node_Pair_list[][2],int Membrane_Triangle_Pair_Nodes[][4],int Membrane_triangle_list[Membrane_num_of_Triangles][3], int Outer_Membrane_num_of_triangles, int Nucleus_Membrane_num_of_triangles, int Membrane_num_of_Triangle_Pairs, int Membrane_num_of_Node_Pairs)//  potential energy
{
    double le0,le1,lmax,lmin;
    double deltax=0, deltay=0, deltaz=0,temp_Node_distance,temp_force;
    double Utemp = 0.0;
    double xpos[3][5],N1[3],N2[3],p3p1[3],p3p2[3],p4p2[3],p4p1[3],Ep2p1[3];// for exmple p3p1 is p3-p1 and ....
    double U=0.0;
    
    
    /// calculate network force:
    int temp_Node_A,temp_Node_B;
    le0=1.15000*Node_radius;
    lmax=1.33000*Node_radius;
    le1=0.85000*Node_radius;
    lmin=0.67000*Node_radius;
    U=0.0;
    
    
    for (int k=0 ; k< Membrane_num_of_Node_Pairs ; k++)
    {
        temp_Node_B=Membrane_Node_Pair_list[k][0];
        temp_Node_A=Membrane_Node_Pair_list[k][1];
        deltax=Membrane_Node_Position[temp_Node_A][0]-Membrane_Node_Position[temp_Node_B][0];
        deltay=Membrane_Node_Position[temp_Node_A][1]-Membrane_Node_Position[temp_Node_B][1];
        deltaz=Membrane_Node_Position[temp_Node_A][2]-Membrane_Node_Position[temp_Node_B][2];
        
        if (Periodic_condtion_status == 1.0) {
            deltax=periodiccondition(deltax);
            deltay=periodiccondition(deltay);
        }
        
        temp_Node_distance=sqrt(deltax*deltax+deltay*deltay+deltaz*deltaz);
        temp_force=0.0;
        double temp_exp_le0=exp(1.0/(le0-temp_Node_distance));
        double temp_exp_le1=exp(1.0/(temp_Node_distance-le1));
        if(temp_Node_distance >le1  & temp_Node_distance < le0 )  //free zone
        {
            Utemp=0 ; // free zone
        }
        
        if(temp_Node_distance > le0  & temp_Node_distance <lmax ) { //bondforce
            
            temp_force = (Membrane_spring_coefficient*temp_exp_le0/(lmax-temp_Node_distance))*( 1.0/(lmax-temp_Node_distance) +  1.0/((le0-temp_Node_distance)*(le0-temp_Node_distance)));
            Utemp= Membrane_spring_coefficient*temp_exp_le0/(lmax-temp_Node_distance);
        }
        if(temp_Node_distance < le1   &  temp_Node_distance > lmin  ) { // repulsive force
            
            temp_force= -(Membrane_spring_coefficient*temp_exp_le1/(temp_Node_distance-lmin))*( 1.0/(temp_Node_distance-lmin) + 1.0/((temp_Node_distance-le1)*(temp_Node_distance-le1)));                 // force on i th from j
            Utemp= Membrane_spring_coefficient*temp_exp_le1/(temp_Node_distance-lmin);
        }
        if(temp_force>965.31  || temp_Node_distance>lmax ) {
            temp_force = 965.31+Membrane_spring_force_cutt_off* ( temp_Node_distance - 1.3280*Node_radius );
            Utemp=   1.81599  + 965.31 * ( temp_Node_distance - 1.3280*Node_radius )+0.5*Membrane_spring_force_cutt_off * ( temp_Node_distance - 1.3280*Node_radius ) * ( temp_Node_distance - 1.3280*Node_radius );
        }
        if(temp_force<-1000.05   ||  temp_Node_distance<lmin ) {
            temp_force =-1000.05-Membrane_spring_force_cutt_off* ( 0.671965*Node_radius - temp_Node_distance );
            Utemp=   1.85038 + 1005.05 * ( 0.671965*Node_radius - temp_Node_distance )+0.5*Membrane_spring_force_cutt_off*( 0.671965*Node_radius - temp_Node_distance )*( 0.671965*Node_radius - temp_Node_distance );
        }
        U+=Utemp;
        
        // implimentation of forces:
        
    }
    /// calculating  triangle-triangle  force (surface force)
    
    for(int i=0 ;i<Membrane_num_of_Triangle_Pairs;i++)
    {
        //PERIODIC
        
        if (Periodic_condtion_status == 1.0 ) {
            for(int i2=0;i2<3;i2++)   // storing coordiates of poses in new array in order to imposing periodic boundary
            {
                
                xpos[i2][1]=Membrane_Node_Position[Membrane_Triangle_Pair_Nodes[i][0]][i2];
                xpos[i2][2]=Membrane_Node_Position[Membrane_Triangle_Pair_Nodes[i][1]][i2];
                xpos[i2][3]=Membrane_Node_Position[Membrane_Triangle_Pair_Nodes[i][2]][i2];
                xpos[i2][4]=Membrane_Node_Position[Membrane_Triangle_Pair_Nodes[i][3]][i2];
            }
            
            for(int i2=0;i2<2;i2++)   // puting pos1 in (0,0,?)  --- only x and y  ,z remains untouch
            {
                
                xpos[i2][2]=xpos[i2][2]-xpos[i2][1];
                xpos[i2][3]=xpos[i2][3]-xpos[i2][1];
                xpos[i2][4]=xpos[i2][4]-xpos[i2][1];
                xpos[i2][1]=0.0;
            }
            if(  xpos[0][2]> 6.0*Node_radius  )
                xpos[0][2]= xpos[0][2]-(Lbox+1.0);
            
            if(  xpos[1][2]> 6.0*Node_radius  )
                xpos[1][2]= xpos[1][2]-(Lbox+1.0);
            
            if(  xpos[0][2]< -6.0*Node_radius  )
                xpos[0][2]= xpos[0][2]+(Lbox+1.0);
            
            if(  xpos[1][2]< -6.0*Node_radius  )
                xpos[1][2]= xpos[1][2]+(Lbox+1.0);
            
            
            
            if(  xpos[0][3]> 6.0*Node_radius  )
                xpos[0][3]= xpos[0][3]-(Lbox+1.0);
            
            if(  xpos[1][3]> 6.0*Node_radius  )
                xpos[1][3]= xpos[1][3]-(Lbox+1.0);
            
            if(  xpos[0][3]< -6.0*Node_radius  )
                xpos[0][3]= +xpos[0][3]+(Lbox+1.0);
            
            if(  xpos[1][3]<- 6.0*Node_radius  )
                xpos[1][3]= +xpos[1][3]+(Lbox+1.0);
            
            
            if(  xpos[0][4]> 6.0*Node_radius  )
                xpos[0][4]= xpos[0][4]-(Lbox+1.0);
            
            if(  xpos[1][4]> 6.0*Node_radius  )
                xpos[1][4]= xpos[1][4]-(Lbox+1.0);
            
            if(  xpos[0][4]< -6.0*Node_radius  )
                xpos[0][4]= +xpos[0][4]+(Lbox+1.0);
            
            if(  xpos[1][4]<- 6.0*Node_radius  )
                xpos[1][4]= +xpos[1][4]+(Lbox+1.0);
        }
        //PERIODIC
        
        for (int l=0; l<3; l++) {
            p3p1[l]=Membrane_Node_Position[Membrane_Triangle_Pair_Nodes[i][2]][l]-Membrane_Node_Position[Membrane_Triangle_Pair_Nodes[i][0]][l];
            p3p2[l]=Membrane_Node_Position[Membrane_Triangle_Pair_Nodes[i][2]][l]-Membrane_Node_Position[Membrane_Triangle_Pair_Nodes[i][1]][l];
            p4p2[l]=Membrane_Node_Position[Membrane_Triangle_Pair_Nodes[i][3]][l]-Membrane_Node_Position[Membrane_Triangle_Pair_Nodes[i][1]][l];
            p4p1[l]=Membrane_Node_Position[Membrane_Triangle_Pair_Nodes[i][3]][l]-Membrane_Node_Position[Membrane_Triangle_Pair_Nodes[i][0]][l];
            Ep2p1[l]=Membrane_Node_Position[Membrane_Triangle_Pair_Nodes[i][1]][l]-Membrane_Node_Position[Membrane_Triangle_Pair_Nodes[i][0]][l];
        }
        
        
        crossvector(N1,p3p1,p3p2);
        crossvector(N2,p4p2,p4p1);
        
        U+=Membrane_bending_coefficient*(1.0 - ( innerproduct(N1,N2)/(vectorlength(N1)*vectorlength(N2)) ) );
        
    }  // for
    
    //I had a very breif (under 2 mins) chat with Prof Ejtehadi and we agreed that we don't need to calculate this during the MC step as the constraint will be applied to the system and timewise it should be and is 'long' compared to the dynamics of the system.
    //**********area
    double Smem=Membrane_surface_area_calculator(Membrane_Node_Position, Membrane_triangle_list, Outer_Membrane_num_of_triangles);
    double Snuc=surfaceareaNucleus(Membrane_Node_Position,Membrane_triangle_list, Outer_Membrane_num_of_triangles);
    double  s0membrane=0.41*Node_radius*Node_radius*Outer_Membrane_num_of_triangles; // initialize the initial surface area
    double  s0nucleus=0.41*Node_radius*Node_radius*Nucleus_Membrane_num_of_triangles; // initialize the initial surface area
    double Uarea=0;
    Uarea=0.5*(K_surfaceConstant_local)*(s0membrane-Smem)*(s0membrane-Smem) + 0.5*(K_surfaceConstant_local)*(s0nucleus-Snuc)*(s0nucleus-Snuc) ;
    //**********area
    
    U=U+Uarea;
    
    return U;
}

double periodiccondition(double dx )
{
    if( dx>Lbox/2.0  )
    {
        return dx-(Lbox+1.0); /// so the new dx is <0
    }
    else if(dx<-Lbox/2.0  )
    {
        return (Lbox+1.0)+dx;  /// so the new dx is >0
    }
    else
    {
        return dx;
    }
}

double Membrane_surface_area_calculator(double  Membrane_Node_Position [][3],  int Membrane_triangle_list[Membrane_num_of_Triangles][3], int Outer_Membrane_num_of_triangles)
{
    double area=0.0;
    int temp_triangle_Node_A,temp_triangle_Node_B,temp_triangle_Node_C;
    double AB[3],AC[3],ABxAC[3];
    
    for (int i=0 ; i<Outer_Membrane_num_of_triangles ; i++  ) /// calculating the sum of area of all triangles :
    {
        temp_triangle_Node_A=Membrane_triangle_list[i][0];
        temp_triangle_Node_B=Membrane_triangle_list[i][1];
        temp_triangle_Node_C=Membrane_triangle_list[i][2];
        for (int l=0; l<3; l++) {
            AB[l]=Membrane_Node_Position[temp_triangle_Node_A][l]-Membrane_Node_Position[temp_triangle_Node_B][l];
            
            AC[l]=Membrane_Node_Position[temp_triangle_Node_A][l]-Membrane_Node_Position[temp_triangle_Node_C][l];
        }
        
        
        if (Periodic_condtion_status == 1.0) {
            AB[0]=periodiccondition( AB[0]);
            AB[1]=periodiccondition( AB[1]);
            AB[2]=periodiccondition( AB[2]);
            
            AC[0]=periodiccondition( AC[0]);
            AC[1]=periodiccondition( AC[1]);
            AC[2]=periodiccondition( AC[2]);
        }
        
        crossvector( ABxAC,AB,AC  );
        area += vectorlength(ABxAC)/2.0;
    }
    
    return area;
}



void ConstantSurfaceForceLocalTriangles(double Membrane_Node_Position[][3],double Membrane_Node_Force[][3],int Membrane_triangle_list[Membrane_num_of_Triangles][3], int Outer_Membrane_num_of_triangles)
{
    
    int temp_node_A,temp_node_B,temp_node_C;
    double AB[3],AC[3],temp_AB_length_squared,temp_AC_length_squared,f0, temp_ABxAC[3];
    
    // cout <<surfacearea(x,tri )<<"   s0="<< s0<<endl;
    double s0_i=0.41*Node_radius*Node_radius; //1.732=3^0.5
    double s_i;
    
    for(  int i=0;i<Membrane_num_of_Triangles;i++)
    {
        
        //*********
        temp_node_A=Membrane_triangle_list[i][0];
        temp_node_B=Membrane_triangle_list[i][1];
        temp_node_C=Membrane_triangle_list[i][2];
        
        
        
        AB[0]=Membrane_Node_Position[temp_node_B][0]-Membrane_Node_Position[temp_node_A][0];
        AB[1]=Membrane_Node_Position[temp_node_B][1]-Membrane_Node_Position[temp_node_A][1];
        AB[2]=Membrane_Node_Position[temp_node_B][2]-Membrane_Node_Position[temp_node_A][2];
        
        AC[0]=Membrane_Node_Position[temp_node_C][0]-Membrane_Node_Position[temp_node_A][0];
        AC[1]=Membrane_Node_Position[temp_node_C][1]-Membrane_Node_Position[temp_node_A][1];
        AC[2]=Membrane_Node_Position[temp_node_C][2]-Membrane_Node_Position[temp_node_A][2];
        
        
        crossvector(temp_ABxAC,AB,AC);
        s_i = vectorlength(temp_ABxAC)/2.0;
        // cout << s_i <<endl;
        
        temp_AB_length_squared=vectorlength(AB)*vectorlength(AB);
        temp_AC_length_squared=vectorlength(AC)*vectorlength(AC);
        
        f0= K_surfaceConstant_local*(s_i  -  s0_i )/2.0*vectorlength(temp_ABxAC);
        
        Membrane_Node_Force[temp_node_A][0] +=  f0 *2.0* ( -temp_AC_length_squared *  AB[0]  -temp_AB_length_squared *  AC[0] + innerproduct(AB,AC)*  ( AB[0]+AC[0] )   ) ;
        Membrane_Node_Force[temp_node_A][1] +=  f0 *2.0* ( -temp_AC_length_squared *  AB[1]  -temp_AB_length_squared *  AC[1] + innerproduct(AB,AC)*  ( AB[1]+AC[1] )   ) ;
        Membrane_Node_Force[temp_node_A][2] +=  f0 *2.0* ( -temp_AC_length_squared *  AB[2]  -temp_AB_length_squared *  AC[2] + innerproduct(AB,AC)*  ( AB[2]+AC[2] )   ) ;
        
        //*********
        temp_node_A=Membrane_triangle_list[i][1];
        temp_node_B=Membrane_triangle_list[i][0];
        temp_node_C=Membrane_triangle_list[i][2];
        
        
        
        AB[0]=Membrane_Node_Position[temp_node_B][0]-Membrane_Node_Position[temp_node_A][0];
        AB[1]=Membrane_Node_Position[temp_node_B][1]-Membrane_Node_Position[temp_node_A][1];
        AB[2]=Membrane_Node_Position[temp_node_B][2]-Membrane_Node_Position[temp_node_A][2];
        
        AC[0]=Membrane_Node_Position[temp_node_C][0]-Membrane_Node_Position[temp_node_A][0];
        AC[1]=Membrane_Node_Position[temp_node_C][1]-Membrane_Node_Position[temp_node_A][1];
        AC[2]=Membrane_Node_Position[temp_node_C][2]-Membrane_Node_Position[temp_node_A][2];
        
        temp_AB_length_squared=vectorlength(AB)*vectorlength(AB);
        temp_AC_length_squared=vectorlength(AC)*vectorlength(AC);
        
        Membrane_Node_Force[temp_node_A][0] +=  f0 *2.0* ( -temp_AC_length_squared *  AB[0]  -temp_AB_length_squared *  AC[0] + innerproduct(AB,AC)*  ( AB[0]+AC[0] )   ) ;
        Membrane_Node_Force[temp_node_A][1] +=  f0 *2.0* ( -temp_AC_length_squared *  AB[1]  -temp_AB_length_squared *  AC[1] + innerproduct(AB,AC)*  ( AB[1]+AC[1] )   ) ;
        Membrane_Node_Force[temp_node_A][2] +=  f0 *2.0* ( -temp_AC_length_squared *  AB[2]  -temp_AB_length_squared *  AC[2] + innerproduct(AB,AC)*  ( AB[2]+AC[2] )   ) ;
        //*********
        
        
        
        //*********
        temp_node_A=Membrane_triangle_list[i][2];
        temp_node_B=Membrane_triangle_list[i][1];
        temp_node_C=Membrane_triangle_list[i][0];
        
        AB[0]=Membrane_Node_Position[temp_node_B][0]-Membrane_Node_Position[temp_node_A][0];
        AB[1]=Membrane_Node_Position[temp_node_B][1]-Membrane_Node_Position[temp_node_A][1];
        AB[2]=Membrane_Node_Position[temp_node_B][2]-Membrane_Node_Position[temp_node_A][2];
        
        AC[0]=Membrane_Node_Position[temp_node_C][0]-Membrane_Node_Position[temp_node_A][0];
        AC[1]=Membrane_Node_Position[temp_node_C][1]-Membrane_Node_Position[temp_node_A][1];
        AC[2]=Membrane_Node_Position[temp_node_C][2]-Membrane_Node_Position[temp_node_A][2];
        
        temp_AB_length_squared=vectorlength(AB)*vectorlength(AB);
        temp_AC_length_squared=vectorlength(AC)*vectorlength(AC);
        
        
        Membrane_Node_Force[temp_node_A][0] +=  f0 *2.0* ( -temp_AC_length_squared *  AB[0]  -temp_AB_length_squared *  AC[0] + innerproduct(AB,AC)*  ( AB[0]+AC[0] )   ) ;
        Membrane_Node_Force[temp_node_A][1] +=  f0 *2.0* ( -temp_AC_length_squared *  AB[1]  -temp_AB_length_squared *  AC[1] + innerproduct(AB,AC)*  ( AB[1]+AC[1] )   ) ;
        Membrane_Node_Force[temp_node_A][2] +=  f0 *2.0* ( -temp_AC_length_squared *  AB[2]  -temp_AB_length_squared *  AC[2] + innerproduct(AB,AC)*  ( AB[2]+AC[2] )   ) ;
        //*********
    }
    
}

double surfaceareaNucleus(double  Membrane_Node_Position [][3],  int Membrane_triangle_list[Membrane_num_of_Triangles][3], int Outer_Membrane_num_of_triangles)
{
    double area;
    int p1,p2,p3;
    double vec1[3],vec2[3],vec3[3];
    
    area=0.0;
    for (int i=Outer_Membrane_num_of_triangles ; i<Membrane_num_of_Triangles ; i++  ) /// calculating the sum of area of all triangles :
    {
        p1=Membrane_triangle_list[i][0];
        p2=Membrane_triangle_list[i][1];
        p3=Membrane_triangle_list[i][2];
        
        vec1[0]=Membrane_Node_Position[p1][0]-Membrane_Node_Position[p2][0];
        vec1[1]=Membrane_Node_Position[p1][1]-Membrane_Node_Position[p2][1];
        vec1[2]=Membrane_Node_Position[p1][2]-Membrane_Node_Position[p2][2];
        
        vec2[0]=Membrane_Node_Position[p1][0]-Membrane_Node_Position[p3][0];
        vec2[1]=Membrane_Node_Position[p1][1]-Membrane_Node_Position[p3][1];
        vec2[2]=Membrane_Node_Position[p1][2]-Membrane_Node_Position[p3][2];
        
        if (Periodic_condtion_status ==1.0) {
            vec1[0]=periodiccondition( vec1[0]);
            vec1[1]=periodiccondition( vec1[1]);
            vec1[2]=periodiccondition( vec1[2]);
            
            vec2[0]=periodiccondition( vec2[0]);
            vec2[1]=periodiccondition( vec2[1]);
            vec2[2]=periodiccondition( vec2[2]);
        }
        
        crossvector( vec3,vec1,vec2  );
        area = area+ vectorlength(vec3)/2.0;
    }
    
    return area;
}




double VolumeNucleus(double Membrane_Node_Position[][3],int Membrane_triangle_list[Membrane_num_of_Triangles][3],int Membrane_Normal_direction[Membrane_num_of_Triangles][2], int Outer_Membrane_num_of_triangles)
{
    int p1,p2,p3;
    double Volume=0.0;
    double vec1[3],vec2[3],vec3[3];
    double x1,x2,x3,y1,y2,y3,z1,z2,z3;
    
    
    
    for (int i =Outer_Membrane_num_of_triangles; i < Membrane_num_of_Triangles  ; i++)  // for each triangle of nucleus /// calculating the sum of Vulume of surface :
        
    {
        p1=Membrane_triangle_list[i][0];
        p2=Membrane_triangle_list[i][1];
        p3=Membrane_triangle_list[i][2];
        
        x1=Membrane_Node_Position[p1][0];
        x2=Membrane_Node_Position[p2][0];
        x3=Membrane_Node_Position[p3][0];
        
        y1=Membrane_Node_Position[p1][1];
        y2=Membrane_Node_Position[p2][1];
        y3=Membrane_Node_Position[p3][1];
        
        z1=Membrane_Node_Position[p1][2];
        z2=Membrane_Node_Position[p2][2];
        z3=Membrane_Node_Position[p3][2];
        
        vec1[0]=Membrane_Node_Position[p2][0]-Membrane_Node_Position[p1][0];
        vec1[1]=Membrane_Node_Position[p2][1]-Membrane_Node_Position[p1][1];
        vec1[2]=Membrane_Node_Position[p2][2]-Membrane_Node_Position[p1][2];
        
        vec2[0]=Membrane_Node_Position[p3][0]-Membrane_Node_Position[p1][0];
        vec2[1]=Membrane_Node_Position[p3][1]-Membrane_Node_Position[p1][1];
        vec2[2]=Membrane_Node_Position[p3][2]-Membrane_Node_Position[p1][2];
        
        crossvector( vec3,vec1,vec2);
        
        Volume=Volume+  (double)  Membrane_Normal_direction [i][1] * ( (x1+x2+x3)*vec3[0]+(y1+y2+y3)*vec3[1]+(z1+z2+z3)*vec3[2]  );
    }
    return Volume/18.0;
}




double VolumeMemerane(double Membrane_Node_Position[][3],int tri[3][Membrane_num_of_Triangles],int Membrane_Normal_direction[Membrane_num_of_Triangles][2], int Outer_Membrane_num_of_triangles)
{
    int p1,p2,p3;
    double Volume=0.0;
    double vec1[3],vec2[3],vec3[3];
    double x1,x2,x3,y1,y2,y3,z1,z2,z3;
    
    for (int i =0; i < Outer_Membrane_num_of_triangles  ; i++)  // for each triangle of nucleus /// calculating the sum of Vulume of surface :
        
    {
        p1=tri[0][i];
        p2=tri[1][i];
        p3=tri[2][i];
        
        x1=Membrane_Node_Position[p1][0];
        x2=Membrane_Node_Position[p2][0];
        x3=Membrane_Node_Position[p3][0];
        
        y1=Membrane_Node_Position[p1][1];
        y2=Membrane_Node_Position[p2][1];
        y3=Membrane_Node_Position[p3][1];
        
        z1=Membrane_Node_Position[p1][2];
        z2=Membrane_Node_Position[p2][2];
        z3=Membrane_Node_Position[p3][2];
        
        vec1[0]=Membrane_Node_Position[p2][0]-Membrane_Node_Position[p1][0];
        vec1[1]=Membrane_Node_Position[p2][1]-Membrane_Node_Position[p1][1];
        vec1[2]=Membrane_Node_Position[p2][2]-Membrane_Node_Position[p1][2];
        
        vec2[0]=Membrane_Node_Position[p3][0]-Membrane_Node_Position[p1][0];
        vec2[1]=Membrane_Node_Position[p3][1]-Membrane_Node_Position[p1][1];
        vec2[2]=Membrane_Node_Position[p3][2]-Membrane_Node_Position[p1][2];
        
        crossvector( vec3,vec1,vec2);
        
        Volume += (double) Membrane_Normal_direction [i][1] * ( (x1+x2+x3)*vec3[0]+(y1+y2+y3)*vec3[1]+(z1+z2+z3)*vec3[2]  );
    }
    return Volume/18.0;
}

//______________________Actin functions








void  Actin_Force_calculator( double  Actin_Node_Position [][3],double  Actin_Node_VelocityRungKuta [][3],double  Actin_Node_Force [][3],double Actin_Node_Pair_List[][3], int Actin_num_of_Bonds, double &Total_Potential_Energy)
{
    
    double deltax,deltay,deltaz ,temp_distance, initial_distance;// defined below in "for loop" in detail
    int node1,node2;
    double  temp_force[3];
    
    /// Spring
    for(int i=0 ;i<Actin_num_of_Bonds ; i++ )  // all beads interaction whit the next one
    {
        node1=(int) Actin_Node_Pair_List[i][0];
        node2=(int) Actin_Node_Pair_List[i][1];
        
        initial_distance=Actin_Node_Pair_List[i][2];
        
        deltax=0.0;
        deltay=0.0;
        deltaz=0.0;
        temp_distance=0.0;
        
        deltax=Actin_Node_Position[node2][0]-Actin_Node_Position[node1][0];// delta x  betwwn i and i+1 th beads
        deltay=Actin_Node_Position[node2][1]-Actin_Node_Position[node1][1];// delta y  betwwn i and i+1 th beads
        deltaz=Actin_Node_Position[node2][2]-Actin_Node_Position[node1][2];// delta z  betwwn i and i+1 th beads
        temp_distance=sqrt(deltax*deltax+deltay*deltay+deltaz*deltaz); // distance btween i th and i+1 th  bead
        
        
        
        if(CytoskeletonNetworkType==0)  {// if PCN is off use simple hookian network
            
            temp_force[0]=-Actin_spring_coefficient*(temp_distance-initial_distance) *  deltax/temp_distance  ;
            temp_force[1]=-Actin_spring_coefficient*(temp_distance-initial_distance) *  deltay/temp_distance  ;
            temp_force[2]=-Actin_spring_coefficient*(temp_distance-initial_distance) *  deltaz/temp_distance  ;
            if (energy_calculation_flag==1.0) {
                Total_Potential_Energy += 0.5*Actin_spring_coefficient*(temp_distance-initial_distance)*(temp_distance-initial_distance);
                
            }
            
            Actin_Node_Force[node1][0] += temp_force[0]  ; // force of springs and dashes
            Actin_Node_Force[node1][1] += temp_force[1]  ;// force of springs and dashes
            Actin_Node_Force[node1][2] += temp_force[2]  ;// force of springs and dashes
            
            Actin_Node_Force[node2][0] += - temp_force[0]  ; // force of springs and dashes
            Actin_Node_Force[node2][1] += - temp_force[1]  ;// force of springs and dashes
            Actin_Node_Force[node2][2] += - temp_force[2]  ;// force of springs and dashes
            
            Actin_Node_Pair_List[i][2] +=  MD_Time_Step*   (  (temp_distance-initial_distance)/temp_distance  )   /Actin_kelvin_damping_coefficient; //We have used the Kelvin model for the Visco elasticity. Kelvin model: A spring and dashpot (in series) are in parralel with a spring.
            
        }  else if(CytoskeletonNetworkType==1  && (temp_distance-initial_distance) > 0.0)  {// passive cable network  on?
            temp_force[0]=-Actin_Passive_Cable_Network_Coefficient*(temp_distance-initial_distance)*  deltax/temp_distance  ;
            temp_force[1]=-Actin_Passive_Cable_Network_Coefficient*(temp_distance-initial_distance)*  deltay/temp_distance  ;
            temp_force[2]=-Actin_Passive_Cable_Network_Coefficient*(temp_distance-initial_distance)*  deltaz/temp_distance  ;
            if (energy_calculation_flag==1.0) {
                Total_Potential_Energy += 0.5*Actin_Passive_Cable_Network_Coefficient*(temp_distance-initial_distance)*(temp_distance-initial_distance);
            }
            
            
            Actin_Node_Force[node1][0] += temp_force[0]  ;// force of springs and dashes
            Actin_Node_Force[node1][1] += temp_force[1]  ;// force of springs and dashes
            Actin_Node_Force[node1][2] += temp_force[2]  ;// force of springs and dashes
            
            Actin_Node_Force[node2][0] += - temp_force[0]  ;// force of springs and dashes
            Actin_Node_Force[node2][1] += - temp_force[1]  ;// force of springs and dashes
            Actin_Node_Force[node2][2] += - temp_force[2]  ;// force of springs and dashes
            
        } else if(CytoskeletonNetworkType==2 )  {//active cable network
            
            if(temp_distance > initial_distance)  { //more detail in the paper: Contractile network models for adherent cells
                
                temp_force[0]=-(KActinACN_EA *(temp_distance-initial_distance) +ACN_TL0)   *  deltax/temp_distance  ;
                temp_force[1]=-(KActinACN_EA *(temp_distance-initial_distance) +ACN_TL0)   *  deltay/temp_distance  ;
                temp_force[2]=-(KActinACN_EA *(temp_distance-initial_distance) +ACN_TL0)   *  deltaz/temp_distance  ;
                if (energy_calculation_flag==1.0) {
                    Total_Potential_Energy += 0.5*KActinACN_EA*(temp_distance-initial_distance)*(temp_distance-initial_distance) + ACN_TL0*(temp_distance-initial_distance);
                }
                
            }
            else if(   (temp_distance < initial_distance)  &&  temp_distance>ACN_LC )   //more detail in the paper: Contractile network models for adherent cells
            {
                temp_force[0]=-(ACN_TL0)   *  deltax/temp_distance  ;
                temp_force[1]=-(ACN_TL0)   *  deltay/temp_distance  ;
                temp_force[2]=-(ACN_TL0)   *  deltaz/temp_distance  ;
                
                if (energy_calculation_flag==1.0) {
                    Total_Potential_Energy += ACN_TL0*(temp_distance-initial_distance);
                }
                
            }
            
            else if(temp_distance < ACN_LC)   //more detail in the paper: Contractile network models for adherent cells
                
            {
                temp_force[0]=-( ACN_TL0 *(temp_distance)/ACN_LC )   *  deltax/temp_distance  ;
                temp_force[1]=-( ACN_TL0 *(temp_distance)/ACN_LC )   *  deltay/temp_distance  ;
                temp_force[2]=-( ACN_TL0 *(temp_distance)/ACN_LC )   *  deltaz/temp_distance  ;
                if (energy_calculation_flag==1.0) {
                    Total_Potential_Energy += 0.5*ACN_TL0*(temp_distance)*(temp_distance)/ACN_LC;
                }
                
            }
            
            
            Actin_Node_Force[node1][0] +=  temp_force[0]  ;// force of springs and dashes
            Actin_Node_Force[node1][1] +=  temp_force[1]  ;// force of springs and dashes
            Actin_Node_Force[node1][2] +=  temp_force[2]  ;// force of springs and dashes
            
            Actin_Node_Force[node2][0] += - temp_force[0]  ;// force of springs and dashes
            Actin_Node_Force[node2][1] += - temp_force[1]  ;// force of springs and dashes
            Actin_Node_Force[node2][2] += - temp_force[2]  ;// force of springs and dashes
            
        }
        
        ///Damping force
        temp_force[0]= Actin_damping_Coefficient*(Actin_Node_VelocityRungKuta[node1][0]-Actin_Node_VelocityRungKuta[node2][0]);
        temp_force[1]= Actin_damping_Coefficient*(Actin_Node_VelocityRungKuta[node1][1]-Actin_Node_VelocityRungKuta[node2][1]);
        temp_force[2]= Actin_damping_Coefficient*(Actin_Node_VelocityRungKuta[node1][2]-Actin_Node_VelocityRungKuta[node2][2]);
        if (energy_calculation_flag==1.0) {
            Total_Potential_Energy -= 0.5*Actin_damping_Coefficient*((Actin_Node_VelocityRungKuta[node1][2]-Actin_Node_VelocityRungKuta[node2][2])*(Actin_Node_VelocityRungKuta[node1][2]-Actin_Node_VelocityRungKuta[node2][2])+(Actin_Node_VelocityRungKuta[node1][1]-Actin_Node_VelocityRungKuta[node2][1])*(Actin_Node_VelocityRungKuta[node1][1]-Actin_Node_VelocityRungKuta[node2][1])+(Actin_Node_VelocityRungKuta[node1][0]-Actin_Node_VelocityRungKuta[node2][0])*(Actin_Node_VelocityRungKuta[node1][0]-Actin_Node_VelocityRungKuta[node2][0]));
        }
        
        
        Actin_Node_Force[node1][0] += temp_force[0]  ; // force of springs and dashes
        Actin_Node_Force[node1][1] += temp_force[1]  ;// force of springs and dashes
        Actin_Node_Force[node1][2] += temp_force[2]  ;// force of springs and dashes
        
        Actin_Node_Force[node2][0] += - temp_force[0]  ; // force of springs and dashes
        Actin_Node_Force[node2][1] += - temp_force[1]  ;// force of springs and dashes
        Actin_Node_Force[node2][2] += - temp_force[2]  ;// force of springs and dashes
        
    }
    
}



void Actin_Membrane_Barrier( double  Actin_Node_Position [][3], double  Actin_Node_Velocity [][3], double Membrane_Node_Position[][3], double Membrane_Node_Velocity[][3],  int Membrane_triangle_list[Membrane_num_of_Triangles][3], int Membrane_Normal_direction[Membrane_num_of_Triangles][2])
{
    double membrane_triangle_COM_position[3]; // coordinates to the centre of mass of the triangles
    double membrane_triangle_COM_velocity[3]; // velocity of the center of mass of the triangles
    double actin_membrane_distance_amplitude;// distance between solvent particle and com of triangle
    double ABxAC[3],AB[3],AC[3]; // normal vector of membrane
    double actin_triangle_distance_vector[3];
    double perpendicular_distance;
    double relevant_velocity[3];
    
    
    for (int i =0; i < Membrane_num_of_Triangles  ; i++)  // for each triangle:
    {
        membrane_triangle_COM_position[0]=(Membrane_Node_Position[ Membrane_triangle_list[i][0]][0] + Membrane_Node_Position[Membrane_triangle_list[i][1]][0] + Membrane_Node_Position[ Membrane_triangle_list[i][2]][0])/3.0;
        membrane_triangle_COM_position[1]=(Membrane_Node_Position[ Membrane_triangle_list[i][0]][1] + Membrane_Node_Position[Membrane_triangle_list[i][1]][1] + Membrane_Node_Position[ Membrane_triangle_list[i][2]][1])/3.0;
        membrane_triangle_COM_position[2]=(Membrane_Node_Position[ Membrane_triangle_list[i][0]][2] + Membrane_Node_Position[Membrane_triangle_list[i][1]][2] + Membrane_Node_Position[ Membrane_triangle_list[i][2]][2])/3.0;
        
        membrane_triangle_COM_velocity[0]=(Membrane_Node_Velocity[ Membrane_triangle_list[i][0]][0] + Membrane_Node_Velocity[Membrane_triangle_list[i][1]][0] + Membrane_Node_Velocity[ Membrane_triangle_list[i][2]][0])/3.0;
        membrane_triangle_COM_velocity[1]=(Membrane_Node_Velocity[ Membrane_triangle_list[i][0]][1] + Membrane_Node_Velocity[Membrane_triangle_list[i][1]][1] + Membrane_Node_Velocity[ Membrane_triangle_list[i][2]][1])/3.0;
        membrane_triangle_COM_velocity[2]=(Membrane_Node_Velocity[ Membrane_triangle_list[i][0]][2] + Membrane_Node_Velocity[Membrane_triangle_list[i][1]][2] + Membrane_Node_Velocity[ Membrane_triangle_list[i][2]][2])/3.0;
        
        for (int ch=Actin_Membrane_shared_num_of_Nodes;ch<Actin_num_of_Nodes;ch++)
        {
            // checking the interaction conditions:
            actin_membrane_distance_amplitude =sqrt( (membrane_triangle_COM_position[0] -Actin_Node_Position[ch][0] ) * (membrane_triangle_COM_position[0] - Actin_Node_Position[ch][0]) + (membrane_triangle_COM_position[1] - Actin_Node_Position[ch][1]) * (membrane_triangle_COM_position[1]- Actin_Node_Position[ch][1]) + (membrane_triangle_COM_position[2] - Actin_Node_Position[ch][2]) * (membrane_triangle_COM_position[2] - Actin_Node_Position[ch][2]));
            if (  actin_membrane_distance_amplitude < sqrt(0.43*Node_radius * 0.43*Node_radius + Actin_Membrane_Radius_of_Hard_Sphere_Interaction*Actin_Membrane_Radius_of_Hard_Sphere_Interaction )  )  // is solvent in proper sphere?
            {
                // cout<<"interactio"<<endl;
                AB[0]=Membrane_Node_Position[ Membrane_triangle_list[i][1]][0]-Membrane_Node_Position[ Membrane_triangle_list[i][0]][0];
                AB[1]=Membrane_Node_Position[ Membrane_triangle_list[i][1]][1]-Membrane_Node_Position[ Membrane_triangle_list[i][0]][1];
                AB[2]=Membrane_Node_Position[ Membrane_triangle_list[i][1]][2]-Membrane_Node_Position[ Membrane_triangle_list[i][0]][2];
                AC[0]=Membrane_Node_Position[ Membrane_triangle_list[i][2]][0]-Membrane_Node_Position[ Membrane_triangle_list[i][0]][0];
                AC[1]=Membrane_Node_Position[ Membrane_triangle_list[i][2]][1]-Membrane_Node_Position[ Membrane_triangle_list[i][0]][1];
                AC[2]=Membrane_Node_Position[ Membrane_triangle_list[i][2]][2]-Membrane_Node_Position[ Membrane_triangle_list[i][0]][2];
                crossvector(ABxAC,AB,AC);
                ABxAC[0]=ABxAC[0]*Membrane_Normal_direction[i][1];
                ABxAC[1]=ABxAC[1]*Membrane_Normal_direction[i][1];
                ABxAC[2]=ABxAC[2]*Membrane_Normal_direction[i][1];
                
                actin_triangle_distance_vector[0]=Actin_Node_Position[ch][0]-membrane_triangle_COM_position[0];
                actin_triangle_distance_vector[1]=Actin_Node_Position[ch][1]-membrane_triangle_COM_position[1];
                actin_triangle_distance_vector[2]=Actin_Node_Position[ch][2]-membrane_triangle_COM_position[2];
                
                perpendicular_distance=innerproduct(actin_triangle_distance_vector,ABxAC)/vectorlength(ABxAC);
                
                if    ( abs( perpendicular_distance )<Actin_Membrane_Radius_of_Hard_Sphere_Interaction )                {
                    
                    relevant_velocity[0] = Actin_Node_Velocity[ch][0]-membrane_triangle_COM_velocity[0];
                    relevant_velocity[1] = Actin_Node_Velocity[ch][1]-membrane_triangle_COM_velocity[1];
                    relevant_velocity[2] = Actin_Node_Velocity[ch][2]-membrane_triangle_COM_velocity[2];
                    
                    if (   innerproduct(relevant_velocity,ABxAC)*Membrane_Normal_direction[i][0]>0) {// is the actin node moving towards the triangle?
                        
                        for (int i1 = 0; i1 < 3; i1++)
                        {
                            Membrane_Node_Velocity[Membrane_triangle_list[i][i1]][0] +=(double) (2*Actin_Node_Mass / (Actin_Node_Mass + 3*Membrane_Node_Mass)) * (Actin_Node_Velocity[ch][0] - membrane_triangle_COM_velocity[0]);
                            Membrane_Node_Velocity[Membrane_triangle_list[i][i1]][1] +=(double) (2*Actin_Node_Mass / (Actin_Node_Mass + 3*Membrane_Node_Mass)) * (Actin_Node_Velocity[ch][1] - membrane_triangle_COM_velocity[1]);
                            Membrane_Node_Velocity[Membrane_triangle_list[i][i1]][2] +=(double) (2*Actin_Node_Mass / (Actin_Node_Mass + 3*Membrane_Node_Mass)) * (Actin_Node_Velocity[ch][2] - membrane_triangle_COM_velocity[2]);
                        }
                        
                        Actin_Node_Velocity[ch][0] +=  - (double) (6*Membrane_Node_Mass / (Actin_Node_Mass + 3*Membrane_Node_Mass)) * (Actin_Node_Velocity[ch][0]- membrane_triangle_COM_velocity[0]);
                        Actin_Node_Velocity[ch][1] +=  - (double) (6*Membrane_Node_Mass / (Actin_Node_Mass + 3*Membrane_Node_Mass)) * (Actin_Node_Velocity[ch][1] - membrane_triangle_COM_velocity[1]);
                        Actin_Node_Velocity[ch][2] += - (double) (6*Membrane_Node_Mass / (Actin_Node_Mass + 3*Membrane_Node_Mass)) * (Actin_Node_Velocity[ch][2] - membrane_triangle_COM_velocity[2]);
                    }
                }
            }
            
        }
        
    }
    
}
void Actin_Membrane_Barrier_2( double  Actin_Node_Position [][3], double  Actin_Node_Velocity [][3], double Membrane_Node_Position[][3], double Membrane_Node_Velocity[][3],  int Membrane_triangle_list[Membrane_num_of_Triangles][3], int Membrane_Normal_direction[Membrane_num_of_Triangles][2])
{
    double membrane_triangle_COM_position[3]; // coordinates to the centre of mass of the triangles
    double membrane_triangle_COM_velocity[3]; // velocity of the center of mass of the triangles
    double actin_membrane_distance_amplitude;// distance between solvent particle and com of triangle
    double ABxAC[3], AB[3], AC[3], ABxAC_unit_vector[3]; // normal vector of membrane
    double actin_triangle_distance_vector[3];
    double perpendicular_distance;
    double relevant_velocity[3];
    
    for (int i =0; i < Membrane_num_of_Triangles  ; i++)
    {
        membrane_triangle_COM_position[0]=(Membrane_Node_Position[ Membrane_triangle_list[i][0]][0] + Membrane_Node_Position[Membrane_triangle_list[i][1]][0] + Membrane_Node_Position[ Membrane_triangle_list[i][2]][0])/3.0;
        membrane_triangle_COM_position[1]=(Membrane_Node_Position[ Membrane_triangle_list[i][0]][1] + Membrane_Node_Position[Membrane_triangle_list[i][1]][1] + Membrane_Node_Position[ Membrane_triangle_list[i][2]][1])/3.0;
        membrane_triangle_COM_position[2]=(Membrane_Node_Position[ Membrane_triangle_list[i][0]][2] + Membrane_Node_Position[Membrane_triangle_list[i][1]][2] + Membrane_Node_Position[ Membrane_triangle_list[i][2]][2])/3.0;
        
        membrane_triangle_COM_velocity[0]=(Membrane_Node_Velocity[ Membrane_triangle_list[i][0]][0] + Membrane_Node_Velocity[Membrane_triangle_list[i][1]][0] + Membrane_Node_Velocity[ Membrane_triangle_list[i][2]][0])/3.0;
        membrane_triangle_COM_velocity[1]=(Membrane_Node_Velocity[ Membrane_triangle_list[i][0]][1] + Membrane_Node_Velocity[Membrane_triangle_list[i][1]][1] + Membrane_Node_Velocity[ Membrane_triangle_list[i][2]][1])/3.0;
        membrane_triangle_COM_velocity[2]=(Membrane_Node_Velocity[ Membrane_triangle_list[i][0]][2] + Membrane_Node_Velocity[Membrane_triangle_list[i][1]][2] + Membrane_Node_Velocity[ Membrane_triangle_list[i][2]][2])/3.0;
        
        for (int actin_counter=Actin_Membrane_shared_num_of_Nodes;actin_counter<Actin_num_of_Nodes;actin_counter++)
        {
            actin_membrane_distance_amplitude = sqrt( (membrane_triangle_COM_position[0] -Actin_Node_Position[actin_counter][0] ) * (membrane_triangle_COM_position[0] - Actin_Node_Position[actin_counter][0]) + (membrane_triangle_COM_position[1] - Actin_Node_Position[actin_counter][1]) * (membrane_triangle_COM_position[1]- Actin_Node_Position[actin_counter][1]) + (membrane_triangle_COM_position[2] - Actin_Node_Position[actin_counter][2]) * (membrane_triangle_COM_position[2] - Actin_Node_Position[actin_counter][2]));
            if (  actin_membrane_distance_amplitude < sqrt(0.43*Node_radius * 0.43*Node_radius + Actin_Membrane_Radius_of_Hard_Sphere_Interaction*Actin_Membrane_Radius_of_Hard_Sphere_Interaction )  )
            {
                AB[0]=Membrane_Node_Position[ Membrane_triangle_list[i][1]][0]-Membrane_Node_Position[ Membrane_triangle_list[i][0]][0];
                AB[1]=Membrane_Node_Position[ Membrane_triangle_list[i][1]][1]-Membrane_Node_Position[ Membrane_triangle_list[i][0]][1];
                AB[2]=Membrane_Node_Position[ Membrane_triangle_list[i][1]][2]-Membrane_Node_Position[ Membrane_triangle_list[i][0]][2];
                AC[0]=Membrane_Node_Position[ Membrane_triangle_list[i][2]][0]-Membrane_Node_Position[ Membrane_triangle_list[i][0]][0];
                AC[1]=Membrane_Node_Position[ Membrane_triangle_list[i][2]][1]-Membrane_Node_Position[ Membrane_triangle_list[i][0]][1];
                AC[2]=Membrane_Node_Position[ Membrane_triangle_list[i][2]][2]-Membrane_Node_Position[ Membrane_triangle_list[i][0]][2];
                crossvector(ABxAC,AB,AC);
                ABxAC[0]=ABxAC[0]*Membrane_Normal_direction[i][1];
                ABxAC[1]=ABxAC[1]*Membrane_Normal_direction[i][1];
                ABxAC[2]=ABxAC[2]*Membrane_Normal_direction[i][1];
                
                ABxAC_unit_vector[0]=ABxAC[0]/vectorlength(ABxAC);
                ABxAC_unit_vector[1]=ABxAC[1]/vectorlength(ABxAC);
                ABxAC_unit_vector[2]=ABxAC[2]/vectorlength(ABxAC);
                
                actin_triangle_distance_vector[0]=Actin_Node_Position[actin_counter][0]-membrane_triangle_COM_position[0];
                actin_triangle_distance_vector[1]=Actin_Node_Position[actin_counter][1]-membrane_triangle_COM_position[1];
                actin_triangle_distance_vector[2]=Actin_Node_Position[actin_counter][2]-membrane_triangle_COM_position[2];
                
                perpendicular_distance=innerproduct(actin_triangle_distance_vector, ABxAC)/vectorlength(ABxAC);
                
                relevant_velocity[0] = Actin_Node_Velocity[actin_counter][0]-membrane_triangle_COM_velocity[0];
                relevant_velocity[1] = Actin_Node_Velocity[actin_counter][1]-membrane_triangle_COM_velocity[1];
                relevant_velocity[2] = Actin_Node_Velocity[actin_counter][2]-membrane_triangle_COM_velocity[2];
                
                if    ( (abs( perpendicular_distance )<Actin_Membrane_Radius_of_Hard_Sphere_Interaction)
                       && (innerproduct(relevant_velocity,ABxAC)*Membrane_Normal_direction[i][0]>0)
                       )
                {
                    double Actin_velocity_N_new, Actin_velocity_N, Membrane_triangle_COM_velocity_N_new, Membrane_triangle_COM_velocity_N;
                    
                    Actin_velocity_N=innerproduct(Actin_Node_Velocity[actin_counter], ABxAC_unit_vector);
                    Membrane_triangle_COM_velocity_N=innerproduct(membrane_triangle_COM_velocity, ABxAC_unit_vector);
                    
                    Actin_velocity_N_new=(Actin_velocity_N*(Actin_Node_Mass-3*Membrane_Node_Mass)+2.0*3*Membrane_Node_Mass*Membrane_triangle_COM_velocity_N)/(Actin_Node_Mass+3*Membrane_Node_Mass);
                    Membrane_triangle_COM_velocity_N_new=(Membrane_triangle_COM_velocity_N*(3*Membrane_Node_Mass-Actin_Node_Mass)+2.0*Actin_Node_Mass*Actin_velocity_N)/(Actin_Node_Mass+3*Membrane_Node_Mass);
                    
                    Actin_Node_Velocity[actin_counter][0]+= (-Actin_velocity_N + Actin_velocity_N_new)*ABxAC_unit_vector[0];
                    Actin_Node_Velocity[actin_counter][1]+= (-Actin_velocity_N + Actin_velocity_N_new)*ABxAC_unit_vector[1];
                    Actin_Node_Velocity[actin_counter][2]+= (-Actin_velocity_N + Actin_velocity_N_new)*ABxAC_unit_vector[2];
                    
                    for (int i1 = 0; i1 < 3; i1++)
                    {
                        Membrane_Node_Velocity[Membrane_triangle_list[i][i1]][0] += (-Membrane_triangle_COM_velocity_N + Membrane_triangle_COM_velocity_N_new)*ABxAC_unit_vector[0];
                        Membrane_Node_Velocity[Membrane_triangle_list[i][i1]][1] += (-Membrane_triangle_COM_velocity_N + Membrane_triangle_COM_velocity_N_new)*ABxAC_unit_vector[1];
                        Membrane_Node_Velocity[Membrane_triangle_list[i][i1]][2] += (-Membrane_triangle_COM_velocity_N + Membrane_triangle_COM_velocity_N_new)*ABxAC_unit_vector[2];
                    }
                }//END OF:  if    ( (abs( perpendicular_distance )<Actin_Membrane_Radius_of_Hard_Sphere_Interaction) &&
            }//END OF: if (  actin_membrane_distance_amplitude < sqrt(0.43*a * 0.43*a +
            
        }//END OF: for (int actin_counter=Actin_Membrane_shared_num_of_Nodes;actin_counter<Actin_num_of_Nodes;
        
    }//END OF: for (int i =0; i < Membrane_num_of_Triangles  ; i++)
    
}

void Nucleus_Membrane_Barrier(int Nucleus_Membrane_list_of_Nodes[], double Membrane_Node_Position[][3], double Membrane_Node_Velocity[][3], int Membrane_triangle_list[Membrane_num_of_Triangles][3], int Membrane_Normal_direction[Membrane_num_of_Triangles][2], int Outer_Membrane_num_of_triangles, int  Outer_Membrane_num_of_Nodes, int Membrane_num_of_Nodes)
{
    double membrane_triangle_COM_position[3]; // place of center of mass of triangle
    double membrane_triangle_COM_velocity[3]; // velocity of center of mass of triangle
    double vcomnew[3]; // velocity of center of mass of triangle
    double Nucleus_membrane_distance_amplitude;// distance between solvent particle and com of triangle
    double ABxAC[3],AB[3],AC[3]; // normal vector of membrane
    double Nucleus_triangle_distance_vector[3];
    double perpendicular_distance;
    double relevant_velocity[3];
    int temp_nucleus_membrane_node;
    
    for (int Membrane_counter =0; Membrane_counter < Outer_Membrane_num_of_triangles  ; Membrane_counter++)  // for each triangle:
    {
        //Calculating the coordinates and velocity of the centre of the Membrane triangle
        membrane_triangle_COM_position[0]=(Membrane_Node_Position[ Membrane_triangle_list[Membrane_counter][0]][0] + Membrane_Node_Position[Membrane_triangle_list[Membrane_counter][1]][0] + Membrane_Node_Position[ Membrane_triangle_list[Membrane_counter][2]][0])/3.0;
        membrane_triangle_COM_position[1]=(Membrane_Node_Position[ Membrane_triangle_list[Membrane_counter][0]][1] + Membrane_Node_Position[Membrane_triangle_list[Membrane_counter][1]][1] + Membrane_Node_Position[ Membrane_triangle_list[Membrane_counter][2]][1])/3.0;
        membrane_triangle_COM_position[2]=(Membrane_Node_Position[ Membrane_triangle_list[Membrane_counter][0]][2] + Membrane_Node_Position[Membrane_triangle_list[Membrane_counter][1]][2] + Membrane_Node_Position[ Membrane_triangle_list[Membrane_counter][2]][2])/3.0;
        
        membrane_triangle_COM_velocity[0]=(Membrane_Node_Velocity[ Membrane_triangle_list[Membrane_counter][0]][0] + Membrane_Node_Velocity[Membrane_triangle_list[Membrane_counter][1]][0] + Membrane_Node_Velocity[ Membrane_triangle_list[Membrane_counter][2]][0])/3.0;
        membrane_triangle_COM_velocity[1]=(Membrane_Node_Velocity[ Membrane_triangle_list[Membrane_counter][0]][1] + Membrane_Node_Velocity[Membrane_triangle_list[Membrane_counter][1]][1] + Membrane_Node_Velocity[ Membrane_triangle_list[Membrane_counter][2]][1])/3.0;
        membrane_triangle_COM_velocity[2]=(Membrane_Node_Velocity[ Membrane_triangle_list[Membrane_counter][0]][2] + Membrane_Node_Velocity[Membrane_triangle_list[Membrane_counter][1]][2] + Membrane_Node_Velocity[ Membrane_triangle_list[Membrane_counter][2]][2])/3.0;
        
        for (int j=0;j<Membrane_num_of_Nodes-Outer_Membrane_num_of_Nodes;j++)
        {
            temp_nucleus_membrane_node=Nucleus_Membrane_list_of_Nodes[j]; // We take a Nucleus node from the list
            
            Nucleus_membrane_distance_amplitude =sqrt( (membrane_triangle_COM_position[0] - Membrane_Node_Position[temp_nucleus_membrane_node][0]) * (membrane_triangle_COM_position[0] - Membrane_Node_Position[temp_nucleus_membrane_node][0]) + (membrane_triangle_COM_position[1] -  Membrane_Node_Position[temp_nucleus_membrane_node][1]) * (membrane_triangle_COM_position[1]- Membrane_Node_Position[temp_nucleus_membrane_node][1]) + (membrane_triangle_COM_position[2] - Membrane_Node_Position[temp_nucleus_membrane_node][2]) * (membrane_triangle_COM_position[2] - Membrane_Node_Position[temp_nucleus_membrane_node][2]));
            
            if ( Nucleus_membrane_distance_amplitude < sqrt(0.43*Node_radius * 0.43*Node_radius + 4*Nucleus_Membrane_Radius_of_Hard_Sphere_Interaction*Nucleus_Membrane_Radius_of_Hard_Sphere_Interaction ) ){
                
                AB[0]=Membrane_Node_Position[ Membrane_triangle_list[Membrane_counter][1]][0]-Membrane_Node_Position[ Membrane_triangle_list[Membrane_counter][0]][0];
                AB[1]=Membrane_Node_Position[ Membrane_triangle_list[Membrane_counter][1]][1]-Membrane_Node_Position[ Membrane_triangle_list[Membrane_counter][0]][1];
                AB[2]=Membrane_Node_Position[ Membrane_triangle_list[Membrane_counter][1]][2]-Membrane_Node_Position[ Membrane_triangle_list[Membrane_counter][0]][2];
                AC[0]=Membrane_Node_Position[ Membrane_triangle_list[Membrane_counter][2]][0]-Membrane_Node_Position[ Membrane_triangle_list[Membrane_counter][0]][0];
                AC[1]=Membrane_Node_Position[ Membrane_triangle_list[Membrane_counter][2]][1]-Membrane_Node_Position[ Membrane_triangle_list[Membrane_counter][0]][1];
                AC[2]=Membrane_Node_Position[ Membrane_triangle_list[Membrane_counter][2]][2]-Membrane_Node_Position[ Membrane_triangle_list[Membrane_counter][0]][2];
                crossvector(ABxAC,AB,AC);
                ABxAC[0]=ABxAC[0]*Membrane_Normal_direction[Membrane_counter][1];
                ABxAC[1]=ABxAC[1]*Membrane_Normal_direction[Membrane_counter][1];
                ABxAC[2]=ABxAC[2]*Membrane_Normal_direction[Membrane_counter][1];
                
                Nucleus_triangle_distance_vector[0]=  Membrane_Node_Position[temp_nucleus_membrane_node][0]-membrane_triangle_COM_position[0];
                Nucleus_triangle_distance_vector[1]=  Membrane_Node_Position[temp_nucleus_membrane_node][1]-membrane_triangle_COM_position[1];
                Nucleus_triangle_distance_vector[2]=  Membrane_Node_Position[temp_nucleus_membrane_node][2]-membrane_triangle_COM_position[2];
                perpendicular_distance=innerproduct(Nucleus_triangle_distance_vector,ABxAC)/vectorlength(ABxAC);
                
                if    (  abs( perpendicular_distance )<2*Nucleus_Membrane_Radius_of_Hard_Sphere_Interaction  ){
                    
                    relevant_velocity[0] =  Membrane_Node_Velocity[temp_nucleus_membrane_node][0]-membrane_triangle_COM_velocity[0];
                    relevant_velocity[1] =  Membrane_Node_Velocity[temp_nucleus_membrane_node][1]-membrane_triangle_COM_velocity[1];
                    relevant_velocity[2] =  Membrane_Node_Velocity[temp_nucleus_membrane_node][2]-membrane_triangle_COM_velocity[2];
                    
                    if (   innerproduct(relevant_velocity,ABxAC)*perpendicular_distance<0    ) // is solvent moving toward triangle?
                    {
                        Membrane_Node_Velocity[temp_nucleus_membrane_node][0] += - (double) (6*Membrane_Node_Mass / (Membrane_Node_Mass + 3*Membrane_Node_Mass)) * (Membrane_Node_Velocity[temp_nucleus_membrane_node][0]-  membrane_triangle_COM_velocity[0]);
                        Membrane_Node_Velocity[temp_nucleus_membrane_node][1] += - (double) (6*Membrane_Node_Mass / (Membrane_Node_Mass + 3*Membrane_Node_Mass)) * (Membrane_Node_Velocity[temp_nucleus_membrane_node][1]- membrane_triangle_COM_velocity[1]);
                        Membrane_Node_Velocity[temp_nucleus_membrane_node][2] += - (double) (6*Membrane_Node_Mass / (Membrane_Node_Mass + 3*Membrane_Node_Mass)) * (Membrane_Node_Velocity[temp_nucleus_membrane_node][2] - membrane_triangle_COM_velocity[2]);
                        
                        vcomnew[0]=  (double) (2*Membrane_Node_Mass / (Membrane_Node_Mass + 3*Membrane_Node_Mass)) * (Membrane_Node_Velocity[temp_nucleus_membrane_node][0] - membrane_triangle_COM_velocity[0]);
                        vcomnew[1]=  (double) (2*Membrane_Node_Mass / (Membrane_Node_Mass + 3*Membrane_Node_Mass)) * (Membrane_Node_Velocity[temp_nucleus_membrane_node][1] - membrane_triangle_COM_velocity[1]);
                        vcomnew[2]=  (double) (2*Membrane_Node_Mass / (Membrane_Node_Mass + 3*Membrane_Node_Mass)) * (Membrane_Node_Velocity[temp_nucleus_membrane_node][2] - membrane_triangle_COM_velocity[2]);
                        
                        for (int i1 = 0; i1 < 3; i1++)
                        {
                            Membrane_Node_Velocity[Membrane_triangle_list[Membrane_counter][i1]][0] += vcomnew[0];
                            Membrane_Node_Velocity[Membrane_triangle_list[Membrane_counter][i1]][1] += vcomnew[1];
                            Membrane_Node_Velocity[Membrane_triangle_list[Membrane_counter][i1]][2] += vcomnew[2];
                        }
                        
                    }
                    
                }
            }
            
            
        }
        
    }
    
}
void Nucleus_Membrane_Barrier_2(int Nucleus_Membrane_list_of_Nodes[], double Membrane_Node_Position[][3], double Membrane_Node_Velocity[][3], int Membrane_triangle_list[Membrane_num_of_Triangles][3], int Membrane_Normal_direction[Membrane_num_of_Triangles][2], int Outer_Membrane_num_of_triangles, int  Outer_Membrane_num_of_Nodes, int Membrane_num_of_Nodes)
{
    double membrane_triangle_COM_position[3]; // place of center of mass of triangle
    double membrane_triangle_COM_velocity[3]; // velocity of center of mass of triangle
    double Nucleus_membrane_distance_amplitude;// distance between solvent particle and com of triangle
    double ABxAC[3], AB[3], AC[3], ABxAC_unit_vector[3]; // normal vector of membrane
    double Nucleus_triangle_distance_vector[3];
    double perpendicular_distance;
    double relevant_velocity[3];
    int temp_nucleus_membrane_node;
    
    for (int Membrane_counter =0; Membrane_counter < Outer_Membrane_num_of_triangles  ; Membrane_counter++)
    {
        //Calculating the coordinates and velocity of the centre of the Membrane triangle
        membrane_triangle_COM_position[0]=(Membrane_Node_Position[ Membrane_triangle_list[Membrane_counter][0]][0] + Membrane_Node_Position[Membrane_triangle_list[Membrane_counter][1]][0] + Membrane_Node_Position[ Membrane_triangle_list[Membrane_counter][2]][0])/3.0;
        membrane_triangle_COM_position[1]=(Membrane_Node_Position[ Membrane_triangle_list[Membrane_counter][0]][1] + Membrane_Node_Position[Membrane_triangle_list[Membrane_counter][1]][1] + Membrane_Node_Position[ Membrane_triangle_list[Membrane_counter][2]][1])/3.0;
        membrane_triangle_COM_position[2]=(Membrane_Node_Position[ Membrane_triangle_list[Membrane_counter][0]][2] + Membrane_Node_Position[Membrane_triangle_list[Membrane_counter][1]][2] + Membrane_Node_Position[ Membrane_triangle_list[Membrane_counter][2]][2])/3.0;
        
        membrane_triangle_COM_velocity[0]=(Membrane_Node_Velocity[ Membrane_triangle_list[Membrane_counter][0]][0] + Membrane_Node_Velocity[Membrane_triangle_list[Membrane_counter][1]][0] + Membrane_Node_Velocity[ Membrane_triangle_list[Membrane_counter][2]][0])/3.0;
        membrane_triangle_COM_velocity[1]=(Membrane_Node_Velocity[ Membrane_triangle_list[Membrane_counter][0]][1] + Membrane_Node_Velocity[Membrane_triangle_list[Membrane_counter][1]][1] + Membrane_Node_Velocity[ Membrane_triangle_list[Membrane_counter][2]][1])/3.0;
        membrane_triangle_COM_velocity[2]=(Membrane_Node_Velocity[ Membrane_triangle_list[Membrane_counter][0]][2] + Membrane_Node_Velocity[Membrane_triangle_list[Membrane_counter][1]][2] + Membrane_Node_Velocity[ Membrane_triangle_list[Membrane_counter][2]][2])/3.0;
        
        for (int j=0;j<Membrane_num_of_Nodes-Outer_Membrane_num_of_Nodes;j++)
        {
            temp_nucleus_membrane_node=Nucleus_Membrane_list_of_Nodes[j]; // We take a Nucleus node from the list
            
            Nucleus_membrane_distance_amplitude =sqrt( (membrane_triangle_COM_position[0] - Membrane_Node_Position[temp_nucleus_membrane_node][0]) * (membrane_triangle_COM_position[0] - Membrane_Node_Position[temp_nucleus_membrane_node][0]) + (membrane_triangle_COM_position[1] -  Membrane_Node_Position[temp_nucleus_membrane_node][1]) * (membrane_triangle_COM_position[1]- Membrane_Node_Position[temp_nucleus_membrane_node][1]) + (membrane_triangle_COM_position[2] - Membrane_Node_Position[temp_nucleus_membrane_node][2]) * (membrane_triangle_COM_position[2] - Membrane_Node_Position[temp_nucleus_membrane_node][2]));
            
            if ( Nucleus_membrane_distance_amplitude < sqrt(0.43*Node_radius * 0.43*Node_radius + 4*Nucleus_Membrane_Radius_of_Hard_Sphere_Interaction*Nucleus_Membrane_Radius_of_Hard_Sphere_Interaction ) )
            {
                
                AB[0]=Membrane_Node_Position[ Membrane_triangle_list[Membrane_counter][1]][0]-Membrane_Node_Position[ Membrane_triangle_list[Membrane_counter][0]][0];
                AB[1]=Membrane_Node_Position[ Membrane_triangle_list[Membrane_counter][1]][1]-Membrane_Node_Position[ Membrane_triangle_list[Membrane_counter][0]][1];
                AB[2]=Membrane_Node_Position[ Membrane_triangle_list[Membrane_counter][1]][2]-Membrane_Node_Position[ Membrane_triangle_list[Membrane_counter][0]][2];
                AC[0]=Membrane_Node_Position[ Membrane_triangle_list[Membrane_counter][2]][0]-Membrane_Node_Position[ Membrane_triangle_list[Membrane_counter][0]][0];
                AC[1]=Membrane_Node_Position[ Membrane_triangle_list[Membrane_counter][2]][1]-Membrane_Node_Position[ Membrane_triangle_list[Membrane_counter][0]][1];
                AC[2]=Membrane_Node_Position[ Membrane_triangle_list[Membrane_counter][2]][2]-Membrane_Node_Position[ Membrane_triangle_list[Membrane_counter][0]][2];
                crossvector(ABxAC,AB,AC);
                ABxAC[0]=ABxAC[0]*Membrane_Normal_direction[Membrane_counter][1];
                ABxAC[1]=ABxAC[1]*Membrane_Normal_direction[Membrane_counter][1];
                ABxAC[2]=ABxAC[2]*Membrane_Normal_direction[Membrane_counter][1];
                
                ABxAC_unit_vector[0]=ABxAC[0]/vectorlength(ABxAC);
                ABxAC_unit_vector[1]=ABxAC[1]/vectorlength(ABxAC);
                ABxAC_unit_vector[2]=ABxAC[2]/vectorlength(ABxAC);
                
                Nucleus_triangle_distance_vector[0]=  Membrane_Node_Position[temp_nucleus_membrane_node][0]-membrane_triangle_COM_position[0];
                Nucleus_triangle_distance_vector[1]=  Membrane_Node_Position[temp_nucleus_membrane_node][1]-membrane_triangle_COM_position[1];
                Nucleus_triangle_distance_vector[2]=  Membrane_Node_Position[temp_nucleus_membrane_node][2]-membrane_triangle_COM_position[2];
                perpendicular_distance=innerproduct(Nucleus_triangle_distance_vector,ABxAC)/vectorlength(ABxAC);
                
                relevant_velocity[0] =  Membrane_Node_Velocity[temp_nucleus_membrane_node][0]-membrane_triangle_COM_velocity[0];
                relevant_velocity[1] =  Membrane_Node_Velocity[temp_nucleus_membrane_node][1]-membrane_triangle_COM_velocity[1];
                relevant_velocity[2] =  Membrane_Node_Velocity[temp_nucleus_membrane_node][2]-membrane_triangle_COM_velocity[2];
                
                if    (  (abs( perpendicular_distance )<2*Nucleus_Membrane_Radius_of_Hard_Sphere_Interaction) && (innerproduct(relevant_velocity,ABxAC)*Membrane_Normal_direction[Membrane_counter][0]>0) )
                {
                    double Nucleus_velocity_N_new, Nucleus_velocity_N, Membrane_triangle_COM_velocity_N_new, Membrane_triangle_COM_velocity_N;
                    
                    Nucleus_velocity_N=innerproduct(Membrane_Node_Velocity[temp_nucleus_membrane_node], ABxAC_unit_vector);
                    Membrane_triangle_COM_velocity_N=innerproduct(membrane_triangle_COM_velocity, ABxAC_unit_vector);
                    
                    Nucleus_velocity_N_new=(Nucleus_velocity_N*(Membrane_Node_Mass-3*Membrane_Node_Mass)+2.0*3*Membrane_Node_Mass*Membrane_triangle_COM_velocity_N)/(Membrane_Node_Mass+3*Membrane_Node_Mass);
                    Membrane_triangle_COM_velocity_N_new=(Membrane_triangle_COM_velocity_N*(3*Membrane_Node_Mass-Membrane_Node_Mass)+2.0*Membrane_Node_Mass*Nucleus_velocity_N)/(Membrane_Node_Mass+3*Membrane_Node_Mass);
                    
                    Membrane_Node_Velocity[temp_nucleus_membrane_node][0]+= (-Nucleus_velocity_N + Nucleus_velocity_N_new)*ABxAC_unit_vector[0];
                    Membrane_Node_Velocity[temp_nucleus_membrane_node][1]+= (-Nucleus_velocity_N + Nucleus_velocity_N_new)*ABxAC_unit_vector[1];
                    Membrane_Node_Velocity[temp_nucleus_membrane_node][2]+= (-Nucleus_velocity_N + Nucleus_velocity_N_new)*ABxAC_unit_vector[2];
                    
                    for (int i1 = 0; i1 < 3; i1++)
                    {
                        Membrane_Node_Velocity[Membrane_triangle_list[Membrane_counter][i1]][0] += (-Membrane_triangle_COM_velocity_N + Membrane_triangle_COM_velocity_N_new)*ABxAC_unit_vector[0];
                        Membrane_Node_Velocity[Membrane_triangle_list[Membrane_counter][i1]][1] += (-Membrane_triangle_COM_velocity_N + Membrane_triangle_COM_velocity_N_new)*ABxAC_unit_vector[1];
                        Membrane_Node_Velocity[Membrane_triangle_list[Membrane_counter][i1]][2] += (-Membrane_triangle_COM_velocity_N + Membrane_triangle_COM_velocity_N_new)*ABxAC_unit_vector[2];
                    }
                }//END OF: if    (  (abs( perpendicular_distance )<2*Nucleus_Membrane_Radius_of_Hard_Sphere_Interaction) &&
            }//END OF: if ( Nucleus_membrane_distance_amplitude < sqrt(0.43*a * 0.43*a + 4*
            
            
        }//END OF: for (int j=0;j<Membrane_num_of_Nodes-Outer_Membrane_num_of_Nodes;j++)
        
    }//END OF: for (int Membrane_counter =0; Membrane_counter < Outer_Membrane_num_of_triangles  ; Membrane_counter++)
    
}


//_______________________________chromatin functions





void Chromatin_Force_calculator(double (&Chromatin_Bead_Position)[Chromatin_num_of_Beads][3],double (&Chromatin_Bead_Velocity)[Chromatin_num_of_Beads][3],double (&Chromatin_Bead_Force)[Chromatin_num_of_Beads][3], double &Total_Potential_Energy)  // calculate force of inside the chromatin beads on each other
{
    double deltax, deltay, deltaz, temp_distance, force;// defined below in "for loop" in detail
    double  temp_vector_1[3],temp_vector_2[3],fi[3],fi1[3],fi2[3],temp_vector_1_length,temp_vector_2_length;
    // Simple spring force for all members of the chain
    for (int temp_chain_counter=0; temp_chain_counter<Chromatin_num_of_chains; temp_chain_counter++)
    {
        for(int i=temp_chain_counter*(Chromatin_num_of_Beads/Chromatin_num_of_chains)  ;i< (temp_chain_counter+1)*(Chromatin_num_of_Beads/Chromatin_num_of_chains) -2; i++ )  //i goes through all the beads on one chain
        {
            deltax=Chromatin_Bead_Position[i][0]-Chromatin_Bead_Position[i+1][0];
            deltay=Chromatin_Bead_Position[i][1]-Chromatin_Bead_Position[i+1][1];
            deltaz=Chromatin_Bead_Position[i][2]-Chromatin_Bead_Position[i+1][2];
            temp_distance=sqrt(deltaz*deltaz+deltay*deltay+deltax*deltax);
            //*******************************************************************************************************
            /*BUG
             |---\   |    |  /---\
             |    |  |    |  |
             |---<   |    |  |  -\
             |    |  |    |  |   |
             |---/   \----/  \---/
             */
            //*******************************************************************************************************
            //***************** Potential BUG: Do not know why the Chromatin_Scaling_Factor should be included in these calculations*********************************
            //*******************************************************************************************************
            force=Chromatin_spring_coefficient*(temp_distance- sigmachromatin*Chromatin_Scaling_Factor);
            //            cout<<"temp_distance- sigmachromatin*Chromatin_Scaling_Factor in ch= "<<temp_distance- sigmachromatin*Chromatin_Scaling_Factor<<endl;
            if (energy_calculation_flag==1.0) {
                Total_Potential_Energy+=0.5*Chromatin_spring_coefficient*(temp_distance- sigmachromatin*Chromatin_Scaling_Factor)*(temp_distance- sigmachromatin*Chromatin_Scaling_Factor);
                //            cout<<"Total_Potential_Energy in ch= "<<Total_Potential_Energy<<endl;
            }
            
            
            Chromatin_Bead_Force[i][0] += +force*deltax/temp_distance ; // force of springs and dashes
            Chromatin_Bead_Force[i][1] += +force*deltay/temp_distance  ;// force of springs and dashes
            Chromatin_Bead_Force[i][2] += +force*deltaz/temp_distance ;// force of springs and dashes
            
            Chromatin_Bead_Force[i+1][0] += -force *  deltax/temp_distance  ;
            Chromatin_Bead_Force[i+1][1] += -force *  deltay/temp_distance  ;
            Chromatin_Bead_Force[i+1][2] += -force *  deltaz/temp_distance  ;
            
            // Bending Force
            //*******************************************************************************************************
            /*BUG
             |---\   |    |  /---\
             |    |  |    |  |
             |---<   |    |  |  -\
             |    |  |    |  |   |
             |---/   \----/  \---/
             */
            //*******************************************************************************************************
            //***************** Potential BUG: The 'Chromatin_bending_coefficient' was set to be zero when I got the code. I am not sure why we do not want the chromatin to bend. It will colaps it it does not bend, right? *********************************
            //*******************************************************************************************************
            //*******************************************************************************************************
            /*BUG
             |---\   |    |  /---\
             |    |  |    |  |
             |---<   |    |  |  -\
             |    |  |    |  |   |
             |---/   \----/  \---/
             */
            //*******************************************************************************************************
            //***************** BUG: Also I will replace the current code with a spring! *********************************
            //*******************************************************************************************************
            temp_vector_1[0] = Chromatin_Bead_Position[i+1][0]-Chromatin_Bead_Position[i][0];
            temp_vector_1[1] = Chromatin_Bead_Position[i+1][1]-Chromatin_Bead_Position[i][1];
            temp_vector_1[2] = Chromatin_Bead_Position[i+1][2]-Chromatin_Bead_Position[i][2];
            temp_vector_1_length=vectorlength(temp_vector_1);
            
            temp_vector_2[0] = Chromatin_Bead_Position[i+2][0]-Chromatin_Bead_Position[i+1][0];
            temp_vector_2[1] = Chromatin_Bead_Position[i+2][1]-Chromatin_Bead_Position[i+1][1];
            temp_vector_2[2] = Chromatin_Bead_Position[i+2][2]-Chromatin_Bead_Position[i+1][2];
            temp_vector_2_length=vectorlength(temp_vector_2);
            
            if (Chromatin_bending_coefficient!=0) {
                fi[0] = Chromatin_bending_coefficient* (-temp_vector_2[0]+innerproduct(temp_vector_1,temp_vector_2) *temp_vector_1[0] /( temp_vector_1_length*temp_vector_1_length )   )/( 2.0*temp_vector_1_length*temp_vector_2_length );
                fi[1] = Chromatin_bending_coefficient* (-temp_vector_2[1]+innerproduct(temp_vector_1,temp_vector_2 )*temp_vector_1[1] /( temp_vector_1_length*temp_vector_1_length )   )/( 2.0*temp_vector_1_length*temp_vector_2_length );
                fi[2] = Chromatin_bending_coefficient* (-temp_vector_2[2]+innerproduct(temp_vector_1,temp_vector_2) *temp_vector_1[2] /( temp_vector_1_length*temp_vector_1_length )   )/( 2.0*temp_vector_1_length*temp_vector_2_length );
                
                fi2[0]=Chromatin_bending_coefficient* (temp_vector_1[0]  -    innerproduct(temp_vector_1,temp_vector_2) *temp_vector_2[0] /( temp_vector_2_length*temp_vector_2_length )   )/( 2.0*temp_vector_1_length*temp_vector_2_length );
                fi2[1]=Chromatin_bending_coefficient* (temp_vector_1[1]  -    innerproduct(temp_vector_1,temp_vector_2 )*temp_vector_2[1] /( temp_vector_2_length*temp_vector_2_length )   )/( 2.0*temp_vector_1_length*temp_vector_2_length );
                fi2[2]=Chromatin_bending_coefficient* (temp_vector_1[2]  -    innerproduct(temp_vector_1,temp_vector_2) *temp_vector_2[2] /( temp_vector_2_length*temp_vector_2_length )   )/( 2.0*temp_vector_1_length*temp_vector_2_length );
                
                fi1[0]=-fi[0]-fi2[0];
                fi1[1]=-fi[1]-fi2[1];
                fi1[2]=-fi[2]-fi2[2];
                
                Chromatin_Bead_Force[i][0] += -  fi[0];
                Chromatin_Bead_Force[i][1] += -  fi[1];
                Chromatin_Bead_Force[i][2] += -  fi[2];
                
                Chromatin_Bead_Force[i+1][0] += -  fi1[0];
                Chromatin_Bead_Force[i+1][1] += -  fi1[1];
                Chromatin_Bead_Force[i+1][2] += -  fi1[2];
                
                Chromatin_Bead_Force[i+2][0] += -  fi2[0];
                Chromatin_Bead_Force[i+2][1] += -  fi2[1];
                Chromatin_Bead_Force[i+2][2] += -  fi2[2];
            }
            
        } //End of 'for' i
        
        
        //Since the previous loop is up to 'i=(temp_chain_counter+1)*(Chromatin_num_of_Beads/Chromatin_num_of_chains)-3' we do add the last spring force calculations of the chain manually here, (for 'i=(temp_chain_counter+1)*(Chromatin_num_of_Beads/Chromatin_num_of_chains)-2')
        deltax=Chromatin_Bead_Position[(temp_chain_counter+1)*(Chromatin_num_of_Beads/Chromatin_num_of_chains)-2][0]-Chromatin_Bead_Position[(temp_chain_counter+1)*(Chromatin_num_of_Beads/Chromatin_num_of_chains) -1][0];
        deltay=Chromatin_Bead_Position[(temp_chain_counter+1)*(Chromatin_num_of_Beads/Chromatin_num_of_chains)-2][1]-Chromatin_Bead_Position[(temp_chain_counter+1)*(Chromatin_num_of_Beads/Chromatin_num_of_chains) -1][1];
        deltaz=Chromatin_Bead_Position[(temp_chain_counter+1)*(Chromatin_num_of_Beads/Chromatin_num_of_chains)-2][2]-Chromatin_Bead_Position[(temp_chain_counter+1)*(Chromatin_num_of_Beads/Chromatin_num_of_chains) -1][2];
        temp_distance=sqrt(deltaz*deltaz+deltay*deltay+deltax*deltax);
        //*******************************************************************************************************
        /*BUG
         |---\   |    |  /---\
         |    |  |    |  |
         |---<   |    |  |  -\
         |    |  |    |  |   |
         |---/   \----/  \---/
         */
        //*******************************************************************************************************
        //***************** Potential BUG: Do not know why the Chromatin_Scaling_Factor should be included in these calculations*********************************
        //*******************************************************************************************************
        force=Chromatin_spring_coefficient*(temp_distance- sigmachromatin*Chromatin_Scaling_Factor);
        if (energy_calculation_flag==1.0) {
            Total_Potential_Energy+=0.5*Chromatin_spring_coefficient*(temp_distance- sigmachromatin*Chromatin_Scaling_Factor)*(temp_distance- sigmachromatin*Chromatin_Scaling_Factor);
            
        }
        
        Chromatin_Bead_Force[(temp_chain_counter+1)*(Chromatin_num_of_Beads/Chromatin_num_of_chains)-2][0] += +force*deltax/temp_distance ; // force of springs and dashes
        Chromatin_Bead_Force[(temp_chain_counter+1)*(Chromatin_num_of_Beads/Chromatin_num_of_chains)-2][1] += +force*deltay/temp_distance  ;// force of springs and dashes
        Chromatin_Bead_Force[(temp_chain_counter+1)*(Chromatin_num_of_Beads/Chromatin_num_of_chains)-2][2] += +force*deltaz/temp_distance ;// force of springs and dashes
        
        Chromatin_Bead_Force[(temp_chain_counter+1)*(Chromatin_num_of_Beads/Chromatin_num_of_chains) -1][0] += -force *  deltax/temp_distance  ;
        Chromatin_Bead_Force[(temp_chain_counter+1)*(Chromatin_num_of_Beads/Chromatin_num_of_chains) -1][1] += -force *  deltay/temp_distance  ;
        Chromatin_Bead_Force[(temp_chain_counter+1)*(Chromatin_num_of_Beads/Chromatin_num_of_chains) -1][2] += -force *  deltaz/temp_distance  ;
    }//end of 'for' temp_chain_counter
    
    
    
    
}


void hardsphereforcechromatin(double Chromatin_Bead_Position[Chromatin_num_of_Beads][3], double Chromatin_Bead_Velocity[Chromatin_num_of_Beads][3], double Chromatin_Bead_Force[Chromatin_num_of_Beads][3], bool chromatin_contant_matrix_calculation_flag, float CmLaststeps[Chromatin_num_of_Beads][Chromatin_num_of_Beads], float Cm[Chromatin_num_of_Beads][Chromatin_num_of_Beads], float Cm1[Chromatin_num_of_Beads][Chromatin_num_of_Beads], float Cm3[Chromatin_num_of_Beads][Chromatin_num_of_Beads], double Dissimilaritycut)  // calculate force of inside the chromatin beads on each other
{
    
    double deltax1,deltay1,deltaz1,distance1 = 0.0;// defined below in "for loop" in detail
    
    double ljmatrix[Chromatin_num_of_Beads][Chromatin_num_of_Beads],flj[3],r2,r6;     // (force amplitute)/r     of l-j will kept in jfmatrix
    
    // update_spatial_contact_matrix
    double contact_coordinates[3];
    double mesh_size= 10.0* double (Nucleus_Membrane_radius) / double (size_spatialcontactmatrix) ;
    int ii,jj,kk;
    
    /// L-J   force :
    if (chromatin_contant_matrix_calculation_flag==true) {
        for(int i=0;i<Chromatin_num_of_Beads ;i++)  // initialize lhmatrix
        {
            for(int j=0;j<Chromatin_num_of_Beads ;j++)
            {
                ljmatrix[i][j]=0.0;
                Cm[i][j]=0;
            }
        }
    }
    
    
    
    
    for (int j=0 ; j< Chromatin_num_of_Beads ; j++)  //calculate the amplitute of force between beads
    {
        for (int i=0 ; i<Chromatin_num_of_Beads ;i++)      // i dnotes rows
        {
            if(i != j   )
            {
                deltax1=Chromatin_Bead_Position[i][0]-Chromatin_Bead_Position[j][0];// delta x  betwwn i and j th beads
                deltay1=Chromatin_Bead_Position[i][1]-Chromatin_Bead_Position[j][1];// delta y  betwwn i and j th beads
                deltaz1=Chromatin_Bead_Position[i][2]-Chromatin_Bead_Position[j][2];// delta z  betwwn i and j th beads
                distance1=sqrt(deltax1*deltax1+deltay1*deltay1+deltaz1*deltaz1); // distance btween i th and i+1 th  bead
                
                if( distance1<    1.122462048*sigmachromatin*Chromatin_Scaling_Factor  && distance1 != 0.0 )
                {
                    r2=distance1*distance1/(sigmachromatin*sigmachromatin*Chromatin_Scaling_Factor*Chromatin_Scaling_Factor);
                    r6=r2*r2*r2;
                    ljmatrix[i][j]=-4.0*(  12.0/(r6*r6)- 6.0/r6   )/distance1;
                    if( ljmatrix[i][j]<-chromatin_force_cut_off)
                    {
                        ljmatrix[i][j]= -chromatin_force_cut_off;
                    }
                }
                if (chromatin_contant_matrix_calculation_flag==true) {
                    if( distance1< Dissimilaritycut )
                    {
                        Cm[i][j]=Cm[i][j]+1;
                        Cm1[i][j]=Cm1[i][j]+1;
                        Cm3[i][j]=Cm3[i][j]+1;
                        CmLaststeps[i][j]=CmLaststeps[i][j]+1;
                        
                    }
                    if(distance1< Dissimilaritycut)
                    {
                        contact_coordinates[0]=(Chromatin_Bead_Position[i][0]+Chromatin_Bead_Position[j][0])/2.0;
                        contact_coordinates[1]=(Chromatin_Bead_Position[i][1]+Chromatin_Bead_Position[j][1])/2.0;
                        contact_coordinates[2]=(Chromatin_Bead_Position[i][2]+Chromatin_Bead_Position[j][2])/2.0;
                        
                        contact_coordinates[0] += -COM_of_nucleus_membrane[0];
                        contact_coordinates[1] += -COM_of_nucleus_membrane[1];
                        contact_coordinates[2] += -COM_of_nucleus_membrane[2];
                        
                        ii=(int)floor( contact_coordinates[0]/mesh_size) + (int)floor(size_spatialcontactmatrix/2);
                        jj=(int)floor( contact_coordinates[1]/mesh_size) + (int)floor(size_spatialcontactmatrix/2);
                        kk=(int)floor( contact_coordinates[2]/mesh_size) + (int)floor(size_spatialcontactmatrix/2);
                        
                        if( abs(contact_coordinates[0]) < effective_thickness_of_spatial_dis )
                        {
                            spatial_contact_matrix_x[jj][kk]=spatial_contact_matrix_x[jj][kk]+1;
                        }
                        
                        if( abs(contact_coordinates[1]) < effective_thickness_of_spatial_dis )
                        {
                            spatial_contact_matrix_y[ii][kk]=spatial_contact_matrix_y[ii][kk]+1;
                        }
                        
                        if( abs(contact_coordinates[2]) < effective_thickness_of_spatial_dis )
                        {
                            spatial_contact_matrix_z[ii][jj]=spatial_contact_matrix_x[ii][jj]+1;
                        }
                        
                    }
                }
                
                //---------------------------------spatial contact
                
                //---------------------------------spatial contact
            }
            else if (chromatin_contant_matrix_calculation_flag==true )
            {
                Cm[i][j]=Cm[i][j]+1; //
                Cm1[i][j]=Cm1[i][j]+1;
                Cm3[i][j]=Cm3[i][j]+1;
                CmLaststeps[i][j]=CmLaststeps[i][j]+1;
                
            }
            
            
            
        }
        
        
    }
    
    
    
    
    
    for (int i=0 ; i<=Chromatin_num_of_Beads-1 ;i++)   // ading LJ force to  calculater forces that effect the i th  bead:
    {
        flj[0]=0.0;
        flj[1]=0.0;
        flj[2]=0.0;
        
        for(int j=0 ;j<=Chromatin_num_of_Beads-1 ;j++)
        {
            flj[0]=flj[0]+ljmatrix[i][j]*(Chromatin_Bead_Position[i][0]-Chromatin_Bead_Position[j][0])/distance1;
            flj[1]=flj[1]+ljmatrix[i][j]*(Chromatin_Bead_Position[i][1]-Chromatin_Bead_Position[j][1])/distance1;
            flj[2]=flj[2]+ljmatrix[i][j]*(Chromatin_Bead_Position[i][2]-Chromatin_Bead_Position[j][2])/distance1;
            
        }
        
        Chromatin_Bead_Force[i][0] += flj[0];
        Chromatin_Bead_Force[i][1] += flj[1];
        Chromatin_Bead_Force[i][2] += flj[2];
    }
    
    /// end of L-J(spring)
    
    
    
}


double potentialenergychromatin(double (&Chromatin_Bead_Position)[Chromatin_num_of_Beads][3])
{
    ///not completed yet!!
    double V,deltax,deltay,deltaz ,distance1;// defined below in "for loop" in detail
    
    V=0.0;
    
    /// Spring
    
    for(int i=0 ;i<Chromatin_num_of_Beads-1 ; i++ )  // all beads interaction whit the next one
    {
        deltax=0.0;
        deltay=0.0;
        deltaz=0.0;
        distance1=0.0;
        
        deltax=Chromatin_Bead_Position[i][0]-Chromatin_Bead_Position[i+1][0];// delta z  betwwn i and i+1 th beads
        deltay=Chromatin_Bead_Position[i][1]-Chromatin_Bead_Position[i+1][1];// delta z  betwwn i and i+1 th beads
        deltaz=Chromatin_Bead_Position[i][2]-Chromatin_Bead_Position[i+1][2];// delta z  betwwn i and i+1 th beads
        distance1=sqrt(deltaz*deltaz+deltay*deltay+deltax*deltax); // distance btween i th and i+1 th  bead
        
        V= V- (Chromatin_spring_coefficient*Rmaxchromatin*Rmaxchromatin/2.0)*log( 1.0-distance1*distance1/(Rmaxchromatin*Rmaxchromatin) );
        
    }
    
    
    double  vec1[3],vec2[3],vec1l,vec2l;
    
    
    for(int i=0 ;i<Chromatin_num_of_Beads-2 ; i++ )  // all beads interaction whit the next two beads
    {
        vec1[0]=Chromatin_Bead_Position[i+1][0]-Chromatin_Bead_Position[i][0];
        vec1[1]=Chromatin_Bead_Position[i+1][1]-Chromatin_Bead_Position[i][1];
        vec1[2]=Chromatin_Bead_Position[i+1][2]-Chromatin_Bead_Position[i][2];
        vec1l=vectorlength(vec1);
        
        vec2[0]=Chromatin_Bead_Position[i+2][0]-Chromatin_Bead_Position[i+1][0];
        vec2[1]=Chromatin_Bead_Position[i+2][1]-Chromatin_Bead_Position[i+1][1];
        vec2[2]=Chromatin_Bead_Position[i+2][2]-Chromatin_Bead_Position[i+1][2];
        vec2l=vectorlength(vec2);
        
        V=V+(Chromatin_bending_coefficient/2.0)*(1.0- innerproduct(vec1,vec2)/(vec1l*vec2l) );
    }
    
    return V;
}

void Chromatin_membrane_Barrier(double Chromatin_Bead_Position[Chromatin_num_of_Beads][3], double Chromatin_Bead_Velocity[Chromatin_num_of_Beads][3], double Chromatin_Bead_Force[Chromatin_num_of_Beads][3], double Membrane_Node_Position[][3], double Membrane_Node_Velocity[][3], int Membrane_triangle_list[Membrane_num_of_Triangles][3], int Membrane_Normal_direction[Membrane_num_of_Triangles][2], int Outer_Membrane_num_of_triangles)
{
    double membrane_triangle_COM_position[3]; // place of center of mass of triangle
    double membrane_triangle_COM_velocity[3]; // velocity of center of mass of triangle
    
    double chromatin_membrane_distance_amplitude;// distance between solvent particle and com of triangle
    double ABxAC[3],AB[3],AC[3]; // normal vector of membrane
    double chromatin_triangle_distance_vector[3];
    double perpendicular_distance;
    double relevant_velocity[3];
    
    for (int i =Outer_Membrane_num_of_triangles; i < Membrane_num_of_Triangles  ; i++)  // for each triangle:
    {
        membrane_triangle_COM_position[0]=(Membrane_Node_Position[ Membrane_triangle_list[i][0]][0] + Membrane_Node_Position[Membrane_triangle_list[i][1]][0] + Membrane_Node_Position[ Membrane_triangle_list[i][2]][0])/3.0;
        membrane_triangle_COM_position[1]=(Membrane_Node_Position[ Membrane_triangle_list[i][0]][1] + Membrane_Node_Position[Membrane_triangle_list[i][1]][1] + Membrane_Node_Position[ Membrane_triangle_list[i][2]][1])/3.0;
        membrane_triangle_COM_position[2]=(Membrane_Node_Position[ Membrane_triangle_list[i][0]][2] + Membrane_Node_Position[Membrane_triangle_list[i][1]][2] + Membrane_Node_Position[ Membrane_triangle_list[i][2]][2])/3.0;
        
        membrane_triangle_COM_velocity[0]=(Membrane_Node_Velocity[ Membrane_triangle_list[i][0]][0] + Membrane_Node_Velocity[Membrane_triangle_list[i][1]][0] + Membrane_Node_Velocity[ Membrane_triangle_list[i][2]][0])/3.0;
        membrane_triangle_COM_velocity[1]=(Membrane_Node_Velocity[ Membrane_triangle_list[i][0]][1] + Membrane_Node_Velocity[Membrane_triangle_list[i][1]][1] + Membrane_Node_Velocity[ Membrane_triangle_list[i][2]][1])/3.0;
        membrane_triangle_COM_velocity[2]=(Membrane_Node_Velocity[ Membrane_triangle_list[i][0]][2] + Membrane_Node_Velocity[Membrane_triangle_list[i][1]][2] + Membrane_Node_Velocity[ Membrane_triangle_list[i][2]][2])/3.0;
        
        for (int ch=0;ch<Chromatin_num_of_Beads;ch++)
        {
            // checking the interaction conditions:
            chromatin_membrane_distance_amplitude =sqrt( (membrane_triangle_COM_position[0] -Chromatin_Bead_Position[ch][0] ) * (membrane_triangle_COM_position[0] - Chromatin_Bead_Position[ch][0]) + (membrane_triangle_COM_position[1] - Chromatin_Bead_Position[ch][1]) * (membrane_triangle_COM_position[1]- Chromatin_Bead_Position[ch][1]) + (membrane_triangle_COM_position[2] - Chromatin_Bead_Position[ch][2]) * (membrane_triangle_COM_position[2] - Chromatin_Bead_Position[ch][2]));
            
            if (  chromatin_membrane_distance_amplitude < sqrt(0.43*Node_radius * 0.43*Node_radius + Chromatin_Membrane_Radius_of_Hard_Sphere_Interaction*Chromatin_Membrane_Radius_of_Hard_Sphere_Interaction )  )
            {
                AB[0]=Membrane_Node_Position[ Membrane_triangle_list[i][1]][0]-Membrane_Node_Position[ Membrane_triangle_list[i][0]][0];
                AB[1]=Membrane_Node_Position[ Membrane_triangle_list[i][1]][1]-Membrane_Node_Position[ Membrane_triangle_list[i][0]][1];
                AB[2]=Membrane_Node_Position[ Membrane_triangle_list[i][1]][2]-Membrane_Node_Position[ Membrane_triangle_list[i][0]][2];
                AC[0]=Membrane_Node_Position[ Membrane_triangle_list[i][2]][0]-Membrane_Node_Position[ Membrane_triangle_list[i][0]][0];
                AC[1]=Membrane_Node_Position[ Membrane_triangle_list[i][2]][1]-Membrane_Node_Position[ Membrane_triangle_list[i][0]][1];
                AC[2]=Membrane_Node_Position[ Membrane_triangle_list[i][2]][2]-Membrane_Node_Position[ Membrane_triangle_list[i][0]][2];
                crossvector(ABxAC,AB,AC);
                ABxAC[0]=ABxAC[0]*Membrane_Normal_direction[i][1];
                ABxAC[1]=ABxAC[1]*Membrane_Normal_direction[i][1];
                ABxAC[2]=ABxAC[2]*Membrane_Normal_direction[i][1];
                
                chromatin_triangle_distance_vector[0]=Chromatin_Bead_Position[ch][0]-membrane_triangle_COM_position[0];
                chromatin_triangle_distance_vector[1]=Chromatin_Bead_Position[ch][1]-membrane_triangle_COM_position[1];
                chromatin_triangle_distance_vector[2]=Chromatin_Bead_Position[ch][2]-membrane_triangle_COM_position[2];
                perpendicular_distance=innerproduct(chromatin_triangle_distance_vector,ABxAC)/vectorlength(ABxAC);
                
                if    (  abs( perpendicular_distance )<Chromatin_Membrane_Radius_of_Hard_Sphere_Interaction  )  // is solvent close enough to surface?
                {
                    
                    relevant_velocity[0] = Chromatin_Bead_Velocity[ch][0]-membrane_triangle_COM_velocity[0];
                    relevant_velocity[1] = Chromatin_Bead_Velocity[ch][1]-membrane_triangle_COM_velocity[1];
                    relevant_velocity[2] = Chromatin_Bead_Velocity[ch][2]-membrane_triangle_COM_velocity[2];
                    
                    if (   innerproduct(relevant_velocity,ABxAC)*perpendicular_distance<0    ) // is solvent moving toward triangle?
                    {
                        for (int i1 = 0; i1 < 3; i1++)
                        {
                            Membrane_Node_Velocity[Membrane_triangle_list[i][i1]][0] += (double) (2*Chromatin_Bead_Mass / (Chromatin_Bead_Mass + 3*Membrane_Node_Mass)) * relevant_velocity[0];
                            Membrane_Node_Velocity[Membrane_triangle_list[i][i1]][1] += (double) (2*Chromatin_Bead_Mass / (Chromatin_Bead_Mass + 3*Membrane_Node_Mass)) * relevant_velocity[1];
                            Membrane_Node_Velocity[Membrane_triangle_list[i][i1]][2] += (double) (2*Chromatin_Bead_Mass / (Chromatin_Bead_Mass + 3*Membrane_Node_Mass)) * relevant_velocity[2];
                        }
                        Chromatin_Bead_Velocity[ch][0] += - (double) (6*Membrane_Node_Mass / (Chromatin_Bead_Mass + 3*Membrane_Node_Mass)) * relevant_velocity[0];
                        Chromatin_Bead_Velocity[ch][1] += - (double) (6*Membrane_Node_Mass / (Chromatin_Bead_Mass + 3*Membrane_Node_Mass)) * relevant_velocity[1];
                        Chromatin_Bead_Velocity[ch][2] += - (double) (6*Membrane_Node_Mass / (Chromatin_Bead_Mass + 3*Membrane_Node_Mass)) * relevant_velocity[2];
                    }
                    
                }
            }
            
        }
        
    }
    
}
void Chromatin_membrane_Barrier_2(double Chromatin_Bead_Position[Chromatin_num_of_Beads][3], double Chromatin_Bead_Velocity[Chromatin_num_of_Beads][3], double Chromatin_Bead_Force[Chromatin_num_of_Beads][3], double Membrane_Node_Position[][3], double Membrane_Node_Velocity[][3], int Membrane_triangle_list[Membrane_num_of_Triangles][3], int Membrane_Normal_direction[Membrane_num_of_Triangles][2], int Outer_Membrane_num_of_triangles)
{
    double membrane_triangle_COM_position[3]; // place of center of mass of triangle
    double membrane_triangle_COM_velocity[3]; // velocity of center of mass of triangle
    double chromatin_membrane_distance_amplitude;// distance between solvent particle and com of triangle
    double ABxAC[3], AB[3], AC[3], ABxAC_unit_vector[3]; // normal vector of membrane
    double chromatin_triangle_distance_vector[3];
    double perpendicular_distance;
    double relevant_velocity[3];
    
    for (int i =Outer_Membrane_num_of_triangles; i < Membrane_num_of_Triangles  ; i++)  // for each triangle:
    {
        membrane_triangle_COM_position[0]=(Membrane_Node_Position[ Membrane_triangle_list[i][0]][0] + Membrane_Node_Position[Membrane_triangle_list[i][1]][0] + Membrane_Node_Position[ Membrane_triangle_list[i][2]][0])/3.0;
        membrane_triangle_COM_position[1]=(Membrane_Node_Position[ Membrane_triangle_list[i][0]][1] + Membrane_Node_Position[Membrane_triangle_list[i][1]][1] + Membrane_Node_Position[ Membrane_triangle_list[i][2]][1])/3.0;
        membrane_triangle_COM_position[2]=(Membrane_Node_Position[ Membrane_triangle_list[i][0]][2] + Membrane_Node_Position[Membrane_triangle_list[i][1]][2] + Membrane_Node_Position[ Membrane_triangle_list[i][2]][2])/3.0;
        
        membrane_triangle_COM_velocity[0]=(Membrane_Node_Velocity[ Membrane_triangle_list[i][0]][0] + Membrane_Node_Velocity[Membrane_triangle_list[i][1]][0] + Membrane_Node_Velocity[ Membrane_triangle_list[i][2]][0])/3.0;
        membrane_triangle_COM_velocity[1]=(Membrane_Node_Velocity[ Membrane_triangle_list[i][0]][1] + Membrane_Node_Velocity[Membrane_triangle_list[i][1]][1] + Membrane_Node_Velocity[ Membrane_triangle_list[i][2]][1])/3.0;
        membrane_triangle_COM_velocity[2]=(Membrane_Node_Velocity[ Membrane_triangle_list[i][0]][2] + Membrane_Node_Velocity[Membrane_triangle_list[i][1]][2] + Membrane_Node_Velocity[ Membrane_triangle_list[i][2]][2])/3.0;
        
        for (int ch=0;ch<Chromatin_num_of_Beads;ch++)
        {
            chromatin_membrane_distance_amplitude =sqrt( (membrane_triangle_COM_position[0] -Chromatin_Bead_Position[ch][0] ) * (membrane_triangle_COM_position[0] - Chromatin_Bead_Position[ch][0]) + (membrane_triangle_COM_position[1] - Chromatin_Bead_Position[ch][1]) * (membrane_triangle_COM_position[1]- Chromatin_Bead_Position[ch][1]) + (membrane_triangle_COM_position[2] - Chromatin_Bead_Position[ch][2]) * (membrane_triangle_COM_position[2] - Chromatin_Bead_Position[ch][2]));
            
            if (  chromatin_membrane_distance_amplitude < sqrt(0.43*Node_radius * 0.43*Node_radius + Chromatin_Membrane_Radius_of_Hard_Sphere_Interaction*Chromatin_Membrane_Radius_of_Hard_Sphere_Interaction )  )
            {
                AB[0]=Membrane_Node_Position[ Membrane_triangle_list[i][1]][0]-Membrane_Node_Position[ Membrane_triangle_list[i][0]][0];
                AB[1]=Membrane_Node_Position[ Membrane_triangle_list[i][1]][1]-Membrane_Node_Position[ Membrane_triangle_list[i][0]][1];
                AB[2]=Membrane_Node_Position[ Membrane_triangle_list[i][1]][2]-Membrane_Node_Position[ Membrane_triangle_list[i][0]][2];
                AC[0]=Membrane_Node_Position[ Membrane_triangle_list[i][2]][0]-Membrane_Node_Position[ Membrane_triangle_list[i][0]][0];
                AC[1]=Membrane_Node_Position[ Membrane_triangle_list[i][2]][1]-Membrane_Node_Position[ Membrane_triangle_list[i][0]][1];
                AC[2]=Membrane_Node_Position[ Membrane_triangle_list[i][2]][2]-Membrane_Node_Position[ Membrane_triangle_list[i][0]][2];
                crossvector(ABxAC,AB,AC);
                
                ABxAC[0]=ABxAC[0]*Membrane_Normal_direction[i][1];
                ABxAC[1]=ABxAC[1]*Membrane_Normal_direction[i][1];
                ABxAC[2]=ABxAC[2]*Membrane_Normal_direction[i][1];
                
                ABxAC_unit_vector[0]=ABxAC[0]/vectorlength(ABxAC);
                ABxAC_unit_vector[1]=ABxAC[1]/vectorlength(ABxAC);
                ABxAC_unit_vector[2]=ABxAC[2]/vectorlength(ABxAC);
                
                chromatin_triangle_distance_vector[0]=Chromatin_Bead_Position[ch][0]-membrane_triangle_COM_position[0];
                chromatin_triangle_distance_vector[1]=Chromatin_Bead_Position[ch][1]-membrane_triangle_COM_position[1];
                chromatin_triangle_distance_vector[2]=Chromatin_Bead_Position[ch][2]-membrane_triangle_COM_position[2];
                perpendicular_distance=innerproduct(chromatin_triangle_distance_vector,ABxAC_unit_vector);
                
                relevant_velocity[0] = Chromatin_Bead_Velocity[ch][0]-membrane_triangle_COM_velocity[0];
                relevant_velocity[1] = Chromatin_Bead_Velocity[ch][1]-membrane_triangle_COM_velocity[1];
                relevant_velocity[2] = Chromatin_Bead_Velocity[ch][2]-membrane_triangle_COM_velocity[2];
                
                if    (  (abs( perpendicular_distance )<Chromatin_Membrane_Radius_of_Hard_Sphere_Interaction) && (innerproduct(relevant_velocity, ABxAC)*perpendicular_distance<0) )
                {
                    double Chromatin_velocity_N_new, Chromatin_velocity_N, Membrane_triangle_COM_velocity_N_new, Membrane_triangle_COM_velocity_N;
                    
                    Chromatin_velocity_N=innerproduct(Chromatin_Bead_Velocity[ch], ABxAC_unit_vector);
                    Membrane_triangle_COM_velocity_N=innerproduct(membrane_triangle_COM_velocity, ABxAC_unit_vector);
                    
                    Chromatin_velocity_N_new=(Chromatin_velocity_N*(Chromatin_Bead_Mass-3*Membrane_Node_Mass)+2.0*3*Membrane_Node_Mass*Membrane_triangle_COM_velocity_N)/(Chromatin_Bead_Mass+3*Membrane_Node_Mass);
                    Membrane_triangle_COM_velocity_N_new=(Membrane_triangle_COM_velocity_N*(3*Membrane_Node_Mass-Actin_Node_Mass)+2.0*Chromatin_Bead_Mass*Chromatin_velocity_N)/(Chromatin_Bead_Mass+3*Membrane_Node_Mass);
                    
                    Chromatin_Bead_Velocity[ch][0]+= (-Chromatin_velocity_N + Chromatin_velocity_N_new)*ABxAC_unit_vector[0];
                    Chromatin_Bead_Velocity[ch][1]+= (-Chromatin_velocity_N + Chromatin_velocity_N_new)*ABxAC_unit_vector[1];
                    Chromatin_Bead_Velocity[ch][2]+= (-Chromatin_velocity_N + Chromatin_velocity_N_new)*ABxAC_unit_vector[2];
                    
                    for (int i1 = 0; i1 < 3; i1++)
                    {
                        Membrane_Node_Velocity[Membrane_triangle_list[i][i1]][0] += (-Membrane_triangle_COM_velocity_N + Membrane_triangle_COM_velocity_N_new)*ABxAC_unit_vector[0];
                        Membrane_Node_Velocity[Membrane_triangle_list[i][i1]][1] += (-Membrane_triangle_COM_velocity_N + Membrane_triangle_COM_velocity_N_new)*ABxAC_unit_vector[1];
                        Membrane_Node_Velocity[Membrane_triangle_list[i][i1]][2] += (-Membrane_triangle_COM_velocity_N + Membrane_triangle_COM_velocity_N_new)*ABxAC_unit_vector[2];
                    }
                } //END OF: if    (  (abs( perpendicular_distance )<Chromatin_Membrane_Radius_of_Hard_Sphere_Interaction) &&
            }//END OF: if (  chromatin_membrane_distance_amplitude < sqrt(0.43*a * 0.43*a +
            
        }//END OF: for (int ch=0;ch<Chromatin_num_of_Beads;ch++)
        
    } //END OF: for (int i =Outer_Membrane_num_of_triangles; i < Membrane_num_of_Triangles  ; i++)
    
}

void update_spatial_contact_matrix(double contact_coordinates[3])
{
    
    double mesh_size= 10.0* double (Nucleus_Membrane_radius) / double (size_spatialcontactmatrix) ;
    int ii,jj,kk;
    
    contact_coordinates[0] += -COM_of_nucleus_membrane[0];
    contact_coordinates[1] += -COM_of_nucleus_membrane[1];
    contact_coordinates[2] += -COM_of_nucleus_membrane[2];
    
    ii=(int)floor( contact_coordinates[0]/mesh_size) + (int)floor(size_spatialcontactmatrix/2);
    jj=(int)floor( contact_coordinates[1]/mesh_size) + (int)floor(size_spatialcontactmatrix/2);
    kk=(int)floor( contact_coordinates[2]/mesh_size) + (int)floor(size_spatialcontactmatrix/2);
    
    if( abs(contact_coordinates[0]) < effective_thickness_of_spatial_dis )
    {
        spatial_contact_matrix_x[jj][kk]=spatial_contact_matrix_x[jj][kk]+1;
    }
    
    if( abs(contact_coordinates[1]) < effective_thickness_of_spatial_dis )
    {
        spatial_contact_matrix_y[ii][kk]=spatial_contact_matrix_y[ii][kk]+1;
    }
    
    if( abs(contact_coordinates[2]) < effective_thickness_of_spatial_dis )
    {
        spatial_contact_matrix_z[ii][jj]=spatial_contact_matrix_x[ii][jj]+1;
    }
}

void COM_of_nucleus_membrane_function(double Membrane_Node_Position[][3],int Nucleus_Membrane_list_of_Nodes[], int  Outer_Membrane_num_of_Nodes, int Membrane_num_of_Nodes)
{
    int lambda;
    COM_of_nucleus_membrane[0]=0.0;
    COM_of_nucleus_membrane[1]=0.0;
    COM_of_nucleus_membrane[2]=0.0;
    
    for(int j=0; j<Membrane_num_of_Nodes- Outer_Membrane_num_of_Nodes;j++) // saving trajectory
    {   lambda=Nucleus_Membrane_list_of_Nodes[j];
        
        COM_of_nucleus_membrane[0] += Membrane_Node_Position[lambda][0];
        COM_of_nucleus_membrane[1] += Membrane_Node_Position[lambda][1];
        COM_of_nucleus_membrane[2] += Membrane_Node_Position[lambda][2];
        
    }
    
    COM_of_nucleus_membrane[0]=  COM_of_nucleus_membrane[0]/( double(Membrane_num_of_Nodes- Outer_Membrane_num_of_Nodes) );
    COM_of_nucleus_membrane[1]=  COM_of_nucleus_membrane[1]/( double(Membrane_num_of_Nodes- Outer_Membrane_num_of_Nodes) );
    COM_of_nucleus_membrane[2]=  COM_of_nucleus_membrane[2]/( double(Membrane_num_of_Nodes- Outer_Membrane_num_of_Nodes) );
    
}


void loop_force_chromatin(double Chromatin_Bead_Position[Chromatin_num_of_Beads][3],double Chromatin_Bead_Force[Chromatin_num_of_Beads][3], double &Total_Potential_Energy)

{
    
    double deltax,deltay,deltaz ,distance1,force;// defined below in "for loop" in detail
    /// loop
    int  loop_lenth =  floor((Chromatin_num_of_Beads/Chromatin_num_of_chains) / number_of_loops_in_each_chain ) - 1;
    int i;
    int i_2;
    
    for (int nchain=0;nchain<Chromatin_num_of_chains;nchain++)
    {
        
        for (int k=1 ; k<= number_of_loops_in_each_chain ; k++)
        {
            
            i=nchain*(Chromatin_num_of_Beads/Chromatin_num_of_chains) + loop_lenth*(k-1) ;
            i_2=i+loop_lenth;
            // cout<<i<<"  "<<i_2<<endl;
            
            // cout<<nchain<<"  "<<i<<endl;
            deltay=0.0;
            deltaz=0.0;
            distance1=0.0;
            force=0.0;
            
            deltax=Chromatin_Bead_Position[i][0]-Chromatin_Bead_Position[i_2][0];// delta z  betwwn i and i+1 th beads
            deltay=Chromatin_Bead_Position[i][1]-Chromatin_Bead_Position[i_2][1];// delta z  betwwn i and i+1 th beads
            deltaz=Chromatin_Bead_Position[i][2]-Chromatin_Bead_Position[i_2][2];// delta z  betwwn i and i+1 th beads
            distance1=sqrt(deltaz*deltaz+deltay*deltay+deltax*deltax); // distance btween i th and i+1 th  bead
            
            
            // force=+(Chromatin_spring_coefficient*Rmaxchromatin*Rmaxchromatin/2.0)*(2.0*distance1/(Rmaxchromatin*Rmaxchromatin))/(1.0-distance1*distance1/(Rmaxchromatin*Rmaxchromatin));
            
            force=Chromatin_spring_coefficient*(distance1- 0.70*sigmachromatin*Chromatin_Scaling_Factor);
            if (energy_calculation_flag==1.0) {
                Total_Potential_Energy+= 0.5*Chromatin_spring_coefficient*(distance1- 0.70*sigmachromatin*Chromatin_Scaling_Factor)*(distance1- 0.70*sigmachromatin*Chromatin_Scaling_Factor);
            }
            
            
            Chromatin_Bead_Force[i][0] += force*deltax/distance1 ; // force of springs and dashes
            Chromatin_Bead_Force[i][1] += force*deltay/distance1  ;// force of springs and dashes
            Chromatin_Bead_Force[i][2] += force*deltaz/distance1 ;// force of springs and dashes
            
            Chromatin_Bead_Force[i_2][0] += -force *  deltax/distance1  ;
            Chromatin_Bead_Force[i_2][1] += -force *  deltay/distance1  ;
            Chromatin_Bead_Force[i_2][2] += -force *  deltaz/distance1  ;
            
            
            
        }
        
        
    }
    
}


//_______________________________ECM functions


void  ECM_rigidity_constructor(double  ECM_Node_Position [][3],double ECM_Node_Pair_List[][3], double ECM_varying_stiffness_coefficient[], int ECM_num_of_Bonds)
{
    int node1,node2;
    
    for(int i=0 ;i<ECM_num_of_Bonds ; i++ )  // all beads interaction whit the next one
    {
        node1=(int) ECM_Node_Pair_List[i][0];
        node2=(int) ECM_Node_Pair_List[i][1];
        
        ECM_varying_stiffness_coefficient[i]=((ECM_Node_Position[node2][0]+ECM_Node_Position[node1][0])/2.0)*(2.0*ECM_Max_Gradient_Coefficient/ECM_Gradient_length)+(ECM_Min_Gradient_Coefficient+ECM_Max_Gradient_Coefficient);
        //        if (((ECM_Node_Position[node2][0]+ECM_Node_Position[node1][0])/2.0) < 0.0) {
        //            ECM_varying_stiffness_coefficient[i]=ECM_Min_Gradient_Coefficient;
        //        } else {
        //            ECM_varying_stiffness_coefficient[i]=ECM_Max_Gradient_Coefficient;
        //        }
        
        //Spining gradient-beginning
        //        ECM_varying_stiffness_coefficient[i]=((ECM_Max_Gradient_Coefficient-ECM_Min_Gradient_Coefficient)/(2.0*pi))*(atan((ECM_Node_Position[node2][2]+ECM_Node_Position[node1][2])/(ECM_Node_Position[node2][0]+ECM_Node_Position[node1][0])));
        //Spining gradient-end
    }
    
    
}



void  ECM_Force( double ECM_Node_Position[][3], double ECM_Node_Velocity[][3], double ECM_Node_Force[][3], double ECM_Node_Pair_List[][3], double ECM_upper_surface_Node_Pairs[], double ECM_varying_stiffness_coefficient[], int ECM_num_of_Bonds, double &Total_Potential_Energy) {
    
    double deltax, deltay, deltaz, temp_distance, initial_distance; // defined below in "for loop" in detail
    int temp_node_1,temp_node_2;
    double  fi[3];
    
    for(int i=0 ;i<ECM_num_of_Bonds ; i++ )  // Pair force, Kelvin spring (The model is very breifly explaind in the 'if' below).
    {
        temp_node_1=(int) ECM_Node_Pair_List[i][0];
        temp_node_2=(int) ECM_Node_Pair_List[i][1];
        initial_distance=ECM_Node_Pair_List[i][2];
        
        deltax=ECM_Node_Position[temp_node_2][0]-ECM_Node_Position[temp_node_1][0];
        deltay=ECM_Node_Position[temp_node_2][1]-ECM_Node_Position[temp_node_1][1];
        deltaz=ECM_Node_Position[temp_node_2][2]-ECM_Node_Position[temp_node_1][2];
        temp_distance=sqrt(deltax*deltax+deltay*deltay+deltaz*deltaz);
        
        fi[0] = -ECM_varying_stiffness_coefficient[i]* ( (temp_distance-initial_distance)/temp_distance ) *  deltax;
        fi[1] = -ECM_varying_stiffness_coefficient[i]* ( (temp_distance-initial_distance)/temp_distance ) *  deltay;
        fi[2] = -ECM_varying_stiffness_coefficient[i]* ( (temp_distance-initial_distance)/temp_distance ) *  deltaz;
        if (energy_calculation_flag==1.0) {
            Total_Potential_Energy+=0.5*ECM_varying_stiffness_coefficient[i]*(temp_distance-initial_distance)*(temp_distance-initial_distance);
        }
        
        
        ECM_Node_Force[temp_node_1][0] += fi[0]  ;
        ECM_Node_Force[temp_node_1][1] += fi[1]  ;
        ECM_Node_Force[temp_node_1][2] += fi[2]  ;
        
        ECM_Node_Force[temp_node_2][0] += - fi[0]  ;
        ECM_Node_Force[temp_node_2][1] += - fi[1]  ;
        ECM_Node_Force[temp_node_2][2] += - fi[2]  ;
        
        
        if( (initial_distance <2.5 && initial_distance>1.5) || (initial_distance >2.5 && temp_distance<initial_distance  )  || (initial_distance <1.5  && temp_distance>initial_distance  ) )
        {
            ECM_Node_Pair_List[i][2] += MD_Time_Step*( (temp_distance-initial_distance)/temp_distance )/ECM_kelvin_damping_coefficient ;
            //We have used the Kelvin model for the Visco elasticity. Kelvin model: A spring and dashpot (in series) are in parralel with a spring.
        }
        
        fi[0]= miuDampECM*(ECM_Node_Velocity[temp_node_1][0]-ECM_Node_Velocity[temp_node_2][0]);
        fi[1]= miuDampECM*(ECM_Node_Velocity[temp_node_1][1]-ECM_Node_Velocity[temp_node_2][1]);
        fi[2]= miuDampECM*(ECM_Node_Velocity[temp_node_1][2]-ECM_Node_Velocity[temp_node_2][2]);
        
        if (energy_calculation_flag==1.0) {
            Total_Potential_Energy-=0.5*miuDampECM*((ECM_Node_Velocity[temp_node_1][0]-ECM_Node_Velocity[temp_node_2][0])*(ECM_Node_Velocity[temp_node_1][0]-ECM_Node_Velocity[temp_node_2][0])+(ECM_Node_Velocity[temp_node_1][0]-ECM_Node_Velocity[temp_node_2][1])*(ECM_Node_Velocity[temp_node_1][1]-ECM_Node_Velocity[temp_node_2][1])*(ECM_Node_Velocity[temp_node_1][1]-ECM_Node_Velocity[temp_node_2][1])+(ECM_Node_Velocity[temp_node_1][2]-ECM_Node_Velocity[temp_node_2][2])*(ECM_Node_Velocity[temp_node_1][2]-ECM_Node_Velocity[temp_node_2][2]));
        }
        
        
        ECM_Node_Force[temp_node_1][0] += fi[0]  ;
        ECM_Node_Force[temp_node_1][1] += fi[1]  ;
        ECM_Node_Force[temp_node_1][2] += fi[2]  ;
        
        ECM_Node_Force[temp_node_2][0] += - fi[0]  ;
        ECM_Node_Force[temp_node_2][1] += - fi[1]  ;
        ECM_Node_Force[temp_node_2][2] += - fi[2]  ;
    }
    
}



void Membrane_ECM_interaction(int istep,double Membrane_Node_Position [][3],double Membrane_Node_Velocity [][3],double Membrane_Node_Force [][3], double  ECM_Node_Position [][3],double  ECM_Node_Velocity [][3],double  ECM_Node_Force [][3],int ECM_surface_triangle_list[ECM_Surface_num_of_Triangles][3],int Membrane_triangle_list[Membrane_num_of_Triangles][3],int Membrane_Normal_direction[Membrane_num_of_Triangles][2],  int ECMandmembrane[][ 1+maximumHowmanyECMCouldInteractWithmembrane] ,double Cellcom[3], int Outer_Membrane_num_of_triangles, int Membrane_num_of_Nodes)
{
    
    
    double ECM_triangle_AB[3], ECM_triangle_AC[3], ECM_triangle_ABxAC_unit_vector[3], ECM_triangle_ABxAC_length;
    double Membrane_triangle_AB[3], Membrane_triangle_AC[3], Membrane_triangle_ABxAC[3]; // normal vector of membrane
    int temp_ECM_node_A, temp_ECM_node_B, temp_ECM_node_C;
    
    double D[3], Membrane_triangle_COM_ECM_distance_vec[3], membrane_tri_COM_ECM_distance_normal_component;
    double a1,a2,a3,a1a2a3,r2,r6;
    double t1,flength,f[3],PD[3];
    double normal2[3], Membrane_triangle_COM[3]; // normal vector of membrane
    
    int index;
    double attachingforce = 0.0;
    int counter1=0;
    //F-net-com=0
    double Fnet__Active[3];  Fnet__Active[0]=0.0; Fnet__Active[1]=0.0; Fnet__Active[2]=0.0;
    double Fnet_Migration_force[3];  Fnet_Migration_force[0]=0.0; Fnet_Migration_force[1]=0.0; Fnet_Migration_force[2]=0.0;
    int Fnet__ActiveList[Outer_Membrane_num_of_triangles],numberFnet__ActiveList=0;
    double fx,fy,fz;
    int directionIsInward,checkmark;
    double comTriVec[3],lvec;
    
    double Migrating___force_amp;
    double Migrating___force_Vec[3][Membrane_num_of_Nodes];
    for (int i=0; i<Membrane_num_of_Nodes; i++) {
        Migrating___force_Vec[0][i]=0;
        Migrating___force_Vec[1][i]=0;
        Migrating___force_Vec[2][i]=0;
    }
    double Total_generated_force=0;
    double alpha_coeff;
    
    
    if(istep % Membrane_ECM_interaction_update_step == 0 )///__new interaction strategy :
    {
        for(int i=0;i<Outer_Membrane_num_of_triangles;i++)
        {
            counter1=0;
            ECMandmembrane[i][0]=0;
            for ( int j=0;j<ECM_Surface_num_of_Triangles;j++)
            {
                Membrane_triangle_AB[0]=Membrane_Node_Position[ Membrane_triangle_list[i][1]][0]-Membrane_Node_Position[ Membrane_triangle_list[i][0]][0];
                Membrane_triangle_AB[1]=Membrane_Node_Position[ Membrane_triangle_list[i][1]][1]-Membrane_Node_Position[ Membrane_triangle_list[i][0]][1];
                Membrane_triangle_AB[2]=Membrane_Node_Position[ Membrane_triangle_list[i][1]][2]-Membrane_Node_Position[ Membrane_triangle_list[i][0]][2];
                Membrane_triangle_AC[0]=Membrane_Node_Position[ Membrane_triangle_list[i][2]][0]-Membrane_Node_Position[ Membrane_triangle_list[i][0]][0];
                Membrane_triangle_AC[1]=Membrane_Node_Position[ Membrane_triangle_list[i][2]][1]-Membrane_Node_Position[ Membrane_triangle_list[i][0]][1];
                Membrane_triangle_AC[2]=Membrane_Node_Position[ Membrane_triangle_list[i][2]][2]-Membrane_Node_Position[ Membrane_triangle_list[i][0]][2];
                crossvector(Membrane_triangle_ABxAC, Membrane_triangle_AB, Membrane_triangle_AC);
                Membrane_triangle_ABxAC[0]*=Membrane_Normal_direction[i][1];
                Membrane_triangle_ABxAC[1]*=Membrane_Normal_direction[i][1];
                Membrane_triangle_ABxAC[2]*=Membrane_Normal_direction[i][1];
                
                temp_ECM_node_A=ECM_surface_triangle_list[j][0];
                temp_ECM_node_B=ECM_surface_triangle_list[j][1];
                temp_ECM_node_C=ECM_surface_triangle_list[j][2];
                
                ECM_triangle_AB[0]=ECM_Node_Position[temp_ECM_node_B][0]-ECM_Node_Position[temp_ECM_node_A][0];
                ECM_triangle_AB[1]=ECM_Node_Position[temp_ECM_node_B][1]-ECM_Node_Position[temp_ECM_node_A][1];
                ECM_triangle_AB[2]=ECM_Node_Position[temp_ECM_node_B][2]-ECM_Node_Position[temp_ECM_node_A][2];
                ECM_triangle_AC[0]=ECM_Node_Position[temp_ECM_node_C][0]-ECM_Node_Position[temp_ECM_node_A][0];
                ECM_triangle_AC[1]=ECM_Node_Position[temp_ECM_node_C][1]-ECM_Node_Position[temp_ECM_node_A][1];
                ECM_triangle_AC[2]=ECM_Node_Position[temp_ECM_node_C][2]-ECM_Node_Position[temp_ECM_node_A][2];
                crossvector(ECM_triangle_ABxAC_unit_vector, ECM_triangle_AB, ECM_triangle_AC);
                ECM_triangle_ABxAC_length=vectorlength(ECM_triangle_ABxAC_unit_vector);
                ECM_triangle_ABxAC_unit_vector[0]/=ECM_triangle_ABxAC_length;
                ECM_triangle_ABxAC_unit_vector[1]/=ECM_triangle_ABxAC_length;
                ECM_triangle_ABxAC_unit_vector[2]/=ECM_triangle_ABxAC_length;  // it is the norml of the ECM element
                
                if(ECMrepultionforceToTrianglesOrnode==0) //0::triangles 1::nodes
                {
                    //com of membrane element
                    Membrane_triangle_COM[0]= ( Membrane_Node_Position[ Membrane_triangle_list[i][0]][0] + Membrane_Node_Position[ Membrane_triangle_list[i][1]][0] +Membrane_Node_Position[ Membrane_triangle_list[i][2]][0] ) /3.0;
                    Membrane_triangle_COM[1]= ( Membrane_Node_Position[ Membrane_triangle_list[i][0]][1] +Membrane_Node_Position[ Membrane_triangle_list[i][1]][1] +Membrane_Node_Position[ Membrane_triangle_list[i][2]][1] ) /3.0;
                    Membrane_triangle_COM[2]= ( Membrane_Node_Position[ Membrane_triangle_list[i][0]][2] +Membrane_Node_Position[ Membrane_triangle_list[i][1]][2] +Membrane_Node_Position[ Membrane_triangle_list[i][2]][2] ) /3.0;
                    
                    
                    // check the condition of  ECM and COM of membrane beig close enough to ech othe
                    Membrane_triangle_COM_ECM_distance_vec[0]= Membrane_triangle_COM[0] -ECM_Node_Position[temp_ECM_node_A][0];
                    Membrane_triangle_COM_ECM_distance_vec[1]= Membrane_triangle_COM[1] -ECM_Node_Position[temp_ECM_node_A][1];
                    Membrane_triangle_COM_ECM_distance_vec[2]= Membrane_triangle_COM[2] -ECM_Node_Position[temp_ECM_node_A][2];
                    membrane_tri_COM_ECM_distance_normal_component=innerproduct(Membrane_triangle_COM_ECM_distance_vec,ECM_triangle_ABxAC_unit_vector);
                    
                    //if they  are close are they in good position?:
                    if( abs(membrane_tri_COM_ECM_distance_normal_component)<  considerTimesHowFarElemets*epsilonECM + epsilonECM*2.0  )
                    {
                        D[0]=Membrane_triangle_COM[0] - membrane_tri_COM_ECM_distance_normal_component * ECM_triangle_ABxAC_unit_vector[0];
                        D[1]=Membrane_triangle_COM[1] - membrane_tri_COM_ECM_distance_normal_component * ECM_triangle_ABxAC_unit_vector[1];
                        D[2]=Membrane_triangle_COM[2] - membrane_tri_COM_ECM_distance_normal_component * ECM_triangle_ABxAC_unit_vector[2];
                        
                        Membrane_ECM_interactionHELPER_newstrategy(ECM_Node_Position, temp_ECM_node_A, temp_ECM_node_B, temp_ECM_node_C, ECM_triangle_ABxAC_unit_vector,D,a1,a2,a3);
                        if(a1>=0 & a2>=0 & a3>=0 )
                        {
                            counter1=counter1+1;
                            ECMandmembrane[i][0]=counter1; // says how many ECM element are interacting with tish elemet
                            ECMandmembrane[i][counter1]=j;  //says what element is interacting !
                            
                            if(counter1 > maximumHowmanyECMCouldInteractWithmembrane )
                            {
                                cout << " error:  maximumHowmanyECMCouldInteractWithmembrane   is NOT enough,increase it! ___the current number is "<< counter1 <<endl;
                            }
                            
                        }
                    }
                }
                
                
                
                
                
                
                
                
                else if(ECMrepultionforceToTrianglesOrnode==1)
                {
                    checkmark=0;
                    
                    
                    //tri 0
                    Membrane_triangle_COM_ECM_distance_vec[0]= Membrane_Node_Position[ Membrane_triangle_list[i][0]][0] -ECM_Node_Position[temp_ECM_node_A][0];
                    Membrane_triangle_COM_ECM_distance_vec[1]= Membrane_Node_Position[ Membrane_triangle_list[i][0]][1] -ECM_Node_Position[temp_ECM_node_A][1];
                    Membrane_triangle_COM_ECM_distance_vec[2]= Membrane_Node_Position[ Membrane_triangle_list[i][0]][2] -ECM_Node_Position[temp_ECM_node_A][2];
                    membrane_tri_COM_ECM_distance_normal_component=innerproduct(Membrane_triangle_COM_ECM_distance_vec, ECM_triangle_ABxAC_unit_vector);
                    if( abs(membrane_tri_COM_ECM_distance_normal_component)<  considerTimesHowFarElemets*epsilonECM + epsilonECM*2.0  )
                    {
                        D[0]=Membrane_Node_Position[ Membrane_triangle_list[i][0]][0] - membrane_tri_COM_ECM_distance_normal_component * ECM_triangle_ABxAC_unit_vector[0];
                        D[1]=Membrane_Node_Position[ Membrane_triangle_list[i][0]][1] - membrane_tri_COM_ECM_distance_normal_component * ECM_triangle_ABxAC_unit_vector[1];
                        D[2]=Membrane_Node_Position[ Membrane_triangle_list[i][0]][2] - membrane_tri_COM_ECM_distance_normal_component * ECM_triangle_ABxAC_unit_vector[2];
                        Membrane_ECM_interactionHELPER_newstrategy(ECM_Node_Position, temp_ECM_node_A, temp_ECM_node_B, temp_ECM_node_C, ECM_triangle_ABxAC_unit_vector,D,a1,a2,a3);
                        if(a1>=0 & a2>=0 & a3>=0 )
                        {
                            checkmark=1;
                        }
                    }
                    
                    //tri 1
                    Membrane_triangle_COM_ECM_distance_vec[0]= Membrane_Node_Position[ Membrane_triangle_list[i][1]][0] -ECM_Node_Position[temp_ECM_node_A][0];
                    Membrane_triangle_COM_ECM_distance_vec[1]= Membrane_Node_Position[ Membrane_triangle_list[i][1]][1] -ECM_Node_Position[temp_ECM_node_A][1];
                    Membrane_triangle_COM_ECM_distance_vec[2]= Membrane_Node_Position[ Membrane_triangle_list[i][1]][2] -ECM_Node_Position[temp_ECM_node_A][2];
                    membrane_tri_COM_ECM_distance_normal_component=innerproduct(Membrane_triangle_COM_ECM_distance_vec, ECM_triangle_ABxAC_unit_vector);
                    if( abs(membrane_tri_COM_ECM_distance_normal_component)<  considerTimesHowFarElemets*epsilonECM + epsilonECM*2.0  )
                    {
                        D[0]=Membrane_Node_Position[ Membrane_triangle_list[i][1]][0] - membrane_tri_COM_ECM_distance_normal_component * ECM_triangle_ABxAC_unit_vector[0];
                        D[1]=Membrane_Node_Position[ Membrane_triangle_list[i][1]][1] - membrane_tri_COM_ECM_distance_normal_component * ECM_triangle_ABxAC_unit_vector[1];
                        D[2]=Membrane_Node_Position[ Membrane_triangle_list[i][1]][2] - membrane_tri_COM_ECM_distance_normal_component * ECM_triangle_ABxAC_unit_vector[2];
                        Membrane_ECM_interactionHELPER_newstrategy(ECM_Node_Position, temp_ECM_node_A, temp_ECM_node_B, temp_ECM_node_C, ECM_triangle_ABxAC_unit_vector,D,a1,a2,a3);
                        if(a1>=0 & a2>=0 & a3>=0 )
                        {
                            checkmark=1;
                        }
                    }
                    
                    
                    
                    //tri 2
                    Membrane_triangle_COM_ECM_distance_vec[0]= Membrane_Node_Position[ Membrane_triangle_list[i][2]][0] -ECM_Node_Position[temp_ECM_node_A][0];
                    Membrane_triangle_COM_ECM_distance_vec[1]= Membrane_Node_Position[ Membrane_triangle_list[i][2]][1] -ECM_Node_Position[temp_ECM_node_A][1];
                    Membrane_triangle_COM_ECM_distance_vec[2]= Membrane_Node_Position[ Membrane_triangle_list[i][2]][2] -ECM_Node_Position[temp_ECM_node_A][2];
                    membrane_tri_COM_ECM_distance_normal_component=innerproduct(Membrane_triangle_COM_ECM_distance_vec, ECM_triangle_ABxAC_unit_vector);
                    if( abs(membrane_tri_COM_ECM_distance_normal_component)<  considerTimesHowFarElemets*epsilonECM + epsilonECM*2.0  )
                    {
                        D[0]=Membrane_Node_Position[ Membrane_triangle_list[i][2]][0] - membrane_tri_COM_ECM_distance_normal_component * ECM_triangle_ABxAC_unit_vector[0];
                        D[1]=Membrane_Node_Position[ Membrane_triangle_list[i][2]][1] - membrane_tri_COM_ECM_distance_normal_component * ECM_triangle_ABxAC_unit_vector[1];
                        D[2]=Membrane_Node_Position[ Membrane_triangle_list[i][2]][2] - membrane_tri_COM_ECM_distance_normal_component * ECM_triangle_ABxAC_unit_vector[2];
                        Membrane_ECM_interactionHELPER_newstrategy(ECM_Node_Position, temp_ECM_node_A, temp_ECM_node_B, temp_ECM_node_C, ECM_triangle_ABxAC_unit_vector,D,a1,a2,a3);
                        if(a1>=0 & a2>=0 & a3>=0 )
                        {
                            checkmark=1;
                        }
                    }
                    
                    
                    
                    if(checkmark==1)
                    {
                        counter1=counter1+1;
                        ECMandmembrane[i][0]=counter1; // says how many ECM element are interacting with tish elemet
                        ECMandmembrane[i][counter1]=j;  //says what element is interacting !
                        
                        if(counter1 > maximumHowmanyECMCouldInteractWithmembrane )
                        {
                            cout << " error:  maximumHowmanyECMCouldInteractWithmembrane   is NOT enough,increase it! ___the current number is "<< counter1 <<endl;
                        }
                    }
                    
                }
                
                
                
                
                
                
            }
        }
        
        
    }
    
    
    
    
    
    
    
    //
    int j=0;
    int triangle_interact_actively_yet;
    for(int i=0;i<Outer_Membrane_num_of_triangles;i++) // only interact with outer membrane elements
    {
        triangle_interact_actively_yet=0;
        for ( int k=0;k< ECMandmembrane[i][0]  ;k++)
        {
            //  cout<<ECMandmembrane[i][0]<<endl;
            Membrane_triangle_AB[0]=Membrane_Node_Position[ Membrane_triangle_list[i][1]][0]-Membrane_Node_Position[ Membrane_triangle_list[i][0]][0];
            Membrane_triangle_AB[1]=Membrane_Node_Position[ Membrane_triangle_list[i][1]][1]-Membrane_Node_Position[ Membrane_triangle_list[i][0]][1];
            Membrane_triangle_AB[2]=Membrane_Node_Position[ Membrane_triangle_list[i][1]][2]-Membrane_Node_Position[ Membrane_triangle_list[i][0]][2];
            Membrane_triangle_AC[0]=Membrane_Node_Position[ Membrane_triangle_list[i][2]][0]-Membrane_Node_Position[ Membrane_triangle_list[i][0]][0];
            Membrane_triangle_AC[1]=Membrane_Node_Position[ Membrane_triangle_list[i][2]][1]-Membrane_Node_Position[ Membrane_triangle_list[i][0]][1];
            Membrane_triangle_AC[2]=Membrane_Node_Position[ Membrane_triangle_list[i][2]][2]-Membrane_Node_Position[ Membrane_triangle_list[i][0]][2];
            crossvector(Membrane_triangle_ABxAC, Membrane_triangle_AB, Membrane_triangle_AC);
            ECM_triangle_ABxAC_length=vectorlength(Membrane_triangle_ABxAC);
            Membrane_triangle_ABxAC[0]=Membrane_triangle_ABxAC[0]*Membrane_Normal_direction[i][1]/ECM_triangle_ABxAC_length;
            Membrane_triangle_ABxAC[1]=Membrane_triangle_ABxAC[1]*Membrane_Normal_direction[i][1]/ECM_triangle_ABxAC_length;
            Membrane_triangle_ABxAC[2]=Membrane_triangle_ABxAC[2]*Membrane_Normal_direction[i][1]/ECM_triangle_ABxAC_length;  // it is the norml of the membrane element
            
            
            //com of membrane element
            Membrane_triangle_COM[0]= ( Membrane_Node_Position[ Membrane_triangle_list[i][0]][0] + Membrane_Node_Position[ Membrane_triangle_list[i][1]][0] + Membrane_Node_Position[ Membrane_triangle_list[i][2]][0]) /3.0;
            Membrane_triangle_COM[1]= ( Membrane_Node_Position[ Membrane_triangle_list[i][0]][1] + Membrane_Node_Position[ Membrane_triangle_list[i][1]][1] + Membrane_Node_Position[ Membrane_triangle_list[i][2]][1]) /3.0;
            Membrane_triangle_COM[2]= ( Membrane_Node_Position[ Membrane_triangle_list[i][0]][2] +Membrane_Node_Position[ Membrane_triangle_list[i][1]][2] + Membrane_Node_Position[ Membrane_triangle_list[i][2]][2]) /3.0;
            
            j=ECMandmembrane[i][k+1];
            temp_ECM_node_A=ECM_surface_triangle_list[j][0];
            temp_ECM_node_B=ECM_surface_triangle_list[j][1];
            temp_ECM_node_C=ECM_surface_triangle_list[j][2];
            
            ECM_triangle_AB[0]=ECM_Node_Position[temp_ECM_node_B][0]-ECM_Node_Position[temp_ECM_node_A][0];
            ECM_triangle_AB[1]=ECM_Node_Position[temp_ECM_node_B][1]-ECM_Node_Position[temp_ECM_node_A][1];
            ECM_triangle_AB[2]=ECM_Node_Position[temp_ECM_node_B][2]-ECM_Node_Position[temp_ECM_node_A][2];
            ECM_triangle_AC[0]=ECM_Node_Position[temp_ECM_node_C][0]-ECM_Node_Position[temp_ECM_node_A][0];
            ECM_triangle_AC[1]=ECM_Node_Position[temp_ECM_node_C][1]-ECM_Node_Position[temp_ECM_node_A][1];
            ECM_triangle_AC[2]=ECM_Node_Position[temp_ECM_node_C][2]-ECM_Node_Position[temp_ECM_node_A][2];
            crossvector(ECM_triangle_ABxAC_unit_vector, ECM_triangle_AB, ECM_triangle_AC);
            ECM_triangle_ABxAC_length=vectorlength(ECM_triangle_ABxAC_unit_vector);
            ECM_triangle_ABxAC_unit_vector[0]=ECM_triangle_ABxAC_unit_vector[0]/ECM_triangle_ABxAC_length;
            ECM_triangle_ABxAC_unit_vector[1]=ECM_triangle_ABxAC_unit_vector[1]/ECM_triangle_ABxAC_length;
            ECM_triangle_ABxAC_unit_vector[2]=ECM_triangle_ABxAC_unit_vector[2]/ECM_triangle_ABxAC_length;
            
            
            // check the condition of  ECM and COM of membrane beig close enough to ech othe
            Membrane_triangle_COM_ECM_distance_vec[0]= Membrane_triangle_COM[0] -ECM_Node_Position[temp_ECM_node_A][0];
            Membrane_triangle_COM_ECM_distance_vec[1]= Membrane_triangle_COM[1] -ECM_Node_Position[temp_ECM_node_A][1];
            Membrane_triangle_COM_ECM_distance_vec[2]= Membrane_triangle_COM[2] -ECM_Node_Position[temp_ECM_node_A][2];
            membrane_tri_COM_ECM_distance_normal_component=innerproduct(Membrane_triangle_COM_ECM_distance_vec, ECM_triangle_ABxAC_unit_vector);
            
            
            
            if( abs(membrane_tri_COM_ECM_distance_normal_component)< epsilonECM*strechingInteractionTill  )
            {
                D[0]=Membrane_triangle_COM[0] - membrane_tri_COM_ECM_distance_normal_component * ECM_triangle_ABxAC_unit_vector[0];
                D[1]=Membrane_triangle_COM[1] - membrane_tri_COM_ECM_distance_normal_component * ECM_triangle_ABxAC_unit_vector[1];
                D[2]=Membrane_triangle_COM[2] - membrane_tri_COM_ECM_distance_normal_component * ECM_triangle_ABxAC_unit_vector[2];
                
                Membrane_ECM_interactionHELPER(ECM_Node_Position, temp_ECM_node_A, temp_ECM_node_B, temp_ECM_node_C, ECM_triangle_ABxAC_unit_vector,D,a1,a2,a3);
                
                if(a1>=0 & a2>=0 & a3>=0 )  // if the com is inside:
                {
                    if( innerproduct( ECM_triangle_ABxAC_unit_vector, Membrane_triangle_ABxAC ) <=0.0  )  // and if the normal of membrne is pointinf toward the ECM
                    {
                        directionIsInward=1;
                    }
                    else
                    {
                        directionIsInward=0;
                    }
                    
                    a1a2a3=a1+a2+a3;
                    
                    PD[0]=D[0]- Membrane_triangle_COM[0];
                    PD[1]=D[1]- Membrane_triangle_COM[1];
                    PD[2]=D[2]- Membrane_triangle_COM[2];// PD is ortogonal to ECM and has orthogonal length between ECM and membrane
                    t1= (epsilonECM/membrane_tri_COM_ECM_distance_normal_component);
                    r2=t1*t1;
                    r6=r2*r2*r2;
                    
                    flength=4*sigmaECM*( 12.0*r6*r6 - 6.0*r6  )/(membrane_tri_COM_ECM_distance_normal_component);  ///???
                    f[0]=flength* PD[0]/membrane_tri_COM_ECM_distance_normal_component;
                    f[1]=flength* PD[1]/membrane_tri_COM_ECM_distance_normal_component;
                    f[2]=flength* PD[2]/membrane_tri_COM_ECM_distance_normal_component;
                    
                    if( (flength > 0) || (flength < 0 && (   (AdhesiveIslandOffOrOn==0) ||   (j<NumberOfAdhesiveIslandNodes) )) )   // attractive Islands? : repulsive in all ECM and attractive only in attractive islands!
                    {
                        
                        if(flength > fECMcuttoff)
                        {
                            flength=fECMcuttoff;
                        }
                        if(flength <- fECMcuttoff)
                        {
                            flength=-fECMcuttoff;
                        }
                        
                        
                        
                        if(flength<0)
                        {
                            attachingforce=2*flength;
                            attachingforce=-ECMstrechingforce;
                            if( flength<-ECMstrechingforce )
                            {
                                attachingforce=-ECMstrechingforce;
                            }
                        }
                        
                        
                        //_____________
                        if (migratingcell ==1)
                        {
                            if(calculateAngleFromCOM ==0)
                            {
                                
                                Migrating___force_amp=-miu_migration*  abs(  pow(abs (innerproduct(Membrane_triangle_ABxAC, cell_polarizatioz_direction_migration)),(etha_migration) ) )   * sign_function(innerproduct(Membrane_triangle_ABxAC,cell_polarizatioz_direction_migration))   ;
                                
                                
                                
                                
                                
                            }
                            else if(calculateAngleFromCOM ==1)
                            {
                                comTriVec[0]=Cellcom[0]-Membrane_triangle_COM[0];
                                comTriVec[1]=Cellcom[1]-Membrane_triangle_COM[1];
                                comTriVec[2]=Cellcom[2]-Membrane_triangle_COM[2];
                                lvec=-vectorlength(comTriVec);
                                comTriVec[0]= comTriVec[0]/lvec;
                                comTriVec[1]= comTriVec[1]/lvec;
                                comTriVec[2]= comTriVec[2]/lvec;
                                
                                Migrating___force_amp=-miu_migration*  abs(  pow(abs (innerproduct(comTriVec,cell_polarizatioz_direction_migration)),(etha_migration) ) )   * sign_function(innerproduct(comTriVec,cell_polarizatioz_direction_migration))   ;
                                
                                
                                
                            }
                            
                        }
                        
                        //normal2[0]=normal[0]-innerproduct(n,normal)*n[0]; // only the paraller component of normal is needed in streching force
                        //normal2[1]=normal[1]-innerproduct(n,normal)*n[1]; // only the paraller component of normal is needed in streching force
                        //normal2[2]=normal[2]-innerproduct(n,normal)*n[2]; // only the paraller component of normal is needed in streching force
                        
                        normal2[0]=Membrane_triangle_ABxAC[0]+ECM_triangle_ABxAC_unit_vector[0]; // only the paraller component of normal is needed in streching force
                        normal2[1]=Membrane_triangle_ABxAC[1]+ECM_triangle_ABxAC_unit_vector[1]; // only the paraller component of normal is needed in streching force
                        normal2[2]=Membrane_triangle_ABxAC[2]+ECM_triangle_ABxAC_unit_vector[2]; // only the paraller component of normal is needed in streching force
                        
                        
                        
                        
                        
                        
                        
                        ///_________________________________Active_force____INTERNAL force_______________///
                        
                        if( flength< 0.0 && abs(membrane_tri_COM_ECM_distance_normal_component)> epsilonECM*strechingInteractionFrom )
                        {
                            if(triangle_interact_actively_yet==0)
                            { triangle_interact_actively_yet=1;
                                //  temp2=sqrt( (Membrane_Node_Position[index][0]-sourcePoint[0])*(Membrane_Node_Position[index][0]-sourcePoint[0])+(Membrane_Node_Position[index][1]-sourcePoint[1])*(Membrane_Node_Position[index][1]-sourcePoint[1])+(Membrane_Node_Position[index][2]-sourcePoint[2])*(Membrane_Node_Position[index][2]-sourcePoint[2]) );
                                //  temp1=asymmetryConst*exp( -temp2*temp2/(pointRange*pointRange) );
                                
                                
                                //Strecing Force
                                
                                fx=attachingforce*normal2[0]/3.0;
                                fy=attachingforce*normal2[1]/3.0;
                                fz=asymmetryConstInZ*attachingforce*normal2[2]/3.0;
                                
                                
                                index= Membrane_triangle_list[i][0];
                                Membrane_Node_Force[index][0] += fx;
                                Membrane_Node_Force[index][1] += fy;
                                Membrane_Node_Force[index][2] += fz;
                                //f-com-net=0
                                Fnet__ActiveList[numberFnet__ActiveList]=i;
                                numberFnet__ActiveList=numberFnet__ActiveList+1;
                                Fnet__Active[0]=Fnet__Active[0]+fx;
                                Fnet__Active[1]=Fnet__Active[1]+fy;
                                Fnet__Active[2]=Fnet__Active[2]+fz;
                                
                                
                                index= Membrane_triangle_list[i][1];
                                Membrane_Node_Force[index][0] += fx;
                                Membrane_Node_Force[index][1] += fy;
                                Membrane_Node_Force[index][2] += fz;
                                //f-com-net=0
                                Fnet__Active[0]=Fnet__Active[0]+fx;
                                Fnet__Active[1]=Fnet__Active[1]+fy;
                                Fnet__Active[2]=Fnet__Active[2]+fz;
                                
                                
                                index= Membrane_triangle_list[i][2];
                                Membrane_Node_Force[index][0] += fx;
                                Membrane_Node_Force[index][1] += fy;
                                Membrane_Node_Force[index][2] += fz;
                                //f-com-net=0
                                Fnet__Active[0]=Fnet__Active[0]+fx;
                                Fnet__Active[1]=Fnet__Active[1]+fy;
                                Fnet__Active[2]=Fnet__Active[2]+fz;
                                
                                
                                //Migration Force  ____ not normalized yet
                                if(migratingcell==1)
                                {
                                    if(Migrating___force_amp<0)
                                    {
                                        fx=Migrating___force_amp*normal2[0]/3.0;
                                        fy=Migrating___force_amp*normal2[1]/3.0;
                                        fz=Migrating___force_amp*normal2[2]/3.0;
                                        
                                        index= Membrane_triangle_list[i][0];
                                        Migrating___force_Vec[0][index]=Migrating___force_Vec[0][index]+fx;
                                        Migrating___force_Vec[1][index]=Migrating___force_Vec[1][index]+fy;
                                        Migrating___force_Vec[2][index]=Migrating___force_Vec[2][index]+fz;
                                        Total_generated_force=Total_generated_force+sqrt(fx*fx+fy*fy+fz*fz);
                                        
                                        index= Membrane_triangle_list[i][1];
                                        Migrating___force_Vec[0][index]=Migrating___force_Vec[0][index]+fx;
                                        Migrating___force_Vec[1][index]=Migrating___force_Vec[1][index]+fy;
                                        Migrating___force_Vec[2][index]=Migrating___force_Vec[2][index]+fz;
                                        Total_generated_force=Total_generated_force+sqrt(fx*fx+fy*fy+fz*fz);
                                        
                                        
                                        index= Membrane_triangle_list[i][2];
                                        Migrating___force_Vec[0][index]=Migrating___force_Vec[0][index]+fx;
                                        Migrating___force_Vec[1][index]=Migrating___force_Vec[1][index]+fy;
                                        Migrating___force_Vec[2][index]=Migrating___force_Vec[2][index]+fz;
                                        Total_generated_force=Total_generated_force+sqrt(fx*fx+fy*fy+fz*fz);
                                    }
                                    
                                    
                                    
                                    if(Migrating___force_amp>0)
                                    {
                                        fx=Back_Contraction_Coefficient*Migrating___force_amp*normal2[0]/3.0;
                                        fy=Back_Contraction_Coefficient*Migrating___force_amp*normal2[1]/3.0;
                                        fz=Back_Contraction_Coefficient*Migrating___force_amp*normal2[2]/3.0;
                                        
                                        index= Membrane_triangle_list[i][0];
                                        Migrating___force_Vec[0][index]=Migrating___force_Vec[0][index]+fx;
                                        Migrating___force_Vec[1][index]=Migrating___force_Vec[1][index]+fy;
                                        Migrating___force_Vec[2][index]=Migrating___force_Vec[2][index]+fz;
                                        Total_generated_force=Total_generated_force+sqrt(fx*fx+fy*fy+fz*fz);
                                        
                                        index= Membrane_triangle_list[i][1];
                                        Migrating___force_Vec[0][index]=Migrating___force_Vec[0][index]+fx;
                                        Migrating___force_Vec[1][index]=Migrating___force_Vec[1][index]+fy;
                                        Migrating___force_Vec[2][index]=Migrating___force_Vec[2][index]+fz;
                                        Total_generated_force=Total_generated_force+sqrt(fx*fx+fy*fy+fz*fz);
                                        
                                        
                                        index= Membrane_triangle_list[i][2];
                                        Migrating___force_Vec[0][index]=Migrating___force_Vec[0][index]+fx;
                                        Migrating___force_Vec[1][index]=Migrating___force_Vec[1][index]+fy;
                                        Migrating___force_Vec[2][index]=Migrating___force_Vec[2][index]+fz;
                                        Total_generated_force=Total_generated_force+sqrt(fx*fx+fy*fy+fz*fz);
                                    }
                                    
                                }
                            }
                        }
                        
                        ///_________________________________Active_force___________________///
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        ///_________________________________ECM-mem-LJforce___________________///
                        
                        if(ECMrepultionforceToTrianglesOrnode==0)
                        {
                            if( abs(membrane_tri_COM_ECM_distance_normal_component)< epsilonECM*2.0  && directionIsInward==1 )
                                
                                if( (flength>0) ||( flength <0 && ECM_LJ_just_Repultion==0 )  )
                                {
                                    {
                                        
                                        index= Membrane_triangle_list[i][0];
                                        Membrane_Node_Force[index][0] += -flength*ECM_triangle_ABxAC_unit_vector[0]/3.0;
                                        Membrane_Node_Force[index][1] += -flength*ECM_triangle_ABxAC_unit_vector[1]/3.0;
                                        Membrane_Node_Force[index][2] += -flength*ECM_triangle_ABxAC_unit_vector[2]/3.0;
                                        
                                        
                                        
                                        index= Membrane_triangle_list[i][1];
                                        Membrane_Node_Force[index][0] += -flength*ECM_triangle_ABxAC_unit_vector[0]/3.0;
                                        Membrane_Node_Force[index][1] += -flength*ECM_triangle_ABxAC_unit_vector[1]/3.0;
                                        Membrane_Node_Force[index][2] += -flength*ECM_triangle_ABxAC_unit_vector[2]/3.0;
                                        
                                        
                                        
                                        index= Membrane_triangle_list[i][2];
                                        Membrane_Node_Force[index][0] += -flength*ECM_triangle_ABxAC_unit_vector[0]/3.0;
                                        Membrane_Node_Force[index][1] += -flength*ECM_triangle_ABxAC_unit_vector[1]/3.0;
                                        Membrane_Node_Force[index][2] += -flength*ECM_triangle_ABxAC_unit_vector[2]/3.0;
                                        
                                        ECM_Node_Force[temp_ECM_node_A][0] += flength*ECM_triangle_ABxAC_unit_vector[0] *a1/a1a2a3;
                                        ECM_Node_Force[temp_ECM_node_A][1] += flength*ECM_triangle_ABxAC_unit_vector[1] *a1/a1a2a3;
                                        ECM_Node_Force[temp_ECM_node_A][2] += flength*ECM_triangle_ABxAC_unit_vector[2] *a1/a1a2a3;
                                        
                                        ECM_Node_Force[temp_ECM_node_B][0] += flength*ECM_triangle_ABxAC_unit_vector[0] *a2/a1a2a3;
                                        ECM_Node_Force[temp_ECM_node_B][1] += flength*ECM_triangle_ABxAC_unit_vector[1] *a2/a1a2a3;
                                        ECM_Node_Force[temp_ECM_node_B][2] += flength*ECM_triangle_ABxAC_unit_vector[2] *a2/a1a2a3;
                                        
                                        ECM_Node_Force[temp_ECM_node_C][0] += flength*ECM_triangle_ABxAC_unit_vector[0] *a3/a1a2a3;
                                        ECM_Node_Force[temp_ECM_node_C][1] += flength*ECM_triangle_ABxAC_unit_vector[1] *a3/a1a2a3;
                                        ECM_Node_Force[temp_ECM_node_C][2] += flength*ECM_triangle_ABxAC_unit_vector[2] *a3/a1a2a3;
                                    }
                                }
                        }
                        ///_________________________________ECM-mem-LJforce___________________///
                        
                        
                        
                    }  //islands!
                    
                    
                    
                    
                }
            }  //...
            
            
            
            if(ECMrepultionforceToTrianglesOrnode==1)
            {
                
                
                //-------->  Tri 0
                index=Membrane_triangle_list[i][0];
                
                Membrane_triangle_COM_ECM_distance_vec[0]= Membrane_Node_Position[index][0] -ECM_Node_Position[temp_ECM_node_A][0];
                Membrane_triangle_COM_ECM_distance_vec[1]= Membrane_Node_Position[index][1] -ECM_Node_Position[temp_ECM_node_A][1];
                Membrane_triangle_COM_ECM_distance_vec[2]= Membrane_Node_Position[index][2] -ECM_Node_Position[temp_ECM_node_A][2];
                membrane_tri_COM_ECM_distance_normal_component=innerproduct(Membrane_triangle_COM_ECM_distance_vec, ECM_triangle_ABxAC_unit_vector);
                
                if( abs(membrane_tri_COM_ECM_distance_normal_component)< epsilonECM*2.0  )
                {
                    D[0]=Membrane_Node_Position[ index][0]  - membrane_tri_COM_ECM_distance_normal_component * ECM_triangle_ABxAC_unit_vector[0];
                    D[1]=Membrane_Node_Position[ index][1]  - membrane_tri_COM_ECM_distance_normal_component * ECM_triangle_ABxAC_unit_vector[1];
                    D[2]=Membrane_Node_Position[ index][2]  - membrane_tri_COM_ECM_distance_normal_component * ECM_triangle_ABxAC_unit_vector[2];
                    
                    Membrane_ECM_interactionHELPER(ECM_Node_Position, temp_ECM_node_A, temp_ECM_node_B, temp_ECM_node_C, ECM_triangle_ABxAC_unit_vector,D,a1,a2,a3);
                    
                    if(a1>=0 & a2>=0 & a3>=0 )  // if the com is inside:
                    {
                        if( innerproduct(ECM_triangle_ABxAC_unit_vector, Membrane_triangle_ABxAC ) <=0.0  )  // and if the normal of membrne is pointinf toward the ECM
                        {
                            a1a2a3=a1+a2+a3;
                            PD[0]=D[0]- Membrane_Node_Position[ index][0];
                            PD[1]=D[1]- Membrane_Node_Position[ index][1];
                            PD[2]=D[2]- Membrane_Node_Position[ index][2];// PD is ortogonal to ECM and has orthogonal length between ECM and membrane
                            t1= (epsilonECM/membrane_tri_COM_ECM_distance_normal_component);
                            r2=t1*t1;
                            r6=r2*r2*r2;
                            
                            flength=4*sigmaECM*( 12.0*r6*r6 - 6.0*r6  )/(membrane_tri_COM_ECM_distance_normal_component);  ///???
                            f[0]=flength* PD[0]/membrane_tri_COM_ECM_distance_normal_component;
                            f[1]=flength* PD[1]/membrane_tri_COM_ECM_distance_normal_component;
                            f[2]=flength* PD[2]/membrane_tri_COM_ECM_distance_normal_component;
                            
                            
                            if(flength > fECMcuttoff)
                            {
                                flength=fECMcuttoff;
                            }
                            if(flength <- fECMcuttoff)
                            {
                                flength=-fECMcuttoff;
                            }
                            
                            if( (flength>0) ||( flength <0 && ECM_LJ_just_Repultion==0 )  )
                            {
                                Membrane_Node_Force[index][0] += -flength*ECM_triangle_ABxAC_unit_vector[0];
                                Membrane_Node_Force[index][1] += -flength*ECM_triangle_ABxAC_unit_vector[1];
                                Membrane_Node_Force[index][2] += -flength*ECM_triangle_ABxAC_unit_vector[2];
                                
                                ECM_Node_Force[temp_ECM_node_A][0] += flength*ECM_triangle_ABxAC_unit_vector[0] *a1/a1a2a3;
                                ECM_Node_Force[temp_ECM_node_A][1] += flength*ECM_triangle_ABxAC_unit_vector[1] *a1/a1a2a3;
                                ECM_Node_Force[temp_ECM_node_A][2] += flength*ECM_triangle_ABxAC_unit_vector[2] *a1/a1a2a3;
                                
                                ECM_Node_Force[temp_ECM_node_B][0] += flength*ECM_triangle_ABxAC_unit_vector[0] *a2/a1a2a3;
                                ECM_Node_Force[temp_ECM_node_B][1] += flength*ECM_triangle_ABxAC_unit_vector[1] *a2/a1a2a3;
                                ECM_Node_Force[temp_ECM_node_B][2] += flength*ECM_triangle_ABxAC_unit_vector[2] *a2/a1a2a3;
                                
                                ECM_Node_Force[temp_ECM_node_C][0] += flength*ECM_triangle_ABxAC_unit_vector[0] *a3/a1a2a3;
                                ECM_Node_Force[temp_ECM_node_C][1] += flength*ECM_triangle_ABxAC_unit_vector[1] *a3/a1a2a3;
                                ECM_Node_Force[temp_ECM_node_C][2] += flength*ECM_triangle_ABxAC_unit_vector[2] *a3/a1a2a3;
                            }
                        }
                        
                        
                        
                        
                        
                        
                        
                        
                        
                    }
                }
                
                
                //-------->  Tri 1
                index=Membrane_triangle_list[i][1];
                
                Membrane_triangle_COM_ECM_distance_vec[0]= Membrane_Node_Position[index][0] -ECM_Node_Position[temp_ECM_node_A][0];
                Membrane_triangle_COM_ECM_distance_vec[1]= Membrane_Node_Position[index][1] -ECM_Node_Position[temp_ECM_node_A][1];
                Membrane_triangle_COM_ECM_distance_vec[2]= Membrane_Node_Position[index][2] -ECM_Node_Position[temp_ECM_node_A][2];
                membrane_tri_COM_ECM_distance_normal_component=innerproduct(Membrane_triangle_COM_ECM_distance_vec, ECM_triangle_ABxAC_unit_vector);
                
                if( abs(membrane_tri_COM_ECM_distance_normal_component)< epsilonECM*2.0  )
                {
                    D[0]=Membrane_Node_Position[ index][0]  - membrane_tri_COM_ECM_distance_normal_component * ECM_triangle_ABxAC_unit_vector[0];
                    D[1]=Membrane_Node_Position[ index][1]  - membrane_tri_COM_ECM_distance_normal_component * ECM_triangle_ABxAC_unit_vector[1];
                    D[2]=Membrane_Node_Position[ index][2]  - membrane_tri_COM_ECM_distance_normal_component * ECM_triangle_ABxAC_unit_vector[2];
                    
                    Membrane_ECM_interactionHELPER(ECM_Node_Position, temp_ECM_node_A, temp_ECM_node_B, temp_ECM_node_C, ECM_triangle_ABxAC_unit_vector,D,a1,a2,a3);
                    
                    if(a1>=0 & a2>=0 & a3>=0 )  // if the com is inside:
                    {
                        if( innerproduct(ECM_triangle_ABxAC_unit_vector, Membrane_triangle_ABxAC ) <=0.0  )  // and if the normal of membrne is pointinf toward the ECM
                        {
                            a1a2a3=a1+a2+a3;
                            PD[0]=D[0]- Membrane_Node_Position[ index][0];
                            PD[1]=D[1]- Membrane_Node_Position[ index][1];
                            PD[2]=D[2]- Membrane_Node_Position[ index][2];// PD is ortogonal to ECM and has orthogonal length between ECM and membrane
                            t1= (epsilonECM/membrane_tri_COM_ECM_distance_normal_component);
                            r2=t1*t1;
                            r6=r2*r2*r2;
                            
                            flength=4*sigmaECM*( 12.0*r6*r6 - 6.0*r6  )/(membrane_tri_COM_ECM_distance_normal_component);  ///???
                            f[0]=flength* PD[0]/membrane_tri_COM_ECM_distance_normal_component;
                            f[1]=flength* PD[1]/membrane_tri_COM_ECM_distance_normal_component;
                            f[2]=flength* PD[2]/membrane_tri_COM_ECM_distance_normal_component;
                            
                            
                            if(flength > fECMcuttoff)
                            {
                                flength=fECMcuttoff;
                            }
                            if(flength <- fECMcuttoff)
                            {
                                flength=-fECMcuttoff;
                            }
                            
                            if( (flength>0) ||( flength <0 && ECM_LJ_just_Repultion==0 )  )
                            {
                                Membrane_Node_Force[index][0] += -flength*ECM_triangle_ABxAC_unit_vector[0];
                                Membrane_Node_Force[index][1] += -flength*ECM_triangle_ABxAC_unit_vector[1];
                                Membrane_Node_Force[index][2] += -flength*ECM_triangle_ABxAC_unit_vector[2];
                                
                                ECM_Node_Force[temp_ECM_node_A][0] += flength*ECM_triangle_ABxAC_unit_vector[0] *a1/a1a2a3;
                                ECM_Node_Force[temp_ECM_node_A][1] += flength*ECM_triangle_ABxAC_unit_vector[1] *a1/a1a2a3;
                                ECM_Node_Force[temp_ECM_node_A][2] += flength*ECM_triangle_ABxAC_unit_vector[2] *a1/a1a2a3;
                                
                                ECM_Node_Force[temp_ECM_node_B][0] += flength*ECM_triangle_ABxAC_unit_vector[0] *a2/a1a2a3;
                                ECM_Node_Force[temp_ECM_node_B][1] += flength*ECM_triangle_ABxAC_unit_vector[1] *a2/a1a2a3;
                                ECM_Node_Force[temp_ECM_node_B][2] += flength*ECM_triangle_ABxAC_unit_vector[2] *a2/a1a2a3;
                                
                                ECM_Node_Force[temp_ECM_node_C][0] += flength*ECM_triangle_ABxAC_unit_vector[0] *a3/a1a2a3;
                                ECM_Node_Force[temp_ECM_node_C][1] += flength*ECM_triangle_ABxAC_unit_vector[1] *a3/a1a2a3;
                                ECM_Node_Force[temp_ECM_node_C][2] += flength*ECM_triangle_ABxAC_unit_vector[2] *a3/a1a2a3;
                            }
                            
                        }
                        
                        
                        
                        
                        
                        
                        
                        
                        
                    }
                }
                //-------->  Tri 2
                index=Membrane_triangle_list[i][2];
                
                Membrane_triangle_COM_ECM_distance_vec[0]= Membrane_Node_Position[index][0] -ECM_Node_Position[temp_ECM_node_A][0];
                Membrane_triangle_COM_ECM_distance_vec[1]= Membrane_Node_Position[index][1] -ECM_Node_Position[temp_ECM_node_A][1];
                Membrane_triangle_COM_ECM_distance_vec[2]= Membrane_Node_Position[index][2] -ECM_Node_Position[temp_ECM_node_A][2];
                membrane_tri_COM_ECM_distance_normal_component=innerproduct(Membrane_triangle_COM_ECM_distance_vec, ECM_triangle_ABxAC_unit_vector);
                
                if( abs(membrane_tri_COM_ECM_distance_normal_component)< epsilonECM*2.0  )
                {
                    D[0]=Membrane_Node_Position[ index][0]  - membrane_tri_COM_ECM_distance_normal_component * ECM_triangle_ABxAC_unit_vector[0];
                    D[1]=Membrane_Node_Position[ index][1]  - membrane_tri_COM_ECM_distance_normal_component * ECM_triangle_ABxAC_unit_vector[1];
                    D[2]=Membrane_Node_Position[ index][2]  - membrane_tri_COM_ECM_distance_normal_component * ECM_triangle_ABxAC_unit_vector[2];
                    
                    Membrane_ECM_interactionHELPER(ECM_Node_Position, temp_ECM_node_A, temp_ECM_node_B, temp_ECM_node_C, ECM_triangle_ABxAC_unit_vector,D,a1,a2,a3);
                    
                    if(a1>=0 & a2>=0 & a3>=0 )  // if the com is inside:
                    {
                        if( innerproduct(ECM_triangle_ABxAC_unit_vector, Membrane_triangle_ABxAC ) <=0.0  )  // and if the normal of membrne is pointinf toward the ECM
                        {
                            a1a2a3=a1+a2+a3;
                            PD[0]=D[0]- Membrane_Node_Position[ index][0];
                            PD[1]=D[1]- Membrane_Node_Position[ index][1];
                            PD[2]=D[2]- Membrane_Node_Position[ index][2];// PD is ortogonal to ECM and has orthogonal length between ECM and membrane
                            t1= (epsilonECM/membrane_tri_COM_ECM_distance_normal_component);
                            r2=t1*t1;
                            r6=r2*r2*r2;
                            
                            flength=4*sigmaECM*( 12.0*r6*r6 - 6.0*r6  )/(membrane_tri_COM_ECM_distance_normal_component);  ///???
                            f[0]=flength* PD[0]/membrane_tri_COM_ECM_distance_normal_component;
                            f[1]=flength* PD[1]/membrane_tri_COM_ECM_distance_normal_component;
                            f[2]=flength* PD[2]/membrane_tri_COM_ECM_distance_normal_component;
                            
                            
                            if(flength > fECMcuttoff)
                            {
                                flength=fECMcuttoff;
                            }
                            if(flength <- fECMcuttoff)
                            {
                                flength=-fECMcuttoff;
                            }
                            
                            if( (flength>0) ||( flength <0 && ECM_LJ_just_Repultion==0 )  )
                            {
                                Membrane_Node_Force[index][0] += -flength*ECM_triangle_ABxAC_unit_vector[0];
                                Membrane_Node_Force[index][1] += -flength*ECM_triangle_ABxAC_unit_vector[1];
                                Membrane_Node_Force[index][2] += -flength*ECM_triangle_ABxAC_unit_vector[2];
                                
                                ECM_Node_Force[temp_ECM_node_A][0] += flength*ECM_triangle_ABxAC_unit_vector[0] *a1/a1a2a3;
                                ECM_Node_Force[temp_ECM_node_A][1] += flength*ECM_triangle_ABxAC_unit_vector[1] *a1/a1a2a3;
                                ECM_Node_Force[temp_ECM_node_A][2] += flength*ECM_triangle_ABxAC_unit_vector[2] *a1/a1a2a3;
                                
                                ECM_Node_Force[temp_ECM_node_B][0] += flength*ECM_triangle_ABxAC_unit_vector[0] *a2/a1a2a3;
                                ECM_Node_Force[temp_ECM_node_B][1] += flength*ECM_triangle_ABxAC_unit_vector[1] *a2/a1a2a3;
                                ECM_Node_Force[temp_ECM_node_B][2] += flength*ECM_triangle_ABxAC_unit_vector[2] *a2/a1a2a3;
                                
                                ECM_Node_Force[temp_ECM_node_C][0] += flength*ECM_triangle_ABxAC_unit_vector[0] *a3/a1a2a3;
                                ECM_Node_Force[temp_ECM_node_C][1] += flength*ECM_triangle_ABxAC_unit_vector[1] *a3/a1a2a3;
                                ECM_Node_Force[temp_ECM_node_C][2] += flength*ECM_triangle_ABxAC_unit_vector[2] *a3/a1a2a3;
                            }
                            
                        }
                        
                        
                        
                        
                        
                        
                        
                        
                        
                    }
                }
                
            }///if(ECMrepultionforceToTrianglesOrnode==1)
            
            
            
            
            
            
            
            
            
            
            
        }
    }
    
    
    
    //Migration Force  ____  normalizing
    if(migratingcell==1)
    {
        if(controll_migration_with_total_trakhing_force==1)
        {
            alpha_coeff=total_tracking_foece / Total_generated_force;
        }
        else
        {
            alpha_coeff=1.0;
        }
        
        total_tracking_foece_now=alpha_coeff*Total_generated_force;
        
        for(int i=0;i<Outer_Membrane_num_of_triangles;i++)
        {
            
            index= Membrane_triangle_list[i][0];
            fx=alpha_coeff* Migrating___force_Vec[0][index];  //calculate normalized force
            fy=alpha_coeff* Migrating___force_Vec[1][index];
            fz=alpha_coeff* Migrating___force_Vec[2][index];
            Membrane_Node_Force[index][0] += fx;   //add force to meme
            Membrane_Node_Force[index][1] += fy;
            Membrane_Node_Force[index][2] += fz;
            Fnet__Active[0]=Fnet__Active[0]+fx;   //add force to F-net_com to make Fnet=0 later
            Fnet__Active[1]=Fnet__Active[1]+fy;
            Fnet__Active[2]=Fnet__Active[2]+fz;
            
            index= Membrane_triangle_list[i][1];
            fx=alpha_coeff* Migrating___force_Vec[0][index];  //calculate normalized force
            fy=alpha_coeff* Migrating___force_Vec[1][index];
            fz=alpha_coeff* Migrating___force_Vec[2][index];
            Membrane_Node_Force[index][0] += fx;   //add force to meme
            Membrane_Node_Force[index][1] += fy;
            Membrane_Node_Force[index][2] += fz;
            Fnet__Active[0]=Fnet__Active[0]+fx;   //add force to F-net_com to make Fnet=0 later
            Fnet__Active[1]=Fnet__Active[1]+fy;
            Fnet__Active[2]=Fnet__Active[2]+fz;
            
            index= Membrane_triangle_list[i][2];
            fx=alpha_coeff* Migrating___force_Vec[0][index];  //calculate normalized force
            fy=alpha_coeff* Migrating___force_Vec[1][index];
            fz=alpha_coeff* Migrating___force_Vec[2][index];
            Membrane_Node_Force[index][0] += fx;   //add force to meme
            Membrane_Node_Force[index][1] += fy;
            Membrane_Node_Force[index][2] += fz;
            Fnet__Active[0]=Fnet__Active[0]+fx;   //add force to F-net_com to make Fnet=0 later
            Fnet__Active[1]=Fnet__Active[1]+fy;
            Fnet__Active[2]=Fnet__Active[2]+fz;
            
            
            
            
        }
        
    }
    
    
    //Migration Force  ____ not normalizing
    
    
    
    
    
    
    
    //=========================F-com-net=0
    
    int i;
    /*
     if( toOuuterNodes+toButtomNodes+toNucleusnodes !=1.0){
     cout  << "Error:       toouternodes+toButtomNodes+toNucleusnodes !=1.0   correct it\n";
     }*/
    
    
    for(int i1=0;i1<numberFnet__ActiveList;i1++)
    {
        i=Fnet__ActiveList[i1];
        index= Membrane_triangle_list[i][0];
        Membrane_Node_Force[index][0] += - (toButtomNodes)* Fnet__Active[0]/( 3.0*(double)numberFnet__ActiveList );
        Membrane_Node_Force[index][1] += - (toButtomNodes)* Fnet__Active[1]/( 3.0*(double)numberFnet__ActiveList );
        Membrane_Node_Force[index][2] += - (toButtomNodes)* Fnet__Active[2]/( 3.0*(double)numberFnet__ActiveList );
        
        index= Membrane_triangle_list[i][1];
        Membrane_Node_Force[index][0] += - (toButtomNodes)*Fnet__Active[0]/( 3.0*(double)numberFnet__ActiveList );
        Membrane_Node_Force[index][1] += - (toButtomNodes)*Fnet__Active[1]/( 3.0*(double)numberFnet__ActiveList );
        Membrane_Node_Force[index][2] += - (toButtomNodes)*Fnet__Active[2]/( 3.0*(double)numberFnet__ActiveList );
        
        index= Membrane_triangle_list[i][2];
        Membrane_Node_Force[index][0] += - (toButtomNodes)*Fnet__Active[0]/( 3.0*(double)numberFnet__ActiveList );
        Membrane_Node_Force[index][1] += - (toButtomNodes)*Fnet__Active[1]/( 3.0*(double)numberFnet__ActiveList );
        Membrane_Node_Force[index][2] += - (toButtomNodes)*Fnet__Active[2]/( 3.0*(double)numberFnet__ActiveList );
    }
    
    for(int i1=Outer_Membrane_num_of_triangles ;i1<Membrane_num_of_Triangles;i1++)
    {
        index= Membrane_triangle_list[i1][0];
        Membrane_Node_Force[index][0] += -(toNucleusnodes)* Fnet__Active[0]/( 3.0*(double)(Membrane_num_of_Triangles-Outer_Membrane_num_of_triangles) );
        Membrane_Node_Force[index][1] += -(toNucleusnodes)* Fnet__Active[1]/( 3.0*(double)(Membrane_num_of_Triangles-Outer_Membrane_num_of_triangles) );
        Membrane_Node_Force[index][2] += -(toNucleusnodes)* Fnet__Active[2]/( 3.0*(double)(Membrane_num_of_Triangles-Outer_Membrane_num_of_triangles) );
        
        index= Membrane_triangle_list[i1][1];
        Membrane_Node_Force[index][0] += -(toNucleusnodes)* Fnet__Active[0]/( 3.0*(double)(Membrane_num_of_Triangles-Outer_Membrane_num_of_triangles) );
        Membrane_Node_Force[index][1] += -(toNucleusnodes)* Fnet__Active[1]/( 3.0*(double)(Membrane_num_of_Triangles-Outer_Membrane_num_of_triangles) );
        Membrane_Node_Force[index][2] += -(toNucleusnodes)* Fnet__Active[2]/( 3.0*(double)(Membrane_num_of_Triangles-Outer_Membrane_num_of_triangles) );
        
        index= Membrane_triangle_list[i1][2];
        Membrane_Node_Force[index][0] += -(toNucleusnodes)* Fnet__Active[0]/( 3.0*(double)(Membrane_num_of_Triangles-Outer_Membrane_num_of_triangles) );
        Membrane_Node_Force[index][1] += -(toNucleusnodes)* Fnet__Active[1]/( 3.0*(double)(Membrane_num_of_Triangles-Outer_Membrane_num_of_triangles) );
        Membrane_Node_Force[index][2] += -(toNucleusnodes)* Fnet__Active[2]/( 3.0*(double)(Membrane_num_of_Triangles-Outer_Membrane_num_of_triangles) );
    }
    
    
    
    for(int i1=0;i1<Outer_Membrane_num_of_triangles;i1++)
    {
        index= Membrane_triangle_list[i1][0];
        Membrane_Node_Force[index][0] += -(toOuuterNodes)* Fnet__Active[0]/( 3.0*(double)Outer_Membrane_num_of_triangles );
        Membrane_Node_Force[index][1] += -(toOuuterNodes)* Fnet__Active[1]/( 3.0*(double)Outer_Membrane_num_of_triangles );
        Membrane_Node_Force[index][2] += -(toOuuterNodes)* Fnet__Active[2]/( 3.0*(double)Outer_Membrane_num_of_triangles );
        
        index= Membrane_triangle_list[i1][1];
        Membrane_Node_Force[index][0] += -(toOuuterNodes)* Fnet__Active[0]/( 3.0*(double)Outer_Membrane_num_of_triangles );
        Membrane_Node_Force[index][1] += -(toOuuterNodes)* Fnet__Active[1]/( 3.0*(double)Outer_Membrane_num_of_triangles );
        Membrane_Node_Force[index][2] += -(toOuuterNodes)* Fnet__Active[2]/( 3.0*(double)Outer_Membrane_num_of_triangles );
        
        index= Membrane_triangle_list[i1][2];
        Membrane_Node_Force[index][0] += -(toOuuterNodes)* Fnet__Active[0]/( 3.0*(double)Outer_Membrane_num_of_triangles );
        Membrane_Node_Force[index][1] += -(toOuuterNodes)* Fnet__Active[1]/( 3.0*(double)Outer_Membrane_num_of_triangles );
        Membrane_Node_Force[index][2] += -(toOuuterNodes)* Fnet__Active[2]/( 3.0*(double)Outer_Membrane_num_of_triangles );
    }
    
    
    
    
    //=========================F-com-net=0
    
    
    
    
    
    
    
}


void Membrane_ECM_interaction_4(int MD_step, double Membrane_Node_Position[][3], double Membrane_Node_Velocity[][3], double Membrane_Node_Force[][3], double ECM_Node_Position[][3], double ECM_Node_Velocity[][3], double ECM_Node_Force[][3], int ECM_surface_triangle_list[ECM_Surface_num_of_Triangles][3], int Membrane_triangle_list[Membrane_num_of_Triangles][3], int Membrane_Normal_direction[Membrane_num_of_Triangles][2], vector<vector<int> > &ECM_Membrane_Trinagle_neighbours_2, double Cellcom[3], int Outer_Membrane_num_of_triangles, int  Outer_Membrane_num_of_Nodes, int Outer_Membrane_list_of_Nodes[], vector<int>  &Membrane_Edge_triangle_and_ECM_neighbours, double &Total_Potential_Energy)
{
    
    double ECM_triangle_AB[3], ECM_triangle_AC[3], ECM_triangle_ABxAC_unit_vector[3], ECM_triangle_BD[3], ECM_triangle_CD[3], ECM_triangle_AD[3], ECM_triangle_BCxBD[3], ECM_triangle_CAxCD[3], ECM_triangle_ABxAD[3], ECM_triangle_BC[3], ECM_triangle_CA[3], ECM_triangle_COM[3];
    
    double Membrane_triangle_AB[3], Membrane_triangle_AC[3], Membrane_triangle_ABxAC[3], Membrane_triangle_ABxAC_unit_vec[3]; // normal vector of membrane
    int temp_ECM_node_A, temp_ECM_node_B, temp_ECM_node_C;
    int temp_Membrane_Node_A, temp_Membrane_Node_B, temp_Membrane_Node_C;
    double Membrane_trinagle_com_shadow_on_ECM_triangle[3], Membrane_triangle_COM_ECM_distance_vec[3], membrane_tri_COM_ECM_distance_normal_component, Membrane_triangle_COM_ECM_distance_vec_length;
    //    int Membrane_ECM_trinagle_num_of_neighbours=0;
    
    double alpha_A, alpha_B, alpha_C, ECM_triangle_area;
    double reduced_radius, reduced_radius_squared, reduced_radius_cubed, reduced_radius_quint, temp_force_magnitude,f[3];
    double Membrane_triangle_COM[3]; // normal vector of membrane
    
    
    //F-net-com=0
    double Fnet__Active[3];  Fnet__Active[0]=0.0; Fnet__Active[1]=0.0; Fnet__Active[2]=0.0;
    double Fnet_Migration_force[3];  Fnet_Migration_force[0]=0.0; Fnet_Migration_force[1]=0.0; Fnet_Migration_force[2]=0.0;
    double random_number;
    
    //Dettachment step
    if(MD_step % Membrane_ECM_triangle_unbinding_update_step == 0 ) //The Membrane_node-ECM_triangle that interact are updated every 'Membrane_ECM_interaction_update_step' step. After update they will continue to initerct until the list is updated again.
    {
        
        for(int membrane_triangle_index=0; membrane_triangle_index<Outer_Membrane_num_of_triangles; membrane_triangle_index++)
        {
            for ( int ecm_neighbour_index=0; ecm_neighbour_index<ECM_Membrane_Trinagle_neighbours_2[membrane_triangle_index].size();
                 ecm_neighbour_index++)
            {
                temp_Membrane_Node_A=Membrane_triangle_list[membrane_triangle_index][0];
                temp_Membrane_Node_B=Membrane_triangle_list[membrane_triangle_index][1];
                temp_Membrane_Node_C=Membrane_triangle_list[membrane_triangle_index][2];
                
                Membrane_triangle_AB[0]=Membrane_Node_Position[ temp_Membrane_Node_B][0]-Membrane_Node_Position[ temp_Membrane_Node_A][0];
                Membrane_triangle_AB[1]=Membrane_Node_Position[ temp_Membrane_Node_B][1]-Membrane_Node_Position[ temp_Membrane_Node_A][1];
                Membrane_triangle_AB[2]=Membrane_Node_Position[ temp_Membrane_Node_B][2]-Membrane_Node_Position[ temp_Membrane_Node_A][2];
                Membrane_triangle_AC[0]=Membrane_Node_Position[ temp_Membrane_Node_C][0]-Membrane_Node_Position[ temp_Membrane_Node_A][0];
                Membrane_triangle_AC[1]=Membrane_Node_Position[ temp_Membrane_Node_C][1]-Membrane_Node_Position[ temp_Membrane_Node_A][1];
                Membrane_triangle_AC[2]=Membrane_Node_Position[ temp_Membrane_Node_C][2]-Membrane_Node_Position[ temp_Membrane_Node_A][2];
                
                crossvector(Membrane_triangle_ABxAC, Membrane_triangle_AB, Membrane_triangle_AC);
                
                Membrane_triangle_ABxAC[0]*=Membrane_Normal_direction[membrane_triangle_index][1];
                Membrane_triangle_ABxAC[1]*=Membrane_Normal_direction[membrane_triangle_index][1];
                Membrane_triangle_ABxAC[2]*=Membrane_Normal_direction[membrane_triangle_index][1];
                
                
                Membrane_triangle_COM[0]= ( Membrane_Node_Position[ temp_Membrane_Node_A][0] + Membrane_Node_Position[ temp_Membrane_Node_B][0] +Membrane_Node_Position[ temp_Membrane_Node_C][0] ) /3.0;
                Membrane_triangle_COM[1]= ( Membrane_Node_Position[ temp_Membrane_Node_A][1] +Membrane_Node_Position[ temp_Membrane_Node_B][1] +Membrane_Node_Position[ temp_Membrane_Node_C][1] ) /3.0;
                Membrane_triangle_COM[2]= ( Membrane_Node_Position[ temp_Membrane_Node_A][2] +Membrane_Node_Position[ temp_Membrane_Node_B][2] +Membrane_Node_Position[ temp_Membrane_Node_C][2] ) /3.0;
                
                Membrane_triangle_ABxAC_unit_vec[0]=Membrane_triangle_ABxAC[0]/vectorlength(Membrane_triangle_ABxAC);
                Membrane_triangle_ABxAC_unit_vec[1]=Membrane_triangle_ABxAC[1]/vectorlength(Membrane_triangle_ABxAC);
                Membrane_triangle_ABxAC_unit_vec[2]=Membrane_triangle_ABxAC[2]/vectorlength(Membrane_triangle_ABxAC);
                
                temp_ECM_node_A=ECM_surface_triangle_list[ECM_Membrane_Trinagle_neighbours_2[membrane_triangle_index][ecm_neighbour_index]][0];
                temp_ECM_node_B=ECM_surface_triangle_list[ECM_Membrane_Trinagle_neighbours_2[membrane_triangle_index][ecm_neighbour_index]][1];
                temp_ECM_node_C=ECM_surface_triangle_list[ECM_Membrane_Trinagle_neighbours_2[membrane_triangle_index][ecm_neighbour_index]][2];
                
                ECM_triangle_AB[0]=ECM_Node_Position[temp_ECM_node_B][0]-ECM_Node_Position[temp_ECM_node_A][0];
                ECM_triangle_AB[1]=ECM_Node_Position[temp_ECM_node_B][1]-ECM_Node_Position[temp_ECM_node_A][1];
                ECM_triangle_AB[2]=ECM_Node_Position[temp_ECM_node_B][2]-ECM_Node_Position[temp_ECM_node_A][2];
                ECM_triangle_AC[0]=ECM_Node_Position[temp_ECM_node_C][0]-ECM_Node_Position[temp_ECM_node_A][0];
                ECM_triangle_AC[1]=ECM_Node_Position[temp_ECM_node_C][1]-ECM_Node_Position[temp_ECM_node_A][1];
                ECM_triangle_AC[2]=ECM_Node_Position[temp_ECM_node_C][2]-ECM_Node_Position[temp_ECM_node_A][2];
                
                crossvector(ECM_triangle_ABxAC_unit_vector, ECM_triangle_AB, ECM_triangle_AC);
                
                ECM_triangle_ABxAC_unit_vector[0]/=vectorlength(ECM_triangle_ABxAC_unit_vector);
                ECM_triangle_ABxAC_unit_vector[1]/=vectorlength(ECM_triangle_ABxAC_unit_vector);
                ECM_triangle_ABxAC_unit_vector[2]/=vectorlength(ECM_triangle_ABxAC_unit_vector);
                
                ECM_triangle_COM[0]= ( ECM_Node_Position[ temp_ECM_node_A][0] + ECM_Node_Position[ temp_ECM_node_B][0] + ECM_Node_Position[ temp_ECM_node_C][0] ) /3.0;
                ECM_triangle_COM[1]= ( ECM_Node_Position[ temp_ECM_node_A][1] + ECM_Node_Position[ temp_ECM_node_B][1] + ECM_Node_Position[ temp_ECM_node_C][1] ) /3.0;
                ECM_triangle_COM[2]= ( ECM_Node_Position[ temp_ECM_node_A][2] + ECM_Node_Position[ temp_ECM_node_B][2] + ECM_Node_Position[ temp_ECM_node_C][2] ) /3.0;
                
                Membrane_triangle_COM_ECM_distance_vec[0]= Membrane_triangle_COM[0] - ECM_triangle_COM[0];
                Membrane_triangle_COM_ECM_distance_vec[1]= Membrane_triangle_COM[1] - ECM_triangle_COM[1];
                Membrane_triangle_COM_ECM_distance_vec[2]= Membrane_triangle_COM[2] - ECM_triangle_COM[2];
                
                Membrane_triangle_COM_ECM_distance_vec_length=vectorlength(Membrane_triangle_COM_ECM_distance_vec);
                
                if ( (Membrane_triangle_COM_ECM_distance_vec_length > 2.5*epsilonECM)
                    || (innerproduct(ECM_triangle_ABxAC_unit_vector, Membrane_triangle_ABxAC_unit_vec) > 0 )
                    )
                {
                    ECM_Membrane_Trinagle_neighbours_2[membrane_triangle_index].erase(ECM_Membrane_Trinagle_neighbours_2[membrane_triangle_index].begin()+ecm_neighbour_index);
                    ecm_neighbour_index--;
                    break;
                } else {
                    random_number = ((double) rand() / (RAND_MAX));
                    if (random_number < Membrane_ECM_UnBinding_prob) {
                        ECM_Membrane_Trinagle_neighbours_2[membrane_triangle_index].erase(ECM_Membrane_Trinagle_neighbours_2[membrane_triangle_index].begin()+ecm_neighbour_index);
                        ecm_neighbour_index--;
                    }
                }
                
            }//END OF: for ( int j=0;j<ECM_Surface_num_of_Triangles;j++)
        } //END OF: for(int membrane_triangle_index=0; membrane_triangle_index<Outer_Membrane_num_of_triangles; membrane_triangle_index++)
        //    exit (EXIT_FAILURE);
        
    } //END OF: if(istep % Membrane_ECM_interaction_update_step == 0 )
    //    int testing_counter;
    //Attachment step
    if(MD_step % Membrane_ECM_triangle_interaction_update_step == 0 ) //The Membrane_node-ECM_triangle that interact are updated every 'Membrane_ECM_interaction_update_step' step. After update they will continue to initerct until the list is updated again.
        
    {
        //streatching-begin
        if (spreading_flag==1.0) {
            Membrane_Edge_triangle_and_ECM_neighbours.clear();
            Membrane_Edge_triangle_and_ECM_neighbours.resize(Outer_Membrane_num_of_triangles);
            for (int i=0; i<Outer_Membrane_num_of_triangles; i++) {
                Membrane_Edge_triangle_and_ECM_neighbours[i]=-1;
            }
        }
        
        //        int spreading_triangle_count=0;
        //streatching-end
        //        testing_counter=0;
        for(int membrane_triangle_index=0; membrane_triangle_index<Outer_Membrane_num_of_triangles; membrane_triangle_index++)
        {
            bool attachmed_list_flag=false;
            for ( int ecm_triangle_index=0; ecm_triangle_index<ECM_Surface_num_of_Triangles; ecm_triangle_index++)
            {
                for (int neighbour_ecm_index=0; neighbour_ecm_index < ECM_Membrane_Trinagle_neighbours_2[membrane_triangle_index].size(); neighbour_ecm_index++)
                {
                    if (ECM_Membrane_Trinagle_neighbours_2[membrane_triangle_index][neighbour_ecm_index] == ecm_triangle_index) {
                        attachmed_list_flag = true;
                        break;
                    }
                }
                if (attachmed_list_flag != true)
                {
                    temp_Membrane_Node_A=Membrane_triangle_list[membrane_triangle_index][0];
                    temp_Membrane_Node_B=Membrane_triangle_list[membrane_triangle_index][1];
                    temp_Membrane_Node_C=Membrane_triangle_list[membrane_triangle_index][2];
                    
                    Membrane_triangle_AB[0]=Membrane_Node_Position[ temp_Membrane_Node_B][0]-Membrane_Node_Position[ temp_Membrane_Node_A][0];
                    Membrane_triangle_AB[1]=Membrane_Node_Position[ temp_Membrane_Node_B][1]-Membrane_Node_Position[ temp_Membrane_Node_A][1];
                    Membrane_triangle_AB[2]=Membrane_Node_Position[ temp_Membrane_Node_B][2]-Membrane_Node_Position[ temp_Membrane_Node_A][2];
                    Membrane_triangle_AC[0]=Membrane_Node_Position[ temp_Membrane_Node_C][0]-Membrane_Node_Position[ temp_Membrane_Node_A][0];
                    Membrane_triangle_AC[1]=Membrane_Node_Position[ temp_Membrane_Node_C][1]-Membrane_Node_Position[ temp_Membrane_Node_A][1];
                    Membrane_triangle_AC[2]=Membrane_Node_Position[ temp_Membrane_Node_C][2]-Membrane_Node_Position[ temp_Membrane_Node_A][2];
                    
                    crossvector(Membrane_triangle_ABxAC, Membrane_triangle_AB, Membrane_triangle_AC);
                    
                    Membrane_triangle_ABxAC[0]*=Membrane_Normal_direction[membrane_triangle_index][1];
                    Membrane_triangle_ABxAC[1]*=Membrane_Normal_direction[membrane_triangle_index][1];
                    Membrane_triangle_ABxAC[2]*=Membrane_Normal_direction[membrane_triangle_index][1];
                    
                    
                    Membrane_triangle_COM[0]= ( Membrane_Node_Position[ temp_Membrane_Node_A][0] + Membrane_Node_Position[ temp_Membrane_Node_B][0] +Membrane_Node_Position[ temp_Membrane_Node_C][0] ) /3.0;
                    Membrane_triangle_COM[1]= ( Membrane_Node_Position[ temp_Membrane_Node_A][1] +Membrane_Node_Position[ temp_Membrane_Node_B][1] +Membrane_Node_Position[ temp_Membrane_Node_C][1] ) /3.0;
                    Membrane_triangle_COM[2]= ( Membrane_Node_Position[ temp_Membrane_Node_A][2] +Membrane_Node_Position[ temp_Membrane_Node_B][2] +Membrane_Node_Position[ temp_Membrane_Node_C][2] ) /3.0;
                    
                    Membrane_triangle_ABxAC_unit_vec[0]=Membrane_triangle_ABxAC[0]/vectorlength(Membrane_triangle_ABxAC);
                    Membrane_triangle_ABxAC_unit_vec[1]=Membrane_triangle_ABxAC[1]/vectorlength(Membrane_triangle_ABxAC);
                    Membrane_triangle_ABxAC_unit_vec[2]=Membrane_triangle_ABxAC[2]/vectorlength(Membrane_triangle_ABxAC);
                    
                    temp_ECM_node_A=ECM_surface_triangle_list[ecm_triangle_index][0];
                    temp_ECM_node_B=ECM_surface_triangle_list[ecm_triangle_index][1];
                    temp_ECM_node_C=ECM_surface_triangle_list[ecm_triangle_index][2];
                    
                    ECM_triangle_AB[0]=ECM_Node_Position[temp_ECM_node_B][0]-ECM_Node_Position[temp_ECM_node_A][0];
                    ECM_triangle_AB[1]=ECM_Node_Position[temp_ECM_node_B][1]-ECM_Node_Position[temp_ECM_node_A][1];
                    ECM_triangle_AB[2]=ECM_Node_Position[temp_ECM_node_B][2]-ECM_Node_Position[temp_ECM_node_A][2];
                    ECM_triangle_AC[0]=ECM_Node_Position[temp_ECM_node_C][0]-ECM_Node_Position[temp_ECM_node_A][0];
                    ECM_triangle_AC[1]=ECM_Node_Position[temp_ECM_node_C][1]-ECM_Node_Position[temp_ECM_node_A][1];
                    ECM_triangle_AC[2]=ECM_Node_Position[temp_ECM_node_C][2]-ECM_Node_Position[temp_ECM_node_A][2];
                    
                    crossvector(ECM_triangle_ABxAC_unit_vector, ECM_triangle_AB, ECM_triangle_AC);
                    
                    ECM_triangle_ABxAC_unit_vector[0]/=vectorlength(ECM_triangle_ABxAC_unit_vector);
                    ECM_triangle_ABxAC_unit_vector[1]/=vectorlength(ECM_triangle_ABxAC_unit_vector);
                    ECM_triangle_ABxAC_unit_vector[2]/=vectorlength(ECM_triangle_ABxAC_unit_vector);
                    
                    ECM_triangle_COM[0]= ( ECM_Node_Position[ temp_ECM_node_A][0] + ECM_Node_Position[ temp_ECM_node_B][0] + ECM_Node_Position[ temp_ECM_node_C][0] ) /3.0;
                    ECM_triangle_COM[1]= ( ECM_Node_Position[ temp_ECM_node_A][1] + ECM_Node_Position[ temp_ECM_node_B][1] + ECM_Node_Position[ temp_ECM_node_C][1] ) /3.0;
                    ECM_triangle_COM[2]= ( ECM_Node_Position[ temp_ECM_node_A][2] + ECM_Node_Position[ temp_ECM_node_B][2] + ECM_Node_Position[ temp_ECM_node_C][2] ) /3.0;
                    
                    Membrane_triangle_COM_ECM_distance_vec[0]= Membrane_triangle_COM[0] - ECM_triangle_COM[0];
                    Membrane_triangle_COM_ECM_distance_vec[1]= Membrane_triangle_COM[1] - ECM_triangle_COM[1];
                    Membrane_triangle_COM_ECM_distance_vec[2]= Membrane_triangle_COM[2] - ECM_triangle_COM[2];
                    
                    membrane_tri_COM_ECM_distance_normal_component=innerproduct(Membrane_triangle_COM_ECM_distance_vec,ECM_triangle_ABxAC_unit_vector);
                    Membrane_triangle_COM_ECM_distance_vec_length=vectorlength(Membrane_triangle_COM_ECM_distance_vec);
                    
                    if ( (Membrane_triangle_COM_ECM_distance_vec_length < spreading_force_max_range*epsilonECM) &&
                        (innerproduct(ECM_triangle_ABxAC_unit_vector, Membrane_triangle_ABxAC_unit_vec) < spreading_force_cos_triangle_interaction_angle))
                    {
                        Membrane_trinagle_com_shadow_on_ECM_triangle[0]=Membrane_triangle_COM[0] - membrane_tri_COM_ECM_distance_normal_component * ECM_triangle_ABxAC_unit_vector[0];
                        Membrane_trinagle_com_shadow_on_ECM_triangle[1]=Membrane_triangle_COM[1] - membrane_tri_COM_ECM_distance_normal_component * ECM_triangle_ABxAC_unit_vector[1];
                        Membrane_trinagle_com_shadow_on_ECM_triangle[2]=Membrane_triangle_COM[2] - membrane_tri_COM_ECM_distance_normal_component * ECM_triangle_ABxAC_unit_vector[2];
                        
                        ECM_triangle_BD[0]=Membrane_trinagle_com_shadow_on_ECM_triangle[0]-ECM_Node_Position[temp_ECM_node_B][0];
                        ECM_triangle_BD[1]=Membrane_trinagle_com_shadow_on_ECM_triangle[1]-ECM_Node_Position[temp_ECM_node_B][1];
                        ECM_triangle_BD[2]=Membrane_trinagle_com_shadow_on_ECM_triangle[2]-ECM_Node_Position[temp_ECM_node_B][2];
                        
                        ECM_triangle_CD[0]=Membrane_trinagle_com_shadow_on_ECM_triangle[0]-ECM_Node_Position[temp_ECM_node_C][0];
                        ECM_triangle_CD[1]=Membrane_trinagle_com_shadow_on_ECM_triangle[1]-ECM_Node_Position[temp_ECM_node_C][1];
                        ECM_triangle_CD[2]=Membrane_trinagle_com_shadow_on_ECM_triangle[2]-ECM_Node_Position[temp_ECM_node_C][2];
                        
                        ECM_triangle_AD[0]=Membrane_trinagle_com_shadow_on_ECM_triangle[0]-ECM_Node_Position[temp_ECM_node_A][0];
                        ECM_triangle_AD[1]=Membrane_trinagle_com_shadow_on_ECM_triangle[1]-ECM_Node_Position[temp_ECM_node_A][1];
                        ECM_triangle_AD[2]=Membrane_trinagle_com_shadow_on_ECM_triangle[2]-ECM_Node_Position[temp_ECM_node_A][2];
                        
                        ECM_triangle_BC[0]=ECM_Node_Position[temp_ECM_node_C][0]-ECM_Node_Position[temp_ECM_node_B][0];
                        ECM_triangle_BC[1]=ECM_Node_Position[temp_ECM_node_C][1]-ECM_Node_Position[temp_ECM_node_B][1];
                        ECM_triangle_BC[2]=ECM_Node_Position[temp_ECM_node_C][2]-ECM_Node_Position[temp_ECM_node_B][2];
                        
                        ECM_triangle_CA[0]=-ECM_triangle_AC[0];
                        ECM_triangle_CA[1]=-ECM_triangle_AC[1];
                        ECM_triangle_CA[2]=-ECM_triangle_AC[2];
                        
                        crossvector(ECM_triangle_BCxBD, ECM_triangle_BC, ECM_triangle_BD);
                        alpha_A=0.5*innerproduct(ECM_triangle_ABxAC_unit_vector, ECM_triangle_BCxBD);
                        
                        crossvector(ECM_triangle_CAxCD, ECM_triangle_CA, ECM_triangle_CD);
                        alpha_B=0.5*innerproduct(ECM_triangle_ABxAC_unit_vector, ECM_triangle_CAxCD);
                        
                        crossvector(ECM_triangle_ABxAD, ECM_triangle_AB, ECM_triangle_AD);
                        alpha_C=0.5*innerproduct(ECM_triangle_ABxAC_unit_vector, ECM_triangle_ABxAD);
                        
                        if(alpha_A>=0 & alpha_B>=0 & alpha_C>=0)
                        {
                            if ((Membrane_triangle_COM_ECM_distance_vec_length < 2.5*epsilonECM) && (innerproduct(ECM_triangle_ABxAC_unit_vector, Membrane_triangle_ABxAC_unit_vec) < -Membrane_ECM_Max_cos_triangle_interaction_angle)) {
                                random_number = ((double) rand() / (RAND_MAX));
                                if (random_number < Membrane_ECM_Binding_prob) {
                                    ECM_Membrane_Trinagle_neighbours_2[membrane_triangle_index].push_back(ecm_triangle_index);
                                }
                            }
                            if ((Membrane_triangle_COM_ECM_distance_vec_length > spreading_force_min_range*epsilonECM) && (innerproduct(ECM_triangle_ABxAC_unit_vector, Membrane_triangle_ABxAC_unit_vec) > -spreading_force_cos_triangle_interaction_angle) && spreading_flag==1.0) {
                                Membrane_Edge_triangle_and_ECM_neighbours[membrane_triangle_index]=ecm_triangle_index;
                                //                                spreading_triangle_count++;
                            }
                        }
                    } //END OF: if ( abs(membrane_tri_COM_ECM_distance_normal_component) && (innerproduct(ECM_triangle_ABxAC_unit_vector, Membrane_triangle_ABxAC_unit_vec) < -Membrane_ECM_Max_cos_triangle_interaction_angle ))
                }//END OF: if (attachment_list_flag != true)
            }//END OF: for ( int j=0;j<ECM_Surface_num_of_Triangles;j++)
        } //END OF: for(int membrane_triangle_index=0; membrane_triangle_index<Outer_Membrane_num_of_triangles; membrane_triangle_index++)
        //    exit (EXIT_FAILURE);
        //        cout<<"spreading_triangle_count:  "<<spreading_triangle_count<<endl;
    } //END OF: if(istep % Membrane_ECM_interaction_update_step == 0 )
    
    //Force calculation
    for(int temp_trinagle_index=0; temp_trinagle_index<Outer_Membrane_num_of_triangles; temp_trinagle_index++)
    {
        if (ECM_Membrane_Trinagle_neighbours_2[temp_trinagle_index].size()!=0) {
            
            temp_Membrane_Node_A=Membrane_triangle_list[temp_trinagle_index][0];
            temp_Membrane_Node_B=Membrane_triangle_list[temp_trinagle_index][1];
            temp_Membrane_Node_C=Membrane_triangle_list[temp_trinagle_index][2];
            
            Membrane_triangle_AB[0]=Membrane_Node_Position[ temp_Membrane_Node_B][0]-Membrane_Node_Position[ temp_Membrane_Node_A][0];
            Membrane_triangle_AB[1]=Membrane_Node_Position[ temp_Membrane_Node_B][1]-Membrane_Node_Position[ temp_Membrane_Node_A][1];
            Membrane_triangle_AB[2]=Membrane_Node_Position[ temp_Membrane_Node_B][2]-Membrane_Node_Position[ temp_Membrane_Node_A][2];
            
            Membrane_triangle_AC[0]=Membrane_Node_Position[ temp_Membrane_Node_C][0]-Membrane_Node_Position[ temp_Membrane_Node_A][0];
            Membrane_triangle_AC[1]=Membrane_Node_Position[ temp_Membrane_Node_C][1]-Membrane_Node_Position[ temp_Membrane_Node_A][1];
            Membrane_triangle_AC[2]=Membrane_Node_Position[ temp_Membrane_Node_C][2]-Membrane_Node_Position[ temp_Membrane_Node_A][2];
            
            crossvector(Membrane_triangle_ABxAC, Membrane_triangle_AB, Membrane_triangle_AC);
            
            Membrane_triangle_ABxAC_unit_vec[0]=Membrane_triangle_ABxAC[0]*Membrane_Normal_direction[temp_trinagle_index][1]/vectorlength(Membrane_triangle_ABxAC);
            Membrane_triangle_ABxAC_unit_vec[1]=Membrane_triangle_ABxAC[1]*Membrane_Normal_direction[temp_trinagle_index][1]/vectorlength(Membrane_triangle_ABxAC);
            Membrane_triangle_ABxAC_unit_vec[2]=Membrane_triangle_ABxAC[2]*Membrane_Normal_direction[temp_trinagle_index][1]/vectorlength(Membrane_triangle_ABxAC);
            
            Membrane_triangle_COM[0]= ( Membrane_Node_Position[ temp_Membrane_Node_A][0] + Membrane_Node_Position[ temp_Membrane_Node_B][0] + Membrane_Node_Position[ temp_Membrane_Node_C][0]) /3.0;
            Membrane_triangle_COM[1]= ( Membrane_Node_Position[ temp_Membrane_Node_A][1] + Membrane_Node_Position[ temp_Membrane_Node_B][1] + Membrane_Node_Position[ temp_Membrane_Node_C][1]) /3.0;
            Membrane_triangle_COM[2]= ( Membrane_Node_Position[ temp_Membrane_Node_A][2] +Membrane_Node_Position[ temp_Membrane_Node_B][2] + Membrane_Node_Position[ temp_Membrane_Node_C][2]) /3.0;
            for (int neighbour_counter=0; neighbour_counter<ECM_Membrane_Trinagle_neighbours_2[temp_trinagle_index].size(); neighbour_counter++)
            {
                temp_ECM_node_A=ECM_surface_triangle_list[ECM_Membrane_Trinagle_neighbours_2[temp_trinagle_index][neighbour_counter]][0];
                temp_ECM_node_B=ECM_surface_triangle_list[ECM_Membrane_Trinagle_neighbours_2[temp_trinagle_index][neighbour_counter]][1];
                temp_ECM_node_C=ECM_surface_triangle_list[ECM_Membrane_Trinagle_neighbours_2[temp_trinagle_index][neighbour_counter]][2];
                
                ECM_triangle_AB[0]=ECM_Node_Position[temp_ECM_node_B][0]-ECM_Node_Position[temp_ECM_node_A][0];
                ECM_triangle_AB[1]=ECM_Node_Position[temp_ECM_node_B][1]-ECM_Node_Position[temp_ECM_node_A][1];
                ECM_triangle_AB[2]=ECM_Node_Position[temp_ECM_node_B][2]-ECM_Node_Position[temp_ECM_node_A][2];
                ECM_triangle_AC[0]=ECM_Node_Position[temp_ECM_node_C][0]-ECM_Node_Position[temp_ECM_node_A][0];
                ECM_triangle_AC[1]=ECM_Node_Position[temp_ECM_node_C][1]-ECM_Node_Position[temp_ECM_node_A][1];
                ECM_triangle_AC[2]=ECM_Node_Position[temp_ECM_node_C][2]-ECM_Node_Position[temp_ECM_node_A][2];
                
                crossvector(ECM_triangle_ABxAC_unit_vector, ECM_triangle_AB, ECM_triangle_AC);
                
                ECM_triangle_ABxAC_unit_vector[0]=ECM_triangle_ABxAC_unit_vector[0]/vectorlength(ECM_triangle_ABxAC_unit_vector);
                ECM_triangle_ABxAC_unit_vector[1]=ECM_triangle_ABxAC_unit_vector[1]/vectorlength(ECM_triangle_ABxAC_unit_vector);
                ECM_triangle_ABxAC_unit_vector[2]=ECM_triangle_ABxAC_unit_vector[2]/vectorlength(ECM_triangle_ABxAC_unit_vector);
                
                ECM_triangle_COM[0]= ( ECM_Node_Position[ temp_ECM_node_A][0] + ECM_Node_Position[ temp_ECM_node_B][0] + ECM_Node_Position[ temp_ECM_node_C][0] ) /3.0;
                ECM_triangle_COM[1]= ( ECM_Node_Position[ temp_ECM_node_A][1] + ECM_Node_Position[ temp_ECM_node_B][1] + ECM_Node_Position[ temp_ECM_node_C][1] ) /3.0;
                ECM_triangle_COM[2]= ( ECM_Node_Position[ temp_ECM_node_A][2] + ECM_Node_Position[ temp_ECM_node_B][2] + ECM_Node_Position[ temp_ECM_node_C][2] ) /3.0;
                
                Membrane_triangle_COM_ECM_distance_vec[0]= Membrane_triangle_COM[0] - ECM_triangle_COM[0];
                Membrane_triangle_COM_ECM_distance_vec[1]= Membrane_triangle_COM[1] - ECM_triangle_COM[1];
                Membrane_triangle_COM_ECM_distance_vec[2]= Membrane_triangle_COM[2] - ECM_triangle_COM[2];
                
                membrane_tri_COM_ECM_distance_normal_component=innerproduct(Membrane_triangle_COM_ECM_distance_vec, ECM_triangle_ABxAC_unit_vector);
                
                Membrane_trinagle_com_shadow_on_ECM_triangle[0]=Membrane_triangle_COM[0] - membrane_tri_COM_ECM_distance_normal_component * ECM_triangle_ABxAC_unit_vector[0];
                Membrane_trinagle_com_shadow_on_ECM_triangle[1]=Membrane_triangle_COM[1] - membrane_tri_COM_ECM_distance_normal_component * ECM_triangle_ABxAC_unit_vector[1];
                Membrane_trinagle_com_shadow_on_ECM_triangle[2]=Membrane_triangle_COM[2] - membrane_tri_COM_ECM_distance_normal_component * ECM_triangle_ABxAC_unit_vector[2];
                
                Membrane_triangle_COM_ECM_distance_vec_length=vectorlength(Membrane_triangle_COM_ECM_distance_vec);
                
                ECM_triangle_BD[0]=Membrane_trinagle_com_shadow_on_ECM_triangle[0]-ECM_Node_Position[temp_ECM_node_B][0];
                ECM_triangle_BD[1]=Membrane_trinagle_com_shadow_on_ECM_triangle[1]-ECM_Node_Position[temp_ECM_node_B][1];
                ECM_triangle_BD[2]=Membrane_trinagle_com_shadow_on_ECM_triangle[2]-ECM_Node_Position[temp_ECM_node_B][2];
                
                ECM_triangle_CD[0]=Membrane_trinagle_com_shadow_on_ECM_triangle[0]-ECM_Node_Position[temp_ECM_node_C][0];
                ECM_triangle_CD[1]=Membrane_trinagle_com_shadow_on_ECM_triangle[1]-ECM_Node_Position[temp_ECM_node_C][1];
                ECM_triangle_CD[2]=Membrane_trinagle_com_shadow_on_ECM_triangle[2]-ECM_Node_Position[temp_ECM_node_C][2];
                
                ECM_triangle_AD[0]=Membrane_trinagle_com_shadow_on_ECM_triangle[0]-ECM_Node_Position[temp_ECM_node_A][0];
                ECM_triangle_AD[1]=Membrane_trinagle_com_shadow_on_ECM_triangle[1]-ECM_Node_Position[temp_ECM_node_A][1];
                ECM_triangle_AD[2]=Membrane_trinagle_com_shadow_on_ECM_triangle[2]-ECM_Node_Position[temp_ECM_node_A][2];
                
                ECM_triangle_BC[0]=ECM_Node_Position[temp_ECM_node_C][0]-ECM_Node_Position[temp_ECM_node_B][0];
                ECM_triangle_BC[1]=ECM_Node_Position[temp_ECM_node_C][1]-ECM_Node_Position[temp_ECM_node_B][1];
                ECM_triangle_BC[2]=ECM_Node_Position[temp_ECM_node_C][2]-ECM_Node_Position[temp_ECM_node_B][2];
                
                ECM_triangle_CA[0]=-ECM_triangle_AC[0];
                ECM_triangle_CA[1]=-ECM_triangle_AC[1];
                ECM_triangle_CA[2]=-ECM_triangle_AC[2];
                
                crossvector(ECM_triangle_BCxBD, ECM_triangle_BC, ECM_triangle_BD);
                alpha_A=0.5*innerproduct(ECM_triangle_ABxAC_unit_vector, ECM_triangle_BCxBD);
                
                crossvector(ECM_triangle_CAxCD, ECM_triangle_CA, ECM_triangle_CD);
                alpha_B=0.5*innerproduct(ECM_triangle_ABxAC_unit_vector, ECM_triangle_CAxCD);
                
                crossvector(ECM_triangle_ABxAD, ECM_triangle_AB, ECM_triangle_AD);
                alpha_C=0.5*innerproduct(ECM_triangle_ABxAC_unit_vector, ECM_triangle_ABxAD);
                
                //            =============================================================
                ECM_triangle_area = alpha_A + alpha_B + alpha_C;
                
                reduced_radius = (1.5*epsilonECM/Membrane_triangle_COM_ECM_distance_vec_length);
                reduced_radius_squared = reduced_radius*reduced_radius;
                reduced_radius_cubed = reduced_radius_squared*reduced_radius;
                reduced_radius_quint = reduced_radius_cubed*reduced_radius_squared;
                
                temp_force_magnitude = 4.0*Membrane_ECM_interaction_strength*( reduced_radius_quint - reduced_radius_cubed )/(1.5*epsilonECM*Membrane_triangle_COM_ECM_distance_vec_length);
                
                if (energy_calculation_flag==1.0) {
                    Total_Potential_Energy+=Membrane_ECM_interaction_strength*(reduced_radius_squared*reduced_radius_squared-2*reduced_radius_squared);
                }
                //                f[0] = temp_force_magnitude * Membrane_triangle_COM_ECM_distance_vec[0];
                f[1] = temp_force_magnitude * Membrane_triangle_COM_ECM_distance_vec[1];
                //                f[2] = temp_force_magnitude * Membrane_triangle_COM_ECM_distance_vec[2];
                f[0]=0.0;
                f[2]=0.0;
                
                for (int abc_counter=0; abc_counter<3; abc_counter++) {
                    Membrane_Node_Force[temp_Membrane_Node_A][abc_counter] += -f[abc_counter]/3.0;
                    Membrane_Node_Force[temp_Membrane_Node_B][abc_counter] += -f[abc_counter]/3.0;
                    Membrane_Node_Force[temp_Membrane_Node_C][abc_counter] += -f[abc_counter]/3.0;
                    
                    ECM_Node_Force[temp_ECM_node_A][abc_counter] += (alpha_A/ECM_triangle_area)*f[abc_counter];
                    ECM_Node_Force[temp_ECM_node_B][abc_counter] += (alpha_B/ECM_triangle_area)*f[abc_counter];
                    ECM_Node_Force[temp_ECM_node_C][abc_counter] += (alpha_C/ECM_triangle_area)*f[abc_counter];
                    
                }
            }//END OF: for (int neighbour_counter=0; neighbour_counter<ECM_Membrane_Trinagle_neighbours_2[i].size(); neighbour_counter++)
            
            //            spreading force-begin
            if (Membrane_Edge_triangle_and_ECM_neighbours[temp_trinagle_index]!=-1 && spreading_flag==1.0) {
                int temp_ecm_spreading_index=Membrane_Edge_triangle_and_ECM_neighbours[temp_trinagle_index];
                temp_ECM_node_A=ECM_surface_triangle_list[temp_ecm_spreading_index][0];
                temp_ECM_node_B=ECM_surface_triangle_list[temp_ecm_spreading_index][1];
                temp_ECM_node_C=ECM_surface_triangle_list[temp_ecm_spreading_index][2];
                
                f[0] = spreading_force_magnitude * Membrane_triangle_COM_ECM_distance_vec[0];
                f[1] = spreading_force_magnitude * Membrane_triangle_COM_ECM_distance_vec[1];
                f[2] = spreading_force_magnitude * Membrane_triangle_COM_ECM_distance_vec[2];
                
                for (int abc_counter=0; abc_counter<3; abc_counter++) {
                    Membrane_Node_Force[temp_Membrane_Node_A][abc_counter] += -f[abc_counter]/3.0;
                    Membrane_Node_Force[temp_Membrane_Node_B][abc_counter] += -f[abc_counter]/3.0;
                    Membrane_Node_Force[temp_Membrane_Node_C][abc_counter] += -f[abc_counter]/3.0;
                    
                    ECM_Node_Force[temp_ECM_node_A][abc_counter] += f[abc_counter]/3.0;
                    ECM_Node_Force[temp_ECM_node_B][abc_counter] += f[abc_counter]/3.0;
                    ECM_Node_Force[temp_ECM_node_C][abc_counter] += f[abc_counter]/3.0;
                    
                }
            }
            //            spreading force-end
            
            
        }// END OF: if (ECM_Membrane_Trinagle_neighbours[i]!=ECM_Surface_num_of_Triangles)
    }//END OF: for(int i=0;i<Outer_Membrane_num_of_triangles;i++)
    
}

void Membrane_ECM_interaction_5(int MD_step, double Membrane_Node_Position[][3], double Membrane_Node_Velocity[][3], double Membrane_Node_Force[][3], double ECM_Node_Position[][3], double ECM_Node_Velocity[][3], double ECM_Node_Force[][3], int ECM_surface_triangle_list[ECM_Surface_num_of_Triangles][3], int Membrane_triangle_list[Membrane_num_of_Triangles][3], int Membrane_Normal_direction[Membrane_num_of_Triangles][2], vector<int> &ECM_Membrane_Trinagle_neighbours_3, double Cellcom[3], int Outer_Membrane_num_of_triangles, int  Outer_Membrane_num_of_Nodes, int Outer_Membrane_list_of_Nodes[], vector<int>  &Membrane_Edge_triangle_and_ECM_neighbours, double &Total_Potential_Energy)
{
    
    double ECM_triangle_AB[3], ECM_triangle_AC[3], ECM_triangle_ABxAC_unit_vector[3], ECM_triangle_BD[3], ECM_triangle_CD[3], ECM_triangle_AD[3], ECM_triangle_BCxBD[3], ECM_triangle_CAxCD[3], ECM_triangle_ABxAD[3], ECM_triangle_BC[3], ECM_triangle_CA[3], ECM_triangle_COM[3];
    
    double Membrane_triangle_AB[3], Membrane_triangle_AC[3], Membrane_triangle_ABxAC[3], Membrane_triangle_ABxAC_unit_vec[3]; // normal vector of membrane
    
    double temp_membrane_ecm_min_distance=1000;
    int temp_membraen_ecm_min_index;
    
    int temp_ECM_node_A, temp_ECM_node_B, temp_ECM_node_C;
    int temp_Membrane_Node_A, temp_Membrane_Node_B, temp_Membrane_Node_C;
    double Membrane_trinagle_com_shadow_on_ECM_triangle[3], Membrane_triangle_COM_ECM_distance_vec[3], membrane_tri_COM_ECM_distance_normal_component, Membrane_triangle_COM_ECM_distance_vec_length;
    //    int Membrane_ECM_trinagle_num_of_neighbours=0;
    
    double alpha_A, alpha_B, alpha_C, ECM_triangle_area;
    double reduced_radius, reduced_radius_squared, reduced_radius_cubed, reduced_radius_quint, temp_force_magnitude,f[3];
    double Membrane_triangle_COM[3]; // normal vector of membrane
    
    
    //F-net-com=0
    double Fnet__Active[3];  Fnet__Active[0]=0.0; Fnet__Active[1]=0.0; Fnet__Active[2]=0.0;
    double Fnet_Migration_force[3];  Fnet_Migration_force[0]=0.0; Fnet_Migration_force[1]=0.0; Fnet_Migration_force[2]=0.0;
    double random_number;
    
    //Dettachment step
    if(MD_step % Membrane_ECM_triangle_unbinding_update_step == 0 ) //The Membrane_node-ECM_triangle that interact are updated every 'Membrane_ECM_interaction_update_step' step. After update they will continue to initerct until the list is updated again.
    {
        
        for(int membrane_triangle_index=0; membrane_triangle_index<Outer_Membrane_num_of_triangles; membrane_triangle_index++)
        {
            //            for ( int ecm_neighbour_index=0; ecm_neighbour_index<ECM_Membrane_Trinagle_neighbours_2[membrane_triangle_index].size();
            //                 ecm_neighbour_index++)
            //            {
            temp_Membrane_Node_A=Membrane_triangle_list[membrane_triangle_index][0];
            temp_Membrane_Node_B=Membrane_triangle_list[membrane_triangle_index][1];
            temp_Membrane_Node_C=Membrane_triangle_list[membrane_triangle_index][2];
            
            Membrane_triangle_AB[0]=Membrane_Node_Position[ temp_Membrane_Node_B][0]-Membrane_Node_Position[ temp_Membrane_Node_A][0];
            Membrane_triangle_AB[1]=Membrane_Node_Position[ temp_Membrane_Node_B][1]-Membrane_Node_Position[ temp_Membrane_Node_A][1];
            Membrane_triangle_AB[2]=Membrane_Node_Position[ temp_Membrane_Node_B][2]-Membrane_Node_Position[ temp_Membrane_Node_A][2];
            Membrane_triangle_AC[0]=Membrane_Node_Position[ temp_Membrane_Node_C][0]-Membrane_Node_Position[ temp_Membrane_Node_A][0];
            Membrane_triangle_AC[1]=Membrane_Node_Position[ temp_Membrane_Node_C][1]-Membrane_Node_Position[ temp_Membrane_Node_A][1];
            Membrane_triangle_AC[2]=Membrane_Node_Position[ temp_Membrane_Node_C][2]-Membrane_Node_Position[ temp_Membrane_Node_A][2];
            
            crossvector(Membrane_triangle_ABxAC, Membrane_triangle_AB, Membrane_triangle_AC);
            
            Membrane_triangle_ABxAC[0]*=Membrane_Normal_direction[membrane_triangle_index][1];
            Membrane_triangle_ABxAC[1]*=Membrane_Normal_direction[membrane_triangle_index][1];
            Membrane_triangle_ABxAC[2]*=Membrane_Normal_direction[membrane_triangle_index][1];
            
            
            Membrane_triangle_COM[0]= ( Membrane_Node_Position[ temp_Membrane_Node_A][0] + Membrane_Node_Position[ temp_Membrane_Node_B][0] +Membrane_Node_Position[ temp_Membrane_Node_C][0] ) /3.0;
            Membrane_triangle_COM[1]= ( Membrane_Node_Position[ temp_Membrane_Node_A][1] +Membrane_Node_Position[ temp_Membrane_Node_B][1] +Membrane_Node_Position[ temp_Membrane_Node_C][1] ) /3.0;
            Membrane_triangle_COM[2]= ( Membrane_Node_Position[ temp_Membrane_Node_A][2] +Membrane_Node_Position[ temp_Membrane_Node_B][2] +Membrane_Node_Position[ temp_Membrane_Node_C][2] ) /3.0;
            
            Membrane_triangle_ABxAC_unit_vec[0]=Membrane_triangle_ABxAC[0]/vectorlength(Membrane_triangle_ABxAC);
            Membrane_triangle_ABxAC_unit_vec[1]=Membrane_triangle_ABxAC[1]/vectorlength(Membrane_triangle_ABxAC);
            Membrane_triangle_ABxAC_unit_vec[2]=Membrane_triangle_ABxAC[2]/vectorlength(Membrane_triangle_ABxAC);
            
            temp_ECM_node_A=ECM_surface_triangle_list[ECM_Membrane_Trinagle_neighbours_3[membrane_triangle_index]][0];
            temp_ECM_node_B=ECM_surface_triangle_list[ECM_Membrane_Trinagle_neighbours_3[membrane_triangle_index]][1];
            temp_ECM_node_C=ECM_surface_triangle_list[ECM_Membrane_Trinagle_neighbours_3[membrane_triangle_index]][2];
            
            ECM_triangle_AB[0]=ECM_Node_Position[temp_ECM_node_B][0]-ECM_Node_Position[temp_ECM_node_A][0];
            ECM_triangle_AB[1]=ECM_Node_Position[temp_ECM_node_B][1]-ECM_Node_Position[temp_ECM_node_A][1];
            ECM_triangle_AB[2]=ECM_Node_Position[temp_ECM_node_B][2]-ECM_Node_Position[temp_ECM_node_A][2];
            ECM_triangle_AC[0]=ECM_Node_Position[temp_ECM_node_C][0]-ECM_Node_Position[temp_ECM_node_A][0];
            ECM_triangle_AC[1]=ECM_Node_Position[temp_ECM_node_C][1]-ECM_Node_Position[temp_ECM_node_A][1];
            ECM_triangle_AC[2]=ECM_Node_Position[temp_ECM_node_C][2]-ECM_Node_Position[temp_ECM_node_A][2];
            
            crossvector(ECM_triangle_ABxAC_unit_vector, ECM_triangle_AB, ECM_triangle_AC);
            
            ECM_triangle_ABxAC_unit_vector[0]/=vectorlength(ECM_triangle_ABxAC_unit_vector);
            ECM_triangle_ABxAC_unit_vector[1]/=vectorlength(ECM_triangle_ABxAC_unit_vector);
            ECM_triangle_ABxAC_unit_vector[2]/=vectorlength(ECM_triangle_ABxAC_unit_vector);
            
            ECM_triangle_COM[0]= ( ECM_Node_Position[ temp_ECM_node_A][0] + ECM_Node_Position[ temp_ECM_node_B][0] + ECM_Node_Position[ temp_ECM_node_C][0] ) /3.0;
            ECM_triangle_COM[1]= ( ECM_Node_Position[ temp_ECM_node_A][1] + ECM_Node_Position[ temp_ECM_node_B][1] + ECM_Node_Position[ temp_ECM_node_C][1] ) /3.0;
            ECM_triangle_COM[2]= ( ECM_Node_Position[ temp_ECM_node_A][2] + ECM_Node_Position[ temp_ECM_node_B][2] + ECM_Node_Position[ temp_ECM_node_C][2] ) /3.0;
            
            Membrane_triangle_COM_ECM_distance_vec[0]= Membrane_triangle_COM[0] - ECM_triangle_COM[0];
            Membrane_triangle_COM_ECM_distance_vec[1]= Membrane_triangle_COM[1] - ECM_triangle_COM[1];
            Membrane_triangle_COM_ECM_distance_vec[2]= Membrane_triangle_COM[2] - ECM_triangle_COM[2];
            
            Membrane_triangle_COM_ECM_distance_vec_length=vectorlength(Membrane_triangle_COM_ECM_distance_vec);
            
            if ( (Membrane_triangle_COM_ECM_distance_vec_length > 2.5*epsilonECM)
                || (innerproduct(ECM_triangle_ABxAC_unit_vector, Membrane_triangle_ABxAC_unit_vec) > 0 )
                )
            {
                ECM_Membrane_Trinagle_neighbours_3[membrane_triangle_index]=-1;
                //                ecm_neighbour_index--;
                break;
            } else {
                random_number = ((double) rand() / (RAND_MAX));
                if (random_number < Membrane_ECM_UnBinding_prob) {
                    ECM_Membrane_Trinagle_neighbours_3[membrane_triangle_index]=-1;
                    //                    ecm_neighbour_index--;
                }
            }
            
            //            }//END OF: for ( int j=0;j<ECM_Surface_num_of_Triangles;j++)
        } //END OF: for(int membrane_triangle_index=0; membrane_triangle_index<Outer_Membrane_num_of_triangles; membrane_triangle_index++)
        //    exit (EXIT_FAILURE);
        
    } //END OF: if(istep % Membrane_ECM_interaction_update_step == 0 )
    //    int testing_counter;
    //Attachment step
    if(MD_step % Membrane_ECM_triangle_interaction_update_step == 0 ) //The Membrane_node-ECM_triangle that interact are updated every 'Membrane_ECM_interaction_update_step' step. After update they will continue to initerct until the list is updated again.
        
    {
        //streatching-begin
        if (spreading_flag==1.0) {
            Membrane_Edge_triangle_and_ECM_neighbours.clear();
            Membrane_Edge_triangle_and_ECM_neighbours.resize(Outer_Membrane_num_of_triangles);
            for (int i=0; i<Outer_Membrane_num_of_triangles; i++) {
                Membrane_Edge_triangle_and_ECM_neighbours[i]=-1;
            }
        }
        
        //        int spreading_triangle_count=0;
        //streatching-end
        //        testing_counter=0;
        for(int membrane_triangle_index=0; membrane_triangle_index<Outer_Membrane_num_of_triangles; membrane_triangle_index++)
        {
            temp_membrane_ecm_min_distance=1000;
            bool attachmed_list_flag=false;
            for ( int ecm_triangle_index=0; ecm_triangle_index<ECM_Surface_num_of_Triangles; ecm_triangle_index++)
            {
                for (int neighbour_ecm_index=0; neighbour_ecm_index < Outer_Membrane_num_of_triangles; neighbour_ecm_index++)
                {
                    if (ECM_Membrane_Trinagle_neighbours_3[neighbour_ecm_index] == ecm_triangle_index) {
                        attachmed_list_flag = true;
                        break;
                    }
                }
                if (attachmed_list_flag != true)
                {
                    temp_Membrane_Node_A=Membrane_triangle_list[membrane_triangle_index][0];
                    temp_Membrane_Node_B=Membrane_triangle_list[membrane_triangle_index][1];
                    temp_Membrane_Node_C=Membrane_triangle_list[membrane_triangle_index][2];
                    
                    Membrane_triangle_AB[0]=Membrane_Node_Position[ temp_Membrane_Node_B][0]-Membrane_Node_Position[ temp_Membrane_Node_A][0];
                    Membrane_triangle_AB[1]=Membrane_Node_Position[ temp_Membrane_Node_B][1]-Membrane_Node_Position[ temp_Membrane_Node_A][1];
                    Membrane_triangle_AB[2]=Membrane_Node_Position[ temp_Membrane_Node_B][2]-Membrane_Node_Position[ temp_Membrane_Node_A][2];
                    Membrane_triangle_AC[0]=Membrane_Node_Position[ temp_Membrane_Node_C][0]-Membrane_Node_Position[ temp_Membrane_Node_A][0];
                    Membrane_triangle_AC[1]=Membrane_Node_Position[ temp_Membrane_Node_C][1]-Membrane_Node_Position[ temp_Membrane_Node_A][1];
                    Membrane_triangle_AC[2]=Membrane_Node_Position[ temp_Membrane_Node_C][2]-Membrane_Node_Position[ temp_Membrane_Node_A][2];
                    
                    crossvector(Membrane_triangle_ABxAC, Membrane_triangle_AB, Membrane_triangle_AC);
                    
                    Membrane_triangle_ABxAC[0]*=Membrane_Normal_direction[membrane_triangle_index][1];
                    Membrane_triangle_ABxAC[1]*=Membrane_Normal_direction[membrane_triangle_index][1];
                    Membrane_triangle_ABxAC[2]*=Membrane_Normal_direction[membrane_triangle_index][1];
                    
                    
                    Membrane_triangle_COM[0]= ( Membrane_Node_Position[ temp_Membrane_Node_A][0] + Membrane_Node_Position[ temp_Membrane_Node_B][0] +Membrane_Node_Position[ temp_Membrane_Node_C][0] ) /3.0;
                    Membrane_triangle_COM[1]= ( Membrane_Node_Position[ temp_Membrane_Node_A][1] +Membrane_Node_Position[ temp_Membrane_Node_B][1] +Membrane_Node_Position[ temp_Membrane_Node_C][1] ) /3.0;
                    Membrane_triangle_COM[2]= ( Membrane_Node_Position[ temp_Membrane_Node_A][2] +Membrane_Node_Position[ temp_Membrane_Node_B][2] +Membrane_Node_Position[ temp_Membrane_Node_C][2] ) /3.0;
                    
                    Membrane_triangle_ABxAC_unit_vec[0]=Membrane_triangle_ABxAC[0]/vectorlength(Membrane_triangle_ABxAC);
                    Membrane_triangle_ABxAC_unit_vec[1]=Membrane_triangle_ABxAC[1]/vectorlength(Membrane_triangle_ABxAC);
                    Membrane_triangle_ABxAC_unit_vec[2]=Membrane_triangle_ABxAC[2]/vectorlength(Membrane_triangle_ABxAC);
                    
                    temp_ECM_node_A=ECM_surface_triangle_list[ecm_triangle_index][0];
                    temp_ECM_node_B=ECM_surface_triangle_list[ecm_triangle_index][1];
                    temp_ECM_node_C=ECM_surface_triangle_list[ecm_triangle_index][2];
                    
                    ECM_triangle_AB[0]=ECM_Node_Position[temp_ECM_node_B][0]-ECM_Node_Position[temp_ECM_node_A][0];
                    ECM_triangle_AB[1]=ECM_Node_Position[temp_ECM_node_B][1]-ECM_Node_Position[temp_ECM_node_A][1];
                    ECM_triangle_AB[2]=ECM_Node_Position[temp_ECM_node_B][2]-ECM_Node_Position[temp_ECM_node_A][2];
                    ECM_triangle_AC[0]=ECM_Node_Position[temp_ECM_node_C][0]-ECM_Node_Position[temp_ECM_node_A][0];
                    ECM_triangle_AC[1]=ECM_Node_Position[temp_ECM_node_C][1]-ECM_Node_Position[temp_ECM_node_A][1];
                    ECM_triangle_AC[2]=ECM_Node_Position[temp_ECM_node_C][2]-ECM_Node_Position[temp_ECM_node_A][2];
                    
                    crossvector(ECM_triangle_ABxAC_unit_vector, ECM_triangle_AB, ECM_triangle_AC);
                    
                    ECM_triangle_ABxAC_unit_vector[0]/=vectorlength(ECM_triangle_ABxAC_unit_vector);
                    ECM_triangle_ABxAC_unit_vector[1]/=vectorlength(ECM_triangle_ABxAC_unit_vector);
                    ECM_triangle_ABxAC_unit_vector[2]/=vectorlength(ECM_triangle_ABxAC_unit_vector);
                    
                    ECM_triangle_COM[0]= ( ECM_Node_Position[ temp_ECM_node_A][0] + ECM_Node_Position[ temp_ECM_node_B][0] + ECM_Node_Position[ temp_ECM_node_C][0] ) /3.0;
                    ECM_triangle_COM[1]= ( ECM_Node_Position[ temp_ECM_node_A][1] + ECM_Node_Position[ temp_ECM_node_B][1] + ECM_Node_Position[ temp_ECM_node_C][1] ) /3.0;
                    ECM_triangle_COM[2]= ( ECM_Node_Position[ temp_ECM_node_A][2] + ECM_Node_Position[ temp_ECM_node_B][2] + ECM_Node_Position[ temp_ECM_node_C][2] ) /3.0;
                    
                    Membrane_triangle_COM_ECM_distance_vec[0]= Membrane_triangle_COM[0] - ECM_triangle_COM[0];
                    Membrane_triangle_COM_ECM_distance_vec[1]= Membrane_triangle_COM[1] - ECM_triangle_COM[1];
                    Membrane_triangle_COM_ECM_distance_vec[2]= Membrane_triangle_COM[2] - ECM_triangle_COM[2];
                    
                    membrane_tri_COM_ECM_distance_normal_component=innerproduct(Membrane_triangle_COM_ECM_distance_vec,ECM_triangle_ABxAC_unit_vector);
                    Membrane_triangle_COM_ECM_distance_vec_length=vectorlength(Membrane_triangle_COM_ECM_distance_vec);
                    
                    if ( (Membrane_triangle_COM_ECM_distance_vec_length < spreading_force_max_range*epsilonECM) &&
                        (innerproduct(ECM_triangle_ABxAC_unit_vector, Membrane_triangle_ABxAC_unit_vec) < spreading_force_cos_triangle_interaction_angle))
                    {
                        Membrane_trinagle_com_shadow_on_ECM_triangle[0]=Membrane_triangle_COM[0] - membrane_tri_COM_ECM_distance_normal_component * ECM_triangle_ABxAC_unit_vector[0];
                        Membrane_trinagle_com_shadow_on_ECM_triangle[1]=Membrane_triangle_COM[1] - membrane_tri_COM_ECM_distance_normal_component * ECM_triangle_ABxAC_unit_vector[1];
                        Membrane_trinagle_com_shadow_on_ECM_triangle[2]=Membrane_triangle_COM[2] - membrane_tri_COM_ECM_distance_normal_component * ECM_triangle_ABxAC_unit_vector[2];
                        
                        ECM_triangle_BD[0]=Membrane_trinagle_com_shadow_on_ECM_triangle[0]-ECM_Node_Position[temp_ECM_node_B][0];
                        ECM_triangle_BD[1]=Membrane_trinagle_com_shadow_on_ECM_triangle[1]-ECM_Node_Position[temp_ECM_node_B][1];
                        ECM_triangle_BD[2]=Membrane_trinagle_com_shadow_on_ECM_triangle[2]-ECM_Node_Position[temp_ECM_node_B][2];
                        
                        ECM_triangle_CD[0]=Membrane_trinagle_com_shadow_on_ECM_triangle[0]-ECM_Node_Position[temp_ECM_node_C][0];
                        ECM_triangle_CD[1]=Membrane_trinagle_com_shadow_on_ECM_triangle[1]-ECM_Node_Position[temp_ECM_node_C][1];
                        ECM_triangle_CD[2]=Membrane_trinagle_com_shadow_on_ECM_triangle[2]-ECM_Node_Position[temp_ECM_node_C][2];
                        
                        ECM_triangle_AD[0]=Membrane_trinagle_com_shadow_on_ECM_triangle[0]-ECM_Node_Position[temp_ECM_node_A][0];
                        ECM_triangle_AD[1]=Membrane_trinagle_com_shadow_on_ECM_triangle[1]-ECM_Node_Position[temp_ECM_node_A][1];
                        ECM_triangle_AD[2]=Membrane_trinagle_com_shadow_on_ECM_triangle[2]-ECM_Node_Position[temp_ECM_node_A][2];
                        
                        ECM_triangle_BC[0]=ECM_Node_Position[temp_ECM_node_C][0]-ECM_Node_Position[temp_ECM_node_B][0];
                        ECM_triangle_BC[1]=ECM_Node_Position[temp_ECM_node_C][1]-ECM_Node_Position[temp_ECM_node_B][1];
                        ECM_triangle_BC[2]=ECM_Node_Position[temp_ECM_node_C][2]-ECM_Node_Position[temp_ECM_node_B][2];
                        
                        ECM_triangle_CA[0]=-ECM_triangle_AC[0];
                        ECM_triangle_CA[1]=-ECM_triangle_AC[1];
                        ECM_triangle_CA[2]=-ECM_triangle_AC[2];
                        
                        crossvector(ECM_triangle_BCxBD, ECM_triangle_BC, ECM_triangle_BD);
                        alpha_A=0.5*innerproduct(ECM_triangle_ABxAC_unit_vector, ECM_triangle_BCxBD);
                        
                        crossvector(ECM_triangle_CAxCD, ECM_triangle_CA, ECM_triangle_CD);
                        alpha_B=0.5*innerproduct(ECM_triangle_ABxAC_unit_vector, ECM_triangle_CAxCD);
                        
                        crossvector(ECM_triangle_ABxAD, ECM_triangle_AB, ECM_triangle_AD);
                        alpha_C=0.5*innerproduct(ECM_triangle_ABxAC_unit_vector, ECM_triangle_ABxAD);
                        
                        if(alpha_A>=0 & alpha_B>=0 & alpha_C>=0)
                        {
                            if ((Membrane_triangle_COM_ECM_distance_vec_length < 2.5*epsilonECM)
                                &&
                                (innerproduct(ECM_triangle_ABxAC_unit_vector, Membrane_triangle_ABxAC_unit_vec) < -Membrane_ECM_Max_cos_triangle_interaction_angle)
                                &&
                                (Membrane_triangle_COM_ECM_distance_vec_length < temp_membrane_ecm_min_distance)) {
                                random_number = ((double) rand() / (RAND_MAX));
                                if (random_number < Membrane_ECM_Binding_prob) {
                                    ECM_Membrane_Trinagle_neighbours_3[membrane_triangle_index]=ecm_triangle_index;
                                    temp_membrane_ecm_min_distance=Membrane_triangle_COM_ECM_distance_vec_length;
                                }
                            }
                            if ((Membrane_triangle_COM_ECM_distance_vec_length > spreading_force_min_range*epsilonECM) && (innerproduct(ECM_triangle_ABxAC_unit_vector, Membrane_triangle_ABxAC_unit_vec) > -spreading_force_cos_triangle_interaction_angle) && spreading_flag==1.0) {
                                Membrane_Edge_triangle_and_ECM_neighbours[membrane_triangle_index]=ecm_triangle_index;
                                //                                spreading_triangle_count++;
                            }
                        }
                    } //END OF: if ( abs(membrane_tri_COM_ECM_distance_normal_component) && (innerproduct(ECM_triangle_ABxAC_unit_vector, Membrane_triangle_ABxAC_unit_vec) < -Membrane_ECM_Max_cos_triangle_interaction_angle ))
                }//END OF: if (attachment_list_flag != true)
            }//END OF: for ( int j=0;j<ECM_Surface_num_of_Triangles;j++)
        } //END OF: for(int membrane_triangle_index=0; membrane_triangle_index<Outer_Membrane_num_of_triangles; membrane_triangle_index++)
        //    exit (EXIT_FAILURE);
        //        cout<<"spreading_triangle_count:  "<<spreading_triangle_count<<endl;
    } //END OF: if(istep % Membrane_ECM_interaction_update_step == 0 )
    
    //Force calculation
    for(int temp_trinagle_index=0; temp_trinagle_index<Outer_Membrane_num_of_triangles; temp_trinagle_index++)
    {
        
        if (ECM_Membrane_Trinagle_neighbours_3[temp_trinagle_index]!=-1) {
            
            temp_Membrane_Node_A=Membrane_triangle_list[temp_trinagle_index][0];
            temp_Membrane_Node_B=Membrane_triangle_list[temp_trinagle_index][1];
            temp_Membrane_Node_C=Membrane_triangle_list[temp_trinagle_index][2];
            
            Membrane_triangle_AB[0]=Membrane_Node_Position[ temp_Membrane_Node_B][0]-Membrane_Node_Position[ temp_Membrane_Node_A][0];
            Membrane_triangle_AB[1]=Membrane_Node_Position[ temp_Membrane_Node_B][1]-Membrane_Node_Position[ temp_Membrane_Node_A][1];
            Membrane_triangle_AB[2]=Membrane_Node_Position[ temp_Membrane_Node_B][2]-Membrane_Node_Position[ temp_Membrane_Node_A][2];
            
            Membrane_triangle_AC[0]=Membrane_Node_Position[ temp_Membrane_Node_C][0]-Membrane_Node_Position[ temp_Membrane_Node_A][0];
            Membrane_triangle_AC[1]=Membrane_Node_Position[ temp_Membrane_Node_C][1]-Membrane_Node_Position[ temp_Membrane_Node_A][1];
            Membrane_triangle_AC[2]=Membrane_Node_Position[ temp_Membrane_Node_C][2]-Membrane_Node_Position[ temp_Membrane_Node_A][2];
            
            crossvector(Membrane_triangle_ABxAC, Membrane_triangle_AB, Membrane_triangle_AC);
            
            Membrane_triangle_ABxAC_unit_vec[0]=Membrane_triangle_ABxAC[0]*Membrane_Normal_direction[temp_trinagle_index][1]/vectorlength(Membrane_triangle_ABxAC);
            Membrane_triangle_ABxAC_unit_vec[1]=Membrane_triangle_ABxAC[1]*Membrane_Normal_direction[temp_trinagle_index][1]/vectorlength(Membrane_triangle_ABxAC);
            Membrane_triangle_ABxAC_unit_vec[2]=Membrane_triangle_ABxAC[2]*Membrane_Normal_direction[temp_trinagle_index][1]/vectorlength(Membrane_triangle_ABxAC);
            
            Membrane_triangle_COM[0]= ( Membrane_Node_Position[ temp_Membrane_Node_A][0] + Membrane_Node_Position[ temp_Membrane_Node_B][0] + Membrane_Node_Position[ temp_Membrane_Node_C][0]) /3.0;
            Membrane_triangle_COM[1]= ( Membrane_Node_Position[ temp_Membrane_Node_A][1] + Membrane_Node_Position[ temp_Membrane_Node_B][1] + Membrane_Node_Position[ temp_Membrane_Node_C][1]) /3.0;
            Membrane_triangle_COM[2]= ( Membrane_Node_Position[ temp_Membrane_Node_A][2] +Membrane_Node_Position[ temp_Membrane_Node_B][2] + Membrane_Node_Position[ temp_Membrane_Node_C][2]) /3.0;
            //            for (int neighbour_counter=0; neighbour_counter<ECM_Membrane_Trinagle_neighbours_2[temp_trinagle_index].size(); neighbour_counter++)
            //            {
            temp_ECM_node_A=ECM_surface_triangle_list[ECM_Membrane_Trinagle_neighbours_3[temp_trinagle_index]][0];
            temp_ECM_node_B=ECM_surface_triangle_list[ECM_Membrane_Trinagle_neighbours_3[temp_trinagle_index]][1];
            temp_ECM_node_C=ECM_surface_triangle_list[ECM_Membrane_Trinagle_neighbours_3[temp_trinagle_index]][2];
            
            ECM_triangle_AB[0]=ECM_Node_Position[temp_ECM_node_B][0]-ECM_Node_Position[temp_ECM_node_A][0];
            ECM_triangle_AB[1]=ECM_Node_Position[temp_ECM_node_B][1]-ECM_Node_Position[temp_ECM_node_A][1];
            ECM_triangle_AB[2]=ECM_Node_Position[temp_ECM_node_B][2]-ECM_Node_Position[temp_ECM_node_A][2];
            ECM_triangle_AC[0]=ECM_Node_Position[temp_ECM_node_C][0]-ECM_Node_Position[temp_ECM_node_A][0];
            ECM_triangle_AC[1]=ECM_Node_Position[temp_ECM_node_C][1]-ECM_Node_Position[temp_ECM_node_A][1];
            ECM_triangle_AC[2]=ECM_Node_Position[temp_ECM_node_C][2]-ECM_Node_Position[temp_ECM_node_A][2];
            
            crossvector(ECM_triangle_ABxAC_unit_vector, ECM_triangle_AB, ECM_triangle_AC);
            
            ECM_triangle_ABxAC_unit_vector[0]=ECM_triangle_ABxAC_unit_vector[0]/vectorlength(ECM_triangle_ABxAC_unit_vector);
            ECM_triangle_ABxAC_unit_vector[1]=ECM_triangle_ABxAC_unit_vector[1]/vectorlength(ECM_triangle_ABxAC_unit_vector);
            ECM_triangle_ABxAC_unit_vector[2]=ECM_triangle_ABxAC_unit_vector[2]/vectorlength(ECM_triangle_ABxAC_unit_vector);
            
            ECM_triangle_COM[0]= ( ECM_Node_Position[ temp_ECM_node_A][0] + ECM_Node_Position[ temp_ECM_node_B][0] + ECM_Node_Position[ temp_ECM_node_C][0] ) /3.0;
            ECM_triangle_COM[1]= ( ECM_Node_Position[ temp_ECM_node_A][1] + ECM_Node_Position[ temp_ECM_node_B][1] + ECM_Node_Position[ temp_ECM_node_C][1] ) /3.0;
            ECM_triangle_COM[2]= ( ECM_Node_Position[ temp_ECM_node_A][2] + ECM_Node_Position[ temp_ECM_node_B][2] + ECM_Node_Position[ temp_ECM_node_C][2] ) /3.0;
            
            Membrane_triangle_COM_ECM_distance_vec[0]= Membrane_triangle_COM[0] - ECM_triangle_COM[0];
            Membrane_triangle_COM_ECM_distance_vec[1]= Membrane_triangle_COM[1] - ECM_triangle_COM[1];
            Membrane_triangle_COM_ECM_distance_vec[2]= Membrane_triangle_COM[2] - ECM_triangle_COM[2];
            
            membrane_tri_COM_ECM_distance_normal_component=innerproduct(Membrane_triangle_COM_ECM_distance_vec, ECM_triangle_ABxAC_unit_vector);
            
            Membrane_trinagle_com_shadow_on_ECM_triangle[0]=Membrane_triangle_COM[0] - membrane_tri_COM_ECM_distance_normal_component * ECM_triangle_ABxAC_unit_vector[0];
            Membrane_trinagle_com_shadow_on_ECM_triangle[1]=Membrane_triangle_COM[1] - membrane_tri_COM_ECM_distance_normal_component * ECM_triangle_ABxAC_unit_vector[1];
            Membrane_trinagle_com_shadow_on_ECM_triangle[2]=Membrane_triangle_COM[2] - membrane_tri_COM_ECM_distance_normal_component * ECM_triangle_ABxAC_unit_vector[2];
            
            Membrane_triangle_COM_ECM_distance_vec_length=vectorlength(Membrane_triangle_COM_ECM_distance_vec);
            
            ECM_triangle_BD[0]=Membrane_trinagle_com_shadow_on_ECM_triangle[0]-ECM_Node_Position[temp_ECM_node_B][0];
            ECM_triangle_BD[1]=Membrane_trinagle_com_shadow_on_ECM_triangle[1]-ECM_Node_Position[temp_ECM_node_B][1];
            ECM_triangle_BD[2]=Membrane_trinagle_com_shadow_on_ECM_triangle[2]-ECM_Node_Position[temp_ECM_node_B][2];
            
            ECM_triangle_CD[0]=Membrane_trinagle_com_shadow_on_ECM_triangle[0]-ECM_Node_Position[temp_ECM_node_C][0];
            ECM_triangle_CD[1]=Membrane_trinagle_com_shadow_on_ECM_triangle[1]-ECM_Node_Position[temp_ECM_node_C][1];
            ECM_triangle_CD[2]=Membrane_trinagle_com_shadow_on_ECM_triangle[2]-ECM_Node_Position[temp_ECM_node_C][2];
            
            ECM_triangle_AD[0]=Membrane_trinagle_com_shadow_on_ECM_triangle[0]-ECM_Node_Position[temp_ECM_node_A][0];
            ECM_triangle_AD[1]=Membrane_trinagle_com_shadow_on_ECM_triangle[1]-ECM_Node_Position[temp_ECM_node_A][1];
            ECM_triangle_AD[2]=Membrane_trinagle_com_shadow_on_ECM_triangle[2]-ECM_Node_Position[temp_ECM_node_A][2];
            
            ECM_triangle_BC[0]=ECM_Node_Position[temp_ECM_node_C][0]-ECM_Node_Position[temp_ECM_node_B][0];
            ECM_triangle_BC[1]=ECM_Node_Position[temp_ECM_node_C][1]-ECM_Node_Position[temp_ECM_node_B][1];
            ECM_triangle_BC[2]=ECM_Node_Position[temp_ECM_node_C][2]-ECM_Node_Position[temp_ECM_node_B][2];
            
            ECM_triangle_CA[0]=-ECM_triangle_AC[0];
            ECM_triangle_CA[1]=-ECM_triangle_AC[1];
            ECM_triangle_CA[2]=-ECM_triangle_AC[2];
            
            crossvector(ECM_triangle_BCxBD, ECM_triangle_BC, ECM_triangle_BD);
            alpha_A=0.5*innerproduct(ECM_triangle_ABxAC_unit_vector, ECM_triangle_BCxBD);
            
            crossvector(ECM_triangle_CAxCD, ECM_triangle_CA, ECM_triangle_CD);
            alpha_B=0.5*innerproduct(ECM_triangle_ABxAC_unit_vector, ECM_triangle_CAxCD);
            
            crossvector(ECM_triangle_ABxAD, ECM_triangle_AB, ECM_triangle_AD);
            alpha_C=0.5*innerproduct(ECM_triangle_ABxAC_unit_vector, ECM_triangle_ABxAD);
            
            //            =============================================================
            ECM_triangle_area = alpha_A + alpha_B + alpha_C;
            
            reduced_radius = (1.5*epsilonECM/Membrane_triangle_COM_ECM_distance_vec_length);
            reduced_radius_squared = reduced_radius*reduced_radius;
            reduced_radius_cubed = reduced_radius_squared*reduced_radius;
            reduced_radius_quint = reduced_radius_cubed*reduced_radius_squared;
            
            temp_force_magnitude = 4.0*Membrane_ECM_interaction_strength*( reduced_radius_quint - reduced_radius_cubed )/(1.5*epsilonECM*Membrane_triangle_COM_ECM_distance_vec_length);
            
            if (energy_calculation_flag==1.0) {
                Total_Potential_Energy+=Membrane_ECM_interaction_strength*(reduced_radius_squared*reduced_radius_squared-2*reduced_radius_squared);
            }
            //                f[0] = temp_force_magnitude * Membrane_triangle_COM_ECM_distance_vec[0];
            f[1] = temp_force_magnitude * Membrane_triangle_COM_ECM_distance_vec[1];
            //                f[2] = temp_force_magnitude * Membrane_triangle_COM_ECM_distance_vec[2];
            f[0]=0.0;
            f[2]=0.0;
            
            for (int abc_counter=0; abc_counter<3; abc_counter++) {
                Membrane_Node_Force[temp_Membrane_Node_A][abc_counter] += -f[abc_counter]/3.0;
                Membrane_Node_Force[temp_Membrane_Node_B][abc_counter] += -f[abc_counter]/3.0;
                Membrane_Node_Force[temp_Membrane_Node_C][abc_counter] += -f[abc_counter]/3.0;
                
                ECM_Node_Force[temp_ECM_node_A][abc_counter] += (alpha_A/ECM_triangle_area)*f[abc_counter];
                ECM_Node_Force[temp_ECM_node_B][abc_counter] += (alpha_B/ECM_triangle_area)*f[abc_counter];
                ECM_Node_Force[temp_ECM_node_C][abc_counter] += (alpha_C/ECM_triangle_area)*f[abc_counter];
                
            }
            //            }//END OF: for (int neighbour_counter=0; neighbour_counter<ECM_Membrane_Trinagle_neighbours_2[i].size(); neighbour_counter++)
            
            //            spreading force-begin
            if (Membrane_Edge_triangle_and_ECM_neighbours[temp_trinagle_index]!=-1 && spreading_flag==1.0) {
                int temp_ecm_spreading_index=Membrane_Edge_triangle_and_ECM_neighbours[temp_trinagle_index];
                temp_ECM_node_A=ECM_surface_triangle_list[temp_ecm_spreading_index][0];
                temp_ECM_node_B=ECM_surface_triangle_list[temp_ecm_spreading_index][1];
                temp_ECM_node_C=ECM_surface_triangle_list[temp_ecm_spreading_index][2];
                
                f[0] = spreading_force_magnitude * Membrane_triangle_COM_ECM_distance_vec[0];
                f[1] = spreading_force_magnitude * Membrane_triangle_COM_ECM_distance_vec[1];
                f[2] = spreading_force_magnitude * Membrane_triangle_COM_ECM_distance_vec[2];
                
                for (int abc_counter=0; abc_counter<3; abc_counter++) {
                    Membrane_Node_Force[temp_Membrane_Node_A][abc_counter] += -f[abc_counter]/3.0;
                    Membrane_Node_Force[temp_Membrane_Node_B][abc_counter] += -f[abc_counter]/3.0;
                    Membrane_Node_Force[temp_Membrane_Node_C][abc_counter] += -f[abc_counter]/3.0;
                    
                    ECM_Node_Force[temp_ECM_node_A][abc_counter] += f[abc_counter]/3.0;
                    ECM_Node_Force[temp_ECM_node_B][abc_counter] += f[abc_counter]/3.0;
                    ECM_Node_Force[temp_ECM_node_C][abc_counter] += f[abc_counter]/3.0;
                    
                }
            }
            //            spreading force-end
            
            
        }// END OF: if (ECM_Membrane_Trinagle_neighbours[i]!=ECM_Surface_num_of_Triangles)
    }//END OF: for(int i=0;i<Outer_Membrane_num_of_triangles;i++)
    
}

void Membrane_ECM_interactionHELPER(double  ECM_Node_Position [][3],int indexA,int indexB,int indexC,double n[3],double D[3],double (&a1),double (&a2),double (&a3))
{
    double AB[3],AD[3],BC[3],BD[3],CA[3],CD[3];
    double temp[3];
    
    double xECMindexA[3],xECMindexC[3],xECMindexB[3],xcom[3];
    double AC[3],ACL,ABL,H[3],HL,H2L;
    
    ///solving the problem of edges
    
    xcom[0]=(ECM_Node_Position[indexA][0] + ECM_Node_Position[indexB][0] + ECM_Node_Position[indexC][0]  )/ 3.0 ;
    xcom[1]=(ECM_Node_Position[indexA][1] + ECM_Node_Position[indexB][1] + ECM_Node_Position[indexC][1]  )/ 3.0 ;
    xcom[2]=(ECM_Node_Position[indexA][2] + ECM_Node_Position[indexB][2] + ECM_Node_Position[indexC][2]  )/ 3.0 ;
    
    //--------------------------------------- indexA
    AC[0]= ECM_Node_Position[indexA][0] -xcom[0]  ;
    AC[1]= ECM_Node_Position[indexA][1] -xcom[1]  ;
    AC[2]= ECM_Node_Position[indexA][2] -xcom[2]  ;
    ACL=vectorlength(AC);
    
    AB[0]= ECM_Node_Position[indexB][0] - ECM_Node_Position[indexA][0] ;
    AB[1]= ECM_Node_Position[indexB][1] - ECM_Node_Position[indexA][1] ;
    AB[2]= ECM_Node_Position[indexB][2] - ECM_Node_Position[indexA][2] ;
    ABL=vectorlength(AB);
    
    H[0]= AC[0]* sqrt(  1.0- innerproduct(AC,AB)* innerproduct(AC,AB) / ( ABL*ABL*ACL*ACL )    )  ;// h=AC*sin(teta)
    H[1]= AC[1]* sqrt(  1.0- innerproduct(AC,AB)* innerproduct(AC,AB) / ( ABL*ABL*ACL*ACL )    )  ;// h=AC*sin(teta
    H[2]= AC[2]* sqrt(  1.0- innerproduct(AC,AB)* innerproduct(AC,AB) / ( ABL*ABL*ACL*ACL )    )  ;// h=AC*sin(teta
    HL=vectorlength(H);
    
    H2L=HL+extengingtriangleECM*epsilonECM;
    
    xECMindexA[0]=xcom[0] + AC[0] * H2L/HL ;
    xECMindexA[1]=xcom[1] + AC[1] * H2L/HL ;
    xECMindexA[2]=xcom[2] + AC[2] * H2L/HL ;
    
    
    
    
    //--------------------------------------- indexB
    AC[0]= ECM_Node_Position[indexB][0] -xcom[0]  ;
    AC[1]= ECM_Node_Position[indexB][1] -xcom[1]  ;
    AC[2]= ECM_Node_Position[indexB][2] -xcom[2]  ;
    ACL=vectorlength(AC);
    
    AB[0]= ECM_Node_Position[indexC][0] - ECM_Node_Position[indexB][0] ;
    AB[1]= ECM_Node_Position[indexC][1] - ECM_Node_Position[indexB][1] ;
    AB[2]= ECM_Node_Position[indexC][2] - ECM_Node_Position[indexB][2] ;
    ABL=vectorlength(AB);
    
    H[0]= AC[0]* sqrt(  1.0- innerproduct(AC,AB)* innerproduct(AC,AB) / ( ABL*ABL*ACL*ACL )    )  ;// h=AC*sin(teta)
    H[1]= AC[1]* sqrt(  1.0- innerproduct(AC,AB)* innerproduct(AC,AB) / ( ABL*ABL*ACL*ACL )    )  ;// h=AC*sin(teta
    H[2]= AC[2]* sqrt(  1.0- innerproduct(AC,AB)* innerproduct(AC,AB) / ( ABL*ABL*ACL*ACL )    )  ;// h=AC*sin(teta
    HL=vectorlength(H);
    
    H2L=HL+extengingtriangleECM*epsilonECM;
    
    xECMindexB[0]=xcom[0] + AC[0] * H2L/HL ;
    xECMindexB[1]=xcom[1] + AC[1] * H2L/HL ;
    xECMindexB[2]=xcom[2] + AC[2] * H2L/HL ;
    
    
    
    
    //--------------------------------------- indexC
    AC[0]= ECM_Node_Position[indexC][0] -xcom[0]  ;
    AC[1]= ECM_Node_Position[indexC][1] -xcom[1]  ;
    AC[2]= ECM_Node_Position[indexC][2] -xcom[2]  ;
    ACL=vectorlength(AC);
    
    AB[0]= ECM_Node_Position[indexA][0] - ECM_Node_Position[indexC][0] ;
    AB[1]= ECM_Node_Position[indexA][1] - ECM_Node_Position[indexC][1] ;
    AB[2]= ECM_Node_Position[indexA][2] - ECM_Node_Position[indexC][2] ;
    ABL=vectorlength(AB);
    
    H[0]= AC[0]* sqrt(  1.0- innerproduct(AC,AB)* innerproduct(AC,AB) / ( ABL*ABL*ACL*ACL )    )  ;// h=AC*sin(teta)
    H[1]= AC[1]* sqrt(  1.0- innerproduct(AC,AB)* innerproduct(AC,AB) / ( ABL*ABL*ACL*ACL )    )  ;// h=AC*sin(teta
    H[2]= AC[2]* sqrt(  1.0- innerproduct(AC,AB)* innerproduct(AC,AB) / ( ABL*ABL*ACL*ACL )    )  ;// h=AC*sin(teta
    HL=vectorlength(H);
    
    H2L=HL+extengingtriangleECM*epsilonECM;
    
    xECMindexC[0]=xcom[0] + AC[0] * H2L/HL ;
    xECMindexC[1]=xcom[1] + AC[1] * H2L/HL ;
    xECMindexC[2]=xcom[2] + AC[2] * H2L/HL ;
    
    
    
    
    
    AB[0]=xECMindexB[0]-xECMindexA[0];
    AB[1]=xECMindexB[1]-xECMindexA[1];
    AB[2]=xECMindexB[2]-xECMindexA[2];
    
    AD[0]=D[0]-xECMindexA[0];
    AD[1]=D[1]-xECMindexA[1];
    AD[2]=D[2]-xECMindexA[2];
    
    BC[0]=xECMindexC[0]-xECMindexB[0];
    BC[1]=xECMindexC[1]-xECMindexB[1];
    BC[2]=xECMindexC[2]-xECMindexB[2];
    
    BD[0]=D[0]-xECMindexB[0];
    BD[1]=D[1]-xECMindexB[1];
    BD[2]=D[2]-xECMindexB[2];
    
    CA[0]=-xECMindexC[0]+xECMindexA[0];
    CA[1]=-xECMindexC[1]+xECMindexA[1];
    CA[2]=-xECMindexC[2]+xECMindexA[2];
    
    CD[0]=D[0]-xECMindexC[0];
    CD[1]=D[1]-xECMindexC[1];
    CD[2]=D[2]-xECMindexC[2];
    
    crossvector(temp,BC,BD);
    a1=innerproduct (n,temp);
    
    crossvector(temp,CA,CD);
    a2=innerproduct (n,temp);
    
    crossvector(temp,AB,AD);
    a3=innerproduct (n,temp);
    
    
    
    
    
    //cout<< a1<<"  "<<a2<<"  "<< a3<<endl;
    
    
}


void Membrane_ECM_interactionHELPER_newstrategy(double ECM_Node_Position [][3], int temp_ECM_node_A, int temp_ECM_node_B, int temp_ECM_node_C, double ECM_triangle_ABxAC_unit_vector[3], double D[3], double (&a1), double (&a2), double (&a3))
{
    double AB[3], AD[3], BC[3], BD[3], CA[3], CD[3];
    double temp[3];
    
    double xECMindexA[3],xECMindexC[3],xECMindexB[3], ECM_triangle_COM[3];
    double AC[3],ACL,ABL,H[3],HL,H2L;
    
    ///solving the problem of edges
    
    ECM_triangle_COM[0]=(ECM_Node_Position[temp_ECM_node_A][0] + ECM_Node_Position[temp_ECM_node_B][0] + ECM_Node_Position[temp_ECM_node_C][0]  )/ 3.0 ;
    ECM_triangle_COM[1]=(ECM_Node_Position[temp_ECM_node_A][1] + ECM_Node_Position[temp_ECM_node_B][1] + ECM_Node_Position[temp_ECM_node_C][1]  )/ 3.0 ;
    ECM_triangle_COM[2]=(ECM_Node_Position[temp_ECM_node_A][2] + ECM_Node_Position[temp_ECM_node_B][2] + ECM_Node_Position[temp_ECM_node_C][2]  )/ 3.0 ;
    
    //--------------------------------------- indexA
    AC[0]= ECM_Node_Position[temp_ECM_node_A][0] -ECM_triangle_COM[0]  ;
    AC[1]= ECM_Node_Position[temp_ECM_node_A][1] -ECM_triangle_COM[1]  ;
    AC[2]= ECM_Node_Position[temp_ECM_node_A][2] -ECM_triangle_COM[2]  ;
    ACL=vectorlength(AC);
    
    AB[0]= ECM_Node_Position[temp_ECM_node_B][0] - ECM_Node_Position[temp_ECM_node_A][0] ;
    AB[1]= ECM_Node_Position[temp_ECM_node_B][1] - ECM_Node_Position[temp_ECM_node_A][1] ;
    AB[2]= ECM_Node_Position[temp_ECM_node_B][2] - ECM_Node_Position[temp_ECM_node_A][2] ;
    ABL=vectorlength(AB);
    
    H[0]= AC[0]* sqrt(  1.0- innerproduct(AC,AB)* innerproduct(AC,AB) / ( ABL*ABL*ACL*ACL )    )  ;// h=AC*sin(teta)
    H[1]= AC[1]* sqrt(  1.0- innerproduct(AC,AB)* innerproduct(AC,AB) / ( ABL*ABL*ACL*ACL )    )  ;// h=AC*sin(teta
    H[2]= AC[2]* sqrt(  1.0- innerproduct(AC,AB)* innerproduct(AC,AB) / ( ABL*ABL*ACL*ACL )    )  ;// h=AC*sin(teta
    HL=vectorlength(H);
    
    H2L=HL+considerTimesHowFarElemets*epsilonECM;
    
    xECMindexA[0]=ECM_triangle_COM[0] + AC[0] * H2L/HL ;
    xECMindexA[1]=ECM_triangle_COM[1] + AC[1] * H2L/HL ;
    xECMindexA[2]=ECM_triangle_COM[2] + AC[2] * H2L/HL ;
    
    //--------------------------------------- indexB
    AC[0]= ECM_Node_Position[temp_ECM_node_B][0] -ECM_triangle_COM[0]  ;
    AC[1]= ECM_Node_Position[temp_ECM_node_B][1] -ECM_triangle_COM[1]  ;
    AC[2]= ECM_Node_Position[temp_ECM_node_B][2] -ECM_triangle_COM[2]  ;
    ACL=vectorlength(AC);
    
    AB[0]= ECM_Node_Position[temp_ECM_node_C][0] - ECM_Node_Position[temp_ECM_node_B][0] ;
    AB[1]= ECM_Node_Position[temp_ECM_node_C][1] - ECM_Node_Position[temp_ECM_node_B][1] ;
    AB[2]= ECM_Node_Position[temp_ECM_node_C][2] - ECM_Node_Position[temp_ECM_node_B][2] ;
    ABL=vectorlength(AB);
    
    H[0]= AC[0]* sqrt(  1.0- innerproduct(AC,AB)* innerproduct(AC,AB) / ( ABL*ABL*ACL*ACL )    )  ;// h=AC*sin(teta)
    H[1]= AC[1]* sqrt(  1.0- innerproduct(AC,AB)* innerproduct(AC,AB) / ( ABL*ABL*ACL*ACL )    )  ;// h=AC*sin(teta
    H[2]= AC[2]* sqrt(  1.0- innerproduct(AC,AB)* innerproduct(AC,AB) / ( ABL*ABL*ACL*ACL )    )  ;// h=AC*sin(teta
    HL=vectorlength(H);
    
    H2L=HL+considerTimesHowFarElemets*epsilonECM;
    
    xECMindexB[0]=ECM_triangle_COM[0] + AC[0] * H2L/HL ;
    xECMindexB[1]=ECM_triangle_COM[1] + AC[1] * H2L/HL ;
    xECMindexB[2]=ECM_triangle_COM[2] + AC[2] * H2L/HL ;
    
    //--------------------------------------- indexC
    AC[0]= ECM_Node_Position[temp_ECM_node_C][0] -ECM_triangle_COM[0]  ;
    AC[1]= ECM_Node_Position[temp_ECM_node_C][1] -ECM_triangle_COM[1]  ;
    AC[2]= ECM_Node_Position[temp_ECM_node_C][2] -ECM_triangle_COM[2]  ;
    ACL=vectorlength(AC);
    
    AB[0]= ECM_Node_Position[temp_ECM_node_A][0] - ECM_Node_Position[temp_ECM_node_C][0] ;
    AB[1]= ECM_Node_Position[temp_ECM_node_A][1] - ECM_Node_Position[temp_ECM_node_C][1] ;
    AB[2]= ECM_Node_Position[temp_ECM_node_A][2] - ECM_Node_Position[temp_ECM_node_C][2] ;
    ABL=vectorlength(AB);
    
    H[0]= AC[0]* sqrt(  1.0- innerproduct(AC,AB)* innerproduct(AC,AB) / ( ABL*ABL*ACL*ACL )    )  ;// h=AC*sin(teta)
    H[1]= AC[1]* sqrt(  1.0- innerproduct(AC,AB)* innerproduct(AC,AB) / ( ABL*ABL*ACL*ACL )    )  ;// h=AC*sin(teta
    H[2]= AC[2]* sqrt(  1.0- innerproduct(AC,AB)* innerproduct(AC,AB) / ( ABL*ABL*ACL*ACL )    )  ;// h=AC*sin(teta
    HL=vectorlength(H);
    
    H2L=HL+considerTimesHowFarElemets*epsilonECM;
    
    xECMindexC[0]=ECM_triangle_COM[0] + AC[0] * H2L/HL ;
    xECMindexC[1]=ECM_triangle_COM[1] + AC[1] * H2L/HL ;
    xECMindexC[2]=ECM_triangle_COM[2] + AC[2] * H2L/HL ;
    
    AB[0]=xECMindexB[0]-xECMindexA[0];
    AB[1]=xECMindexB[1]-xECMindexA[1];
    AB[2]=xECMindexB[2]-xECMindexA[2];
    
    AD[0]=D[0]-xECMindexA[0];
    AD[1]=D[1]-xECMindexA[1];
    AD[2]=D[2]-xECMindexA[2];
    
    BC[0]=xECMindexC[0]-xECMindexB[0];
    BC[1]=xECMindexC[1]-xECMindexB[1];
    BC[2]=xECMindexC[2]-xECMindexB[2];
    
    BD[0]=D[0]-xECMindexB[0];
    BD[1]=D[1]-xECMindexB[1];
    BD[2]=D[2]-xECMindexB[2];
    
    CA[0]=-xECMindexC[0]+xECMindexA[0];
    CA[1]=-xECMindexC[1]+xECMindexA[1];
    CA[2]=-xECMindexC[2]+xECMindexA[2];
    
    CD[0]=D[0]-xECMindexC[0];
    CD[1]=D[1]-xECMindexC[1];
    CD[2]=D[2]-xECMindexC[2];
    
    crossvector(temp,BC,BD);
    a1=innerproduct (ECM_triangle_ABxAC_unit_vector,temp);
    
    crossvector(temp,CA,CD);
    a2=innerproduct (ECM_triangle_ABxAC_unit_vector,temp);
    
    crossvector(temp,AB,AD);
    a3=innerproduct (ECM_triangle_ABxAC_unit_vector,temp);
    
}

void cellshift( double  Membrane_Node_Position [][3],double Actin_Node_Position[Actin_num_of_Nodes][3],double Chromatin_Bead_Position[Chromatin_num_of_Beads][3],double  ECM_Node_Position [][3] ,double  Membrane_Node_Velocity [][3],double Actin_Node_Velocity[Actin_num_of_Nodes][3],double Chromatin_Bead_Velocity[Chromatin_num_of_Beads][3], int Membrane_num_of_Nodes)
{
    for (int i=0; i<Membrane_num_of_Nodes ; i++)
    {
        Membrane_Node_Position[i][0] += membraneshiftinXdirection;
        Membrane_Node_Position[i][2] += membraneshiftinZdirection;
    }
    for (int i=0; i<Actin_num_of_Nodes ; i++)
    {
        Actin_Node_Position[i][0]= Actin_Node_Position[i][0]+membraneshiftinXdirection;
        Actin_Node_Position[i][2]= Actin_Node_Position[i][2]+membraneshiftinZdirection;
    }
    
    for (int i=0; i<Chromatin_num_of_Beads ; i++)
    {
        Chromatin_Bead_Position[i][0] += membraneshiftinXdirection;
        Chromatin_Bead_Position[i][2] += membraneshiftinZdirection;
    }
    
    for (int i=0; i<Membrane_num_of_Nodes ; i++)
    {
        Membrane_Node_Velocity[i][1] += cell_downward_speed;
        
    }
    for (int i=0; i<Actin_num_of_Nodes ; i++)
    {
        Actin_Node_Velocity[i][1]= Actin_Node_Velocity[i][1]+cell_downward_speed;
        
    }
    
    for (int i=0; i<Chromatin_num_of_Beads ; i++)
    {
        Chromatin_Bead_Velocity[i][1] += cell_downward_speed;
        
    }
    
}


//___________________migration

void CellCOM( double com[3],double  Membrane_Node_Position [][3],double  Actin_Node_Position [][3] ,double Chromatin_Bead_Position[Chromatin_num_of_Beads][3], int Membrane_num_of_Nodes)
{
    com[0]=0;
    com[1]=0;
    com[2]=0;
    double sigmaM=Membrane_num_of_Nodes*Membrane_Node_Mass+Actin_num_of_Nodes*Actin_Node_Mass+Chromatin_num_of_Beads*Chromatin_Bead_Mass;
    //----------------------membrane---------------------
    for(int i=0;i<Membrane_num_of_Nodes;i++)
    {
        com [0] += Membrane_Node_Mass*Membrane_Node_Position [i][0];
        com [1] += Membrane_Node_Mass*Membrane_Node_Position [i][1];
        com [2] += Membrane_Node_Mass*Membrane_Node_Position [i][2];
//        sigmaM += Membrane_Node_Mass;
    }
    //----------------------actin---------------------
    
    for(int i=0;i<Actin_num_of_Nodes;i++)
        //        for(int i=0;i<Actin_num_of_Nodes;i++)
    {
        com [0] += Actin_Node_Mass* Actin_Node_Position [i] [0];
        com [1] += Actin_Node_Mass* Actin_Node_Position [i] [1];
        com [2] += Actin_Node_Mass* Actin_Node_Position [i] [2];
//        sigmaM += Actin_Node_Mass;
        
    }
    
    //----------------------chromatin---------------------
    for(int i=0;i<Chromatin_num_of_Beads;i++)
    {
        com [0] += Chromatin_Bead_Mass*Chromatin_Bead_Position [i][0];
        com [1] += Chromatin_Bead_Mass*Chromatin_Bead_Position [i][1];
        com [2] += Chromatin_Bead_Mass*Chromatin_Bead_Position [i][2];
//        sigmaM += Chromatin_Bead_Mass;
    }
    
    //----------------------------------------
    com [0]=com[0]/sigmaM;
    com [1]=com[1]/sigmaM;
    com [2]=com[2]/sigmaM;
}



void Vector_transformation (double MV[3],double  M[3][3] ,double V[3])
{
    MV[0]= M[0][0] * V [0] + M[0][1] * V [1]  +M[0][2] * V [2] ;
    MV[1]= M[1][0] * V [0] + M[1][1] * V [1]  +M[1][2] * V [2] ;
    MV[2]= M[2][0] * V [0] + M[2][1] * V [1]  +M[2][2] * V [2] ;
}

//______________________Solvent functions

void initialzesolvent()
{
    
    if(BoundaryType ==1 || BoundaryType==2  || BoundaryType==5 || BoundaryType==6)
    {
        for (int i = 0; i < nsolvent; i++) // random place and velocity
        {
            xsolvent[0][i] = (double) (rand() % (Lx * 1000)) / 1000 - (Lx/2);
            xsolvent[1][i] = (double) (rand() % (Ly * 1000)) / 1000 - (Ly/2);
            xsolvent[2][i] = (double) (rand() % (Lz * 1000)) / 1000 - (Lz/2);
            
            vsolvent[2][i] =(double) (rand() % (int) (VmaxSolvent * 100000)) / 100000 - ((double)VmaxSolvent/2);
            vsolvent[1][i] =(double) (rand() % (int)  (VmaxSolvent * 100000)) / 100000 - ((double) VmaxSolvent/2);
            vsolvent[0][i] =(double) (rand() % (int) (VmaxSolvent * 100000)) / 100000 - ((double)VmaxSolvent/2);
            
        }
    }
    if(BoundaryType ==3 || BoundaryType==4)
    {
        for (int i = 0; i < nsolvent; i++) // random place and velocity
        {
            xsolvent[0][i] = (double) (rand() % (Lx * 1000)) / 1000 - (Lx/2);
            xsolvent[1][i] = (double) (rand() % (Ly * 1000)) / 1000 - (Ly/2);
            xsolvent[2][i] = (double) (rand() % (Lz * 1000)) / 1000 - (Lz/2);
            
            vsolvent[2][i] = (double) (rand() % (int) (VmaxSolvent * 100000)) / 100000 - ((double)VmaxSolvent/2);
            vsolvent[1][i] = (double) (rand() % (int)  (VmaxSolvent * 100000)) / 100000 - ((double)VmaxSolvent/2);
            vsolvent[0][i] =Vinitialflow;
            
        }
    }
}

void solventBoundaries()
{
    
    if(BoundaryType==1  || BoundaryType==3)
    {
        for (int j = 0; j < nsolvent; j++)
        {
            // ___________________ Boundary Condition Solvent
            //x:
            if (xsolvent[0][j] > Lx/2.0)
                xsolvent[0][j] = xsolvent[0][j] - Lx ;
            
            if (xsolvent[0][j] < -Lx/2.0)
                xsolvent[0][j] = xsolvent[0][j] + Lx ;
            
            // slip condition
            //y and z :
            double ly =Ly/2.0-wallthickness;
            double lz =Lz/2.0-wallthickness;
            if (   (xsolvent[1][j]>ly & vsolvent[1][j] >0.0)   ||   (xsolvent[1][j]<-ly & vsolvent[1][j] <0.0)  )
            {
                vsolvent[1][j] = -vsolvent[1][j];
            }
            
            else  if (   (xsolvent[2][j]>lz & vsolvent[2][j] >0.0)   ||   (xsolvent[2][j]<-lz & vsolvent[2][j] <0.0)  )
            {
                vsolvent[2][j] = -vsolvent[2][j];
            }
        }
    }
    
    else if(BoundaryType==2 ||  BoundaryType==4  || BoundaryType==5)
    {
        for (int j = 0; j < nsolvent; j++)
        {
            // ___________________ Boundary Condition Solvent
            //x:
            if (xsolvent[0][j] > Lx/2.0)
                xsolvent[0][j] = xsolvent[0][j] - Lx ;
            
            if (xsolvent[0][j] < -Lx/2.0)
                xsolvent[0][j] = xsolvent[0][j] + Lx ;
            
            // no-slip condition
            //y and z :
            double ly =Ly/2.0-wallthickness;
            double lz =Lz/2.0-wallthickness;
            if (   (xsolvent[1][j]>ly & vsolvent[1][j] >0.0)   ||   (xsolvent[1][j]<-ly & vsolvent[1][j] <0.0)  )
            {
                vsolvent[0][j] = -vsolvent[0][j];
                vsolvent[1][j] = -vsolvent[1][j];
                vsolvent[2][j] = -vsolvent[2][j];
            }
            
            else    if (   (xsolvent[2][j]>lz & vsolvent[2][j] >0.0)   ||   (xsolvent[2][j]<-lz & vsolvent[2][j] <0.0)  )
            {
                vsolvent[0][j] = -vsolvent[0][j];
                vsolvent[1][j] = -vsolvent[1][j];
                vsolvent[2][j] = -vsolvent[2][j];
            }
        }
    }
    
    
    else if(BoundaryType==6)
    {
        for (int j = 0; j < nsolvent; j++)
        {
            // ___________________ Boundary Condition Solvent
            //x:
            if (xsolvent[0][j] > Lx/2.0)
                xsolvent[0][j] = xsolvent[0][j] - Lx ;
            
            if (xsolvent[0][j] < -Lx/2.0)
                xsolvent[0][j] = xsolvent[0][j] + Lx ;
            
            // shear flow in z plane and slip in y plane
            //y and z :
            double ly =Ly/2.0-wallthickness;
            double lz =Lz/2.0-wallthickness;
            
            if (   (xsolvent[1][j]>ly & vsolvent[1][j] >0.0)   ||   (xsolvent[1][j]<-ly & vsolvent[1][j] <0.0)  )
            {
                vsolvent[1][j] = -vsolvent[1][j];
            }
            else   if (   (xsolvent[2][j]>lz & vsolvent[2][j] >0.0)   ||   (xsolvent[2][j]<-lz & vsolvent[2][j] <0.0)  )
            {
                vsolvent[2][j] = -vsolvent[2][j];
            }
            
            
            if (  xsolvent[1][j]>lz     )
            {
                vsolvent[0][j] = +Vshear;
            }
            else   if (  xsolvent[1][j]<-lz     )
            {
                vsolvent[0][j] = -Vshear;
            }
            
        }
        
    }
    
    
    
}


void collision()
{
    double u[3][nb];
    
    double n1, n2, n3; // arbitrary vector in space around which particles rotate
    
    double x, y, z;
    
    
    for (int i = 0; i < nb; i++)
    {
        u[0][i] = 0;
        u[1][i] = 0;
        u[2][i] = 0;
    }
    
    updateHowmanySolventsineachBox();
    
    // calculate velocity of center of mass for each box
    for (int i = 0; i < nb; i++)
    {
        for (int j = 0; j < howmanySolventsineachBox[i]; j++)
        {
            u[0][i] += vsolvent[0][eachSolventinwichBox[i][j]];
            u[1][i] += vsolvent[1][eachSolventinwichBox[i][j]];
            u[2][i] += vsolvent[2][eachSolventinwichBox[i][j]];
        }
        
        if (howmanySolventsineachBox[i] != 0)
        {
            u[0][i] = (double) u[0][i] / howmanySolventsineachBox[i];
            u[1][i] = (double) u[1][i] / howmanySolventsineachBox[i];
            u[2][i] = (double) u[2][i] / howmanySolventsineachBox[i];
        }
        if (howmanySolventsineachBox[i] == 0)
        {
            u[0][i] = 0.0;
            u[1][i] = 0.0;
            u[2][i] = 0.0;
        }
    }
    
    // rotate velocities in boxes
    for (int i = 0; i < nb; i++)
    {
        n1 = (double) (rand() % 10000) / 10000 - 0.5;  // create random vector  Random coordinates -0.5 +0.5
        n2 = (double) (rand() % 10000) / 10000 - 0.5;
        n3 = (double) (rand() % 10000) / 10000 - 0.5;
        
        double r = (double) sqrt(n1 * n1 + n2 * n2 + n3 * n3);
        
        n1 = (double) n1 / r;
        n2 = (double) n2 / r;
        n3 = (double) n3 / r;
        
        for (int j = 0; j < howmanySolventsineachBox[i]; j++)
        {
            if (eachSolventinwichBox[i][j] >= 0 && eachSolventinwichBox[i][j] < nsolvent)
            {
                
                x =-u[0][i]+ vsolvent[0][eachSolventinwichBox[i][j]];
                y =-u[1][i]+ vsolvent[1][eachSolventinwichBox[i][j]];
                z =-u[2][i]+ vsolvent[2][eachSolventinwichBox[i][j]];
                
                vsolvent[0][eachSolventinwichBox[i][j]]=u[0][i] +  n1*(n1*x+n2*y+n3*z) + (-n3*y + n2*z );
                vsolvent[1][eachSolventinwichBox[i][j]]=u[1][i] +  n2*(n1*x+n2*y+n3*z) + ( n3*x - n1*z );
                vsolvent[2][eachSolventinwichBox[i][j]]=u[2][i] +  n3*(n1*x+n2*y+n3*z) + (-n2*x + n1*y );
                
            }
        }
    }
    
}


void updateHowmanySolventsineachBox()
{
    
    //--------------------------------+UPDATE  howmanySolventsineachBox-----------------
    int x, y, z;
    int index;
    for (int i = 0; i < nb; i++)
    {
        howmanySolventsineachBox[i] = 0;
    }
    // specify the box of each particle
    for (int i = 0; i < nsolvent; i++)
    {
        if ( xsolvent[0][i] > -Lx/2 && xsolvent[0][i] < Lx/2 && xsolvent[1][i] > -Ly/2 && xsolvent[1][i] < Ly/2 && xsolvent[2][i] > -Lz/2 && xsolvent[2][i] < Lz/2)
        {
            x =floor (xsolvent[0][i] + (Lx/2));// x y z will lable the corresponding box of each particle
            y =floor (xsolvent[1][i] + (Ly/2));
            z =floor (xsolvent[2][i] + (Lz/2));
            
            index = (x *  Ly * Lz) + (z * Ly) + y;  // the lable of each box in another way
            eachSolventinwichBox[index][howmanySolventsineachBox[index]] = i;                // says that i th particle is here!
            howmanySolventsineachBox[index]++;                           //says that  this box have how many particles!
        }
        
    }
    
    //------------------------------------------------------------------
    
    
    
}




void findcom(double vcom[3][Membrane_num_of_Triangles],double pcom[3][Membrane_num_of_Triangles],int Membrane_triangle_list[Membrane_num_of_Triangles][3], double (&Membrane_Node_Position)[][3], double Membrane_Node_Velocity[][3])
{
    double p1, p2, p3, v1, v2, v3;
    for (int i = 0; i < Membrane_num_of_Triangles; i++)
    {
        // find center of mass of triangle (place)
        p1 = Membrane_Node_Position[Membrane_triangle_list[i][0]][0] + Membrane_Node_Position[Membrane_triangle_list[i][1]][0] + Membrane_Node_Position[Membrane_triangle_list[i][2]][0];
        p2 = Membrane_Node_Position[Membrane_triangle_list[i][0]][1] + Membrane_Node_Position[Membrane_triangle_list[i][1]][1] + Membrane_Node_Position[Membrane_triangle_list[i][2]][1];
        p3 = Membrane_Node_Position[Membrane_triangle_list[i][0]][2] + Membrane_Node_Position[Membrane_triangle_list[i][1]][2] + Membrane_Node_Position[Membrane_triangle_list[i][2]][2];
        
        pcom[0][i] = (double) p1 / 3.0;
        pcom[1][i] = (double) p2 / 3.0;
        pcom[2][i] = (double) p3 / 3.0;
        
        // center of mass velocity
        v1 = Membrane_Node_Velocity[Membrane_triangle_list[i][0]][0] + Membrane_Node_Velocity[Membrane_triangle_list[i][1]][0] + Membrane_Node_Velocity[Membrane_triangle_list[i][2]][0];
        v2 = Membrane_Node_Velocity[Membrane_triangle_list[i][0]][1] + Membrane_Node_Velocity[Membrane_triangle_list[i][1]][1] + Membrane_Node_Velocity[Membrane_triangle_list[i][2]][1];
        v3 = Membrane_Node_Velocity[Membrane_triangle_list[i][0]][2] + Membrane_Node_Velocity[Membrane_triangle_list[i][1]][2] + Membrane_Node_Velocity[Membrane_triangle_list[i][2]][2];
        
        vcom[0][i] = (double) v1 / 3.0;
        vcom[1][i] = (double) v2 / 3.0;
        vcom[2][i] = (double) v3 / 3.0;
    }
}


void interaction(double Membrane_Node_Position[][3], double Membrane_Node_Velocity[][3], int Membrane_triangle_list[Membrane_num_of_Triangles][3], int Membrane_Normal_direction[Membrane_num_of_Triangles][2])
{
    
    double pcom[3]; // place of center of mass of triangle
    double vcom[3]; // velocity of center of mass of triangle
    double vcomnew[3]; // velocity of center of mass of triangle
    int neiborBoxes[27];
    int part ; // index of paticle
    int index;
    int x, y, z;
    double dis;// distance between solvent particle and com of triangle
    double normal[3],a1a2[3],a1a3[3]; // normal vector of membrane
    double TriangleSolvent[3];
    double inORout;
    double vsvt[3];
    double v1, v2, v3;
    
    for (int i = 0; i < Membrane_num_of_Triangles; i++)  // for each triangle:
    {
        pcom[0]=(Membrane_Node_Position[ Membrane_triangle_list[i][0]][0] + Membrane_Node_Position[Membrane_triangle_list[i][1]][0]+Membrane_Node_Position[ Membrane_triangle_list[i][2]][0])/3.0;
        pcom[1]=(Membrane_Node_Position[ Membrane_triangle_list[i][0]][1] + Membrane_Node_Position[Membrane_triangle_list[i][1]][1]+Membrane_Node_Position [ Membrane_triangle_list[i][2]][1])/3.0;
        pcom[2]=(Membrane_Node_Position[ Membrane_triangle_list[i][0]][2] + Membrane_Node_Position[Membrane_triangle_list[i][1]][2]+Membrane_Node_Position [ Membrane_triangle_list[i][2]][2])/3.0;
        
        vcom[0]=(Membrane_Node_Velocity[ Membrane_triangle_list[i][0]][0] + Membrane_Node_Velocity[Membrane_triangle_list[i][1]][0] + Membrane_Node_Velocity[ Membrane_triangle_list[i][2]][0])/3.0;
        vcom[1]=(Membrane_Node_Velocity[ Membrane_triangle_list[i][0]][1] + Membrane_Node_Velocity[Membrane_triangle_list[i][1]][1] + Membrane_Node_Velocity[ Membrane_triangle_list[i][2]][1])/3.0;
        vcom[2]=(Membrane_Node_Velocity[ Membrane_triangle_list[i][0]][2] + Membrane_Node_Velocity[Membrane_triangle_list[i][1]][2] + Membrane_Node_Velocity[ Membrane_triangle_list[i][2]][2])/3.0;
        
        x =floor (pcom[0] + (Lx/2.0));
        y =floor (pcom[1] + (Ly/2.0));
        z =floor (pcom[2] + (Lz/2.0));
        
        // Main box and its neigbors:
        neiborBoxes[0]=index  = (x *  Ly * Lz) + (z * Ly) + y;
        
        //  cout<< x<<"  " <<y<<"  "<<z<<"   pcom:"<< pcom[0][i]<<"  " <<pcom[1][i]<<"  "<<pcom[2][i]<<(x *  Ly * Lz) + (z * Ly) + y<<endl;
        
        neiborBoxes[1]  = ((x+1) *  Ly * Lz) + (z * Ly) + y;
        neiborBoxes[2]  = (x *  Ly * Lz) + (z * Ly) + (y+1);
        neiborBoxes[3]  = (x *  Ly * Lz) + ((z+1) * Ly) + y;
        neiborBoxes[4]  = ((x-1) *  Ly * Lz) + (z * Ly) + y;
        neiborBoxes[5]  = (x *  Ly * Lz) + (z * Ly) + (y-1);
        neiborBoxes[6]  = (x *  Ly * Lz) + ((z-1) * Ly) + y;
        
        neiborBoxes[7]  = ((x+1) *  Ly * Lz) + (z * Ly) + (y+1);
        neiborBoxes[8]  = ((x+1) *  Ly * Lz) + (z * Ly) + (y-1);
        neiborBoxes[9]  = ((x+1) *  Ly * Lz) + ((z+1) * Ly) + y;
        neiborBoxes[10] = ((x+1) *  Ly * Lz) + ((z-1) * Ly) + y;
        neiborBoxes[11]  = ((x-1) *  Ly * Lz) + (z * Ly) + (y+1);
        neiborBoxes[12]  = ((x-1) *  Ly * Lz) + (z * Ly) + (y-1);
        neiborBoxes[13]  = ((x-1) *  Ly * Lz) + ((z+1) * Ly) + y;
        neiborBoxes[14] = ((x-1) *  Ly * Lz) + ((z-1) * Ly) + y;
        neiborBoxes[15] = (x *  Ly * Lz) + ((z+1) * Ly) + (y+1);
        neiborBoxes[16] = (x *  Ly * Lz) + ((z-1) * Ly) + (y+1);
        neiborBoxes[17] = (x *  Ly * Lz) + ((z+1) * Ly) + (y-1);
        neiborBoxes[18] = (x *  Ly * Lz) + ((z-1) * Ly) + (y-1);
        
        neiborBoxes[19]  = ((x+1) *  Ly * Lz) + ((z+1) * Ly) + (y+1);
        neiborBoxes[20]  = ((x+1) *  Ly * Lz) + ((z-1) * Ly) + (y+1);
        neiborBoxes[21]  = ((x+1) *  Ly * Lz) + ((z+1) * Ly) + (y-1);
        neiborBoxes[22]  = ((x+1) *  Ly * Lz) + ((z-1) * Ly) + (y-1);
        neiborBoxes[23]  = ((x-1) *  Ly * Lz) + ((z+1) * Ly) + (y+1);
        neiborBoxes[24]  = ((x-1) *  Ly * Lz) + ((z-1) * Ly) + (y+1);
        neiborBoxes[25]  = ((x-1) *  Ly * Lz) + ((z+1) * Ly) + (y-1);
        neiborBoxes[26]  = ((x-1) *  Ly * Lz) + ((z-1) * Ly) + (y-1);
        
        for (int box=0;box<27;box++)
        {
            if(neiborBoxes[box]>-1 & neiborBoxes[box] < nb )
            {
                
                //  cout<< neiborBoxes[box] <<"   "<<howmanySolventsineachBox[neiborBoxes[box] ]<<endl;
                for (int j = 0; j < howmanySolventsineachBox[neiborBoxes[box]  ]; j++)
                {
                    // checking the interaction conditions:
                    part = eachSolventinwichBox[neiborBoxes[box]][j];
                    dis =sqrt( (pcom[0] - xsolvent[0][part]) * (pcom[0] - xsolvent[0][part]) + (pcom[1] - xsolvent[1][part]) * (pcom[1]- xsolvent[1][part]) + (pcom[2] - xsolvent[2][part]) * (pcom[2] - xsolvent[2][part]));
                    
                    if (  dis < sqrt(0.43*Node_radius * 0.43*Node_radius + lbs*lbs )  )  // is silvent in proper sphere?
                    {
                        // cout<<"interactio"<<endl;
                        a1a2[0]=Membrane_Node_Position[ Membrane_triangle_list[i][1]][0]-Membrane_Node_Position[ Membrane_triangle_list[i][0]][0];
                        a1a2[1]=Membrane_Node_Position[ Membrane_triangle_list[i][1]][1]-Membrane_Node_Position[ Membrane_triangle_list[i][0]][1];
                        a1a2[2]=Membrane_Node_Position[ Membrane_triangle_list[i][1]][2]-Membrane_Node_Position[ Membrane_triangle_list[i][0]][2];
                        a1a3[0]=Membrane_Node_Position[ Membrane_triangle_list[i][2]][0]-Membrane_Node_Position[ Membrane_triangle_list[i][0]][0];
                        a1a3[1]=Membrane_Node_Position[ Membrane_triangle_list[i][2]][1]-Membrane_Node_Position[ Membrane_triangle_list[i][0]][1];
                        a1a3[2]=Membrane_Node_Position[ Membrane_triangle_list[i][2]][2]-Membrane_Node_Position[ Membrane_triangle_list[i][0]][2];
                        crossvector(normal,a1a2,a1a3);
                        normal[0]=normal[0]*Membrane_Normal_direction[i][1];
                        normal[1]=normal[1]*Membrane_Normal_direction[i][1];
                        normal[2]=normal[2]*Membrane_Normal_direction[i][1];
                        
                        TriangleSolvent[0]=xsolvent[0][part]-pcom[0];
                        TriangleSolvent[1]=xsolvent[1][part]-pcom[1];
                        TriangleSolvent[2]=xsolvent[2][part]-pcom[2];
                        inORout=innerproduct(TriangleSolvent,normal)/vectorlength(normal);
                        
                        if    (  abs( inORout )<lbs  )  // is solvent close enough to surface?
                        {
                            
                            vsvt[0] = vsolvent[0][part]-vcom[0];
                            vsvt[1] = vsolvent[1][part]-vcom[1];
                            vsvt[2] = vsolvent[2][part]-vcom[2];
                            
                            v1 = vsolvent[0][part];
                            v2 = vsolvent[1][part];
                            v3 = vsolvent[2][part];
                            
                            if (   innerproduct(vsvt,normal)*inORout<0    ) // is solvent moving toward triangle?
                            {
                                vsolvent[0][part] = vsolvent[0][part] - (double) (6*Membrane_Node_Mass / (msolvent + 3*Membrane_Node_Mass)) * (vsolvent[0][part] - vcom[0]);
                                vsolvent[1][part] = vsolvent[1][part] - (double) (6*Membrane_Node_Mass / (msolvent + 3*Membrane_Node_Mass)) * (vsolvent[1][part] - vcom[1]);
                                vsolvent[2][part] = vsolvent[2][part] - (double) (6*Membrane_Node_Mass / (msolvent + 3*Membrane_Node_Mass)) * (vsolvent[2][part] - vcom[2]);
                                
                                vcomnew[0]= vcom[0] + (double) (2*msolvent / (msolvent + 3*Membrane_Node_Mass)) * (v1 - vcom[0]);
                                vcomnew[1]= vcom[1] + (double) (2*msolvent / (msolvent + 3*Membrane_Node_Mass)) * (v2 - vcom[1]);
                                vcomnew[2]= vcom[2] + (double) (2*msolvent / (msolvent + 3*Membrane_Node_Mass)) * (v3 - vcom[2]);
                                
                                for (int i1 = 0; i1 < 3; i1++)
                                {
                                    Membrane_Node_Velocity[Membrane_triangle_list[i][i1]][0] +=- vcom[0]+ vcomnew[0];
                                    Membrane_Node_Velocity[Membrane_triangle_list[i][i1]][1] +=- vcom[1]+ vcomnew[1];
                                    Membrane_Node_Velocity[Membrane_triangle_list[i][i1]][2] +=- vcom[2]+ vcomnew[2];
                                }
                                
                            }
                        }
                    }
                }
            }
        }
    }
    
    
}








//-----------------------------thermostat------------------
double kineticenergymembrane (double Membrane_Node_Velocity[][3], int Membrane_num_of_Nodes)
{
    double vcom[3];
    vcom[0]=0.0;
    vcom[1]=0.0;
    vcom[2]=0.0;
    
    for ( int i=0 ; i<=Membrane_num_of_Nodes ; i++)
    {
        vcom[0]=vcom[0]+ Membrane_Node_Velocity[i][0];
        vcom[1]=vcom[1]+ Membrane_Node_Velocity[i][1];
        vcom[2]=vcom[2]+ Membrane_Node_Velocity[i][2];
    }
    
    vcom[0]=vcom[0]/((double) Membrane_num_of_Nodes);
    vcom[1]=vcom[1]/((double) Membrane_num_of_Nodes);
    vcom[2]=vcom[2]/((double) Membrane_num_of_Nodes);
    
    
    double K;
    K=0;
    for ( int i=0 ; i<=Membrane_num_of_Nodes ; i++)
    {
        K+=(Membrane_Node_Velocity[i][0]-vcom[0]) *(Membrane_Node_Velocity[i][0]-vcom[0]) +
        (Membrane_Node_Velocity[i][1]-vcom[1]) *(Membrane_Node_Velocity[i][1]-vcom[1]) +
        (Membrane_Node_Velocity[i][2]-vcom[2]) *(Membrane_Node_Velocity[i][2]-vcom[2]) ;
        
    }
    
    K=Membrane_Node_Mass*K/(2.0);
    
    
    return K;
    
}

double kineticenergyactin (double  Actin_Node_Velocity [][3])
{
    double vcom[3];
    vcom[0]=0.0;
    vcom[1]=0.0;
    vcom[2]=0.0;
    for ( int i=Actin_Membrane_shared_num_of_Nodes ; i<=Actin_num_of_Nodes ; i++)
    {
        vcom[0]= vcom[0]+ Actin_Node_Velocity[i][0];
        vcom[1]= vcom[1]+ Actin_Node_Velocity[i][1];
        vcom[2]= vcom[2]+ Actin_Node_Velocity[i][2];
        
    }
    
    vcom[0]= vcom[0]/((double)Actin_num_of_Nodes);
    vcom[1]= vcom[1]/((double)Actin_num_of_Nodes);
    vcom[2]= vcom[2]/((double)Actin_num_of_Nodes);
    
    double K;
    K=0;
    for ( int i=Actin_Membrane_shared_num_of_Nodes ; i<=Actin_num_of_Nodes ; i++)
    {
        K=K+ (Actin_Node_Velocity[i][0]-vcom[0])*(Actin_Node_Velocity[i][0]-vcom[0])+
        (Actin_Node_Velocity[i][1]-vcom[1])*(Actin_Node_Velocity[i][1]-vcom[1])+
        (Actin_Node_Velocity[i][2]-vcom[2])*(Actin_Node_Velocity[i][2]-vcom[2]);
        
    }
    
    
    
    K=Actin_Node_Mass*K/(2.0);
    
    return K;
    
}

double kineticenergyECM (double  ECM_Node_Velocity [][3])
{
    
    
    double K;
    K=0;
    for ( int i=0 ; i<=ECM_num_of_Nodes ; i++)
    {
        K=K+ ECM_Node_Velocity[i][0]*ECM_Node_Velocity[i][0]+
        ECM_Node_Velocity[i][1]*ECM_Node_Velocity[i][1]+
        ECM_Node_Velocity[i][2]*ECM_Node_Velocity[i][2];
        
    }
    
    
    
    K=ECM_Node_Mass*K/(2.0);
    
    return K;
    
}



double kineticenergysolvent ()
{
    double K;
    K=0;
    for ( int i=0 ; i<=nsolvent ; i++)
    {
        K=K+ vsolvent[0][i]*vsolvent[0][i]+
        vsolvent[1][i]*vsolvent[1][i]+
        vsolvent[2][i]*vsolvent[2][i];
        
    }
    
    
    
    K=msolvent*K/(2.0);
    
    
    return K;
    
}


double kineticenergychromatin (double Chromatin_Bead_Velocity[Chromatin_num_of_Beads][3])
{
    double vcom[3];
    vcom[0]=0.0;
    vcom[1]=0.0;
    vcom[2]=0.0;
    for ( int i=0 ; i<=Chromatin_num_of_Beads ; i++)
    {
        vcom[0]= vcom[0]+ Chromatin_Bead_Velocity[i][0];
        vcom[1]= vcom[1]+ Chromatin_Bead_Velocity[i][1];
        vcom[2]= vcom[2]+ Chromatin_Bead_Velocity[i][2];
    }
    vcom[0]= vcom[0]/((double)Chromatin_num_of_Beads);
    vcom[1]= vcom[1]/((double)Chromatin_num_of_Beads);
    vcom[2]= vcom[2]/((double)Chromatin_num_of_Beads);
    
    double K;
    K=0;
    for ( int i=0 ; i<=Chromatin_num_of_Beads ; i++)
    {
        K+=(Chromatin_Bead_Velocity[i][0]-vcom[0])*(Chromatin_Bead_Velocity[i][0]-vcom[0])+
        (Chromatin_Bead_Velocity[i][1]-vcom[1])*(Chromatin_Bead_Velocity[i][1]-vcom[1])+
        (Chromatin_Bead_Velocity[i][2]-vcom[2])*(Chromatin_Bead_Velocity[i][2]-vcom[2]);
    }
    
    K=Chromatin_Bead_Mass*K/(2.0);
    
    
    return K;
    
}


void Thermostat(int istep, double Membrane_Node_Velocity[][3], double Actin_Node_Velocity[][3], double Chromatin_Bead_Velocity[Chromatin_num_of_Beads][3], double ECM_Node_Velocity[][3], int Membrane_num_of_Nodes)
{
    double alpha;
    //----------------------membrane---------------------
    
    alpha=sqrt(      (3* Membrane_num_of_Nodes  *KT) / kineticenergymembrane( Membrane_Node_Velocity, Membrane_num_of_Nodes )         );/// NOTE THAT THERMOSTATE IS FOR MEMBRANE YET. IN ORDER TO
    /// UPDATE IT FOR BOTH THE MEMBRANE AND NUCLEI WE HAVE TO
    /// WRITE  alpha=sqrt(      (2*3*Membrane_num_of_Nodes*KT) / kineticenergy( Membrane_Node_Velocity,vnuclei
    for(int i=0;i<Membrane_num_of_Nodes;i++)
    {
        Membrane_Node_Velocity [i][0]*=alpha;
        Membrane_Node_Velocity [i][1]*=alpha;
        Membrane_Node_Velocity [i][2]*=alpha;
    }
    
    
    
    
    //----------------------actin---------------------
    alpha=sqrt(      (3* Actin_num_of_Nodes  *KT) / kineticenergyactin( Actin_Node_Velocity )         );/// NOTE THAT THERMOSTATE IS FOR MEMBRANE YET. IN ORDER TO
    /// UPDATE IT FOR BOTH THE MEMBRANE AND NUCLEI WE HAVE TO
    /// WRITE  alpha=sqrt(      (2*3*Membrane_num_of_Nodes*KT) / kineticenergy( Membrane_Node_Velocity,vnuclei
    for(int i=Actin_Membrane_shared_num_of_Nodes;i<Actin_num_of_Nodes;i++)
    {
        Actin_Node_Velocity [i] [0]= Actin_Node_Velocity [i] [0]*alpha;
        Actin_Node_Velocity [i] [1]= Actin_Node_Velocity [i] [1]*alpha;
        Actin_Node_Velocity [i] [2]= Actin_Node_Velocity [i] [2]*alpha;
    }
    
    
    
    //----------------------solvent---------------------
    
    alpha=sqrt(      (3* nsolvent  *KT) / kineticenergysolvent(  )         );/// NOTE THAT THERMOSTATE IS FOR MEMBRANE YET. IN ORDER TO
    /// UPDATE IT FOR BOTH THE MEMBRANE AND NUCLEI WE HAVE TO
    /// WRITE  alpha=sqrt(      (2*3*Membrane_num_of_Nodes*KT) / kineticenergy( Membrane_Node_Velocity,vnuclei
    for(int i=0;i<nsolvent;i++)
    {
        vsolvent [0][i]=vsolvent [0][i]*alpha;
        vsolvent [1][i]=vsolvent [1][i]*alpha;
        vsolvent [2][i]=vsolvent [2][i]*alpha;
    }
    
    
    //----------------------chromatin---------------------
    
    if(istep%ThermostatOnChromatin==0)
    {
        alpha=sqrt(      (3*Chromatin_num_of_Beads*KT) / kineticenergychromatin( Chromatin_Bead_Velocity )         );
        for(int i=0;i<Chromatin_num_of_Beads;i++)
        {
            Chromatin_Bead_Velocity [i][0] *= alpha;
            Chromatin_Bead_Velocity [i][1] *= alpha;
            Chromatin_Bead_Velocity [i][2] *= alpha;
        }
    }
    
}

//-----------------------------thermostat------------------






//---------------------------restart-----------------------
void restartsave( double Membrane_Node_Position[][3], double Membrane_Node_Velocity [][3], int Membrane_triangle_list[Membrane_num_of_Triangles][3], int Membrane_Normal_direction[Membrane_num_of_Triangles][2], int Membrane_Triangle_Pair_Nodes[][4], int bondslist[][2], double Actin_Node_Position[Actin_num_of_Nodes][3], double Actin_Node_Velocity[Actin_num_of_Nodes][3], double Actin_Node_Pair_List[][3], double Chromatin_Bead_Position[Chromatin_num_of_Beads][3], double Chromatin_Bead_Velocity[Chromatin_num_of_Beads][3], double ECM_Node_Position[][3], double ECM_Node_Velocity [][3], int Membrane_num_of_Triangle_Pairs, int Actin_num_of_Bonds, int Membrane_num_of_Nodes)
{
    ofstream restart;
    restart.open("restart_1.txt");
    restart << std:: fixed;
    
    for (int i=0; i<Membrane_num_of_Nodes ; i++)
    {
        for (int j=0; j<3 ; j++)
        {
            restart<< Membrane_Node_Position[i][j]<<endl;
        }
    }
    
    for (int i=0; i<Membrane_num_of_Nodes ; i++)
    {
        for (int j=0; j<3 ; j++)
        {
            restart<< Membrane_Node_Velocity[i][j]<<endl;
        }
    }
    
    
    for (int i=0; i<Membrane_num_of_Triangles ; i++)
    {
        for (int j=0; j<3 ; j++)
        {
            restart<< Membrane_triangle_list[i][j]<<endl;
        }
    }
    
    for (int i=0; i<Membrane_num_of_Triangles ; i++)
    {
        for (int j=0; j<2 ; j++)
        {
            restart<< Membrane_Normal_direction[i][j]<<endl;
        }
    }
    
    
    for (int i=0; i<Membrane_num_of_Triangle_Pairs ; i++)
    {
        for (int j=0; j<4 ; j++)
        {
            restart<< Membrane_Triangle_Pair_Nodes[i][j]<<endl;
        }
    }
    
    for (int i=0; i<Membrane_num_of_Triangle_Pairs ; i++)
    {
        for (int j=0; j<2 ; j++)
        {
            restart<< bondslist[i][j]<<endl;
        }
    }
    
    for (int i=0; i<Actin_num_of_Nodes ; i++)
    {
        for (int j=0; j<3 ; j++)
        {
            restart<< Actin_Node_Position[i][j]<<endl;
        }
    }
    for (int i=0; i<Actin_num_of_Nodes ; i++)
    {
        for (int j=0; j<3 ; j++)
        {
            restart<< Actin_Node_Velocity[i][j]<<endl;
        }
    }
    
    
    for (int i=0; i<Actin_num_of_Bonds ; i++)
    {
        for (int j=0; j<3 ; j++)
        {
            restart<< Actin_Node_Pair_List[i][j]<<endl;
        }
    }
    
    for (int i=0; i<nsolvent ; i++)
    {
        for (int j=0; j<3 ; j++)
        {
            restart<< xsolvent[j][i]<<endl;
        }
    }
    
    
    for (int i=0; i<nsolvent ; i++)
    {
        for (int j=0; j<3 ; j++)
        {
            restart<< vsolvent[j][i]<<endl;
        }
    }
    
    
    
    for (int i=0; i<Chromatin_num_of_Beads ; i++)
    {
        for (int j=0; j<3 ; j++)
        {
            restart<< Chromatin_Bead_Position[i][j]<<endl;
        }
    }
    
    
    for (int i=0; i<Chromatin_num_of_Beads ; i++)
    {
        for (int j=0; j<3 ; j++)
        {
            restart<< Chromatin_Bead_Velocity[i][j]<<endl;
        }
    }
    
    
    for (int i=0; i<ECM_num_of_Nodes ; i++)
    {
        for (int j=0; j<3 ; j++)
        {
            restart<< ECM_Node_Position[i][j]<<endl;
        }
    }
    
    
    
    for (int i=0; i<ECM_num_of_Nodes ; i++)
    {
        for (int j=0; j<3 ; j++)
        {
            restart<< ECM_Node_Velocity[i][j]<<endl;
        }
    }
    
    
    
    
}


void restartread(double Membrane_Node_Position [][3], double Membrane_Node_Velocity [][3], int Membrane_triangle_list[Membrane_num_of_Triangles][3], int Membrane_Normal_direction[Membrane_num_of_Triangles][2], int Membrane_Triangle_Pair_Nodes[][4], int bondslist[][2], double Actin_Node_Position[Actin_num_of_Nodes][3], double Actin_Node_Velocity[Actin_num_of_Nodes][3],double Actin_Node_Pair_List[][3], double Chromatin_Bead_Position[Chromatin_num_of_Beads][3], double Chromatin_Bead_Velocity[Chromatin_num_of_Beads][3], double ECM_Node_Position [][3], double ECM_Node_Velocity [][3], int Membrane_num_of_Triangle_Pairs, int Actin_num_of_Bonds, int Membrane_num_of_Nodes)
{
    
    
    ifstream restart;
    //    restart.open("restart-backup.txt");
    restart.open("restart_s_2.txt");
    
    for (int i=0; i<Membrane_num_of_Nodes ; i++)
    {
        for (int j=0; j<3 ; j++)
        {
            restart>> Membrane_Node_Position[i][j];
        }
    }
    
    for (int i=0; i<Membrane_num_of_Nodes ; i++)
    {
        for (int j=0; j<3 ; j++)
        {
            restart>> Membrane_Node_Velocity[i][j];
        }
    }
    
    
    for (int i=0; i<Membrane_num_of_Triangles ; i++)
    {
        for (int j=0; j<3 ; j++)
        {
            restart>> Membrane_triangle_list[i][j];
        }
    }
    
    for (int i=0; i<Membrane_num_of_Triangles ; i++)
    {
        for (int j=0; j<2 ; j++)
        {
            restart>> Membrane_Normal_direction[i][j];
        }
    }
    
    
    for (int i=0; i<Membrane_num_of_Triangle_Pairs ; i++)
    {
        for (int j=0; j<4 ; j++)
        {
            restart>> Membrane_Triangle_Pair_Nodes[i][j];
        }
    }
    
    
    
    for (int i=0; i<Membrane_num_of_Triangle_Pairs ; i++)
    {
        for (int j=0; j<2 ; j++)
        {
            restart>> bondslist[i][j];
        }
    }
    
    
    for (int i=0; i<Actin_num_of_Nodes ; i++)
    {
        for (int j=0; j<3 ; j++)
        {
            restart>> Actin_Node_Position[i][j];
        }
    }
    
    
    
    for (int i=0; i<Actin_num_of_Nodes ; i++)
    {
        for (int j=0; j<3 ; j++)
        {
            restart>> Actin_Node_Velocity[i][j];
        }
    }
    
    
    for (int i=0; i<Actin_num_of_Bonds ; i++)
    {
        for (int j=0; j<3 ; j++)
        {
            restart>> Actin_Node_Pair_List[i][j];
        }
    }
    
    for (int i=0; i<nsolvent ; i++)
    {
        for (int j=0; j<3 ; j++)
        {
            restart>> xsolvent[j][i];
        }
    }
    
    
    for (int i=0; i<nsolvent ; i++)
    {
        for (int j=0; j<3 ; j++)
        {
            restart>> vsolvent[j][i];
        }
    }
    
    
    
    for (int i=0; i<Chromatin_num_of_Beads ; i++)
    {
        for (int j=0; j<3 ; j++)
        {
            restart>> Chromatin_Bead_Position[i][j];
        }
    }
    
    
    for (int i=0; i<Chromatin_num_of_Beads ; i++)
    {
        for (int j=0; j<3 ; j++)
        {
            restart>> Chromatin_Bead_Velocity[i][j];
        }
    }
    
    for (int i=0; i<ECM_num_of_Nodes ; i++)
    {
        for (int j=0; j<3 ; j++)
        {
            restart>> ECM_Node_Position[i][j];
        }
    }
    
    
    
    for (int i=0; i<ECM_num_of_Nodes ; i++)
    {
        for (int j=0; j<3 ; j++)
        {
            restart>> ECM_Node_Velocity[i][j];
        }
    }
    
    
}





//---------------------------restart-----------------------



void povray_output_creator(int currentStep, double Membrane_Node_Position[][3], int  Membrane_triangle_list [Membrane_num_of_Triangles][3], int Membrane_Normal_direction[Membrane_num_of_Triangles][2], int Membrane_Node_Pair_list[][2], double Actin_Node_Position[Actin_num_of_Nodes][3], double Actin_Node_Pair_List[][3], double Chromatin_Bead_Position[Chromatin_num_of_Beads][3], double ECM_Node_Position[][3], double ECM_Node_Pair_List[][3], int ECM_surface_triangle_list[ECM_Surface_num_of_Triangles][3], int Outer_Membrane_num_of_triangles, int Membrane_num_of_Node_Pairs, int Outer_Membrane_num_of_Node_Pairs, int Actin_num_of_Bonds, int ECM_num_of_Bonds, int Membrane_num_of_Nodes)
{
    
    ///____________________________POV
    ///____________________________POV
    ///____________________________POV
    ///____________________________POV
    ///____________________________POV
    
    
    
    double normalvectorofnodes[Membrane_num_of_Nodes][3];
    for(int i=0;i<Membrane_num_of_Nodes;i++)
    {
        normalvectorofnodes[i][0]=0.0;
        normalvectorofnodes[i][1]=0.0;
        normalvectorofnodes[i][2]=0.0;
    }
    double normal[3]; // normal vector of membrane
    double a1a2[3],a1a3[3]; // normal vector of membrane
    for(int i=0;i<Membrane_num_of_Triangles;i++)  //this loop caclulate the avrage normal vector of each node
    {
        a1a2[0]=Membrane_Node_Position[ Membrane_triangle_list[i][1]][0]-Membrane_Node_Position[ Membrane_triangle_list[i][0]][0];
        a1a2[1]=Membrane_Node_Position[ Membrane_triangle_list[i][1]][1]-Membrane_Node_Position[ Membrane_triangle_list[i][0]][1];
        a1a2[2]=Membrane_Node_Position[ Membrane_triangle_list[i][1]][2]-Membrane_Node_Position[ Membrane_triangle_list[i][0]][2];
        a1a3[0]=Membrane_Node_Position[ Membrane_triangle_list[i][2]][0]-Membrane_Node_Position[ Membrane_triangle_list[i][0]][0];
        a1a3[1]=Membrane_Node_Position[ Membrane_triangle_list[i][2]][1]-Membrane_Node_Position[ Membrane_triangle_list[i][0]][1];
        a1a3[2]=Membrane_Node_Position[ Membrane_triangle_list[i][2]][2]-Membrane_Node_Position[ Membrane_triangle_list[i][0]][2];
        crossvector(normal,a1a2,a1a3);
        normal[0]=normal[0]*Membrane_Normal_direction[i][1]/vectorlength(normal);
        normal[1]=normal[1]*Membrane_Normal_direction[i][1]/vectorlength(normal);
        normal[2]=normal[2]*Membrane_Normal_direction[i][1]/vectorlength(normal);  // it is the norml of the membrane element
        
        normalvectorofnodes [ Membrane_triangle_list[i][0] ][0]= normalvectorofnodes [ Membrane_triangle_list[i][0] ][0] + normal[0];
        normalvectorofnodes [ Membrane_triangle_list[i][0] ][1]= normalvectorofnodes [ Membrane_triangle_list[i][0] ][1] + normal[1];
        normalvectorofnodes [ Membrane_triangle_list[i][0] ][2]= normalvectorofnodes [ Membrane_triangle_list[i][0] ][2] + normal[2];
        
        
        normalvectorofnodes [ Membrane_triangle_list[i][1] ][0]= normalvectorofnodes [ Membrane_triangle_list[i][1] ][0] + normal[0];
        normalvectorofnodes [ Membrane_triangle_list[i][1] ][1]= normalvectorofnodes [ Membrane_triangle_list[i][1] ][1] + normal[1];
        normalvectorofnodes [ Membrane_triangle_list[i][1] ][2]= normalvectorofnodes [ Membrane_triangle_list[i][1] ][2] + normal[2];
        
        normalvectorofnodes [ Membrane_triangle_list[i][2] ][0]= normalvectorofnodes [ Membrane_triangle_list[i][2] ][0] + normal[0];
        normalvectorofnodes [ Membrane_triangle_list[i][2] ][1]= normalvectorofnodes [ Membrane_triangle_list[i][2] ][1] + normal[1];
        normalvectorofnodes [ Membrane_triangle_list[i][2] ][2]= normalvectorofnodes [ Membrane_triangle_list[i][2] ][2] + normal[2];
    }
    
    
    
    
    
    ///=============GENRAL===============
    int w1,w2;
    //////________
    string pov_file_name;
    pov_file_name="results/zpov_of_step_";
    pov_file_name=pov_file_name+to_string((currentStep));
    pov_file_name=pov_file_name+".pov";
    
    ofstream pov;
    pov.open(pov_file_name.c_str() );
    pov << std:: fixed;
    
    /////_________
    
    
    
    pov<< "global_settings {   ambient_light rgb<1, 1, 1>   photons {  spacing 0.2  autostop 1  media 60  max_trace_level 6 } } "<<endl;
    pov<< " #include \"colors.inc\"   " <<endl;
    pov<< " #include \"textures.inc\"   " <<endl;
    pov<< " #include \"metals.inc\"   " <<endl;
    pov<< " #include \"glass.inc\"   " <<endl;
    pov<< " #default{ finish{ ambient 0.1 diffuse 0.9 }}   " <<endl;
    
    pov<< "  camera{  location<50,50,50>  look_at <-4,-20,0>    right 0.5*4/3*x     up 0.5*y}"<<endl;
    pov << "background {White} light_source { <000, 100, 000> color White }   light_source { <100, 100, 100> color White }   light_source { <-100, 100, 100> color White } light_source { <00, 100, 100> color White } light_source { <00, 200, 00> color White }"<<endl;
    
    
    
    pov<< "#declare membrane = texture {    pigment{color rgbt<1.0,0.1,0.10,0.3>   }   finish {     diffuse 1.2 brilliance  70.0   ambient 0.30} } " <<endl;
    
    pov<< " #declare nucleus =  texture {   pigment{color rgbt<0.50,0.1,0.10,0.81>   }   finish { diffuse 0.6    ambient 3.0 } } " <<endl;
    
    pov<< "#declare chromatins = texture {    pigment{color rgbt<0.0,0.10,1.0,0.0>   }    finish {   diffuse 10.0  phong 2.1 phong_size 100  ambient 1.8     } }" <<endl;
    
    
    pov<< " #declare ECM = texture { pigment{color rgbt<10/255,10/255,10/255,0>  }   }" <<endl;
    
    
    pov<< " #declare chromatins =texture {    pigment{color rgbt<1.0,0.1,0.10,0.03>   }   finish {     diffuse 1.2 brilliance  70.0   ambient 0.30} }  " <<endl;
    
    pov<< " #declare membranebonds =  texture {    pigment{color rgbt<1.0,0.1,0.10,0.03>   }   finish {     diffuse 1.2 brilliance  70.0   ambient 0.30} }" <<endl;
    pov<< " #declare nucleusbonds =  texture {    pigment{color rgbt<1.0,0.1,0.10,0.03>   }   finish {     diffuse 1.2 brilliance  70.0   ambient 0.30} }" <<endl;
    pov<< "  #declare ECMbonds = texture {    pigment{color rgb<50/255,50/255,50/255>  }}  " <<endl;
    
    pov<< "#declare  radiusChromatin=" << Chromatin_Scaling_Factor*sigmachromatin*0.5 <<";"<<endl;
    pov<< "#declare  radiBondMeM=0.03"  <<";"<<endl;
    pov<< "#declare  radiBondECM=0.1"  <<";"<<endl;
    
    pov << " " <<endl;
    pov << " " <<endl;
    pov << " " <<endl;
    pov << " " <<endl;
    pov << " " <<endl;
    pov << " " <<endl;
    pov << " " <<endl;
    pov << " " <<endl;
    pov << " " <<endl;
    pov << " " <<endl;
    pov << " " <<endl;
    pov << " " <<endl;
    pov << " " <<endl;
    pov << " " <<endl;
    pov << " " <<endl;
    
    
    ///=============MEMRANE===============
    pov << "union {   // MEMBRANE____________________________________________________________________" <<endl;
    for(int i=0;i<Outer_Membrane_num_of_triangles;i++)
    {
        
        pov<< " triangle { < "<<Membrane_Node_Position[Membrane_triangle_list[i][0]][0]<<","<<Membrane_Node_Position[Membrane_triangle_list[i][0]][1] << ","<<Membrane_Node_Position[Membrane_triangle_list[i][0]][2]<< ">,";
        // pov<< " <"<<normalvectorofnodes[Membrane_triangle_list[0][i]][0]<<","<<normalvectorofnodes[Membrane_triangle_list[0][i]][1] << ","<<normalvectorofnodes[Membrane_triangle_list[0][i]][2]<< ">,";
        pov<< " <"<<Membrane_Node_Position[Membrane_triangle_list[i][1]][0]<<","<<Membrane_Node_Position[Membrane_triangle_list[i][1]][1] << ","<<Membrane_Node_Position[Membrane_triangle_list[i][1]][2]<< ">,";
        // pov<< " <"<<normalvectorofnodes[Membrane_triangle_list[1][i]][0]<<","<<normalvectorofnodes[Membrane_triangle_list[1][i]][1] << ","<<normalvectorofnodes[Membrane_triangle_list[1][i]][2]<< ">,";
        pov<< " <"<<Membrane_Node_Position[Membrane_triangle_list[i][2]][0]<<","<<Membrane_Node_Position[Membrane_triangle_list[i][2]][1] << ","<<Membrane_Node_Position[Membrane_triangle_list[i][2]][2]<< ">";
        // pov<< " <"<<normalvectorofnodes[Membrane_triangle_list[2][i]][0]<<","<<normalvectorofnodes[Membrane_triangle_list[2][i]][1] << ","<<normalvectorofnodes[Membrane_triangle_list[2][i]][2]<< ">";
        pov<< " } " <<endl;
    }
    pov<< " texture {membrane} no_shadow " <<endl;
    pov << " }" <<endl;
    pov << " " <<endl;
    pov << " " <<endl;
    pov << " " <<endl;
    pov << " " <<endl;
    pov << " " <<endl;
    
    
    
    
    
    
    pov << "union {   // nucleus___________________________________________________________________" <<endl;
    for(int i=Outer_Membrane_num_of_triangles;i<Membrane_num_of_Triangles;i++)
    {
        pov<< " smooth_triangle { < "<<Membrane_Node_Position[Membrane_triangle_list[i][0]][0]<<","<<Membrane_Node_Position[Membrane_triangle_list[i][0]][1] << ","<<Membrane_Node_Position[Membrane_triangle_list[i][0]][2]<< ">,";
        pov<< " <"<<normalvectorofnodes[Membrane_triangle_list[i][0]][0]<<","<<normalvectorofnodes[Membrane_triangle_list[i][0]][1] << ","<<normalvectorofnodes[Membrane_triangle_list[i][0]][2]<< ">,";
        pov<< " <"<<Membrane_Node_Position[Membrane_triangle_list[i][1]][0]<<","<<Membrane_Node_Position[Membrane_triangle_list[i][1]][1] << ","<<Membrane_Node_Position[Membrane_triangle_list[i][1]][2]<< ">,";
        pov<< " <"<<normalvectorofnodes[Membrane_triangle_list[i][1]][0]<<","<<normalvectorofnodes[Membrane_triangle_list[i][1]][1] << ","<<normalvectorofnodes[Membrane_triangle_list[i][1]][2]<< ">,";
        pov<< " <"<<Membrane_Node_Position[Membrane_triangle_list[i][2]][0]<<","<<Membrane_Node_Position[Membrane_triangle_list[i][2]][1] << ","<<Membrane_Node_Position[Membrane_triangle_list[i][2]][2]<< ">,";
        pov<< " <"<<normalvectorofnodes[Membrane_triangle_list[i][2]][0]<<","<<normalvectorofnodes[Membrane_triangle_list[i][2]][1] << ","<<normalvectorofnodes[Membrane_triangle_list[i][2]][2]<< ">";
        pov<< " } " <<endl;
    }
    pov<< " texture {nucleus} no_shadow " <<endl;
    pov << " }" <<endl;
    pov << " " <<endl;
    pov << " " <<endl;
    pov << " " <<endl;
    pov << " " <<endl;
    pov << " " <<endl;
    
    
    pov << "union { ///=========================================membrane bonds======================================="<<endl;
    for(int i=0;i<Outer_Membrane_num_of_Node_Pairs;i++)
    {
        
        w1=(int)  Membrane_Node_Pair_list[i][0];
        w2=(int)  Membrane_Node_Pair_list[i][1];
        
        pov<< "cylinder { <" << Membrane_Node_Position[w1][0]<< ","<< Membrane_Node_Position[w1][1]<< "," << Membrane_Node_Position[w1][2]<< "> ," ;
        pov<< "  <" << Membrane_Node_Position[w2][0]<< ","<< Membrane_Node_Position[w2][1]<< ","<< Membrane_Node_Position[w2][2]<<"> ," <<"radiBondMeM";
        pov<< "} " <<endl;
        
        
    }
    
    pov<< " texture {membranebonds} no_shadow " <<endl;
    pov << " }" <<endl;
    
    pov << "union { ///=========================================nucleus bonds======================================="<<endl;
    for(int i=Outer_Membrane_num_of_Node_Pairs;i<Membrane_num_of_Node_Pairs;i++)
    {
        
        w1=(int)  Membrane_Node_Pair_list[i][0];
        w2=(int)  Membrane_Node_Pair_list[i][1];
        
        pov<< "cylinder { <" << Membrane_Node_Position[w1][0]<< ","<< Membrane_Node_Position[w1][1]<< "," << Membrane_Node_Position[w1][2]<< "> ," ;
        pov<< "  <" << Membrane_Node_Position[w2][0]<< ","<< Membrane_Node_Position[w2][1]<< ","<< Membrane_Node_Position[w2][2]<<"> ," <<"radiBondMeM" ;
        pov<< "} " <<endl;
        
        
    }
    
    pov<< " texture {nucleusbonds} no_shadow " <<endl;
    pov << " }" <<endl;
    
    
    ///============membrane bonds================
    
    pov << "union {   // ECM___________________________________________________________________" <<endl;
    for(int i=0;i<ECM_Surface_num_of_Triangles;i++)
    {
        pov << " triangle { <"<< ECM_Node_Position[ ECM_surface_triangle_list[i][0] ] [0]  << " , "<< ECM_Node_Position[ ECM_surface_triangle_list[i][0] ] [1] << " , "<< ECM_Node_Position[ ECM_surface_triangle_list[i][0] ] [2]  << ">,";
        pov << " <"<< ECM_Node_Position[ ECM_surface_triangle_list[i][1] ] [0]  << " , "<< ECM_Node_Position[ ECM_surface_triangle_list[i][1] ] [1] << " , "<< ECM_Node_Position[ ECM_surface_triangle_list[i][1] ] [2]  << ">,";
        pov << " <"<< ECM_Node_Position[ ECM_surface_triangle_list[i][2] ] [0]  << " , "<< ECM_Node_Position[ ECM_surface_triangle_list[i][2] ] [1] << " , "<< ECM_Node_Position[ ECM_surface_triangle_list[i][2] ] [2]  << ">" <<endl;
        pov<< "} " <<endl;
    }
    pov<< " texture {ECM}  " <<endl;
    
    pov << " }" <<endl;
    
    
    
    pov << "union { ///=========================================chromatins======================================="<<endl;
    for(int i=0;i<Chromatin_num_of_Beads;i++)
    {
        pov<< "sphere { <" << Chromatin_Bead_Position[i][0]<< ","<< Chromatin_Bead_Position[i][1]<< ","<< Chromatin_Bead_Position[i][2]<< "> ," << "radiusChromatin";
        pov<< "} " <<endl;
        
    }
    
    for (int nchain=0;nchain<Chromatin_num_of_chains;nchain++)
    {
        for(int i=nchain*(Chromatin_num_of_Beads/Chromatin_num_of_chains)  ;i< (nchain+1)*(Chromatin_num_of_Beads/Chromatin_num_of_chains) -1; i++ )  // all beads interaction whit the next one
        {
            pov<< "cylinder { <" << Chromatin_Bead_Position[i][0]<< ","<< Chromatin_Bead_Position[i][1]<< ","<< Chromatin_Bead_Position[i][2]<< "> ," ;
            pov<< "  <" << Chromatin_Bead_Position[i+1][0]<< ","<< Chromatin_Bead_Position[i+1][1]<< ","<< Chromatin_Bead_Position[i+1][2]<< "> ," <<"radiusChromatin";
            pov<< "} " <<endl;
        }
    }
    
    pov<< " texture {chromatins}  no_shadow" <<endl;
    pov << " }" <<endl;
    
    ///=========================================chromatins=======================================
    
    
    ///=========================================actin=================================
    
    pov << "union { ///=========================================actin======================================="<<endl;
    for(int i=0;i<Actin_num_of_Bonds;i++)
    {
        
        w1=(int)  Actin_Node_Pair_List[i][0];
        w2=(int)  Actin_Node_Pair_List[i][1];
        
        if( (w1<Actin_Membrane_shared_num_of_Nodes & w2>Actin_Membrane_shared_num_of_Nodes) ||  (w2<Actin_Membrane_shared_num_of_Nodes & w1>Actin_Membrane_shared_num_of_Nodes)  )
        {
            pov<< "cylinder { <" << Actin_Node_Position[w1][0]<< ","<< Actin_Node_Position[w1][1]<< ","<< Actin_Node_Position[w1][2]<< "> ," ;
            pov<< "  <" << Actin_Node_Position[w2][0]<< ","<<Actin_Node_Position[w2][1]<< ","<<Actin_Node_Position[w2][2]<< "> ," <<0.05 ;
            pov<< "} " <<endl;
        }
        
    }
    
    pov<< " pigment{color rgb<50/255,160/255,160/255>  transmit 1.0 }  " <<endl;
    pov << " finish {   diffuse 0.9   }no_shadow" <<endl;
    pov << " }" <<endl;
    ///=========================================actin=================================
    
    ///=========================================ECM network=================================
    
    pov << "union { ///=========================================ECM bonds======================================="<<endl;
    for(int i=0;i<ECM_num_of_Bonds;i++)
    {
        
        w1=(int)  ECM_Node_Pair_List[i][0];
        w2=(int)  ECM_Node_Pair_List[i][1];
        
        pov<< "cylinder { <" << ECM_Node_Position[w1][0]<< ","<< ECM_Node_Position[w1][1]<< ","<< ECM_Node_Position[w1][2]<< "> ," ;
        pov<< "  <" << ECM_Node_Position[w2][0]<< ","<<ECM_Node_Position[w2][1]<< ","<<ECM_Node_Position[w2][2]<< "> ," << "radiBondECM";
        pov<< "} " <<endl;
        
        
    }
    
    pov<< " texture {ECMbonds} no_shadow " <<endl;
    pov << " }" <<endl;
    
    ///=========================================ECM network=================================
    
    ///=============finalize===============
    ///____________________________POV
    ///____________________________POV
    ///____________________________POV
    ///____________________________POV
    ///____________________________POV
    
}

void generate_initial_condition_report (string initial_condition_file_name, int Membrane_num_of_Nodes){
    ofstream write_initial_condition;
    write_initial_condition.open(initial_condition_file_name.c_str() );
    
    write_initial_condition<<"MD_num_of_steps"<<"=\t"<<MD_num_of_steps<<endl;
    write_initial_condition<<"savingstep"<<"=\t"<<savingstep<<endl;
    write_initial_condition<<"MD_Time_Step"<<"=\t"<<MD_Time_Step<<endl;
    write_initial_condition<<"KT"<<"=\t"<<KT<<endl;
    write_initial_condition<<"RunThermostatePerstep"<<"=\t"<<RunThermostatePerstep<<endl;
    write_initial_condition<<"\n\n";
    write_initial_condition<<"Membrane_num_of_Nodes"<<"=\t"<<Membrane_num_of_Nodes<<endl;
    write_initial_condition<<"Membrane_num_of_Triangles"<<"=\t"<<Membrane_num_of_Triangles<<endl;
    write_initial_condition<<"Membrane_Radius"<<"=\t"<<Membrane_Radius<<endl;
    write_initial_condition<<"Nucleus_Membrane_radius"<<"=\t"<<Nucleus_Membrane_radius<<endl;
    write_initial_condition<<"K_surfaceConstant_local"<<"=\t"<<K_surfaceConstant_local<<endl;
    write_initial_condition<<"Membrane_spring_coefficient"<<"=\t"<<Membrane_spring_coefficient<<endl;
    write_initial_condition<<"Membrane_bending_coefficient"<<"=\t"<<Membrane_bending_coefficient<<endl;
    write_initial_condition<<"membrane_damping_coefficient"<<"=\t"<<membrane_damping_coefficient<<endl;
    write_initial_condition<<"Membrane_Node_Mass"<<"=\t"<<Membrane_Node_Mass<<endl;
    write_initial_condition<<"fluidity"<<"=\t"<<fluidity<<endl;
    write_initial_condition<<"Nucleus_Membrane_Radius_of_Hard_Sphere_Interaction"<<"=\t"<<Nucleus_Membrane_Radius_of_Hard_Sphere_Interaction<<endl;
    write_initial_condition<<"\n\n";
    write_initial_condition<<"k_actine_membrane"<<"=\t"<<k_actine_membrane<<endl;
    write_initial_condition<<"Actin_num_of_Nodes"<<"=\t"<<Actin_num_of_Nodes<<endl;
    write_initial_condition<<"Actin_Membrane_shared_num_of_Nodes"<<"=\t"<<Actin_Membrane_shared_num_of_Nodes<<endl;
    write_initial_condition<<"Actin_spring_coefficient"<<"=\t"<<Actin_spring_coefficient<<endl;
    write_initial_condition<<"Membrane_barrier_calculation_rate"<<"=\t"<<Membrane_barrier_calculation_rate<<endl;
    write_initial_condition<<"CytoskeletonNetworkType"<<"=\t"<<CytoskeletonNetworkType<<endl;
    write_initial_condition<<"Actin_Passive_Cable_Network_Coefficient"<<"=\t"<<Actin_Passive_Cable_Network_Coefficient<<endl;
    write_initial_condition<<"KActinACN_EA"<<"=\t"<<KActinACN_EA<<endl;
    write_initial_condition<<"ACN_TL0"<<"=\t"<<ACN_TL0<<endl;
    write_initial_condition<<"ACN_LC"<<"=\t"<<ACN_LC<<endl;
    write_initial_condition<<"Actin_kelvin_damping_coefficient"<<"=\t"<<Actin_kelvin_damping_coefficient<<endl;
    write_initial_condition<<"Actin_damping_Coefficient"<<"=\t"<<Actin_damping_Coefficient<<endl;
    write_initial_condition<<"Actin_Node_Mass"<<"=\t"<<Actin_Node_Mass<<endl;
    write_initial_condition<<"minimumlengthActin"<<"=\t"<<minimumlengthActin<<endl;
    write_initial_condition<<"maximumlengthActin"<<"=\t"<<maximumlengthActin<<endl;
    write_initial_condition<<"Actin_Membrane_Radius_of_Hard_Sphere_Interaction"<<"=\t"<<Actin_Membrane_Radius_of_Hard_Sphere_Interaction<<endl;
    write_initial_condition<<"\n\n";
    write_initial_condition<<"Chromatin_num_of_Beads"<<"=\t"<<Chromatin_num_of_Beads<<endl;
    write_initial_condition<<"Chromatin_num_of_chains"<<"=\t"<<Chromatin_num_of_chains<<endl;
    write_initial_condition<<"Chromatin_Bead_Mass"<<"=\t"<<Chromatin_Bead_Mass<<endl;
    write_initial_condition<<"Chromatin_spring_coefficient"<<"=\t"<<Chromatin_spring_coefficient<<endl;
    write_initial_condition<<"Chromatin_bending_coefficient"<<"=\t"<<Chromatin_bending_coefficient<<endl;
    write_initial_condition<<"sigmachromatin"<<"=\t"<<sigmachromatin<<endl;
    write_initial_condition<<"Rmaxchromatin"<<"=\t"<<Rmaxchromatin<<endl;
    write_initial_condition<<"Chromatin_Membrane_Radius_of_Hard_Sphere_Interaction"<<"=\t"<<Chromatin_Membrane_Radius_of_Hard_Sphere_Interaction<<endl;
    write_initial_condition<<"Chromatin_Scaling_Factor"<<"=\t"<<Chromatin_Scaling_Factor<<endl;
    write_initial_condition<<"ThermostatOnChromatin"<<"=\t"<<ThermostatOnChromatin<<endl;
    write_initial_condition<<"chromatin_force_cut_off"<<"=\t"<<chromatin_force_cut_off<<endl;
    write_initial_condition<<"chromatin_force_cut_off"<<"=\t"<<chromatin_force_cut_off<<endl;
    write_initial_condition<<"\n\n";
    write_initial_condition<<"ECM_num_of_Nodes"<<"=\t"<<ECM_num_of_Nodes<<endl;
    write_initial_condition<<"ECM_Surface_num_of_Triangles"<<"=\t"<<ECM_Surface_num_of_Triangles<<endl;
    write_initial_condition<<"AdhesiveIslandOffOrOn"<<"=\t"<<AdhesiveIslandOffOrOn<<endl;
    write_initial_condition<<"NumberOfAdhesiveIslandNodes"<<"=\t"<<NumberOfAdhesiveIslandNodes<<endl;
    write_initial_condition<<"ECM_Thickness"<<"=\t"<<ECM_Thickness<<endl;
    write_initial_condition<<"ECM_Flexibility"<<"=\t"<<ECM_Flexibility<<endl;
    write_initial_condition<<"sigmaECM"<<"=\t"<<sigmaECM<<endl;
    write_initial_condition<<"ECM_LJ_just_Repultion"<<"=\t"<<ECM_LJ_just_Repultion<<endl;
    write_initial_condition<<"epsilonECM"<<"=\t"<<epsilonECM<<endl;
    write_initial_condition<<"ECM_kelvin_damping_coefficient"<<"=\t"<<ECM_kelvin_damping_coefficient<<endl;
    write_initial_condition<<"miuDampECM"<<"=\t"<<miuDampECM<<endl;
    write_initial_condition<<"fECMcuttoff"<<"=\t"<<fECMcuttoff<<endl;
    write_initial_condition<<"ECM_Node_Mass"<<"=\t"<<ECM_Node_Mass<<endl;
    write_initial_condition<<"ECM_Min_Gradient_Coefficient"<<"=\t"<<ECM_Min_Gradient_Coefficient<<endl;
    write_initial_condition<<"ECM_Max_Gradient_Coefficient"<<"=\t"<<ECM_Max_Gradient_Coefficient<<endl;
    write_initial_condition<<"ECM_Gradient_length"<<"=\t"<<ECM_Gradient_length<<endl;
    write_initial_condition<<"spreading_force_magnitude"<<"=\t"<<spreading_force_magnitude<<endl;
    write_initial_condition<<"spreading_force_min_range"<<"=\t"<<spreading_force_min_range<<endl;
    write_initial_condition<<"spreading_force_max_range"<<"=\t"<<spreading_force_max_range<<endl;
    write_initial_condition<<"spreading_flag"<<"=\t"<<spreading_flag<<endl;
    write_initial_condition<<"spreading_force_cos_triangle_interaction_angle"<<"=\t"<<spreading_force_cos_triangle_interaction_angle<<endl;
    write_initial_condition<<"\n\n";
    write_initial_condition<<"Actin_membrane_stiff_spring_coefficient"<<"=\t"<<Actin_membrane_stiff_spring_coefficient<<endl;
    write_initial_condition<<"Actin_membrane_damping_coefficient"<<"=\t"<<Actin_membrane_damping_coefficient<<endl;
    
}





