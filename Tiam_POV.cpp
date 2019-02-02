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

//Povray------
# define export_povray_step_distance 10000  // clear!
void povray_output_creator(int currentStep, double Membrane_Node_Position[][3], int  Membrane_triangle_list[Membrane_num_of_Triangles][3], int Membrane_Normal_direction[Membrane_num_of_Triangles][2], int Membrane_Node_Pair_list[][2], double Actin_Node_Position[Actin_num_of_Nodes][3], double Actin_Node_Pair_List[][3], double Chromatin_Bead_Position[Chromatin_num_of_Beads][3], double  ECM_Node_Position[][3], double ECM_Node_Pair_List[][3], int ECM_surface_triangle_list[ECM_Surface_num_of_Triangles][3], int Outer_Membrane_num_of_triangles, int Membrane_num_of_Node_Pairs, int Outer_Membrane_num_of_Node_Pairs, int Actin_num_of_Bonds, int ECM_num_of_Bonds, int Membrane_num_of_Nodes);

int Outer_Membrane_num_of_Node_Pairs ; //Ali's Comment: In previous versions of the code, the membrane and the membrane nucleus where introduced using a single mesh file. So in the initialisation of the code we would seperate the nodes into the outer membrane and the nucleus. Tiam's comment: usefull in POV ray- modulate in void sortingbonds(int bondslist[][2],int tri[3][Membrane_num_of_Triangles])

int main{
    
    cout<<"Beginning the MD loop\n";
    for(int MD_Step=0 ;MD_Step<=MD_num_of_steps ; MD_Step++)
    {

        if (MD_Step%(export_povray_step_distance)==0    ) // for last steps cacculate contact matrix
        {
            povray_output_creator(MD_Step, Membrane_Node_Position, Membrane_triangle_list, Membrane_Normal_direction,  Membrane_Node_Pair_list, Actin_Node_Position, Actin_Node_Pair_List, Chromatin_Bead_Position,   ECM_Node_Position, ECM_Node_Pair_List, ECM_surface_triangle_list, Outer_Membrane_num_of_triangles, Membrane_num_of_Node_Pairs, Outer_Membrane_num_of_Node_Pairs, Actin_num_of_Bonds, ECM_num_of_Bonds, Membrane_num_of_Nodes);
        }
    }
    
    
    
    
    
    
    
    
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////
    
}
    //pov
    povray_output_creator( MD_num_of_steps, Membrane_Node_Position, Membrane_triangle_list,
                          Membrane_Normal_direction, Membrane_Node_Pair_list, Actin_Node_Position,Actin_Node_Pair_List, Chromatin_Bead_Position, ECM_Node_Position, ECM_Node_Pair_List, ECM_surface_triangle_list, Outer_Membrane_num_of_triangles, Membrane_num_of_Node_Pairs, Outer_Membrane_num_of_Node_Pairs, Actin_num_of_Bonds, ECM_num_of_Bonds, Membrane_num_of_Nodes);
    //pov
    


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
