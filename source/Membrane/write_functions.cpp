//
//  write_functions.cpp
//  Mem
//
//  Created by Ali Farnudi on 19/08/2018.
//  Copyright © 2018 Ali Farnudi. All rights reserved.
//

#include "Membrane.h"
#include "General_constants.h"
#include "General_functions.hpp"

using std::ofstream;


void Membrane:: write_pov_traj(std::string traj_name, std::string label, int currentStep){
     ///=============GENRAL===============
    int w1,w2;
    //////________
    std::string pov_file_name;
    pov_file_name=traj_name;
    pov_file_name=pov_file_name+"_of_step_";
    pov_file_name=pov_file_name+std::to_string((currentStep));
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
    
    //pov<< " #declare nucleus =  texture {   pigment{color rgbt<0.50,0.1,0.10,0.81>   }   finish { diffuse 0.6    ambient 3.0 } } " <<endl;
    
   // pov<< "#declare chromatins = texture {    pigment{color rgbt<0.0,0.10,1.0,0.0>   }    finish {   diffuse 10.0  phong 2.1 phong_size 100  ambient 1.8     } }" <<endl;
    
    
   // pov<< " #declare ECM = texture { pigment{color rgbt<10/255,10/255,10/255,0>  }   }" <<endl;
    
    
    //pov<< " #declare chromatins =texture {    pigment{color rgbt<1.0,0.1,0.10,0.03>   }   finish {     diffuse 1.2 brilliance  70.0   ambient 0.30} }  " <<endl;
    
    pov<< " #declare membranebonds =  texture {    pigment{color rgbt<1.0,0.1,0.10,0.03>   }   finish {     diffuse 1.2 brilliance  70.0   ambient 0.30} }" <<endl;
   // pov<< " #declare nucleusbonds =  texture {    pigment{color rgbt<1.0,0.1,0.10,0.03>   }   finish {     diffuse 1.2 brilliance  70.0   ambient 0.30} }" <<endl;
  //  pov<< "  #declare ECMbonds = texture {    pigment{color rgb<50/255,50/255,50/255>  }}  " <<endl;
    
   // pov<< "#declare  radiusChromatin=" << Chromatin_Scaling_Factor*sigmachromatin*0.5 <<";"<<endl;
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
    
    ///=============MEMBRANE===============
    
    pov << "union {   // MEMBRANE____________________________________________________________________" <<endl;
    
    for(int i=0;i<Num_of_Triangles;i++)
    {
        
        pov<< " triangle { < "<<Node_Position[Triangle_list[i][0]][0]<<","<<Node_Position[Triangle_list[i][0]][1] << ","<<Node_Position[Triangle_list[i][0]][2]<< ">,";
        pov<< " <"<<Node_Position[Triangle_list[i][1]][0]<<","<<Node_Position[Triangle_list[i][1]][1] << ","<<Node_Position[Triangle_list[i][1]][2]<< ">,";
        pov<< " <"<<Node_Position[Triangle_list[i][2]][0]<<","<<Node_Position[Triangle_list[i][2]][1] << ","<<Node_Position[Triangle_list[i][2]][2]<< ">";
        pov<< " } " <<endl;
    }
    pov<< " texture {membrane} no_shadow " <<endl;
    pov << " }" <<endl;
    pov << " " <<endl;
    pov << " " <<endl;
    pov << " " <<endl;
    pov << " " <<endl;
    pov << " " <<endl;
    pov << "union { ///=========================================membrane bonds======================================="<<endl;
    for(int i=0;i<Num_of_Node_Pairs;i++)
    {
        
        w1=(int)  Node_Bond_list[i][0];
        w2=(int)  Node_Bond_list[i][1];
        
        pov<< "cylinder { <" << Node_Position[w1][0]<< ","<< Node_Position[w1][1]<< "," << Node_Position[w1][2]<< "> ," ;
        pov<< "  <" << Node_Position[w2][0]<< ","<< Node_Position[w2][1]<< ","<< Node_Position[w2][2]<<"> ," <<"radiBondMeM";
        pov<< "} " <<endl;
        
        
    }
    
    pov<< " texture {membranebonds} no_shadow " <<endl;
    pov << " }" <<endl;
}



void Membrane::export_for_resume(int MD_step){
    ofstream write_resume_file;
    std::string resume_file_name="Results/Resumes/Resume_Membrane_"+std::to_string(index)+"_";
    resume_file_name+=file_time;
    resume_file_name+=".txt";
    write_resume_file.open(resume_file_name.c_str());
    
    write_resume_file<<MD_step<<endl;
    write_resume_file<<Num_of_Nodes<<endl;
    for (int i=0; i<Num_of_Nodes; i++) {
        write_resume_file<<Node_Position[i][0]<<"\t"<<Node_Position[i][1]<<"\t"<<Node_Position[i][2]<<"\n";
        write_resume_file<<Node_Velocity[i][0]<<"\t"<<Node_Velocity[i][1]<<"\t"<<Node_Velocity[i][2]<<"\n";
    }
    
    write_resume_file<<Num_of_Triangles<<endl;
    for (int i=0; i<Num_of_Triangles; i++) {
        write_resume_file<<Triangle_list[i][0]<<"\t"<<Triangle_list[i][1]<<"\t"<<Triangle_list[i][2]<<"\n";
    }
    
    write_resume_file<<Num_of_Node_Pairs<<endl;
    for (int i=0; i<Num_of_Node_Pairs; i++) {
        write_resume_file<<Node_Bond_list[i][0]<<"\t"<<Node_Bond_list[i][1]<<"\n";
    }
    //In the import function we should call the neighbour list constructor
    write_resume_file<<Num_of_Triangle_Pairs<<endl;
    for (int i=0; i<Num_of_Triangle_Pairs; i++) {
        write_resume_file<<Triangle_pair_list[i][0]<<"\t"<<Triangle_pair_list[i][1]<<"\n";
        write_resume_file<<Triangle_Pair_Nodes[i][0]<<"\t"<<Triangle_Pair_Nodes[i][1]<<"\t"<<Triangle_Pair_Nodes[i][2]<<"\t"<<Triangle_Pair_Nodes[i][3]<<"\t"<<"\n";
    }
    write_resume_file<<Max_node_pair_length<<"\t"<<Min_node_pair_length<<"\t"<<Average_node_pair_length<<endl;
}

void Membrane::export_for_resume(int MD_step, MyAtomInfo atoms[], int atom_count){
    ofstream write_resume_file;
    std::string resume_file_name="Results/Resumes/Resume_Membrane_"+std::to_string(index)+"_";
    resume_file_name+=file_time;
    resume_file_name+=".txt";
    write_resume_file.open(resume_file_name.c_str());
    
    write_resume_file<<MD_step<<endl;
    write_resume_file<<Num_of_Nodes<<endl;
    
    for (int i=atom_count; i<atom_count+Num_of_Nodes; i++) {
        Node_Position[i-atom_count][0] = atoms[i].posInNm[0];
        Node_Position[i-atom_count][1] = atoms[i].posInNm[1];
        Node_Position[i-atom_count][2] = atoms[i].posInNm[2];
        
        Node_Velocity[i-atom_count][0] = atoms[i].velocityInNmperPs[0];
        Node_Velocity[i-atom_count][1] = atoms[i].velocityInNmperPs[1];
        Node_Velocity[i-atom_count][2] = atoms[i].velocityInNmperPs[2];
    }
    
    for (int i=0; i<Num_of_Nodes; i++) {
        write_resume_file<<Node_Position[i][0]<<"\t"<<Node_Position[i][1]<<"\t"<<Node_Position[i][2]<<"\n";
        write_resume_file<<Node_Velocity[i][0]<<"\t"<<Node_Velocity[i][1]<<"\t"<<Node_Velocity[i][2]<<"\n";
    }
    
    write_resume_file<<Num_of_Triangles<<endl;
    for (int i=0; i<Num_of_Triangles; i++) {
        write_resume_file<<Triangle_list[i][0]<<"\t"<<Triangle_list[i][1]<<"\t"<<Triangle_list[i][2]<<"\n";
    }
    
    write_resume_file<<Num_of_Node_Pairs<<endl;
    for (int i=0; i<Num_of_Node_Pairs; i++) {
        write_resume_file<<Node_Bond_list[i][0]<<"\t"<<Node_Bond_list[i][1]<<"\n";
    }
    //In the import function we should call the neighbour list constructor
    write_resume_file<<Num_of_Triangle_Pairs<<endl;
    for (int i=0; i<Num_of_Triangle_Pairs; i++) {
        write_resume_file<<Triangle_pair_list[i][0]<<"\t"<<Triangle_pair_list[i][1]<<"\n";
        write_resume_file<<Triangle_Pair_Nodes[i][0]<<"\t"<<Triangle_Pair_Nodes[i][1]<<"\t"<<Triangle_Pair_Nodes[i][2]<<"\t"<<Triangle_Pair_Nodes[i][3]<<"\t"<<"\n";
    }
    write_resume_file<<Max_node_pair_length<<"\t"<<Min_node_pair_length<<"\t"<<Average_node_pair_length<<endl;
}


void Membrane::generate_report()
{
    std::string Report_file_name;
    Report_file_name= "Results/Reports/Report_Membrane_"+std::to_string(index)+"_";
    Report_file_name+=file_time;
    Report_file_name+=".txt";
    
    ofstream Report;
    Report.open(Report_file_name.c_str());
    
    if (!Report.is_open()){
        cout<<"Couldn't generate report. Please check that the required directories exist. The required directory should be ./bin/Results/Reports or ./Results/Reports"<<endl;
        Report.close();
        cout<<"Generating report in an alternative directory:\nbin/Results/Reports/"<<endl;
        Report_file_name= "bin/Results/Reports/Report_Membrane_"+std::to_string(index)+"_";
        Report_file_name+=file_time;
        Report_file_name+=".txt";
        Report.open(Report_file_name.c_str());
    }
    
    Report<< std:: fixed;
    Report<<"MD Params in the general configuration file:\n---------------\n";
    Report<<"Simulation_Time_In_Ps\t"<<GenConst::Simulation_Time_In_Ps<<endl;
    Report<<"Step_Size_In_Fs\t"<<GenConst::Step_Size_In_Fs<<endl;
    Report<<"Report_Interval_In_Fs\t"<<GenConst::Report_Interval_In_Fs<<endl;
    Report<<"temperature (kelvin)\t"<<GenConst::temperature<<endl;
    Report<<"MC_step\t"<<GenConst::MC_step<<endl;
    Report<<"Mem_fluidity\t"<<GenConst::Mem_fluidity<<endl;
//    Report<<"Lbox\t"<<GenConst::Lbox<<endl;
    Report<<"Periodic_box\t"<<GenConst::Periodic_box<<endl;
    Report<<"Num_of_Membranes\t"<<GenConst::Num_of_Membranes<<endl;
    Report<<"trajectory_file_name\t"<<GenConst::trajectory_file_name<<endl;
    
    
    Report<<"\nMembrane Params:\n---------------\n";
    Report<<"Mesh_file_name\t"<<Mesh_file_name<<endl;
    Report<<"Node_Mass\t"<<Node_Mass<<endl;
    Report<<"Node_radius\t"<<Node_radius<<endl;
    Report<<"Membrane Radius\t"<<Radius<<endl;
    Report<<"spring_model\t"<<spring_model;
    if (spring_model==1){
        Report<<"\tFENE logarithmic barrier"<<endl;
    } else if (spring_model==2) {
        Report<<"\tHookian"<<endl;
    } else if (spring_model==3) {
        Report<<"\tcustom"<<endl;
    }
    Report<<"Spring coefficient\t"<<Spring_coefficient<<endl;
    Report<<"Bending coefficient\t"<<Bending_coefficient<<endl;
    Report<<"Damping coefficient\t"<<Damping_coefficient<<endl;
    Report<<"Shift_in_X_direction\t"<<Shift_in_X_direction<<endl;
    Report<<"Shift_in_Y_direction\t"<<Shift_in_Y_direction<<endl;
    Report<<"Shift_in_Z_direction\t"<<Shift_in_Z_direction<<endl;
    Report<<"x_speed\t"<<x_speed<<endl;
    Report<<"y_speed\t"<<y_speed<<endl;
    Report<<"z_speed\t"<<z_speed<<endl;
    Report<<"X_in_mem\t"<<X_in<<endl;
    Report<<"Y_in_mem\t"<<Y_in<<endl;
    Report<<"Z_in_mem\t"<<Z_in<<endl;
    Report<<"position_scale_x\t"<<X_scale<<endl;
    Report<<"position_scale_y\t"<<Y_scale<<endl;
    Report<<"position_scale_z\t"<<Z_scale<<endl;
    Report<<"rescale_factor\t"<<rescale_factor<<endl;
    Report<<"sigma_LJ_12_6\t"<<sigma_LJ_12_6<<endl;
    Report<<"epsilon_LJ_12_6\t"<<epsilon_LJ_12_6<<endl;
    
    Report<<"Minimum node pair length\t"<<Min_node_pair_length<<endl;
    Report<<"Maximum node pair length\t"<<Max_node_pair_length<<endl;
    Report<<"Average node pair length\t"<<Average_node_pair_length<<endl;
    Report<<"# of Nodes\t"<<get_num_of_nodes()<<endl;
    Report<<"# of Triangles\t"<<get_num_of_triangle()<<endl;
    Report<<"# of triangle pairs="<<Num_of_Triangle_Pairs<<endl;
    
    Report.close();
    
    
    
}
     
