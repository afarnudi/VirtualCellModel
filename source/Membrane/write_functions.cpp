//
//  write_functions.cpp
//  Mem
//
//  Created by Ali Farnudi on 19/08/2018.
//  Copyright © 2018 Ali Farnudi. All rights reserved.
//

#include "Membrane.h"
#include "General_constants.h"

void Membrane::write_traj (string traj_name, string label){
    ofstream Trajectory;
    Trajectory.open(traj_name.c_str(), ios::app);
    Trajectory << std:: fixed;
    for(int j=0; j< Num_of_Nodes;j++) // saving trajectory
    {
        Trajectory << label <<setprecision(5)<< setw(20)<<Node_Position[j][0]<< setw(20)<<Node_Position[j][1]<< setw(20)<<Node_Position[j][2]<<endl;
    }
}
void Membrane:: write_pov_traj(string traj_name, string label, int currentStep){
     ///=============GENRAL===============
    int w1,w2;
    //////________
    string pov_file_name;
    pov_file_name=traj_name;
    pov_file_name=pov_file_name+"_of_step_";
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
    
   ;
    
    
    ///============membrane bonds================
}
void Membrane::relaxation_traj (void)
{
    string energy_file_name;
    string traj_file_name;
    
    traj_file_name="Results/Relaxation/Relaxation_"+GenConst::trajectory_file_name+"Membrane_"+to_string(index)+"_"+file_time+".xyz";
    //trajectory:
    
    ofstream Trajectory;
    
    Trajectory.open(traj_file_name.c_str(), ios::app);
    Trajectory << std:: fixed;
    Trajectory <<Num_of_Nodes<<endl;
    Trajectory << " nodes  "<<endl;
    for(int j=0; j< Num_of_Nodes; j++) // saving trajectory
    {
        Trajectory << "mem" <<setprecision(5)<< setw(20)<<Node_Position[j][0]<< setw(20)<<Node_Position[j][1]<< setw(20)<<Node_Position[j][2]<<endl;
    }
    
}


void Membrane::export_for_resume(int MD_step){
    ofstream write_resume_file;
    string resume_file_name="Results/Resumes/Resume_Membrane_"+to_string(index)+"_";
    resume_file_name+=file_time;
    resume_file_name+=".txt";
    write_resume_file.open(resume_file_name.c_str());
    
    write_resume_file<<MD_step<<endl;
    write_resume_file<<Num_of_Nodes<<endl;
    for (int i=0; i<Num_of_Nodes; i++) {
        write_resume_file<<Node_Position[i][0]<<"\t"<<Node_Position[i][1]<<"\t"<<Node_Position[i][2]<<"\n";
        write_resume_file<<Node_Velocity[i][0]<<"\t"<<Node_Velocity[i][1]<<"\t"<<Node_Velocity[i][2]<<"\n";
        //Node_force=0
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
    //run_check_ for max and min node distances
}


void Membrane::generate_report()
{
    string Report_file_name;
    Report_file_name= "Results/Reports/Report_Membrane_"+to_string(index)+"_";
    Report_file_name+=file_time;
    Report_file_name+=".txt";
    
    ofstream Report;
    Report.open(Report_file_name.c_str());
    Report<< std:: fixed;
    Report<<"General MD Params:\n---------------\n";
    Report<<"MD_num_of_steps"<<setw(20)<<GenConst::MD_num_of_steps<<endl;
    Report<<"MD_traj_save_step"<<setw(20)<<GenConst::MD_traj_save_step<<endl;
    Report<<"MD_Time_Step"<<setw(20)<<GenConst::MD_Time_Step<<endl;
    Report<<"MD_T"<<setw(20)<<GenConst::MD_T<<endl;
    Report<<"MD_thrmo_step"<<setw(20)<<GenConst::MD_thrmo_step<<endl;
    Report<<"Bussi_tau"<<setw(20)<<GenConst::Bussi_tau<<endl;
    Report<<"MC_step"<<setw(20)<<GenConst::MC_step<<endl;
    Report<<"Mem_fluidity"<<setw(20)<<GenConst::Mem_fluidity<<endl;
    Report<<"Lbox"<<setw(20)<<GenConst::Lbox<<endl;
    Report<<"Periodic_condtion_status"<<setw(20)<<GenConst::Periodic_condtion_status<<endl;
    Report<<"Num_of_Membranes"<<setw(20)<<GenConst::Num_of_Membranes<<endl;
    Report<<"trajectory_file_name"<<setw(20)<<GenConst::trajectory_file_name<<endl;
    
    
    Report<<"Membrane Params:\n---------------\n";
    Report<<"Node Mass"<< setw(20)<<Node_Mass<<endl;
    Report<<"Radius"<< setw(20)<<Radius<<endl;
    Report<<"spring_model"<< setw(20)<<spring_model<<endl;
    if (spring_model==1){
        Report<<"Membrane Spring Model:"<< setw(20)<<"logarithmic barrier"<<endl;
    } else if (spring_model==2) {
        Report<<"Membrane Spring Model:"<< setw(20)<<"Hookian"<<endl;
    } else if (spring_model==3) {
        Report<<"Membrane Spring Model:"<< setw(20)<<"FENE"<<endl;
    }
    Report<<"Spring coefficient"<< setw(20)<<Spring_coefficient<<endl;
    Report<<"Bending coefficient"<< setw(20)<<Bending_coefficient<<endl;
    Report<<"Damping coefficient"<< setw(20)<<Damping_coefficient<<endl;
    Report<<"K_surfaceConstant_local"<< setw(20)<<K_surfaceConstant_local<<endl;
    Report<<"Shift_in_X_direction"<< setw(20)<<Shift_in_X_direction<<endl;
    Report<<"Shift_in_Y_direction"<< setw(20)<<Shift_in_Y_direction<<endl;
    Report<<"Shift_in_Z_direction"<< setw(20)<<Shift_in_Z_direction<<endl;
    Report<<"Downward_speed"<< setw(20)<<Downward_speed<<endl;
    Report<<"X_in_mem"<< setw(20)<<X_in<<endl;
    Report<<"Y_in_mem"<< setw(20)<<Y_in<<endl;
    Report<<"Z_in_mem"<< setw(20)<<Z_in<<endl;
    Report<<"position_scale_x"<< setw(20)<<X_scale<<endl;
    Report<<"position_scale_y"<< setw(20)<<Y_scale<<endl;
    Report<<"position_scale_z"<< setw(20)<<Z_scale<<endl;
    
    Report<<"Minimum node pair length"<< setw(20)<<Min_node_pair_length<<endl;
    Report<<"Maximum node pair length"<< setw(20)<<Max_node_pair_length<<endl;
    Report<<"Average node pair length"<< setw(20)<<Average_node_pair_length<<endl;
    Report<<"# of Nodes "<< setw(20)<<return_num_of_nodes()<<endl;
    Report<<"# of Triangles "<< setw(20)<<return_num_of_triangle()<<endl;
    
    
    
    
    
}

void Membrane::write_parameters(int MD_Step){
    //    string energy_file_name;
    string traj_file_name;
//    omega_calculator_2();
    double a[3]={Omega[0],Omega[1],Omega[2]};
    double Omega_len=vector_length(a);
    
    traj_file_name="Results/Param_"+GenConst::trajectory_file_name+"Membrane_"+to_string(index)+"_"+file_time+".txt";
    ofstream Trajectory;
    
    Trajectory.open(traj_file_name.c_str(), ios::app);
    Trajectory << std:: fixed;
    if (GenConst::File_header==false) {
        Trajectory<<"MD Step\t"<<"Total Kinetic Energy\t"
        <<"omega_x\t"<<"omega_y\t"<<"omega_z\t"<<"omega\n";
        GenConst::File_header=true;
    }
    Trajectory<<MD_Step<<"\t"<<Total_Kinetic_Energy
    <<"\t"<<Omega[0]<<"\t"<<Omega[1]<<"\t"<<Omega[2]<<"\t"<<Omega_len
    <<endl;
}


void Membrane::export_relaxed(int MD_step){
    ofstream write_resume_file;
    string resume_file_name="Results/Relaxation/Resume_Membrane_"+to_string(index)+"_";
    resume_file_name+=file_time;
    resume_file_name+=".txt";
    write_resume_file.open(resume_file_name.c_str());
    
    write_resume_file<<MD_step<<endl;
    write_resume_file<<Num_of_Nodes<<endl;
    for (int i=0; i<Num_of_Nodes; i++) {
        write_resume_file<<Node_Position[i][0]<<"\t"<<Node_Position[i][1]<<"\t"<<Node_Position[i][2]<<"\n";
        write_resume_file<<Node_Velocity[i][0]<<"\t"<<Node_Velocity[i][1]<<"\t"<<Node_Velocity[i][2]<<"\n";
        //Node_force=0
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
    //run_check_ for max and min node distances
}
