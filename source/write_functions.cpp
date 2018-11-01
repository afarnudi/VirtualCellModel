//
//  write_functions.cpp
//  Mem
//
//  Created by Ali Farnudi on 19/08/2018.
//  Copyright © 2018 Ali Farnudi. All rights reserved.
//

#include "write_functions.hpp"

void Results (ECM ecm, string label, char* buffer)
{
    //output file names:
    
    string energy_file_name;
    string traj_file_name;
    ecm.output_file_neme+=buffer;
    traj_file_name="Results_";
    traj_file_name+=ecm.output_file_neme;
    traj_file_name+=".xyz";
    energy_file_name="Results_Potential_Energy";
    energy_file_name+=ecm.output_file_neme;
    
    //trajectory:
    
    ofstream Trajectory;
    
    Trajectory.open(traj_file_name.c_str(), ios::app);
    Trajectory << std:: fixed;
    Trajectory << ecm.return_num_of_nodes()<<endl;
    Trajectory << " nodes  "<<endl;
    for(int j=0; j< ecm.return_num_of_nodes();j++) // saving trajectory
    {
        Trajectory << label <<setprecision(5)<< setw(20)<<ecm.Node_Position[j][0]<< setw(20)<<ecm.Node_Position[j][1]<< setw(20)<<ecm.Node_Position[j][2]<<endl;
    }
    
    //    //Energy:
    //    ofstream Potential_Energy;
    //    Potential_Energy.open(energy_file_name.c_str(), ios::app);
    //    Potential_Energy<<membrane.Total_Potential_Energy<<endl;
}

void Results (Membrane membrane, string label, char* buffer)
{
    //output file names:
    
    string energy_file_name;
    string traj_file_name;

    membrane.output_file_neme+=buffer;
    traj_file_name="Results_";
    traj_file_name+=membrane.output_file_neme;
    traj_file_name+=".xyz";
    energy_file_name="Results_Potential_Energy";
    energy_file_name+=membrane.output_file_neme;
	 
     
    
    //trajectory:
    
    ofstream Trajectory;
	
    Trajectory.open(traj_file_name.c_str(), ios::app);
    Trajectory << std:: fixed;
    Trajectory << membrane.return_num_of_nodes()<<endl;
    Trajectory << " nodes  "<<endl;
    for(int j=0; j< membrane.return_num_of_nodes();j++) // saving trajectory
    {
        Trajectory << label <<setprecision(5)<< setw(20)<<membrane.Node_Position[j][0]<< setw(20)<<membrane.Node_Position[j][1]<< setw(20)<<membrane.Node_Position[j][2]<<endl;
    }
	
	
    
	
    //    //Energy:
    //    ofstream Potential_Energy;
    //    Potential_Energy.open(energy_file_name.c_str(), ios::app);
    //    Potential_Energy<<membrane.Total_Potential_Energy<<endl;
}

void MembraneGeneratingReport (char* buffer, Membrane membrane )
{
	string Report_file_name;
	Report_file_name= "Membrane_Report_";
	Report_file_name+=buffer;
    Report_file_name+=".txt";
	
	ofstream Report;
	Report.open(Report_file_name.c_str());
	Report<< std:: fixed;
	Report<<"***General Constants***"<<endl;
//    Report<<"Number of MD steps"<< setw(20)<<MD_num_of_steps<<endl;
//    Report<<"MD time step"<< setw(20)<<MD_Time_Step<<endl;
//    Report<<"KT"<< setw(20)<<MD_KT<<endl;
//    Report<<"MD time step"<< setw(20)<<MD_Time_Step<<endl;
	Report<<"***Membrane Properties***"<<endl;
	Report<<"Membrane Node Mass"<< setw(20)<<membrane.Node_Mass<<endl;
	Report<<"Membrane Radius"<< setw(20)<<membrane.Radius<<endl;
	Report<<"Minimum node pair length"<< setw(20)<<membrane.Min_node_pair_length<<endl;
	Report<<"Maximum node pair length"<< setw(20)<<membrane.Max_node_pair_length<<endl;
	Report<<"Average node pair length"<< setw(20)<<membrane.Average_node_pair_length<<endl;
	Report<<"number of Membrane's Nodes "<< setw(20)<<membrane.return_num_of_nodes()<<endl;
	Report<<"number of Membrane's Triangles "<< setw(20)<<membrane.return_num_of_triangle()<<endl;
	Report<<"number of Membrane's Nodes"<< setw(20)<<membrane.return_num_of_nodes()<<endl;
	Report<<"Membrane Spring coefficient"<< setw(20)<<membrane.Spring_coefficient<<endl;
	Report<<"Membrane Bending coefficient"<< setw(20)<<membrane.Bending_coefficient<<endl;
	Report<<"Membrane Damping coefficient"<< setw(20)<<membrane.Damping_coefficient<<endl;
	if (membrane.spring_model==1)
	{Report<<"Membrane Spring Model:"<< setw(20)<<"Standard Model"<<endl;}
	if (membrane.spring_model==2)
	{Report<<"Membrane Spring Model:"<< setw(20)<<"Houkian"<<endl;}
	
}

void ParticleGeneratingReport (char* buffer, Membrane particle )
{
	string Report_file_name;
	Report_file_name= "Particle_Report_";
	Report_file_name+=buffer;
    Report_file_name+=".txt";
	
	ofstream Report;
	Report.open(Report_file_name.c_str());
	Report<< std:: fixed;
	Report<<"***General Constants***"<<endl;
//    Report<<"Number of MD steps"<< setw(20)<<MD_num_of_steps<<endl;
//    Report<<"MD time step"<< setw(20)<<MD_Time_Step<<endl;
//    Report<<"KT"<< setw(20)<<MD_KT<<endl;
//    Report<<"MD time step"<< setw(20)<<MD_Time_Step<<endl;
	Report<<"***Particle Properties***"<<endl;
	Report<<"Particle Node Mass"<< setw(20)<<particle.Node_Mass<<endl;
	Report<<"Particle Radius"<< setw(20)<<particle.Radius<<endl;
	Report<<"number of Particle's Nodes "<< setw(20)<<particle.return_num_of_nodes()<<endl;
	Report<<"number of Particle's Triangles "<< setw(20)<<particle.return_num_of_triangle()<<endl;
	Report<<"number of Particle's Nodes"<< setw(20)<<particle.return_num_of_nodes()<<endl;
	Report<<"Particle Spring coefficient"<< setw(20)<<particle.Spring_coefficient<<endl;
	Report<<"Particle Bending coefficient"<< setw(20)<<particle.Bending_coefficient<<endl;
	Report<<"Particle Damping coefficient"<< setw(20)<<particle.Damping_coefficient<<endl;
	if (particle.spring_model==1)
	{Report<<"Particle Spring Model:"<< setw(20)<<"Standard Model"<<endl;}
	if (particle.spring_model==2)
	{Report<<"Membrane Spring Model:"<< setw(20)<<"Houkian"<<endl;}
}


void EcmGeneratingReport (char* buffer, ECM ecm )
{
	string Report_file_name;
	Report_file_name= "ECM_Report_";
    Report_file_name+=buffer;
    Report_file_name+=".txt";
	
	ofstream Report;
	Report.open(Report_file_name.c_str());
	Report<< std:: fixed;
	Report<<"***General Constants***"<<endl;
//    Report<<"Number of MD steps"<< setw(20)<<MD_num_of_steps<<endl;
//    Report<<"MD time step"<< setw(20)<<MD_Time_Step<<endl;
//    Report<<"KT"<< setw(20)<<MD_KT<<endl;
//    Report<<"MD time step"<< setw(20)<<MD_Time_Step<<endl;
	Report<<"***ECM Properties***"<<endl;
	// 
	//Report<<"number of ECM's Nodes "<< setw(20)<<ecm.Num_of_Nodes<<endl;
	//Report<<"number of ECM's Triangles "<< setw(20)<<ecm.Num_of_Triangles<<endl;
}

void checkingForce (Membrane membrane, int MD_Step, char* buffer)
{
	string check_force_name="check_force_";
	check_force_name+=buffer;
	ofstream check_force;
	check_force.open(check_force_name.c_str(), ios::app);
	check_force<< "MD step = "<<MD_Step<<endl;
	for(int j=0; j< membrane.return_num_of_nodes();j++) // saving Forces
    {
        check_force  <<setprecision(5)<< setw(20)<<membrane.Node_Force[j][0]<< setw(20)<<membrane.Node_Force[j][1]<< setw(20)<<membrane.Node_Force[j][2]<<endl;
    }
	check_force<< "_________________________________________________________"<<endl;
}

void Membrane::export_for_resume(char* buffer, int MD_step){
    ofstream write_resume_file;
    string resume_file_name="Resume_";
    resume_file_name+=buffer;
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
