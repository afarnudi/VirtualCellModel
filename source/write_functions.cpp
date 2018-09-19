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

void generatingReport (char* buffer, Membrane membrane )
{
	string Report_file_name;
	Report_file_name= "Report_";
	Report_file_name+=buffer;
	
	ofstream Report;
	Report.open(Report_file_name.c_str());
	Report<< std:: fixed;
	Report<<"***General Constants***"<<endl;
	Report<<"Number of MD steps"<< setw(20)<<MD_num_of_steps<<endl;
	Report<<"MD time step"<< setw(20)<<MD_Time_Step<<endl;
	Report<<"KT"<< setw(20)<<KT<<endl;
	Report<<"MD time step"<< setw(20)<<MD_Time_Step<<endl;
	Report<<"***Membrane Properties***"<<endl;
	Report<<"Membrane Node Mass"<< setw(20)<<membrane.Node_Mass<<endl;
	Report<<"Membrane Radius"<< setw(20)<<membrane.Radius<<endl;
	Report<<"number of Membrane's Nodes "<< setw(20)<<membrane.Num_of_Nodes<<endl;
	Report<<"number of Membrane's Triangles "<< setw(20)<<membrane.Num_of_Triangles<<endl;
	Report<<"number of Membrane's Nodes"<< setw(20)<<membrane.Num_of_Nodes<<endl;
	Report<<"Membrane Spring coefficient"<< setw(20)<<membrane.Spring_coefficient<<endl;
	Report<<"Membrane Bending coefficient"<< setw(20)<<membrane.Bending_coefficient<<endl;
	Report<<"Membrane Damping coefficient"<< setw(20)<<membrane.Damping_coefficient<<endl;
	if (membrane.spring_model==1)
	{Report<<"Membrane Spring Model:"<< setw(20)<<"Standard Model"<<endl;}
	if (membrane.spring_model==2)
	{Report<<"Membrane Spring Model:"<< setw(20)<<"Houkian"<<endl;}
	
}

void generatingReport (char* buffer, Membrane membrane, Membrane particle )
{
	string Report_file_name;
	Report_file_name= "Report_";
	Report_file_name+=buffer;
	
	ofstream Report;
	Report.open(Report_file_name.c_str());
	Report<< std:: fixed;
	Report<<"***General Constants***"<<endl;
	Report<<"Number of MD steps"<< setw(20)<<MD_num_of_steps<<endl;
	Report<<"MD time step"<< setw(20)<<MD_Time_Step<<endl;
	Report<<"KT"<< setw(20)<<KT<<endl;
	Report<<"MD time step"<< setw(20)<<MD_Time_Step<<endl;
	Report<<"***Membrane Properties***"<<endl;
	Report<<"Membrane Node Mass"<< setw(20)<<membrane.Node_Mass<<endl;
	Report<<"Membrane Radius"<< setw(20)<<membrane.Radius<<endl;
	Report<<"number of Membrane's Nodes "<< setw(20)<<membrane.Num_of_Nodes<<endl;
	Report<<"number of Membrane's Triangles "<< setw(20)<<membrane.Num_of_Triangles<<endl;
	Report<<"number of Membrane's Nodes"<< setw(20)<<membrane.Num_of_Nodes<<endl;
	Report<<"Membrane Spring coefficient"<< setw(20)<<membrane.Spring_coefficient<<endl;
	Report<<"Membrane Bending coefficient"<< setw(20)<<membrane.Bending_coefficient<<endl;
	Report<<"Membrane Damping coefficient"<< setw(20)<<membrane.Damping_coefficient<<endl;
	if (membrane.spring_model==1)
	{Report<<"Membrane Spring Model:"<< setw(20)<<"Standard Model"<<endl;}
	if (membrane.spring_model==2)
	{Report<<"Membrane Spring Model:"<< setw(20)<<"Houkian"<<endl;}
	Report<<"***Particle Properties***"<<endl;
	Report<<"Particle Node Mass"<< setw(20)<<particle.Node_Mass<<endl;
	Report<<"Particle Radius"<< setw(20)<<particle.Radius<<endl;
	Report<<"number of Particle's Nodes "<< setw(20)<<particle.Num_of_Nodes<<endl;
	Report<<"number of Particle's Triangles "<< setw(20)<<particle.Num_of_Triangles<<endl;
	Report<<"number of Particle's Nodes"<< setw(20)<<particle.Num_of_Nodes<<endl;
	Report<<"Particle Spring coefficient"<< setw(20)<<particle.Spring_coefficient<<endl;
	Report<<"Particle Bending coefficient"<< setw(20)<<particle.Bending_coefficient<<endl;
	Report<<"Particle Damping coefficient"<< setw(20)<<particle.Damping_coefficient<<endl;
	if (particle.spring_model==1)
	{Report<<"Particle Spring Model:"<< setw(20)<<"Standard Model"<<endl;}
	if (particle.spring_model==2)
	{Report<<"Membrane Spring Model:"<< setw(20)<<"Houkian"<<endl;}
}


void generatingReport (char* buffer, Membrane membrane, ECM ecm )
{
	string Report_file_name;
	Report_file_name= "Report_";
	Report_file_name+=buffer;
	
	ofstream Report;
	Report.open(Report_file_name.c_str());
	Report<< std:: fixed;
	Report<<"***General Constants***"<<endl;
	Report<<"Number of MD steps"<< setw(20)<<MD_num_of_steps<<endl;
	Report<<"MD time step"<< setw(20)<<MD_Time_Step<<endl;
	Report<<"KT"<< setw(20)<<KT<<endl;
	Report<<"MD time step"<< setw(20)<<MD_Time_Step<<endl;
	Report<<"***Membrane Properties***"<<endl;
	Report<<"Membrane Node Mass"<< setw(20)<<membrane.Node_Mass<<endl;
	Report<<"Membrane Radius"<< setw(20)<<membrane.Radius<<endl;
	Report<<"number of Membrane's Nodes "<< setw(20)<<membrane.Num_of_Nodes<<endl;
	Report<<"number of Membrane's Triangles "<< setw(20)<<membrane.Num_of_Triangles<<endl;
	Report<<"number of Membrane's Nodes"<< setw(20)<<membrane.Num_of_Nodes<<endl;
	Report<<"Membrane Spring coefficient"<< setw(20)<<membrane.Spring_coefficient<<endl;
	Report<<"Membrane Bending coefficient"<< setw(20)<<membrane.Bending_coefficient<<endl;
	Report<<"Membrane Damping coefficient"<< setw(20)<<membrane.Damping_coefficient<<endl;
	if (membrane.spring_model==1)
	{Report<<"Membrane Spring Model:"<< setw(20)<<"Standard Model"<<endl;}
	if (membrane.spring_model==2)
	{Report<<"Membrane Spring Model:"<< setw(20)<<"Houkian"<<endl;}
	Report<<"***ECM Properties***"<<endl;
	// 
	//Report<<"number of ECM's Nodes "<< setw(20)<<ecm.Num_of_Nodes<<endl;
	//Report<<"number of ECM's Triangles "<< setw(20)<<ecm.Num_of_Triangles<<endl;
}