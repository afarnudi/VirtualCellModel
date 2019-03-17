//
//  write_functions.cpp
//  Mem
//
//  Created by Ali Farnudi on 19/08/2018.
//  Copyright © 2018 Ali Farnudi. All rights reserved.
//

#include "Actin.h"
#include "General_constants.h"

void Actin::write_traj (string traj_name, string label){
    ofstream Trajectory;
    Trajectory.open(traj_name.c_str(), ios::app);
    Trajectory << std:: fixed;
    for(int j=0; j< Num_of_Nodes;j++) // saving trajectory
    {
        Trajectory << label <<setprecision(5)<< setw(20)<<Node_Position[j][0]<< setw(20)<<Node_Position[j][1]<< setw(20)<<Node_Position[j][2]<<endl;
    }
}





void Actin::generate_report()
{
    string Report_file_name;
    Report_file_name= "Results/Reports/Report_Actin_"+to_string(index)+"_";
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
    Report<<"Num_of_Actins"<<setw(20)<<GenConst::Num_of_Actins<<endl;
    Report<<"trajectory_file_name"<<setw(20)<<GenConst::trajectory_file_name<<endl;
    
    
    Report<<"Actin Params:\n---------------\n";
    Report<<"# of Nodes"<< setw(20) <<Num_of_Nodes<<endl;
    Report<<"Node Mass"<< setw(20)<<Node_Mass<<endl;
    Report<<"Node Radius"<< setw(20)<<Node_radius<<endl;
    Report<<"Number of Node Pairs"<< setw(20)<<Num_of_Node_Pairs<<endl;
    Report<<"spring_model"<< setw(20)<<spring_model<<endl;
    if (spring_model==0)
    {
        Report<<"Actin Spring Model:"<< setw(20)<<"Hookian"<<endl;
    }
    if (spring_model==2)
    {
        Report<<"Actin Spring Model:"<< setw(20)<<"Kelvin"<<endl;
    }
    Report<<"Spring coefficient"<< setw(20)<<Spring_coefficient<<endl;
    Report<<"Kelvin Damping Coefficient"<< setw(20)<<Kelvin_Damping_Coefficient<<endl;
    Report<<"Shift_in_X_direction"<< setw(20)<<Shift_in_X_direction<<endl;
    Report<<"Shift_in_Y_direction"<< setw(20)<<Shift_in_Y_direction<<endl;
    Report<<"Shift_in_Z_direction"<< setw(20)<<Shift_in_Z_direction<<endl;
    Report<<"Downward_speed"<< setw(20)<<Downward_speed<<endl;
    
}
