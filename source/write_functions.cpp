//
//  write_functions.cpp
//  Mem
//
//  Created by Ali Farnudi on 19/08/2018.
//  Copyright © 2018 Ali Farnudi. All rights reserved.
//

#include "write_functions.hpp"



//void Results (ECM ecm, string label, char* buffer)
//{
//    //output file names:
//    string energy_file_name;
//    string traj_file_name;
//    ecm.output_file_neme+=buffer;
//    traj_file_name="Results/Results_";
//    traj_file_name+=ecm.output_file_neme;
//    traj_file_name+=".xyz";
//    energy_file_name="Results_Potential_Energy";
//    energy_file_name+=ecm.output_file_neme;
//
//    //trajectory:
//
//    ofstream Trajectory;
//
//    Trajectory.open(traj_file_name.c_str(), ios::app);
//    Trajectory << std:: fixed;
//    Trajectory << ecm.return_num_of_nodes()<<endl;
//    Trajectory << " nodes  "<<endl;
//    for(int j=0; j< ecm.return_num_of_nodes();j++) // saving trajectory
//    {
//        Trajectory << label <<setprecision(5)<< setw(20)<<ecm.Node_Position[j][0]<< setw(20)<<ecm.Node_Position[j][1]<< setw(20)<<ecm.Node_Position[j][2]<<endl;
//    }
//
//    //    //Energy:
//    //    ofstream Potential_Energy;
//    //    Potential_Energy.open(energy_file_name.c_str(), ios::app);
//    //    Potential_Energy<<membrane.Total_Potential_Energy<<endl;
//}

//void Results (Membrane membrane, string label, char* buffer)
//{
//    //output file names:
//
//    string energy_file_name;
//    string traj_file_name;
//
//    membrane.output_file_neme+=buffer;
//    traj_file_name="Results/Results_";
//    traj_file_name+=membrane.output_file_neme;
//    traj_file_name+=".xyz";
//    energy_file_name="Results/Results_Potential_Energy";
//    energy_file_name+=membrane.output_file_neme;
//    //trajectory:
//
//    ofstream Trajectory;
//
//    Trajectory.open(traj_file_name.c_str(), ios::app);
//    Trajectory << std:: fixed;
//    Trajectory << membrane.return_num_of_nodes()<<endl;
//    Trajectory << " nodes  "<<endl;
//    for(int j=0; j< membrane.return_num_of_nodes();j++) // saving trajectory
//    {
//        Trajectory << label <<setprecision(5)<< setw(20)<<membrane.Node_Position[j][0]<< setw(20)<<membrane.Node_Position[j][1]<< setw(20)<<membrane.Node_Position[j][2]<<endl;
//    }
//
//
//
//
//    //    //Energy:
//    //    ofstream Potential_Energy;
//    //    Potential_Energy.open(energy_file_name.c_str(), ios::app);
//    //    Potential_Energy<<membrane.Total_Potential_Energy<<endl;
//}



//
//void checkingForce (Membrane membrane, int MD_Step, char* buffer)
//{
//    string check_force_name="Results/check_force_";
//    check_force_name+=buffer;
//    ofstream check_force;
//    check_force.open(check_force_name.c_str(), ios::app);
//    check_force<< "MD step = "<<MD_Step<<endl;
//    for(int j=0; j< membrane.return_num_of_nodes();j++) // saving Forces
//    {
//        check_force  <<setprecision(5)<< setw(20)<<membrane.Node_Force[j][0]<< setw(20)<<membrane.Node_Force[j][1]<< setw(20)<<membrane.Node_Force[j][2]<<endl;
//    }
//    check_force<< "_________________________________________________________"<<endl;
//}



