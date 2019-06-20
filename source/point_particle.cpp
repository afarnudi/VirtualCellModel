#include "point_particle.h"
#include "General_constants.h"
#include <sstream>

using namespace std;

void point_particle::initialize_point_particle(double x, double y, double z)
{
    Node_Position.resize(3); 
    Node_Velocity.resize(3);
    Neighbour.resize(3);
    Node_Force.resize(3);
    Node_Position[0]=x;
    Node_Position[1]=y;
    Node_Position[2]=z;
    for (int i=0; i<3; i++){
        Node_Velocity[i]=0;
        Node_Force[i]=0;

    }
    cout<<"point particle class is initiated"<<endl;
}




 void point_particle::MD_Evolution_beginning(double MD_Time_Step){

            Node_Position[0] += Node_Velocity[0]*MD_Time_Step + Node_Force[0]*MD_Time_Step*MD_Time_Step/(Node_Mass*2.0);
            Node_Position[1] += Node_Velocity[1]*MD_Time_Step + Node_Force[1]*MD_Time_Step*MD_Time_Step/(Node_Mass*2.0);
            Node_Position[2] += Node_Velocity[2]*MD_Time_Step + Node_Force[2]*MD_Time_Step*MD_Time_Step/(Node_Mass*2.0);


            Node_Velocity[0] += Node_Force[0]*MD_Time_Step/(Node_Mass*2.0);
            Node_Velocity[1] += Node_Force[1]*MD_Time_Step/(Node_Mass*2.0);
            Node_Velocity[2] += Node_Force[2]*MD_Time_Step/(Node_Mass*2.0);


            Node_Force[0]=0.0;
            Node_Force[1]=0.0;
            Node_Force[2]=0.0; 
        }

 

void point_particle::MD_Evolution_end(double MD_Time_Step){
    
        
        Node_Velocity[0] += Node_Force[0]*MD_Time_Step/(Node_Mass*2.0);
        Node_Velocity[1] += Node_Force[1]*MD_Time_Step/(Node_Mass*2.0);
        Node_Velocity[2] += Node_Force[2]*MD_Time_Step/(Node_Mass*2.0);
        
    
    
}

void point_particle::write_traj (string traj_name){
    ofstream Trajectory;
    Trajectory.open(traj_name.c_str(), ios::app);
    Trajectory << std:: fixed;

        Trajectory << "pointparticle" <<setprecision(5)<< setw(20)<<Node_Position[0]<< setw(20)<<Node_Position[1]<< setw(20)<<Node_Position[2]<<endl;
}

void point_particle::generate_report(){
    {
    string Report_file_name;
    Report_file_name= "Results/Reports/Report_pointparticle_"+to_string(index)+"_";
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

    Report<<"Num_of_point particles"<<setw(20)<<GenConst::Num_of_pointparticles<<endl;
    Report<<"trajectory_file_name"<<setw(20)<<GenConst::trajectory_file_name<<endl;
    
    
    Report<<"point particle Params:\n---------------\n";
    Report<<"Node Mass"<< setw(20)<<Node_Mass<<endl;
    
    Report<<"X_position"<< setw(20)<<Node_Position[0]<<endl;
    Report<<"Y_position"<< setw(20)<<Node_Position[1]<<endl;
    Report<<"Z_position"<< setw(20)<<Node_Position[2]<<endl;
    
    Report<<"point particle- point particle interaction Params:\n---------------\n";
    Report<<"Sigma"<<return_P_P_sigma()<<endl;
    Report<<"epsilon"<<return_P_P_epsilon()<<endl;
    Report<<"cut off"<<return_P_P_cut_off()<<endl;

    Report<<"point particle- Membrane interaction Params:\n---------------\n";
    Report<<"Sigma"<<return_P_Membrane_sigma()<<endl;
    Report<<"epsilon"<<return_P_Membrane_epsilon()<<endl;
    Report<<"cut off"<<return_P_Membrane_cut_off()<<endl;

    
}
}

void point_particle::import_config(string config_file_name){
    
    map<string, double>::iterator it;
    
    
    //string resume_file_name, Mesh_file_name="non";
    
    
    ifstream read_config_file(config_file_name.c_str());

    if (read_config_file.is_open()) {
        cout<<"'"<<config_file_name<<"' file opened successfully.\n";
        string line;
        int line_num=0;
        string comment="//";
        
        //        char delimiter=' ';
        while(getline(read_config_file, line)){
            line_num++;
            if(line.empty()){
                continue;
            }
            
            istringstream iss(line);
            vector<string> split(istream_iterator<string>{iss}, istream_iterator<string>());
            
            if (split[0] == comment || (split[0][0]=='/' && split[0][1]=='/')) {
                break;
            }
            
            param_map[split[0]]=stod(split[1]);
            set_map_parameter(split[0], param_map[split[0]]);
            
                //                break;
                //                cout<<split[i]<<"\t";
            
        }//End of while(getline(read_config_file, line))
        cout<<"here"<<endl;
        
    } else {//End of if (read_config_file.is_open())
        cout<<"Couldn't open the '"<<config_file_name<<"' file.\n";
        exit(EXIT_FAILURE);
    }
    
      
    initialize_point_particle(X_0,Y_0, Z_0);  
    
}

void point_particle::set_map_parameter(string param_name, double param_value){

    if (param_name=="Node_Mass") {
        Node_Mass=param_value;
    } else if (param_name=="on_or_off_MD_evolution"){
       on_or_off_MD_evolution=param_value;
    } else if (param_name=="X_0"){
        X_0=param_value;
    } else if (param_name=="Y_0"){
        Y_0=param_value;
    }  else if (param_name=="Z_0"){
       Z_0=param_value;
    } else if (param_name=="P_P_sigma"){
        P_P_sigma=param_value;
    } else if (param_name=="P_P_epsilon"){
        P_P_epsilon=param_value;
    } else if (param_name=="P_P_cut_off"){
        P_P_cut_off=param_value;
    } else if (param_name=="P_Membrane_sigma"){
        P_Membrane_sigma=param_value;
    } else if (param_name=="P_Membrane_epsilon"){
        P_Membrane_epsilon=param_value;
    } else if (param_name=="P_Membrane_cut_off"){
        P_Membrane_cut_off=param_value;
    }
    
}

//point_particle::~point_particle()
//{
//}

