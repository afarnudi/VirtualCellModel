//
//  main.cpp
//  Test
//
//  Created by Sajjad Sadeghizadeh on 6/24/19.
//  Copyright © 2019 Sajjad Sadeghizadeh. All rights reserved.
//

#include <iostream>
//#include <fstream>
//#include <stdlib.h>
//#include <string>
//#include <iostream>
//#include <fstream>
#include <vector>
//#include <math.h>
//#include <map>
//#include <iomanip>
//#include <iterator>


using namespace std ;



int main(){
    int type=0;
    cout<<"Welcome to VCM config file generator. Please Choose a method. "<<'\n'<<"1- simple" << '\n' << "2- complete"<<'\n'<< "3- example"<<'\n';
    cin>> type;
    
    FILE* pFile;
    pFile = fopen ("general_config_file.txt","w");
    
    int num_mem , num_act , num_ecm , num_chrom ;
    int total_steps , traj_steps;
    int load_checkpoint = 0;
    int integrator = 0 ;
    double step_size;
    double friction=0;
    double temperature = 300;
    //string mesh_path , old_mesh_path;
    vector<string> mesh_path;
    double d_input;
    int i_input;
    string s_input;
    
    if (type==1)
    {
        cout<<"enter MD total number of steps"<<'\n';
        cin>>total_steps;
        cout<<"enter step size in fs"<<'\n';
        cin>>step_size;
        cout<<"enter saving rate in number of steps"<<'\n';
        cin>>traj_steps;
        cout<<"choose integrator type \n";
        cout<<"verlet --> 0 \n" << "brownian --> 1 \n" << "langevin --> 2 \n";
        cin>>integrator;
        
        if (integrator>0)
        {
            cout<<"enter friction in inverse ps"<<'\n';
            cin>>friction;
            cout<<"enter temperature in Kelvin"<<'\n';
            cin>>temperature;
            
        }
        
        
        fprintf(pFile,"//general_config_file \n");
        fprintf(pFile,"MD_num_of_steps  %d \n" , total_steps);
        fprintf(pFile,"MD_traj_save_step  %d \n" , traj_steps);
        fprintf(pFile, "Step_Size_In_Fs %6.3f \n" , step_size);
        fprintf(pFile,"Load_from_checkpoint  %d \n" , load_checkpoint);
        fprintf(pFile,"Integrator_type  %d \n" , integrator);
        if(integrator>0)
        {
            fprintf(pFile, "frictionInPs %6.2f \n" , friction);
            fprintf(pFile, "temperature %6.2f \n" , temperature);
        }
        
        
        vector<double> mass, radius, stiffness, bending, shift_x, shift_y , shift_z, vel_x, vel_y, vel_z, sigma_LJ, epsilon_LJ, f_x, f_y, f_z, kelvin_damp, contractile_k1, contractile_k2, contractile_force, contractile_hill_co, contractile_rmin, contractile_rmax, actin_r0, abp_spring_coefficient, abp_k1, abp_k2, abp_force, abp_hill_co, abp_rmin, abp_rmax, abp_r0 , stiffness_gradient_x, stiffness_gradient_y, stiffness_gradient_z, receptor_density, receptor_gradient_x, receptor_gradient_y, receptor_gradient_z, receptor_center_x, receptor_center_y, receptor_center_z ;
        
        vector<int> spring_model, ext_force,class_type, contractile_model, abp_spring_model, abp_model, receptor_type;
        //int spring_model = 2;
        // vector<double> stiffness;
        //vector<double> bending;
        //double stiffness = 1,bending = 1;
        //vector<double> shift_x;
        //vector<double> shift_y;
        //vector<double> shift_z;
        //double shift_x = 0,shift_y = 0,shift_z =0;
        //double x_speed = 0,y_speed = 0,z_speed = 0;
        //double sigma_LJ = 1, epsilon_LJ = 0;
        //double f_x = 0,f_y = 0,f_z =0;
        // int ext_force=0;
        //int class_type=1;
        //double kelvin_damp=0;
        //double contractile_k1=0, contractile_k2=0;
        //double contractile_force=0,contractile_hill_co=0;
        //int contractile_model=1;
        // double contractile_rmin=0 ,contractile_rmax=10;
        //double actin_r0 = 1;
        // int abp_spring_model=2;
        //        double abp_spring_coefficient=0;
        //        double abp_k1=0, abp_k2=0;
        //        double abp_force=0,abp_hill_co=0;
        //        int abp_model=1;
        //        double abp_rmin=0 ,abp_rmax=10;
        //        double abp_r0 = 1;
        
        
        
        
        
        
        cout<<"How many membranes do you want to initialise? \n";
        cin>>num_mem;
        
        
        for(int i=0; i<num_mem ; i++)
        {
            cout<<"for membrane number "<<i << '\n';
            
            if (i>0)
            {
                cout<<"enter the path (relative to binary) to the mesh file. if you want to initialize a same membrane, type repeat. \n ";
                cin>>s_input;
                mesh_path.push_back(s_input);
                // cin>>mesh_path;
            }
            else
            {
                cout<<"enter the path (relative to binary) to the mesh file.  \n ";
                cin>>s_input;
                mesh_path.push_back(s_input);
                //cin>>mesh_path;
            }
            
            
            if (mesh_path[i].compare("repeat") != 0)
            {
                //old_mesh_path = mesh_path;
                cout<<"enter node mass \n";
                cin>>d_input;
                mass.push_back(d_input);
                cout<<"enter node radius \n";
                cin>>d_input;
                radius.push_back(d_input);
                cout<<"enter spring model \n";
                cin>>i_input;
                spring_model.push_back(i_input);
                cout<<"enter spring stiffness coefficient \n";
                cin>>d_input;
                stiffness.push_back(d_input);
                
                if(spring_model[i]==4)
                {
                    cout<<"enter Kelvin damping coefficient";
                    cin>>d_input;
                    kelvin_damp.push_back(d_input);
                }
                else
                {
                    kelvin_damp.push_back(0);
                }
                
                cout<<"enter bending coefficient \n";
                cin>>d_input;
                bending.push_back(d_input);
                cout<<"enter shift in X in nm \n";
                cin>>d_input;
                shift_x.push_back(d_input);
                cout<<"enter shift in Y in nm \n";
                cin>>d_input;
                shift_y.push_back(d_input);
                cout<<"enter shift in Z in nm \n";
                cin>>d_input;
                shift_z.push_back(d_input);
                cout<<"enter initial velocity in X nm/ps \n";
                cin>>d_input;
                vel_x.push_back(d_input);
                cout<<"enter initial velocity in Y nm/ps \n";
                cin>>d_input;
                vel_y.push_back(d_input);
                cout<<"enter initial velocity in Z nm/ps \n";
                cin>>d_input;
                vel_z.push_back(d_input);
                cout<<"enter sigma for LJ interaction \n";
                cin>>d_input;
                sigma_LJ.push_back(d_input);
                cout<<"enter epsilon for LJ interaction \n";
                cin>>d_input;
                epsilon_LJ.push_back(d_input);
                cout<<"do you want an external force on this object? \n"<<"No --> 0 \n"<<"Yes --> 1 \n";
                cin>>i_input;
                ext_force.push_back(i_input);
                if (ext_force[i]>0)
                {
                    cout<<"enter external force on this object in X \n";
                    cin>>d_input;
                    f_x.push_back(d_input);
                    cout<<"enter external force on this object in X \n";
                    cin>>d_input;
                    f_y.push_back(d_input);
                    cout<<"enter external force on this object in X \n";
                    cin>>d_input;
                    f_z.push_back(d_input);
                }
                else
                {
                    f_x.push_back(0);
                    f_y.push_back(0);
                    f_z.push_back(0);
                }
            }
            
            else
            {
                cout<<"enter the membrane number you want to copy \n";
                int num_copy;
                cin>>num_copy;
                mesh_path[i] = mesh_path[num_copy];
                
                mass.push_back(mass[num_copy]);
                radius.push_back(radius[num_copy]);
                spring_model.push_back(spring_model[num_copy]);
                stiffness.push_back(stiffness[num_copy]);
                
                kelvin_damp.push_back(kelvin_damp[num_copy]);
                
                
                
                bending.push_back(bending[num_copy]);
                
                
                vel_x.push_back(vel_x[num_copy]);
                
                vel_y.push_back(vel_y[num_copy]);
                
                vel_z.push_back(vel_z[num_copy]);
                
                sigma_LJ.push_back(sigma_LJ[num_copy]);
                
                epsilon_LJ.push_back(epsilon_LJ[num_copy]);
                
                ext_force.push_back(ext_force[num_copy]);
                
                
                f_x.push_back(f_x[num_copy]);
                
                f_y.push_back(f_y[num_copy]);
                
                f_z.push_back(f_z[num_copy]);
                
                
                cout<<"enter shift in X in nm \n";
                cin>>d_input;
                shift_x.push_back(d_input);
                cout<<"enter shift in Y in nm \n";
                cin>>d_input;
                shift_y.push_back(d_input);
                cout<<"enter shift in Z in nm \n";
                cin>>d_input;
                shift_z.push_back(d_input);
                
            }
            
            
            //string mem_config_name = "mem" + to_string(i) + "_config_file.txt" ;
            
            
            //            FILE* pFile;
            //            pFile = fopen (mem_config_name.c_str(),"w");
            
            fprintf(pFile,"//membrane_config_file %d ",i);
            //fprintf(pFile, "%s", to_string(i).c_str());
            fprintf(pFile, "\n");
            fprintf(pFile,"Mesh_file_name  %d  ",1);
            fprintf(pFile,"%s", mesh_path[i].c_str());
            fprintf(pFile, "\n");
            fprintf(pFile, "Node_mass %6.2f \n" , mass[i]);
            fprintf(pFile, "Node_radius %6.2f \n" , radius[i]);
            fprintf(pFile,"spring_model  %d  ",spring_model[i]);
            fprintf(pFile, "Spring_coefficient %6.2f \n" , stiffness[i]);
            
            if(spring_model[i]==4)
            {
                fprintf(pFile, "Kelvin_Damping_Coefficient %6.2f \n" , kelvin_damp[i]);
            }
            
            fprintf(pFile, "Bending_coefficient %6.2f \n" , bending[i]);
            fprintf(pFile, "Shift_in_X_direction %6.2f \n" , shift_x[i]);
            fprintf(pFile, "Shift_in_Y_direction %6.2f \n" , shift_y[i]);
            fprintf(pFile, "Shift_in_Z_direction %6.2f \n" , shift_z[i]);
            fprintf(pFile, "x_speed %6.2f \n" , vel_x[i]);
            fprintf(pFile, "y_speed %6.2f \n" , vel_y[i]);
            fprintf(pFile, "z_speed %6.2f \n" , vel_z[i]);
            fprintf(pFile, "sigma_LJ_12_6 %6.2f \n" , sigma_LJ[i]);
            fprintf(pFile, "epsilon_LJ_12_6 %6.2f \n" , epsilon_LJ[i]);
            fprintf(pFile,"ext_force  %d  \n",ext_force[i]);
            fprintf(pFile, "x_force_constant %6.2f \n" , f_x[i]);
            fprintf(pFile, "y_force_constant %6.2f \n" , f_y[i]);
            fprintf(pFile, "z_force_constant %6.2f \n" , f_z[i]);
            
        }
        
        
        
        
        mass.clear();
        radius.clear();
        stiffness.clear();
        bending.clear();
        shift_x.clear();
        shift_y.clear();
        shift_z.clear();
        vel_x.clear();
        vel_y.clear();
        vel_z.clear();
        sigma_LJ.clear();
        epsilon_LJ.clear();
        f_x.clear();
        f_y.clear();
        f_z.clear();
        kelvin_damp.clear();
        spring_model.clear();
        ext_force.clear();
        class_type.clear();
        mesh_path.clear();
        
        
        
        cout<<"How many actins do you want to initialise? \n";
        cin>>num_act;
        
        for(int i=0; i<num_act ; ++i)
        {
            cout<<"for act number "<<i << '\n';
            
            if (i>0)
            {
                cout<<"enter the path (relative to binary) to the mesh file. if you want to initialize a same actin, type repeat. \n ";
                
                cin>>s_input;
                mesh_path.push_back(s_input);
            }
            else
            {
                cout<<"enter the path (relative to binary) to the mesh file.  \n ";
                
                cin>>s_input;
                mesh_path.push_back(s_input);
            }
            
            
            if (mesh_path[i].compare("repeat") != 0)
            {
                cout<<"choose your actin type \n"<<"gmsh mesh --> 1 \n"<<"filaments --> 2 \n";
                cin>>i_input;
                class_type.push_back(i_input);
                cout<<"enter node mass \n";
                cin>>d_input;
                mass.push_back(d_input);
                cout<<"enter node radius \n";
                cin>>d_input;
                radius.push_back(d_input);
                
                cout<<"enter shift in X in nm \n";
                cin>>d_input;
                shift_x.push_back(d_input);
                cout<<"enter shift in Y in nm \n";
                cin>>d_input;
                shift_y.push_back(d_input);
                cout<<"enter shift in Z in nm \n";
                cin>>d_input;
                shift_z.push_back(d_input);
                
                
                cout<<"enter sigma for LJ interaction \n";
                cin>>d_input;
                sigma_LJ.push_back(d_input);
                cout<<"enter epsilon for LJ interaction \n";
                cin>>d_input;
                epsilon_LJ.push_back(d_input);
                cout<<"do you want an external force on this object? \n"<<"No --> 0 \n"<<"Yes --> 1 \n";
                cin>>i_input;
                ext_force.push_back(i_input);
                if (ext_force[i]>0)
                {
                    cout<<"enter external force on this object in X \n";
                    cin>>d_input;
                    f_x.push_back(d_input);
                    cout<<"enter external force on this object in X \n";
                    cin>>d_input;
                    f_y.push_back(d_input);
                    cout<<"enter external force on this object in X \n";
                    cin>>d_input;
                    f_z.push_back(d_input);
                }
                else
                {
                    f_x.push_back(0);
                    f_y.push_back(0);
                    f_z.push_back(0);
                }
                
                if(class_type[i]==1)
                {
                    cout<<"enter spring model \n";
                    cin>>i_input;
                    spring_model.push_back(i_input);
                    cout<<"enter spring stiffness coefficient \n";
                    cin>>d_input;
                    stiffness.push_back(d_input);
                    
                    if(spring_model[i]==4)
                    {
                        cout<<"enter Kelvin damping coefficient";
                        cin>>d_input;
                        kelvin_damp.push_back(d_input);
                    }
                    else
                    {
                        kelvin_damp.push_back(0);
                    }
                    
                    
                    contractile_model.push_back(1);
                    contractile_k1.push_back(0);
                    contractile_k2.push_back(0);
                    contractile_force.push_back(0);
                    contractile_hill_co.push_back(0);
                    contractile_rmin.push_back(0);
                    contractile_rmax.push_back(10);
                    actin_r0.push_back(1);
                    
                    abp_spring_model.push_back(2);
                    abp_model.push_back(1);
                    abp_k1.push_back(0);
                    abp_k2.push_back(0);
                    abp_force.push_back(0);
                    abp_hill_co.push_back(0);
                    abp_rmin.push_back(0);
                    abp_rmax.push_back(10);
                    abp_r0.push_back(1);
                    
                    
                }
                
                
                if(class_type[i]==2)
                {
                    cout<<"setup actin filaments.each filament consists of three springs and a contractile element with constant force parallel";
                    
                    cout<<"enter actin filament main spring model \n";
                    cin>>i_input;
                    spring_model.push_back(i_input);
                    cout<<"enter actin filament spring stiffness coefficient \n";
                    cin>>d_input;
                    stiffness.push_back(d_input);
                    
                    cout<<"enter actin filament second spring stiffness coefficient. this spring works when filament length is less than nominal length. \n ";
                    cin>>d_input;
                    contractile_k1.push_back(d_input);
                    cout<<"enter actin filament third spring stiffness coefficient. this spring works when filament length is more than nominal length. \n ";
                    cin>>d_input;
                    contractile_k2.push_back(d_input);
                    cout<<"enter contractile force constant. \n";
                    cin>>d_input;
                    contractile_force.push_back(d_input);
                    
                    cout<<"enter contractile force type. \n" << "constant --> 1 \n"<<"hill model --> 2 \n" ;
                    cin>>i_input;
                    contractile_model.push_back(i_input);
                    if(contractile_model[i]==2)
                    {
                        cout<<"hill coefficient \n";
                        cin>>d_input;
                        contractile_hill_co.push_back(d_input);
                    }
                    else
                    {
                        contractile_hill_co.push_back(0);
                    }
                    cout<<"enter rmin factor \n";
                    cin>>d_input;
                    contractile_rmin.push_back(d_input);
                    cout<<"enter rmax factor \n";
                    cin>>d_input;
                    contractile_rmax.push_back(d_input);
                    cout<<"enter r0 factor \n";
                    cin>>d_input;
                    actin_r0.push_back(d_input);
                    
                    
                    
                    cout<<"setup abp filaments.each filament consists of three springs and a contractile element with constant force parallel";
                    
                    cout<<"enter abp filament main spring model \n";
                    cin>>i_input;
                    abp_spring_model.push_back(i_input);
                    cout<<"enter abp filament spring stiffness coefficient \n";
                    cin>>d_input;
                    abp_spring_coefficient.push_back(d_input);
                    
                    cout<<"enter abp filament second spring stiffness coefficient. this spring works when filament length is less than nominal length. \n ";
                    cin>>d_input;
                    abp_k1.push_back(d_input);
                    cout<<"enter abp filament third spring stiffness coefficient. this spring works when filament length is more than nominal length. \n ";
                    cin>>d_input;
                    abp_k2.push_back(d_input);
                    cout<<"enter contractile force constant. \n";
                    cin>>d_input;
                    abp_force.push_back(d_input);
                    
                    cout<<"enter abp force type. \n" << "constant --> 1 \n"<<"hill model --> 2 \n" ;
                    cin>>i_input;
                    abp_model.push_back(i_input);
                    if(abp_model[i]==2)
                    {
                        cout<<"hill coefficient \n";
                        cin>>d_input;
                        abp_hill_co.push_back(d_input);
                    }
                    else
                    {
                        abp_hill_co.push_back(0);
                    }
                    cout<<"enter rmin factor \n";
                    cin>>d_input;
                    abp_rmin.push_back(d_input);
                    cout<<"enter rmax factor \n";
                    cin>>d_input;
                    abp_rmax.push_back(d_input);
                    cout<<"enter r0 factor \n";
                    cin>>d_input;
                    abp_r0.push_back(d_input);
                    
                    
                    if((spring_model[i]==4) || (abp_spring_model[i]==4))
                    {
                        cout<<"enter Kelvin damping coefficient";
                        cin>>d_input;
                        kelvin_damp.push_back(d_input);
                    }
                    else
                    {
                        kelvin_damp.push_back(0);
                    }
                    
                }
                
                
            }
            
            else
            {
                cout<<"enter the actin number you want to copy \n";
                int num_copy;
                cin>>num_copy;
                mesh_path[i] = mesh_path[num_copy];
                
                class_type.push_back(class_type[num_copy]);
                mass.push_back(mass[num_copy]);
                radius.push_back(radius[num_copy]);
                spring_model.push_back(spring_model[num_copy]);
                stiffness.push_back(stiffness[num_copy]);
                
                kelvin_damp.push_back(kelvin_damp[num_copy]);
                
                //                vel_x.push_back(vel_x[num_copy]);
                //
                //                vel_y.push_back(vel_y[num_copy]);
                //
                //                vel_z.push_back(vel_z[num_copy]);
                
                sigma_LJ.push_back(sigma_LJ[num_copy]);
                
                epsilon_LJ.push_back(epsilon_LJ[num_copy]);
                
                ext_force.push_back(ext_force[num_copy]);
                
                
                f_x.push_back(f_x[num_copy]);
                
                f_y.push_back(f_y[num_copy]);
                
                f_z.push_back(f_z[num_copy]);
                
                
                contractile_model.push_back(contractile_model[num_copy]);
                contractile_k1.push_back(contractile_k1[num_copy]);
                contractile_k2.push_back(contractile_k2[num_copy]);
                contractile_force.push_back(contractile_force[num_copy]);
                contractile_hill_co.push_back(contractile_hill_co[num_copy]);
                contractile_rmin.push_back(contractile_rmin[num_copy]);
                contractile_rmax.push_back(contractile_rmax[num_copy]);
                actin_r0.push_back(actin_r0[num_copy]);
                
                abp_spring_model.push_back(abp_spring_model[num_copy]);
                abp_model.push_back(abp_model[num_copy]);
                abp_k1.push_back(abp_k1[num_copy]);
                abp_k2.push_back(abp_k2[num_copy]);
                abp_force.push_back(abp_force[num_copy]);
                abp_hill_co.push_back(abp_hill_co[num_copy]);
                abp_rmin.push_back(abp_rmin[num_copy]);
                abp_rmax.push_back(abp_rmax[num_copy]);
                abp_r0.push_back(abp_r0[num_copy]);
                
                cout<<"enter shift in X in nm \n";
                cin>>d_input;
                shift_x.push_back(d_input);
                cout<<"enter shift in Y in nm \n";
                cin>>d_input;
                shift_y.push_back(d_input);
                cout<<"enter shift in Z in nm \n";
                cin>>d_input;
                shift_z.push_back(d_input);
            }
            
            
            
            
            fprintf(pFile,"//actin_config_file %d ",i);
            //fprintf(pFile, "%s", to_string(i).c_str());
            fprintf(pFile, "\n");
            fprintf(pFile,"Mesh_file_name  %d  ",class_type[i]);
            fprintf(pFile,"%s", mesh_path[i].c_str());
            fprintf(pFile, "\n");
            fprintf(pFile, "Node_mass %6.2f \n" , mass[i]);
            fprintf(pFile, "Node_radius %6.2f \n" , radius[i]);
            
            fprintf(pFile,"spring_model  %d  ",spring_model[i]);
            fprintf(pFile, "Spring_coefficient %6.2f \n" , stiffness[i]);
            
            if(class_type[i]==2)
            {
                fprintf(pFile,"Contractile_model  %d  ",contractile_model[i]);
                fprintf(pFile, "Contractile_hill_co %6.2f \n" , contractile_hill_co[i]);
                fprintf(pFile, "Contractile_force %6.2f \n" , contractile_force[i]);
                fprintf(pFile, "Contractile_k1 %6.2f \n" , contractile_k1[i]);
                fprintf(pFile, "Contractile_k2 %6.2f \n" , contractile_k2[i]);
                fprintf(pFile, "Contractile_rmin_factor %6.2f \n" , contractile_rmin[i]);
                fprintf(pFile, "Contractile_rmax_factor %6.2f \n" , contractile_rmax[i]);
                fprintf(pFile, "actin_r0factor %6.2f \n" , actin_r0[i]);
                
                fprintf(pFile,"abp_spring_model  %d  ",abp_spring_model[i]);
                fprintf(pFile, "abp_Spring_coefficient %6.2f \n" , abp_spring_coefficient[i]);
                fprintf(pFile,"abp_model  %d  ",abp_model[i]);
                fprintf(pFile, "abp_hill_co %6.2f \n" , abp_hill_co[i]);
                fprintf(pFile, "abp_force %6.2f \n" , abp_force[i]);
                fprintf(pFile, "abp_k1 %6.2f \n" , abp_k1[i]);
                fprintf(pFile, "abp_k2 %6.2f \n" , abp_k2[i]);
                fprintf(pFile, "abp_rmin_factor %6.2f \n" , abp_rmin[i]);
                fprintf(pFile, "abp_rmax_factor %6.2f \n" , abp_rmax[i]);
                fprintf(pFile, "abp_r0factor %6.2f \n" , abp_r0[i]);
            }
            
            
            
            if((spring_model[i]==4) || (abp_spring_model[i]==4))
            {
                fprintf(pFile, "Kelvin_Damping_Coefficient %6.2f \n" , kelvin_damp[i]);
            }
            
            
            
            //                    fprintf(pFile, "Bending_coefficient %6.2f \n" , bending);
            fprintf(pFile, "Shift_in_X_direction %6.2f \n" , shift_x[i]);
            fprintf(pFile, "Shift_in_Y_direction %6.2f \n" , shift_y[i]);
            fprintf(pFile, "Shift_in_Z_direction %6.2f \n" , shift_z[i]);
            //                    fprintf(pFile, "x_speed %6.2f \n" , x_speed);
            //                    fprintf(pFile, "y_speed %6.2f \n" , y_speed);
            //                    fprintf(pFile, "z_speed %6.2f \n" , z_speed);
            fprintf(pFile, "sigma_LJ_12_6 %6.2f \n" , sigma_LJ[i]);
            fprintf(pFile, "epsilon_LJ_12_6 %6.2f \n" , epsilon_LJ[i]);
            fprintf(pFile,"ext_force  %d  \n",ext_force[i]);
            fprintf(pFile, "x_force_constant %6.2f \n" , f_x[i]);
            fprintf(pFile, "y_force_constant %6.2f \n" , f_y[i]);
            fprintf(pFile, "z_force_constant %6.2f \n" , f_z[i]);
            
        }
        
        
        
        mass.clear();
        radius.clear();
        stiffness.clear();
        bending.clear();
        shift_x.clear();
        shift_y.clear();
        shift_z.clear();
        vel_x.clear();
        vel_y.clear();
        vel_z.clear();
        sigma_LJ.clear();
        epsilon_LJ.clear();
        f_x.clear();
        f_y.clear();
        f_z.clear();
        kelvin_damp.clear();
        spring_model.clear();
        ext_force.clear();
        class_type.clear();
        mesh_path.clear();
        
        contractile_model.clear();
        contractile_k1.clear();
        contractile_k2.clear();
        contractile_force.clear();
        contractile_hill_co.clear();
        contractile_rmin.clear();
        contractile_rmax.clear();
        actin_r0.clear();
        
        abp_spring_model.clear();
        abp_model.clear();
        abp_k1.clear();
        abp_k2.clear();
        abp_force.clear();
        abp_hill_co.clear();
        abp_rmin.clear();
        abp_rmax.clear();
        abp_r0.clear();
        
        
        
        
        
        
        
        
        
        cout<<"How many ECMs do you want to initialise? \n";
        cin>>num_ecm;
        
        for(int i=0; i<num_ecm ; i++)
        {
            cout<<"for ECM number "<<i << '\n';
            
            if (i>0)
            {
                cout<<"enter the path (relative to binary) to the mesh file. if you want to initialize a same ECM, type repeat. \n ";
                cin>>s_input;
                mesh_path.push_back(s_input);
                // cin>>mesh_path;
            }
            else
            {
                cout<<"enter the path (relative to binary) to the mesh file.  \n ";
                cin>>s_input;
                mesh_path.push_back(s_input);
                //cin>>mesh_path;
            }
            
            
            if (mesh_path[i].compare("repeat") != 0)
            {
                //old_mesh_path = mesh_path;
                cout<<"enter ECM dimension \n";
                cin>>i_input;
                class_type.push_back(i_input);
                cout<<"enter node mass \n";
                cin>>d_input;
                mass.push_back(d_input);
                cout<<"enter node radius \n";
                cin>>d_input;
                radius.push_back(d_input);
                cout<<"enter spring model \n";
                cin>>i_input;
                spring_model.push_back(i_input);
                cout<<"enter spring stiffness coefficient \n";
                cin>>d_input;
                stiffness.push_back(d_input);
                
                if(spring_model[i]==4)
                {
                    cout<<"enter Kelvin damping coefficient";
                    cin>>d_input;
                    kelvin_damp.push_back(d_input);
                }
                else
                {
                    kelvin_damp.push_back(0);
                }
                
                //                cout<<"enter bending coefficient \n";
                //                cin>>d_input;
                //                bending.push_back(d_input);
                cout<<"enter shift in X in nm \n";
                cin>>d_input;
                shift_x.push_back(d_input);
                cout<<"enter shift in Y in nm \n";
                cin>>d_input;
                shift_y.push_back(d_input);
                cout<<"enter shift in Z in nm \n";
                cin>>d_input;
                shift_z.push_back(d_input);
                cout<<"enter initial velocity in X nm/ps \n";
                cin>>d_input;
                vel_x.push_back(d_input);
                cout<<"enter initial velocity in Y nm/ps \n";
                cin>>d_input;
                vel_y.push_back(d_input);
                cout<<"enter initial velocity in Z nm/ps \n";
                cin>>d_input;
                vel_z.push_back(d_input);
                cout<<"enter sigma for LJ interaction \n";
                cin>>d_input;
                sigma_LJ.push_back(d_input);
                cout<<"enter epsilon for LJ interaction \n";
                cin>>d_input;
                epsilon_LJ.push_back(d_input);
                cout<<"do you want an external force on this object? \n"<<"No --> 0 \n"<<"Yes --> 1 \n";
                cin>>i_input;
                ext_force.push_back(i_input);
                if (ext_force[i]>0)
                {
                    cout<<"enter external force on this object in X \n";
                    cin>>d_input;
                    f_x.push_back(d_input);
                    cout<<"enter external force on this object in X \n";
                    cin>>d_input;
                    f_y.push_back(d_input);
                    cout<<"enter external force on this object in X \n";
                    cin>>d_input;
                    f_z.push_back(d_input);
                }
                else
                {
                    f_x.push_back(0);
                    f_y.push_back(0);
                    f_z.push_back(0);
                }
                cout<<"enter stiffness gradient in X direction \n";
                cin>>d_input;
                stiffness_gradient_x.push_back(d_input);
                cout<<"enter stiffness gradient in Y direction \n";
                cin>>d_input;
                stiffness_gradient_y.push_back(d_input);
                cout<<"enter stiffness gradient in Z direction \n";
                cin>>d_input;
                stiffness_gradient_z.push_back(d_input);
                
                
                cout<<"enter receptor type on ECM \n"<<"linear --> 1 \n" << "tape --> 2 \n";
                cin>>i_input;
                receptor_type.push_back(i_input);
                
                cout<<"enter receptor density at receptor center. you'll set receptor center in next step. \n";
                cin>>d_input;
                receptor_density.push_back(d_input);
                
                cout<<"enter receptor center in X direction \n";
                cin>>d_input;
                receptor_center_x.push_back(d_input);
                cout<<"enter receptor center in Y direction \n";
                cin>>d_input;
                receptor_center_y.push_back(d_input);
                cout<<"enter receptor center in Z direction \n";
                cin>>d_input;
                receptor_center_z.push_back(d_input);
                if(receptor_type[i]==1)
                {
                    cout<<"enter receptor density gradient in X direction \n";
                    cin>>d_input;
                    receptor_gradient_x.push_back(d_input);
                    cout<<"enter receptor density gradient in Y direction \n";
                    cin>>d_input;
                    receptor_gradient_y.push_back(d_input);
                    cout<<"enter receptor density gradient in Z direction \n";
                    cin>>d_input;
                    receptor_gradient_z.push_back(d_input);
                }
                else
                {
                    receptor_gradient_x.push_back(0);
                    receptor_gradient_y.push_back(0);
                    receptor_gradient_z.push_back(0);
                }
                
                
                
            }
            
            else
            {
                cout<<"enter the ECM number you want to copy \n";
                int num_copy;
                cin>>num_copy;
                mesh_path[i] = mesh_path[num_copy];
                
                class_type.push_back(class_type[num_copy]);
                mass.push_back(mass[num_copy]);
                radius.push_back(radius[num_copy]);
                spring_model.push_back(spring_model[num_copy]);
                stiffness.push_back(stiffness[num_copy]);
                
                kelvin_damp.push_back(kelvin_damp[num_copy]);
                
                
                
                // bending.push_back(bending[num_copy]);
                
                
                vel_x.push_back(vel_x[num_copy]);
                
                vel_y.push_back(vel_y[num_copy]);
                
                vel_z.push_back(vel_z[num_copy]);
                
                sigma_LJ.push_back(sigma_LJ[num_copy]);
                
                epsilon_LJ.push_back(epsilon_LJ[num_copy]);
                
                ext_force.push_back(ext_force[num_copy]);
                
                
                f_x.push_back(f_x[num_copy]);
                
                f_y.push_back(f_y[num_copy]);
                
                f_z.push_back(f_z[num_copy]);
                
                stiffness_gradient_x.push_back(stiffness_gradient_x[num_copy]);
                stiffness_gradient_y.push_back(stiffness_gradient_y[num_copy]);
                stiffness_gradient_z.push_back(stiffness_gradient_z[num_copy]);
                
                
                receptor_type.push_back(receptor_type[num_copy]);
                receptor_density.push_back(receptor_density[num_copy]);
                receptor_center_x.push_back(receptor_center_x[num_copy]);
                receptor_center_y.push_back(receptor_center_y[num_copy]);
                receptor_center_z.push_back(receptor_center_z[num_copy]);
                
                receptor_gradient_x.push_back(receptor_gradient_x[num_copy]);
                receptor_gradient_y.push_back(receptor_gradient_y[num_copy]);
                receptor_gradient_z.push_back(receptor_gradient_z[num_copy]);
                
                cout<<"enter shift in X in nm \n";
                cin>>d_input;
                shift_x.push_back(d_input);
                cout<<"enter shift in Y in nm \n";
                cin>>d_input;
                shift_y.push_back(d_input);
                cout<<"enter shift in Z in nm \n";
                cin>>d_input;
                shift_z.push_back(d_input);
                
            }
            
            
            fprintf(pFile,"//ECM_config_file %d ",i);
            //fprintf(pFile, "%s", to_string(i).c_str());
            fprintf(pFile, "\n");
            fprintf(pFile,"Mesh_file_name  %d  ",class_type[i]);
            fprintf(pFile,"%s", mesh_path[i].c_str());
            fprintf(pFile, "\n");
            fprintf(pFile, "Node_mass %6.2f \n" , mass[i]);
            fprintf(pFile, "Node_radius %6.2f \n" , radius[i]);
            fprintf(pFile,"spring_model  %d  ",spring_model[i]);
            fprintf(pFile, "Spring_coefficient %6.2f \n" , stiffness[i]);
            
            if(spring_model[i]==4)
            {
                fprintf(pFile, "Kelvin_Damping_Coefficient %6.2f \n" , kelvin_damp[i]);
            }
            
            //fprintf(pFile, "Bending_coefficient %6.2f \n" , bending[i]);
            fprintf(pFile, "Shift_in_X_direction %6.2f \n" , shift_x[i]);
            fprintf(pFile, "Shift_in_Y_direction %6.2f \n" , shift_y[i]);
            fprintf(pFile, "Shift_in_Z_direction %6.2f \n" , shift_z[i]);
            fprintf(pFile, "x_speed %6.2f \n" , vel_x[i]);
            fprintf(pFile, "y_speed %6.2f \n" , vel_y[i]);
            fprintf(pFile, "z_speed %6.2f \n" , vel_z[i]);
            fprintf(pFile, "sigma_LJ_12_6 %6.2f \n" , sigma_LJ[i]);
            fprintf(pFile, "epsilon_LJ_12_6 %6.2f \n" , epsilon_LJ[i]);
            fprintf(pFile,"ext_force  %d  \n",ext_force[i]);
            fprintf(pFile, "x_force_constant %6.2f \n" , f_x[i]);
            fprintf(pFile, "y_force_constant %6.2f \n" , f_y[i]);
            fprintf(pFile, "z_force_constant %6.2f \n" , f_z[i]);
            
            fprintf(pFile, "stiffness_gradient_x %6.2f \n" , stiffness_gradient_x[i]);
            fprintf(pFile, "stiffness_gradient_y %6.2f \n" , stiffness_gradient_y[i]);
            fprintf(pFile, "stiffness_gradient_z %6.2f \n" , stiffness_gradient_z[i]);
            
            fprintf(pFile,"receptor_type  %d  \n",receptor_type[i]);
            fprintf(pFile, "receptor_density %6.2f \n" , receptor_density[i]);
            fprintf(pFile, "receptor_center_x %6.2f \n" , receptor_center_x[i]);
            fprintf(pFile, "receptor_center_y %6.2f \n" , receptor_center_y[i]);
            fprintf(pFile, "receptor_center_z %6.2f \n" , receptor_center_z[i]);
            fprintf(pFile, "receptor_gradient_x %6.2f \n" , receptor_gradient_x[i]);
            fprintf(pFile, "receptor_gradient_y %6.2f \n" , receptor_gradient_y[i]);
            fprintf(pFile, "receptor_gradient_z %6.2f \n" , receptor_gradient_z[i]);
            
        }
        
        mass.clear();
        radius.clear();
        stiffness.clear();
        bending.clear();
        shift_x.clear();
        shift_y.clear();
        shift_z.clear();
        vel_x.clear();
        vel_y.clear();
        vel_z.clear();
        sigma_LJ.clear();
        epsilon_LJ.clear();
        f_x.clear();
        f_y.clear();
        f_z.clear();
        kelvin_damp.clear();
        spring_model.clear();
        ext_force.clear();
        class_type.clear();
        mesh_path.clear();
        
//        contractile_model.clear();
//        contractile_k1.clear();
//        contractile_k2.clear();
//        contractile_force.clear();
//        contractile_hill_co.clear();
//        contractile_rmin.clear();
//        contractile_rmax.clear();
//        actin_r0.clear();
//
//        abp_spring_model.clear();
//        abp_model.clear();
//        abp_k1.clear();
//        abp_k2.clear();
//        abp_force.clear();
//        abp_hill_co.clear();
//        abp_rmin.clear();
//        abp_rmax.clear();
//        abp_r0.clear();
        
        stiffness_gradient_x.clear();
        stiffness_gradient_y.clear();
        stiffness_gradient_z.clear();
        receptor_type.clear();
        receptor_density.clear();
        receptor_center_x.clear();
        receptor_center_y.clear();
        receptor_center_z.clear();
        receptor_gradient_x.clear();
        receptor_gradient_y.clear();
        receptor_gradient_z.clear();
        
        
        //clear
        
        cout<<"How many chromatins do you want to initialise? \n";
        cin>>num_chrom;
        
        for(int i=0; i<num_chrom ; i++)
        {
            
        }
        
        
        
    }
    
    else if (type==2)
    {
        
    }
    
    else if (type==3)
    {
        
    }
    
    else
    {
        cout<<"incorrect input"<<'\n';
    }
    
    
    
    fclose(pFile);
}
