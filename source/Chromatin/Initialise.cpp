//
//  Initialise.cpp
//  Membrae
//
//  Created by Ali Farnudi on 14/10/2018.
//  Copyright © 2018 Ali Farnudi. All rights reserved.
//

#include "Chromatin.h"
#include <algorithm>
//void build_random_chain(void);

using std::cout;
using std::endl;

void Chromatin::initialise(void){
    vector<vector<double> > zero_vec(Num_of_Nodes, vector<double> (Num_of_Nodes));
    Contact_Matrix=zero_vec;
    
    AB_index.resize(Num_of_Nodes);
    generate(AB_index.begin(), AB_index.end(), []() {
        return rand() % 2;
    });
    
    cout<<"\nInitialising the Chromatin Class..."<<endl;
    build_random_chain();
//    Pack(7);
    shift_node_positions();
    cout<<"Chromatin class initiated.\n\n";
    
}

void Chromatin::initialise(double min_radius){
    vector<vector<double> > zero_vec(Num_of_Nodes, vector<double> (Num_of_Nodes));
    Contact_Matrix=zero_vec;
    
    AB_index.resize(Num_of_Nodes);
    generate(AB_index.begin(), AB_index.end(), []() {
        return rand() % 2;
    });
    
    cout<<"\nInitialising the Chromatin Class..."<<endl;
    build_random_chain();
    Pack(min_radius-2);
    shift_node_positions();
    cout<<"Chromatin class initiated.\n";
    
}


void Chromatin::build_random_chain(void){
    
    Node_Position.resize(Num_of_Nodes);
    Node_Velocity.resize(Num_of_Nodes);
    Node_Force.resize(Num_of_Nodes);
    
    
    double velocity_COM[3]={0};
    
    for(int i=0;i<Num_of_Nodes;i++)
    {
        Node_Velocity[i].resize(3,0);
        Node_Force[i].resize(3,0);
        Node_Position[i].resize(3,0);
    }
    const double PI  =3.141592653589793238463;
    int attempt_counter=0;
    
    for (int i=1; i<Num_of_Nodes; i++) {
        double theta=((double)rand()/(double)RAND_MAX)*PI;
        double phi=((double)rand()/(double)RAND_MAX)*2*PI;
        double temp_x=Node_Position[i-1][0]+Node_radius*2*sin(theta)*cos(phi);
        double temp_y=Node_Position[i-1][1]+Node_radius*2*sin(theta)*sin(phi);
        double temp_z=Node_Position[i-1][2]+Node_radius*2*cos(theta);
        
        bool accept_random_step=true;
        for (int k=0; k<i-1; k++) {
            double temp_dist=sqrt( (temp_x-Node_Position[k][0])*(temp_x-Node_Position[k][0]) + (temp_y-Node_Position[k][1])*(temp_y-Node_Position[k][1]) + (temp_z-Node_Position[k][2])*(temp_z-Node_Position[k][2])  );
            if (temp_dist<Node_radius*2.1){
                i--;
                accept_random_step=false;
                attempt_counter++;
                if (attempt_counter>1000) {
                    cout<<"The chromatin has reached a dead end (chain length = "<<i<<" ) before reaching the set node number"<<Num_of_Nodes<<"\n";
                    exit(EXIT_FAILURE);
                }
                break;
            }
        } // for (int k=0; k<i; k++)
        
        if (accept_random_step) {
            attempt_counter=0;
            
            Node_Position[i][0] = temp_x;
            Node_Position[i][1] = temp_y;
            Node_Position[i][2] = temp_z;
            
            Node_Velocity[i][0]=((double)rand()/(double)RAND_MAX)*2-1;
            Node_Velocity[i][1]=((double)rand()/(double)RAND_MAX)*2-1;
            Node_Velocity[i][2]=((double)rand()/(double)RAND_MAX)*2-1;
            velocity_COM[0] += Node_Velocity[i][0];
            velocity_COM[1] += Node_Velocity[i][1];
            velocity_COM[2] += Node_Velocity[i][2];
        }
    } // for (int i=1; i<Num_of_Nodes; i++)
    
    velocity_COM[0]/=Num_of_Nodes;
    velocity_COM[1]/=Num_of_Nodes;
    velocity_COM[2]/=Num_of_Nodes;
    
    for (int i=0; i<Num_of_Nodes; i++) {
        for (int j=0; j<3; j++) {
            Node_Velocity[i][j] -= velocity_COM[j];
        }
    }
}

void Chromatin::Pack(double min_radius){
    cout<<"\nBeginnig the Chromatin packing\nProgress:\nmin radius = "<<min_radius<<endl;
    int progress=0;
    int MD_packing_Steps=4000;
    double Sphere_Radius=0;
    
    double Max_dist=chromatin_prepack();
    
    if (Max_dist > min_radius) {
        double slope=(min_radius-Max_dist)/MD_packing_Steps;
        
        
        for(int MD_Step=0 ;MD_Step<=MD_packing_Steps ; MD_Step++){
            
            Sphere_Radius=slope*MD_Step+Max_dist;
            
            for (int i=0; i<100; i++) {
                MD_Evolution_beginning(GenConst::Step_Size_In_Fs);
                packing_potential(Sphere_Radius);
                Strong_spring();
                hard_sphere();
                MD_Evolution_end(GenConst::Step_Size_In_Fs);
            }
            if (MD_Step!=0 && MD_Step%10==0) {
                packing_traj();
                if (MD_Step%100==0) {
                    rescale_velocities(0.95);
                }
                
            }
            if (int(100*MD_Step/MD_packing_Steps) > progress){
                cout<<"[ "<<progress<<"% ]\t step: "<<MD_Step*100<<"\r" << std::flush;
                progress+=5;
            }
        }
        
        double velocity_COM[3]={0};
        //    double position_COM[3]={0};
        
        for (int i=0; i<Num_of_Nodes; i++) {
            for (int j=0; j<3; j++) {
                Node_Velocity[i][j]=((double)rand()/(double)RAND_MAX)*2-1;
                velocity_COM[j] += Node_Velocity[i][j];
                //            position_COM[j] += Node_Position[i][j];
                Node_Force[i][j]=0;
            }
        }
        
        velocity_COM[0]/=Num_of_Nodes;
        velocity_COM[1]/=Num_of_Nodes;
        velocity_COM[2]/=Num_of_Nodes;
        //        cout<<"Chromatin com velocity after packing:\n";
        //        cout<<velocity_COM[0]<<"\t"<<velocity_COM[1]<<"\t"<<velocity_COM[2]<<endl;
        //    position_COM[0]/=Num_of_Nodes;
        //    position_COM[1]/=Num_of_Nodes;
        //    position_COM[2]/=Num_of_Nodes;
        
        for (int i=0; i<Num_of_Nodes; i++) {
            for (int j=0; j<3; j++) {
                Node_Velocity[i][j] += -velocity_COM[j];
                //            Node_Position[i][j] += -position_COM[j];
            }
        }
        export_pack(0);
    }
    
}

void Chromatin::packing_potential(double Sphere_Radius){
    double le1,lmin;
    double deltax, deltay, deltaz, Node_distance, force;
    double interaction_strength=0.05;
    
    lmin=Node_radius;
    le1=lmin*4;
    
    for (int i=0; i<Num_of_Nodes; i++) {
        
        deltax=Node_Position[i][0];
        deltay=Node_Position[i][1];
        deltaz=Node_Position[i][2];
        
        double a[3]={deltax, deltay, deltaz};
        double pos=vector_length(a);
        Node_distance=Sphere_Radius-pos;
        if (Node_distance<0) {
            cout<<"Chromatin got out of the confinement."<<endl;
            cout<<"Chromatin index = "<<i<<endl;
            cout<<"Node distance from origin = "<<pos<<endl;
            cout<<"Wall (Sphere_Radius) = "<<Sphere_Radius<<endl;
            cout<<"Node distance from wall = "<<Node_distance<<endl;
            exit(EXIT_FAILURE);
        } else if(Node_distance < le1   &&  Node_distance > lmin  ) {
            double exp_le1=exp(1.0/(Node_distance-le1));
            force = ( (interaction_strength*exp_le1)/(Node_distance-lmin) )*( 1.0/(Node_distance-lmin)+1.0/( (Node_distance-le1)*(Node_distance-le1) ) );
            
            force /= pos;
            
            Node_Force[i][0] += -force*deltax;
            Node_Force[i][1] += -force*deltay;
            Node_Force[i][2] += -force*deltaz;
        }
    }
}

void Chromatin::reset_com_velocity(void){
    double velocity_com[3]={0};
    
    for (int i=0; i<Num_of_Nodes; i++) {
        for (int j=0; j<3; j++) {
            velocity_com[j] += Node_Velocity[i][j];
        }
    }
    velocity_com[0]/=Num_of_Nodes;
    velocity_com[1]/=Num_of_Nodes;
    velocity_com[2]/=Num_of_Nodes;
    
    for (int i=0; i<Num_of_Nodes; i++) {
        for (int j=0; j<3; j++) {
            Node_Velocity[i][j] -= 0.5*velocity_com[j];
        }
    }
}

void Chromatin::rescale_velocities(double scale){
    for (int i=0; i<Num_of_Nodes; i++) {
        for (int j=0; j<3; j++) {
            Node_Velocity[i][j] *= scale;
        }
    }
}

double Chromatin::chromatin_prepack(void){
    
    double Max_dist=0;
    double pos_com[3]={0};
    
    for (int i=0; i<Num_of_Nodes; i++) {
        pos_com[0] += Node_Position[i][0];
        pos_com[1] += Node_Position[i][1];
        pos_com[2] += Node_Position[i][2];
    }
    
    pos_com[0] /= Num_of_Nodes;
    pos_com[1] /= Num_of_Nodes;
    pos_com[2] /= Num_of_Nodes;
    for (int i=0; i<Num_of_Nodes; i++) {
        Node_Position[i][0] -= pos_com[0];
        Node_Position[i][1] -= pos_com[1];
        Node_Position[i][2] -= pos_com[2];
        double a[3]={Node_Position[i][0], Node_Position[i][1], Node_Position[i][2]};
        double temp_dist=vector_length(a);
        if (temp_dist > Max_dist) {
            Max_dist=temp_dist;
        }
    }
    cout<<"Max_dist + Node_radius = "<<Max_dist+2*Node_radius<<endl;
    return Max_dist+2*Node_radius;
}

void Chromatin::export_pack(int MD_step){
    std::ofstream write_resume_file;
    string resume_file_name="Results/Relaxation/Resume_Chromatin_"+std::to_string(index)+"_";
    resume_file_name+=file_time;
    resume_file_name+=".txt";
    write_resume_file.open(resume_file_name.c_str());
    
    write_resume_file<<MD_step<<endl;
    write_resume_file<<Num_of_Nodes<<endl;
    for (int i=0; i<Num_of_Nodes; i++) {
        write_resume_file<<Node_Position[i][0]<<"\t"<<Node_Position[i][1]<<"\t"<<Node_Position[i][2]<<"\n";
        write_resume_file<<Node_Velocity[i][0]<<"\t"<<Node_Velocity[i][1]<<"\t"<<Node_Velocity[i][2]<<"\n";
        write_resume_file<<AB_index[i]<<"\n";
        //Node_force=0
    }
}
