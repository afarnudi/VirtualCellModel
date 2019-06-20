//
//  check.cpp
//  Mem
//
//  Created by Ali Farnudi on 19/08/2018.
//  Copyright © 2018 Ali Farnudi. All rights reserved.
//

#include "Membrane.h"

void Membrane::Relax_1(void){
    check();
    if (Relaxation){
        cout<<"\nBeginnig the Relaxation\nProgress:\n";
        double temp_Bending_coefficient=Bending_coefficient;
        Bending_coefficient=350*GenConst::MD_T*GenConst::K;
        int Relaxation_progress=0;
        if (Relaxation_Process_Model==1){
            double relax_temp= GenConst::MD_T;
            for(int MD_Step=0 ;MD_Step<MD_num_of_Relaxation_steps ; MD_Step++){
                MD_Evolution_beginning(GenConst::MD_Time_Step);
                Hookian();
                Bending_potetial_2(0);
                MD_Evolution_end(GenConst::MD_Time_Step);
                if(MD_Step%GenConst::MD_traj_save_step==0){
                    relaxation_traj();}
                if(MD_Step%50000==0){
                    Thermostat_2(relax_temp);
                    relax_temp=relax_temp*0.8;
                }
                if (int(100*MD_Step/MD_num_of_Relaxation_steps)>Relaxation_progress){
                    cout<<"[ "<<Relaxation_progress<<"% ]\t step: "<<MD_Step<<"\r" << std::flush;
                    Relaxation_progress+=5;
                }
            } //End of for(int MD_Step=0 ;MD_Step<=GenConst::MD_num_of_Relaxation_steps ; MD_Step++)
            calculate_mesh_properties();
        }// End of if (Membranes[i].return_Relaxation_Prosses_Model()==1)
        if (Relaxation_Process_Model==2){
            node_distance_correction();
            cout<<"Level2"<<endl;
            int Relaxation_progress=get_correction_progress();
            double relax_temp= GenConst::MD_T;
            for(int MD_Step=0 ;MD_Step<=MD_num_of_Relaxation_steps ; MD_Step++){
                MD_Evolution_beginning(GenConst::MD_Time_Step);
                Relaxation_potential();
                Bending_potetial_2(0);
                MD_Evolution_end(GenConst::MD_Time_Step);
                if(MD_Step%GenConst::MD_traj_save_step==0){
                    relaxation_traj();
                }
                if(MD_Step%50000==0){
                    Thermostat_2(relax_temp);
                    relax_temp=relax_temp*0.8;
                    
                    if (int(100*MD_Step/(MD_num_of_Relaxation_steps + MD_correction_steps * 100))>Relaxation_progress){
                        cout<<"[ "<<Relaxation_progress<<"% ]\t step: "<<MD_Step<<"\r" << std::flush;
                        Relaxation_progress+=5;
                    }
                } // End of for(int MD_Step=0 ;MD_Step<=GenConst::MD_num_of_Relaxation_steps ; MD_Step++)
                calculate_mesh_properties();
            } // End of if (Membranes[i].return_Relaxation_Prosses_Model()==2)
            Bending_coefficient=temp_Bending_coefficient;
            export_relaxed(MD_num_of_Relaxation_steps);
        }
    } //End if (Relaxation)
    
} // End Relax_1(void)



void Membrane::node_distance_correction(void){
    correction_progress=0;
    cout<<"level1"<<endl;
    // double MD_relax_Steps_1=2000;//2000;
    double slope=(Max_node_pair_length/1.8-Min_node_pair_length)/MD_correction_steps, min=Min_node_pair_length;
    // double temp_Bending_coefficient=Bending_coefficient;
    //Bending_coefficient=70*GenConst::MD_T*GenConst::K;
    for(int MD_Step=0 ;MD_Step<=MD_correction_steps ; MD_Step++){
        //Setting the min angle of triangles to 20 dgrees or pi/9
        Min_node_pair_length=slope*MD_Step+min;
        for (int i=0; i<100; i++) {
            MD_Evolution_beginning(GenConst::MD_Time_Step);
            Relaxation_potential();
            Bending_potetial_2(0);
            MD_Evolution_end(GenConst::MD_Time_Step);
        }
        if (int(100*MD_Step/(MD_correction_steps*100 + MD_num_of_Relaxation_steps))>correction_progress){
            cout<<"Hi"<<endl;
            cout<<"[ "<<correction_progress<<"% ]\t step: "<<MD_Step<<"\r" << std::flush;
            correction_progress+=5;
        }
    } //End of for (int MD_Step=0 ;MD_Step<=MD_num_of_steps ; MD_Step++)
    //  Bending_coefficient=temp_Bending_coefficient;
}//end


