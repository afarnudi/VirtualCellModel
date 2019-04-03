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
            int Relaxation_progress=0; 
           if (Relaxation_Process_Model==1){
            double relax_temp= GenConst::MD_T;  
            for(int MD_Step=0 ;MD_Step<=GenConst::MD_num_of_Relaxation_steps ; MD_Step++){
                MD_Evolution_beginning(GenConst::MD_Time_Step);
                Elastic_Force_Calculator(0);
                MD_Evolution_end(GenConst::MD_Time_Step);
                if(MD_Step%100==0){
                    Thermostat_2(relax_temp);
                    relax_temp=relax_temp/2;}
                if (int(100*MD_Step/GenConst::MD_num_of_Relaxation_steps)>Relaxation_progress){
            cout<<"[ "<<Relaxation_progress<<"% ]\t step: "<<MD_Step<<"\r" << std::flush;
            Relaxation_progress+=5;}
            // I am not sure that this will work fine in the case of having more than one membrane.
              //  if (MD_Step%GenConst::MD_traj_save_step == 0){
                //Trajectory << num_of_elements<<endl;
                //Trajectory << " nodes  "<<endl;
                //string label="Membrane_"+to_string(i);
               // write_traj(traj_file_name, label);
                }//End of if (MD_Step%GenConst::MD_traj_save_step == 0)
            
           } //End of for(int MD_Step=0 ;MD_Step<=GenConst::MD_num_of_Relaxation_steps ; MD_Step++)
           calculate_mesh_properties();
           }// End of if (Membranes[i].return_Relaxation_Prosses_Model()==1)
           if (Relaxation_Process_Model==2){
               node_distance_correction();
               int Relaxation_progress=return_correction_progress();
                double relax_temp= GenConst::MD_T;
                for(int MD_Step=0 ;MD_Step<=GenConst::MD_num_of_Relaxation_steps ; MD_Step++){
                MD_Evolution_beginning(GenConst::MD_Time_Step);
                Elastic_Force_Calculator(0);
                MD_Evolution_end(GenConst::MD_Time_Step);
                if(MD_Step%100==0){
                    Thermostat_2(relax_temp);
                    relax_temp=relax_temp/2;}
                if (int(100*MD_Step/(GenConst::MD_num_of_Relaxation_steps + GenConst::MD_correction_steps * 100))>Relaxation_progress){
            cout<<"[ "<<Relaxation_progress<<"% ]\t step: "<<MD_Step<<"\r" << std::flush;
            Relaxation_progress+=5;}
        } // End of for(int MD_Step=0 ;MD_Step<=GenConst::MD_num_of_Relaxation_steps ; MD_Step++)
        calculate_mesh_properties();
        } // End of if (Membranes[i].return_Relaxation_Prosses_Model()==2)
    } 
    

void Membrane::Relax_2(void){
    node_distance_correction();
    calculate_mesh_properties();
}
