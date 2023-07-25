#include "Membrane.h"

void Membrane::Elastic_Force_Calculator(double theta_0)
{
	
    Total_Potential_Energy=0.0;
	if (spring_model==potentialModelIndex.Model["FENE"]) {FENE_log();}
    if (spring_model==3) {custom();}
    
    if (dihedral_bending_coefficient!=0) {
        Bending_potetial_2(0);
    }
    
}

void Membrane::Mechanical_Energy_calculator()
{
    
    
    Total_Bond_Energy=0;
    Total_Bending_Energy=0;
    
    if (spring_model==potentialModelIndex.Model["FENE"]) {
        string errorMessage = TWARN;
        errorMessage +="Error: Membrane mechanical energy calculator:  FENE potential is underdevelopment and is not supported at the moment.\n";
        errorMessage +=TRESET;
        throw std::runtime_error(errorMessage);
    } else if (spring_model==potentialModelIndex.Model["Harmonic"]) {
        harmonic_potential_calculator();
    } else if (spring_model==potentialModelIndex.Model["None"]){
        Total_Bond_Energy=0;
    }
    
    if (dihedral_bending_model==potentialModelIndex.Model["Dihedral"]) {
        calculating_dihedral_energy();
    } else if (dihedral_bending_model==potentialModelIndex.Model["None"]){
        Total_Bending_Energy=0;
    }
    
    Total_Potential_Energy = Total_Bending_Energy + Total_Bond_Energy;
//    cout<< "Total_Bending_Energy "<<Total_Bending_Energy<<"  Total_Bond_Energy "<<Total_Bond_Energy<<"  Total_Potential_Energy "<<Total_Potential_Energy<<endl;
}
