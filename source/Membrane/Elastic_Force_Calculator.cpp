#include "Membrane.h"

void Membrane::Elastic_Force_Calculator(double theta_0)
{
	
    Total_Potential_Energy=0.0;
	if (spring_model==GenConst::potential.Model["FENE"]) {FENE_log();}
    if (spring_model==3) {custom();}
    
    if (Bending_coefficient!=0) {
        Bending_potetial_2(0);
    }
    
}
