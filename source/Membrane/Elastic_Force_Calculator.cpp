#include "Membrane.h"

void Membrane::Elastic_Force_Calculator(double theta_0)
{
	
    Total_Potential_Energy=0.0;
	if (spring_model==1) {log_barrier();}
    if (spring_model==2) {Hookian();}
    if (spring_model==3) {FENE();}
    if (spring_model==4) {Relaxation_potential();}
    
    if (Bending_coefficient!=0) {Bending_potetial_2(0);}
    
}
