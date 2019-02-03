#include "Membrane.h"

void Membrane::Elastic_Force_Calculator(double theta_0)
{
	
    Total_Potential_Energy=0.0;
	if (spring_model==1) {potential_1();}
    if (spring_model==2) {potential_2();}//Spring
    if (spring_model==3) {FENE();}

    Bending_potetial_2(0);
}
