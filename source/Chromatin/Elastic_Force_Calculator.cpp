#include "Chromatin.h"

void Chromatin::Elastic_Force_Calculator(void)
{
    Total_Potential_Energy=0.0;
	if (spring_model==0) {potential_1();}
    if (spring_model==1) {FENE();}
    if (spring_model==2) {Strong_spring();}
    hard_sphere();
}
