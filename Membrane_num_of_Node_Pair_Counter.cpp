#include "Membrane.h"

int Membrane::Membrane_num_of_Node_Pair_Counter()
{
    int bondslist[10*Membrane_num_of_Nodes][2];
    for (int j=0 ; j< Membrane_num_of_Nodes*10 ; j++)
    {
        bondslist[j][0]=-1;
        bondslist[j][1]=-1;
    }
    
    int temp_Membrane_num_of_Node_Pairs=0;
    int temp_Membrane_triangle_Node_A, temp_Membrane_triangle_Node_B, temp_Membrane_triangle_Node_C;
    
    int repeatednumber1=0;
    int repeatednumber2=0;
    int repeatednumber3=0;
    
    for(int i=0;i<Membrane_triangle_list.size();i++)
    {
        temp_Membrane_triangle_Node_A= Membrane_triangle_list[i][0];
        temp_Membrane_triangle_Node_B= Membrane_triangle_list[i][1];
        temp_Membrane_triangle_Node_C= Membrane_triangle_list[i][2];
        
        for(int j=0;j<10*Membrane_num_of_Nodes;j++)
        {
            if(  ( bondslist[j][0]==temp_Membrane_triangle_Node_A &  bondslist[j][1]==temp_Membrane_triangle_Node_B )  || ( bondslist[j][0]==temp_Membrane_triangle_Node_B &  bondslist[j][1]==temp_Membrane_triangle_Node_A )    )
            {
                repeatednumber1=1;
            }
            
            if(  ( bondslist[j][0]==temp_Membrane_triangle_Node_B &  bondslist[j][1]==temp_Membrane_triangle_Node_C )  || ( bondslist[j][0]==temp_Membrane_triangle_Node_C &  bondslist[j][1]==temp_Membrane_triangle_Node_B )    )
            {
                repeatednumber2=1;
            }
            
            if(  ( bondslist[j][0]==temp_Membrane_triangle_Node_A &  bondslist[j][1]==temp_Membrane_triangle_Node_C )  || ( bondslist[j][0]==temp_Membrane_triangle_Node_C &  bondslist[j][1]==temp_Membrane_triangle_Node_A )    )
            {
                repeatednumber3=1;
            }
        }
        
        if(repeatednumber1==0)
        {
            bondslist[temp_Membrane_num_of_Node_Pairs][0]=temp_Membrane_triangle_Node_A;       //note that first node store in i and second in i+numofbonds  ---Membrane_Node_Pair_list[2*numofbonds]
            bondslist[temp_Membrane_num_of_Node_Pairs][1]=temp_Membrane_triangle_Node_B;
            temp_Membrane_num_of_Node_Pairs++;
            
        }
        
        if(repeatednumber2==0)
        {
            bondslist[temp_Membrane_num_of_Node_Pairs][0]=temp_Membrane_triangle_Node_B;       //note that first node store in i and second in i+numofbonds  ---Membrane_Node_Pair_list[2*numofbonds]
            bondslist[temp_Membrane_num_of_Node_Pairs][1]=temp_Membrane_triangle_Node_C;
            temp_Membrane_num_of_Node_Pairs++;
            
        }
        
        if(repeatednumber3==0)
        {
            bondslist[temp_Membrane_num_of_Node_Pairs][0]=temp_Membrane_triangle_Node_A;       //note that first node store in i and second in i+numofbonds  ---Membrane_Node_Pair_list[2*numofbonds]
            bondslist[temp_Membrane_num_of_Node_Pairs][1]=temp_Membrane_triangle_Node_C;
            temp_Membrane_num_of_Node_Pairs++;
            
        }
        
        repeatednumber1=0;
        repeatednumber2=0;
        repeatednumber3=0;
    }
    
    //    cout<<"Outer_Membrane_num_of_Node_Pairs_test=\t"<<Outer_Membrane_num_of_Node_Pairs<<endl;
    //    exit (EXIT_FAILURE);
    //*******************************************************************************************************
    /*BUG
     |---\   |    |  /---\
     |    |  |    |  |
     |---<   |    |  |  -\
     |    |  |    |  |   |
     |---/   \----/  \---/
     */
    //*******************************************************************************************************
    //***************** Potential BUG: This counter gives exactley the same number as the *******************
    //***************** 'Membrane_triangle_pair_counter'. I think we can use that number, *******************
    //***************** but we have to check it first *******************************************************
    //*******************************************************************************************************
	return temp_Membrane_num_of_Node_Pairs;
	
    
}
