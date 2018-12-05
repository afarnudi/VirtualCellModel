//#include "Membrane.h"
//
//void Membrane::Membrane_num_of_Node_Pair_Counter_2()
//
//{
//    Membrane_Node_Pair_list.resize(Membrane_num_of_Node_Pairs);
//    for(int i=0 ; i<Membrane_num_of_Node_Pairs ; i++)
//    {
//        Membrane_Node_Pair_list[i].resize(2);
//    }
//    for (int j=0 ; j< Membrane_num_of_Node_Pairs ; j++)
//    {
//        Membrane_Node_Pair_list[j][0]=-1;
//        Membrane_Node_Pair_list[j][1]=-1;
//    }
//    
//    int temp_Membrane_num_of_Node_Pairs=0;
//    int temp_Membrane_triangle_Node_A, temp_Membrane_triangle_Node_B, temp_Membrane_triangle_Node_C;
//    
//    int repeatednumber1=0;
//    int repeatednumber2=0;
//    int repeatednumber3=0;
//    for(int i=0;i<Membrane_triangle_list.size();i++)
//    {
//        temp_Membrane_triangle_Node_A= Membrane_triangle_list[i][0];
//        temp_Membrane_triangle_Node_B= Membrane_triangle_list[i][1];
//        temp_Membrane_triangle_Node_C= Membrane_triangle_list[i][2];
//        
//        for(int j=0;j<Membrane_num_of_Node_Pairs;j++)
//        {
//            if(  ( Membrane_Node_Pair_list[j][0]==temp_Membrane_triangle_Node_A &  Membrane_Node_Pair_list[j][1]==temp_Membrane_triangle_Node_B )  || ( Membrane_Node_Pair_list[j][0]==temp_Membrane_triangle_Node_B &  Membrane_Node_Pair_list[j][1]==temp_Membrane_triangle_Node_A )    )
//            {
//                repeatednumber1=1;
//            }
//            
//            if(  ( Membrane_Node_Pair_list[j][0]==temp_Membrane_triangle_Node_B &  Membrane_Node_Pair_list[j][1]==temp_Membrane_triangle_Node_C )  || ( Membrane_Node_Pair_list[j][0]==temp_Membrane_triangle_Node_C &  Membrane_Node_Pair_list[j][1]==temp_Membrane_triangle_Node_B )    )
//            {
//                repeatednumber2=1;
//            }
//            
//            if(  ( Membrane_Node_Pair_list[j][0]==temp_Membrane_triangle_Node_A &  Membrane_Node_Pair_list[j][1]==temp_Membrane_triangle_Node_C )  || ( Membrane_Node_Pair_list[j][0]==temp_Membrane_triangle_Node_C &  Membrane_Node_Pair_list[j][1]==temp_Membrane_triangle_Node_A )    )
//            {
//                repeatednumber3=1;
//            }
//        }
//        
//        if(repeatednumber1==0)
//        {
//            Membrane_Node_Pair_list[temp_Membrane_num_of_Node_Pairs][0]=temp_Membrane_triangle_Node_A;       //note that first node store in i and second in i+numofbonds  ---Membrane_Node_Pair_list[2*numofbonds]
//            Membrane_Node_Pair_list[temp_Membrane_num_of_Node_Pairs][1]=temp_Membrane_triangle_Node_B;
//            
//            temp_Membrane_num_of_Node_Pairs++;
//            
//        }
//        
//        if(repeatednumber2==0)
//        {
//            Membrane_Node_Pair_list[temp_Membrane_num_of_Node_Pairs][0]=temp_Membrane_triangle_Node_B;       //note that first node store in i and second in i+numofbonds  ---Membrane_Node_Pair_list[2*numofbonds]
//            Membrane_Node_Pair_list[temp_Membrane_num_of_Node_Pairs][1]=temp_Membrane_triangle_Node_C;
//            
//            temp_Membrane_num_of_Node_Pairs++;
//            
//        }
//        
//        if(repeatednumber3==0)
//        {
//            Membrane_Node_Pair_list[temp_Membrane_num_of_Node_Pairs][0]=temp_Membrane_triangle_Node_A;       //note that first node store in i and second in i+numofbonds  ---Membrane_Node_Pair_list[2*numofbonds]
//            Membrane_Node_Pair_list[temp_Membrane_num_of_Node_Pairs][1]=temp_Membrane_triangle_Node_C;
//            
//            temp_Membrane_num_of_Node_Pairs++;
//            
//        }
//        
//        repeatednumber1=0;
//        repeatednumber2=0;
//        repeatednumber3=0;
//    }
//    
//}
//
//
