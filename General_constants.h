//
//  General_constants.h
//  Cell-Durotaxis
//
//  Created by Ali Farnudi on 26/08/2017.
//  Copyright © 2017 Ali Farnudi. All rights reserved.
//

#ifndef General_constants_h
#define General_constants_h

#define MD_num_of_steps  5000//35000// number of MD stps
#define savingstep    1000//The step on which the trajector of the membrane is saved.
#define MD_Time_Step     0.001 // time length of steps in MD
#define KT     1.0  // KT the quanta of energy
#define pi     3.141592 // clear !
#define RunThermostatePerstep   100 //
#define Node_radius 1.0 // Interacion range of the nodes: use to be 1.0 !
#define mcstep 1

#define mcstep    1
#define fluidity       0.002 //Used in the MC step

#define Lbox  1000.0///    (size of square periodic box-1)
#define Periodic_condtion_status 0.0 //status 0.0 = off (The Periodic update will not be executed in the 'Main MD' loop). status = 1.0 = on

#endif /* General_constants_h */
