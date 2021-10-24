//
//  General_constants.h
//  Cell-Durotaxis
//
//  Created by Ali Farnudi on 26/08/2017.
//  Copyright © 2017 Ali Farnudi. All rights reserved.
//

#ifndef OpenMM_WCA_hpp
#define OpenMM_WCA_hpp

#include "OpenMM_structs.h"
#include "Interaction_table.hpp"


void init_WCAFC_interaction(vector<OpenMM::CustomNonbondedForce*> &WCAFCs,
                            const MyAtomInfo                      atoms[],
                            set<int>                              set_1,
                            set<int>                              set_2,
                            string                                set_1_name,
                            string                                set_2_name,
                            bool                                  use_max_radius);

void init_WCAFC_interaction(vector<OpenMM::CustomNonbondedForce*> &WCAFCs,
                            const MyAtomInfo                      atoms[],
                            set<int>                              set_1,
                            set<int>                              set_2,
                            int                                   set_1_index,
                            int                                   set_2_index,
                            string                                set_1_name,
                            string                                set_2_name,
                            bool                                  use_max_radius);
void init_WCAFC_interaction(vector<OpenMM::CustomNonbondedForce*> &WCAFCs,
                          const MyAtomInfo                      atoms[],
                          vector<set<int> >                     set_1,
                          vector<set<int> >                     set_2,
                          int                                   set_1_index,
                          int                                   set_2_index,
                          string                                set_1_name,
                          string                                set_2_name);
void init_WCAFC_interaction(vector<OpenMM::CustomNonbondedForce*> &WCAFCs,
                          const MyAtomInfo                       atoms[],
                          vector<vector<set<int> > >             set_1,
                          vector<set<int> >                      set_2,
                          int                                    set_1_index,
                          int                                    set_2_index,
                          string                                 set_1_name,
                          string                                 set_2_name,
                            bool                                   use_max_radius);
void init_WCAFC_interaction(vector<OpenMM::CustomNonbondedForce*> &WCAFCs,
                          const MyAtomInfo                       atoms[],
                          vector<vector<std::set<int> > >             set_1,
                          vector<vector<std::set<int> > >             set_2,
                          int                                    set_1_index,
                          int                                    set_2_index,
                          string                                 set_1_name,
                          string                                 set_2_name);








void init_WCA_interaction(vector<OpenMM::CustomNonbondedForce*> &WCAs,
                          const MyAtomInfo                      atoms[],
                          set<int>                              set_1,
                          set<int>                              set_2,
                          string                                set_1_name,
                          string                                set_2_name,
                          bool                                  use_max_radius);

void init_WCA_interaction(vector<OpenMM::CustomNonbondedForce*> &WCAs,
                          const MyAtomInfo                      atoms[],
                          set<int>                              set_1,
                          set<int>                              set_2,
                          int                                   set_1_index,
                          int                                   set_2_index,
                          string                                set_1_name,
                          string                                set_2_name,
                          bool                                  use_max_radius);

void init_WCA_interaction(vector<OpenMM::CustomNonbondedForce*> &WCAs,
                          const MyAtomInfo                      atoms[],
                          vector<set<int> >                     set_1,
                          vector<set<int> >                     set_2,
                          int                                   set_1_index,
                          int                                   set_2_index,
                          string                                set_1_name,
                          string                                set_2_name);

void init_WCA_interaction(vector<OpenMM::CustomNonbondedForce*> &WCAs,
                          const MyAtomInfo                      atoms[],
                          vector<vector<std::set<int> > >       set_1,
                          vector<std::set<int> >                set_2,
                          int                                   set_1_index,
                          int                                   set_2_index,
                          string                                set_1_name,
                          string                                set_2_name,
                          bool                                  use_max_radius);

void init_WCA_interaction(vector<OpenMM::CustomNonbondedForce*> &WCAs,
                          const MyAtomInfo                       atoms[],
                          vector<vector<std::set<int> > >             set_1,
                          vector<vector<std::set<int> > >             set_2,
                          int                                    set_1_index,
                          int                                    set_2_index,
                          string                                 set_1_name,
                          string                                 set_2_name);



double get_WCA_CutoffDistance(double    radius_1,
                              double    radius_2,
                              bool      use_max_radius);

double get_WCA_sigma(double    radius_1,
                     double    radius_2,
                     bool      use_max_radius);

#endif
