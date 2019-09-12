#ifndef GENERAL_CLASS_FUNCTIONS_H
#define GENERAL_CLASS_FUNCTIONS_H

#include <vector>

#include "Membrane.h"
#include "Actin.h"
#include "ECM.h"
#include "Chromatin.h"

void Export_classes_for_resume(std::vector<Membrane>  &membranes,
                               std::vector<Actin>     &actins,
                               std::vector<ECM>       &ecms,
                               std::vector<Chromatin> &chromatins,
                               double                  time,
                               MyAtomInfo*             all_atoms);


using std::set;

bool check_for_membrane_update(vector<Membrane>    &membranes,
                               double               time);

void updateOpenMMforces(vector<Membrane>       &membranes,
                        MyOpenMMData*          omm,
                        double                 time,
                        MyAtomInfo             atoms[],
                        Bonds*                 bonds,
                        vector<set<int> >      &membrane_set,
                        vector<set<int> >      &actin_set,
                        vector<set<int> >      &ecm_set,
                        vector<set<int> >      &chromatin_set,
                        vector<vector<int> >   interaction_map);

#endif // GENERAL_CLASS_FUNCTIONS_H
