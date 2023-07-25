#ifndef GENERAL_CLASS_FUNCTIONS_H
#define GENERAL_CLASS_FUNCTIONS_H

#include <vector>

#include "Membrane.h"
#include "Actin.h"
#include "ECM.h"
#include "Chromatin.h"
#include "Interaction_table.hpp"

void Update_classes(std::vector<Membrane>  &membranes,
                    std::vector<Actin>     &actins,
                    std::vector<ECM>       &ecms,
                    std::vector<Chromatin> &chromatins,
                    double                  time,
                    MyAtomInfo*             all_atoms);


using std::set;

bool check_for_membrane_update(vector<Membrane>    &membranes,
                               double               time);

void updateOpenMMforces(vector<Membrane>                &membranes,
                        vector<Chromatin>                chromos,
                        MyOpenMMData*                    omm,
                        double                           time,
                        MyAtomInfo                       atoms[],
                        Bonds*                           bonds,
                        vector<set<int> >               &membrane_set,
                        NonBondInteractionMap            interaction_map);

void expand(vector<Chromatin>                chromos,
            MyOpenMMData*                    omm);

#endif // GENERAL_CLASS_FUNCTIONS_H
