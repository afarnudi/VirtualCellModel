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
                               MyAtomInfo*             all_atoms,
                               int                    &atom_count,
                               int                    &bond_count,
                               int                    &dihe_count);

#endif // GENERAL_CLASS_FUNCTIONS_H
