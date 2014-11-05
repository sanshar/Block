#ifndef __SORT_NPDM_ARRAY_H
#define __SORT_NPDM_ARRAY_H

#include <string>
#include "molcas_types.h"

void sort1pdm (FORTINT N_act, FORTINT iRoot, FORTINT jRoot);
void sort2pdm (FORTINT N_act, FORTINT iRoot, FORTINT jRoot);
void sort3pdm (FORTINT N_act, FORTINT iRoot, FORTINT jRoot);

#endif // __SORT_NPDM_ARRAY_H
