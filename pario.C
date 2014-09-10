#include <stdio.h>
#include <iostream>
#include "pario.h"

#ifdef MOLPRO_I8
#include "global/CxOutputStream.h"
blockout Bout(&xout);
blockerr Berr(&xerr);
#else
blockout Bout;
blockerr Berr;
#endif

std::ostream &bout = *(Bout.outstream);
std::ostream &berr = *(Berr.errstream);
