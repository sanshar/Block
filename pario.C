#include <stdio.h>
#include <iostream>
#include "pario.h"

#ifdef MOLPRO
#include "global/CxOutputStream.h"
blockout Bout(0,&xout);
blockerr Berr(0,&xerr);
#else
blockout Bout;
blockerr Berr;
#endif

std::ostream &bout = *(Bout.outstream);
std::ostream &berr = *(Berr.errstream);
