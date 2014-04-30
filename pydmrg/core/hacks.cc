/*
 * File: hacks.cc
 *
 * Some hack functions to cheat cython compiler and Block global variables
 */

#include <string>
#include "config.h"
#include "global.h"
#include "Symmetry.h"
#include "input.h"
#include "orbstring.h"
#include "hacks.h"

using namespace SpinAdapted;


void init_dmrginp(char *conf)
{
    v_1.rhf = true;
    v_2.rhf = true;
    std::string configFile(conf);
    dmrginp = Input(configFile);
    //dmrginp.initialize_defaults();
    dmrginp.initCumulTimer();
    Orbstring::init(dmrginp.slater_size());
}

int get_last_site_id()
{
    return dmrginp.last_site() - 1;
}

void initialize_default_dmrginp(char *fcidump, std::string& prefix, std::string& inpsym)
{
    dmrginp.initCumulTimer();
    v_1.rhf = true;
    v_2.rhf = true;
    //dmrginp.initialize_defaults();

    // TODO: remove this, use more natrual way to handle Block's IO
    std::string *save_prefix = const_cast<std::string *>(&dmrginp.save_prefix());
    std::string *load_prefix = const_cast<std::string *>(&dmrginp.load_prefix());
    *save_prefix = prefix;
    *load_prefix = prefix;

    sym = inpsym; // FIXME: the arg inpsym of InitialiseTable has no effects
    if (inpsym != "c1") {
        Symmetry::InitialiseTable(inpsym);
    }

    std::string orbfile(fcidump);
    dmrginp.readorbitalsfile(orbfile, v_1, v_2);

    // class Slater and OrbString store its size in global scope in OrbString.C
    Orbstring::init(dmrginp.slater_size());
}
