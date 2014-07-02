/*
 * File:
 */
#include <string.h>
#include "config.h"
#include "StateInfo.h"

using namespace SpinAdapted;

int save_stateinfo(char *filesi, StateInfo *si)
{
    std::ofstream ofs(filesi, std::ios::binary);
    boost::archive::binary_oarchive save_state(ofs);

    save_state << *si;
    ofs.close();
    return 0;
}

int load_stateinfo(char *filesi, StateInfo *si)
{
    std::ifstream ifs(filesi, std::ios::binary);
    boost::archive::binary_iarchive load_state(ifs);

    load_state >> *si;
    ifs.close();
    return 0;
}

std::vector<int> *x_StateInfo_quantaMap(StateInfo *s, int lquanta_id,
                                        int rquanta_id)
{
    return &s->quantaMap(lquanta_id, rquanta_id);
}

char *x_StateInfo_allowedQuanta(StateInfo *s, int lquanta_id,
                                int rquanta_id)
{
    return &s->allowedQuanta(lquanta_id, rquanta_id);
}

// memcpy here since we don't expose the ObjectMatrix class
int get_whole_StateInfo_allowedQuanta(StateInfo *s, char *tftab)
{
    std::vector<char>& rep = s->allowedQuanta.rep;
    memcpy(tftab, &rep[0], rep.size()*sizeof(char));
    return 0;
}

// part of StatInfo::CollectQuanta
// a = a+b
void union_StateInfo_quanta(StateInfo *a, StateInfo *b)
{
    std::vector<SpinQuantum> duplicateQuanta = b->quanta;
    sort (duplicateQuanta.begin (), duplicateQuanta.end ());
    unique_copy (duplicateQuanta.begin (), duplicateQuanta.end (),
                 back_inserter (a->quanta));
    sort (a->quanta.begin (), a->quanta.end ());

    a->quantaStates.resize (a->quanta.size ());
    a->oldToNewState.resize (a->quanta.size ());

    for (int i = 0; i < a->quanta.size (); ++i) {
        for (int j = 0; j < b->quanta.size (); ++j) {
            if (b->quanta [j] == a->quanta [i]) {
                a->quantaStates [i] += b->quantaStates [j];
                // a->oldToNewState may have problem.  oldToNewState
                // saves the index for unMap, new.unMap[oldToNew[i]]
                // see wavefunction.C, BaseOperator.C
                a->oldToNewState [i].push_back (j);
            }
        }
    }

    a->totalStates = accumulate(a->quantaStates.begin (),
                                a->quantaStates.end (), 0);

    // unBlockedIndex is only assigned in CollectQuanta
    a->UnBlockIndex();
}
