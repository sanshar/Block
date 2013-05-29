#include <string>
#include "GAInput.h"
using namespace std;

genetic::GAInput::GAInput(void)
{
  max_community  = 32;
  max_generation = 10000;
  max_cells      = 0;
  thre_cloning   = 0.90;
  thre_mutation  = 0.10;
  max_elite      = 1;
  scaling        = 1.0;
  exponent       = 2.0;
  select         = GAUSS;
  random_seed    = 1414; //setting default seed
  fiedler        = 1;
}

genetic::GAInput::GAInput(ifstream& config)
{
  Configure(config);
}

void genetic::GAInput::Configure(ifstream& config)
{
  string entry;
  while(config >> entry)
  {
    if(entry == "maxcomm")  config >> max_community;
    if(entry == "maxgen")   config >> max_generation;
    if(entry == "maxcell")  config >> max_cells;
    if(entry == "cloning")  config >> thre_cloning;
    if(entry == "mutation") config >> thre_mutation;
    if(entry == "elite")    config >> max_elite;
    if(entry == "scale")    config >> scaling;
    if(entry == "exponent") config >> exponent;
    if(entry == "seed")     config >> random_seed;
    if(entry == "fiedler")  config >> fiedler;
    if(entry == "method")
    {
      config >> entry;
      if(entry == "gauss")    select = GAUSS;
      if(entry == "boltzmann") select = BOLTZMANN;
      if(entry == "uniform")  select = UNIFORM;
    }
  }
}
