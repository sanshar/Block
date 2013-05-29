#ifndef GA_INPUT_H
#define GA_INPUT_H

#include <fstream>
#include <cstring>
#include <boost/serialization/serialization.hpp>
using namespace std;

namespace genetic
{
  enum SELECT_TYPE { GAUSS, BOLTZMANN, UNIFORM };

  struct GAInput
  {
  private:
    friend class boost::serialization::access;
    template <class Archive>
    void serialize(Archive& ar, const unsigned int version)
    {
      ar & max_community;
      ar & max_generation;
      ar & max_cells;
      ar & thre_cloning;
      ar & thre_mutation;
      ar & max_elite;
      ar & scaling;
      ar & exponent;
      ar & random_seed;
      ar & select;
      ar & fiedler;
    }

  public:
    int max_community;
    int max_generation;
    int max_cells;
    double thre_cloning;
    double thre_mutation;
    int max_elite;
    double scaling;
    double exponent;
    unsigned int random_seed;
    int fiedler;
    SELECT_TYPE select;

    GAInput(void);
    GAInput(ifstream& config);
    void Configure(ifstream& config);
  };
};

#endif
