#ifndef GENE_H
#define GENE_H

#include <iostream>
#include <ostream>
#include <iomanip>
#include <vector>
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/vector.hpp>
#include "genetic_utils.h"
using namespace std;

namespace genetic
{
  class Gene
  {
  private:
    friend class boost::serialization::access;
    template <class Archive>
    void serialize(Archive& ar, const unsigned int version)
    {
      ar & m_sequence;
    }

  private:
    vector<int> m_sequence;
    static int m_length;

  public:
    static int& Length(void) { return m_length; }
    const static Gene Random(void)
    {
      Gene gen(RandomSequence(m_length));
      return gen;
    }

    const vector<int>& Sequence(void) const { return m_sequence; }
    Gene(void) { }
   ~Gene(void) { }
    Gene(const vector<int>& sequence)
    {
      if(m_length == 0) m_length = sequence.size();
      m_sequence = sequence;
    }
    Gene(const Gene& other)
    {
      m_sequence = other.m_sequence;
    }
    Gene& operator= (const Gene& other)
    {
      m_sequence = other.m_sequence;
      return *this;
    }
    Gene  operator* (const Gene& other) const
    {
      return Gene(CrossOver(m_sequence, other.m_sequence));
    }
    void pMutate(void) // point mutation
    {
      m_sequence = PointMutation(m_sequence);
    }
    void gMutate(void) // global mutation
    {
      m_sequence = GlobalMutation(m_sequence);
    }
    friend ostream& operator<< (ostream& ost, const Gene& gene)
    {
      ost << "Gene::Sequence = ";
      for(int i = 0; i < m_length; ++i) ost << setw(3) << gene.m_sequence[i] + 1 << ","; ost << flush;
      return ost;
    }
  };

};

#endif
