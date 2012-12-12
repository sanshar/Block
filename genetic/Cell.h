#ifndef CELL_H
#define CELL_H

#include <iostream>
#include <boost/function.hpp>
#include <boost/serialization/serialization.hpp>
#include "Gene.h"
using namespace std;

namespace genetic
{
  // GA Individual defined by class Cell
  class Cell
  {
  private:
    friend class boost::serialization::access;
    template <class Archive>
    void serialize(Archive& ar, const unsigned int version)
    {
      ar & m_gen;
      ar & m_fit;
    }

  private:
    Gene   m_gen;
    double m_fit;

  public:
    static boost::function<double(const Gene&)> Evaluate;

    Cell(void) { }
   ~Cell(void) { }

    Cell(const Gene& gen) { Create(gen); }

    Cell(const Cell& other)
    {
      m_gen = other.m_gen;
      m_fit = other.m_fit;
    }

    Cell& operator= (const Cell& other)
    {
      m_gen = other.m_gen;
      m_fit = other.m_fit;
      return *this;
    }

    void Create(const Gene& gen)
    {
      m_gen = gen;
      m_fit = Evaluate(gen);
    }

    inline bool operator== (const Cell& other) const { return m_fit == other.m_fit; }
    inline bool operator!= (const Cell& other) const { return m_fit != other.m_fit; }
    inline bool operator<  (const Cell& other) const { return m_fit <  other.m_fit; }
    inline bool operator>  (const Cell& other) const { return m_fit >  other.m_fit; }

    const Gene& Gen(void) const { return m_gen; }
    const double& Fitness(void) const { return m_fit; }

    friend ostream& operator<< (ostream& ost, const Cell& cell)
    {
      ost << cell.m_gen << " / f = " << scientific << cell.m_fit << flush;
      return ost;
    }
  };

};

//#ifndef EVALUATE_H
//#include "Evaluate.h"
//#endif

#endif
