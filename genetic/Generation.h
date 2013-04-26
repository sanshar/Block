#ifndef GENERATION_H
#define GENERATION_H

#include <iostream>
#include <vector>
#include <algorithm>
#include "Cell.h"
#include "GAInput.h"
using namespace std;

namespace genetic
{
  extern GAInput gainput;

  // Generation
  class Generation
  {
  private:
    vector<Cell> m_cells;
    vector<double> m_prob;

    void ComputeProbability(void);

    void BoltzmannProb(void);
    void GaussProb(void);
    void UniformProb(void);

  public:
    Generation(void);
    Generation(const Generation& other);

    // God Makes Life
    void Initialize(void);

    const Cell& Select(void) const;
          Cell& Select(void);

    inline const Cell& Select(const int& i) const { return m_cells[i]; }
    inline       Cell& Select(const int& i)       { return m_cells[i]; }

    inline const Cell& Min(void) const { return *min_element(m_cells.begin(), m_cells.end()); }
    inline       Cell& Min(void)       { return *min_element(m_cells.begin(), m_cells.end()); }

    inline const Cell& Max(void) const { return *max_element(m_cells.begin(), m_cells.end()); }
    inline       Cell& Max(void)       { return *max_element(m_cells.begin(), m_cells.end()); }

    void Generate(const Generation& ancestor);

    void AddFiedler(std::vector<int> fiedlerorder);

    friend ostream& operator<< (ostream& ost, const Generation& g){
      using std::setw;
      using std::endl;
      for(int i = 0; i < g.m_cells.size(); ++i)
        ost << setw(4) << i << ": " << g.m_cells[i] << " ( " << g.m_prob[i] << " ) " << endl;
      return ost;
    }
  };
};

#endif
