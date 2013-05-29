#include <iostream>
#include <cmath>
#include "Generation.h"
#include "genetic_utils.h"
using namespace std;

   void genetic::Generation::AddFiedler(std::vector<int> fiedlerorder){
      m_cells[0].Create(fiedlerorder);
      ComputeProbability();
   }

void genetic::Generation::ComputeProbability(void)
{
  // Sort by Fitness
  sort(m_cells.begin(), m_cells.end());
       if(gainput.select == genetic::GAUSS)    GaussProb();
  else if(gainput.select == genetic::BOLTZMANN) BoltzmannProb();
  else                                         UniformProb();
}

void genetic::Generation::GaussProb(void)
{
  int nCells = m_cells.size();
  m_prob = vector< double >(nCells, 0.0);
  for(int i = 0; i < nCells; ++i) m_prob[i] = sqrt(m_cells[i].Fitness());
  double minFit = m_prob[0];
  for(int i = 0; i < nCells; ++i)
  {
    double value = m_prob[i] - minFit;
    m_prob[i] = exp(-gainput.scaling * value * value);
  }
}

void genetic::Generation::BoltzmannProb(void)
{
  int nCells = m_cells.size();
  m_prob = vector< double >(nCells, 0.0);
  for(int i = 0; i < nCells; ++i) m_prob[i] = sqrt(m_cells[i].Fitness());
  double minFit = m_prob[0];
  for(int i = 0; i < nCells; ++i) m_prob[i] = exp(-gainput.scaling * (m_prob[i] - minFit));
}

void genetic::Generation::UniformProb(void)
{
  int nCells = m_cells.size();
  m_prob = vector< double >(nCells, 0.0);
  for(int i = 0; i < nCells; ++i) m_prob[i] = 1.0;
}

genetic::Generation::Generation() { Initialize(); }

genetic::Generation::Generation(const genetic::Generation& other) { m_cells = other.m_cells; }

void genetic::Generation::Initialize(void)
{
  int& nCells = gainput.max_cells;
  m_cells = vector<Cell>(nCells, Cell());
  for(int i = 0; i < nCells; ++i) {
     m_cells[i].Create(Gene::Random());
  }
  ComputeProbability();
}

const genetic::Cell& genetic::Generation::Select(void) const
{
  int& nCells = gainput.max_cells;
  int iCell = 0;
  while(1)
  {
    iCell = irand(nCells);
    if(m_prob[iCell] > drand(1.0)) break;
  }
  return m_cells[iCell];
}

genetic::Cell& genetic::Generation::Select(void)
{
  int& nCells = gainput.max_cells;
  int iCell = 0;
  while(1)
  {
    iCell = irand(nCells);
    if(m_prob[iCell] > drand(1.0)) break;
  }
  return m_cells[iCell];
}

void genetic::Generation::Generate(const genetic::Generation& ancestor)
{
  int& nCells = gainput.max_cells;
  m_cells = vector<Cell>(nCells, Cell());

  int iCell = 0;

  // Chose Elites
  for(; iCell < gainput.max_elite; ++iCell) m_cells[iCell] = ancestor.m_cells[iCell];

  // Selection to Next Generation
  while(iCell < nCells)
  {
    Gene child;

    if(drand(1.0) > gainput.thre_cloning)
    {
      child = ancestor.Select().Gen();
    }
    else
    {
      Cell father(ancestor.Select());
      Cell mother(ancestor.Select());
//    if(father.Fitness() > mother.Fitness())
        child = father.Gen() * mother.Gen();
//    else
//      child = mother.Gen() * father.Gen();
    }
    // Mutation
    if(drand(1.0) < gainput.thre_mutation) child.pMutate();
    if(drand(1.0) < gainput.thre_mutation) child.gMutate();

    m_cells[iCell++] = Cell(child);
  }

  ComputeProbability();
}
