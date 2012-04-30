#ifndef IRREP_SPACE_HEADER
#define IRREP_SPACE_HEADER
#include "Symmetry.h"

//In nonabelian symmetries each irrep could have many basis vectors, this class
// represents the irrep. Each basis vector in a irrep is represented by a 
// different class IrrepVector
namespace SpinAdapted{
class IrrepSpace
{
 private:
  friend class boost::serialization::access;
  template<class Archive> void serialize(Archive & ar, const unsigned int version)
  {
    ar & irrep;
  }
  int irrep;
 public:
  IrrepSpace() : irrep(0) {}
  explicit IrrepSpace(int ir) : irrep(ir) {}

  std::vector<IrrepSpace> operator+=(IrrepSpace rhs);
  bool operator==(IrrepSpace rhs) const;
  bool operator!=(IrrepSpace rhs) const;
  bool operator<(IrrepSpace rhs) const;
  friend std::vector<IrrepSpace> operator+(IrrepSpace lhs, IrrepSpace rhs);
  void Save(std::ofstream &ofs);
  void Load(std::ifstream &ifs);
  friend ostream& operator<<(ostream& os, const IrrepSpace s);
  int getirrep() const ;
};
}
#endif
