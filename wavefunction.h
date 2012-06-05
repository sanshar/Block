/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is free software: you can redistribute it and/or modify         
it under the terms of the GNU General Public License as published by         
the Free Software Foundation, either version 3 of the License, or            
(at your option) any later version.                                          
                                                                             
This program is distributed in the hope that it will be useful,              
but WITHOUT ANY WARRANTY; without even the implied warranty of               
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                
GNU General Public License for more details.                                 
                                                                             
You should have received a copy of the GNU General Public License            
along with this program.  If not, see <http://www.gnu.org/licenses/>.        
*/


#ifndef SPIN_WAVEFUNCTION_HEADER
#define SPIN_WAVEFUNCTION_HEADER
#include "BaseOperator.h"

namespace SpinAdapted{
class Wavefunction : public SpinAdapted::SparseMatrix
{
private:
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version)
  {
    ar & boost::serialization::base_object<SparseMatrix>(*this);
    ar & onedot;
  } 
  bool onedot;
public:
 Wavefunction() : onedot(false), SparseMatrix(){}
  Wavefunction(const Wavefunction& wf) : SparseMatrix(wf), onedot(wf.onedot) {}
  Wavefunction(const SpinQuantum dQ, const SpinBlock* b, const bool onedot_) : onedot(onedot_) { initialise(dQ, b, onedot_); }
  void initialise(const SpinQuantum dQ, const SpinBlock* b, const bool &onedot_);
  const bool &get_onedot() const {return onedot;}
  void set_onedot(bool p_onedot) {onedot = p_onedot;}

  boost::shared_ptr<SparseMatrix> getworkingrepresentation(const SpinBlock* block) {return boost::shared_ptr<SparseMatrix> (this);}
  void build(const SpinBlock& b) {}
  double redMatrixElement(Csf c1, vector<Csf>& ladder, const SpinBlock* b=0) {return 0.0;}
  void CollectFrom(const RowVector& C);
  void FlattenInto(Matrix& C);

  void LoadWavefunctionInfo (StateInfo &waveInfo, const std::vector<int>& sites, const int wave_num);
  void SaveWavefunctionInfo (const StateInfo &waveInfo, const std::vector<int>& sites, const int wave_num);

  void UnCollectQuantaAlongRows(const StateInfo& sRow, const StateInfo& sCol);
  void AllowQuantaFor(const StateInfo& sRow, const StateInfo& sCol, const SpinQuantum q);
  void CollectQuantaAlongRows(const StateInfo& sRow, const StateInfo& sCol);
  void CollectQuantaAlongColumns(const StateInfo& sRow, const StateInfo& sCol);
  void UnCollectQuantaAlongColumns(const StateInfo& sRow, const StateInfo& sCol);
  Wavefunction& operator+=(const Wavefunction& other);
};
}
#endif
