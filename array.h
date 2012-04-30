#include <vector>

namespace SpinAdapted{

  template <class T> class Vector : public vector<T>
  {
   public:
    Vector() : vector<T>(){}
    Vector(const Vector<T>& x) : vector<T>(x->vector<T>) {}
    
      
