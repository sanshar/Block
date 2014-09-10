#ifndef NPDM_SPIN_TRANSFORMATION_H
#define NPDM_SPIN_TRANSFORMATION_H

#include <cmath>
#include <limits>
#include <iterator>
#include <boost/spirit/include/qi.hpp>
#ifndef BOOST_1_56_0
#include <boost/spirit/home/phoenix/container.hpp>
#include <boost/spirit/home/phoenix/operator.hpp>
#include <boost/spirit/home/phoenix/object/construct.hpp>
#include <boost/spirit/home/phoenix/function/function.hpp>
#include <boost/spirit/home/phoenix/bind/bind_function.hpp>
#else
#include <boost/phoenix/stl/container.hpp>
#include <boost/phoenix/operator.hpp>
#include <boost/phoenix/object/construct.hpp>
#include <boost/phoenix/function/function.hpp>
#include <boost/phoenix/bind/bind_function.hpp>
#endif
//FIXME
#include "npdm_tensor_operator.h"
//#include "tensor_operator.h"

namespace SpinAdapted {
namespace Npdm {
namespace expression {

//===========================================================================================================================================================

  void appendTensorOp(std::vector<TensorOp>& invec, int index, int sign) 
  {
    TensorOp top(index, sign);
    std::vector<TensorOp> outvec;
    if (invec.size() == 0) {invec.push_back(top); return;}

    if (sign < 0)
      invec.push_back(top);
    else
      invec.insert(invec.begin(),top);
  }

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

  void addTensorOptoVec(std::vector<TensorOp>& invec, int index, int sign) 
  {
    TensorOp top(index, sign);
    std::vector<TensorOp> outvec;
    if (invec.size() == 0) {invec.push_back(top); return;}
    for (int i =0; i<invec.size(); i++) {
      TensorOp& op = invec[i];
      TensorOp opcopy(op);
      for (int spin = abs(op.Spin-top.Spin); spin<= op.Spin+top.Spin; spin+=2) {
	opcopy.product(top, spin, 0);
	outvec.push_back(opcopy);
	opcopy = op;
      }
    }
    invec = outvec;
  }

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

  void addTensorOpVectoVec(std::vector<TensorOp>& invec1, std::vector<TensorOp>& invec2) 
  {
    std::vector<TensorOp> outvec;
    if (invec2.size() == 0) return;
    for (int i2 = 0; i2<invec2.size(); i2++) {
      TensorOp& top = invec2[i2];
      
      if (invec1.size() == 0) {outvec.push_back(top);}
      else {
	for (int i =0; i<invec1.size(); i++) {
	  TensorOp& op = invec1[i];
	  TensorOp opcopy(op);
	  for (int spin = abs(op.Spin-top.Spin); spin<= op.Spin+top.Spin; spin+=2) {
	    opcopy.product(top, spin, 0);
	    outvec.push_back(opcopy);
	    opcopy = op;
	  }
	}
      }
      
    }
    invec1 = outvec;
  }

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

  template <typename FPT, typename Iterator> struct grammar : boost::spirit::qi::grammar <Iterator, FPT(), boost::spirit::ascii::space_type>
  {

      struct optype_ : boost::spirit::qi::symbols<char, int>
      {
	optype_()
	  {
	    add
	      ("C" , 1)
	      ("D" , -1)
	      ;
	  }
      };
      
      optype_ optype;

    
    boost::spirit::qi::rule <Iterator, FPT(), boost::spirit::ascii::space_type> expression, term;
    
    grammar() : grammar::base_type(expression)
    {
      
      using boost::spirit::qi::_1;
      using boost::spirit::qi::_2;
      using boost::spirit::qi::int_;
      using boost::spirit::qi::char_;
      using boost::spirit::qi::_val;
      using boost::spirit::qi::eps;


      //expression =  *( (optype>>int_)[boost::phoenix::push_back(_val, boost::phoenix::construct<TensorOp>(_2, _1))] ); //get all tensorops
      //expression =  *( (optype>>int_)[boost::phoenix::bind(appendTensorOp, _val, _2, _1)] );  //get tensorops with C before D
      

      expression =  +( ('(' >> expression >> ')')[boost::phoenix::bind(addTensorOpVectoVec, _val, _1)] |  //single term in a bracket
		      term[boost::phoenix::bind(addTensorOpVectoVec, _val, _1)]  //single term product of operators		      
		      );

      term =  +( (optype>>int_)[boost::phoenix::bind(addTensorOptoVec, _val, _2, _1)] 
		      );

    }
  };
  
//-----------------------------------------------------------------------------------------------------------------------------------------------------------

  template <typename FPT, typename Iterator>
    bool parse(Iterator &iter,
	       Iterator end,
	       const grammar<FPT,Iterator> &g,
	       FPT &result)
  {
    return boost::spirit::qi::phrase_parse(
					   iter, end, g, boost::spirit::ascii::space, result);
  }
  

//===========================================================================================================================================================

}
}
}

#endif
