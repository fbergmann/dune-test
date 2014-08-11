#ifndef REACTIONADAPTER_H
#define REACTIONADAPTER_H

#include <dune/copasi/utilities/sbmlhelper.hh>

//! Adapter for system of equations
/*!
  \tparam M model type
  This class is used for adapting system of ODE's (right side)
  to system of PDE's (as source term)
  What need to be implemented is only the evaluate function
*/
template<typename RF>
class ReactionAdapter
{
public:

    //! constructor stores reference to the model
    ReactionAdapter(const Dune::ParameterTree & param_)
        : param(param_)
	   , v1(param.sub("Reaction").template get<RF>("v1"))
	   , v2(param.sub("Reaction").template get<RF>("v2"))
	   , V(param.sub("Reaction").template get<RF>("V"))
	   , Km(param.sub("Reaction").template get<RF>("Km"))
	   , k1(param.sub("Reaction").template get<RF>("k1"))
	   , compartment_1(param.sub("Reaction").template get<RF>("compartment_1"))
	   , c2(param.sub("Reaction").template get<RF>("c2"))
	   , Membrane2(param.sub("Reaction").template get<RF>("Membrane2"))

    {

    }

    void preStep(RF time,RF dt,int stages)
    {

    }

    //! evaluate model in entity eg with local values x
    template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
    void evaluate (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r)
    {
        RF r1 = (1.00000000000000000E+000 * k1 * x(lfsu,1) * x(lfsu,0) / compartment_1 + 1.00000000000000000E+000 * v1 / compartment_1) + ( - 1.00000000000000000E+000) * V * x(lfsu,0) * 1.00000000000000000E+000 / (Km + x(lfsu,0)) / compartment_1;
        RF r2 = ( - 1.00000000000000000E+000) * k1 * x(lfsu,1) * x(lfsu,0) / compartment_1 + 1.00000000000000000E+000 * v2 / compartment_1;

        r.accumulate(lfsv,0,-r1*eg.geometry().volume());
        r.accumulate(lfsv,1,-r2*eg.geometry().volume());


    }

private:
    const Dune::ParameterTree & param;
    const RF v1;
    const RF v2;
    const RF V;
    const RF Km;
    const RF k1;
    const RF compartment_1;
    const RF c2;
    const RF Membrane2;

};

#endif // REACTIONADAPTER_H
