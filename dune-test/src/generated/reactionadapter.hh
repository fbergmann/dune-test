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
template<typename RF, int N>
class ReactionAdapter
{
public:

    //! constructor stores reference to the model
    ReactionAdapter(const Dune::ParameterTree & param_)
        : param(param_)
	   , J0_k_1(param.sub("Reaction").template get<RF>("J0_k_1"))
	   , J1_k_2(param.sub("Reaction").template get<RF>("J1_k_2"))
	   , J2_k_3(param.sub("Reaction").template get<RF>("J2_k_3"))
	   , J3_k_4(param.sub("Reaction").template get<RF>("J3_k_4"))
	   , A(param.sub("Reaction").template get<RF>("A"))
	   , B(param.sub("Reaction").template get<RF>("B"))
	   , D(param.sub("Reaction").template get<RF>("D"))
	   , E(param.sub("Reaction").template get<RF>("E"))
	   , cell(param.sub("Reaction").template get<RF>("cell"))
	   , EC(param.sub("Reaction").template get<RF>("EC"))
	   , Membrane0(param.sub("Reaction").template get<RF>("Membrane0"))

    {

    }

    void preStep(RF time,RF dt,int stages)
    {

    }

    //! evaluate model in entity eg with local values x
    template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
    void evaluate (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r)
    {
        RF r1 = 1.00000000000000000E+000 * J0_k_1 * A / cell + 1.00000000000000000E+000 * J1_k_2 * x(lfsu,0) * x(lfsu,0) * x(lfsu,1) / cell +  - 1.00000000000000000E+000 * J2_k_3 * x(lfsu,0) * B / cell +  - 1.00000000000000000E+000 * J3_k_4 * x(lfsu,0) / cell;
        RF r2 =  - 1.00000000000000000E+000 * J1_k_2 * x(lfsu,0) * x(lfsu,0) * x(lfsu,1) / cell + 1.00000000000000000E+000 * J2_k_3 * x(lfsu,0) * B / cell;

        r.accumulate(lfsv,0,-r1*eg.geometry().volume());
        r.accumulate(lfsv,1,-r2*eg.geometry().volume());


//        RF r1 = -k1*x(lfsu,0)*x(lfsu,1)+v2;
//        RF r2 = k1*x(lfsu,0)*x(lfsu,1)+v1-V*x(lfsu,1)/(Km+x(lfsu,1));
//
//        r.accumulate(lfsv,0,-r1*eg.geometry().volume());
//        r.accumulate(lfsv,1,-r2*eg.geometry().volume());
    }

private:
    const Dune::ParameterTree & param;
    const RF J0_k_1;
    const RF J1_k_2;
    const RF J2_k_3;
    const RF J3_k_4;
    const RF A;
    const RF B;
    const RF D;
    const RF E;
    const RF cell;
    const RF EC;
    const RF Membrane0;

};

#endif // REACTIONADAPTER_H
