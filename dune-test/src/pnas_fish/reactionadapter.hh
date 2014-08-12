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
	   , KMOLE(param.sub("Reaction").template get<RF>("KMOLE"))
	   , re0_cu(param.sub("Reaction").template get<RF>("re0_cu"))
	   , re0_c3(param.sub("Reaction").template get<RF>("re0_c3"))
	   , re0_c2(param.sub("Reaction").template get<RF>("re0_c2"))
	   , re0_c1(param.sub("Reaction").template get<RF>("re0_c1"))
	   , re0_U(param.sub("Reaction").template get<RF>("re0_U"))
	   , re1_cv(param.sub("Reaction").template get<RF>("re1_cv"))
	   , re1_c6(param.sub("Reaction").template get<RF>("re1_c6"))
	   , re1_c5(param.sub("Reaction").template get<RF>("re1_c5"))
	   , re1_c4(param.sub("Reaction").template get<RF>("re1_c4"))
	   , re1_V(param.sub("Reaction").template get<RF>("re1_V"))
	   , re2_cw(param.sub("Reaction").template get<RF>("re2_cw"))
	   , re2_c9(param.sub("Reaction").template get<RF>("re2_c9"))
	   , re2_c8(param.sub("Reaction").template get<RF>("re2_c8"))
	   , re2_c7(param.sub("Reaction").template get<RF>("re2_c7"))
	   , re2_W(param.sub("Reaction").template get<RF>("re2_W"))
	   , outside(param.sub("Reaction").template get<RF>("outside"))
	   , eye(param.sub("Reaction").template get<RF>("eye"))
	   , fish(param.sub("Reaction").template get<RF>("fish"))
	   , Membrane0(param.sub("Reaction").template get<RF>("Membrane0"))
	   , Membrane1(param.sub("Reaction").template get<RF>("Membrane1"))

    {

    }

    void preStep(RF time,RF dt,int stages)
    {

    }

    //! evaluate model in entity eg with local values x
    template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
    void evaluate (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r)
    {
        RF r1 = 1.00000000000000000E+000 * ((((((re0_c1 * x(lfsu,2) + re0_c2 * x(lfsu,3)) + re0_c3) >= re0_U) ? re0_U : 0.00000000000000000E+000) + (((((re0_c1 * x(lfsu,2) + re0_c2 * x(lfsu,3)) + re0_c3) > 0.00000000000000000E+000 && ((re0_c1 * x(lfsu,2) + re0_c2 * x(lfsu,3)) + re0_c3) < re0_U)) ? (re0_c1 * x(lfsu,2) + re0_c2 * x(lfsu,3)) + re0_c3 : 0.00000000000000000E+000)) + ( - re0_cu * x(lfsu,0))) / fish;
        RF r2 = (( - 1.00000000000000000E+000) * ((((((re0_c1 * x(lfsu,2) + re0_c2 * x(lfsu,3)) + re0_c3) >= re0_U) ? re0_U : 0.00000000000000000E+000) + (((((re0_c1 * x(lfsu,2) + re0_c2 * x(lfsu,3)) + re0_c3) > 0.00000000000000000E+000 && ((re0_c1 * x(lfsu,2) + re0_c2 * x(lfsu,3)) + re0_c3) < re0_U)) ? (re0_c1 * x(lfsu,2) + re0_c2 * x(lfsu,3)) + re0_c3 : 0.00000000000000000E+000)) + ( - re0_cu * x(lfsu,0))) / fish + ( - 1.00000000000000000E+000) * ((((((re1_c4 * x(lfsu,0) + re1_c5 * x(lfsu,3)) + re1_c6) >= re1_V) ? re1_V : 0.00000000000000000E+000) + (((((re1_c4 * x(lfsu,0) + re1_c5 * x(lfsu,3)) + re1_c6) > 0.00000000000000000E+000 && ((re1_c4 * x(lfsu,0) + re1_c5 * x(lfsu,3)) + re1_c6) < re1_V)) ? (re1_c4 * x(lfsu,0) + re1_c5 * x(lfsu,3)) + re1_c6 : 0.00000000000000000E+000)) + ( - re1_cv * x(lfsu,2))) / fish) + ( - 1.00000000000000000E+000) * ((((((re2_c7 * x(lfsu,0) + re2_c8 * x(lfsu,2)) + re2_c9) >= re2_W) ? re2_W : 0.00000000000000000E+000) + (((((re2_c7 * x(lfsu,0) + re2_c8 * x(lfsu,2)) + re2_c9) > 0.00000000000000000E+000 && ((re2_c7 * x(lfsu,0) + re2_c8 * x(lfsu,2)) + re2_c9) < re2_W)) ? (re2_c7 * x(lfsu,0) + re2_c8 * x(lfsu,2)) + re2_c9 : 0.00000000000000000E+000)) + ( - re2_cw * x(lfsu,3))) / fish;
        RF r3 = 1.00000000000000000E+000 * ((((((re1_c4 * x(lfsu,0) + re1_c5 * x(lfsu,3)) + re1_c6) >= re1_V) ? re1_V : 0.00000000000000000E+000) + (((((re1_c4 * x(lfsu,0) + re1_c5 * x(lfsu,3)) + re1_c6) > 0.00000000000000000E+000 && ((re1_c4 * x(lfsu,0) + re1_c5 * x(lfsu,3)) + re1_c6) < re1_V)) ? (re1_c4 * x(lfsu,0) + re1_c5 * x(lfsu,3)) + re1_c6 : 0.00000000000000000E+000)) + ( - re1_cv * x(lfsu,2))) / fish;
        RF r4 = 1.00000000000000000E+000 * ((((((re2_c7 * x(lfsu,0) + re2_c8 * x(lfsu,2)) + re2_c9) >= re2_W) ? re2_W : 0.00000000000000000E+000) + (((((re2_c7 * x(lfsu,0) + re2_c8 * x(lfsu,2)) + re2_c9) > 0.00000000000000000E+000 && ((re2_c7 * x(lfsu,0) + re2_c8 * x(lfsu,2)) + re2_c9) < re2_W)) ? (re2_c7 * x(lfsu,0) + re2_c8 * x(lfsu,2)) + re2_c9 : 0.00000000000000000E+000)) + ( - re2_cw * x(lfsu,3))) / fish;

        r.accumulate(lfsv,0,-r1*eg.geometry().volume());
        r.accumulate(lfsv,1,-r2*eg.geometry().volume());
        r.accumulate(lfsv,2,-r3*eg.geometry().volume());
        r.accumulate(lfsv,3,-r4*eg.geometry().volume());


    }

private:
    const Dune::ParameterTree & param;
    const RF KMOLE;
    const RF re0_cu;
    const RF re0_c3;
    const RF re0_c2;
    const RF re0_c1;
    const RF re0_U;
    const RF re1_cv;
    const RF re1_c6;
    const RF re1_c5;
    const RF re1_c4;
    const RF re1_V;
    const RF re2_cw;
    const RF re2_c9;
    const RF re2_c8;
    const RF re2_c7;
    const RF re2_W;
    const RF outside;
    const RF eye;
    const RF fish;
    const RF Membrane0;
    const RF Membrane1;

};

#endif // REACTIONADAPTER_H
