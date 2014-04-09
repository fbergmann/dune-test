#ifndef REACTIONADAPTER_H
#define REACTIONADAPTER_H

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
       , v1(param.sub("Reaction").template get<RF>("v1"))
       , v2(param.sub("Reaction").template get<RF>("v2"))
       , V(param.sub("Reaction").template get<RF>("V"))
       , k1(param.sub("Reaction").template get<RF>("k1"))
       , Km(param.sub("Reaction").template get<RF>("Km")
       )
    {

    }

    void preStep(RF time,RF dt,int stages)
    {

    }

    //! evaluate model in entity eg with local values x
    template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
    void evaluate (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r)
    {
        RF r1 = -k1*x(lfsu,0)*x(lfsu,1)+v2;
        RF r2 = k1*x(lfsu,0)*x(lfsu,1)+v1-V*x(lfsu,1)/(Km+x(lfsu,1));

        r.accumulate(lfsv,0,-r1*eg.geometry().volume());
        r.accumulate(lfsv,1,-r2*eg.geometry().volume());
    }

private:
    const Dune::ParameterTree & param;
    const RF v1;
    const RF v2;
    const RF V;
    const RF k1;
    const RF Km;

};

#endif // REACTIONADAPTER_H
