#ifndef INITIAL_CONDITIONS_H
#define INITIAL_CONDITIONS_H

/** \brief A function for initial values of u_0
 */
template<typename GV, typename RF>
class U0Initial
        : public Dune::PDELab::GridFunctionBase<Dune::PDELab::
        GridFunctionTraits<GV,RF,1,Dune::FieldVector<RF,1> >, U0Initial<GV,RF> >
{
    const GV& gv;

public:
    typedef Dune::PDELab::GridFunctionTraits<GV,RF,1,Dune::FieldVector<RF,1> > Traits;

    U0Initial (const GV& gv_)
      : gv(gv_)
      , param()
    {}

    //! construct from grid view
    U0Initial (const GV& gv_, const Dune::ParameterTree & param_)
        : gv(gv_)
        , param(param_)
	   , v1(param.sub("Reaction").template get<RF>("v1"))
	   , v2(param.sub("Reaction").template get<RF>("v2"))
	   , V(param.sub("Reaction").template get<RF>("V"))
	   , Km(param.sub("Reaction").template get<RF>("Km"))
	   , k1(param.sub("Reaction").template get<RF>("k1"))
	   , compartment_1(param.sub("Reaction").template get<RF>("compartment_1"))
	   , c2(param.sub("Reaction").template get<RF>("c2"))
	   , Membrane2(param.sub("Reaction").template get<RF>("Membrane2"))

    {}

    //! evaluate extended function on element
    inline void evaluate (const typename Traits::ElementType& e,
                          const typename Traits::DomainType& xlocal,
                          typename Traits::RangeType& __initial) const
    {
    
        const int dim = Traits::GridViewType::Grid::dimension;
        typedef typename Traits::GridViewType::Grid::ctype ctype;
        Dune::FieldVector<ctype,dim> x = e.geometry().global(xlocal);

        __initial = 1.00000000000000010E-001;

        
    }

    //! get a reference to the grid view
    inline const GV& getGridView () {return gv;}
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



/** \brief A function for initial values of u_1
 */
template<typename GV, typename RF>
class U1Initial
        : public Dune::PDELab::GridFunctionBase<Dune::PDELab::
        GridFunctionTraits<GV,RF,1,Dune::FieldVector<RF,1> >, U1Initial<GV,RF> >
{
    const GV& gv;

public:
    typedef Dune::PDELab::GridFunctionTraits<GV,RF,1,Dune::FieldVector<RF,1> > Traits;

    U1Initial (const GV& gv_)
      : gv(gv_)
      , param()
    {}

    //! construct from grid view
    U1Initial (const GV& gv_, const Dune::ParameterTree & param_)
        : gv(gv_)
        , param(param_)
	   , v1(param.sub("Reaction").template get<RF>("v1"))
	   , v2(param.sub("Reaction").template get<RF>("v2"))
	   , V(param.sub("Reaction").template get<RF>("V"))
	   , Km(param.sub("Reaction").template get<RF>("Km"))
	   , k1(param.sub("Reaction").template get<RF>("k1"))
	   , compartment_1(param.sub("Reaction").template get<RF>("compartment_1"))
	   , c2(param.sub("Reaction").template get<RF>("c2"))
	   , Membrane2(param.sub("Reaction").template get<RF>("Membrane2"))

    {}

    //! evaluate extended function on element
    inline void evaluate (const typename Traits::ElementType& e,
                          const typename Traits::DomainType& xlocal,
                          typename Traits::RangeType& __initial) const
    {
    
        const int dim = Traits::GridViewType::Grid::dimension;
        typedef typename Traits::GridViewType::Grid::ctype ctype;
        Dune::FieldVector<ctype,dim> x = e.geometry().global(xlocal);

        __initial = ((((x[0] > 2.00000000000000000E+001 && x[0] < 2.50000000000000000E+001) && (x[1] > 2.00000000000000000E+001 && x[1] < 2.50000000000000000E+001))) ? 2.00000000000000000E+001 : 5.00000000000000000E-001);

        
    }

    //! get a reference to the grid view
    inline const GV& getGridView () {return gv;}
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


#endif // INITIAL_CONDITIONS_H
