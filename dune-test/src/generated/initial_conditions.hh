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

    //! construct from grid view
    U0Initial (const GV& gv_)
        : gv(gv_)
    {}

    //! evaluate extended function on element
    inline void evaluate (const typename Traits::ElementType& e,
                          const typename Traits::DomainType& xlocal,
                          typename Traits::RangeType& y) const
    {
        const int dim = Traits::GridViewType::Grid::dimension;
        typedef typename Traits::GridViewType::Grid::ctype ctype;
        Dune::FieldVector<ctype,dim> x = e.geometry().global(xlocal);

        //new initial condition

        if (x[0] > 30. -1.e-06 &&  x[0] < 50. +1e-06 &&
                x[1] > 45. -1.e-06  &&  x[1] < 55. +1.e-06)
            y = 20.;
        else
            y= 0.5;

        return;
    }

    //! get a reference to the grid view
    inline const GV& getGridView () {return gv;}
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

    //! construct from grid view
    U1Initial (const GV& gv_)
        : gv(gv_)
    {}

    //! evaluate extended function on element
    inline void evaluate (const typename Traits::ElementType& e,
                          const typename Traits::DomainType& xlocal,
                          typename Traits::RangeType& y) const
    {
        //const int dim = Traits::GridViewType::Grid::dimension;
        //typedef typename Traits::GridViewType::Grid::ctype ctype;
        //Dune::FieldVector<ctype,dim> x = e.geometry().global(xlocal);


        //new initial condition

        y= 0.1;



        return;
    }

    //! get a reference to the grid view
    inline const GV& getGridView () {return gv;}
};


#endif // INITIAL_CONDITIONS_H
