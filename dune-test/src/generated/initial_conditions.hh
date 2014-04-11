#ifndef INITIAL_CONDITIONS_H
#define INITIAL_CONDITIONS_H

#include <dune/copasi/utilities/datahelper.hh>


/** \brief A function for initial values of u_0
 */
template<typename GV, typename RF>
class U0Initial
        : public Dune::PDELab::GridFunctionBase<Dune::PDELab::
        GridFunctionTraits<GV,RF,1,Dune::FieldVector<RF,1> >, U0Initial<GV,RF> >
{
    const GV& gv;
    const DataHelper data;

public:
    typedef Dune::PDELab::GridFunctionTraits<GV,RF,1,Dune::FieldVector<RF,1> > Traits;

    //! construct from grid view
    U0Initial (const GV& gv_)
        : gv(gv_)
        , data()
    {}
    U0Initial (const GV& gv_, const std::string& fileName = "")
        : gv(gv_)
        , data(fileName)
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
        if (data.isValid())
        {
            y = data((double)x[0], (double)x[1]);
            return;
        }

        auto r = 0.5*sqrt((x[0]-45.)*(x[0]-45.)+(x[1]-50.)*(x[1]-50.)+0.6*(x[0]-45.)*(x[1]-50.));
        y = 1.5*1./(1.+0.1*r*r);

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
    const DataHelper data;

public:
    typedef Dune::PDELab::GridFunctionTraits<GV,RF,1,Dune::FieldVector<RF,1> > Traits;

    //! construct from grid view
    U1Initial (const GV& gv_)
        : gv(gv_)
        , data()
    {}
    U1Initial (const GV& gv_, const std::string& fileName = "")
        : gv(gv_)
        , data(fileName)
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
        if (data.isValid())
        {
            y = data((double)x[0], (double)x[1]);
            return;
        }

        auto r = 0.5*sqrt((x[0]-45.)*(x[0]-45.)+(x[1]-50.)*(x[1]-50.)+0.6*(x[0]-45.)*(x[1]-50.));
        y = 2.*(1.-1./(1.+0.05*r*r));

    }

    //! get a reference to the grid view
    inline const GV& getGridView () {return gv;}
};


#endif // INITIAL_CONDITIONS_H
