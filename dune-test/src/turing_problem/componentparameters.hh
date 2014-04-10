// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set ts=4 sw=2 et sts=2:

#include<dune/copasi/utilities/componentparameters.hh>
#include"local_operator.hh"

//! Transport in water phase
template<typename GV, typename RF>
class DiffusionParameter :
        public Dune::PDELab::DiffusionMulticomponentInterface<Dune::PDELab::DiffusionParameterTraits<GV,RF>,
        DiffusionParameter<GV,RF> >
{
    enum {dim=GV::Grid::dimension};

public:
    typedef Dune::PDELab::DiffusionParameterTraits<GV,RF> Traits;

    DiffusionParameter(const Dune::ParameterTree & param, const std::string cname)
        : time(0.)
        , Dt(param.sub(cname).template get<RF>("D"))
    {}


    //! tensor permeability
    typename Traits::RangeFieldType
    D (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
    {
        return Dt;
    }

    //! source/reaction term

    typename Traits::RangeFieldType
    q (const typename Traits::ElementType& e, const typename Traits::DomainType& ) const
    {
        return 0.0;
    }

    //! boundary condition type function
    // 0 means Neumann
    // 1 means Dirichlet
    int
    bc (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x) const
    {
        return 0; // only neumann
    }

    //! Dirichlet boundary condition value
    typename Traits::RangeFieldType
    g (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x) const
    {
        return 0; // Dirichlet is zero
    }

    //! Neumann boundary condition
    typename Traits::RangeFieldType
    j (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x) const
    {
        return 0.0;
    }

    //! set time for subsequent evaluation
    void setTime (RF t)
    {
        time = t;
    }

    void setTimeTarget(RF time_, RF dt_)
    {
        tend = time_;
    }

    //! to be called once before each time step
    void preStep (RF time_, RF dt_, int stages)
    {

    }


private:
    RF time, tend, dt;
    const RF Dt;
};
