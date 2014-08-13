// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set ts=4 sw=2 et sts=2:
#ifndef DUNE_COPASI_DIFFUSIONPARAMETERS_HH
#define DUNE_COPASI_DIFFUSIONPARAMETERS_HH

#include<dune/copasi/utilities/componentparameters.hh>
#include "local_operator.hh"

//! Transport in water phase
template<typename GV, typename RF>
class DiffusionParameter :
  public Dune::PDELab::DiffusionMulticomponentInterface<Dune::PDELab::DiffusionParameterTraits<GV,RF>,
                                                        DiffusionParameter<GV,RF> >
{
  enum {dim=GV::Grid::dimension};

public:
  typedef Dune::PDELab::DiffusionParameterTraits<GV,RF> Traits;
  typedef typename Traits::BCType BCType;

  DiffusionParameter(const Dune::ParameterTree & param, const std::string cname)
    : time(0.)
    , Dt(param.sub(cname).template get<RF>("D"))
    , Xmin(param.sub(cname).template get<RF>("Xmin", 0))
    , Xmax(param.sub(cname).template get<RF>("Xmax", 0))
    , Ymin(param.sub(cname).template get<RF>("Ymin", 0))
    , Ymax(param.sub(cname).template get<RF>("Ymax", 0))
    , width(param.sub("Domain").template get<int>("width", 0))
    , height(param.sub("Domain").template get<int>("height", 0))
    , boundarytype(BCType::Neumann)
  {
    //   param.sub(cname).template get<BCType>("BCType");
    int bc = param.sub(cname).template get<int>("BCType");
    switch(bc) {
    case 1:
      boundarytype=BCType::Dirichlet;
      break;
    case -1:
      boundarytype=BCType::Neumann;
      break;
    case -2:
      boundarytype=BCType::Outflow;
      break;
    case -3:
      boundarytype=BCType::None;
      break;
    }
  }


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
  BCType
  bctype (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x) const
  {
    return boundarytype;
  }

  //! Dirichlet boundary condition value
  typename Traits::RangeFieldType
  g (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x_) const
  {
    //  if (!is.boundary() || boundarytype == BCType::Dirichlet)
    //       return 0; // Dirichlet is zero
    typename Traits::DomainType x = is.geometry().global(x_);
    if (x[0] < 1e-6)
      return Xmin;
    if (x[0] > width - 1e-6)
      return Xmax;
    if (x[1] < 1e-6)
      return Ymin;
    if (x[1] > height - 1e-6)
      return Ymax;
    return 0.0;

  }

  //! Neumann boundary condition
  typename Traits::RangeFieldType
  j (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x_) const
  {
    //   if (!is.boundary() || boundarytype!=BCType::Neumann)
    //    return 0.0;
    typename Traits::DomainType x = is.geometry().global(x_);
    if (x[0] < 1e-6)
      return Xmin;
    if (x[0] > width - 1e-6)
      return Xmax;
    if (x[1] < 1e-6)
      return Ymin;
    if (x[1] > height - 1e-6)
      return Ymax;
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
  const RF Xmin;
  const RF Xmax;
  const RF Ymin;
  const RF Ymax;
  const RF width;
  const RF height;
  BCType boundarytype;
};


#endif
