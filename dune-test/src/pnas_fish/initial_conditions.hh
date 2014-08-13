// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set ts=4 sw=2 et sts=2:
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

  {}

  //! evaluate extended function on element
  inline void evaluate (const typename Traits::ElementType& e,
                        const typename Traits::DomainType& xlocal,
                        typename Traits::RangeType& __initial) const
  {

    const int dim = Traits::GridViewType::Grid::dimension;
    typedef typename Traits::GridViewType::Grid::ctype ctype;
    Dune::FieldVector<ctype,dim> x = e.geometry().global(xlocal);

    __initial = (((x[1] > 1.41000000000000000E+002 && x[1] < 1.53000000000000000E+002)) ? 4.00000000000000000E+000 : 0.00000000000000000E+000);


  }

  //! get a reference to the grid view
  inline const GV& getGridView () {return gv;}
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

  {}

  //! evaluate extended function on element
  inline void evaluate (const typename Traits::ElementType& e,
                        const typename Traits::DomainType& xlocal,
                        typename Traits::RangeType& __initial) const
  {

    __initial = 0.00000000000000000E+000;


  }

  //! get a reference to the grid view
  inline const GV& getGridView () {return gv;}
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



/** \brief A function for initial values of u_2
 */
template<typename GV, typename RF>
class U2Initial
  : public Dune::PDELab::GridFunctionBase<Dune::PDELab::
                                          GridFunctionTraits<GV,RF,1,Dune::FieldVector<RF,1> >, U2Initial<GV,RF> >
{
  const GV& gv;

public:
  typedef Dune::PDELab::GridFunctionTraits<GV,RF,1,Dune::FieldVector<RF,1> > Traits;

  U2Initial (const GV& gv_)
    : gv(gv_)
    , param()
  {}

  //! construct from grid view
  U2Initial (const GV& gv_, const Dune::ParameterTree & param_)
    : gv(gv_)
    , param(param_)
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

  {}

  //! evaluate extended function on element
  inline void evaluate (const typename Traits::ElementType& e,
                        const typename Traits::DomainType& xlocal,
                        typename Traits::RangeType& __initial) const
  {

    const int dim = Traits::GridViewType::Grid::dimension;
    typedef typename Traits::GridViewType::Grid::ctype ctype;
    Dune::FieldVector<ctype,dim> x = e.geometry().global(xlocal);

    __initial = (((x[1] > 1.41000000000000000E+002 && x[1] < 1.53000000000000000E+002)) ? 4.00000000000000000E+000 : 0.00000000000000000E+000);


  }

  //! get a reference to the grid view
  inline const GV& getGridView () {return gv;}
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



/** \brief A function for initial values of u_3
 */
template<typename GV, typename RF>
class U3Initial
  : public Dune::PDELab::GridFunctionBase<Dune::PDELab::
                                          GridFunctionTraits<GV,RF,1,Dune::FieldVector<RF,1> >, U3Initial<GV,RF> >
{
  const GV& gv;

public:
  typedef Dune::PDELab::GridFunctionTraits<GV,RF,1,Dune::FieldVector<RF,1> > Traits;

  U3Initial (const GV& gv_)
    : gv(gv_)
    , param()
  {}

  //! construct from grid view
  U3Initial (const GV& gv_, const Dune::ParameterTree & param_)
    : gv(gv_)
    , param(param_)
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

  {}

  //! evaluate extended function on element
  inline void evaluate (const typename Traits::ElementType& e,
                        const typename Traits::DomainType& xlocal,
                        typename Traits::RangeType& __initial) const
  {

    const int dim = Traits::GridViewType::Grid::dimension;
    typedef typename Traits::GridViewType::Grid::ctype ctype;
    Dune::FieldVector<ctype,dim> x = e.geometry().global(xlocal);

    __initial = (((x[1] > 1.41000000000000000E+002 && x[1] < 1.53000000000000000E+002)) ? 4.00000000000000000E+000 : 0.00000000000000000E+000);


  }

  //! get a reference to the grid view
  inline const GV& getGridView () {return gv;}
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


#endif // INITIAL_CONDITIONS_H
