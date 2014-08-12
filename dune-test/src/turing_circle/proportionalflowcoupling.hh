#ifndef PROPORTIONALFLOWCOUPLING_HH
#define PROPORTIONALFLOWCOUPLING_HH

#include <dune/geometry/quadraturerules.hh>

#include <dune/pdelab/multidomain/couplingutilities.hh>
#include <dune/pdelab/localoperator/pattern.hh>
#include <dune/pdelab/localoperator/idefault.hh>

class ProportionalFlowCoupling :
  public Dune::PDELab::MultiDomain::CouplingOperatorDefaultFlags,
  public Dune::PDELab::MultiDomain::NumericalJacobianCoupling<ProportionalFlowCoupling>,
  public Dune::PDELab::MultiDomain::NumericalJacobianApplyCoupling<ProportionalFlowCoupling>,
  public Dune::PDELab::MultiDomain::FullCouplingPattern,
  public Dune::PDELab::InstationaryLocalOperatorDefaultMethods<double>
{

public:

  ProportionalFlowCoupling(double intensity)
    : _intensity(intensity)
  {}

  static const bool doAlphaCoupling = true;
  static const bool doPatternCoupling = true;

  template<typename IG, typename LFSU1, typename LFSU2, typename X, typename LFSV1, typename LFSV2,
           typename R>
  void alpha_coupling
  ( const IG& ig,
    const LFSU1& lfsu_s, const X& x_s, const LFSV1& lfsv_s,
    const LFSU2& lfsu_n, const X& x_n, const LFSV2& lfsv_n,
    R& r_s, R& r_n) const
  {
    // domain and range field type
    typedef typename LFSU1::Traits::FiniteElementType::
      Traits::LocalBasisType::Traits::DomainFieldType DF;
    typedef typename LFSU1::Traits::FiniteElementType::
      Traits::LocalBasisType::Traits::RangeFieldType RF;
    typedef typename LFSU1::Traits::FiniteElementType::
      Traits::LocalBasisType::Traits::RangeType RangeType;

    const int intorder = 4;

    typedef typename LFSU1::Traits::SizeType size_type;

    // dimensions
    const int dim = IG::dimension;

    // select quadrature rule
    Dune::GeometryType gtface = ig.geometryInInside().type();
    const Dune::QuadratureRule<DF,dim-1>& rule = Dune::QuadratureRules<DF,dim-1>::rule(gtface,intorder);

    // loop over quadrature points and integrate normal flux
    for (typename Dune::QuadratureRule<DF,dim-1>::const_iterator it=rule.begin(); it!=rule.end(); ++it)
      {
        // position of quadrature point in local coordinates of element
        Dune::FieldVector<DF,dim> local1 = ig.geometryInInside().global(it->position());
        Dune::FieldVector<DF,dim> local2 = ig.geometryInOutside().global(it->position());

        // evaluate ansatz shape functions (assume Galerkin for now)
        std::vector<RangeType> phi1(lfsv_s.size());
        lfsv_s.finiteElement().localBasis().evaluateFunction(local1,phi1);

        std::vector<RangeType> phi2(lfsv_n.size());
        lfsv_n.finiteElement().localBasis().evaluateFunction(local2,phi2);

        RF u_s(0.0);
        for (size_t i=0; i<lfsu_s.size(); i++)
          u_s += x_s(lfsu_s,i) * phi1[i];

        RF u_n(0.0);
        for (size_t i=0; i<lfsu_n.size(); i++)
          u_n += x_n(lfsu_n,i) * phi2[i];

        RF u_diff = _intensity*(u_s - u_n);

        // integrate J
        RF factor = it->weight()*ig.geometry().integrationElement(it->position());
        for (size_type i=0; i<lfsv_s.size(); i++)
          r_s.accumulate(lfsv_s,i,u_diff*phi1[i]*factor);
        for (size_type i=0; i<lfsv_n.size(); i++)
          r_n.accumulate(lfsv_n,i,-u_diff*phi2[i]*factor);

      }
  }

private:
  const double _intensity;

};


template<typename TReal>
class ContinuousValueContinuousFlowCoupling
  : public Dune::PDELab::MultiDomain::CouplingOperatorDefaultFlags
  , public Dune::PDELab::MultiDomain::NumericalJacobianCoupling<ContinuousValueContinuousFlowCoupling<TReal> >
  , public Dune::PDELab::MultiDomain::NumericalJacobianApplyCoupling<ContinuousValueContinuousFlowCoupling<TReal> >
  , public Dune::PDELab::MultiDomain::FullCouplingPattern
  , public Dune::PDELab::InstationaryLocalOperatorDefaultMethods<TReal>
{

public:

  ContinuousValueContinuousFlowCoupling(int intorder, double intensity = 1.0)
    : _intorder(intorder)
    , _intensity(intensity)
  {}

  static const bool doAlphaCoupling = true;
  static const bool doPatternCoupling = true;

  template<typename IG, typename LFSU1, typename LFSU2, typename X, typename LFSV1, typename LFSV2,
           typename R>
  void alpha_coupling
  ( const IG& ig,
    const LFSU1& lfsu1, const X& x1, const LFSV1& lfsv1,
    const LFSU2& lfsu2, const X& x2, const LFSV2& lfsv2,
    R& r1, R& r2) const
  {
    // domain and range field type
    typedef typename LFSU1::Traits::FiniteElementType::
      Traits::LocalBasisType::Traits::DomainFieldType DF;
    typedef typename LFSU1::Traits::FiniteElementType::
      Traits::LocalBasisType::Traits::RangeFieldType RF;
    typedef typename LFSU1::Traits::FiniteElementType::
      Traits::LocalBasisType::Traits::RangeType RangeType;
    typedef typename LFSU1::Traits::FiniteElementType::
      Traits::LocalBasisType::Traits::JacobianType JacobianType;

    typedef typename LFSU1::Traits::SizeType size_type;
    typedef typename IG::Element Element;
    typedef typename Element::Geometry::Jacobian GeometryJacobian;

    const double h_F = (ig.geometry().corner(0) - ig.geometry().corner(1)).two_norm();

    // dimensions
    const int dim = IG::dimension;

    // select quadrature rule
    Dune::GeometryType gtface = ig.geometryInInside().type();
    const Dune::QuadratureRule<DF,dim-1>& rule = Dune::QuadratureRules<DF,dim-1>::rule(gtface,_intorder);

    // save entity pointers to adjacent elements
    Element e1 = ig.insideElement();
    Element e2 = ig.outsideElement();

    // loop over quadrature points and integrate normal flux
    for (typename Dune::QuadratureRule<DF,dim-1>::const_iterator it=rule.begin(); it!=rule.end(); ++it)
      {
        // position of quadrature point in local coordinates of element
        Dune::FieldVector<DF,dim> local1 = ig.geometryInInside().global(it->position());
        Dune::FieldVector<DF,dim> local2 = ig.geometryInOutside().global(it->position());

        // evaluate ansatz shape functions (assume Galerkin for now)
        std::vector<RangeType> phi1(lfsv1.size());
        lfsv1.finiteElement().localBasis().evaluateFunction(local1,phi1);

        std::vector<RangeType> phi2(lfsv2.size());
        lfsv2.finiteElement().localBasis().evaluateFunction(local2,phi2);

        // evaluate gradient of shape functions
        std::vector<JacobianType> js1(lfsu1.size());
        lfsu1.finiteElement().localBasis().evaluateJacobian(local1,js1);
        std::vector<JacobianType> js2(lfsu2.size());
        lfsu2.finiteElement().localBasis().evaluateJacobian(local2,js2);

        // transform gradient to real element
        const GeometryJacobian& jac1 = e1.geometry().jacobianInverseTransposed(local1);
        const GeometryJacobian& jac2 = e2.geometry().jacobianInverseTransposed(local2);

        std::vector<Dune::FieldVector<RF,dim> > gradphi1(lfsu1.size());
        for (size_t i=0; i<lfsu1.size(); i++)
          {
            gradphi1[i] = 0.0;
            jac1.umv(js1[i][0],gradphi1[i]);
          }
        std::vector<Dune::FieldVector<RF,dim> > gradphi2(lfsu2.size());
        for (size_t i=0; i<lfsu2.size(); i++)
          {
            gradphi2[i] = 0.0;
            jac2.umv(js2[i][0],gradphi2[i]);
          }

        // compute gradient of u1
        Dune::FieldVector<RF,dim> gradu1(0.0);
        for (size_t i=0; i<lfsu1.size(); i++)
          gradu1.axpy(x1(lfsu1,i),gradphi1[i]);
        // compute gradient of u2
        Dune::FieldVector<RF,dim> gradu2(0.0);
        for (size_t i=0; i<lfsu2.size(); i++)
          gradu2.axpy(x2(lfsu2,i),gradphi2[i]);


        RF u1(0.0);
        for (size_t i=0; i<lfsu1.size(); i++)
          u1 += x1(lfsu1,i) * phi1[i];

        RF u2(0.0);
        for (size_t i=0; i<lfsu2.size(); i++)
          u2 += x2(lfsu2,i) * phi2[i];

        Dune::FieldVector<DF,dim> normal1 = ig.unitOuterNormal(it->position());
        Dune::FieldVector<DF,dim> normal2 = ig.unitOuterNormal(it->position());
        normal2 *= -1;
        //const RF mean_u = 0.5*(u1 + u2);
        const RF jump_u1 = u1 - u2;
        const RF jump_u2 = u2 - u1;
        Dune::FieldVector<RF,dim> mean_gradu = gradu1 + gradu2;
        mean_gradu *= 0.5;
        const double theta = 1;
        const double alpha = _intensity;

        // integrate
        const RF factor = it->weight()*ig.geometry().integrationElement(it->position());
        for (size_type i=0; i<lfsv1.size(); i++)
          {
            RF value = 0.0;
            value -= /* epsilon * */ (mean_gradu * normal1) * phi1[i];
            value -= theta * /* epsilon * */ 0.5 * (gradphi1[i] * normal1) * jump_u1;
            value += alpha / h_F * jump_u1 * phi1[i];
            r1.accumulate(lfsv1,i,factor * value);
          }
        for (size_type i=0; i<lfsv2.size(); i++)
          {
            RF value = 0.0;
            value -= /* epsilon * */ (mean_gradu * normal2) * phi2[i];
            value -= theta * /* epsilon * */ 0.5 * (gradphi2[i] * normal2) * jump_u2;
            value += alpha / h_F * jump_u2 * phi2[i];
            r2.accumulate(lfsv2,i,factor * value);
          }
      }
  }

private:
  const int _intorder;
  const double _intensity;

};


#endif // PROPORTIONALFLOWCOUPLING_HH
