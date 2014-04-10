#ifndef DUNE_TEST_LOCAL_OPERATOR_HH
#define DUNE_TEST_LOCAL_OPERATOR_HH

#include<dune/geometry/referenceelements.hh>
#include <dune/geometry/quadraturerules.hh>
#include <dune/pdelab/common/geometrywrapper.hh>
#include <dune/pdelab/localoperator/pattern.hh>
#include <dune/pdelab/localoperator/flags.hh>
#include <dune/pdelab/localoperator/idefault.hh>
#include <dune/pdelab/localoperator/defaultimp.hh>

namespace Dune {
  namespace PDELab {

    // reaction base class
    // do not compute anything
    class ReactionBaseAdapter
    {
    public:
      //! constructor stores reference to the model
      ReactionBaseAdapter()
      {}

      //! do one step
      template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
      void evaluate (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r)
      {

      }

      template<typename RF>
      void preStep(RF time,RF dt,int stages)
      {

      }

    };

    //! traits class for two phase parameter class
    template<typename GV, typename RF>
    struct DiffusionParameterTraits
    {
      //! \brief the grid view
      typedef GV GridViewType;

      //! \brief Enum for domain dimension
      enum {
        //! \brief dimension of the domain
        dimDomain = GV::dimension
      };

      //! \brief Export type for domain field
      typedef typename GV::Grid::ctype DomainFieldType;

      //! \brief domain type
      typedef Dune::FieldVector<DomainFieldType,dimDomain> DomainType;

      //! \brief domain type
      typedef Dune::FieldVector<DomainFieldType,dimDomain-1> IntersectionDomainType;

      //! \brief Export type for range field
      typedef RF RangeFieldType;

      //! \brief range type
      typedef Dune::FieldVector<RF,GV::dimensionworld> RangeType;

      //! grid types
      typedef typename GV::Traits::template Codim<0>::Entity ElementType;
      typedef typename GV::Intersection IntersectionType;
    };


    //! base class for parameter class
    template<class T, class Imp>
    class DiffusionMulticomponentInterface
    {
    public:
      typedef T Traits;

      //! scalar diffusion coefficient
      typename Traits::RangeFieldType
      D (const typename Traits::ElementType& e, const typename Traits::DomainType& x, std::size_t i) const
      {
        return asImp().component(i).D(e,x);
      }

      //! source term
      typename Traits::RangeFieldType
      q (const typename Traits::ElementType& e, const typename Traits::DomainType& x, std::size_t i) const
      {
        return asImp().component(i).q(e,x);
      }

      int
      bc (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x, std::size_t i) const
      {
        return asImp().component(i).bc(is,x);
      }

      //! Dirichlet boundary condition on inflow
      typename Traits::RangeFieldType
      g (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x, std::size_t i) const
      {
        return asImp().component(i).g(is,x);
      }

      //! Neumann boundary condition
      typename Traits::RangeFieldType
      j  (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x, std::size_t i) const
      {
        return asImp().component(i).j(is,x);
      }

      template<std::size_t i>
      void setTime(typename Traits::RangeFieldType t)
      {
        asImp().setTime(t);
      }

    private:
      Imp& asImp () {return static_cast<Imp &> (*this);}
      const Imp& asImp () const {return static_cast<const Imp &>(*this);}
    };


    /** A local operator for solving the two-component reaction-diffusion
     *  system of Fitzhugh-Nagumo type. See also
     *  http://en.wikipedia.org/wiki/Reaction-diffusion
     *
     *   - d_0 \Delta u_0 - f_0 = 0   in \Omega
     *   - d_1 \Delta u_1 - f_1 = 0   in \Omega
     *
     * with natural boundary conditions
     *   \nabla u_0 \cdot n = 0   on \partial\Omega
     *   \nabla u_1 \cdot v = 0   on \partial\Omega
     *
     * with conforming finite elements on all types of grids in any dimension
     */
    template<typename TP, typename RA = ReactionBaseAdapter>
    class Example05LocalOperator :
      public Dune::PDELab::NumericalJacobianApplyVolume<Example05LocalOperator<TP,RA> >,
      public Dune::PDELab::NumericalJacobianVolume<Example05LocalOperator<TP,RA> >,
      public Dune::PDELab::FullVolumePattern,
      public Dune::PDELab::LocalOperatorDefaultFlags,
      public Dune::PDELab::InstationaryLocalOperatorDefaultMethods<typename TP::Traits::RangeFieldType>
    {
    public:
      // pattern assembly flags
      enum { doPatternVolume = true };

      // residual assembly flags
      enum { doAlphaVolume = true };

      // constructor stores parameters
      Example05LocalOperator (TP& tp_, unsigned int intorder_=2):
        tp(tp_),
        ra(raDefault()),
        intorder(intorder_) //new
      {}

      // constructor stores parameters
      Example05LocalOperator (TP& tp_, RA& ra_, unsigned int intorder_=2):
        tp(tp_),
        ra(ra_),
        intorder(intorder_) //new
      {}

      // volume integral depending on test and ansatz functions
      template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
      void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
      {
        // select the two components (assume Galerkin scheme U=V)
        typedef typename LFSU::template Child<0>::Type LFSU0;       // extract components
        const LFSU0& lfsu0 = lfsu.template child<0>();           // with template magic
        typedef typename LFSU::template Child<1>::Type LFSU1;
        const LFSU1& lfsu1 = lfsu.template child<1>();

        // domain and range field type (assume both components have same RF)
        typedef typename LFSU0::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::DomainFieldType DF;
        typedef typename LFSU0::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeFieldType RF;
        typedef typename LFSU0::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::JacobianType JacobianType;
        typedef typename LFSU0::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeType RangeType;
        typedef typename LFSU::Traits::SizeType size_type;

        // dimensions
        const int dim = EG::Geometry::dimension;
        const int dimw = EG::Geometry::dimensionworld;

        // cell geometry
        const Dune::FieldVector<DF,dim>&
          cell_center_local = Dune::ReferenceElements<DF,dim>::general(eg.geometry().type()).position(0,0);

        // select quadrature rule
        Dune::GeometryType gt = eg.geometry().type();
        const Dune::QuadratureRule<DF,dim>& rule = Dune::QuadratureRules<DF,dim>::rule(gt,intorder);

        // loop over quadrature points
        for (typename Dune::QuadratureRule<DF,dim>::const_iterator
               it=rule.begin(); it!=rule.end(); ++it)
          {
            // evaluate basis functions on reference element
            std::vector<RangeType> phi0(lfsu0.size());
            lfsu0.finiteElement().localBasis().evaluateFunction(it->position(),phi0);
            std::vector<RangeType> phi1(lfsu1.size());
            lfsu1.finiteElement().localBasis().evaluateFunction(it->position(),phi1);
            // compute u_0, u_1 at integration point
            RF u_0=0.0;
            for (size_type i=0; i<lfsu0.size(); i++)
              u_0 += x(lfsu0,i)*phi0[i];                // localIndex() maps dof within
            RF u_1=0.0;                                             // leaf space to all dofs
            for (size_type i=0; i<lfsu1.size(); i++)                // within given element
              u_1 += x(lfsu1,i)*phi1[i];


            // evaluate gradient of basis functions on reference element
            std::vector<JacobianType> js0(lfsu0.size());
            lfsu0.finiteElement().localBasis().evaluateJacobian(it->position(),js0);
            std::vector<JacobianType> js1(lfsu1.size());
            lfsu1.finiteElement().localBasis().evaluateJacobian(it->position(),js1);


            // transform gradients from reference element to real element
            const Dune::FieldMatrix<DF,dimw,dim>
              jac = eg.geometry().jacobianInverseTransposed(it->position());
            std::vector<Dune::FieldVector<RF,dim> > gradphi0(lfsu0.size());
            for (size_type i=0; i<lfsu0.size(); i++)
              jac.mv(js0[i][0],gradphi0[i]);
            std::vector<Dune::FieldVector<RF,dim> > gradphi1(lfsu1.size());
            for (size_type i=0; i<lfsu1.size(); i++)
              jac.mv(js1[i][0],gradphi1[i]);

            // compute gradient of u_0, u_1
            Dune::FieldVector<RF,dim> gradu0(0.0);
            for (size_type i=0; i<lfsu0.size(); i++)
              gradu0.axpy(x(lfsu0,i),gradphi0[i]);
            Dune::FieldVector<RF,dim> gradu1(0.0);
            for (size_type i=0; i<lfsu1.size(); i++)
              gradu1.axpy(x(lfsu1,i),gradphi1[i]);


            // integrate both components
            RF factor = it->weight()*eg.geometry().integrationElement(it->position());
            // diffusion coefficient is assumed to be constant on element
            RF d_0 = tp.D(eg.entity(),cell_center_local,0);
            // eq. 0: - d_0 \Delta u_0 - (-k1*u_0*u_1 +v2) = 0
            for (size_type i=0; i<lfsu0.size(); i++)
              r.accumulate(lfsu0,i,(d_0*(gradu0*gradphi0[i])
                                    -ra.evaluate(0,u_0,u_1)*phi0[i])*factor);

            // diffusion coefficient is assumed to be constant on element
            RF d_1 = tp.D(eg.entity(),cell_center_local,1);
            // eq. 1: - d_1 \Delta u_1 - (k1*u_0*u_1 +v1 - v*u_1*1/(km+u_1)) = 0
            for (size_type i=0; i<lfsu1.size(); i++)
              r.accumulate(lfsu1,i,(d_1*(gradu1*gradphi1[i])
                                    -ra.evaluate(1,u_0,u_1)*phi1[i])*factor);

          }
      }

    private:
      static RA & raDefault()
      {
        static RA ra;
        return ra;
      }

      TP& tp;
      RA& ra;
      unsigned int intorder;
    };


    /** a local operator for the mass operator (L_2 integral) in the system
     *
     * u_0
     * u_1
     */

    template<class TP>
    class Example05TimeLocalOperator
      :public Dune::PDELab::FullVolumePattern,
       public Dune::PDELab::LocalOperatorDefaultFlags,
       public Dune::PDELab::InstationaryLocalOperatorDefaultMethods<typename TP::Traits::RangeFieldType>
    {
    public:
      // pattern assembly flags
      enum { doPatternVolume = true };

      // residual assembly flags
      enum { doAlphaVolume = true };

      // constructor remembers parameters
      Example05TimeLocalOperator (TP& tp_, unsigned int intorder_=2)
        : tp(tp_), intorder(intorder_) {}

      // volume integral depending on test and ansatz functions
      template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
      void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
      {
        // select the two components (assume Galerkin scheme U=V)
        typedef typename LFSU::template Child<0>::Type LFSU0;       // extract components
        const LFSU0& lfsu0 = lfsu.template child<0>();           // with template magic
        typedef typename LFSU::template Child<1>::Type LFSU1;
        const LFSU1& lfsu1 = lfsu.template child<1>();

        // domain and range field type (assume both components have same RF)
        typedef typename LFSU0::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::DomainFieldType DF;
        typedef typename LFSU0::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeFieldType RF;
        typedef typename LFSU0::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeType RangeType;
        typedef typename LFSU::Traits::SizeType size_type;

        // dimensions
        const int dim = EG::Geometry::dimension;

        // select quadrature rule
        Dune::GeometryType gt = eg.geometry().type();
        const Dune::QuadratureRule<DF,dim>& rule = Dune::QuadratureRules<DF,dim>::rule(gt,intorder);

        // loop over quadrature points
        for (typename Dune::QuadratureRule<DF,dim>::const_iterator it=rule.begin(); it!=rule.end(); ++it)
          {
            // evaluate basis functions on reference element
            std::vector<RangeType> phi0(lfsu0.size());
            lfsu0.finiteElement().localBasis().evaluateFunction(it->position(),phi0);
            std::vector<RangeType> phi1(lfsu1.size());
            lfsu1.finiteElement().localBasis().evaluateFunction(it->position(),phi1);

            // compute u_0, u_1 at integration point
            RF u_0=0.0;
            for (size_type i=0; i<lfsu0.size(); i++) u_0 += x(lfsu0,i)*phi0[i];
            RF u_1=0.0;
            for (size_type i=0; i<lfsu1.size(); i++) u_1 += x(lfsu1,i)*phi1[i];

            // integration
            RF factor = it->weight() * eg.geometry().integrationElement(it->position());
            for (size_type i=0; i<lfsu0.size(); i++)
              r.accumulate(lfsu0,i,u_0*phi0[i]*factor);
            for (size_type i=0; i<lfsu1.size(); i++)
              r.accumulate(lfsu1,i,u_1*phi1[i]*factor);
          }
      }


      // jacobian of volume term
      template<typename EG, typename LFSU, typename X, typename LFSV, typename M>
      void jacobian_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv,
                            M& mat) const
      {
        // select the two components (assume Galerkin scheme U=V)
        typedef typename LFSU::template Child<0>::Type LFSU0;       // extract components
        const LFSU0& lfsu0 = lfsu.template child<0>();           // with template magic
        typedef typename LFSU::template Child<1>::Type LFSU1;
        const LFSU1& lfsu1 = lfsu.template child<1>();
        // domain and range field type
        typedef typename LFSU0::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::DomainFieldType DF;
        typedef typename LFSU0::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeFieldType RF;
        typedef typename LFSU0::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeType RangeType;
        typedef typename LFSU0::Traits::SizeType size_type;

        // dimensions
        const int dim = EG::Geometry::dimension;

        // select quadrature rule
        Dune::GeometryType gt = eg.geometry().type();
        const Dune::QuadratureRule<DF,dim>& rule = Dune::QuadratureRules<DF,dim>::rule(gt,intorder);

        // loop over quadrature points
        for (typename Dune::QuadratureRule<DF,dim>::const_iterator it=rule.begin(); it!=rule.end(); ++it)
          {


            // evaluate basis functions
            std::vector<RangeType> phi0(lfsu0.size());
            lfsu0.finiteElement().localBasis().evaluateFunction(it->position(),phi0);

            // evaluate basis functions
            std::vector<RangeType> phi1(lfsu1.size());
            lfsu1.finiteElement().localBasis().evaluateFunction(it->position(),phi1);

            RF factor = it->weight() * eg.geometry().integrationElement(it->position());
            for (size_type j=0; j<lfsu0.size(); j++)
              for (size_type i=0; i<lfsu0.size(); i++)
                mat.accumulate(lfsu0,i,lfsu0,j,phi0[i]*phi0[j]*factor);

            for (size_type j=0; j<lfsu1.size(); j++)
              for (size_type i=0; i<lfsu1.size(); i++)
                mat.accumulate(lfsu1,i,lfsu1,j,phi1[i]*phi1[j]*factor);


          }
      }

    private:
      TP& tp;
      unsigned int intorder;
    };
  }
}

#endif //DUNE_TEST_LOCAL_OPERATOR_HH
