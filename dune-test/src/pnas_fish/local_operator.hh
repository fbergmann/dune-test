// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set ts=4 sw=2 et sts=2:
#ifndef DUNE_PDELAB_MULTICOMPONENTTRANSPORTOP_HH
#define DUNE_PDELAB_MULTICOMPONENTTRANSPORTOP_HH

#include<dune/common/exceptions.hh>
#include<dune/common/fvector.hh>
#include<dune/common/static_assert.hh>

#include <dune/geometry/referenceelements.hh>

#include<dune/pdelab/common/geometrywrapper.hh>
#include<dune/pdelab/common/function.hh>
#include<dune/pdelab/localoperator/defaultimp.hh>
#include<dune/pdelab/localoperator/pattern.hh>
#include<dune/pdelab/localoperator/flags.hh>
#include<dune/pdelab/localoperator/idefault.hh>

#define ANALYTICAL_JACOBIAN FALSE

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




/** a local operator for a cell-centered finite folume scheme for
        the diffusion-reaction equation

        \nabla \cdot \{v u- D \nabla u \} = q in \Omega
        u = g on \Gamma_D
        \{v u - D \nabla u \} \cdot \nu = j on \Gamma_N
        outflow on \Gamma_O

        Modified version for the case

        d_t (c(x,t)u(x,t)) + \nabla \cdot \{v u - D \nabla u \} = q in \Omega

        where c(x,t) may become zero. We assume that the following holds:

        c(x,t+dt) <= eps  ==>  c(x,t) <= eps



        \tparam TP  parameter class implementing ComponentDiffusionParameterInterface
    */
template<typename TP, typename RA = ReactionBaseAdapter>
class MulticomponentCCFVSpatialDiffusionOperator :
        public NumericalJacobianVolume<MulticomponentCCFVSpatialDiffusionOperator<TP,RA> >,
        public NumericalJacobianApplyVolume<MulticomponentCCFVSpatialDiffusionOperator<TP,RA> >,
        public NumericalJacobianSkeleton<MulticomponentCCFVSpatialDiffusionOperator<TP,RA> >,
        public NumericalJacobianBoundary<MulticomponentCCFVSpatialDiffusionOperator<TP,RA> >,
        public NumericalJacobianApplySkeleton<MulticomponentCCFVSpatialDiffusionOperator<TP,RA> >,
        public NumericalJacobianApplyBoundary<MulticomponentCCFVSpatialDiffusionOperator<TP,RA> >,

        public FullSkeletonPattern,
        public FullVolumePattern,
        public LocalOperatorDefaultFlags,
        public InstationaryLocalOperatorDefaultMethods<typename TP::Traits::RangeFieldType>
{

public:
    // pattern assembly flags
    enum { doPatternVolume = true };
    enum { doPatternSkeleton = true };

    // residual assembly flags
    enum { doAlphaVolume  = true };
    enum { doAlphaSkeleton  = true };
    enum { doAlphaBoundary  = true };

    enum { dim = TP::Traits::GridViewType::dimension };

    MulticomponentCCFVSpatialDiffusionOperator (TP& tp_)
        : tp(tp_), ra(raDefault())
    {
    }

    MulticomponentCCFVSpatialDiffusionOperator (TP& tp_, RA& ra_)
        : tp(tp_), ra(ra_)
    {
    }


    // volume integral depending on test and ansatz functions
    template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
    void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
    {

        typedef typename LFSV::template Child<0>::Type Space;

        // domain and range field type
        typedef typename Space::Traits::FiniteElementType::
                Traits::LocalBasisType::Traits::DomainFieldType DF;

        // dimensions
        const int dim = EG::Geometry::dimension;

        // cell center
        const Dune::FieldVector<DF,dim>&
                inside_local = Dune::ReferenceElements<DF,dim>::general(eg.entity().type()).position(0,0);


        // here could described the source (right hand side of the sytem of equations)
        // for each equation (big ammount of code)
        for (std::size_t k = 0;k<TP::COMPONENTS;k++)
        {
            // evaluate source term
            typename TP::Traits::RangeFieldType q = tp.q(eg.entity(),inside_local,k);
            r.accumulate(lfsv,k,-q*eg.geometry().volume());
        }

        ra.evaluate(eg.entity(),lfsu,x,lfsv,r);

    }



    // skeleton integral depending on test and ansatz functions
    // We put the Dirchlet evaluation also in the alpha term to save some geometry evaluations
    template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
    void alpha_boundary (const IG& ig,
                         const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                         R& r_s) const
    {
        typedef typename LFSV::template Child<0>::Type Space;

        // domain and range field type
        typedef typename Space::Traits::FiniteElementType::
                Traits::LocalBasisType::Traits::DomainFieldType DF;
        typedef typename Space::Traits::FiniteElementType::
                Traits::LocalBasisType::Traits::RangeFieldType RF;

        // face geometry
        const Dune::FieldVector<DF,IG::dimension-1>&
                face_local = Dune::ReferenceElements<DF,IG::dimension-1>::general(ig.geometry().type()).position(0,0);
        RF face_volume = ig.geometry().volume();

        const Dune::FieldVector<DF,IG::dimension>&
                inside_local = Dune::ReferenceElements<DF,IG::dimension>::general(ig.inside()->type()).position(0,0);



        Dune::FieldVector<DF,IG::dimension>
                inside_global = ig.inside()->geometry().center();
        Dune::FieldVector<DF,IG::dimension>
                outside_global = ig.geometry().center();
        inside_global -= outside_global;
        RF distance = inside_global.two_norm();

        for (std::size_t k = 0;k<TP::COMPONENTS;k++)
        {
            typename TP::Traits::RangeFieldType D_inside = tp.D(*(ig.inside()),inside_local,k);
            // evaluate boundary condition type
            int bc = tp.bc(ig.intersection(),face_local,k);

            // do things depending on boundary condition type
            if (bc==0) // Neumann boundary
            {
                typename TP::Traits::RangeFieldType j = tp.j(ig.intersection(),face_local,k);
                r_s.accumulate(lfsu_s,k,j*face_volume);
            }


            if (bc==1) // Dirichlet boundary
            {
                typename TP::Traits::RangeFieldType g;
                g=tp.g(ig.intersection(),face_local,k);
                r_s.accumulate(lfsu_s,k,( - D_inside*(g-x_s(lfsu_s,k))/distance)*face_volume);
            }
        }

    }

#ifdef ANALYTICAL_JACOBIAN
    // jacobian of boundary term
    template<typename IG, typename LFSU, typename X, typename LFSV, typename M>
    void jacobian_boundary (const IG& ig,
                            const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                            M& mat_ss) const
    {
        typedef typename LFSV::template Child<0>::Type Space;

        // domain and range field type
        typedef typename Space::Traits::FiniteElementType::
                Traits::LocalBasisType::Traits::DomainFieldType DF;
        typedef typename Space::Traits::FiniteElementType::
                Traits::LocalBasisType::Traits::RangeFieldType RF;

        // face geometry
        const Dune::FieldVector<DF,IG::dimension-1>&
                face_local = Dune::ReferenceElements<DF,IG::dimension-1>::general(ig.geometry().type()).position(0,0);
        RF face_volume = ig.geometry().volume();

        const Dune::FieldVector<DF,IG::dimension>&
                inside_local = Dune::ReferenceElements<DF,IG::dimension>::general(ig.inside()->type()).position(0,0);
        Dune::FieldVector<DF,IG::dimension>
                inside_global = ig.inside()->geometry().center();
        Dune::FieldVector<DF,IG::dimension>
                outside_global = ig.geometry().center();
        inside_global -= outside_global;
        RF distance = inside_global.two_norm();

        for (std::size_t k = 0;k<TP::COMPONENTS;k++)
        {
            // evaluate boundary condition type
            int bc = tp.bc(ig.intersection(),face_local,k);

            // do things depending on boundary condition type
            if (bc==0) // Neumann boundary
            {
                return;
            }

            if (bc==1) // Dirichlet boundary
            {
                typename TP::Traits::RangeFieldType D_inside = tp.D(*(ig.inside()),inside_local,k);
                mat_ss.accumulate(lfsu_s,k,lfsu_s,k,D_inside/distance*face_volume);
                return;
            }
        }
    }
#endif //ANALYTICAL_JACOBIAN


    // skeleton integral depending on test and ansatz functions
    // We put the Dirchlet evaluation also in the alpha term to save some geometry evaluations
    template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
    void alpha_skeleton (const IG& ig,
                         const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                         const LFSU& lfsu_n, const X& x_n, const LFSV& lfsv_n,
                         R& r_s, R& r_n) const
    {

        typedef typename LFSV::template Child<0>::Type Space;

        // domain and range field type
        typedef typename Space::Traits::FiniteElementType::
                Traits::LocalBasisType::Traits::DomainFieldType DF;
        typedef typename Space::Traits::FiniteElementType::
                Traits::LocalBasisType::Traits::RangeFieldType RF;

        // face geometry
        RF face_volume = ig.geometry().volume();
        const Dune::FieldVector<DF,IG::dimension>&
                inside_local = Dune::ReferenceElements<DF,IG::dimension>::general(ig.inside()->type()).position(0,0);
        const Dune::FieldVector<DF,IG::dimension>&
                outside_local = Dune::ReferenceElements<DF,IG::dimension>::general(ig.outside()->type()).position(0,0);

        // distance between cell centers in global coordinates
        Dune::FieldVector<DF,IG::dimension>
                inside_global = ig.inside()->geometry().center();
        Dune::FieldVector<DF,IG::dimension>
                outside_global = ig.outside()->geometry().center();
        inside_global -= outside_global;
        RF distance = inside_global.two_norm();

        for (std::size_t k = 0;k<TP::COMPONENTS;k++)
        {
            // evaluate diffusion coefficients
            typename TP::Traits::RangeFieldType D_inside = tp.D(*(ig.inside()),inside_local,k);
            typename TP::Traits::RangeFieldType D_outside = tp.D(*(ig.outside()),outside_local,k);
            typename TP::Traits::RangeFieldType D_avg = 2.0/(1.0/(D_inside+1E-40) + 1.0/(D_outside+1E-40));

            // diffusive flux
            r_s.accumulate(lfsu_s,k,-(D_avg*(x_n(lfsu_n,k)-x_s(lfsu_s,k))/distance)*face_volume);
            r_n.accumulate(lfsu_n,k,(D_avg*(x_n(lfsu_n,k)-x_s(lfsu_s,k))/distance)*face_volume);
        }

    }

#ifdef ANALYTICAL_JACOBIAN
    // jacobian of skeleton term
    template<typename IG, typename LFSU, typename X, typename LFSV, typename M>
    void jacobian_skeleton (const IG& ig,
                            const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                            const LFSU& lfsu_n, const X& x_n, const LFSV& lfsv_n,
                            M& mat_ss, M& mat_sn,
                            M& mat_ns, M& mat_nn) const
    {
        typedef typename LFSV::template Child<0>::Type Space;

        // domain and range field type
        typedef typename Space::Traits::FiniteElementType::
                Traits::LocalBasisType::Traits::DomainFieldType DF;
        typedef typename Space::Traits::FiniteElementType::
                Traits::LocalBasisType::Traits::RangeFieldType RF;

        // face geometry
        RF face_volume = ig.geometry().volume();
        const Dune::FieldVector<DF,IG::dimension>&
                inside_local = Dune::ReferenceElements<DF,IG::dimension>::general(ig.inside()->type()).position(0,0);
        const Dune::FieldVector<DF,IG::dimension>&
                outside_local = Dune::ReferenceElements<DF,IG::dimension>::general(ig.outside()->type()).position(0,0);


        // distance between cell centers in global coordinates
        Dune::FieldVector<DF,IG::dimension>
                inside_global = ig.inside()->geometry().center();
        Dune::FieldVector<DF,IG::dimension>
                outside_global = ig.outside()->geometry().center();
        inside_global -= outside_global;
        RF distance = inside_global.two_norm();

        for (std::size_t k = 0;k<TP::COMPONENTS;k++)
        {

            // evaluate diffusion coefficients
            typename TP::Traits::RangeFieldType D_inside = tp.D(*(ig.inside()),inside_local,k);
            typename TP::Traits::RangeFieldType D_outside = tp.D(*(ig.outside()),outside_local,k);
            typename TP::Traits::RangeFieldType D_avg = 2.0/(1.0/(D_inside+1E-40) + 1.0/(D_outside+1E-40));

            // diffusive flux
            mat_ss.accumulate(lfsu_s,k,lfsu_s,k,D_avg/distance*face_volume);
            mat_sn.accumulate(lfsu_n,k,lfsu_n,k,-D_avg/distance*face_volume);

            mat_ns.accumulate(lfsu_s,k,lfsu_s,k,-D_avg/distance*face_volume);
            mat_nn.accumulate(lfsu_n,k,lfsu_n,k,D_avg/distance*face_volume);

        }

    }
#endif //ANALYTICAL_JACOBIAN

    //! to be called once before each time step
    void preStep (typename TP::Traits::RangeFieldType time, typename TP::Traits::RangeFieldType dt,
                  int stages)
    {
        tp.preStep(time,dt,stages);
        ra.preStep(time,dt,stages);
    }

    //! to be called once before each stage
    void preStage (typename TP::Traits::RangeFieldType time, int r)
    {

    }

    //! to be called once at the end of each stage
    void postStage ()
    {
    }

private:
    static RA & raDefault()
    {
        static RA ra;
        return ra;
    }


    TP& tp;
    RA& ra;
};


/** a local operator for the storage operator
     *
     * \f{align*}{
     \int_\Omega c(x,t) uv dx
     * \f}
     *
     * version where c(x,t) may become zero.
     */
template<class TP>
class MulticomponentCCFVTemporalOperator
      : public NumericalJacobianVolume<MulticomponentCCFVTemporalOperator<TP> >,
        public NumericalJacobianApplyVolume<MulticomponentCCFVTemporalOperator<TP> >,
        public FullVolumePattern,
        public LocalOperatorDefaultFlags,
        public InstationaryLocalOperatorDefaultMethods<typename TP::Traits::RangeFieldType>
{
public:
    // pattern assembly flags
    enum { doPatternVolume = true };

    // residual assembly flags
    enum { doAlphaVolume = true };

    MulticomponentCCFVTemporalOperator (TP& tp_)
        : tp(tp_)
    {
    }

    // volume integral depending on test and ansatz functions
    template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
    void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
    {
        for (std::size_t k = 0;k<TP::COMPONENTS;k++)
        {
            // residual contribution
            r.accumulate(lfsu,k,x(lfsu,k)*eg.geometry().volume());
        }
    }


#ifdef ANALYTICAL_JACOBIAN
    // jacobian of volume term
    template<typename EG, typename LFSU, typename X, typename LFSV, typename M>
    void jacobian_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv,
                          M& mat) const
    {
        for (std::size_t k = 0;k<TP::COMPONENTS;k++)
        {
            // residual contribution
            mat.accumulate(lfsu,k,lfsu,k,eg.geometry().volume());
        }
    }
#endif //ANALYTICAL_JACOBIAN

    //! to be called once before each time step
    void preStep (typename TP::Traits::RangeFieldType time, typename TP::Traits::RangeFieldType dt,
                  int stages)
    {
        tp.preStep(time,dt,stages);
        tp.setTimeTarget(time,dt);
    }

    //! to be called once before each stage
    void preStage (typename TP::Traits::RangeFieldType time, int r)
    {
    }

    //! to be called once at the end of each stage
    void postStage ()
    {
    }


    //suggest time step, asked after first stage
    typename TP::Traits::RangeFieldType suggestTimestep (typename TP::Traits::RangeFieldType dt) const
    {
        return std::numeric_limits<typename TP::Traits::RangeFieldType>::max(); //initial value should be big enough
    }

private:
    TP& tp;
    typename TP::Traits::RangeFieldType time;
};

}
}

#endif
