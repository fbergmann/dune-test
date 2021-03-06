// -*- tab-width: 4; indent-tabs-mode: nil -*-
/** \file

    \brief Solve two-component diffusion-reaction system
    with finite volume scheme. This class tests the
    multi domain grid.
*/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <math.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <stdlib.h>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/common/static_assert.hh>
#include <dune/common/timer.hh>
#include <dune/common/parametertreeparser.hh>
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
#include <dune/grid/io/file/gmshreader.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/grid/multidomaingrid/indexsets.hh>
#include <dune/grid/multidomaingrid/subdomainset.hh>
#include <dune/grid/multidomaingrid.hh>
#include <dune/grid/common/gridinfo.hh>

#if HAVE_ALBERTA
#include <dune/grid/albertagrid.hh>
#include <dune/grid/albertagrid/dgfparser.hh>
#endif

#if HAVE_UG
#include <dune/grid/uggrid.hh>
#endif

#if HAVE_ALUGRID
#include <dune/grid/alugrid.hh>
#include <dune/grid/io/file/dgfparser/dgfalu.hh>
#include <dune/grid/io/file/dgfparser/dgfparser.hh>
#endif

#include <dune/istl/bvector.hh>
#include <dune/istl/operators.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/io.hh>
#include <dune/istl/superlu.hh>

#include <dune/copasi/utilities/newton.hh>
#include <dune/copasi/utilities/newtonutilities.hh>
#include <dune/copasi/utilities/timemanager.hh>
#include <dune/copasi/utilities/sbmlhelper.hh>
#include <dune/copasi/utilities/solutioncontrol.hh>
#include <dune/copasi/utilities/mdpvdwriter.hh>

#include <dune/pdelab/common/function.hh>
#include <dune/pdelab/common/vtkexport.hh>


#include <dune/pdelab/backend/istl/tags.hh>
#include <dune/pdelab/backend/istl/descriptors.hh>
#include <dune/pdelab/backend/istl/bcrsmatrixbackend.hh>

#include <dune/pdelab/backend/istlvectorbackend.hh>
#include <dune/pdelab/backend/istlmatrixbackend.hh>
#include <dune/pdelab/backend/istlsolverbackend.hh>
#include <dune/pdelab/backend/seqistlsolverbackend.hh>

#include <dune/pdelab/constraints/conforming.hh>
#include <dune/pdelab/constraints/constraints.hh>
#include <dune/pdelab/ordering/orderingbase.hh>

#include <dune/pdelab/finiteelementmap/p0fem.hh>
#include <dune/pdelab/constraints/p0.hh>
#include <dune/pdelab/constraints/constraints.hh>
#include <dune/pdelab/gridoperator/gridoperator.hh>
#include <dune/pdelab/gridoperator/onestep.hh>
#include <dune/pdelab/stationary/linearproblem.hh>
#include <dune/pdelab/instationary/onestep.hh>
#include <dune/pdelab/gridfunctionspace/subspace.hh>

#include <dune/pdelab/multidomain/multidomaingridfunctionspace.hh>
#include <dune/pdelab/multidomain/coupling.hh>
#include <dune/pdelab/multidomain/subproblemlocalfunctionspace.hh>
#include <dune/pdelab/multidomain/gridoperator.hh>
#include <dune/pdelab/multidomain/subproblem.hh>
#include <dune/pdelab/constraints/conforming.hh>
#include <dune/pdelab/multidomain/constraints.hh>
#include <dune/pdelab/multidomain/interpolate.hh>
#include <dune/pdelab/multidomain/vtk.hh>

#include "componentparameters.hh"
#include "local_operator.hh"
#include "initial_conditions.hh"
#include "reactionadapter.hh"
#include "proportionalflowcoupling.hh"

/** \brief Control time step after reaction.

    If some concentration is negative, then it returns false.
    Otherwise it returns true

    \tparam V           Vector backend (containing DOF)
*/
template<typename V>
bool controlReactionTimeStep (V &v)
{
    for (auto it=v.begin(); it!=v.end();++it)
        if (*it<0) return false;

    return true;
}

template<class MDGrid, class GV>
void run (MDGrid& grid, const GV& gv0, const GV& gv1, Dune::ParameterTree & param)
{
    const bool verbosity = param.sub("Verbosity").get<bool>("verbosity", false);
    if (verbosity)
    Dune::gridinfo(gv0.grid());
    typedef typename GV::Grid::ctype DF;
    typedef double RF;
    const int dim = GV::dimension;
    Dune::Timer watch;


    // output to vtk file (true/false)
    const bool graphics = param.get<bool>("writeVTK", false);

    ConvergenceAdaptiveTimeStepper<RF> timemanager(param);

    // for each component parameters different classes
    typedef DiffusionParameter<typename MDGrid::LeafGridView,RF> DP;
    DP dp1(param, "species_2");
    DP dp2(param, "species_1");


    typedef Dune::PDELab::MulticomponentDiffusion<typename MDGrid::LeafGridView,RF,DP,2> VCT;
    VCT vct(dp1,dp2);

    // <<<2>>> Make grid function space
    typedef Dune::PDELab::P0LocalFiniteElementMap<DF,RF,dim> FEM;
    FEM fem(Dune::GeometryType(Dune::GeometryType::cube,dim));
    typedef Dune::PDELab::P0ParallelConstraints CON;
    typedef Dune::PDELab::ISTLVectorBackend<> VBE0;
    typedef Dune::PDELab::GridFunctionSpace<GV,FEM,CON,VBE0> GFS;

    typedef Dune::PDELab::ISTLVectorBackend<> VBE;
    //            <Dune::PDELab::ISTLParameters::static_blocking,VCT::COMPONENTS> VBE;
    typedef Dune::PDELab::PowerGridFunctionSpace<GFS,VCT::COMPONENTS,VBE,
            Dune::PDELab::EntityBlockedOrderingTag> TPGFS;
    watch.reset();
    CON con;
    GFS gfs0(gv0,fem,con);
    TPGFS tpgfs0(gfs0);
    tpgfs0.child(0).name("species0_0");
    tpgfs0.child(1).name("species0_1");
    GFS gfs1(gv1,fem,con);
    TPGFS tpgfs1(gfs1);
    tpgfs1.child(0).name("species1_0");
    tpgfs1.child(1).name("species1_1");

    typedef Dune::PDELab::MultiDomain::MultiDomainGridFunctionSpace<
      MDGrid,
      VBE,
      Dune::PDELab::LexicographicOrderingTag,
      TPGFS,
      TPGFS
      > MultiGFS;

    MultiGFS multigfs(grid,tpgfs0,tpgfs1);

    if (verbosity) std::cout << "=== function space setup " <<  watch.elapsed() << " s" << std::endl;

    // some informations
    //multigfs.update();
    //if (verbosity) std::cout << "number of DOF =" << multigfs.globalSize() << std::endl;


    // <<<2b>>> make subspaces for visualization
    //typedef Dune::PDELab::GridFunctionSubSpace<MultiGFS,Dune::TypeTree::TreePath<0> > U0SUB;
    //U0SUB u0sub(multigfs);
    //typedef Dune::PDELab::GridFunctionSubSpace<MultiGFS,Dune::TypeTree::TreePath<1> > U1SUB;
    //U1SUB u1sub(multigfs);

    typedef ReactionAdapter<RF> RA;
    RA ra(param);

    typedef Dune::PDELab::MulticomponentCCFVSpatialDiffusionOperator<VCT,RA> LOP;
    LOP lop(vct,ra); //local operator including reaction adapter
    typedef Dune::PDELab::MulticomponentCCFVTemporalOperator<VCT> TLOP;
    TLOP tlop(vct);

    typedef Dune::PDELab::MultiDomain::SubDomainEqualityCondition<MDGrid> Condition;

    Condition c0(0);
    Condition c1(1);

    typedef Dune::PDELab::MultiDomain::SubProblem<MultiGFS,MultiGFS,LOP,Condition, 0> LeftSubProblem_dt0;
    typedef Dune::PDELab::MultiDomain::SubProblem<MultiGFS,MultiGFS,TLOP,Condition, 0> LeftSubProblem_dt1;

    LeftSubProblem_dt0 left_sp_dt0(lop,c0);
    LeftSubProblem_dt1 left_sp_dt1(tlop,c0);

    typedef Dune::PDELab::MultiDomain::SubProblem<MultiGFS,MultiGFS,LOP,Condition,1> RightSubProblem_dt0;
    typedef Dune::PDELab::MultiDomain::SubProblem<MultiGFS,MultiGFS,TLOP,Condition, 1> RightSubProblem_dt1;

    RightSubProblem_dt0 right_sp_dt0(lop,c1);
    RightSubProblem_dt1 right_sp_dt1(tlop,c1);

    ContinuousValueContinuousFlowCoupling<RF> proportionalFlowCoupling(4,0.1);

    typedef Dune::PDELab::MultiDomain::Coupling<LeftSubProblem_dt0,RightSubProblem_dt0,ContinuousValueContinuousFlowCoupling<RF> > Coupling;
    Coupling coupling(left_sp_dt0,right_sp_dt0,proportionalFlowCoupling);

    typedef typename MultiGFS::template ConstraintsContainer<RF>::Type CC;
    CC cg;

    // need to get bounds out of the DiffusionParameter
    auto constraints = Dune::PDELab::MultiDomain::constraints<RF>
    (
       multigfs,
       Dune::PDELab::MultiDomain::constrainSubProblem(left_sp_dt0,vct),
       Dune::PDELab::MultiDomain::constrainSubProblem(right_sp_dt0,vct)
    );

    constraints.assemble(cg);

    typedef Dune::PDELab::istl::BCRSMatrixBackend<> MBE;
    MBE mbe(27); // Maximal number of nonzeroes per row can be cross-checked by patternStatistics().

    typedef Dune::PDELab::MultiDomain::GridOperator<
      MultiGFS,MultiGFS,
      MBE,RF,RF,RF,CC,CC,
      LeftSubProblem_dt0,
      RightSubProblem_dt0,
      Coupling
      > GridOperator_dt0;

    typedef Dune::PDELab::MultiDomain::GridOperator<
      MultiGFS,MultiGFS,
      MBE,RF,RF,RF,CC,CC,
      LeftSubProblem_dt1,
      RightSubProblem_dt1
      > GridOperator_dt1;

    typedef Dune::PDELab::OneStepGridOperator<GridOperator_dt0,GridOperator_dt1> GridOperator;

    // <<<7>>> make vector for old time step and initialize
    typedef typename GridOperator::Traits::Domain V;
    V uold(multigfs,0.0);
    V unew(multigfs,0.0);

    // initial conditions
    typedef U0Initial<typename MDGrid::LeafGridView,RF> U0InitialType;
    typedef U1Initial<typename MDGrid::LeafGridView,RF> U1InitialType;
    U0InitialType u0initial(grid.leafGridView(),param);
    U1InitialType u1initial(grid.leafGridView(),param);


    //typedef Dune::PDELab::CompositeGridFunction<U0InitialType,U1InitialType> UInitialType; //new
    typedef Dune::PDELab::PowerGridFunction<U0InitialType,2> UInitialType;
    UInitialType uinitial(u0initial); //new
    //Dune::PDELab::interpolate(uinitial,multigfs,uold);

    Dune::PDELab::MultiDomain::interpolateOnTrialSpace
            (multigfs,uold,uinitial,left_sp_dt0,uinitial,right_sp_dt0);

    //Dune::PDELab::MultiDomain::interpolateOnTrialSpace
    //        (multigfs,uold,uinitial,left_sp_dt0,uinitial,right_sp_dt0);


    unew = uold;

    GridOperator_dt0 go_dt_0(multigfs,multigfs,
                             cg,cg,
                             mbe,
                             left_sp_dt0,
                             right_sp_dt0,
                             coupling);

    GridOperator_dt1 go_dt_1(multigfs,multigfs,
                             cg,cg,
                             mbe,
                             left_sp_dt1,
                             right_sp_dt1);

    GridOperator igo(go_dt_0,go_dt_1);


  // <<<5>>> Select a linear solver backend
#if HAVE_MPI
#if HAVE_SUPERLU
  typedef Dune::PDELab::ISTLBackend_BCGS_AMG_SSOR<GridOperator> LS;
  LS ls(multigfs,param.sub("Newton").get<int>("LSMaxIterations", 100),param.sub("Newton").get<int>("LinearVerbosity", 0));
#else
  typedef Dune::PDELab::ISTLBackend_OVLP_BCGS_SSORk<GFS,CC> LS;
  LS ls(multigfs,cc,5000,5,param.sub("Newton").get<int>("LinearVerbosity", 0));
#endif
#else //!parallel
#if HAVE_SUPERLU
    typedef Dune::PDELab::ISTLBackend_SEQ_SuperLU LS;
    LS ls (param.sub("Newton").get<int>("LinearVerbosity", 0));
#else
    typedef Dune::PDELab::ISTLBackend_SEQ_BCGS_SSOR LS;
  LS ls(5000,param.sub("Newton").get<int>("LinearVerbosity", 0));
#endif
#endif

    // <<<6>>> Solver for non-linear problem per stage
    typedef Dune::PDELab::Newton<GridOperator,LS,V> PDESOLVER;
    PDESOLVER pdesolver(igo,ls);
    pdesolver.setLineSearchStrategy(PDESOLVER::hackbuschReuskenAcceptBest); //strategy for linesearch
    typedef Dune::PDELab::NewtonParameters NewtonParameters;
    NewtonParameters newtonparameters(param.sub("Newton"));
    newtonparameters.set(pdesolver);

    // <<<7>>> time-stepper
    Dune::PDELab::Alexander2Parameter<RF> method;
    Dune::PDELab::OneStepMethod<RF,GridOperator,PDESOLVER,V,V> osm(method,igo,pdesolver);
    Dune::PDELab::TimeSteppingMethods<RF> tsmethods;
    tsmethods.setTimestepMethod(osm,param.get<std::string>("timesolver"));

    osm.setVerbosityLevel(param.sub("Verbosity").get<int>("Instationary", 0));

    char basename[255];
    sprintf(basename,"%s-%01d",param.get<std::string>("VTKname","").c_str(),param.sub("Domain").get<int>("refine"));

    typedef Dune::MdPVDWriter<MDGrid, MultiGFS, 2> PVDWriter;
    PVDWriter pvdwriter(grid, multigfs, basename, Dune::VTK::conforming);

    // <<<9>>> time loop

    if (graphics)
    {
        pvdwriter.write(0, uold);
    }


    double dt = timemanager.getTimeStepSize();
    while (!timemanager.finalize())
    {

        dt = timemanager.getTimeStepSize();

        if (gv0.comm().rank() == 0)
        {
            if (verbosity)
            {
                std::cout << "======= solve reaction problem ======= from " << timemanager.getTime() << "\n";
                std::cout << "time " << timemanager.getTime() << " dt " << timemanager.getTimeStepSize() << std::endl;
            }
        }

        watch.reset();

        try
        {
            // do time step
            osm.apply(timemanager.getTime(),dt,uold,unew);

            if (!controlReactionTimeStepParallel(gv0, unew))
            {
                timemanager.notifyFailure();
                unew = uold;
                continue;

            }
            uold = unew;

            if (gv0.comm().rank() == 0 && verbosity)
                std::cout << "... done\n";
        }
        // newton linear search error
        catch (Dune::PDELab::NewtonLineSearchError) {
            if (gv0.comm().rank() == 0 && verbosity)
                std::cout << "Newton Linesearch Error" << std::endl;
            timemanager.notifyFailure();
            unew = uold;
            if (!verbosity)
            {
              std::cout << "x" << std::flush;
            }
            continue;
        }
        catch (Dune::PDELab::NewtonNotConverged) {
            if (gv0.comm().rank() == 0 && verbosity)
                std::cout << "Newton Convergence Error" << std::endl;
            timemanager.notifyFailure();
            unew = uold;
            if (!verbosity)
            {
              std::cout << "x" << std::flush;
            }
            continue;
        }
        catch (Dune::PDELab::NewtonLinearSolverError) {
            if (gv0.comm().rank() == 0 && verbosity )
                std::cout << "Newton Linear Solver Error" << std::endl;
            timemanager.notifyFailure();
            unew = uold;
            if (!verbosity)
            {
              std::cout << "x" << std::flush;
            }
            continue;
        }
        catch (Dune::ISTLError) {
            if (gv0.comm().rank() == 0 && verbosity )
                std::cout << "ISTL Error" << std::endl;
            timemanager.notifyFailure();
            unew = uold;
            if (!verbosity)
            {
              std::cout << "x" << std::flush;
            }
            continue;
        }

        if (gv0.comm().rank() == 0 && verbosity)
            std::cout << "... took : " << watch.elapsed() << " s" << std::endl;

        // notify success in this timestep
        timemanager.notifySuccess(5);

        if (graphics && timemanager.isTimeForOutput())
        {
          pvdwriter.write(timemanager.getTime(), uold);
          if (!verbosity)
            std::cout << "o" << std::flush;
        }
        else if (!verbosity)
        {
          std::cout << "." << std::flush;
        }
    }

    std::cout << std::endl;
}


bool isInside(const Dune::FieldVector<double, 2>& point, const Dune::FieldVector<double, 2>& dimension)
{

    const auto& x = point[0];
    const auto& y = point[1];
    const auto& width = dimension[0];
    const auto& height = dimension[1];

    bool inside=((((pow(x - 2.50000000000000000E+001, 2.00000000000000000E+000)) + (pow(y - 2.50000000000000000E+001, 2.00000000000000000E+000))) < 2.00000000000000000E+001 * 2.00000000000000000E+001) ? 1.00000000000000000E+000 : 0.00000000000000000E+000);
    return inside;


}

//===============================================================
// Main program with grid setup
//===============================================================

int main(int argc, char** argv)
{
    try{
        //Maybe initialize Mpi
        Dune::MPIHelper& helper = Dune::MPIHelper::instance(argc, argv);

        std::string configfile = argv[0]; configfile += ".conf";

        // parse cmd line parameters
        if (argc > 1 && argv[1][0] != '-')
            configfile = argv[1];
        for (int i = 1; i < argc; i++)
        {
            if (std::string(argv[i]) == "--help" || std::string(argv[i]) == "-h")
            {
                if(helper.rank()==0)
                    std::cout << "usage: ./turing_circle <configfile> [OPTIONS]" << std::endl;
                return 0;
            }
        }

        // read parameters from file
        Dune::ParameterTree param;
        Dune::ParameterTreeParser parser;
        parser.readINITree(configfile, param);

        int dim=param.sub("Domain").get<int>("dim", 2);

        if (dim==2)
        {
            // make grid
            Dune::FieldVector<double,2> L;
            L[0] = param.get<double>("Domain.height");
            L[1] = param.get<double>("Domain.width");

            std::array<int,2> N;
            N[0] = param.get<int>("Domain.nx");
            N[1] = param.get<int>("Domain.ny");

            std::bitset<2> periodic(false);
            int overlap = param.get<int>("overlap", 0);

            typedef Dune::YaspGrid<2> HostGrid;
            HostGrid hostgrid(helper.getCommunicator(),L,N,periodic,overlap);
            hostgrid.globalRefine(param.get<int>("Domain.refine",0));

            typedef Dune::MultiDomainGrid<HostGrid,Dune::mdgrid::ArrayBasedTraits<2,8,8> > MDGrid;
            MDGrid mdgrid(hostgrid,true);
            typedef typename MDGrid::LeafGridView MDGV;
            typedef typename MDGV::template Codim<0>::Iterator Iterator;
            typedef typename MDGrid::SubDomainIndexType SubDomainIndexType;

            MDGV mdgv = mdgrid.leafGridView();
            mdgrid.startSubDomainMarking();
            for(Iterator it = mdgv.template begin<0>(); it != mdgv.template end<0>(); ++it)
            {
                SubDomainIndexType subdomain = 0;
                mdgrid.removeFromAllSubDomains(*it);

                // figure out, whether it is part of the
                // region, and if so add it like so:

                const auto& center = it->geometry().center();
                if (isInside(center, L))
                {
                    mdgrid.addToSubDomain(subdomain,*it);
                }
                else
                {
                    mdgrid.addToSubDomain(1,*it);
                }
            }

            mdgrid.preUpdateSubDomains();
            mdgrid.updateSubDomains();
            mdgrid.postUpdateSubDomains();

            typedef MDGrid::SubDomainGrid SDGrid;
            const SDGrid& sdgrid0 = mdgrid.subDomain(0);
            const SDGrid& sdgrid1 = mdgrid.subDomain(1);
            SDGrid::LeafGridView sdgv0 = sdgrid0.leafGridView();
            SDGrid::LeafGridView sdgv1 = sdgrid1.leafGridView();

            // solve problem
            run(mdgrid,sdgv0,sdgv1,param);
        }
    }
    catch (Dune::Exception &e){
        std::cerr << "Dune reported error: " << e << std::endl;
        return 1;
    }
    //catch (...){
    //  std::cerr << "Unknown exception thrown!" << std::endl;
    //  return 1;
    //}
}
