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

#include <dune/pdelab/common/function.hh>
#include <dune/pdelab/common/vtkexport.hh>
#include <dune/pdelab/instationary/pvdwriter.hh>

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

#include "componentparameters.hh"
#include "local_operator.hh"
#include "initial_conditions.hh"
#include "reactionadapter.hh"


/** \brief Control time step after reaction.

    If some concentration is negative, then it returns false.
    Otherwise it returns true

    \tparam V           Vector backend (containing DOF)
*/
template<typename V>
bool controlReactionTimeStep (V &v)
{
    //for (auto it=v.begin(); it!=v.end();++it)
    //    if (*it<0) return false;

    return true;
}

template<class GV>
void run (const GV& gv, Dune::ParameterTree & param)
{
    const bool verbosity = param.sub("Verbosity").get<bool>("verbosity", false);
    if (verbosity)
    Dune::gridinfo(gv.grid());
    typedef typename GV::Grid::ctype DF;
    typedef double RF;
    const int dim = GV::dimension;
    Dune::Timer watch;


    // output to vtk file (true/false)
    const bool graphics = param.get<bool>("writeVTK", false);

    ConvergenceAdaptiveTimeStepper<RF> timemanager(param);

    // for each component parameters different classes
    typedef DiffusionParameter<GV,RF> DP;
    DP dp1(param, "u");
    DP dp2(param, "X");
    DP dp3(param, "v");
    DP dp4(param, "w");


    typedef Dune::PDELab::MulticomponentDiffusion<GV,RF,DP,4> VCT;
    VCT vct(dp1,dp2,dp3,dp4);

    // <<<2>>> Make grid function space
    typedef Dune::PDELab::P0LocalFiniteElementMap<DF,RF,dim> FEM;
    FEM fem(Dune::GeometryType(Dune::GeometryType::cube,dim));
    typedef Dune::PDELab::P0ParallelConstraints CON;
    typedef Dune::PDELab::ISTLVectorBackend<> VBE0;
    typedef Dune::PDELab::GridFunctionSpace<GV,FEM,CON,VBE0> GFS;

    typedef Dune::PDELab::ISTLVectorBackend
            <Dune::PDELab::ISTLParameters::static_blocking,VCT::COMPONENTS> VBE;
    typedef Dune::PDELab::PowerGridFunctionSpace<GFS,VCT::COMPONENTS,VBE,
            Dune::PDELab::EntityBlockedOrderingTag> TPGFS;
    watch.reset();
    CON con;
    GFS gfs(gv,fem,con);
    TPGFS tpgfs(gfs);
    if (verbosity) std::cout << "=== function space setup " <<  watch.elapsed() << " s" << std::endl;

    // some informations
    gfs.update();
    if (verbosity) std::cout << "number of DOF =" << gfs.globalSize() << std::endl;


    // <<<2b>>> make subspaces for visualization
    typedef Dune::PDELab::GridFunctionSubSpace<TPGFS,Dune::TypeTree::TreePath<0> > U0SUB;
    U0SUB u0sub(tpgfs);
    typedef Dune::PDELab::GridFunctionSubSpace<TPGFS,Dune::TypeTree::TreePath<1> > U1SUB;
    U1SUB u1sub(tpgfs);
    typedef Dune::PDELab::GridFunctionSubSpace<TPGFS,Dune::TypeTree::TreePath<2> > U2SUB;
    U2SUB u2sub(tpgfs);
    typedef Dune::PDELab::GridFunctionSubSpace<TPGFS,Dune::TypeTree::TreePath<3> > U3SUB;
    U3SUB u3sub(tpgfs);


    // make constraints map and initialize it from a function (boundary conditions)
    typedef typename TPGFS::template ConstraintsContainer<RF>::Type CC;
    CC cg;
    cg.clear();
    Dune::PDELab::constraints(tpgfs,cg,false);

    typedef ReactionAdapter<RF> RA;
    RA ra(param);

    typedef Dune::PDELab::MulticomponentCCFVSpatialDiffusionOperator<VCT,RA> LOP;
    LOP lop(vct,ra); //local operator including reaction adapter
    typedef Dune::PDELab::MulticomponentCCFVTemporalOperator<VCT> TLOP;
    TLOP tlop(vct);
    typedef Dune::PDELab::istl::BCRSMatrixBackend<> MBE;
    MBE mbe(5); // Maximal number of nonzeroes per row can be cross-checked by patternStatistics().
    // grid operators
    typedef Dune::PDELab::GridOperator<TPGFS,TPGFS,LOP,MBE,RF,RF,RF,CC,CC> CGO0;
    CGO0 cgo0(tpgfs,cg,tpgfs,cg,lop,mbe);
    typedef Dune::PDELab::GridOperator<TPGFS,TPGFS,TLOP,MBE,RF,RF,RF,CC,CC> CGO1;
    CGO1 cgo1(tpgfs,cg,tpgfs,cg,tlop,mbe);
    // one step grid operator for instationary problems

    typedef Dune::PDELab::OneStepGridOperator<CGO0,CGO1,true> IGO;
    IGO igo(cgo0,cgo1);


    // <<<7>>> make vector for old time step and initialize
    typedef typename IGO::Traits::Domain V;
    V uold(tpgfs,0.0);
    V unew(tpgfs,0.0);

    // initial conditions
    typedef U0Initial<GV,RF> U0InitialType;
    U0InitialType u0initial(gv,param);
    typedef U1Initial<GV,RF> U1InitialType;
    U1InitialType u1initial(gv,param);
    typedef U2Initial<GV,RF> U2InitialType;
    U2InitialType u2initial(gv,param);
    typedef U3Initial<GV,RF> U3InitialType;
    U3InitialType u3initial(gv,param);


    typedef Dune::PDELab::CompositeGridFunction<U0InitialType,U1InitialType,U2InitialType,U3InitialType> UInitialType; //new
    UInitialType uinitial(u0initial,u1initial,u2initial,u3initial); //new
    Dune::PDELab::interpolate(uinitial,tpgfs,uold);
    unew = uold;

  // <<<5>>> Select a linear solver backend
#if HAVE_MPI
#if HAVE_SUPERLU
  typedef Dune::PDELab::ISTLBackend_BCGS_AMG_SSOR<IGO> LS;
  LS ls(tpgfs,param.sub("Newton").get<int>("LSMaxIterations", 100),param.sub("Newton").get<int>("LinearVerbosity", 0));
#else
  typedef Dune::PDELab::ISTLBackend_OVLP_BCGS_SSORk<GFS,CC> LS;
  LS ls(tpgfs,cc,5000,5,param.sub("Newton").get<int>("LinearVerbosity", 0));
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
    typedef Dune::PDELab::Newton<IGO,LS,V> PDESOLVER;
    PDESOLVER pdesolver(igo,ls);
    pdesolver.setLineSearchStrategy(PDESOLVER::hackbuschReuskenAcceptBest); //strategy for linesearch
    typedef Dune::PDELab::NewtonParameters NewtonParameters;
    NewtonParameters newtonparameters(param.sub("Newton"));
    newtonparameters.set(pdesolver);

    // <<<7>>> time-stepper
    Dune::PDELab::Alexander2Parameter<RF> method;
    Dune::PDELab::OneStepMethod<RF,IGO,PDESOLVER,V,V> osm(method,igo,pdesolver);
    Dune::PDELab::TimeSteppingMethods<RF> tsmethods;
    tsmethods.setTimestepMethod(osm,param.get<std::string>("timesolver"));

    osm.setVerbosityLevel(param.sub("Verbosity").get<int>("Instationary", 0));

    char basename[255];
    sprintf(basename,"%s-%01d",param.get<std::string>("VTKname","").c_str(),param.sub("Domain").get<int>("refine"));

    typedef Dune::PVDWriter<GV> PVDWriter;
    PVDWriter pvdwriter(gv,basename,Dune::VTK::conforming);
    // discrete grid functions
    typedef Dune::PDELab::DiscreteGridFunction<U0SUB,V> U0DGF;
    U0DGF u0dgf(u0sub,uold);
    pvdwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<U0DGF>(u0dgf,"u"));
    typedef Dune::PDELab::DiscreteGridFunction<U1SUB,V> U1DGF;
    U1DGF u1dgf(u1sub,uold);
    pvdwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<U1DGF>(u1dgf,"X"));
    typedef Dune::PDELab::DiscreteGridFunction<U2SUB,V> U2DGF;
    U2DGF u2dgf(u2sub,uold);
    pvdwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<U2DGF>(u2dgf,"v"));
    typedef Dune::PDELab::DiscreteGridFunction<U3SUB,V> U3DGF;
    U3DGF u3dgf(u3sub,uold);
    pvdwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<U3DGF>(u3dgf,"w"));



    // <<<9>>> time loop

    if (graphics)
        pvdwriter.write(0);


    double dt = timemanager.getTimeStepSize();
    while (!timemanager.finalize())
    {
        dt = timemanager.getTimeStepSize();
        if (gv.comm().rank() == 0)
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

            if (!controlReactionTimeStep(gv, unew))
            {
                timemanager.notifyFailure();
                unew = uold;
                continue;

            }
            uold = unew;

            if (gv.comm().rank() == 0 && verbosity)
                std::cout << "... done\n";
        }
        // newton linear search error
        catch (Dune::PDELab::NewtonLineSearchError) {
            if (gv.comm().rank() == 0 && verbosity)
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
            if (gv.comm().rank() == 0 && verbosity)
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
            if (gv.comm().rank() == 0 && verbosity )
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
            if (gv.comm().rank() == 0 && verbosity )
                std::cout << "ISTL Error" << std::endl;
            timemanager.notifyFailure();
            unew = uold;
            if (!verbosity)
            {
              std::cout << "x" << std::flush;
            }
            continue;
        }

        if (gv.comm().rank() == 0 && verbosity)
            std::cout << "... took : " << watch.elapsed() << " s" << std::endl;

        // notify success in this timestep
        timemanager.notifySuccess(5);

        if (graphics && timemanager.isTimeForOutput())
        {
          pvdwriter.write(timemanager.getTime());
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

    bool inside=(((((pow((x - 3.00000000000000000E+002 * 4.19999999999999980E-001), 2.00000000000000000E+000) / pow(3.00000000000000000E+002 * 3.70000000000000000E-001, 2.00000000000000000E+000) + pow((y - 3.00000000000000000E+002 / 2.00000000000000000E+000), 2.00000000000000000E+000) / pow(3.00000000000000000E+002 * 2.50000000000000000E-001, 2.00000000000000000E+000)) < 1.00000000000000000E+000 || (y < (( - 3.00000000000000000E+002 * 1.00000000000000010E-001) + x) && y > (3.00000000000000000E+002 * 1.10000000000000010E+000 + ( - x)) && x < 3.00000000000000000E+002 * 9.00000000000000020E-001)) && !((pow((x - 3.00000000000000000E+002 * 2.50000000000000000E-001), 2.00000000000000000E+000) / pow(8.00000000000000020E-002 * 3.00000000000000000E+002, 2.00000000000000000E+000) + pow((y - 3.00000000000000000E+002 * 4.50000000000000010E-001), 2.00000000000000000E+000) / pow(8.00000000000000020E-002 * 3.00000000000000000E+002, 2.00000000000000000E+000)) < 1.00000000000000000E+000))) ? 1.00000000000000000E+000 : 0.00000000000000000E+000);
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
                    std::cout << "usage: ./pnas_fish <configfile> [OPTIONS]" << std::endl;
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
            }

            mdgrid.preUpdateSubDomains();
            mdgrid.updateSubDomains();
            mdgrid.postUpdateSubDomains();

            typedef MDGrid::SubDomainGrid SDGrid;
            const SDGrid& sdgrid = mdgrid.subDomain(0);
            SDGrid::LeafGridView sdgv = sdgrid.leafGridView();

            // solve problem
            run(sdgv,param);
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

