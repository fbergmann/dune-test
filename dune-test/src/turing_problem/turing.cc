// -*- tab-width: 4; indent-tabs-mode: nil -*-
/** \file

    \brief Solve two-component diffusion-reaction system
    with conforming finite elements
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
#include <dune/common/mpihelper.hh>
#include <dune/common/parametertreeparser.hh>

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

#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
#include <dune/grid/io/file/gmshreader.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/grid/common/gridinfo.hh>

#include <dune/istl/bvector.hh>
#include <dune/istl/operators.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/io.hh>
#include <dune/istl/superlu.hh>


#include <dune/copasi/utilities/newton.hh>
#include <dune/copasi/utilities/newtonutilities.hh>
#include <dune/copasi/utilities/timemanager.hh>

#include<dune/pdelab/newton/newton.hh>
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

#include <dune/pdelab/gridoperator/onestep.hh>
#include <dune/pdelab/stationary/linearproblem.hh>
#include <dune/pdelab/instationary/onestep.hh>

#include <dune/pdelab/finiteelementmap/pkfem.hh>
#include <dune/pdelab/finiteelementmap/qkfem.hh>
#include <dune/pdelab/finiteelementmap/rannacher_turek2dfem.hh>

#include <dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include <dune/pdelab/gridfunctionspace/genericdatahandle.hh>
#include <dune/pdelab/gridfunctionspace/interpolate.hh>
#include <dune/pdelab/gridfunctionspace/vtk.hh>
#include <dune/pdelab/gridfunctionspace/subspace.hh>

#include <dune/pdelab/gridoperator/gridoperator.hh>
#include <dune/pdelab/gridoperator/onestep.hh>

#include <dune/pdelab/stationary/linearproblem.hh>

#include <dune/pdelab/instationary/onestep.hh>



#include"local_operator.hh"
#include"turing_initial.hh"


/** \brief Control time step after reaction.

    If some concentration is negative, then it returns false.
    Otherwise it returns true

    \tparam GV          The grid view
    \tparam V           Vector backend (containing DOF)
*/
template<class V>
bool controlReactionTimeStep (V &v)
{
  for (auto it=v.begin();it!=v.end();++it)
    if (*it < 0.)
      {
        std::cout << "concentration negative: " << *it << std::endl;
        return false;
      }
  std::cout << "concentration control was successfull" << std::endl;
  return true;
}


template<class GV, int k = 2>
void run (const GV& gv, Dune::ParameterTree & param)
{
  Dune::gridinfo(gv.grid());
  // <<<1>>> Choose domain and range field type
  typedef typename GV::Grid::ctype DF;
  typedef double RF;
  Dune::Timer watch;

  // output to vtk file (true/false)
  const bool graphics = param.get<bool>("writeVTK", false);

  ConvergenceAdaptiveTimeStepper<RF> timemanager(param);

  // <<<2>>> Make grid function space for the system (Q1 Elements)

  typedef Dune::PDELab::QkLocalFiniteElementMap<GV,DF,RF,k> FEM0;
  FEM0 fem0(gv);
  typedef Dune::PDELab::NoConstraints CON;
  
  typedef Dune::PDELab::ISTLVectorBackend<> VBE0;
  typedef Dune::PDELab::GridFunctionSpace<GV,FEM0,CON,VBE0> GFS0;
  GFS0 gfs0(gv,fem0);


  typedef Dune::PDELab::ISTLVectorBackend <Dune::PDELab::ISTLParameters::static_blocking,2> VBE;               // block size 2
  typedef Dune::PDELab::PowerGridFunctionSpace<GFS0,2,VBE,
                                               Dune::PDELab::EntityBlockedOrderingTag> GFS;
  GFS gfs(gfs0);
  typedef typename GFS::template ConstraintsContainer<RF>::Type CC;

  // some informations
  gfs.update();
  std::cout << "number of DOF =" << gfs.globalSize() << std::endl;

  typedef Dune::PDELab::GridFunctionSubSpace
    <GFS,Dune::TypeTree::TreePath<0> > U0SUB;
  U0SUB u0sub(gfs);
  typedef Dune::PDELab::GridFunctionSubSpace
    <GFS,Dune::TypeTree::TreePath<1> > U1SUB;
  U1SUB u1sub(gfs);

  // <<<3>>> Make FE function


  // reaction terms for local operator
  const RF d_0 = param.sub("Reaction").get<RF>("d0");
  const RF d_1 = param.sub("Reaction").get<RF>("d1");

  const RF k1 = param.sub("Reaction").get<RF>("k1");
  const RF Km = param.sub("Reaction").get<RF>("Km");
  const RF v1 = param.sub("Reaction").get<RF>("v1");
  const RF v2 = param.sub("Reaction").get<RF>("v2");
  const RF V = param.sub("Reaction").get<RF>("V");
  const RF tau = param.sub("Reaction").get<RF>("tau");

  const int integrationorder = param.get<int>("integrationorder");

  typedef Example05LocalOperator LOP;
  LOP lop(d_0,d_1,
          v1, v2, V,k1, Km,
          integrationorder); //(integration order 2 should be high enough
  typedef Example05TimeLocalOperator TLOP;
  TLOP tlop(tau,2*k);
  typedef Dune::PDELab::istl::BCRSMatrixBackend<> MBE;
  MBE mbe(integrationorder == 1 ? 9 : 25);
  typedef Dune::PDELab::GridOperator<GFS,GFS,LOP,MBE,RF,RF,RF,CC,CC> GO0;
  GO0 go0(gfs,gfs,lop,mbe);
  typedef Dune::PDELab::GridOperator<GFS,GFS,TLOP,MBE,RF,RF,RF,CC,CC> GO1;
  GO1 go1(gfs,gfs,tlop,mbe);

  typedef Dune::PDELab::OneStepGridOperator<GO0,GO1> IGO;
  IGO igo(go0,go1);


  // <<<4>>> initial value
  // typedef typename Dune::PDELab::BackendVectorSelector<GFS,RF>::Type U;
  typedef typename IGO::Traits::Domain U;

  U uold(gfs,0.0);
  typedef U0Initial<GV,RF> U0InitialType;
  U0InitialType u0initial(gv);
  typedef U1Initial<GV,RF> U1InitialType;
  U1InitialType u1initial(gv);
  typedef Dune::PDELab::CompositeGridFunction<U0InitialType,U1InitialType> UInitialType;
  UInitialType uinitial(u0initial,u1initial);
  Dune::PDELab::interpolate(uinitial,gfs,uold);

 
  // <<<5>>> Select a linear solver backend
  // now I'm using direct solver superlu (good for testing)
  // but for big problems, iterative solvers should be suitable
#if HAVE_SUPERLU
  typedef Dune::PDELab::ISTLBackend_SEQ_SuperLU LS;
  LS ls (param.sub("Newton").get<int>("LinearVerbosity", 0));
  //LS ls(true);
#else
  typedef Dune::PDELab::ISTLBackend_SEQ_BCGS_SSOR LS;
  LS ls(5000,true);
#endif
  //LS ls(gfs,param.sub("Newton").get<int>("LSMaxIterations", 1000),param.sub("Newton").get<int>("lineark", 5),param.sub("Newton").get<int>("LinearVerbosity", 0));

  // <<<6>>> Solver for non-linear problem per stage
  typedef Dune::PDELab::NewtonReaction<IGO,LS,U> PDESOLVER;
  PDESOLVER pdesolver(igo,ls);
  pdesolver.setLineSearchStrategy(PDESOLVER::hackbuschReuskenAcceptBest); //strategy for linesearch
  typedef Dune::PDELab::NewtonParameters NewtonParameters;
  NewtonParameters newtonparameters(param.sub("Newton"));
  newtonparameters.set(pdesolver);

  // <<<7>>> time-stepper
  //Dune::PDELab::Alexander2Parameter<RF> method;
  Dune::PDELab::ImplicitEulerParameter<RF> method;
  Dune::PDELab::OneStepMethod<RF,IGO,PDESOLVER,U,U> osm(method,igo,pdesolver);
  osm.setVerbosityLevel(param.sub("Verbosity").get<int>("Instationary", 0));

  // discrete grid functions
  typedef Dune::PDELab::DiscreteGridFunction<U0SUB,U> U0DGF;
  U0DGF u0dgf(u0sub,uold);
  typedef Dune::PDELab::DiscreteGridFunction<U1SUB,U> U1DGF;
  U1DGF u1dgf(u1sub,uold);

  // <<<8>>> graphics for initial guess
  std::stringstream basename;
  basename << "turing_" << k << "Q" << k;
  typedef Dune::PVDWriter<GV> PVDWriter;
  PVDWriter pvdwriter(gv,basename.str());
  pvdwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<U0DGF>(u0dgf,"u0"));
  pvdwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<U1DGF>(u1dgf,"u1"));

  // <<<9>>> time loop
  U unew(gfs,0.0);
  unew = uold;


  if (graphics)
    {
      pvdwriter.write(0);
    }


  double dt = timemanager.getTimeStepSize();
  while (!timemanager.finalize())
    {

      dt = timemanager.getTimeStepSize();
      if (gv.comm().rank() == 0)
        {
          std::cout << "======= solve reaction problem =======\n";
          std::cout << "time " << timemanager.getTime() << " dt " << timemanager.getTimeStepSize() << std::endl;
          watch.reset();
        }
      watch.reset();
      try {
        // do time step
        osm.apply(timemanager.getTime(),dt,uold,unew);
        if (!controlReactionTimeStep(unew))
          {
            timemanager.notifyFailure();
            unew = uold;
            continue;

          }
        uold = unew;

        if (gv.comm().rank() == 0)
          std::cout << "... done\n";
      }
      // newton linear search error
      catch (Dune::PDELab::NewtonLineSearchError) {
        if (gv.comm().rank() == 0)
          std::cout << "Newton Linesearch Error" << std::endl;
        timemanager.notifyFailure();
        unew = uold;
        continue;
      }
      catch (Dune::PDELab::NewtonNotConverged) {
        if (gv.comm().rank() == 0)
          std::cout << "Newton Convergence Error" << std::endl;
        timemanager.notifyFailure();
        unew = uold;
        continue;
      }
      catch (Dune::PDELab::NewtonLinearSolverError) {
        if (gv.comm().rank() == 0)
          std::cout << "Newton Linear Solver Error" << std::endl;
        timemanager.notifyFailure();
        unew = uold;
        continue;
      }
      catch (Dune::ISTLError) {
        if (gv.comm().rank() == 0)
          std::cout << "ISTL Error" << std::endl;
        timemanager.notifyFailure();
        unew = uold;
        continue;
      }

      if (gv.comm().rank() == 0)
        std::cout << "... took : " << watch.elapsed() << " s" << std::endl;

      // notify success in this timestep
      timemanager.notifySuccess(5);

      if (graphics && timemanager.isTimeForOutput())
        {
          pvdwriter.write(timemanager.getTime());
        }

    }




}


//===============================================================
// Main program with grid setup
//===============================================================

int main(int argc, char** argv)
{
  try{
    //Maybe initialize Mpi
    Dune::MPIHelper& helper = Dune::MPIHelper::instance(argc, argv);
    if(Dune::MPIHelper::isFake)
      std::cout<< "This is a sequential program." << std::endl;
    else
      {
        if(helper.rank()==0)
          std::cout << "parallel run on " << helper.size() << " process(es)" << std::endl;
      }

    std::string configfile = argv[0]; configfile += ".conf";

    // parse cmd line parameters
    if (argc > 1 && argv[1][0] != '-')
      configfile = argv[1];
    for (int i = 1; i < argc; i++)
      {
        if (std::string(argv[i]) == "--help" || std::string(argv[i]) == "-h")
          {
            if(helper.rank()==0)
              std::cout << "usage: ./turing <configfile> [OPTIONS]" << std::endl;
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
        Dune::array<int,2> N(Dune::fill_array<int,2>(1));
        N[0] = param.get<int>("Domain.nx");
        N[1] = param.get<int>("Domain.ny");
        std::bitset<2> periodic(false);
        int overlap = param.get<int>("overlap", 0);;
        Dune::YaspGrid<2> grid(L,N,periodic,overlap);
        grid.globalRefine(param.get<int>("Domain.refine",0));

        typedef Dune::YaspGrid<2>::LeafGridView GV;
        const GV& gv=grid.leafGridView();
        // solve problem
        run<GV,1>(gv,param);
      }
  }
  catch (Dune::Exception &e){
    std::cerr << "Dune reported error: " << e << std::endl;
    return 1;
  }
  /*catch (...){
    std::cerr << "Unknown exception thrown!" << std::endl;
    return 1;
    }*/
}
