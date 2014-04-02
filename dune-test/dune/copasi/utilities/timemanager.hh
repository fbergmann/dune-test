// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set ts=4 sw=2 et sts=2:
#ifndef DUNE_DYCAP_TIME_STEP_MANAGER_HH
#define DUNE_DYCAP_TIME_STEP_MANAGER_HH

#include <memory>
#include <dune/common/timer.hh>
#include <dune/common/parametertreeparser.hh>
#include <dune/pdelab/gridoperator/common/timesteppingparameterinterface.hh>
#include <dune/pdelab/instationary/onestep.hh>

/**
   \brief Manager object to provide time step size according to
   convergence rate of the Newton solver within a given range.

   A more sophisticated version which derives from this implementation
   is the MultiSpanConvergenceAdaptiveTimeStepper which allows to
   define multiple intervals with different ranges for valid step
   sizes.

   \tparam Real A real floating point type to represent the time
   values.

   \tparam The type of the coefficient vectors as expected by the
   driven one step method.
*/
template <typename Real>
class ConvergenceAdaptiveTimeStepper
{

protected:
  Real dx;

protected:
  Real final_time;
  Real dt;
  Real time;
  Real dt_min;
  Real dt_eps;
  Real dt_max;
  Real dt_plot;
  Real next_plot_time;
  Real backup_dt;
  int increase_rate;
  Dune::Timer total_time;
  Dune::Timer step_time;
  long int total_steps;
  bool notified;
  unsigned int verbosity_level;

  //! \brief Constructor for convenient use in derived classes
  ConvergenceAdaptiveTimeStepper()
    : backup_dt(0), total_steps(0), notified(true), verbosity_level(0)
  {
      if (verbosity_level)
          std::cout << "ConvergenceAdaptiveTimeStepper" << std::endl;
  }

public:
  typedef Real TimeType;

  /**
     \brief Constructor

     The ParameterTree object is required to provide the following
     keys:

     <ul>

     <li> propagation.time: The final time to which to integrate
     (beginning with zero).

     <li> propagation.dt: The initial time step.

     <li> propagation.dt_min: The minimal time step.

     <li> propagation.dt_max: The maximal time step.

     <li>propagation.dt_eps: An numerical epsilon value appropriate
     for the current time step choice. It should be big enough, that
     your problem may be solved for dt == dt_eps without running into
     numerical difficulties.

     <li> propagation.dt_increase_rate: The average convergence rate
     below which the time step is increased

     <li> propagation.dt_plot: The time step at which to provide
     visualization output.

     </ul>

     \param[in] configuration A parameter tree which is expected to
     provide the keys given above.

  */
  ConvergenceAdaptiveTimeStepper(const Dune::ParameterTree & configuration)
    : final_time(configuration.get<Real>("Timeloop.time")),
      dt(configuration.get<Real>("Timeloop.dt")),
      time(0.0),
      dt_min(configuration.get<Real>("Timeloop.dt_min")),
      dt_eps(configuration.get<Real>("Timeloop.dt_eps",1.e-6)),
      dt_max(configuration.get<Real>("Timeloop.dt_max")),
      dt_plot(configuration.get<Real>("Timeloop.dt_plot")),
      next_plot_time(time),
      backup_dt(0),
      increase_rate(configuration.get<Real>("Timeloop.increase_rate")),
      total_steps(0),
      verbosity_level(configuration.get<int>("Verbosity.verbosity"))
  {}

  /**
     \brief Provides the current time step size.

  */
  Real getTimeStepSize()
  {
    // If the success or failure of last step was not notified, then
    // just return the old step.
    if(!notified)
      {
        if (verbosity_level)
          std::cout << "not notified, dt is " << dt << std::endl;
        return dt;
      }

    notified = false;

    // Assert that the step limits are satisfied
    if(dt > dt_max) dt = dt_max;
    if(dt < dt_min) dt = dt_min;

    // Check whether final time is achieved.
    if(time >= final_time - dt_eps){
      if (verbosity_level)
        {
          const Real sec_elapsed = total_time.elapsed();
          const Real hours_elapsed = sec_elapsed / 3600.0;
          std::cout << "Finished time integration at time " << time << std::endl;
          std::cout << "Total computation time: " << hours_elapsed << " h ( or " << sec_elapsed << " s )" << std::endl;
        }
      return 0;
    }

    // Plot if time for plotting is passed and set next plot time
    if(time > next_plot_time - dt_eps){
      next_plot_time = time + dt_plot;
    }
    if(dt_plot == 0.0)
      next_plot_time = time + dt;

    // Adapt time step size with regard to the next plot time
    if(std::abs(time + dt - next_plot_time) < dt_eps
       || time + dt > next_plot_time){

      // Only set backup time if the change is significant
      if(next_plot_time - time < 0.5 * dt)
        {
          backup_dt = dt;
          if (verbosity_level)
            std::cout << "set backup " << dt << std::endl;
        }
      dt = next_plot_time - time;
      if (verbosity_level)
        std::cout << "Adapt time step to reach next output time: dt = " << dt << std::endl;
    }

    // Adapt time step size with regard to the final time
    if(std::abs(time + dt - final_time) < dt_eps
       || time + dt > final_time){
      dt = final_time - time;
      if (verbosity_level)
        std::cout << "Adapt time step to reach final time: dt = " << dt << std::endl;
    }

    if (verbosity_level)
      std::cout << "normal return dt is " << dt << std::endl;
    return dt;
  }

  //! \brief Check if an output is due according to the setting of the
  //! dt_plot parameter.
  bool isTimeForOutput()
  {
    //this variant avoids next timespan changes
    return (time > next_plot_time - dt_eps);
  }

  //! \brief Return current time.
  Real getTime(){
    return time;
  }

  //! \brief Return dt_eps.
  Real getDtEps(){
    return dt_eps;
  }

  //! \brief Notify the stepper about a failed time step.
  void notifyFailure()
  {
    increase_rate = 0;
    notified = true;
    backup_dt = 0.;
    // Begin with time integration step
    dt *= 0.5;
    if (dt<dt_min) {
      if (verbosity_level)
        std::cerr << "Time step too small. Exiting..." << std::endl;
      exit(1);
    }
  }

  //! \brief Notify the stepper about a succeeded time step and
  //! provide the average convergence rate of the applied nonlinear
  //! solver.
  //! The impl. now: timestep is bigger only if previous
  //! Newton-iteration converged succesfully
  void notifySuccess(const Real convergence_rate)
  {
    notified = true;

    // Postprocess time step size with regard to Newton convergence
    // rate
    total_steps++;

    time += dt;

    if(backup_dt){
      dt = backup_dt;

      if (verbosity_level)
        std::cout << "Output time reached. Reset dt = " << dt << std::endl;
    }
    // The convergence dependend time step adaption is not applied
    // if an output time was approached exactly, as the convergence
    // rate might be less meaningful.
    else{
      if(convergence_rate < increase_rate){
        if(dt * 2.0 <= dt_max){
          dt *= 2.0;
          if (verbosity_level)
            std::cout << "Increasing time step to " << dt << std::endl;
        }
        else if(dt < dt_max){
          dt = dt_max;
          if (verbosity_level)
            std::cout << "Increasing time step to " << dt << std::endl;
        }
      }
    }
    if (verbosity_level)
      std::cout << "Time step computed in " << step_time.elapsed() << " s" << std::endl;

    //if (increase_rate==0)
    increase_rate=1000;
    backup_dt = 0;
  }

  bool finalize()
  {
    return time>=final_time-1.e-6;
  }

};


namespace Dune {
  namespace PDELab {

    template<typename RF>
    class TimeSteppingMethods
    {
    public:
      TimeSteppingMethods()
      {
        methods.clear();
        // explicit methods
        methods.insert(std::make_pair("ExplicitEuler",std::shared_ptr<TimeSteppingParameterInterface<RF> >(new Dune::PDELab::ExplicitEulerParameter<RF>())));
        methods.insert(std::make_pair("Heun",std::shared_ptr<TimeSteppingParameterInterface<RF> >(new Dune::PDELab::HeunParameter<RF>())));
        methods.insert(std::make_pair("Shu3",std::shared_ptr<TimeSteppingParameterInterface<RF> >(new Dune::PDELab::Shu3Parameter<RF>())));
        methods.insert(std::make_pair("RK4",std::shared_ptr<TimeSteppingParameterInterface<RF> >(new Dune::PDELab::RK4Parameter<RF>())));
        // implicit methods
        methods.insert(std::make_pair("ImplicitEuler",std::shared_ptr<TimeSteppingParameterInterface<RF> >(new Dune::PDELab::ImplicitEulerParameter<RF>())));
        methods.insert(std::make_pair("OneStepTheta",std::shared_ptr<TimeSteppingParameterInterface<RF> >(new Dune::PDELab::OneStepThetaParameter<RF>(0.5))));
        methods.insert(std::make_pair("Alexander2",std::shared_ptr<TimeSteppingParameterInterface<RF> >(new Dune::PDELab::Alexander2Parameter<RF>())));
        methods.insert(std::make_pair("FractionalStep",std::shared_ptr<TimeSteppingParameterInterface<RF> >(new Dune::PDELab::FractionalStepParameter<RF>())));
        methods.insert(std::make_pair("Alexander3",std::shared_ptr<TimeSteppingParameterInterface<RF> >(new Dune::PDELab::Alexander3Parameter<RF>())));
      }

      template<typename OSM>
      void setTimestepMethod(OSM& osm, std::string method_name)
      {

        if (methods.count(method_name) )
          {
            std::shared_ptr<TimeSteppingParameterInterface<RF> > method = methods.find(method_name)->second;
            osm.setMethod(*method);
          }
        else  DUNE_THROW(Exception,"method " << method_name << "can NOT be found");
      }

    private:
      std::map<std::string, std::shared_ptr<TimeSteppingParameterInterface<RF> > > methods;
    };

  }
}


#endif



