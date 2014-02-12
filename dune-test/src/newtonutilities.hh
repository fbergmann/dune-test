// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set ts=4 sw=2 et sts=2:
#ifndef DUNE_DYCAP_NEWTON_UTILS
#define DUNE_DYCAP_NEWTON_UTILS

#include<dune/common/parametertree.hh>

namespace Dune {
namespace PDELab {

//! interpret a parameter tree as a set of options for the newton solver
/**

       example configuration:

       \code
       [NewtonParameters]

       ReassembleThreshold = 0.1
       LineSearchMaxIterations = 10
       MaxIterations = 7
       AbsoluteLimit = 1e-6
       Reduction = 1e-4
       LinearReduction = 1e-3
       LineSearchDamping  = 0.9
       Verbosity = 2
       \endcode
    */
struct NewtonParameters : public Dune::ParameterTree
{
public:
    NewtonParameters(const Dune::ParameterTree p) :
        Dune::ParameterTree(p)
    {
    }
    NewtonParameters & operator = (const Dune::ParameterTree & p)
    {
        static_cast<ParameterTree &>(*this) = p;
        return *this;
    }

    //! apply the parameter set to an instance of Dune::PDELab::Newton
    template<typename N>
    void set(N & newton) const
    {
        if (this->hasKey("ReassembleThreshold"))
            newton.setReassembleThreshold(
                        this->get<double>("ReassembleThreshold"));
        if (this->hasKey("LineSearchMaxIterations"))
            newton.setLineSearchMaxIterations(
                        this->get<int>("LineSearchMaxIterations"));
        if (this->hasKey("MaxIterations"))
            newton.setMaxIterations(
                        this->get<int>("MaxIterations"));
        if (this->hasKey("AbsoluteLimit"))
            newton.setAbsoluteLimit(
                        this->get<double>("AbsoluteLimit"));
        if (this->hasKey("Reduction"))
            newton.setReduction(
                        this->get<double>("Reduction"));
        if (this->hasKey("LinearReduction"))
            newton.setMinLinearReduction(
                        this->get<double>("LinearReduction"));
        if (this->hasKey("LineSearchDamping"))
            newton.setLineSearchDampingFactor(
                        this->get<double>("LineSearchDamping"));
        if (this->hasKey("Verbosity"))
            newton.setVerbosityLevel(
                        this->get<int>("Verbosity"));
    }
};

}
}
#endif
