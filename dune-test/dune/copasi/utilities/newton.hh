// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=8 sw=4 sts=4:
#ifndef DUNE_PDELAB_MODIFIED_NEWTON_HH
#define DUNE_PDELAB_MODIFIED_NEWTON_HH

#include <dune/pdelab/newton/newton.hh>

namespace Dune
{
    namespace PDELab
    {

        /** \brief Modified Newton method
            
            The problem with chemical reactions and it's models
            is that current concentration (if it is not treated in the model)
            can become negative, which is not required at all!

            Solution for this problem can be in each linear search iteration
            in Newton method control, if the corresponding DOF are < 0. 
            If yes, the update will be modified to set negative values to zero.

            This is first vanilla implementation (change of Newton from PDELab).
            Newton will be in PDELab changed in the future that users will be able
            to change its functionality.

            The second solution would be to control solution after each timestep
            and if the solution is negative, set time to previous time Told and solution
            to previous solution (e.g. with timestepmanager) and solve the problem 
            from Told with smaller timestep.

         */
        template<class GOS, class TrlV, class TstV>
        class NewtonLineSearchReaction : public virtual NewtonBase<GOS,TrlV,TstV>
        {
            typedef GOS GridOperator;
            typedef TrlV TrialVector;
            typedef TstV TestVector;

            typedef typename TestVector::ElementType RFType;

        public:
            enum Strategy { noLineSearch,
                            hackbuschReusken,
                            hackbuschReuskenAcceptBest };

            NewtonLineSearchReaction(GridOperator& go, TrialVector& u_)
                : NewtonBase<GOS,TrlV,TstV>(go,u_)
                , strategy(hackbuschReusken)
                , maxit(10)
                , damping_factor(0.5)
            {}

            NewtonLineSearchReaction(GridOperator& go)
                : NewtonBase<GOS,TrlV,TstV>(go)
                , strategy(hackbuschReusken)
                , maxit(10)
                , damping_factor(0.5)
            {}

            void setLineSearchStrategy(Strategy strategy_)
            {
                strategy = strategy_;
            }

            void setLineSearchMaxIterations(unsigned int maxit_)
            {
                maxit = maxit_;
            }

            void setLineSearchDampingFactor(RFType damping_factor_)
            {
                damping_factor = damping_factor_;
            }

            template <typename V>
            void solutioncontrol(V& v)
            {
                for (auto it=(v.begin()); it!=(v.end()); it++)
                    if (*it<0)
                        {
                            if (this->verbosity_level >= 4)
                                std::cout << "solution " << *it << "set to zero" << std::endl;
                            *it = 0;
                        }
            }


            virtual void line_search(TrialVector& z, TestVector& r)
            {
                if (strategy == noLineSearch)
                    {
                        this->u->axpy(-1.0, z);                     // TODO: vector interface
                        this->defect(r);
                        return;
                    }

                if (this->verbosity_level >= 4)
                    std::cout << "      Performing line search..." << std::endl;
                RFType lambda = 1.0;
                RFType best_lambda = 0.0;
                RFType best_defect = this->res.defect;
                TrialVector prev_u(*this->u);  // TODO: vector interface
                unsigned int i = 0;
                ios_base_all_saver restorer(std::cout); // store old ios flags


                TrialVector curr_u(*this->u);
                while (1)
                    {
                        if (this->verbosity_level >= 4)
                            std::cout << "          trying line search damping factor:   "
                                      << std::setw(12) << std::setprecision(4) << std::scientific
                                      << lambda;
                        
                        curr_u.axpy(-lambda, z);

                        *this->u = curr_u;
                        solutioncontrol(*this->u);


                        try {
                            this->defect(r);
                            if (this->verbosity_level >= 4)
                                std::cout << " new defect " << this->res.defect << std::endl;
                        }
                        catch (NewtonDefectError)
                            {
                                if (this->verbosity_level >= 4)
                                    std::cout << "          Nans detected" << std::endl;
                            }       // ignore NaNs and try again with lower lambda

                        if (this->res.defect <= (1.0 - lambda/4) * this->prev_defect)
                            {
                                if (this->verbosity_level >= 4)
                                    std::cout << "          line search converged" << std::endl;
                                break;
                            }

                        if (this->res.defect < best_defect)
                            {
                                best_defect = this->res.defect;
                                best_lambda = lambda;
                            }

                        if (++i >= maxit)
                            {
                                if (this->verbosity_level >= 4)
                                    std::cout << "          max line search iterations exceeded" << std::endl;
                                switch (strategy)
                                    {
                                    case hackbuschReusken:

                                        *this->u = prev_u;
                                        solutioncontrol(*this->u);
                                        this->defect(r);
                                        DUNE_THROW(NewtonLineSearchError,
                                                   "NewtonLineSearch::line_search(): line search failed, "
                                                   "max iteration count reached, "
                                                   "defect did not improve enough");
                                    case hackbuschReuskenAcceptBest:
                                        if (best_lambda == 0.0)
                                            {
                                                *this->u = prev_u;
                                                solutioncontrol(*this->u);
                                                this->defect(r);
                                                DUNE_THROW(NewtonLineSearchError,
                                                           "NewtonLineSearch::line_search(): line search failed, "
                                                           "max iteration count reached, "
                                                           "defect did not improve in any of the iterations");
                                            }
                                        if (best_lambda != lambda)
                                            {
                                                *this->u = prev_u;
                                                this->u->axpy(-best_lambda, z);
                                                solutioncontrol(*this->u);
                                                this->defect(r);
                                            }
                                        break;
                                    case noLineSearch:
                                        break;
                                    }
                                break;
                            }

                        lambda *= damping_factor;
                        curr_u = prev_u;                          // TODO: vector interface
                    }
                if (this->verbosity_level >= 4)
                    std::cout << "          line search damping factor:   "
                              << std::setw(12) << std::setprecision(4) << std::scientific
                              << lambda << std::endl;
            }



        protected:
            Strategy strategy;
            unsigned int maxit;
            RFType damping_factor;
        };




        template<class GOS, class S, class TrlV, class TstV = TrlV>
        class NewtonReaction : public NewtonSolver<GOS,S,TrlV,TstV>
                     , public NewtonTerminate<GOS,TrlV,TstV>
                     , public NewtonLineSearchReaction<GOS,TrlV,TstV>
                     , public NewtonPrepareStep<GOS,TrlV,TstV>
        {
            typedef GOS GridOperator;
            typedef S Solver;
            typedef TrlV TrialVector;
            typedef TstV TestVector;

            typedef typename TestVector::ElementType RFType;

        public:
            NewtonReaction(GridOperator& go, TrialVector& u_, Solver& solver_)
                : NewtonBase<GOS,TrlV,TstV>(go,u_)
                , NewtonSolver<GOS,S,TrlV,TstV>(go,u_,solver_)
                , NewtonTerminate<GOS,TrlV,TstV>(go,u_)
                , NewtonLineSearchReaction<GOS,TrlV,TstV>(go,u_)
                , NewtonPrepareStep<GOS,TrlV,TstV>(go,u_)
            {}
            NewtonReaction(GridOperator& go, Solver& solver_)
                : NewtonBase<GOS,TrlV,TstV>(go)
                , NewtonSolver<GOS,S,TrlV,TstV>(go,solver_)
                , NewtonTerminate<GOS,TrlV,TstV>(go)
                , NewtonLineSearchReaction<GOS,TrlV,TstV>(go)
                , NewtonPrepareStep<GOS,TrlV,TstV>(go)
            {}


        };


    }

}

#endif // DUNE_PDELAB_NEWTON_HH
