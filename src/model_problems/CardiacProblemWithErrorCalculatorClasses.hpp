#ifndef CARDIACPROBLEMSWITHERRORCALCULATORCLASSES_HPP_
#define CARDIACPROBLEMSWITHERRORCALCULATORCLASSES_HPP_


#include "MonodomainProblem.hpp"
#include "BidomainProblem.hpp"
#include "BidomainWithBathProblem.hpp"
#include "PetscVecTools.hpp"

#include "L2ErrorSquaredCalculator.hpp"         //defined in this project
#include "H1SemiNormErrorSquaredCalculator.hpp" //defined in this project


/*
 *  A class for calculating the L2 and H1 errors of the solution at a given time, taking in the solution vector of the
 *  monodomain (solution = voltage) or bidomain/bidomain-with-bath (solution = voltage and phi_e) problems at that time.
 *  The current values of the L-infinity(L2) and L2(H1) norms of the error, for both variables, are then updated accordingly.
 */
template<unsigned DIM>
class ErrorCalculatorForCardiacProblems
{
public:
    AbstractScalarFunction<DIM>* mpVoltageExactSolution;
    AbstractScalarFunction<DIM>* mpExtracellularPotentialExactSolution;

    double mVoltageLinfL2Error;
    double mVoltageL2H1Error;
    double mExtracellularPotentialLinfL2Error;
    double mExtracellularPotentialL2H1Error;

    Vec mBidomainVoltageSolution;
    Vec mBidomainExtracellularPotentialSolution;

public:
    ErrorCalculatorForCardiacProblems(AbstractScalarFunction<DIM>* pVoltageExactSolution,
                    AbstractScalarFunction<DIM>* pExtracellularPotentialExactSolution = NULL)
        : mpVoltageExactSolution(pVoltageExactSolution),
          mpExtracellularPotentialExactSolution(pExtracellularPotentialExactSolution),
          mVoltageLinfL2Error(0.0),
          mVoltageL2H1Error(0.0),
          mExtracellularPotentialLinfL2Error(0.0),
          mExtracellularPotentialL2H1Error(0.0),
          mBidomainVoltageSolution(NULL)
    {
    }

    ~ErrorCalculatorForCardiacProblems()
    {
        if(mBidomainVoltageSolution)
        {
            PetscTools::Destroy(mBidomainVoltageSolution);
            PetscTools::Destroy(mBidomainExtracellularPotentialSolution);
        }
    }

    void DoErrorCalculation(double time, AbstractTetrahedralMesh<DIM,DIM>& rMesh, Vec solution, bool notMonodomain)
    {
        Vec& r_voltage = solution;

        if(notMonodomain)
        {
            if(mBidomainVoltageSolution==NULL)
            {
                mBidomainVoltageSolution = PetscTools::CreateAndSetVec(rMesh.GetNumNodes(), 0.0);
                mBidomainExtracellularPotentialSolution  = PetscTools::CreateAndSetVec(rMesh.GetNumNodes(), 0.0);
            }


            // Get the voltage and phi_e out of the full solution (=[V0 phi0 V1 phi1 .. VN phiN]) as two Vecs.
            // Parallel comment: sometimes fails in parallel, probably because the ownership isn't guaranteed to be the same
            // between the three Vecs?
            int lo, hi;
            VecGetOwnershipRange(mBidomainVoltageSolution, &lo, &hi);
            for(int index = lo; index < hi; index++)
            {
                PetscVecTools::SetElement(mBidomainVoltageSolution, index, PetscVecTools::GetElement(solution, 2*index));
                PetscVecTools::SetElement(mBidomainExtracellularPotentialSolution, index, PetscVecTools::GetElement(solution, 2*index+1));
            }

            r_voltage = mBidomainVoltageSolution;
        }

        // voltage
        DoSingleErrorCalculation(mpVoltageExactSolution, time, rMesh, r_voltage, mVoltageLinfL2Error, mVoltageL2H1Error, true);

        // Phi_e
        if(notMonodomain)
        {
            // Note: phi_e at initial condition is just set to zeros (really, the solver should solve the second bidomain equation
            // given initial V to determine what initial phi_e is, but the bidomain solver doesn't do this, and phi_e is just initialised to zero)
            // so we have to not include it. Recall phi_e initial value doesn't affect anything.
            if(fabs(time)>1e-12) // ie if t!=0
            {
                DoSingleErrorCalculation(mpExtracellularPotentialExactSolution, time, rMesh, mBidomainExtracellularPotentialSolution, mExtracellularPotentialLinfL2Error, mExtracellularPotentialL2H1Error, false);
            }
        }

        std::cout << time << "\n";
    }

    void DoSingleErrorCalculation(AbstractScalarFunction<DIM>* pExactSolution, double time, AbstractTetrahedralMesh<DIM,DIM>& rMesh, Vec solution,
                                  double& rLinfL2Error, double& rL2H1Error, bool ignoreBath)
    {
        L2ErrorSquaredCalculator<DIM> l2_calc(pExactSolution,time,ignoreBath);
        double l2_error_sqd = l2_calc.Calculate(rMesh,solution);
        double l2_error = sqrt(l2_error_sqd);
        rLinfL2Error = std::max(rLinfL2Error, l2_error);

        H1SemiNormErrorSquaredCalculator<DIM> h1_calc(pExactSolution,time,ignoreBath);
        double h1_seminorm_error_sqd = h1_calc.Calculate(rMesh,solution);
        double h1_error = sqrt(l2_error_sqd + h1_seminorm_error_sqd);

        double dt = HeartConfig::Instance()->GetPrintingTimeStep();
        double end_time = HeartConfig::Instance()->GetSimulationDuration();

        double factor = ( (fabs(time)<1e-12 || fabs(time-end_time)<1e-12) ? 0.5 : 1.0 );
        rL2H1Error += h1_error*factor*dt;
    }
};



/*
 *  This class is basically MonodomainProblem, but also inherits from the ErrorCalculatorForCardiacProblems defined above, and uses it
 *  to compute the contribution to the L-inf(L2) and L2(H1) errors at the end of every timestep (and beginning of first timestep).
 */
template<unsigned DIM>
class MonodomainProblemWithErrorCalculator : public MonodomainProblem<DIM>, public ErrorCalculatorForCardiacProblems<DIM>
{
public:
    MonodomainProblemWithErrorCalculator(AbstractCardiacCellFactory<DIM>* pCellFactory,
                                         AbstractScalarFunction<DIM>* pExactSolution)
        : MonodomainProblem<DIM>(pCellFactory),
          ErrorCalculatorForCardiacProblems<DIM>(pExactSolution)
    {
    }

    // overloaded method
    void AtBeginningOfTimestep(double time)
    {
        MonodomainProblem<DIM>::AtBeginningOfTimestep(time);
        if(fabs(time)<1e-12) // ie if t=0
        {
            this->mSolution = this->CreateInitialCondition(); // make sure mSolution is set up
            this->DoErrorCalculation(time,*(this->mpMesh),this->GetSolution(),false/*ie no phi_e*/);
            PetscTools::Destroy(this->mSolution);
            this->mSolution = NULL;
        }
    }

    // overloaded method
    void OnEndOfTimestep(double time)
    {
        MonodomainProblem<DIM>::OnEndOfTimestep(time);
        DoErrorCalculation(time,*(this->mpMesh),this->GetSolution(),false/*ie no phi_e*/);
    }
};


/*
 *  Same as MonodomainProblemWithErrorCalculator but for BidomainProblem
 */
template<unsigned DIM>
class BidomainProblemWithErrorCalculator : public BidomainProblem<DIM>, public ErrorCalculatorForCardiacProblems<DIM>
{
public:
    BidomainProblemWithErrorCalculator(AbstractCardiacCellFactory<DIM>* pCellFactory,
                                       AbstractScalarFunction<DIM>* pVoltageExactSolution,
                                       AbstractScalarFunction<DIM>* pExtracellularPotentialExactSolution)
        : BidomainProblem<DIM>(pCellFactory),
          ErrorCalculatorForCardiacProblems<DIM>(pVoltageExactSolution,pExtracellularPotentialExactSolution)
    {
    }

    void AtBeginningOfTimestep(double time)
    {
        BidomainProblem<DIM>::AtBeginningOfTimestep(time);
        if(fabs(time)<1e-12)
        {
            this->mSolution = this->CreateInitialCondition();
            this->DoErrorCalculation(time,*(this->mpMesh),this->GetSolution(),true/*different to above*/);
            PetscTools::Destroy(this->mSolution);
            this->mSolution = NULL;
        }
    }

    void OnEndOfTimestep(double time)
    {
        BidomainProblem<DIM>::OnEndOfTimestep(time);
        DoErrorCalculation(time,*(this->mpMesh),this->GetSolution(),true/*different to monodomain version*/);
    }
};

/*
 *  Same as BidomainProblemWithErrorCalculator but for BidomainWithBathProblem
 */
template<unsigned DIM>
class BidomainWithBathProblemWithErrorCalculator : public BidomainWithBathProblem<DIM>, public ErrorCalculatorForCardiacProblems<DIM>
{
public:
    BidomainWithBathProblemWithErrorCalculator(AbstractCardiacCellFactory<DIM>* pCellFactory,
                                               AbstractScalarFunction<DIM>* pVoltageExactSolution,
                                               AbstractScalarFunction<DIM>* pExtracellularPotentialExactSolution)
        : BidomainWithBathProblem<DIM>(pCellFactory),
          ErrorCalculatorForCardiacProblems<DIM>(pVoltageExactSolution,pExtracellularPotentialExactSolution)
    {
    }

    void AtBeginningOfTimestep(double time)
    {
        BidomainWithBathProblem<DIM>::AtBeginningOfTimestep(time);
        if(fabs(time)<1e-12)
        {
            this->mSolution = this->CreateInitialCondition();
            this->DoErrorCalculation(time,*(this->mpMesh),this->GetSolution(),true);
            PetscTools::Destroy(this->mSolution);
            this->mSolution = NULL;
        }
    }


    void OnEndOfTimestep(double time)
    {
        BidomainWithBathProblem<DIM>::OnEndOfTimestep(time);
        DoErrorCalculation(time,*(this->mpMesh),this->GetSolution(),true);
    }
};


#endif // CARDIACPROBLEMSWITHERRORCALCULATORCLASSES_HPP_

