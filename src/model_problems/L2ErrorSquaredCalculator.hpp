#ifndef L2ERRORSQUAREDCALCULATOR_HPP_
#define L2ERRORSQUAREDCALCULATOR_HPP_


#include "AbstractFunctionalCalculator.hpp"


/*
 * This class computes the integral of (V_computed - V_exact)^2 (or the same with extracellular cellular potential),
 * ie the L2 error, squared, of the computed solution, at the given time.
 *
 * Inherits from the AbstractFunctionalCalculator, which allows the user to compute integrals over the domain by
 * just specifying the integrand.
 *
 */
template<unsigned DIM>
class L2ErrorSquaredCalculator : public AbstractFunctionalCalculator<DIM,DIM,1>
{
private:
    AbstractScalarFunction<DIM>* mpExactSolution; // class which knows how to compute the exact solution
    double mTime;

    // For bath problems, there are dummy voltage variables at bath nodes, these nodes must be ignored in
    // the integral (only integrate over tissue)
    bool mIgnoreBath;

    double GetIntegrand(ChastePoint<DIM>& rX,
                        c_vector<double,1>& rU, // what is in rU depends on whether it is the voltage or phi_e that is given by the calling class
                        c_matrix<double,1,DIM>& rGradU)
    {
        double soln_exact = mpExactSolution->Compute(mTime,rX); // replace with zero to compute q for Section 4.3.1 of paper
        return (rU(0) - soln_exact)*(rU(0) - soln_exact);
    }

public:
    L2ErrorSquaredCalculator(AbstractScalarFunction<DIM>* pExactSolution, double time, bool ignoreBath)
      : mpExactSolution(pExactSolution),
        mTime(time),
        mIgnoreBath(ignoreBath)
    {
    }

    bool ShouldSkipThisElement(Element<DIM,DIM>& rElement)
    {
        if( mIgnoreBath && rElement.GetAttribute()==HeartRegionCode::GetValidBathId() )
        {
            return true;
        }
        return false;
    }
};


#endif //L2ERRORSQUAREDCALCULATOR_HPP_
