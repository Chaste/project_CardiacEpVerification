#ifndef H1SEMINORMERRORSQUAREDCALCULATOR_HPP_
#define H1SEMINORMERRORSQUAREDCALCULATOR_HPP_


#include "AbstractFunctionalCalculator.hpp"

/*
 * Calculates the integral of (\grad V_computed - \grad V_exact)^2 over the domain
 * (or the same with extracellular cellular potential), ie the H1 semi-norm (second part of the H1 norm) of the
 * error, squared, for the computed solution at the given time.
 *
 * Inherits from the AbstractFunctionalCalculator, which allows the user to compute integrals over the domain by
 * just specifying the integrand.
 *
 */
template<unsigned DIM>
class H1SemiNormErrorSquaredCalculator : public AbstractFunctionalCalculator<DIM,DIM,1>
{
private:
    AbstractScalarFunction<DIM>* mpExactSolution;  // class which knows how to compute the exact solution
    double mTime;

    // For bath problems, there are dummy voltage variables at bath nodes, these nodes must be ignored in
    // the integral (only integrate over tissue)
    bool mIgnoreBath;

    double GetIntegrand(ChastePoint<DIM>& rX,
                        c_vector<double,1>& rU,  // what is in rU depends on whether it is the voltage or phi_e that is given by the calling class
                        c_matrix<double,1,DIM>& rGradU)
    {
        c_vector<double,DIM> grad_soln_exact = mpExactSolution->ComputeGradient(mTime,rX);
        c_vector<double,DIM> diff;
        double ret;
        if(DIM==1)
        {
            diff(0) = grad_soln_exact(0) - rGradU(0,0);
            ret = diff(0)*diff(0);
        }
        else if(DIM==2)
        {
            diff(0) = grad_soln_exact(0) - rGradU(0,0);
            diff(1) = grad_soln_exact(1) - rGradU(0,1);
            ret = diff(0)*diff(0) + diff(1)*diff(1);
        }
        else
        {
            diff(0) = grad_soln_exact(0) - rGradU(0,0);
            diff(1) = grad_soln_exact(1) - rGradU(0,1);
            diff(2) = grad_soln_exact(2) - rGradU(0,2);
            ret = diff(0)*diff(0) + diff(1)*diff(1) + diff(2)*diff(2);
        }
        return ret;
    }

public:
    H1SemiNormErrorSquaredCalculator(AbstractScalarFunction<DIM>* pExactSolution, double time, bool ignoreBath)
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


#endif //H1SEMINORMERRORSQUAREDCALCULATOR_HPP_
