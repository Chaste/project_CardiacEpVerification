#ifndef MODELPROBLEMEXACTSOLUTIONCLASSES_HPP_
#define MODELPROBLEMEXACTSOLUTIONCLASSES_HPP_

#include <cxxtest/TestSuite.h>


#include "UblasCustomFunctions.hpp" // to include c_vector etc
#include "ChastePoint.hpp"


/*
 * This file just defines various classes that represent the exact solutions of model problems.
 */


/* Useful abstract base class */
template<unsigned DIM>
class AbstractScalarFunction
{
public:
    AbstractScalarFunction()
    {}
    virtual ~AbstractScalarFunction()
    {}

    virtual double Compute(double t, ChastePoint<DIM>& rX)=0;
    virtual c_vector<double,DIM> ComputeGradient(double t, ChastePoint<DIM>& rX)=0;
};




/* Function which returns (1+t)^(1/2) * F(x,y,z) = (1+t)^(1/2) * cos(m1 pi x) * cos(m2 pi y) * cos(m3 pi z) */
template<unsigned DIM>
double TimeScaledFunctionF(ChastePoint<DIM>& rX, double t, unsigned m1, unsigned m2, unsigned m3)
{
    double ret = sqrt(1+t) * cos(m1*M_PI*rX[0]);
    if(DIM>1)
    {
        ret *= cos(m2*M_PI*rX[1]);
    }
    if(DIM>2)
    {
        ret *= cos(m3*M_PI*rX[2]);
    }
    return ret;
}


template<unsigned DIM>
c_vector<double,DIM> TimeScaledGradientOfF(ChastePoint<DIM>& rX, double time, unsigned m1, unsigned m2, unsigned m3)
{
    c_vector<double,DIM> ret;

    if(DIM==1)
    {
        ret(0) = -1.0 * m1*M_PI*sqrt(1+time)*sin(m1*M_PI*rX[0]); // note just using -m1 leads to junk as m1 is unsigned
    }
    else if(DIM==2)
    {
        ret(0) = -1.0 * m1*M_PI*sqrt(1+time)*sin(m1*M_PI*rX[0])*cos(m2*M_PI*rX[1]);
        ret(1) = -1.0 * m2*M_PI*sqrt(1+time)*cos(m1*M_PI*rX[0])*sin(m2*M_PI*rX[1]);
    }
    else
    {
        ret(0) = -1.0 * m1*M_PI*sqrt(1+time)*sin(m1*M_PI*rX[0])*cos(m2*M_PI*rX[1])*cos(m3*M_PI*rX[2]);
        ret(1) = -1.0 * m2*M_PI*sqrt(1+time)*cos(m1*M_PI*rX[0])*sin(m2*M_PI*rX[1])*cos(m3*M_PI*rX[2]);
        ret(2) = -1.0 * m3*M_PI*sqrt(1+time)*cos(m1*M_PI*rX[0])*cos(m2*M_PI*rX[1])*sin(m3*M_PI*rX[2]);
    }
    return ret;
}


/* Function which returns 1 + x y^2 z^3, or 1D or 2D variants */
template<unsigned DIM>
double FunctionG(ChastePoint<DIM>& rX)
{
    double ret;
    if(DIM==1)
    {
        ret = 1.0 + rX[0];
    }
    else if(DIM==2)
    {
        ret = 1.0 + rX[0]*pow(rX[1],2);
    }
    else
    {
        ret = 1.0 + rX[0]*pow(rX[1],2)*pow(rX[2],3);
    }
    return ret;
}

/* Exact solution for the voltage in the monodomain and bidomain model problems */
template<unsigned DIM>
class VoltageExactSolution : public AbstractScalarFunction<DIM>
{
private:
    unsigned mM1;
    unsigned mM2;
    unsigned mM3;

public:
    VoltageExactSolution(unsigned m1, unsigned m2=0, unsigned m3=0)
       : AbstractScalarFunction<DIM>(),
         mM1(m1),
         mM2(m2),
         mM3(m3)
    {
    }

    double Compute(double time, ChastePoint<DIM>& rX)
    {
        return TimeScaledFunctionF(rX,time,mM1,mM2,mM3);
    }

    c_vector<double,DIM> ComputeGradient(double time, ChastePoint<DIM>& rX)
    {
        return TimeScaledGradientOfF(rX,time,mM1,mM2,mM3);
    }
};

/* Exact solution for phi_e in the monodomain and bidomain model problems */
template<unsigned DIM>
class ExtracellularPotentialExactSolution : public AbstractScalarFunction<DIM>
{
private:
    unsigned mM1;
    unsigned mM2;
    unsigned mM3;
    double mK;

public:
    ExtracellularPotentialExactSolution(double k, unsigned m1, unsigned m2=0, unsigned m3=0)
       : AbstractScalarFunction<DIM>(),
         mM1(m1),
         mM2(m2),
         mM3(m3),
         mK(k)
    {
    }

    double Compute(double time, ChastePoint<DIM>& rX)
    {
        return -mK*TimeScaledFunctionF(rX,time,mM1,mM2,mM3);
    }

    c_vector<double,DIM> ComputeGradient(double time, ChastePoint<DIM>& rX)
    {
        return -mK*TimeScaledGradientOfF(rX,time,mM1,mM2,mM3);
    }
};


/* Exact solution for the voltage in bidomain-with-bath model problem */
template<unsigned DIM>
class VoltageExactSolutionBath : public AbstractScalarFunction<DIM>
{
private:
    unsigned mM1;
    double mAlpha;
    double mExtracellularConductivity;

public:
    VoltageExactSolutionBath(unsigned m1, double alpha, double extracellularConductivity)
       : AbstractScalarFunction<DIM>(),
         mM1(m1),
         mAlpha(alpha),
         mExtracellularConductivity(extracellularConductivity)
    {
    }


    double Compute(double t, ChastePoint<DIM>& rX)
    {
        double x = rX[0];
        if(x<0 || x>1)
        {
            NEVER_REACHED;
        }

        return TimeScaledFunctionF<DIM>(rX,t,mM1,0,0) - mAlpha*rX[0]/mExtracellularConductivity;
    }


    c_vector<double,DIM> ComputeGradient(double t, ChastePoint<DIM>& rX)
    {
        c_vector<double,DIM> ret = zero_vector<double>(DIM);
        double x = rX[0];

        if(x<0 || x>1)
        {
            NEVER_REACHED;
        }

        ret = TimeScaledGradientOfF<DIM>(rX,t,mM1,0,0);
        ret(0) -= mAlpha/mExtracellularConductivity;

        return ret;
    }
};


/* Exact solution for phi_e in bidomain-with-bath model problem */
template<unsigned DIM>
class ExtracellularPotentialExactSolutionBath : public AbstractScalarFunction<DIM>
{
private:
    unsigned mM1;
    double mK;
    double mAlpha;
    double mExtracellularConductivity;
    double mBathConductivity;

public:
    ExtracellularPotentialExactSolutionBath(unsigned m1, double k, double alpha, double extracellularConductivity, double bathConductivity)
       : AbstractScalarFunction<DIM>(),
         mM1(m1),
         mK(k),
         mAlpha(alpha),
         mExtracellularConductivity(extracellularConductivity),
         mBathConductivity(bathConductivity)
    {
    }

    double Compute(double t, ChastePoint<DIM>& rX)
    {
        double x = rX[0];
        double constant = -mAlpha*0.5/mExtracellularConductivity;

        double ret;
        if(x<0)
        {
            ChastePoint<DIM> temp(0,0,0);
            ret = -mK*TimeScaledFunctionF(temp, t, mM1, 0, 0) + mAlpha*x/mBathConductivity + constant;
        }
        else if(x<=1+1e-12)
        {
            ret = -mK*TimeScaledFunctionF(rX, t, mM1, 0, 0) + mAlpha*x/mExtracellularConductivity + constant;
        }
        else
        {
            ChastePoint<DIM> temp(1.0,0,0);
            ret = -mK*TimeScaledFunctionF(temp, t, mM1, 0, 0) + mAlpha*1/mExtracellularConductivity + mAlpha*(x-1)/mBathConductivity + constant;
        }
        return ret;
    }


    c_vector<double,DIM> ComputeGradient(double t, ChastePoint<DIM>& rX)
    {
        c_vector<double,DIM> ret = zero_vector<double>(DIM);
        double x = rX[0];
        if(x<0)
        {
            ChastePoint<DIM> temp(0,0,0);
            ret = -mK*TimeScaledGradientOfF(temp, t, mM1, 0, 0);
            ret(0) += mAlpha/mBathConductivity;
        }
        else if(x<=1+1e-12)
        {
            ret = -mK*TimeScaledGradientOfF(rX, t, mM1, 0, 0);
            ret(0) += mAlpha/mExtracellularConductivity;
        }
        else
        {
            ChastePoint<DIM> temp(1.0,0,0);
            ret = -mK*TimeScaledGradientOfF(temp, t, mM1, 0, 0);
            ret(0) += mAlpha/mBathConductivity;
        }
        return ret;
    }
};



#endif // MODELPROBLEMEXACTSOLUTIONCLASSES_HPP_
