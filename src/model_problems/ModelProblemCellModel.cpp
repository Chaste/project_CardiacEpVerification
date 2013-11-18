
#include "ModelProblemCellModel.hpp"
#include "OdeSystemInformation.hpp"
#include "HeartConfig.hpp"
#include <cmath>


/* Constructor */
ModelProblemCellModel::ModelProblemCellModel(
        boost::shared_ptr<AbstractIvpOdeSolver> pOdeSolver,
        boost::shared_ptr<AbstractStimulusFunction> pIntracellularStimulus,
        double beta)
    : AbstractCardiacCell(pOdeSolver, 4, 0, pIntracellularStimulus),
      mBeta(beta)
{
    mpSystemInfo = OdeSystemInformation<ModelProblemCellModel>::Instance();
    Init();
}

/* Destructor */
ModelProblemCellModel::~ModelProblemCellModel(void)
{
}

/* Derivatives function */
void ModelProblemCellModel::EvaluateYDerivatives(double time, const std::vector<double> &rY, std::vector<double>& rDY)
{
    double V = rY[0];
    double u1 = rY[1];
    double u2 = rY[2];
    double u3 = rY[3];

    rDY[0] = 0.0; // normally it would be -iionic here (as rY[0]=voltage), but cell model will never be run outside of tissue simulation so this value will be ignored

    /* Define functions prescribed in paper */
    rDY[1] = (u1+u3-V)*(u1+u3-V)*u2*u2 + (u1+u3-V)*u2*u2*(V-u3)/2.0;
    rDY[2] = -(u1+u3-V)*u2*u2*u2;
    rDY[3] = 0.0;
}

/* Define Iionic function prescribed in paper */
double ModelProblemCellModel::GetIIonic(const std::vector<double>* pStateVariables)
{
    if (!pStateVariables)
    {
        pStateVariables = &mStateVariables; // (cell model GetIionic should be set up to use internal state if no state variables are passed in)
    }

    double V = (*pStateVariables)[0];
    double u1 = (*pStateVariables)[1];
    double u2 = (*pStateVariables)[2];
    double u3 = (*pStateVariables)[3];

    double Cm = HeartConfig::Instance()->GetCapacitance();
    double chi = HeartConfig::Instance()->GetSurfaceAreaToVolumeRatio();

    return -Cm*(u1+u3-V)*u2*u2*(V-u3)/2.0 + mBeta*(V-u3)/chi;
}



/* This is required but these initial values will be overwritten by the cell factory */
template<>
void OdeSystemInformation<ModelProblemCellModel>::Initialise(void)
{
    /*
     * State variables
     */
    this->mVariableNames.push_back("V");
    this->mVariableUnits.push_back("mV");
    this->mInitialConditions.push_back(0.0);

    this->mVariableNames.push_back("u1");
    this->mVariableUnits.push_back("");
    this->mInitialConditions.push_back(0.0);

    this->mVariableNames.push_back("u2");
    this->mVariableUnits.push_back("");
    this->mInitialConditions.push_back(0.0);

    this->mVariableNames.push_back("u3");
    this->mVariableUnits.push_back("");
    this->mInitialConditions.push_back(0.0);

    this->mInitialised = true;
}
