
#ifndef _MODELPROBLEMCELLMODEL_HPP_
#define _MODELPROBLEMCELLMODEL_HPP_

#include "AbstractCardiacCell.hpp"
#include "AbstractStimulusFunction.hpp"
#include <vector>


/* A custom cell model for the model problems, defined as in the paper
 * Note: cell models are usually auto-generated from cellml files.
 */
class ModelProblemCellModel : public AbstractCardiacCell
{
private:
    /* The free constant beta */
    double mBeta;

public:
    /* Constructor, takes in a solver and stimulus (so that interface matches other cell models), and the constant beta. */
    ModelProblemCellModel(boost::shared_ptr<AbstractIvpOdeSolver> pOdeSolver,
                          boost::shared_ptr<AbstractStimulusFunction> pIntracellularStimulus,
                          double beta);

    /* Destructor */
    ~ModelProblemCellModel();

    /* The function saying what the derivatives of the state variables are */
    void EvaluateYDerivatives(double time, const std::vector<double> &rY, std::vector<double>& rDY);

    /* The ionic current function */
    double GetIIonic(const std::vector<double>* pStateVariables=NULL);
};

#endif //_MODELPROBLEMCELLMODEL_HPP_
