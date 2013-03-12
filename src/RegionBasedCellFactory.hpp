#ifndef REGIONBASEDCELLFACTORY_HPP_
#define REGIONBASEDCELLFACTORY_HPP_

#include "AbstractCardiacCellFactory.hpp"
#include "SimpleStimulus.hpp"
#include "MultiStimulus.hpp"
#include "AbstractCardiacCell.hpp"
#include "ChasteEllipsoid.hpp"



// A concrete version of AbstractCardiacCellFactory() which allows the user
// to set any number of 'stimulated regions', defined as a cuboid or sphere
// together with a time of stimulus. (Magnitude of stimulus is taken in in
// the constructor. Also, this class allows the user to set different initial
// conditions on the cell model.
template<class CELL, unsigned DIM>
class RegionBasedCellFactory : public AbstractCardiacCellFactory<DIM>
{
private:
    std::vector<double> mCellInitialConditions;
    std::vector<AbstractChasteRegion<DIM>*> mStimulatedRegions;
    std::vector<boost::shared_ptr<SimpleStimulus> > mStimuli;
    double mStimulusMagnitude;
    double mStimulusDuration;

public:

    RegionBasedCellFactory(double stimulusMagnitude, double stimulusDuration = 0.5)
        : AbstractCardiacCellFactory<DIM>(),
          mStimulusMagnitude(stimulusMagnitude),
          mStimulusDuration(stimulusDuration)
    {
    }

    virtual ~RegionBasedCellFactory()
    {
        for(unsigned i=0; i<mStimulatedRegions.size(); i++)
        {
            delete mStimulatedRegions[i];
        }
    }


    void AddStimulatedSphere(c_vector<double,DIM> sphereCentre, double radius, double stimulusTime)
    {
        ChastePoint<DIM> centre(sphereCentre);
        ChastePoint<DIM> radii(radius, radius, radius);
        ChasteEllipsoid<DIM>* p_sphere = new ChasteEllipsoid<DIM>(centre, radii);

        mStimulatedRegions.push_back(p_sphere);
        mStimuli.push_back(boost::shared_ptr<SimpleStimulus>(new SimpleStimulus(mStimulusMagnitude, mStimulusDuration, stimulusTime)));
    }

    void AddStimulatedCuboid(c_vector<double,DIM> pointA, c_vector<double,DIM> pointB, double stimulusTime)
    {
        ChastePoint<DIM> point_a(pointA);
        ChastePoint<DIM> point_b(pointB);
        ChasteCuboid<DIM>* p_cuboid = new ChasteCuboid<DIM>(point_a, point_b);
        mStimulatedRegions.push_back(p_cuboid);
        mStimuli.push_back(boost::shared_ptr<SimpleStimulus>(new SimpleStimulus(mStimulusMagnitude, mStimulusDuration, stimulusTime)));

    }

    void SetCellModelInitialConditions(std::vector<double> initialConditions)
    {
        mCellInitialConditions = initialConditions;
        CELL cell(this->mpSolver, this->mpZeroStimulus);
        if(cell.GetNumberOfStateVariables()!=initialConditions.size())
        {
            EXCEPTION("Initial conditions size doesn't match cell model");
        }
    }

    AbstractCardiacCell* CreateCardiacCellForTissueNode(unsigned node)
    {
        CELL* p_cell;
        c_vector<double,DIM> location = this->GetMesh()->GetNode(node)->rGetLocation();
        
        
        std::vector<bool> contained_in_regions(mStimulatedRegions.size(), false);
        bool contained_in_none = true;

        for(unsigned i=0; i<mStimulatedRegions.size(); i++)
        {
            if(mStimulatedRegions[i]->DoesContain(this->GetMesh()->GetNode(node)->rGetLocation()))
            {
                contained_in_regions[i] = true;
                contained_in_none = false;
            }
        }

        boost::shared_ptr<MultiStimulus> p_multi_stim;

        if(!contained_in_none)
        {
            p_multi_stim = boost::shared_ptr<MultiStimulus>(new MultiStimulus);
            for(unsigned i=0; i<contained_in_regions.size(); i++)
            {
                if(contained_in_regions[i])
                {
                    p_multi_stim->AddStimulus(mStimuli[i]);
                }
            }
            p_cell = new CELL(this->mpSolver, p_multi_stim);
        }
        else
        {
            p_cell = new CELL(this->mpSolver, this->mpZeroStimulus);
        }

        if(mCellInitialConditions.size()>0)
        {
            for(unsigned i=0; i<p_cell->GetNumberOfStateVariables(); i++)
            {
                p_cell->rGetStateVariables()[i] = mCellInitialConditions[i];
            }
        }

        return p_cell;
    }
};


#endif /* REGIONBASEDCELLFACTORY_HPP_ */
