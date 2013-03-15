
#ifndef TESTCONDUCTIONVELOCITYCASESTUDYLITERATEPAPER_HPP_
#define TESTCONDUCTIONVELOCITYCASESTUDYLITERATEPAPER_HPP_

/*
 *  = Determining the accuracy of the conduction velocity in 1D =
 *
 *  This is the main code for the first calculation verification case study. The simulation just involves a simple
 *  extension of a standard monodomain simulation.
 *
 *  The following are all standard includes:
 */
#include <cxxtest/TestSuite.h>
#include "MonodomainProblem.hpp"
#include "DistributedTetrahedralMesh.hpp"
#include "AbstractCardiacCellFactory.hpp"
#include "LuoRudy1991.hpp"
#include "PetscSetupAndFinalize.hpp"

/* Define a simple cell factory which creates Luo-Rudy cells, and stimulates the given region.
 */
class SimpleCellFactory : public AbstractCardiacCellFactory<1>
{
private:
    boost::shared_ptr<SimpleStimulus> mpStimulus;
    double mStimWidth;

public:
    SimpleCellFactory(double stimWidth)
        : AbstractCardiacCellFactory<1>(),
          mpStimulus(new SimpleStimulus(-5e5, 0.5)),
          mStimWidth(stimWidth)
    {
    }

    AbstractCardiacCell* CreateCardiacCellForTissueNode(unsigned nodeIndex)
    {
        double x = this->GetMesh()->GetNode(nodeIndex)->rGetLocation()[0];
        if (x<=mStimWidth)
        {
            return new CellLuoRudy1991FromCellML(mpSolver, mpStimulus);
        }
        else
        {
            return new CellLuoRudy1991FromCellML(mpSolver, mpZeroStimulus);
        }
    }
};

/* This class inherits from `MonodomainProblem` but does some extra work at the end of every (printing) timestep.
 * (Note: the printing timestep will be set to be the same as the pde timestep).
 */
class MonodomainProblemWithCvComputer1d : public MonodomainProblem<1>
{
private:
    double mQuarterNodeActivationTime;          // to be computed: activation time of node at x=0.25
    double mThreeQuarterNodeActivationTime;     // to be computed: activation time of node at x=0.75
    ReplicatableVector* mpVoltageLastTimestep;  // We will need to save the voltages at the last timestep
    double mLastTime;

public:
    MonodomainProblemWithCvComputer1d(AbstractCardiacCellFactory<1>* pCellFactory)
        : MonodomainProblem<1>(pCellFactory),
          mQuarterNodeActivationTime(DBL_MAX),
          mThreeQuarterNodeActivationTime(DBL_MAX),
          mpVoltageLastTimestep(NULL),
          mLastTime(DBL_MAX)
    {
    }

    ~MonodomainProblemWithCvComputer1d()
    {
        delete mpVoltageLastTimestep;
    }

    /* At the end of every timestep, determine if the voltage for the nodes at x=0.25 and x=0.75 have just become positive. If so,
     * use the current value of the voltage and the last value of the voltage, and linear interpolation over time, to get the time the
     * voltage became positive.
     */
    void OnEndOfTimestep(double time)
    {
        unsigned quarter_index = (this->mpMesh->GetNumNodes()-1)/4;
        unsigned three_quarter_index = 3*(this->mpMesh->GetNumNodes()-1)/4;

        double width = this->mpMesh->GetWidth(0/*ie x-direction*/);

        // a quick check that the nodes are in the expected places
        if((fabs(this->mpMesh->GetNode(quarter_index)->rGetLocation()[0]-0.25*width)>1e-12) || (fabs(this->mpMesh->GetNode(three_quarter_index)->rGetLocation()[0]-0.75*width)>1e-12))
        {
            std::cout << this->mpMesh->GetNode(quarter_index)->rGetLocation()[0] << " " << 0.25*width << " "
                      << this->mpMesh->GetNode(three_quarter_index)->rGetLocation()[0] << "\n";
            NEVER_REACHED;
        }

        ReplicatableVector voltage_repl(this->mSolution);
        if(mpVoltageLastTimestep==NULL)
        {
            mpVoltageLastTimestep = new ReplicatableVector(this->mSolution);
        }

        // determine activation times
        if(mQuarterNodeActivationTime==DBL_MAX && voltage_repl[quarter_index]>0)
        {
            double V2 = voltage_repl[quarter_index];
            double V1 = (*mpVoltageLastTimestep)[quarter_index];
            double t2 = time;
            double t1 = mLastTime;

            mQuarterNodeActivationTime = (V2*t1 - V1*t2)/(V2-V1);
        }

        if(mThreeQuarterNodeActivationTime==DBL_MAX && voltage_repl[three_quarter_index]>0)
        {
            double V2 = voltage_repl[three_quarter_index];
            double V1 = (*mpVoltageLastTimestep)[three_quarter_index];
            double t2 = time;
            double t1 = mLastTime;

            mThreeQuarterNodeActivationTime = (V2*t1 - V1*t2)/(V2-V1);
        }

        mLastTime = time;
    }

    /* Get the conduction velocity from the activation times. */
    double GetConductionVelocity()
    {
        unsigned quarter_index = (this->mpMesh->GetNumNodes()-1)/4;
        unsigned three_quarter_index = 3*(this->mpMesh->GetNumNodes()-1)/4;

        if( mQuarterNodeActivationTime==DBL_MAX || mThreeQuarterNodeActivationTime==DBL_MAX)
        {
            EXCEPTION("One of the nodes was not stimulated");
        }

        double dx = this->mpMesh->GetNode(three_quarter_index)->rGetLocation()[0] - this->mpMesh->GetNode(quarter_index)->rGetLocation()[0];
        return dx/(mThreeQuarterNodeActivationTime-mQuarterNodeActivationTime);
    }
};




class TestConductionVelocityCaseStudyLiteratePaper : public CxxTest::TestSuite
{
private:
    /* The main simulation function: */
    void Run(double parametersScaleFactor/*how much to scale h and dt*/, bool doTest=false /*see later*/)
    {
        /* Define some initial parameters: */
        double width = 1.0;      //cm
        double stim_width = 0.1; //cm
        double end_time = 10.0;  //ms

        /* Define h and dt. Note h is proportional to dt. */
        double init_h = 0.05; // cm, ie 500 um
        double h  = init_h*parametersScaleFactor;
        double dt = 0.01*parametersScaleFactor;

        double dt_ode = dt;
        double dt_pde = dt;
        double printing_dt = dt;

        /* Run a standard monodomain simulation, except use our class `MonodomainProblemWithCvComputer1d`: */
        DistributedTetrahedralMesh<1,1> mesh;
        mesh.ConstructRegularSlabMesh(h, width);

        HeartConfig::Instance()->SetSimulationDuration(end_time); //ms

        std::stringstream output_dir;
        output_dir << "CalculationVerification1d_" << parametersScaleFactor;
        HeartConfig::Instance()->SetOutputDirectory(output_dir.str());
        HeartConfig::Instance()->SetOutputFilenamePrefix("results");
        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(dt_ode, dt_pde, printing_dt);

        SimpleCellFactory cell_factory(stim_width);

        MonodomainProblemWithCvComputer1d monodomain_problem( &cell_factory );
        monodomain_problem.SetMesh(&mesh);

        // only output at one node
        std::vector<unsigned> nodes_to_be_output;
        unsigned three_quarter_index = 3*(mesh.GetNumNodes()-1)/4;
        nodes_to_be_output.push_back(three_quarter_index);
        monodomain_problem.SetOutputNodes(nodes_to_be_output);

        //monodomain_problem.SetWriteInfo();
        monodomain_problem.Initialise();
        monodomain_problem.Solve();

        /* Print results: */
        std::cout << std::setprecision(9);
        double cv = monodomain_problem.GetConductionVelocity();
        std::cout << h << ", " << dt << ", " << cv << "\n";


        /* If in 'testing mode', which only applies if the coarsest mesh is being used, we do a quick
         * test that nothing has changed:
         */
        if(doTest)
        {
            if(parametersScaleFactor!=1.0)
            {
               EXCEPTION("Test with factor=1");
            }
            TS_ASSERT_DELTA(cv, 0.0736377, 1e-4);
        }
    }

/* The code which runs the above. The Richardson extrapolation of the results in done outside of Chaste. See the folder named 'other' in
 * this project of a text file of the results. */
public:
    void TestRunTest() throw (Exception)
    {
        Run(1.0, true);
    }

    void TestRun1dCv() throw (Exception)
    {
        unsigned num_sims = 5; // change to 9 for all paper results

        HeartConfig::Instance()->SetUseAbsoluteTolerance(1e-9);
        for(unsigned N=0; N<num_sims; N++)
        {
            double factor = 1.0/pow(2,N);
            Run(factor);
        }
    }
};




#endif //TESTCONDUCTIONVELOCITYCASESTUDYLITERATEPAPER_HPP_


