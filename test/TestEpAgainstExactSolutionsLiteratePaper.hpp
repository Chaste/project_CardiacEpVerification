
#ifndef TESTEPAGAINSTEXACTSOLUTIONS_HPP_
#define TESTEPAGAINSTEXACTSOLUTIONS_HPP_


/* = Testing the cardiac solvers against exact solutions =
 *
 * This is the main code for solving the monodomain, bidomain and bidomain-with-bath model problems.
 *
 * The following are all standard includes: */
#include <cxxtest/TestSuite.h>
#include "AbstractCardiacCellFactory.hpp"
#include "DistributedTetrahedralMesh.hpp"
#include "EulerIvpOdeSolver.hpp"
#include "ZeroStimulus.hpp"
#include "PetscSetupAndFinalize.hpp"
/* Next, the includes for classes that are defined in this project.
 * `ModelProblemCellModel` defines the cell model for the model problems: */
#include "ModelProblemCellModel.hpp"
/* This file contains several classes which represent the exact solutions of the model problems: */
#include "ModelProblemExactSolutionClasses.hpp"
/* A class for calculating the L2 error (squared) at given time using a computed solution and exact solution: */
#include "L2ErrorSquaredCalculator.hpp"
/* A class for calculating the second part of the H1 error (squared) at given time using a computed solution and exact solution: */
#include "H1SemiNormErrorSquaredCalculator.hpp"
/* This next file defines three classes which inherit from `MonodomainProblem`, `BidomainProblem` and `BidomainWithBathProblem`
 * but also do error calculations (using the above).
 */
#include "CardiacProblemWithErrorCalculatorClasses.hpp"


/* We have to define cell factories to create the cell models as always. This cell factory
 * creates the custom `ModelProblemCellModel` for each node, gives each a zero stimulus,
 * passes in some parameters, and sets up the initial conditions used in the monodomain
 * and bidomain model problems.
 */
template<unsigned DIM>
class NonBathModelProblemCellFactory : public AbstractCardiacCellFactory<DIM>
{
private:
    double mBeta;
    unsigned mM1;
    unsigned mM2;
    unsigned mM3;

public:
    NonBathModelProblemCellFactory(double beta, unsigned m1, unsigned m2=0, unsigned m3 = 0)
        : AbstractCardiacCellFactory<DIM>(boost::shared_ptr<AbstractIvpOdeSolver>(new EulerIvpOdeSolver)),
          mBeta(beta),
          mM1(m1),
          mM2(m2),
          mM3(m3)
    {
    }

    virtual ~NonBathModelProblemCellFactory()
    {
    }

    AbstractCardiacCell* CreateCardiacCellForTissueNode(unsigned node)
    {
        ModelProblemCellModel* p_cell = new ModelProblemCellModel(this->mpSolver, this->mpZeroStimulus, mBeta);
        ChastePoint<DIM> point = this->GetMesh()->GetNode(node)->GetPoint();

        double F = TimeScaledFunctionF(point, 0.0, mM1, mM2, mM3); // just (1+t)^(1/2) * F(x,y,z), defined in ModelProblemExactSolutionClasses
        double G = FunctionG(point);                               // defined in ModelProblemExactSolutionClasses

        p_cell->rGetStateVariables()[0] = F;
        p_cell->rGetStateVariables()[1] = G + F;
        p_cell->rGetStateVariables()[2] = sqrt(1.0/G);

        return p_cell;
    }
};


/* This second cell factory is the same as the above except it passes in the initial conditions
 * defined in the bidomain-with-bath model problem.
 */
template<unsigned DIM>
class BathModelProblemCellFactory : public AbstractCardiacCellFactory<DIM>
{
private:
    double mBeta; // beta as defined in the paper
    unsigned mM1; // the m1 in the definition of F
    double mAlpha;// alpha as in the paper
    double mExtracellularConductivity; // the variable s_e in the paper

public:
    BathModelProblemCellFactory(double beta, unsigned m1, double alpha, double extracellularConductivity)
        : AbstractCardiacCellFactory<DIM>(boost::shared_ptr<AbstractIvpOdeSolver>(new EulerIvpOdeSolver)),
          mBeta(beta),
          mM1(m1),
          mAlpha(alpha),
          mExtracellularConductivity(extracellularConductivity)
    {
    }

    virtual ~BathModelProblemCellFactory()
    {
    }

    AbstractCardiacCell* CreateCardiacCellForTissueNode(unsigned node)
    {
        ModelProblemCellModel* p_cell = new ModelProblemCellModel(this->mpSolver, this->mpZeroStimulus, mBeta);
        ChastePoint<DIM> point = this->GetMesh()->GetNode(node)->GetPoint();

        double F = TimeScaledFunctionF(point, 0.0, mM1, 0.0, 0.0); // just (1+t)^(1/2) * F(x,y,z), defined in ModelProblemExactSolutionClasses
        double G = FunctionG(point);                               // defined in ModelProblemExactSolutionClasses
        double x = point[0];

        p_cell->rGetStateVariables()[0] = F - mAlpha*x/mExtracellularConductivity;
        p_cell->rGetStateVariables()[1] = G + F;
        p_cell->rGetStateVariables()[2] = sqrt(1.0/G);
        p_cell->rGetStateVariables()[3] = -mAlpha*x/mExtracellularConductivity;

        return p_cell;
    }
};



/*
 *  The code for running the model problems are defined in 'tests', as with all code in Chaste (see main documentation).
 */
class TestEpAgainstExactSolutionsLiteratePaper : public CxxTest::TestSuite
{
private:
    /* [The following function can be mostly ignored as an alternative is used]:
     *
     * A function for computing errors, taking in the output directory of a cardiac problem and a class
     * saying how to calculate the exact solution. Note: since the output directory must be written for
     * this function to be used, and with small dt this will lead to a huge datafile (if printing timestep
     * is dt), we no longer use this method for calculating errors. However the code has been kept as it is
     * convenient. The errors are now calculated by computing the contributions to the error at the end of
     * each timestep, as the simulation is progressing. See `MonodomainProblemWithErrorCalculator`
     * and related classes.
     */
    template<unsigned DIM,unsigned PROBLEM_DIM>
    void ComputeErrors(std::string outputDirectory, AbstractScalarFunction<DIM>* pExactSolution,
                       DistributedTetrahedralMesh<DIM,DIM>& rMesh, double printingTimestep,
                       std::string variable     /* Should be "V" or "Phi_e" */,
                       double& rReturnedLinfL2Error,
                       double& rReturnedL2H1Error,
                       bool testingMode = false /* if this is true, than the returned 'errors' are actually just the norms of the exact solution */)
    {
        if(variable != "V" && variable !="Phi_e")
        {
            EXCEPTION("Bad input");
        }

        // get the data of the data file
        Hdf5DataReader reader(outputDirectory,"results");
        unsigned num_timesteps = reader.GetUnlimitedDimensionValues().size();

        DistributedVectorFactory factory(rMesh.GetNumNodes());
        Vec solution = factory.CreateVec();

        rReturnedLinfL2Error = 0.0;
        rReturnedL2H1Error = 0.0;

        // The phi_e written to file for initial condition is just zeros, so needs to be ignored
        // (Really, the solver should solve the second bidomain equation given initial V to determine what
        // initial phi_e is, but the bidomain solver doesn't do this, and phi_e is just initialised to zero)
        unsigned first = (variable == "V" ? 0 : 1);

        for(unsigned timestep=first; timestep<num_timesteps; timestep++)
        {
            reader.GetVariableOverNodes(solution, variable, timestep);
            if(testingMode)
            {
                PetscVecTools::Zero(solution);
            }

            double time = printingTimestep*timestep;

            L2ErrorSquaredCalculator<DIM> l2_calc(pExactSolution,time,(variable=="V"));
            double l2_error_sqd = l2_calc.Calculate(rMesh,solution);
            double l2_error = sqrt(l2_error_sqd);

            H1SemiNormErrorSquaredCalculator<DIM> h1_calc(pExactSolution,time,(variable=="V"));
            double h1_seminorm_error_sqd = h1_calc.Calculate(rMesh,solution);
            double h1_error = sqrt(l2_error_sqd + h1_seminorm_error_sqd);

            rReturnedLinfL2Error = std::max(rReturnedLinfL2Error, l2_error);

            double factor = (timestep==0 || timestep+1==num_timesteps ? 0.5 : 1.0);
            rReturnedL2H1Error += h1_error*factor*printingTimestep;
        }
    }


    /*
     * == Monodomain model problem ==
     *
     * This function is the code for running the monodomain model problem, which we will walk through.
     */
    template<unsigned DIM>
    void RunMonodomainProblem(double parametersScaleFactor /*how much to scale h and dt*/, bool doTest=false /*see later*/)
    {
        /* Define h and dt based on the `parametersScaleFactor`. Note dt is proportional to h^2^ as required. */
        double init_h = 0.1;
        double h      = init_h*parametersScaleFactor; // note: everything dimensionless
        double dt_ode = 0.1*parametersScaleFactor*parametersScaleFactor;
        double dt_pde = 0.1*parametersScaleFactor*parametersScaleFactor;
        double dt_printing = dt_pde;

        /* Define conductivities in each direction as specified in the paper. */
        double s1 = 1.1/(M_PI*M_PI);
        double s2 = 1.2/(M_PI*M_PI);
        double s3 = 0.3/(M_PI*M_PI);

        /* Define the integers m,,1,,, m,,2,, m,,3,, that go into F, as specified in the paper. */
        unsigned m1 = 1;
        unsigned m2 = 2;
        unsigned m3 = 3;

        /* Set up the mesh to be the unit line/square/cube. Also compute the constant c. */
        DistributedTetrahedralMesh<DIM,DIM> mesh;
        double c;

        if(DIM==1)
        {
            c = -s1*m1*m1*M_PI*M_PI;
            mesh.ConstructRegularSlabMesh(h, 1.0);
        }
        else if(DIM==2)
        {
            c = - (s1*m1*m1*M_PI*M_PI + s2*m2*m2*M_PI*M_PI);
            mesh.ConstructRegularSlabMesh(h, 1.0, 1.0);
        }
        else
        {
            c = - (s1*m1*m1*M_PI*M_PI + s2*m2*m2*M_PI*M_PI + s3*m3*m3*M_PI*M_PI);
            mesh.ConstructRegularSlabMesh(h, 1.0, 1.0, 1.0);
        }

        /* End time is T = 1. */
        double end_time = 1.0;
        HeartConfig::Instance()->SetSimulationDuration(end_time);

        /* Define an output directory, but see comments below. */
        std::stringstream output_dir;
        output_dir << "MonodomainExactSolution_" << DIM << "D_" << parametersScaleFactor;
        HeartConfig::Instance()->SetOutputDirectory(output_dir.str());
        HeartConfig::Instance()->SetOutputFilenamePrefix("results");

        /* Set the timesteps, define the conductivity tensor, the capacitance and the surface-to-volume ratio. */
        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(dt_ode, dt_pde, dt_printing);
        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(s1,s2,s3));
        HeartConfig::Instance()->SetCapacitance(2.0);
        HeartConfig::Instance()->SetSurfaceAreaToVolumeRatio(3.0);

        /* The following cell factory creates cell models for the model problem cell model defined in the paper, and
         * with the required initial values of (V,u,,1,,,u,,2,,,u,,3,,).
         */
        NonBathModelProblemCellFactory<DIM> cell_factory(c, m1, m2, m3);

        /* This class can be used to get the exact solution for the voltage: (1+t)^1/2^F(x). */
        VoltageExactSolution<DIM> voltage_soln(m1,m2,m3);

        /* Define the monodomain problem class. `MonodomainProblemWithErrorCalculator` (defined in this project) is just
         * `MonodomainProblem` but also does error computing at the end of each timestep.
         */
        MonodomainProblemWithErrorCalculator<DIM> monodomain_problem( &cell_factory, &voltage_soln );

        /* Set the mesh, '''don't''' write output unless in testing mode (since printing_dt = pde_dt, a lot of output would be
         * written to file (printing dt is small so that error computations can be carried out each dt), initialise
         * and solve.
         */
        monodomain_problem.SetMesh(&mesh);
        monodomain_problem.PrintOutput(doTest);

        //monodomain_problem.SetWriteInfo();
        monodomain_problem.Initialise();
        monodomain_problem.Solve();

        /* Print the errors to screen (the commas and semi-colon are for easy copy & paste into Matlab). */
        std::cout << std::setprecision(10);
        std::cout << h << ", " << monodomain_problem.mVoltageLinfL2Error << ", " << monodomain_problem.mVoltageL2H1Error << ";\n";


        /* There is also some test code. To check that the code (in particular the error calculators) is doing what it should be, if DIM=1 and the coarsest mesh
         * is being used, we have the following tests.
         */
        if(doTest)
        {
            if(DIM!=1 || parametersScaleFactor!=1.0)
            {
                EXCEPTION("Test mode is only for 1d with factor=1");
            }

            /* Check nothing has changed: */
            TS_ASSERT_DELTA(monodomain_problem.mVoltageLinfL2Error, 0.0472926, 1e-5);
            TS_ASSERT_DELTA(monodomain_problem.mVoltageL2H1Error, 0.255933, 1e-5);

            /* Check the two ways of computing the error give the same results: */
            double l_inf_l2;
            double l2_h1;
            ComputeErrors<DIM,1>(output_dir.str(), &voltage_soln, mesh, dt_printing, "V", l_inf_l2, l2_h1);
            TS_ASSERT_DELTA(l_inf_l2, monodomain_problem.mVoltageLinfL2Error, 1e-6);
            TS_ASSERT_DELTA(l2_h1, monodomain_problem.mVoltageL2H1Error, 1e-6);

            /* Finally, check the second way of computing the error: when `ComputeErrors()` is called with `true` as the last parameter it
             * ignores the numerical solution, ie just calculates the norm of the exact solution, which we can also
             * calculate analytically.
             */
            double linf_l2_norm_V;
            double l2_h1_norm_V;
            ComputeErrors<DIM,1>(output_dir.str(), &voltage_soln, mesh,  dt_printing, "V", linf_l2_norm_V, l2_h1_norm_V, true);
            TS_ASSERT_DELTA(linf_l2_norm_V, 1.0, 1e-4); // ||V(t)||_2 = sqrt ((1+t)/2) so max = 1
            TS_ASSERT_DELTA(l2_h1_norm_V, 2.855, 1e-1); // sqrt (3(1+pi^2)/2)
        }
    }

    /* == Bidomain model problem ==
     *
     *
     * Next, the main function for solving the bidomain model problem. This is basically the same as the monodomain code, except has
     * an extracellular conductivity, and gets the errors for both voltage and extracellular potential. */
    template<unsigned DIM>
    void RunBidomainProblem(double parametersScaleFactor, bool doTest=false)
    {
    	/* Code same as monodomain version: */
        double init_h = 0.1;
        double h      = init_h*parametersScaleFactor; // note: everything dimensionless
        double dt_ode = 0.1*parametersScaleFactor*parametersScaleFactor;
        double dt_pde = 0.1*parametersScaleFactor*parametersScaleFactor;
        double dt_printing = dt_pde;

        double s1 = 1.1/(M_PI*M_PI);
        double s2 = 1.2/(M_PI*M_PI);
        double s3 = 0.3/(M_PI*M_PI);

        unsigned m1 = 1.0;
        unsigned m2 = 2.0;
        unsigned m3 = 3.0;

        DistributedTetrahedralMesh<DIM,DIM> mesh;
        double c;

        if(DIM==1)
        {
            c = -s1*m1*m1*M_PI*M_PI;
            mesh.ConstructRegularSlabMesh(h, 1.0);
        }
        else if(DIM==2)
        {
            c = - (s1*m1*m1*M_PI*M_PI + s2*m2*m2*M_PI*M_PI);
            mesh.ConstructRegularSlabMesh(h, 1.0, 1.0);
        }
        else
        {
            c = - (s1*m1*m1*M_PI*M_PI + s2*m2*m2*M_PI*M_PI + s3*m3*m3*M_PI*M_PI);
            mesh.ConstructRegularSlabMesh(h, 1.0, 1.0, 1.0);
        }

        /* Define the constant k: */
        double k = 1.0/sqrt(2);
        double sigma_e_factor = (1.0-k)/k;

    	/* Code similar to monodomain version: */
        double end_time = 1.0;
        HeartConfig::Instance()->SetSimulationDuration(end_time);

        std::stringstream output_dir;
        output_dir << "BidomainExactSolution_" << DIM << "D_" << parametersScaleFactor;
        HeartConfig::Instance()->SetOutputDirectory(output_dir.str());
        HeartConfig::Instance()->SetOutputFilenamePrefix("results");

        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(dt_ode, dt_pde, dt_printing);

        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(s1,s2,s3));
        HeartConfig::Instance()->SetExtracellularConductivities(Create_c_vector(s1*sigma_e_factor,s2*sigma_e_factor,s3*sigma_e_factor));
        HeartConfig::Instance()->SetCapacitance(2.0);
        HeartConfig::Instance()->SetSurfaceAreaToVolumeRatio(3.0);

        NonBathModelProblemCellFactory<DIM> cell_factory(c*(1-k), m1, m2, m3);

        /* Classes for returning V=(1+t)^1/2^F(x) and phi_e = -k(1+t)^1/2^F(x) and their derivatives: */
        VoltageExactSolution<DIM> voltage_soln(m1,m2,m3);
        ExtracellularPotentialExactSolution<DIM> phi_e_soln(k,m1,m2,m3);

        /* Solve and print errors. */
        BidomainProblemWithErrorCalculator<DIM> bidomain_problem( &cell_factory, &voltage_soln, &phi_e_soln );
        bidomain_problem.PrintOutput(doTest);

        bidomain_problem.SetMesh(&mesh);
        //bidomain_problem.SetWriteInfo();
        bidomain_problem.Initialise();
        bidomain_problem.Solve();

        std::cout << std::setprecision(10);
        std::cout << h << ", " << bidomain_problem.mVoltageLinfL2Error << ", " << bidomain_problem.mVoltageL2H1Error << ", " << bidomain_problem.mExtracellularPotentialLinfL2Error << ", " << bidomain_problem.mExtracellularPotentialL2H1Error << ";\n";

        /* Test code similar to monodomain version: */
        if(doTest)
        {
            if(DIM!=1 || parametersScaleFactor!=1.0)
            {
                EXCEPTION("Test mode is only for 1d with factor=1");
            }

            TS_ASSERT_DELTA(bidomain_problem.mVoltageLinfL2Error, 0.0426506, 1e-5);
            TS_ASSERT_DELTA(bidomain_problem.mVoltageL2H1Error, 0.253785, 1e-5);
            TS_ASSERT_DELTA(bidomain_problem.mExtracellularPotentialLinfL2Error, 0.0298811, 1e-5);
            TS_ASSERT_DELTA(bidomain_problem.mExtracellularPotentialL2H1Error, 0.172321, 1e-5);

            double l_inf_l2_V;
            double l2_h1_V;
            double l_inf_l2_phi_e;
            double l2_h1_phi_e;
            ComputeErrors<DIM,2>(output_dir.str(), &voltage_soln, mesh, dt_printing, "V", l_inf_l2_V, l2_h1_V);
            ComputeErrors<DIM,2>(output_dir.str(), &phi_e_soln,   mesh, dt_printing, "Phi_e", l_inf_l2_phi_e, l2_h1_phi_e);
            TS_ASSERT_DELTA(l_inf_l2_V, bidomain_problem.mVoltageLinfL2Error, 1e-6);
            TS_ASSERT_DELTA(l2_h1_V, bidomain_problem.mVoltageL2H1Error, 1e-6);
            TS_ASSERT_DELTA(l_inf_l2_phi_e, bidomain_problem.mExtracellularPotentialLinfL2Error, 1e-6);
            TS_ASSERT_DELTA(l2_h1_phi_e, bidomain_problem.mExtracellularPotentialL2H1Error, 1e-6);

            double linf_l2_norm_V;
            double l2_h1_norm_V;
            double linf_l2_norm_phi_e;
            double l2_h1_norm_phi_e;
            ComputeErrors<DIM,1>(output_dir.str(), &voltage_soln, mesh,  dt_printing, "V", linf_l2_norm_V, l2_h1_norm_V, true);
            TS_ASSERT_DELTA(linf_l2_norm_V, 1.0, 1e-4); // ||V(t)||_2 = sqrt ((1+t)/2) so max = 1
            TS_ASSERT_DELTA(l2_h1_norm_V, 2.855, 1e-1); // sqrt (3(1+pi^2)/2)

            ComputeErrors<DIM,1>(output_dir.str(), &phi_e_soln, mesh,  dt_printing, "Phi_e", linf_l2_norm_phi_e, l2_h1_norm_phi_e, true);
            TS_ASSERT_DELTA(linf_l2_norm_phi_e, 1.0/sqrt(2), 1e-4); // phi_e = -k*V
            TS_ASSERT_DELTA(l2_h1_norm_phi_e, 2.855/sqrt(2), 1e-1); //
        }
    }

    /* == Bidomain-with-bath model problem ==
     *
     * Finally, the bidomain-with-bath-model problem: */
    template<unsigned DIM>
    void RunBidomainWithBathProblem(double parametersScaleFactor, bool doTest=false)
    {
        /* All this is as before: */
        double init_h = 0.1;
        double h      = init_h*parametersScaleFactor; // dimensionless
        double dt_ode = 0.1*parametersScaleFactor*parametersScaleFactor;
        double dt_pde = 0.1*parametersScaleFactor*parametersScaleFactor;
        double dt_printing = dt_pde;

        double s1 = 1.1/(M_PI*M_PI);
        double s2 = 1.2/(M_PI*M_PI);
        double s3 = 0.3/(M_PI*M_PI);

        /* Here, F(x,y,z) = cos(pi*x), whatever the dimension. */
        unsigned m1 = 1;

        DistributedTetrahedralMesh<DIM,DIM> mesh;
        double c = -s1*m1*m1*M_PI*M_PI;

        /* Set up the domain: x in `[-1,2]`; y,z in `[0,1]`. Note the translation at the end. */
        c_vector<double,DIM> disp = zero_vector<double>(DIM);
        disp(0) = -1.0;
        if(DIM==1)
        {
            mesh.ConstructRegularSlabMesh(h, 3.0);
        }
        else if(DIM==2)
        {
            mesh.ConstructRegularSlabMesh(h, 3.0, 1.0);
        }
        else
        {
            mesh.ConstructRegularSlabMesh(h, 3.0, 1.0, 1.0);
        }
        mesh.Translate(disp);

        /* Set appropriate elements as bath: */
        for(unsigned i=0; i<mesh.GetNumElements(); i++)
        {
            double x = mesh.GetElement(i)->CalculateCentroid()[0];
            if( (x<0) || (x>1) )
            {
                mesh.GetElement(i)->SetAttribute(HeartRegionCode::GetValidBathId());
            }
        }

        /* As before: */
        double k = 1.0/sqrt(2);
        double sigma_e_factor = (1.0-k)/k;

        /* Set up the bath conductivity and the electrodes, I,,E,,=-alpha on x=-1, and I,,E,,=alpha on x=2: */
        double alpha = 0.01;
        double extracellular_conductivity = s1*sigma_e_factor;
        double bath_conductivity = extracellular_conductivity/2.0;
        // the following says provide a stimulus of -alpha on the x=min surface (the second parameter '0' indicates x).
        // The code then computes that the stimulus on the opposite surface should be alpha for conservation of current. The false says no ground electrode
        HeartConfig::Instance()->SetElectrodeParameters(false, 0, -alpha, -1.0/*switch on time*/, 1000/*switch off time*/);
        HeartConfig::Instance()->SetBathConductivity(bath_conductivity);

        /* As before: */
        double end_time = 1.0;
        HeartConfig::Instance()->SetSimulationDuration(end_time);

        std::stringstream output_dir;
        output_dir << "BidomainWithBathExactSolution_" << DIM << "D_" << parametersScaleFactor;
        HeartConfig::Instance()->SetOutputDirectory(output_dir.str());
        HeartConfig::Instance()->SetOutputFilenamePrefix("results");
        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(dt_ode, dt_pde, dt_printing);

        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(s1,s2,s3));
        HeartConfig::Instance()->SetExtracellularConductivities(Create_c_vector(s1*sigma_e_factor,s2*sigma_e_factor,s3*sigma_e_factor));
        HeartConfig::Instance()->SetCapacitance(2.0);
        HeartConfig::Instance()->SetSurfaceAreaToVolumeRatio(3.0);

        /* Use the bath version of the cell factory (different initial conditions): */
        BathModelProblemCellFactory<DIM> cell_factory(c*(1-k), m1, alpha, extracellular_conductivity);

        /* Solve and output errors: */
        VoltageExactSolutionBath<DIM> voltage_soln(m1,alpha,extracellular_conductivity);
        ExtracellularPotentialExactSolutionBath<DIM> phi_e_soln(m1,k,alpha,extracellular_conductivity,bath_conductivity);

        BidomainWithBathProblemWithErrorCalculator<DIM> bidomain_problem( &cell_factory, &voltage_soln, &phi_e_soln);
        bidomain_problem.PrintOutput(doTest);

        bidomain_problem.SetMesh(&mesh);
        //bidomain_problem.SetWriteInfo();
        bidomain_problem.Initialise();
        bidomain_problem.Solve();

        std::cout << std::setprecision(10);
        std::cout << h << ", " << bidomain_problem.mVoltageLinfL2Error << ", " << bidomain_problem.mVoltageL2H1Error  << ", " << bidomain_problem.mExtracellularPotentialLinfL2Error << ", " << bidomain_problem.mExtracellularPotentialL2H1Error << ";\n";


        /* Similar to before: */
        if(doTest)
        {
            if(DIM!=1 || parametersScaleFactor!=1.0)
            {
               EXCEPTION("Test mode is only for 1d with factor=1");
            }
            TS_ASSERT_DELTA(bidomain_problem.mVoltageLinfL2Error, 0.0426506, 1e-4);
            TS_ASSERT_DELTA(bidomain_problem.mVoltageL2H1Error, 0.253785, 1e-4);
            TS_ASSERT_DELTA(bidomain_problem.mExtracellularPotentialLinfL2Error, 0.0562449, 1e-3);
            TS_ASSERT_DELTA(bidomain_problem.mExtracellularPotentialL2H1Error, 0.174176, 1e-3);

            double l_inf_l2_V;
            double l2_h1_V;
            double l_inf_l2_phi_e;
            double l2_h1_phi_e;
            ComputeErrors<DIM,2>(output_dir.str(), &voltage_soln, mesh, dt_printing, "V", l_inf_l2_V, l2_h1_V);
            ComputeErrors<DIM,2>(output_dir.str(), &phi_e_soln,   mesh, dt_printing, "Phi_e", l_inf_l2_phi_e, l2_h1_phi_e);
            TS_ASSERT_DELTA(l_inf_l2_V, bidomain_problem.mVoltageLinfL2Error, 1e-6);
            TS_ASSERT_DELTA(l2_h1_V, bidomain_problem.mVoltageL2H1Error, 1e-6);
            TS_ASSERT_DELTA(l_inf_l2_phi_e, bidomain_problem.mExtracellularPotentialLinfL2Error, 1e-6);
            TS_ASSERT_DELTA(l2_h1_phi_e, bidomain_problem.mExtracellularPotentialL2H1Error, 1e-6);
        }
    }

    /* == Main test ==
     *
     * Finally, we have the public 'tests', which actually run the simulations. Note that as integrals over the domain have to be
     * computed each timestep, virtually all the time is spent in the error calculations, certainly for the larger N (ie smaller h). Also,
     * dt is proportional to h^2, so for 3D, the work increases by at least a factor of 32 (2^3^ 2^2^) as N doubles. N=3 for
     * 3D monodomain and bidomain were run on a cluster.
     */
public:
    void TestRunTests() throw (Exception)
    {
        RunMonodomainProblem<1>(1.0,true);
        RunBidomainProblem<1>(1.0,true);
        RunBidomainWithBathProblem<1>(1.0,true);
    }

    void doneTestMonodomain1d() throw (Exception)
    {
        for(unsigned N=0; N<5; N++)
        {
            double factor = 1.0/pow(2,N);
            RunMonodomainProblem<1>(factor);
        }
    }

    void doneTestMonodomain2d() throw (Exception)
    {
        for(unsigned N=0; N<5; N++)
        {
            double factor = 1.0/pow(2,N);
            RunMonodomainProblem<2>(factor);
        }
    }


    void doneTestMonodomain3d() throw (Exception)
    {
        HeartConfig::Instance()->SetUseAbsoluteTolerance(1e-9);
        for(unsigned N=0; N<4; N++) // note: N=3 takes very long
        {
            double factor = 1.0/pow(2,N);
            RunMonodomainProblem<3>(factor);
        }
    }

    void doneTestBidomain1d() throw (Exception)
    {
        HeartConfig::Instance()->SetUseAbsoluteTolerance(1e-9);
        for(unsigned N=0; N<5; N++)
        {
            double factor = 1.0/pow(2,N);
            RunBidomainProblem<1>(factor);
        }
    }

    void doneTestBidomain2d() throw (Exception)
    {
        HeartConfig::Instance()->SetUseAbsoluteTolerance(1e-9);
        for(unsigned N=0; N<5; N++)
        {
            double factor = 1.0/pow(2,N);
            RunBidomainProblem<2>(factor);
        }
    }

    void todoTestBidomain3d() throw (Exception)
    {
        HeartConfig::Instance()->SetUseAbsoluteTolerance(1e-11);
        for(unsigned N=0; N<4; N++) // note: N=3 takes very long
        {
            double factor = 1.0/pow(2,N);
            RunBidomainProblem<3>(factor);
        }
    }

    void doneTestBidomainWithBath2d() throw (Exception)
    {
        HeartConfig::Instance()->SetUseAbsoluteTolerance(1e-9);
        for(unsigned N=0; N<5; N++)
        {
            double factor = 1.0/pow(2,N);
            RunBidomainWithBathProblem<2>(factor);
        }
    }


};

#endif //TESTEPAGAINSTEXACTSOLUTIONS_HPP_

/*
Entering TestBidomain3d
0.1, 0.1167727055, 1.99370425, 0.08253275914, 1.350677342;
0.05, 0.03257310981, 1.037230142, 0.02301724628, 0.7258621174;
0.025, 0.008372181866, 0.5243032404, 0.005915685911, 0.3697861914;

*/

