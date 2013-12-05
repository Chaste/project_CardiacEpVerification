#ifndef TESTREENTRYONRABBITMESHLITERATEPAPER_HPP_
#define TESTREENTRYONRABBITMESHLITERATEPAPER_HPP_

/* == Reentry simulations on realistic rabbit geometry ==
 *
 * This file provides the code used to run the simulations on the realistic rabbit geometry in the second calculation
 * verification case study.
 *
 * See main cardiac tutorials for more detailed descriptions of cardiac simulations.
 *
 */

/* First, we have some standard includes. */
#include <cxxtest/TestSuite.h>
#include <boost/assign.hpp>
#include "CardiacSimulationArchiver.hpp"
#include "MonodomainProblem.hpp"
#include "PetscSetupAndFinalize.hpp"
/* The cell model we use is the Mahajan2008 cell model, which is included with Chaste. */
#include "Mahajan2008BackwardEuler.hpp"
/* The following cell factory is defined in this project. It allows the user to specify
 * stimuli for chosen spheres and cuboids.
 */
#include "RegionBasedCellFactory.hpp"

/* A simple enumeration for which mesh to use: */
typedef enum GeometryOption_
{
    COARSERES_ISOTROPIC = 0,
    MEDIUMRES_ISOTROPIC,
    FULLRES_ISOTROPIC
} GeometryOption;


/* This cell factory inherits from the `RegionBasedCellFactory` but adds a stimulated spherical
 * region at the apex of this geometry.
 */
template<class CELL>
class RegionBasedCellFactoryWithApexS1 : public RegionBasedCellFactory<CELL,3>
{
public:
    RegionBasedCellFactoryWithApexS1()
        : RegionBasedCellFactory<CELL,3>(-500000.0/*stimulus magnitude*/)
    {
        c_vector<double,3> apex_region_centre;
        apex_region_centre(0) = 0.3952;
        apex_region_centre(1) = 1.1142;
        apex_region_centre(2) = 0.2093;
        double apex_region_radius = 0.075;

        this->AddStimulatedSphere(apex_region_centre, apex_region_radius, 0.0 /*stim time*/);
    }
};


/* The main test class: */
class TestReentryOnRabbitMeshLiteratePaper : public CxxTest::TestSuite
{
private:
    /* This method sets some numerical options: */
    void SetHeartConfigForTest()
    {
        HeartConfig::Instance()->Reset();
        HeartConfig::Instance()->SetKSPSolver("cg");
        HeartConfig::Instance()->SetKSPPreconditioner("bjacobi");
        HeartConfig::Instance()->SetUseAbsoluteTolerance(1e-5);
    }

public:
    /* The main simulation method: */
    void TestReentryS1andS2() throw(Exception)
    {
        /* The parts that the user can easily change are all listed here: */
        GeometryOption geometry = COARSERES_ISOTROPIC;  // Other options found in enumeration above
        double s2_time = 170;                           // Time of S2 stimulus in ms. Set this to 'DBL_MAX' for there to be no S2.
        double end_time = 1000;                         // in ms
        double printing_time = 10.0;                    // in ms
        bool write_archive = false;                     // Whether to write an archive at the end of the simulation (so can simulation can be reloaded and run).
        std::string notes = "";                         // Anything here is added to output directory name (see below).


        /* Initial set up: */
        SetHeartConfigForTest();
        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.1 /*ode dt*/, 0.1/*pde dt*/, printing_time);
        HeartConfig::Instance()->SetSimulationDuration(end_time);

        /* The conductivity used is 0.9333, chosen arbitrarily so that reentry is sustained for several rotations. This number is
         * the monodomain conductivity corresponding the Clerc 1976 intra- and extra-cellular conductivities, scaled by (almost exactly)
         * 70%.
         */
        double conductivity = 14.0/15.0;
        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(conductivity,conductivity,conductivity));

        /* '''The coarse mesh is in the repository, but there others are not - you will need to obtain the medium/fine meshes and provide their locations here'''
         * (see project intro page for details on where to get the meshes). */
        switch(geometry)
        {
            case COARSERES_ISOTROPIC:
            {
             	HeartConfig::Instance()->SetMeshFileName("apps/texttest/weekly/Propagation3d/OxfordRabbitHeart_482um");
                break;
            }
            case MEDIUMRES_ISOTROPIC:
            {
                // ** User should provide correct location **
                HeartConfig::Instance()->SetMeshFileName("../Data/OxfordRabbitHeart/medium_res/OxfordRabbitHeart_251um");
                break;
            }
            case FULLRES_ISOTROPIC:
            {
                // ** User should provide correct location **
                HeartConfig::Instance()->SetMeshFileName("../Data/OxfordRabbitHeart/full_res/OxfordRabbitHeart_binary"); //125um
                break;
            }
            default:
            {
                NEVER_REACHED;
                break;
            }
        };


        /* Declare the cell factory. Pass in initial conditions for the cell models (allowing the user to do this is useful additional functionality
         * implemented in `RegionBasedCellFactory`), and add the S2 stimulus. */

        RegionBasedCellFactoryWithApexS1<CellMahajan2008FromCellMLBackwardEuler> cell_factory;

        // These initial conditions where obtained by pacing a Mahajan single cell model every 220ms, for which APD is about 160.
        std::vector<double> init_conds = boost::assign::list_of(-85.44857927)(0.001433209884)(0.9859815937)(0.9702682509)(2.306476493e-05)(0.822170133)(0.03754955151)(5.610314802e-05)(0.1172144035)(0.02298433849)(0.02854004733)(0.1243973603)(0.1668792438)(0.004199973424)(0.1383989158)(0.004096766759)(0.938898475)(94.11365627)(0.0180643094)(17.35116329)(2.005202555)(0.655139027)(0.771004794)(105.6033538)(42.27624508)(40.1459155);
        cell_factory.SetCellModelInitialConditions(init_conds);

        if(s2_time != DBL_MAX)
        {
            c_vector<double,3> s2_stim_centre;
            s2_stim_centre(0) = 0.2438;
            s2_stim_centre(1) = 1.113;
            s2_stim_centre(2) = 1.225;
            double s2_stim_radius = 0.9;
            cell_factory.AddStimulatedSphere(s2_stim_centre, s2_stim_radius, s2_time);
        }

        /* The following code just sets up an output directory name based on the chosen options. */
        std::stringstream ss;
        switch(geometry)
        {
            case COARSERES_ISOTROPIC:
            {
                ss << "CoarseResIso";
                break;
            }
            case MEDIUMRES_ISOTROPIC:
            {
                ss << "MediumResIso";
                break;
            }
            case FULLRES_ISOTROPIC:
            {
                ss << "FullResIso";
                break;
            }
            default:
            {
                NEVER_REACHED;
                break;
            }
        };
        if(s2_time != DBL_MAX)
        {
            ss << "_s2_" << s2_time;
        }
        ss << "_End_" << end_time;
        if(notes.length()>0)
        {
            ss << "_" << notes;
        }
        HeartConfig::Instance()->SetOutputDirectory(ss.str());



        /* Run the simulation: */
        MonodomainProblem<3> cardiac_problem(&cell_factory);

        // cardiac_problem.SetWriteInfo();
        cardiac_problem.Initialise();
        cardiac_problem.Solve();


        /* Write the archive if required and print out timings. */
        if(write_archive)
        {
            CardiacSimulationArchiver<MonodomainProblem<3> >::Save(cardiac_problem, "archived_" + ss.str());
        }

        HeartEventHandler::Headings();
        HeartEventHandler::Report();
    }
};


#endif /*TESTREENTRYONRABBITMESHLITERATEPAPER_HPP_*/

