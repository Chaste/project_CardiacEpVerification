#ifndef TESTREENTRYRABBIT_HPP_
#define TESTREENTRYRABBIT_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/assign.hpp>
#include "CardiacSimulationArchiver.hpp"
#include "UblasCustomFunctions.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "RegionBasedCellFactory.hpp"
#include "LuoRudy1991.hpp"
#include "Mahajan2008BackwardEuler.hpp"
#include "Shannon2004.hpp"
#include "MonodomainProblem.hpp"
#include "TrianglesMeshReader.hpp"
#include "DistributedTetrahedralMesh.hpp"
#include "Hdf5ToMeshalyzerConverter.hpp"

typedef enum GeometryOption_
{
    COARSERES_ISOTROPIC = 0,
    MEDIUMRES_ISOTROPIC,
    FULLRES_ISOTROPIC
} GeometryOption;


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



class TestReentryRabbit : public CxxTest::TestSuite
{
private:
    void SetHeartConfigForTest()
    {
        HeartConfig::Instance()->Reset();
        HeartConfig::Instance()->SetKSPSolver("cg");
        HeartConfig::Instance()->SetKSPPreconditioner("bjacobi");
        HeartConfig::Instance()->SetUseAbsoluteTolerance(1e-5);
    }

public:
    void TestReentryS1andS2() throw(Exception)
    {
        double conductivity_scaled_by = 1.5;

        GeometryOption geometry = COARSERES_ISOTROPIC;
        double s2_time = 170; // set to DBL_MAX for no s2 to be added (and no s2 info added to output dir name)
        double end_time = 1000;
        double printing_time = 10.0;
        bool write_archive = false;
        std::string notes = ""; // added to output dir name

        ////////////////////////////////////////////////
        // Initial set up, timesteps, end time
        ////////////////////////////////////////////////
        SetHeartConfigForTest();
        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.1, 0.1, printing_time);
        HeartConfig::Instance()->SetSimulationDuration(end_time);

        ////////////////////////////////////////////////
        // Geometry and conductivities
        //
        // todo: where does 1.334, 0.176 come from
        ////////////////////////////////////////////////

        double base_fibre_conductivity = 1.4; //   = 1.75*7/(1.75+7)
        double base_trans_conductivity = DBL_MAX;

        switch(geometry)
        {
            case COARSERES_ISOTROPIC:
            case MEDIUMRES_ISOTROPIC:
            case FULLRES_ISOTROPIC:
            {
                HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(base_fibre_conductivity/conductivity_scaled_by,base_fibre_conductivity/conductivity_scaled_by, base_fibre_conductivity/conductivity_scaled_by));
                break;
            }
            case FULLRES_ANISOTROPIC:
            {
                NEVER_REACHED; // set base_trans_conductivity above
                HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(base_fibre_conductivity/conductivity_scaled_by,base_trans_conductivity/conductivity_scaled_by, base_trans_conductivity/conductivity_scaled_by));
                break;
            }
            default:
            {
                NEVER_REACHED;
                break;
            }
        };

        switch(geometry)
        {
            case COARSERES_ISOTROPIC:
            {
             	HeartConfig::Instance()->SetMeshFileName("apps/texttest/weekly/Propagation3d/heart_chaste2_renum_i_triangles");
                break;
            }
            case MEDIUMRES_ISOTROPIC:
            {
                HeartConfig::Instance()->SetMeshFileName("../../data/low_res_downsampled_oxford_rabbit/heart_chaste_renum_i_triangles");
                break;
            }
            case FULLRES_ISOTROPIC:
            {
                HeartConfig::Instance()->SetMeshFileName("../Data/OxfordRabbitHeart/fullres/OxfordRabbitHeart_binary");
                break;
            }
            default:
            {
                NEVER_REACHED;
                break;
            }
        };

        //////////////////////////////////////////////////
        // Cell factory: cell initial conditons, s1, s2
        //////////////////////////////////////////////////

        RegionBasedCellFactoryWithApexS1<CellMahajan2008FromCellMLBackwardEuler> cell_factory;

        // obtained by pacing mahajan every 220ms (down in increments from 300ms), so APD is about 160.
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

        ////////////////////////////////////////////////
        // Output directory name
        ////////////////////////////////////////////////
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
        ss << "_CondScale_" << conductivity_scaled_by;
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



        ////////////////////////////////////////////////
        // Run simulation
        ////////////////////////////////////////////////
        MonodomainProblem<3> cardiac_problem(&cell_factory);

        //// cluster sometimes has issues with communication
        // cardiac_problem.SetWriteInfo();

        cardiac_problem.Initialise();
        cardiac_problem.Solve();


        ////////////////////////////////////////////////
        // Post-solve
        ////////////////////////////////////////////////
        if(write_archive)
        {
            CardiacSimulationArchiver<MonodomainProblem<3> >::Save(cardiac_problem, "archived_" + ss.str());
        }

        HeartEventHandler::Headings();
        HeartEventHandler::Report();
    }
};


#endif /*TESTREENTRYRABBIT_HPP_*/

