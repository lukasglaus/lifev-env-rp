//@HEADER
/*
*******************************************************************************

    Copyright (C) 2004, 2005, 2007 EPFL, Politecnico di Milano, INRIA
    Copyright (C) 2010 EPFL, Politecnico di Milano, Emory University

    This file is part of LifeV.

    LifeV is free software; you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    LifeV is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with LifeV.  If not, see <http://www.gnu.org/licenses/>.

*******************************************************************************
*/
//@HEADER


/*!
    @file
    @brief Simple test using the electromechanical passive constitutive laws


    @date 12-2014
    @author Simone Rossi <simone.rossi@epfl.ch>

    @contributor
    @mantainer Simone Rossi <simone.rossi@epfl.ch>
 */


#include <lifev/core/LifeV.hpp>
#include <lifev/electrophysiology/solver/ElectroETAMonodomainSolver.hpp>
#include <lifev/electrophysiology/solver/IonicModels/IonicMinimalModel.hpp>
#include <lifev/electrophysiology/util/HeartUtility.hpp>

#include <lifev/structure/solver/StructuralConstitutiveLawData.hpp>

#include <lifev/structure/solver/StructuralConstitutiveLaw.hpp>
#include <lifev/structure/solver/StructuralOperator.hpp>

#include <lifev/em/solver/mechanics/EMStructuralOperator.hpp>
#include <lifev/em/solver/mechanics/EMStructuralConstitutiveLaw.hpp>
#include <lifev/em/solver/EMETAFunctors.hpp>

#include <lifev/core/filter/ExporterEnsight.hpp>
#ifdef HAVE_HDF5
#include <lifev/core/filter/ExporterHDF5.hpp>
#endif
#include <lifev/core/filter/ExporterEmpty.hpp>

#include <lifev/eta/fem/ETFESpace.hpp>
#include <lifev/eta/expression/Integrate.hpp>

//#include <lifev/em/solver/mechanics/materials/EMMaterial.hpp>


#include <lifev/bc_interface/3D/bc/BCInterface3D.hpp>
#include <lifev/structure/solver/NeoHookeanMaterialNonLinear.hpp>

using namespace LifeV;



int main (int argc, char** argv)
{

#ifdef HAVE_MPI
    MPI_Init ( &argc, &argv );
#endif

    //start with the communicator
    boost::shared_ptr<Epetra_Comm>  comm ( new Epetra_MpiComm (MPI_COMM_WORLD) );
    if ( comm->MyPID() == 0 )
    {
        cout << "% using MPI" << std::endl;
    }

    //===========================================================
    //===========================================================
    //              READ DATAFILE and CREATE OUTPUT FOLDER
    //===========================================================
    //===========================================================
    GetPot command_line (argc, argv);
    const string data_file_name = command_line.follow ("data", 2, "-f", "--file");
    GetPot dataFile (data_file_name);

    EMData emdata;
    emdata.setup (dataFile);

    //When launching the executable use the flag:  -o OutputFolderName
    std::string problemFolder = EMUtility::createOutputFolder (command_line, *comm);


    //===========================================================
    //===========================================================
    //              LOAD MESH
    //===========================================================
    //===========================================================
    if ( comm->MyPID() == 0 )
    {
        std::cout << "Reading Mesh Name and Path...\n";
    }

    std::string meshName = dataFile ( "solid/space_discretization/mesh_file", "" );
    std::string meshPath = dataFile ( "solid/space_discretization/mesh_dir", "./" );

    typedef RegionMesh<LinearTetra>                         mesh_Type;
    typedef boost::shared_ptr<mesh_Type>                    meshPtr_Type;

    meshPtr_Type localSolidMesh ( new mesh_Type ( comm ) );
    meshPtr_Type fullSolidMesh ( new mesh_Type ( comm ) );

    MeshUtility::loadMesh (localSolidMesh, fullSolidMesh, meshName, meshPath);

    if ( comm->MyPID() == 0 )
    {
        std::cout << " Done!" << std::endl;
    }

    //===========================================================
    //===========================================================
    //              FINITE ELEMENT SPACES
    //===========================================================
    //===========================================================

    //Define the finite element space for bc and exporter
    // and the ET finite element space for assembly
    if ( comm->MyPID() == 0 )
    {
        std::cout << "setup spaces ... ";
    }

    typedef FESpace< RegionMesh<LinearTetra>, MapEpetra >               solidFESpace_Type;
    typedef boost::shared_ptr<solidFESpace_Type>                        solidFESpacePtr_Type;

    typedef ETFESpace< RegionMesh<LinearTetra>, MapEpetra, 3, 3 >       solidETFESpace_Type;
    typedef boost::shared_ptr<solidETFESpace_Type>                      solidETFESpacePtr_Type;



    std::string dOrder =  dataFile ( "solid/space_discretization/order", "P1");
    solidFESpacePtr_Type dFESpace ( new solidFESpace_Type (localSolidMesh, dOrder, 3, comm) );
    solidETFESpacePtr_Type dETFESpace ( new solidETFESpace_Type (localSolidMesh, & (dFESpace->refFE() ), & (dFESpace->fe().geoMap() ), comm) );

    if ( comm->MyPID() == 0 )
    {
        std::cout << " Done!" << std::endl;
    }

    //===========================================================
    //===========================================================
    //              BOUNDARY CONDITIONS
    //===========================================================
    //===========================================================
    if ( comm->MyPID() == 0 )
    {
        std::cout << "setup bc ... ";
    }

    typedef BCHandler                                          bc_Type;
    typedef StructuralOperator< RegionMesh<LinearTetra> >      physicalSolver_Type;
    typedef BCInterface3D< bc_Type, physicalSolver_Type >      bcInterface_Type;
    typedef boost::shared_ptr< bcInterface_Type >              bcInterfacePtr_Type;

    bcInterfacePtr_Type                     solidBC ( new bcInterface_Type() );
    solidBC->createHandler();
    solidBC->fillHandler ( data_file_name, "solid" );
    solidBC->handler()->bcUpdate ( *dFESpace->mesh(), dFESpace->feBd(), dFESpace->dof() );

    if ( comm->MyPID() == 0 )
    {
        std::cout << " Done!" << std::endl;
    }

    //===========================================================
    //===========================================================
    //              SOLID MECHANICS
    //===========================================================
    //===========================================================

    if ( comm->MyPID() == 0 )
    {
        std::cout << "setup constitutive law data ... ";
    }

    //Setup the data of the constitutive law
    boost::shared_ptr<StructuralConstitutiveLawData> dataStructure (new StructuralConstitutiveLawData( ) );
    dataStructure->setup (dataFile);

    if ( comm->MyPID() == 0 )
    {
        std::cout << " Done!" << std::endl;
    }

    //Setup the structural operator
    if ( comm->MyPID() == 0 )
    {
        std::cout << "setup structural operator ... " << std::endl;
    }
    //! 1. Constructor of the structuralSolver
    EMStructuralOperator< RegionMesh<LinearTetra> > solid;

    // solidBC is copied inside the StructuralOperator
    // any changes to it don't affect the M_BCh inside the solver
    // after a change in solidBC add it in the solver!
    solid.setup ( dataStructure, dFESpace, dETFESpace, solidBC -> handler(), comm);
    solid.setDataFromGetPot (dataFile);
    solid.EMMaterial()->setParameters(emdata);

    std::string passiveMaterialType  =  dataFile ( "solid/physics/EMPassiveMaterialType", "");
    if ( passiveMaterialType == "PHO" )
    {
        // load fibers and sheets fields
        std::string fiberFileName  =  dataFile ( "solid/space_discretization/fiber_file_name", "fiber");
        std::string sheetFileName  =  dataFile ( "solid/space_discretization/sheet_file_name", "sheet");
        std::string fiberFieldName =  dataFile ( "solid/space_discretization/fiber_field_name", "fiber");
        std::string sheetFieldName =  dataFile ( "solid/space_discretization/sheet_field_name", "sheet");
        std::string fiberDir       =  dataFile ( "solid/space_discretization/fiber_dir", "./");
        std::string sheetDir       =  dataFile ( "solid/space_discretization/sheet_dir", "./");

        if ( comm->MyPID() == 0 )
        {
            std::cout << "Importing fibers field\n";
            std::cout << "Fibers file name: " << fiberFileName  << "\n";
            std::cout << "Fibers field name: " << fiberFieldName << "\n";
            std::cout << "Fibers dir name: " << fiberDir       << "\n";
        }
        ElectrophysiologyUtility::importVectorField (  solid.EMMaterial()->fiberVectorPtr(),
                                                       fiberFileName,
                                                       fiberFieldName,
                                                       localSolidMesh,
                                                       fiberDir,
                                                       dOrder );
        if ( comm->MyPID() == 0 )
        {
            std::cout << "Importing sheets field\n";
            std::cout << "Sheets file name: " << fiberFileName  << "\n";
            std::cout << "Sheets field name: " << fiberFieldName << "\n";
            std::cout << "Sheets dir name: " << fiberDir       << "\n";
        }
        ElectrophysiologyUtility::importVectorField (  solid.EMMaterial()->sheetVectorPtr(),
                                                       sheetFileName,
                                                       sheetFieldName,
                                                       localSolidMesh,
                                                       sheetDir,
                                                       dOrder );
    }
    else
    {
        solid.EMMaterial() -> setupFiberVector( 1.0, 0.0, 0.0);
        solid.EMMaterial() -> setupSheetVector( 0.0, 1.0, 0.0);
    }
// solid.EMMaterial() -> setFiberVector ( solid.EMMaterial()->fiberVectorPtr() );
// solid.EMMaterial() -> setSheetVector ( solid.EMMaterial()->sheetVectorPtr() );

// solid.EMMaterial() -> setupFiberVector( 1.0, 0.0, 0.0);
// solid.EMMaterial() -> setupSheetVector( 0.0, 1.0, 0.0);

    solid.setNewtonParameters(dataFile);
    solid.buildSystem (1.0);

    if ( comm->MyPID() == 0 )
    {
        std::cout << " Done!" << std::endl;
    }

    //===========================================================
    //===========================================================
    //              CREATE EXPORTER and SAVE INITIAL SOLUTION
    //===========================================================
    //===========================================================

    boost::shared_ptr< Exporter<RegionMesh<LinearTetra> > > exporter;
    exporter.reset ( new ExporterHDF5<RegionMesh<LinearTetra> > ( dataFile, "structure" ) );
    exporter->setPostDir ( problemFolder );
    exporter->setMeshProcId ( localSolidMesh, comm->MyPID() );
    exporter->addVariable ( ExporterData<RegionMesh<LinearTetra> >::VectorField, "displacement", dFESpace, solid.displacementPtr(), UInt (0) );
    exporter->postProcess ( 0 );

    boost::shared_ptr< Exporter<RegionMesh<LinearTetra> > > fibers_exporter;
    fibers_exporter.reset ( new ExporterHDF5<RegionMesh<LinearTetra> > ( dataFile, "fibers" ) );
    fibers_exporter->setPostDir ( problemFolder );
    fibers_exporter->setMeshProcId ( localSolidMesh, comm->MyPID() );
    fibers_exporter->addVariable ( ExporterData<RegionMesh<LinearTetra> >::VectorField, "fibers", dFESpace, solid.EMMaterial()->fiberVectorPtr(), UInt (0) );
    fibers_exporter->addVariable ( ExporterData<RegionMesh<LinearTetra> >::VectorField, "sheets", dFESpace, solid.EMMaterial()->sheetVectorPtr(), UInt (0) );
    fibers_exporter->postProcess ( 0 );

    std::cout << "Displacement vector size: " << solid.displacementPtr()->size() << std::endl;
    std::cout << "Fiber vector size: " << solid.EMMaterial()->fiberVectorPtr()->size() << std::endl;
    //===========================================================
    //===========================================================
    //         SOLVE
    //===========================================================
    //===========================================================
    Real dt =  dataFile ( "solid/time_discretization/timestep", 0.1);
    Real endTime =  dataFile ( "solid/time_discretization/endtime", 1.0);
    ID LVFlag =  dataFile ( "solid/boundary_conditions/LV_flag", 0);
    Real LVPreloadPressure =  dataFile ( "solid/boundary_conditions/LV_preload_pressure", 0.0);
    bool deformedPressure =  dataFile ( "solid/boundary_conditions/deformed_pressure", 1 );

    if ( deformedPressure)
    {
        std::cout << "Setting pressure in the deformed configuration\n";
    }
    else
    {
        std::cout << "Setting pressure in the reference configuration\n";
    }
    solid.setBCFlag( LVFlag );

    for (Real time (0.0); time < endTime;)
    {
        time += dt;
        std::cout << "\nTime: " << time << std::endl;
        solid.data() -> dataTime() -> updateTime();
        std::cout << "----- Preload Time: " << time << std::endl;

        solidBC -> updatePhysicalSolverVariables();

        solid.bcH() = solidBC -> handler();

        solid.setLVPressureBC( -time*LVPreloadPressure );

        solid.iterate ( solidBC -> handler() , deformedPressure );

        exporter->postProcess ( time );
    }

    //===========================================================
    //===========================================================
    //              CLOSE EXPORTER
    //===========================================================
    //===========================================================
    exporter -> closeFile();
    fibers_exporter->closeFile();

#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    return 0;
}
