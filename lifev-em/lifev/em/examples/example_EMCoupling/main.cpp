//============================================//
// Includes
//============================================//

#include <lifev/core/LifeV.hpp>

// Passive material
#include <lifev/electrophysiology/solver/ElectroETAMonodomainSolver.hpp>
#include <lifev/electrophysiology/util/HeartUtility.hpp>
#include <lifev/structure/solver/StructuralConstitutiveLawData.hpp>
#include <lifev/em/solver/mechanics/EMStructuralOperator.hpp>
#include <lifev/em/solver/mechanics/EMStructuralConstitutiveLaw.hpp>
#include <lifev/bc_interface/3D/bc/BCInterface3D.hpp>

// Solver
#include <lifev/em/solver/EMSolver.hpp>

// Exporter
#include <lifev/core/filter/ExporterEnsight.hpp>
#include <lifev/core/filter/ExporterVTK.hpp>
#ifdef HAVE_HDF5
#include <lifev/core/filter/ExporterHDF5.hpp>
#endif
#include <lifev/core/filter/ExporterEmpty.hpp>

// Resize mesh
#include <lifev/core/mesh/MeshUtility.hpp>

// B.C. modification
#include <lifev/core/fem/BCVector.hpp>

// Circulation
#include <lifev/em/solver/circulation/Circulation.hpp>

// Volume computation
#include <lifev/em/solver/circulation/CirculationVolumeIntegrator.hpp>

//#include <fenv.h>



// Namespaces
using namespace LifeV;


//============================================//
// Functions
//============================================//

Real Iapp (const Real& t, const Real&  X, const Real& Y, const Real& Z, const ID& /*i*/)
{
    return ( Y > 1.5 && Y < 3 /*std::abs(X) < 1 && std::abs(Z-3) < 1 && Y < 0*/ /*&& Y < 0.25 && Z < 0.25 */ && t < 7 && t > 5 ? 30 : 0 );
    // setAppliedCurrent in electrophys. module.
}

Real potentialMultiplyerFcn (const Real& t, const Real&  X, const Real& Y, const Real& Z, const ID& /*i*/)
{
    return 0; // ( Y < 2.5 && Y > 0.5 ? 1.0 : 0.0 );
}


int main (int argc, char** argv)
{

//    feenableexcept(FE_INVALID | FE_OVERFLOW);
    
    
    //============================================//
    // Typedefs
    //============================================//
    
    typedef RegionMesh<LinearTetra>                         mesh_Type;
    typedef boost::shared_ptr<mesh_Type>                    meshPtr_Type;
    
    typedef boost::function < Real (const Real & t,
                                    const Real &   x,
                                    const Real &   y,
                                    const Real & z,
                                    const ID&   /*i*/ ) >   function_Type;
    
    typedef VectorEpetra                                    vector_Type;
    typedef boost::shared_ptr<vector_Type>                  vectorPtr_Type;
    
    typedef BCVector                                        bcVector_Type;
    typedef boost::shared_ptr<bcVector_Type>                bcVectorPtr_Type;
    
    typedef EMMonodomainSolver<mesh_Type>                   monodomain_Type;

    typedef FESpace< mesh_Type, MapEpetra >                 solidFESpace_Type;
    typedef boost::shared_ptr<solidFESpace_Type>            solidFESpacePtr_Type;
    
    typedef ETFESpace< mesh_Type, MapEpetra, 3, 3 >         solidETFESpace_Type;
    typedef boost::shared_ptr<solidETFESpace_Type>          solidETFESpacePtr_Type;
    
    typedef StructuralConstitutiveLawData                   constitutiveLaw_Type;
    typedef boost::shared_ptr<constitutiveLaw_Type>         constitutiveLawPtr_Type;
    
    typedef BCHandler                                       bc_Type;
    typedef StructuralOperator< mesh_Type >                 physicalSolver_Type;
    typedef BCInterface3D< bc_Type, physicalSolver_Type >   bcInterface_Type;
    typedef boost::shared_ptr< bcInterface_Type >           bcInterfacePtr_Type;
    

    //============================================//
    // Declare communicator and solver
    //============================================//
    
#ifdef HAVE_MPI
    MPI_Init ( &argc, &argv );
#endif

    boost::shared_ptr<Epetra_Comm>  comm ( new Epetra_MpiComm (MPI_COMM_WORLD) );
    if ( comm->MyPID() == 0 )
    {
        cout << "% using MPI" << std::endl;
    }
    //displayer->leaderPrint ("% using MPI");

    EMSolver<mesh_Type, monodomain_Type> solver(comm);

    
    //============================================//
    // Read data file and create output folder
    //============================================//

    GetPot command_line (argc, argv);
    const std::string data_file_name = command_line.follow ("data", 2, "-f", "--file");
    GetPot dataFile (data_file_name);
    std::string problemFolder = EMUtility::createOutputFolder (command_line, *comm);

    
    //============================================//
    // Setup material data
    //============================================//
    
    EMData emdata;
    emdata.setup (dataFile);
    
    
    //============================================//
    // Load mesh
    //============================================//
    
    if ( comm->MyPID() == 0 )
    {
        std::cout << "Load mesh...\n";
    }
    
    std::string meshName = dataFile("solid/space_discretization/mesh_file", "cube4.mesh");
    std::string meshPath = dataFile("solid/space_discretization/mesh_dir", "./");

    solver.loadMesh (meshName, meshPath);
    
    if ( comm->MyPID() == 0 )
    {
        std::cout << " Done!" << endl;
    }
    
    
    //============================================//
    // Resize mesh
    //============================================//
    
    if ( comm->MyPID() == 0 )
    {
        std::cout << "Resizing mesh..." << endl;
    }
    
    std::vector<Real> scale (3, dataFile("solid/space_discretization/mesh_scaling", 1.0));
    std::vector<Real> rotate { dataFile("solid/space_discretization/mesh_rotation_0", 0.0) , dataFile("solid/space_discretization/mesh_rotation_1", 0.0) , dataFile("solid/space_discretization/mesh_rotation_2", 0.0) };
    std::vector<Real> translate { dataFile("solid/space_discretization/mesh_translation_0", 0.0) , dataFile("solid/space_discretization/mesh_translation_1", 0.0) , dataFile("solid/space_discretization/mesh_translation_2", 0.0) };
    
    MeshUtility::MeshTransformer<mesh_Type> transformerFull (* (solver.fullMeshPtr() ) );
    MeshUtility::MeshTransformer<mesh_Type> transformerLocal (* (solver.localMeshPtr() ) );
    
    transformerFull.transformMesh (scale, rotate, translate);
    transformerLocal.transformMesh (scale, rotate, translate);
    
    if ( comm->MyPID() == 0 )
    {
        std::cout << " Done!" << endl;
    }
    
    
    //============================================//
    // Setup solver (including fe-spaces & b.c.)
    //============================================//
    
    if( 0 == comm->MyPID() )
    {
    	std::cout << "Setting up EM solver ... ";
    }
    
    solver.setup (dataFile);
    
    if( 0 == comm->MyPID() )
    {
    	std::cout << " done!" << std::endl;
    }

    
    //============================================//
    // Setup anisotropy vectors
    //============================================//
    
    if( 0 == comm->MyPID() )
    {
    	std::cout << "Setting up anisotropy vectors ...";
    }
 
    bool anisotropy = dataFile ( "solid/space_discretization/anisotropic", false );

    if ( anisotropy )
    {
        std::string fiberFileName  =  dataFile ( "solid/space_discretization/fiber_name", "FiberDirection");
        std::string sheetFileName  =  dataFile ( "solid/space_discretization/sheet_name", "SheetsDirection");
        std::string fiberFieldName =  dataFile ( "solid/space_discretization/fiber_fieldname", "fibers");
        std::string sheetFieldName =  dataFile ( "solid/space_discretization/sheet_fieldname", "sheets");
        std::string fiberDir       =  dataFile ( "solid/space_discretization/fiber_dir", "./");
        std::string sheetDir       =  dataFile ( "solid/space_discretization/sheet_dir", "./");
        std::string elementOrder   =  dataFile ( "solid/space_discretization/order", "P1");

        solver.setupFiberVector ( fiberFileName, fiberFieldName, fiberDir, elementOrder );
        solver.setupMechanicalSheetVector ( sheetFileName, sheetFieldName, sheetDir, elementOrder );
    }
    else
    {
        solver.setupFiberVector (1., 0., 0.);
        solver.setupSheetVector (0., 1., 0.);
    }
    
    if( 0 == comm->MyPID() )
    {
    	std::cout << " done!" << std::endl;
    }

    
    //============================================//
    // Initialize electrophysiology
    //============================================//
    
    if( 0 == comm->MyPID() )
    {
    	std::cout << "Initialize electrophysiology ... ";
    }
    
    solver.initialize();
    
    // Set potential on certain flags
    UInt lvendo = dataFile( "electrophysiology/flags/lvendo", 36 );
    //UInt rvendo = dataFile( "electrophysiology/flags/rvendo", 37 );
    //UInt rvseptum = dataFile( "electrophysiology/flags/rvseptum", 38 );
    ElectrophysiologyUtility::setValueOnBoundary ( * (solver.electroSolverPtr()->potentialPtr() ), solver.fullMeshPtr(), 1.0, lvendo );
    //ElectrophysiologyUtility::setValueOnBoundary ( * (solver.electroSolverPtr()->potentialPtr() ), solver.fullMeshPtr(), 1.0, rvendo );
    //ElectrophysiologyUtility::setValueOnBoundary ( * (solver.electroSolverPtr()->potentialPtr() ), solver.fullMeshPtr(), 1.0, rvseptum);
    
    // Restrict the potential set by a function
    vectorPtr_Type potentialMultiplyer ( new vector_Type ( solver.electroSolverPtr()->potentialPtr()->map() ) ); // or: vectorPtr_Type potentialMultiplyer ( new vector_Type ( *solver.electroSolverPtr()->potentialPtr() ) );
    function_Type potMult = &potentialMultiplyerFcn;
    solver.electroSolverPtr()->feSpacePtr()->interpolate( potMult, *potentialMultiplyer, 0 );
    *solver.electroSolverPtr()->potentialPtr() *= *potentialMultiplyer;
    
    if( 0 == comm->MyPID() )
    {
    	std::cout << " done!" << std::endl;
    }

    
    //============================================//
    // Building Matrices
    //============================================//

    if( 0 == comm->MyPID() )
    {
    	std::cout << "Building matrices ... ";
    }
    
    solver.twoWayCoupling();
    solver.structuralOperatorPtr()->setNewtonParameters(dataFile);
    solver.buildSystem();
    
    if( 0 == comm->MyPID() )
    {
    	std::cout << " done!" << std::endl;
    }

    
    //============================================//
    // Setup exporter for EMSolver
    //============================================//
    
    if ( 0 == comm->MyPID() )
    {
        std::cout << "Setting up exporters .. " << std::endl;
    }

    solver.setupExporters (problemFolder);
    
    if ( 0 == comm->MyPID() )
    {
        std::cout << " done!" << std::endl;
    }
    
    
    //============================================//
    // Setup vector & exporter for activation time
    //============================================//
    
    vectorPtr_Type activationTimeVector ( new vector_Type ( solver.electroSolverPtr()->potentialPtr() -> map() ) );
    *activationTimeVector = -1.0;
    
    ExporterHDF5< RegionMesh <LinearTetra> > activationTimeExporter;
    activationTimeExporter.setMeshProcId (solver.localMeshPtr(), solver.comm()->MyPID() );
    activationTimeExporter.addVariable (ExporterData<mesh_Type>::ScalarField, "Activation Time",
                                        solver.electroSolverPtr()->feSpacePtr(), activationTimeVector, UInt (0) );
    activationTimeExporter.setPrefix ("ActivationTime");
    activationTimeExporter.setPostDir (problemFolder);
    
    
    //============================================//
    // Electric stimulus function
    //============================================//
    
    function_Type stim = &Iapp;
    
    
    //============================================//
    // Body circulation
    //============================================//
    
    const std::string circulationInputFile = command_line.follow ("inputfile", 2, "-cif", "--cifile");
    const std::string circulationOutputFile = command_line.follow ("solution.txt", 2, "-cof", "--cofile");

    Circulation circulationSolver( circulationInputFile );
    if ( 0 == comm->MyPID() ) circulationSolver.exportSolution( circulationOutputFile );

    std::vector<std::vector<std::string> > bcNames { { "lv" , "p" } , { "rv" , "p" } };
    std::vector<double> bcValues ( 2 , 5 );
  
    
    // Flow rate between two vertices
    auto Q = [] (Circulation& c, const std::string& N1, const std::string& N2) { return c.solution ( std::vector<std::string> {N1, N2} ); };
    auto p = [] (Circulation& c, const std::string& N1) { return c.solution ( N1 ); };

    //============================================//
    // Kept-normal boundary conditions
    //============================================//

    // Get b.c. flags
    // ID LvFlag =  dataFile ( "solid/boundary_conditions/LvFlag", 0);
    
    // Boundary vector normal in deformed configuration
    // solver.structuralOperatorPtr() -> setBCFlag( LvFlag );
    
    
    //============================================//
    // Modifiable-value boundary condition
    //============================================//

    ID LVFlag =  dataFile ( "solid/boundary_conditions/VariableBoundaryConditions/LVFlag", 0 );
    ID RVFlag =  dataFile ( "solid/boundary_conditions/VariableBoundaryConditions/RVFlag", 0 );

    vectorPtr_Type endoLvVectorPtr( new vector_Type ( solver.structuralOperatorPtr() -> displacement().map(), Repeated ) );
    bcVectorPtr_Type pLvBCVectorPtr( new bcVector_Type( *endoLvVectorPtr, solver.structuralOperatorPtr() -> dispFESpacePtr() -> dof().numTotalDof(), 1 ) );
    
    // Todo: Normal boundary condition!!
    
    //solver.bcInterfacePtr() -> handler() -> addBC("LvPressure", LVFlag, Natural, Full, *pLvBCVectorPtr, 3); // BC for using function which keeps bc normal
    solver.bcInterfacePtr() -> handler() -> addBC("LvPressure", LVFlag, Natural, Normal, *pLvBCVectorPtr); // BC which is only at time zero normal
    solver.bcInterfacePtr() -> handler() -> bcUpdate( *solver.structuralOperatorPtr() -> dispFESpacePtr() -> mesh(), solver.structuralOperatorPtr() -> dispFESpacePtr() -> feBd(), solver.structuralOperatorPtr() -> dispFESpacePtr() -> dof() );
    //if ( 0 == comm->MyPID() ) solver.bcInterfacePtr() -> handler() -> showMe();


    //============================================//
    // Volume integrators
    //============================================//
    
    auto& disp = solver.structuralOperatorPtr() -> displacement();
    auto FESpace = solver.structuralOperatorPtr() -> dispFESpacePtr();
    auto dETFESpace = solver.electroSolverPtr() -> displacementETFESpacePtr();
    auto ETFESpace = solver.electroSolverPtr() -> ETFESpacePtr();
    
    VolumeIntegrator LV (std::vector<int> {23}, "Left Ventricle", solver.fullMeshPtr(), solver.localMeshPtr(), ETFESpace, FESpace);
    //VolumeIntegrator RV (std::vector<int> {37, 38}, "Right Ventricle", solver.fullMeshPtr(), solver.localMeshPtr(), ETFESpace);

    Real LVVolume = LV.volume(disp, dETFESpace, - 1);
    //Real RVVolume = RV.volume(disp, dETFESpace, 1);

    solver.saveSolution (-1.0);

    std::vector<double> VCirc { LV.volume(disp, dETFESpace, - 1) };
    std::vector<double> VCircNew (VCirc);
    std::vector<double> VFe (VCirc);
    std::vector<double> VFeNew (VCirc);

    
    //============================================//
    // Preload
    //============================================//
    
    Real pPreloadLvBC = dataFile ( "solid/boundary_conditions/LvPreloadPressure", 0.0);
    int preloadSteps = dataFile ( "solid/boundary_conditions/numPreloadSteps", 0);
    
    for (int i (1); i <= preloadSteps; i++)
    {
        if ( 0 == comm->MyPID() )
        {
            std::cout << "\n*********************";
            std::cout << "\nPreload step: " << i << " / " << preloadSteps;
            std::cout << "\n*********************\n";
        }
        
        // Update pressure b.c.
        *endoLvVectorPtr = - p(circulationSolver, "lv") * 0.001333224 * i / preloadSteps;
        pLvBCVectorPtr.reset ( ( new bcVector_Type (*endoLvVectorPtr, solver.structuralOperatorPtr() -> dispFESpacePtr() -> dof().numTotalDof(), 1) ) );
        solver.bcInterfacePtr() -> handler() -> modifyBC(LVFlag, *pLvBCVectorPtr);
        
        // Solve mechanics
        solver.bcInterfacePtr() -> updatePhysicalSolverVariables();
        solver.solveMechanics();
    }
    
    
    //============================================//
    // Time loop
    //============================================//
    
    //UInt couplingFeJacobianIter (
    
    Real dt_activation = solver.data().electroParameter<Real>("timestep");
    Real dt_mechanics = solver.data().solidParameter<Real>("timestep");
    Real endtime = solver.data().electroParameter<Real>("endtime");
    UInt saveIter = static_cast<UInt>( dt_mechanics / dt_activation );
    UInt maxiter = static_cast<UInt>( endtime / dt_activation ) ;
    Real t = 0;
    solver.saveSolution (0.0);
    
    VFeNew.at(0) = LVVolume;

    for (int k (1); k <= maxiter; k++)
    {
        solver.electroSolverPtr() -> registerActivationTime (*activationTimeVector, t, 0.9);
        
        if ( 0 == comm->MyPID() )
        {
            std::cout << "\n*********************";
            std::cout << "\nTIME = " << t+dt_activation;
            std::cout << "\n*********************\n";
        }

        t = t + dt_activation;
        unsigned int iter (0);

        solver.solveElectrophysiology (stim, t);
        solver.solveActivation (dt_activation);
        
        if ( k % saveIter == 0 )
        {
            //============================================//
            // Solve circlation
            //============================================//
            circulationSolver.iterate(dt_mechanics/1000, bcNames, bcValues, iter);
            VCircNew.at(0) = VCirc.at(0) + dt_mechanics/1000 * ( Q(circulationSolver, "la", "lv") - Q(circulationSolver, "lv", "sa") );
            
            //============================================//
            // Solve mechanics
            //============================================//
            *endoLvVectorPtr = - bcValues.at(0) * 0.001333224;
            pLvBCVectorPtr.reset ( ( new bcVector_Type (*endoLvVectorPtr, solver.structuralOperatorPtr() -> dispFESpacePtr() -> dof().numTotalDof(), 1) ) );
            solver.bcInterfacePtr() -> handler() -> modifyBC(LVFlag, *pLvBCVectorPtr);
            solver.structuralOperatorPtr() -> data() -> dataTime() -> setTime(t);
            solver.bcInterfacePtr() -> updatePhysicalSolverVariables();
            solver.solveMechanics();
            VFeNew.at(0) = LV.volume(disp, dETFESpace, - 1);

            if ( 0 == comm->MyPID() )
            {
                std::cout << "\n*********************";
                std::cout << "\nCoupling: Newton iteration nr. " << iter;
                std::cout << "\nPressure: " << bcValues.at(0);
                std::cout << "\nFE-Volume: " << VFeNew.at(0);
                std::cout << "\nCirculation-Volume: " << VCircNew.at(0);
                std::cout << "\nResidual = " << std::abs(VFeNew.at(0) - VCircNew.at(0));
                std::cout << "\n*********************\n";
            }
            
            if ( ! circulationSolver.coupling().converged(VFeNew, VCircNew, 1e-6) )
            {
                const double dpFactor ( 1e-3 );
                
                //============================================//
                // Jacobian fe
                //============================================//
                std::vector<double> VFePert (VFe);
                
                *endoLvVectorPtr = - bcValues.at(0) * 0.001333224 * (1.0 + 1e-3);
                pLvBCVectorPtr.reset ( ( new bcVector_Type (*endoLvVectorPtr, solver.structuralOperatorPtr() -> dispFESpacePtr() -> dof().numTotalDof(), 1) ) );
                solver.bcInterfacePtr() -> handler() -> modifyBC(LVFlag, *pLvBCVectorPtr);
                solver.structuralOperatorPtr() -> data() -> dataTime() -> setTime(t);
                solver.bcInterfacePtr() -> updatePhysicalSolverVariables();
                solver.solveMechanics();
                VFePert.at(0) = LV.volume(disp, dETFESpace, - 1);
                
                const double Jfe = ( VFePert.at(0) - VFeNew.at(0) ) / ( bcValues.at(0) * dpFactor );
                
                //============================================//
                // Newton iterations
                //============================================//
                while ( ! circulationSolver.coupling().converged(VFeNew, VCircNew, 1e-6) )
                {
                    ++iter;
                    
                    if ( 0 == comm->MyPID() )
                    {
                        std::cout << "\n*********************************************************************************";
                        std::cout << "\nCoupling: Newton iteration nr. " << iter;
                        std::cout << "\n*********************************************************************************\n";
                    }
                    
                    //============================================//
                    // Jacobian circulation
                    //============================================//
                    std::vector<double> bcValuesPert { bcValues.at(0) * (1.0 + dpFactor) , bcValues.at(1) };
                    std::vector<double> VCircPert (1);
                    circulationSolver.iterate(dt_mechanics/1000, bcNames, bcValuesPert, iter);
                    VCircPert.at(0) = VCirc.at(0) + dt_mechanics/1000 * ( Q(circulationSolver, "la", "lv") - Q(circulationSolver, "lv", "sa") );

                    const double Jcirc = ( VCircPert.at(0) - VCircNew.at(0) ) / ( bcValues.at(0) * dpFactor );
                    
                    //============================================//
                    // Jacobian fe
                    //============================================//
                    if ( ! iter % 5 )
                    {
                        std::vector<double> VFePert (VFe);
                    
                        *endoLvVectorPtr = - bcValues.at(0) * 0.001333224 * (1.0 + dpFactor);
                        pLvBCVectorPtr.reset ( ( new bcVector_Type (*endoLvVectorPtr, solver.structuralOperatorPtr() -> dispFESpacePtr() -> dof().numTotalDof(), 1) ) );
                        solver.bcInterfacePtr() -> handler() -> modifyBC(LVFlag, *pLvBCVectorPtr);
                        solver.structuralOperatorPtr() -> data() -> dataTime() -> setTime(t);
                        solver.bcInterfacePtr() -> updatePhysicalSolverVariables();
                        solver.solveMechanics();
                        VFePert.at(0) = LV.volume(disp, dETFESpace, - 1);
                    
                        const double Jfe = ( VFePert.at(0) - VFeNew.at(0) ) / ( bcValues.at(0) * dpFactor );
                    }
                    
                    //============================================//
                    // Update pressure b.c.
                    //============================================//
                    const double J = Jfe - Jcirc;
                    bcValues.at(0) += - std::pow(J, -1) * (VFeNew.at(0) - VCircNew.at(0));
                    bcValues.at(1) = std::max( bcValues.at(0) / 6 , 5.0 );

                    //============================================//
                    // Solve circulation
                    //============================================//
                    circulationSolver.iterate(dt_mechanics/1000, bcNames, bcValues, iter);
                    VCircNew.at(0) = VCirc.at(0) + dt_mechanics/1000 * ( Q(circulationSolver, "la", "lv") - Q(circulationSolver, "lv", "sa") );

                    //============================================//
                    // Solve mechanics
                    //============================================//
                    *endoLvVectorPtr = - bcValues.at(0) * 0.001333224;
                    pLvBCVectorPtr.reset ( ( new bcVector_Type (*endoLvVectorPtr, solver.structuralOperatorPtr() -> dispFESpacePtr() -> dof().numTotalDof(), 1) ) );
                    solver.bcInterfacePtr() -> handler() -> modifyBC(LVFlag, *pLvBCVectorPtr);
                    solver.structuralOperatorPtr() -> data() -> dataTime() -> setTime(t);
                    solver.bcInterfacePtr() -> updatePhysicalSolverVariables();
                    solver.solveMechanics();
                    VFeNew.at(0) = LV.volume(disp, dETFESpace, - 1);
                    
                    if ( 0 == comm->MyPID() )
                    {
                        std::cout << "\n*********************";
                        std::cout << "\nCoupling: Newton iteration nr. " << iter;
                        std::cout << "\nPressure: " << bcValues.at(0);
                        std::cout << "\nFE-Volume: " << VFeNew.at(0);
                        std::cout << "\nFE-Volume perturbed: " << VFePert.at(0);
                        std::cout << "\nFE-Jacobian: " << Jfe;
                        std::cout << "\nCirculation-Volume: " << VCircNew.at(0);
                        std::cout << "\nCirculation-Volume perturbed: " << VCircPert.at(0);
                        std::cout << "\nCirculation-Jacobian: " << Jcirc;
                        std::cout << "\nResidual = " << std::abs(VFeNew.at(0) - VCircNew.at(0));
                        std::cout << "\n*********************\n";
                    }
                }
            }
            
            if ( 0 == comm->MyPID() )
            {
                std::cout << "\n*********************";
                std::cout << "\nCoupling converged after " << iter << " iterations";
                std::cout << "\n*********************\n";
            }
            
            //============================================//
            // Set volumes equal to new ones
            //============================================//
            VCirc = VCircNew;
            VFe = VFeNew;
            
            //============================================//
            // Export circulation solution
            //============================================//
            if ( 0 == comm->MyPID() ) circulationSolver.exportSolution( circulationOutputFile );
        }
        
        //============================================//
        // Export FE-solution
        //============================================//
        solver.saveSolution(t);
        activationTimeExporter.postProcess(t); // usually stored after time loop
    }

    
    //============================================//
    // Close all exporters
    //============================================//
    
    solver.closeExporters();
    activationTimeExporter.closeFile();



#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    return 0;
}
