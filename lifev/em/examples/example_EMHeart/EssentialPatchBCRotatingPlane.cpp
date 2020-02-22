/*

 * EssentialPatchBCRotatingPlane.cpp
 *
 *  Created on: Feb 5, 2020
 *      Author: lglaus
 */

#ifndef EssentialPatchBCRotatingPlane_hpp
#define EssentialPatchBCRotatingPlane_hpp

//#include "EssentialPatchBCRotatingPlane.h"
#include <stdio.h>
#include <lifev/em/examples/example_EMHeart/EssentialPatchBC.hpp>
#include <lifev/em/examples/example_EMHeart/GenericFactory.hpp>
//#include <string>
#include <cmath>
#include <vector>
#define PI 3.14159265359

namespace LifeV
{

class EssentialPatchBCRotatingPlane : public EssentialPatchBC {
public:

    vectorPtr_Type m_p2currentPositionVector;
    Real m_maxDisplacement;
    Vector3D normal_vector;
    Vector3D axis_perp_t;
    double angleOfTime;
    Vector3D starting_point;
    Vector3D direction_to_axis; //2020.02.10 lg
    Vector3D axis_direction; //2020.02.10 lg
    Real distance_to_axis; //2020.02.10 lg
    double maximum_angle; //2020.02.10 lg
    double minimum_angle; //2020.02.10 lg
    Vector3D pointOnHeart;
    Real m_tduration;
    Real m_tmax;
    int rotation_direction;
    unsigned int nodeOnPatchCounterone;
    unsigned int nodeOnPatchCountertwo;
    
void setup(const GetPot& dataFile, const std::string& name)
{
    super::setup(dataFile, name);

    //Import starting point on heart
    for (UInt i (0); i < 3 ;++i )
    {
    pointOnHeart[i] = dataFile( ("solid/boundary_conditions/" + m_Name + "/pointOnHeart").c_str(), 0.0, i );
    }
    
    //Import direction from starting point on heart to point in axis
    for (UInt j (0); j < 3; ++j)
    {
    direction_to_axis[j] = dataFile ( ("solid/boundary_conditions/" + m_Name + "/direction_to_axis").c_str() , 0.0 , j );
    }
    direction_to_axis = normalize_vector(direction_to_axis);
    
    //Import distance between starting point on heart and point in axis
    distance_to_axis = dataFile ( ("solid/boundary_conditions/" + m_Name + "/distance_to_axis").c_str(), 1.0 );
    
    //Import direction of rotation axis
    for (UInt j (0); j < 3; ++j)
    {
    axis_direction[j] = dataFile ( ("solid/boundary_conditions/" + m_Name + "/axis_direction").c_str() , 0.0 , j );
    }
    axis_direction = normalize_vector(axis_direction);
    
    starting_point=calculate_pAxis(pointOnHeart,direction_to_axis,distance_to_axis);//starting_point already defined in EssentialPatchBC.hpp
    //std::cout<<"\nsetup:starting_point = ("<<starting_point[0]<<","<<starting_point[1]<<","<<starting_point[2]<<")"<<endl;
    
    //Import the initial opening angle of the patches
    maximum_angle = dataFile ( ("solid/boundary_conditions/" + m_Name + "/maximum_angle").c_str(), 1.0 );
    maximum_angle = (maximum_angle * PI)/180;
    //std::cout<<"\nsetup:initial maximum_angle = "<<maximum_angle*180/PI<<"° degree"<<endl;
    rotation_direction = dataFile ( ("solid/boundary_conditions/" + m_Name + "/rotation_direction").c_str(), 1.0 );
    
    //Import the final (=smallest) opening angle of the patches
    minimum_angle = dataFile ( ("solid/boundary_conditions/" + m_Name + "/minimum_angle").c_str(), 1.0 );
    minimum_angle = (minimum_angle * 3.141)/180;
    //std::cout<<"seteup:initial minimum_angle = "<<minimum_angle*180/PI<<"° degree"<<endl;
    //In order to have no translation of the patches
    m_maxDisplacement=0;
    
    // Temporal activation parameter
    m_tmax = dataFile ( "solid/patches/tmax", 0. );
    m_tduration = dataFile ( "solid/patches/tduration", 0. );
    //std::cout<<"\nsetup:initial m_tduration = "<<m_tduration;
    
    //initial normal vector for applyPatchBC
    angleOfTime=calculate_angleOfTime(0.0);
    normal_vector=createNormalVector (0.0);
    m_patchDirection=normal_vector;
    
    //if ( solver.comm()->MyPID() == 0 ) std::cout<<"setup completed";
    //std::cout<<"\nsetup:setup completed"<<endl;
    
    //displayImportantVariables(0);
}

//Normalizes a vector
Vector3D normalize_vector (Vector3D vector)
    {
    Real abs=std::sqrt(vector[0]*vector[0]+vector[1]*vector[1]+vector[2]*vector[2]);
        if (abs !=0)
        {
            vector[0]=vector[0]/abs;
            vector[1]=vector[1]/abs;
            vector[2]=vector[2]/abs;
        }
        if (abs!=1){std::cout<<"\nCareful:During import a vector had to be normalized"<<endl;}
        if (abs=0){std::cout<<"\nCareful:Absolute value of a supposed to be normalized vector is zero"<<endl;}
    return vector;
    }

//Calculate the location of pAxis, a point within the rotating axis
Vector3D calculate_pAxis (const Vector3D pointOnHeart,const Vector3D direction_to_axis,const Real distance_to_axis)
    {
        Vector3D p_axis;
        p_axis[0]=pointOnHeart[0]+direction_to_axis[0]*distance_to_axis;
        p_axis[1]=pointOnHeart[1]+direction_to_axis[1]*distance_to_axis;
        p_axis[2]=pointOnHeart[2]+direction_to_axis[2]*distance_to_axis;
        
        return p_axis;
    }

//Calculate the opening angle(Degree) in function of time
double calculate_angleOfTime (Real time)
    {
        double angle;
        
        if (std::fmod(time,800.)/m_tduration < 0.5)
            {
                angle=maximum_angle/2 - (maximum_angle/2 - minimum_angle/2)*(std::fmod(time,800.)/(m_tduration/2));
            }
        if (std::fmod(time,800.)/m_tduration >= 0.5 && std::fmod(time,800.)/m_tduration<1)
            {
                angle=minimum_angle/2 + (maximum_angle/2 - minimum_angle/2)*(std::fmod(time,800.)-(m_tduration/2))/(m_tduration/2);
            }
        if (std::fmod(time,800.)/m_tduration >= 1)
            {
                angle=maximum_angle/2;
            }
        //std::cout << "\nangle= " <<angle*180/PI<<"degree"<<endl;
        angle=angle*rotation_direction;
        //std::cout << "\nangle after multiplication with rotation_direction= " <<angle*180/PI<<"degree"<<endl;
        
        return angle;
    }
    
/* old: out of duty since 18.02.2020
//Calculate the opening angle(Degree) in function of time
double calculate_angleOfTime (Real time)
    {
        double angle;
        
        if (std::fmod(time,m_tduration)/m_tduration < 0.5)
            {
                angle=maximum_angle/2 - (maximum_angle/2 - minimum_angle/2)*(std::fmod(time,(m_tduration/2)))/(m_tduration/2);
               
            }
            else
            {
                angle=minimum_angle/2 + (maximum_angle/2 - minimum_angle/2)*(std::fmod(time,(m_tduration/2)))/(m_tduration/2);
            }
        std::cout << "\nangle= " <<angle*180/PI<<"degree"<<endl;
        angle=angle*rotation_direction;
        std::cout << "\nangle after multiplication with rotation_direction= " <<angle*180/PI<<"degree"<<endl;
        
        return angle;
    }
*/
    
    
    void displayImportantVariables (Real time){
        std::cout<<"\n//============================================";
        std::cout<<"\n// displayImportantVariables";
        std::cout<<"\n//============================================";
        std::cout<<"\nTime= "<< time;
        std::cout<<"\nstarting_point = ("<<starting_point[0]<<","<<starting_point[1]<<","<<starting_point[2]<<")";
        std::cout<<"\nsetup:normal_vector = ("<<normal_vector[0]<<","<<normal_vector[1]<<","<<normal_vector[2]<<")";
        std::cout<<"\nangleOfTime= "<< angleOfTime*180/PI<<" degree";
        std::cout<<"\ninitial maximum_angle = "<<maximum_angle*180/PI<<" degree";
        std::cout<<"\ninitial minimum_angle = "<<minimum_angle*180/PI<<" degree";
        std::cout<<"\ninitial m_tduration = "<<m_tduration;
    }

//cos in degree or radian?
//Changes direction and rotates "direction_to_axis" around "axis_direction" with the angle "angleOfTime"
Vector3D rotateVectorAroundAxis (double angleOfTime)
    {
    Vector3D normalOfPatch;
    Vector3D rotatedVector;
                    
    Real product1=(axis_direction[0]*axis_direction[0]*(1-std::cos(angleOfTime))+std::cos(angleOfTime))*direction_to_axis[0];
    
    Real product2= (axis_direction[0]*axis_direction[1]*(1-std::cos(angleOfTime))-axis_direction[2]*std::sin(angleOfTime))*direction_to_axis[1];
    
    Real product3= (axis_direction[0]*axis_direction[2]*(1-std::cos(angleOfTime))+axis_direction[1]*std::sin(angleOfTime))*direction_to_axis[2];
    
    rotatedVector[0]= -1 * (product1 + product2 + product3);
                    
    Real product4= (axis_direction[0]*axis_direction[1]*(1-std::cos(angleOfTime))+axis_direction[2]*std::sin(angleOfTime))*direction_to_axis[0];
     
    Real product5=(axis_direction[1]*axis_direction[1]*(1-std::cos(angleOfTime))+std::cos(angleOfTime))*direction_to_axis[1];
     
    Real product6= (axis_direction[1]*axis_direction[2]*(1-std::cos(angleOfTime))-axis_direction[0]*std::sin(angleOfTime))*direction_to_axis[2];
     
    rotatedVector[1]= -1 * (product4 + product5 + product6);
                    
    Real product7= (axis_direction[0]*axis_direction[2]*(1-std::cos(angleOfTime))-axis_direction[1]*std::sin(angleOfTime))*direction_to_axis[0];
        
    Real product8= (axis_direction[1]*axis_direction[2]*(1-std::cos(angleOfTime))+axis_direction[0]*std::sin(angleOfTime))*direction_to_axis[1];
        
    Real product9= (axis_direction[2]*axis_direction[2]*(1-std::cos(angleOfTime))+std::cos(angleOfTime))*direction_to_axis[2];
        
    rotatedVector[2]= -1 * (product7 + product8 + product9);
    
    //std::cout<<"\nrotatedVector= ("<<rotatedVector[0]<<","<<rotatedVector[1]<<","<<rotatedVector[2]<<")";
    return rotatedVector;
    }
 
 //Cross product between two vectors creates a vector normal to them: c = a x b
 Vector3D createNormalVector (Real time)
    {
    //std::cout<<"\ncreateNormalVector: time= "<<time;
    angleOfTime = calculate_angleOfTime(time);
    //std::cout<<"\ncreateNormalVector: angle= "<<angle*180/PI;
    axis_perp_t = rotateVectorAroundAxis(angle);
    //std::cout<<"\ncreateNormalVector: axis_perp_t= ("<<axis_perp_t[0]<<","<<axis_perp_t[1]<<","<<axis_perp_t[2]<<")";
    //std::cout<<"\ncreateNormalVector: axis_direction= ("<<axis_direction[0]<<","<<axis_direction[1]<<","<<axis_direction[2]<<")";
    normal_vector=axis_perp_t.cross(axis_direction);
    //std::cout<<"\ncreateNormalVector: axis_perp_t cross axis_direction=normalToPlane= ("<<normalToPlane[0]<<","<<normalToPlane[1]<<","<<normalToPlane[2]<<")";
    /*
    normalToPlane[0]=axis_direction[1]*axis_perp_t[2]-axis_direction[2]*axis_perp_t[1];
    normalToPlane[1]=axis_direction[2]*axis_perp_t[0]-axis_direction[0]*axis_perp_t[2];
    normalToPlane[2]=axis_direction[0]*axis_perp_t[1]-axis_direction[1]*axis_perp_t[0];
    */
        
    normal_vector.normalize();
    //std::cout<<"createNormalVector: normalToPlane (normalized)= ("<<normalToPlane[0]<<","<<normalToPlane[1]<<","<<normalToPlane[2]<<")";
    normal_vector=normalToPlane*rotation_direction;
    //std::cout<<"createNormalVector: rotation_direction = "<<rotation_direction;
    //std::cout<<"createNormalVector: normalToPlane (normalized+multiplied with rotation_direction)= ("<<normalToPlane[0]<<","<<normalToPlane[1]<<","<<normalToPlane[2]<<")";
    return normal_vector;
    }
                                                                                   
//ggf hier eine Zeitabhängigkeit in normal_vector und starting_point einfügen, m_maxDisplacement ggf durch activationFunction(time)) ersetzen
const bool nodeOnPatch(const Vector3D& coord, const Real& time)
{

    bool nodeInArea = false;

        //m_maxDisplacement = dataFile ( ("solid/boundary_conditions/" + m_Name + "/displacement").c_str(), 1.0 );

        if((normal_vector[0]*coord[0] + normal_vector[1]*coord[1] + normal_vector[2]*coord[2] - normal_vector[0]*starting_point[0]- normal_vector[1]*starting_point[1]-normal_vector[2]*starting_point[2]) <= 0)
            {
                nodeInArea = true;
            }
            else
            {
                nodeInArea = false;
            }


            return nodeInArea;
}
                                                                                   
const bool nodeOnPatchCurrent(const Vector3D& coord, const Real& time)
{
    bool nodeInArea = 0;

    //as shift we had + 0.35
        if((normal_vector[0]*coord[0] + normal_vector[1]*coord[1] + normal_vector[2]*coord[2] - normal_vector[0]*starting_point[0]- normal_vector[1]*starting_point[1]-normal_vector[2]*starting_point[2]) <= 0)
        {
            nodeInArea = true;
        }
        else
        {
            nodeInArea = false;
        }


        return nodeInArea;

}

    void nodeOnPatchdisplayer(const Vector3D& coord, const Real& time)
    {
        
        std::cout<<"\nnodeOnPatchCurrent:axis_perp_t= ("<<axis_perp_t[0]<<","<<axis_perp_t[1]<<","<<axis_perp_t[2]<<")";
        std::cout<<"\nnodeOnPatchCurrent:normal_vector= ("<<normal_vector[0]<<","<<normal_vector[1]<<","<<normal_vector[2]<<")";
        std::cout<<"\nnodeOnPatchCurrent:coord= ("<<coord[0]<<","<<coord[1]<<","<<coord[2]<<")";
        std::cout<<"\nnodeOnPatchCurrent:starting_point= ("<<starting_point[0]<<","<<starting_point[1]<<","<<starting_point[2]<<")";
        std::cout<<"\nnodeOnPatchCurrent:Skalarprodukt="<<(normal_vector[0]*coord[0] + normal_vector[1]*coord[1] + normal_vector[2]*coord[2] - normal_vector[0]*starting_point[0]- normal_vector[1]*starting_point[1]-normal_vector[2]*starting_point[2]);
        
    }

void modifyPatchArea(EMSolver<RegionMesh<LinearTetra>, EMMonodomainSolver<RegionMesh<LinearTetra> > >& solver,const int& newFlag, const Real& time)
{
    if ( solver.comm()->MyPID() == 0 ) std::cout << "\nWE ARE IN MODIFY PATCH AREA " << std::endl;

            auto p2FeSpace = solver.electroSolverPtr()->feSpacePtr();
            auto p2dFeSpace = solver.structuralOperatorPtr()->dispFESpacePtr();
            FESpace<RegionMesh<LinearTetra>, MapEpetra > p1FESpace (p2FeSpace->mesh(), "P1", 1, p2FeSpace->mesh()->comm());

            //create an epetra vector to set it equal to one where it is part of patch
            VectorEpetra p1ScalarFieldFaces (p1FESpace.map());

            p1ScalarFieldFaces *= 0.0;

            Int p1ScalarFieldFacesDof = p1ScalarFieldFaces.epetraVector().MyLength();

            int globalIdArray[p1ScalarFieldFacesDof];

            p1ScalarFieldFaces.blockMap().MyGlobalElements(globalIdArray);

            m_patchFlag = newFlag;

            //std::cout << "This is patchFlag in modify Patch Area: " << m_patchFlag << std::endl;
            const auto& mesh = solver.localMeshPtr(); // variable mesh which we use later for for loop; we assign a local Mesh pointer to it

            const auto& meshfull = solver.fullMeshPtr();
            //auto numPoints = meshfull->numPoints();

            //getPatchRegion(solver, m_patchFlag, time);

                // Create patches by changing the markerID (flag) locally
                nodeOnPatchCounterone=0; //here we just initalise an unsigned integer variable
                nodeOnPatchCountertwo=0;
    
                for (int j(0); j < mesh->numBoundaryFacets(); j++) //returns number of boundary facets
                        {
                            auto& face = mesh->boundaryFacet(j);
                            auto faceFlag = face.markerID();
                            //std::cout << "This is face marker ID: " << face.markerID() << std::endl;
                            //if (faceFlag == m_PrevFlag)
                            //{
                            int numPointsOnFace(0);

                            for (int k(0); k < 3; ++k) //k < 3 was before; this is just a test
                                 {
                                     //auto coord = face.point(k).coordinates();
                                     ID pointGlobalId = face.point(k).id();
                                     auto coord = face.point(k).coordinates();
                                     auto pointInPatch = nodeOnPatchCurrent(coord, time);
                                     
                                    if ( solver.comm()->MyPID() == 0 ) std::cout<<"\n\nnodedisplayer for: "<<m_Name;
                                     if ( solver.comm()->MyPID() == 0 ) std::cout<<"\nnode number "<<j<<" of "<<mesh->numBoundaryFacets();
                                     if ( solver.comm()->MyPID() == 0 ) nodeOnPatchdisplayer(coord, time);
                                     
                                     if(pointInPatch == true)
                                         {
                                             nodeOnPatchCounterone++;
                                             ++numPointsOnFace;
                                             for(int n = 0; n < p1ScalarFieldFacesDof; n++)
                                                 {
                                                     if(pointGlobalId == globalIdArray[n])
                                                         {
                                                         //++numPointsOnFace;
                                                         p1ScalarFieldFaces[pointGlobalId] = 1.0;
                                                         }
                                                 }
                                         }

                                 }
                
                             if (numPointsOnFace >= 1) // if there are more than two points on face we execute the if statement; not completly sure here
                             {
                                     //std::cout << "We are now changing the faceID" << std::endl;
                                     //std::cout << "" << std::endl;
                                     face.setMarkerID(m_patchFlag);
                                     //std::cout << "This is the set face flag: " ;
                                     //face.Marker::showMe(std::cout);
                                    nodeOnPatchCountertwo++;
                             }
                
                           //}
                        }

                if ( solver.comm()->MyPID() == 0 ) std::cout<<"\nOn patch "<<m_Name<<" "<<nodeOnPatchCounterone<<" nodes had to be moved"<<nodeOnPatchCountertwo;
        
                m_patchFacesLocationPtr.reset (new vector_Type (p2FeSpace->map() ));
                *m_patchFacesLocationPtr = p2FeSpace->feToFEInterpolate(p1FESpace, p1ScalarFieldFaces);
                //*m_patchFacesLocationPtr = p1ScalarFieldFaces;

 }

//this is directional vectorfield for p2 elements

void modifyPatchBC(EMSolver<RegionMesh<LinearTetra>, EMMonodomainSolver<RegionMesh<LinearTetra> > >& solver, const Real& time, int& PatchFlag)
{
    angleOfTime=calculate_angleOfTime(time);
    normal_vector=createNormalVector (time);
    m_patchDirection=normal_vector;
    
    if ( solver.comm()->MyPID() == 0 )std::cout<<"\nin modifyPatchBC: "<<m_Name;
    if ( solver.comm()->MyPID() == 0 )displayImportantVariables(time);
    
    //std::cout << "This is value of time variable: "<< time << std::endl;
    //int adder = 12;
        //const int constantPatchFlag = PatchFlag;
        //const int constantPatchFlag;
    //std::string patchNameAdder = std::to_string(adder); //converts double variable time to string
    //m_Name = m_Name + patchNameAdder;

    const int currentPatchFlag = PatchFlag;

    /*
    if(PatchFlag == 900 && time != 0)
    {
        m_flagIncreaserOne += 10;
        currentPatchFlag = m_flagIncreaserOne;
    }

    if(PatchFlag == 901 && time != 0)
    {
        m_flagIncreaserTwo += 10;
        currentPatchFlag = m_flagIncreaserTwo;
    }
    */
    //std::cout << "This is modified PatchName: " << m_Name << std::endl;

    //std::cout << "This is patchFlag in modifyPatchBC which we give modifyPatchArea: " << constantPatchFlag << std::endl;


        auto dFeSpace = solver.structuralOperatorPtr()->dispFESpacePtr();
        
        modifyPatchArea(solver, currentPatchFlag, time);

        Real currentPatchDisp = activationFunction(time);
        if ( 0 == solver.comm()->MyPID() ) std::cout << "\nEssentialPatchBC: " << m_Name << " rotated to angle " << (angleOfTime*180)/PI << " degree of ["<<(minimum_angle*180*rotation_direction)/(2*PI)<<","<<(maximum_angle*180*rotation_direction)/(2*PI)<<"]";

        m_patchDispPtr = directionalVectorField(solver,dFeSpace, m_patchDirection, currentPatchDisp, time);

        m_patchDispBCPtr.reset( new bcVector_Type( *m_patchDispPtr, dFeSpace->dof().numTotalDof(), 1 ) );
    /*
    if (49.99  <= time && time  <= 50.02)
    {
         solver.bcInterfacePtr() -> handler()->addBC (m_Name, currentPatchFlag,  Essential, Component, *m_patchDispBCPtr, m_patchComponent);
         solver.bcInterfacePtr()->handler()->modifyBC(currentPatchFlag, *m_patchDispBCPtr);
        if ( 0 == solver.comm()->MyPID() ) solver.bcInterfacePtr() -> handler() -> showMe();
    }

    if (time > 51)
    {
        std::cout << "We are now modifing the BC which we inserted later" << std::endl;
        solver.bcInterfacePtr()->handler()->modifyBC(currentPatchFlag, *m_patchDispBCPtr);
    }
    */
    //solver.bcInterfacePtr() -> handler()->addBC (m_Name, currentPatchFlag,  Essential, Component, *m_patchDispBCPtr, m_patchComponent);
    //solver.bcInterfacePtr() -> handler()->addBC (m_Name, m_patchFlag,  Essential, Component, *m_patchDispBCPtr, m_patchComponent);
    // if ( 0 == solver.comm()->MyPID() ) solver.bcInterfacePtr() -> handler() -> showMe();
     solver.bcInterfacePtr()->handler()->modifyBC(currentPatchFlag, *m_patchDispBCPtr); //This was the version how it worked
        //solver.bcInterfacePtr()->handler()->modifyBC(m_patchFlag, *m_patchDispBCPtr); //this is old version
       //solver.bcInterfacePtr() -> handler()->addBC (m_Name, m_patchFlag,  Essential, Component, *m_patchDispBCPtr, m_patchComponent);//idea is now that we add everytime a new BC
}

vectorPtr_Type directionalVectorField(EMSolver<RegionMesh<LinearTetra>, EMMonodomainSolver<RegionMesh<LinearTetra> > >& solver,const boost::shared_ptr<FESpace<RegionMesh<LinearTetra>, MapEpetra >> dFeSpace, Vector3D& direction, const Real& disp, const Real& time)
{
    Vector3D current_point_on_plane=starting_point;
    Real distance;

        //m_p2currentPositionVector = vectorPtr_Type (new VectorEpetra( dFeSpace->map(), Repeated ));
    //std::cout << "NOW WE ARE IN DIRECTIONAL VECTOR FIELD" << std::endl;


            auto p2PositionVector = p2PositionVectorInitial(dFeSpace, solver);
    /*
        if(time == 0.0)
        {
        m_p2currentPositionVector = vectorPtr_Type (new VectorEpetra( dFeSpace->map(), Repeated ));
        *m_p2currentPositionVector = p2PositionVectorInitial(dFeSpace, solver);
        }
        */

            vectorPtr_Type p2PatchDisplacement (new VectorEpetra( dFeSpace->map(), Repeated ));


            auto nCompLocalDof = p2PatchDisplacement->epetraVector().MyLength() / 3;


            /*
            bool sameMap = m_p2currentPositionVector->blockMap().SameAs(p2PatchDisplacement->blockMap());

            if(sameMap == true)
            {
                std::cout << "The maps are the same" << std::endl;
            }
            else
            {
                std::cout << "The maps are not the same" << std::endl;
            }

            bool samePoint = m_p2currentPositionVector->blockMap().PointSameAs(p2PatchDisplacement->blockMap());
            if(samePoint == true)
            {
                    std::cout << "The points are the same" << std::endl;
            }
            else
            {
                     std::cout << "The points are not the same" << std::endl;
            }
            */
    
            direction.normalize(); //need to be careful; direction and normal_vector aren't the same anymore; after that direction is the normalised normal_vector
    /*
            if(normal_vector[2] != 0.0)
            {
                //In thoughtdofdofdofdofs we set cooridnates x and y equal to zero and solve for z coordinate and store it in current_point_on_plane[0]
                //here we just added the max value of distance vector (3.5155), let's see how it works
                current_point_on_plane[2] = (normal_vector[0]*starting_point[0] + normal_vector[1]*starting_point[1] + normal_vector[2]*starting_point[2] +activationFunction(time))/normal_vector[2];
                current_point_on_plane[1] = 0;
                current_point_on_plane[0] = 0;

                //std::cout << "This is coordinate of current point on plane" << current_point_on_plane[2] << std::endl;

             }
            else if (normal_vector[2] == 0.0 && normal_vector[1] != 0.0)
            {
                
                current_point_on_plane[1] = (normal_vector[0]*starting_point[0] + normal_vector[1]*starting_point[1] + normal_vector[2]*starting_point[2]  + activationFunction(time))/normal_vector[1];
                current_point_on_plane[2] = 0;
                current_point_on_plane[0] = 0;
            }
            else if (normal_vector[2] == 0.0 && normal_vector[1] == 0.0 && normal_vector[0] != 0.0)
            {
                current_point_on_plane[0] = (normal_vector[0]*starting_point[0] + normal_vector[1]*starting_point[1] + normal_vector[2]*starting_point[2]  +activationFunction(time))/normal_vector[0];
                current_point_on_plane[1] = 0;
                current_point_on_plane[2] = 0;
            }
            else
             {
                std::cout << "A normal  vector in the data file of (0, 0 , 0) doesn't make sense" << std::endl;
            }
    */

        //std::cout << "THIS IS CURRENT POINT ON PLANE "  << current_point_on_plane[0] << "       " << current_point_on_plane[1] << "         " << current_point_on_plane[2] << std::endl;
        // std::cout << current_point_on_plane[1] << std::endl;
        // std::cout << current_point_on_plane[2] << std::endl;



            //std::cout << "This is Length of EpetraDisplacementVector in DirectionalVectorfield: " << p2PatchDisplacement->epetraVector().MyLength() << std::endl;

            for (int j (0); j < nCompLocalDof; ++j)
            {
                        // Get coordinates

                        UInt iGID = p2PatchDisplacement->blockMap().GID (j);
                        UInt jGID = p2PatchDisplacement->blockMap().GID (j + nCompLocalDof);
                        UInt kGID = p2PatchDisplacement->blockMap().GID (j + 2 * nCompLocalDof);

                        /*
                        UInt iGID = p2PositionVector.blockMap().GID (j);
                        UInt jGID = p2PositionVector.blockMap().GID (j + nCompLocalDof);
                        UInt kGID = p2PositionVector.blockMap().GID (j + 2 * nCompLocalDof);
                        */

                        Vector3D coordinates;
                
                        coordinates(0) = p2PositionVector[iGID];
                        coordinates(1) = p2PositionVector[jGID];
                        coordinates(2) = p2PositionVector[kGID];

                        //if ( solver.comm()->MyPID() == 0 ) std::cout << "\n\nIteration number" << j << "of" <<nCompLocalDof-1 << "Coordinate vector" << coordinates(0) << "," << coordinates(1) << "," << coordinates(2) << ",";
                
                        /*
                        coordinates(0) = (*m_p2currentPositionVector)[iGID];
                        coordinates(1) = (*m_p2currentPositionVector)[jGID];
                        coordinates(2) = (*m_p2currentPositionVector)[kGID];
                        */

                        Vector3D QP; //define here the vector that goes from Q (point on plane) to point P

                        QP = coordinates - current_point_on_plane;

                        //std::cout << "THESE ARE COORDINATES OF QP: " << QP[0] << "       " << QP[1] << "       " << QP[2] << std::endl;


                     //   QP[0] = coordinates[0] - current_point_on_plane[0];
                     //   QP[1] = coordinates[1] - current_point_on_plane[1];
                     //   QP[2] = coordinates[2] - current_point_on_plane[2];


                        //next we do the projection of the vector QP ond the normalvector to the plane; we look at the sign which it has in the if statement

                        //if(QP[0]*direction[0] + QP[1]*direction[1] + QP[2]*direction[2] <=0) //we use here normalised normal vector
                        if(QP.dot(direction) <= 0) //here i have change to > 0
                        {
                            //if the dot product is smaller equal zero, then we want to apply a displacement; for that we calculate the distance from point P to plane and then say this is the displacement we want
            
                            distance = abs(QP.dot(direction));
                            //distance = abs(QP[0]*direction[0]+ QP[1]*direction[1]+QP[2]*direction[2]);
                            //std::cout << " THIS IS CALCULATED DISTANCE IN VECTORFIELD: " << distance << std::endl;


                            //////////////////Here we write distance to file
                            auto currentprocessor = dFeSpace->mesh()->comm()->MyPID();
                            std::ostringstream oss;//this is to convert int to string
                            oss << currentprocessor;

                            std::string path = "/cluster/home/lglaus/LIFE5/lifev-env-rp/lifev-em-install-debug/lifev/em/examples/example_EMHeart/distancefiles/distances_" + oss.str() + ".dat"; //2020.02.05 lg
                            //std::string path = "/cluster/home/pamstad/LIFE5/lifev-env/lifev-em-build/lifev/em/examples/example_EMHeart/distancefiles/distances_" + oss.str() + ".dat"; //2020.02.05 lg
                            std::ofstream writer(path.c_str(), std::ios_base::app);
                            if(writer.is_open()==false)
                            {
                                std::cout << "error occured while opening the file" << std::endl;
                            }
                            writer << coordinates(0) << "\t\t" << coordinates(1) << "\t\t" << coordinates(2) << "\t\t"  << distance << std::endl;
                            writer.close();
                            //////////////////Here the writing to the file ends


                        }
                        else
                        {
                            distance = 0.0;
                        }

                        Vector3D displacement_vector;
                        displacement_vector[0] = distance*direction[0]; // distance*direction[0];
                        displacement_vector[1] = distance*direction[1]; //  distance*direction[1];
                        displacement_vector[2] = distance*direction[2]; // distance*direction[2];

                //std::cout << "This is displacmeent VEctor: " <<  displacement_vector(0) << "          " << displacement_vector(1) << "         " << displacement_vector(2) << std::endl;

                        (*p2PatchDisplacement)[iGID] = displacement_vector[0];
                        (*p2PatchDisplacement)[jGID] = displacement_vector[1];
                        (*p2PatchDisplacement)[kGID] = displacement_vector[2];

                        /*
                        (*m_p2currentPositionVector)[iGID] = (*m_p2currentPositionVector)[iGID] + (*p2PatchDisplacement)[iGID];
                        (*m_p2currentPositionVector)[jGID] = (*m_p2currentPositionVector)[jGID] + (*p2PatchDisplacement)[jGID];
                        (*m_p2currentPositionVector)[kGID] = (*m_p2currentPositionVector)[kGID] + (*p2PatchDisplacement)[kGID];
                        */
                        /*
                        auto currentprocessor = dFeSpace->mesh()->comm()->MyPID();
                        std::ostringstream oss;//this is to convert int to string
                        oss << currentprocessor;

                        std::string path_two = "/cluster/home/pamstad/LIFE5/lifev-env/lifev-em-build/lifev/em/examples/example_EMHeart/currentPositionVector/currentPositionVector_" + oss.str() + ".dat";
                        std::ofstream writer_two(path_two.c_str(), std::ios_base::app);
                        if(writer_two.is_open()==false)
                        {
                                      std::cout << "error occured while opening the file" << std::endl;
                        }
                        writer_two << (*m_p2currentPositionVector)[iGID] << "\t\t" << (*m_p2currentPositionVector)[jGID] << "\t\t" << (*m_p2currentPositionVector)[kGID] << std::endl;
                        writer_two.close();
                        */


            }
            //p2PositionVector += *p2PatchDisplacement;
            //*m_p2currentPositionVector += *p2PatchDisplacement;
        
        if(time == 0)
        {
            *p2PatchDisplacement *= 0.0;
        }

            return p2PatchDisplacement;

}

    vectorPtr_Type patch_attacher(EMSolver<RegionMesh<LinearTetra>, EMMonodomainSolver<RegionMesh<LinearTetra> > >& solver,const boost::shared_ptr<FESpace<RegionMesh<LinearTetra>, MapEpetra >> dFeSpace, Vector3D& direction, const Real& disp, const Real& time)
    {
        Vector3D current_point_on_plane=starting_point;
        Real distance;

            //m_p2currentPositionVector = vectorPtr_Type (new VectorEpetra( dFeSpace->map(), Repeated ));
        //std::cout << "NOW WE ARE IN DIRECTIONAL VECTOR FIELD" << std::endl;


                auto p2PositionVector = p2PositionVectorInitial(dFeSpace, solver);
        /*
            if(time == 0.0)
            {
            m_p2currentPositionVector = vectorPtr_Type (new VectorEpetra( dFeSpace->map(), Repeated ));
            *m_p2currentPositionVector = p2PositionVectorInitial(dFeSpace, solver);
            }
            */

                vectorPtr_Type p2PatchDisplacement (new VectorEpetra( dFeSpace->map(), Repeated ));


                auto nCompLocalDof = p2PatchDisplacement->epetraVector().MyLength() / 3;


                /*
                bool sameMap = m_p2currentPositionVector->blockMap().SameAs(p2PatchDisplacement->blockMap());

                if(sameMap == true)
                {
                    std::cout << "The maps are the same" << std::endl;
                }
                else
                {
                    std::cout << "The maps are not the same" << std::endl;
                }

                bool samePoint = m_p2currentPositionVector->blockMap().PointSameAs(p2PatchDisplacement->blockMap());
                if(samePoint == true)
                {
                        std::cout << "The points are the same" << std::endl;
                }
                else
                {
                         std::cout << "The points are not the same" << std::endl;
                }
                */
        
                direction.normalize(); //need to be careful; direction and normal_vector aren't the same anymore; after that direction is the normalised normal_vector
        /*
                if(normal_vector[2] != 0.0)
                {
                    //In thoughtdofdofdofdofs we set cooridnates x and y equal to zero and solve for z coordinate and store it in current_point_on_plane[0]
                    //here we just added the max value of distance vector (3.5155), let's see how it works
                    current_point_on_plane[2] = (normal_vector[0]*starting_point[0] + normal_vector[1]*starting_point[1] + normal_vector[2]*starting_point[2] +activationFunction(time))/normal_vector[2];
                    current_point_on_plane[1] = 0;
                    current_point_on_plane[0] = 0;

                    //std::cout << "This is coordinate of current point on plane" << current_point_on_plane[2] << std::endl;

                 }
                else if (normal_vector[2] == 0.0 && normal_vector[1] != 0.0)
                {
                    
                    current_point_on_plane[1] = (normal_vector[0]*starting_point[0] + normal_vector[1]*starting_point[1] + normal_vector[2]*starting_point[2]  + activationFunction(time))/normal_vector[1];
                    current_point_on_plane[2] = 0;
                    current_point_on_plane[0] = 0;
                }
                else if (normal_vector[2] == 0.0 && normal_vector[1] == 0.0 && normal_vector[0] != 0.0)
                {
                    current_point_on_plane[0] = (normal_vector[0]*starting_point[0] + normal_vector[1]*starting_point[1] + normal_vector[2]*starting_point[2]  +activationFunction(time))/normal_vector[0];
                    current_point_on_plane[1] = 0;
                    current_point_on_plane[2] = 0;
                }
                else
                 {
                    std::cout << "A normal  vector in the data file of (0, 0 , 0) doesn't make sense" << std::endl;
                }
        */

            //std::cout << "THIS IS CURRENT POINT ON PLANE "  << current_point_on_plane[0] << "       " << current_point_on_plane[1] << "         " << current_point_on_plane[2] << std::endl;
            // std::cout << current_point_on_plane[1] << std::endl;
            // std::cout << current_point_on_plane[2] << std::endl;



                //std::cout << "This is Length of EpetraDisplacementVector in DirectionalVectorfield: " << p2PatchDisplacement->epetraVector().MyLength() << std::endl;

                for (int j (0); j < nCompLocalDof; ++j)
                {
                            // Get coordinates

                            UInt iGID = p2PatchDisplacement->blockMap().GID (j);
                            UInt jGID = p2PatchDisplacement->blockMap().GID (j + nCompLocalDof);
                            UInt kGID = p2PatchDisplacement->blockMap().GID (j + 2 * nCompLocalDof);

                            /*
                            UInt iGID = p2PositionVector.blockMap().GID (j);
                            UInt jGID = p2PositionVector.blockMap().GID (j + nCompLocalDof);
                            UInt kGID = p2PositionVector.blockMap().GID (j + 2 * nCompLocalDof);
                            */

                            Vector3D coordinates;
                    
                            coordinates(0) = p2PositionVector[iGID];
                            coordinates(1) = p2PositionVector[jGID];
                            coordinates(2) = p2PositionVector[kGID];

                            //if ( solver.comm()->MyPID() == 0 ) std::cout << "\n\nIteration number" << j << "of" <<nCompLocalDof-1 << "Coordinate vector" << coordinates(0) << "," << coordinates(1) << "," << coordinates(2) << ",";
                    
                            /*
                            coordinates(0) = (*m_p2currentPositionVector)[iGID];
                            coordinates(1) = (*m_p2currentPositionVector)[jGID];
                            coordinates(2) = (*m_p2currentPositionVector)[kGID];
                            */

                            Vector3D QP; //define here the vector that goes from Q (point on plane) to point P

                            QP = coordinates - current_point_on_plane;

                            //std::cout << "THESE ARE COORDINATES OF QP: " << QP[0] << "       " << QP[1] << "       " << QP[2] << std::endl;


                         //   QP[0] = coordinates[0] - current_point_on_plane[0];
                         //   QP[1] = coordinates[1] - current_point_on_plane[1];
                         //   QP[2] = coordinates[2] - current_point_on_plane[2];


                            //next we do the projection of the vector QP ond the normalvector to the plane; we look at the sign which it has in the if statement

                            //if(QP[0]*direction[0] + QP[1]*direction[1] + QP[2]*direction[2] <=0) //we use here normalised normal vector
                            if(QP.dot(direction) <= 0) //here i have change to > 0
                            {
                                //if the dot product is smaller equal zero, then we want to apply a displacement; for that we calculate the distance from point P to plane and then say this is the displacement we want
                
                                distance = abs(QP.dot(direction));
                                //distance = abs(QP[0]*direction[0]+ QP[1]*direction[1]+QP[2]*direction[2]);
                                //std::cout << " THIS IS CALCULATED DISTANCE IN VECTORFIELD: " << distance << std::endl;


                                //////////////////Here we write distance to file
                                auto currentprocessor = dFeSpace->mesh()->comm()->MyPID();
                                std::ostringstream oss;//this is to convert int to string
                                oss << currentprocessor;

                                std::string path = "/cluster/home/lglaus/LIFE5/lifev-env-rp/lifev-em-install-debug/lifev/em/examples/example_EMHeart/distancefiles/distances_" + oss.str() + ".dat"; //2020.02.05 lg
                                //std::string path = "/cluster/home/pamstad/LIFE5/lifev-env/lifev-em-build/lifev/em/examples/example_EMHeart/distancefiles/distances_" + oss.str() + ".dat"; //2020.02.05 lg
                                std::ofstream writer(path.c_str(), std::ios_base::app);
                                if(writer.is_open()==false)
                                {
                                    std::cout << "error occured while opening the file" << std::endl;
                                }
                                writer << coordinates(0) << "\t\t" << coordinates(1) << "\t\t" << coordinates(2) << "\t\t"  << distance << std::endl;
                                writer.close();
                                //////////////////Here the writing to the file ends


                            }
                            else
                            {
                                distance = 0.0;
                            }

                            Vector3D displacement_vector;
                            displacement_vector[0] = distance*direction[0]; // distance*direction[0];
                            displacement_vector[1] = distance*direction[1]; //  distance*direction[1];
                            displacement_vector[2] = distance*direction[2]; // distance*direction[2];

                    //std::cout << "This is displacmeent VEctor: " <<  displacement_vector(0) << "          " << displacement_vector(1) << "         " << displacement_vector(2) << std::endl;

                            (*p2PatchDisplacement)[iGID] = displacement_vector[0];
                            (*p2PatchDisplacement)[jGID] = displacement_vector[1];
                            (*p2PatchDisplacement)[kGID] = displacement_vector[2];

                            /*
                            (*m_p2currentPositionVector)[iGID] = (*m_p2currentPositionVector)[iGID] + (*p2PatchDisplacement)[iGID];
                            (*m_p2currentPositionVector)[jGID] = (*m_p2currentPositionVector)[jGID] + (*p2PatchDisplacement)[jGID];
                            (*m_p2currentPositionVector)[kGID] = (*m_p2currentPositionVector)[kGID] + (*p2PatchDisplacement)[kGID];
                            */
                            /*
                            auto currentprocessor = dFeSpace->mesh()->comm()->MyPID();
                            std::ostringstream oss;//this is to convert int to string
                            oss << currentprocessor;

                            std::string path_two = "/cluster/home/pamstad/LIFE5/lifev-env/lifev-em-build/lifev/em/examples/example_EMHeart/currentPositionVector/currentPositionVector_" + oss.str() + ".dat";
                            std::ofstream writer_two(path_two.c_str(), std::ios_base::app);
                            if(writer_two.is_open()==false)
                            {
                                          std::cout << "error occured while opening the file" << std::endl;
                            }
                            writer_two << (*m_p2currentPositionVector)[iGID] << "\t\t" << (*m_p2currentPositionVector)[jGID] << "\t\t" << (*m_p2currentPositionVector)[kGID] << std::endl;
                            writer_two.close();
                            */


                }
                //p2PositionVector += *p2PatchDisplacement;
                //*m_p2currentPositionVector += *p2PatchDisplacement;
            
            if(time == 0)
            {
                *p2PatchDisplacement *= 0.0;
            }

                return p2PatchDisplacement;

    }
    
vector_Type displayDirectionalVectorField(EMSolver<RegionMesh<LinearTetra>, EMMonodomainSolver<RegionMesh<LinearTetra> > >& solver, const Real& time)
{
    Real distance;
    
    Vector3D normalVector=normal_vector;
    Vector3D startingPoint=starting_point;
    
    //Vector3D normalVector; //already in class defined?
    //Vector3D startingPoint;//already in class defined?
    
    //direction.normalize(); // 2020.02.06 lg

    if ( solver.comm()->MyPID() == 0 ){std::cout<<"\nstartingPoint in displayDirectionalVectorfield= "<<startingPoint[0]<<","<<startingPoint[1]<<","<<startingPoint[2]<<")";}
    
    //startingPoint[0] = -3.76487;
    //startingPoint[1] = -10.6687;
    //startingPoint[2] = -0.36572;
    
    
    if ( solver.comm()->MyPID() == 0 ){std::cout<<"\nnormalVector in displayDirectionalVectorfield= "<<normalVector[0]<<","<<normalVector[1]<<","<<normalVector[2]<<")";}
    
    //normalVector[0] = 0.665647;
    //normalVector[1] = 0.695607;
    //normalVector[2] = -0.270367;
    
    if ( solver.comm()->MyPID() == 0 ){std::cout<<"\nangleOfTime in displayDirectionalVectorfield= "<<angleOfTime;}
    
    Vector3D direction = normalVector; // 2020.02.08 lg
    direction.normalize(); // 2020.02.08 lg
    
    Vector3D current_point_on_plane=starting_point;
    
    if ( solver.comm()->MyPID() == 0 ) std::cout << "\n\nvector direction=( " << direction[0] << ", " << direction[1] << ", " << direction[2] << ")\n\n"; // 2020.02.06 lg
    
    //first we want to set up the initial vector
    auto p2dFeSpace = solver.structuralOperatorPtr()->dispFESpacePtr();
        const auto& meshFull = solver.fullMeshPtr();
            FESpace<RegionMesh<LinearTetra>, MapEpetra > p1dFESpace (p2dFeSpace->mesh(), "P1", 3, p2dFeSpace->mesh()->comm());

                    VectorEpetra p1PositionVector (p1dFESpace.map());
                    p1PositionVector *= 0.0;

                    Int p1nLocalPoints = p1PositionVector.epetraVector().MyLength() / 3;

                    //std::cout << "This is length of positionVector: " << p1PositionVector.epetraVector().MyLength() << std::endl;

                    for (int j (0); j < p1nLocalPoints; j++)
                    {
                         UInt iGID = p1PositionVector.blockMap().GID (j);
                         UInt jGID = p1PositionVector.blockMap().GID (j + p1nLocalPoints);
                         UInt kGID = p1PositionVector.blockMap().GID (j + 2 * p1nLocalPoints);


                         //Vector3D coord = meshFull->point(iGID).coordinates();

                         p1PositionVector[iGID] = meshFull->point(iGID).x();
                         p1PositionVector[jGID] = meshFull->point(iGID).y();
                         p1PositionVector[kGID] = meshFull->point(iGID).z();



                    }

                    VectorEpetra p2PositionVector ( m_dispPtr->map() );
                    p2PositionVector = p2dFeSpace->feToFEInterpolate(p1dFESpace, p1PositionVector);

                    //now we have the positionvector


                    if(time == 0.0)
                    {
                        m_currentPositionVector = vectorPtr_Type (new VectorEpetra( p1dFESpace.map(), Repeated ));
                        //*m_currentPositionVector = p2PositionVector;
                    }


                    vectorPtr_Type p2PatchDisplacement (new VectorEpetra( p1dFESpace.map(), Repeated ));
                    auto nCompLocalDof = p2PatchDisplacement->epetraVector().MyLength() / 3;



                    /*
                    if(normalVector[2] != 0.0)
                                {
                                    //In thoughts we set cooridnates x and y equal to zero and solve for z coordinate and store it in current_point_on_plane[0]
                                    //here we just added the max value of distance vector (3.5155), let's see how it works
                                    current_point_on_plane[2] = (normalVector[0]*startingPoint[0] + normalVector[1]*startingPoint[1] + normalVector[2]*startingPoint[2] +activationFunction(time))/normalVector[2];
                                    current_point_on_plane[1] = 0;
                                    current_point_on_plane[0] = 0;

                                    //std::cout << "This is coordinate of current point on plane" << current_point_on_plane[2] << std::endl;

                                }
                     */
                     
                    /*
                                else if (normal_vector[2] == 0.0 && normal_vector[1] != 0.0)
                                {

                                    current_point_on_plane[1] = (normal_vector[0]*starting_point[0] + normal_vector[1]*starting_point[1] + normal_vector[2]*starting_point[2]  + activationFunction(time))/normal_vector[1];
                                    current_point_on_plane[2] = 0;
                                    current_point_on_plane[0] = 0;
                                }
                                else if (normal_vector[2] == 0.0 && normal_vector[1] == 0.0 && normal_vector[0] != 0.0)
                                {
                                    current_point_on_plane[0] = (normal_vector[0]*starting_point[0] + normal_vector[1]*starting_point[1] + normal_vector[2]*starting_point[2]  +activationFunction(time))/normal_vector[0];
                                    current_point_on_plane[1] = 0;
                                    current_point_on_plane[2] = 0;
                                }
                                else
                                {
                                    std::cout << "A normal  vector in the data file of (0, 0 , 0) doesn't make sense" << std::endl;
                                }
                        */


                     for (int j (0); j < nCompLocalDof; ++j)
                                {
                                            // Get coordinates

                                            UInt iGID = p2PatchDisplacement->blockMap().GID (j);
                                            UInt jGID = p2PatchDisplacement->blockMap().GID (j + nCompLocalDof);
                                            UInt kGID = p2PatchDisplacement->blockMap().GID (j + 2 * nCompLocalDof);

                                            /*
                                            UInt iGID = p2PositionVector.blockMap().GID (j);
                                            UInt jGID = p2PositionVector.blockMap().GID (j + nCompLocalDof);
                                            UInt kGID = p2PositionVector.blockMap().GID (j + 2 * nCompLocalDof);
                                            */

                                            Vector3D coordinates;

                                            /*
                                            coordinates(0) = (*m_currentPositionVector)[iGID];
                                            coordinates(1) = (*m_currentPositionVector)[jGID];
                                            coordinates(2) = (*m_currentPositionVector)[kGID];
                                            */

                                            coordinates(0) = p2PositionVector[iGID];
                                            coordinates(1) = p2PositionVector[jGID];
                                            coordinates(2) = p2PositionVector[kGID];

                                            Vector3D QP; //define here the vector that goes from Q (point on plane) to point P

                                            QP = coordinates - current_point_on_plane;

                                            if ( solver.comm()->MyPID() == 0 )
                                            {
                                                //std::cout << "These are coordinates: " << coordinates(0) << "        " << coordinates(1) << "       " << coordinates(2) << std::endl;
                                                //std::cout << "This is current point on plane " << current_point_on_plane(0) << "         " << current_point_on_plane(1) << "        " << current_point_on_plane(2) << std::endl;
                                            }


                                            if(QP.dot(direction) <= 0) //here i have change to > 0
                                            {
                                                        //distance = 1.0;
                                                        distance = abs(QP.dot(direction));
                                            }
                                            else
                                              {
                                                        distance = 0.0;
                                              }

                                            Vector3D displacement_vector;
                                            displacement_vector[0] = distance*direction[0];
                                            displacement_vector[1] = distance*direction[1];
                                            displacement_vector[2] = distance*direction[2];

                                            if ( solver.comm()->MyPID() == 0 )
                                            {
                                                //std::cout << "This is dispalcement vector: " << displacement_vector[0] << "       " << displacement_vector[1] << "          " << displacement_vector[2] << std::endl;
                                            }

                                            (*p2PatchDisplacement)[iGID] = displacement_vector[0];
                                            (*p2PatchDisplacement)[jGID] = displacement_vector[1];
                                            (*p2PatchDisplacement)[kGID] = displacement_vector[2];

                                                                    (*m_currentPositionVector)[iGID] = (*m_currentPositionVector)[iGID] + (*p2PatchDisplacement)[iGID];
                                                                    (*m_currentPositionVector)[jGID] = (*m_currentPositionVector)[jGID] + (*p2PatchDisplacement)[jGID];
                                                                    (*m_currentPositionVector)[kGID] = (*m_currentPositionVector)[kGID] + (*p2PatchDisplacement)[kGID];

                                }

                     return *p2PatchDisplacement;
}

//THIS IS THE CORRECTED VERSION OF THE DIRECTIONAL VECTORFIELD; LETS SEE
/*
vectorPtr_Type EssentialPatchBCMovingPlane::directionalVectorField(const boost::shared_ptr<FESpace<RegionMesh<LinearTetra>, MapEpetra >> dFeSpace, Vector3D& direction, const Real& disp, const Real& time) const
{
    Vector3D current_point_on_plane;
    Real distance;

    // auto p2PositionVector = p2PositionVectorDisplaced(dFeSpace);
            auto p2PositionVector = p2PositionVectorInitial(dFeSpace);

            auto nCompLocalDof_two = p2PositionVector.epetraVector().MyLength() / 3;

            //Here we allocate memory for the epetra vector
            vectorPtr_Type p2PatchDisplacement (new VectorEpetra( dFeSpace->map(), Repeated ));

            auto nCompLocalDof = p2PatchDisplacement->epetraVector().MyLength() / 3;

            direction.normalize(); //need to be careful; direction and normal_vector aren't the same anymore; after that, direction is the normalised normal_vector

            //the next steps we do is to get a point on the plane; we call that point Q; with that point we can calculate a vector QP to any point P of the heart
            //that vector QP we can later project onto the normalised normalvector of the plane; with that we can do two things:
            // 1) can determine if left or right to the plane ; 2) can also determine the distance from point to plane which is distance for displacement

            if(normal_vector[2] != 0)
            {
                //In thoughts we set cooridnates x and y equal to zero and solve for z coordinate and store it in current_point_on_plane[0]
                //here we just added the max value of distance vector (3.5155), let's see how it works
                current_point_on_plane[2] = (normal_vector[0]*starting_point[0] + normal_vector[1]*starting_point[1] + normal_vector[2]*starting_point[2]-0.5 +activationFunction(time))/normal_vector[2];
                current_point_on_plane[1] = 0;
                current_point_on_plane[0] = 0;

                //std::cout << "This is coordinate of current point on plane" << current_point_on_plane[2] << std::endl;

             }
            else if (normal_vector[2] == 0 && normal_vector[1] != 0)
            {

                current_point_on_plane[1] = (normal_vector[0]*starting_point[0] + normal_vector[1]*starting_point[1] + normal_vector[2]*starting_point[2] - 0.5 + activationFunction(time))/normal_vector[1];
                current_point_on_plane[2] = 0;
                current_point_on_plane[0] = 0;
            }
            else if (normal_vector[2] == 0 && normal_vector[1] == 0 && normal_vector[0] != 0)
            {
                current_point_on_plane[0] = (normal_vector[0]*starting_point[0] + normal_vector[1]*starting_point[1] + normal_vector[2]*starting_point[2] - 0.5  +activationFunction(time))/normal_vector[0];
                current_point_on_plane[1] = 0;
                current_point_on_plane[2] = 0;
            }
            else
             {
                std::cout << "A normal  vector in the data file of (0, 0 , 0) doesn't make sense" << std::endl;
            }

            //now we have a point on our plane which is stored in the variable current_point_on_Plane

            //In the next steps we want to get the coordinates of the points of the heart

            for (int j (0); j < nCompLocalDof; ++j)
            {
                        // Get coordinates
                        UInt iGID = p2PatchDisplacement->blockMap().GID (j);
                        UInt jGID = p2PatchDisplacement->blockMap().GID (j + nCompLocalDof);
                        UInt kGID = p2PatchDisplacement->blockMap().GID (j + 2 * nCompLocalDof);

                        //I think the coordinates are correct; what will be the bigger problem is to get right GID to fill the epetra p2PatchDisplacement EpetraVector
                        Vector3D coordinates;
                        coordinates(0) = p2PositionVector[iGID];
                        coordinates(1) = p2PositionVector[jGID];
                        coordinates(2) = p2PositionVector[kGID];


                        //std::cout << "THESE ARE Coordinates OF NODES: " << coordinates(0) << "\t" << coordinates(1) << "\t" << coordinates(2) << std::endl;
                        //we have here now the coordinates of of localDOF j; we call this point P

                        Vector3D QP; //define here the vector that goes from Q (point on plane) to point P


                        QP = coordinates - current_point_on_plane;

                     //   QP[0] = coordinates[0] - current_point_on_plane[0];
                     //   QP[1] = coordinates[1] - current_point_on_plane[1];
                     //   QP[2] = coordinates[2] - current_point_on_plane[2];


                        //next we do the projection of the vector QP ond the normalvector to the plane; we look at the sign which it has in the if statement

                        //if(QP[0]*direction[0] + QP[1]*direction[1] + QP[2]*direction[2] <=0) //we use here normalised normal vector
                        if(QP.dot(direction) < 0) //here i have change <= to <
                        {
                            //if the dot product is smaller equal zero, then we want to apply a displacement; for that we calculate the distance from point P to plane and then say this is the displacement we want

                            distance = abs(QP[0]*direction[0]+ QP[1]*direction[1]+QP[2]*direction[2]);
                            //std::cout << " THIS IS CALCULATED DISTANCE IN VECTORFIELD: " << distance << std::endl;


                            //////////////////Here we write distance to file
                            auto currentprocessor = dFeSpace->mesh()->comm()->MyPID();
                            std::ostringstream oss;//this is to convert int to string
                            oss << currentprocessor;

                            std::string path = "/cluster/home/pamstad/LIFE5/lifev-env/lifev-em-build/lifev/em/examples/example_EMHeart/distancefiles/distances_" + oss.str() + ".dat";
                            std::ofstream writer(path.c_str(), std::ios_base::app);
                            if(writer.is_open()==false)
                            {
                                std::cout << "error occured while opening the file" << std::endl;
                            }
                            writer << distance << std::endl;
                            writer.close();
                            //////////////////Here the writing to the file ends


                        }
                        else
                        {
                            distance = 0;
                        }

                        Vector3D displacement_vector;
                        displacement_vector[0] = distance*direction[0];
                        displacement_vector[1] = distance*direction[1];
                        displacement_vector[2] = distance*direction[2];


                        (*p2PatchDisplacement)[iGID] = displacement_vector[0];
                        (*p2PatchDisplacement)[jGID] = displacement_vector[1];
                        (*p2PatchDisplacement)[kGID] = displacement_vector[2];

            }


            return p2PatchDisplacement;

}
*/

//REGISTER(EssentialPatchBC, EssentialPatchBCMovingPlane);
//}//this is Klammer von LifeV namespace

};

REGISTER(EssentialPatchBC, EssentialPatchBCRotatingPlane);

}
#endif /* EssentialPatchBCRotatingPlane_hpp */
