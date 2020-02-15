/*
 * EssentialPatchBCRotatingPlane.h
 *
 *  Created on: Feb 5, 2020
 *      Author: lglaus
 */

#ifndef LIFEV_EM_EXAMPLES_EXAMPLE_EMHEART_ESSENTIALPATCHBCROTATINGPLANE_H_
#define LIFEV_EM_EXAMPLES_EXAMPLE_EMHEART_ESSENTIALPATCHBCROTATINGPLANE_H_

#include <stdio.h>
#include <lifev/em/examples/example_EMHeart/EssentialPatchBC.hpp>
#include <lifev/em/examples/example_EMHeart/GenericFactory.hpp>
#include <vector>

#define PI 3.14159265359

namespace LifeV
{

class EssentialPatchBCRotatingPlane : public EssentialPatchBC {
public:

    int rotation_direction;
    Real angleOfTime;
    

    virtual void setup(const GetPot& dataFile, const std::string& name,EMSolver<RegionMesh<LinearTetra>, EMMonodomainSolver<RegionMesh<LinearTetra> > >& solver); //In the setup function the basic things come in it, like as name of patch, flag, direction vector, displacement vector and so on

    virtual const bool nodeOnPatch(const Vector3D& coord, const Real& time);

    virtual void modifyPatchArea(EMSolver<RegionMesh<LinearTetra>, EMMonodomainSolver<RegionMesh<LinearTetra> > >& solver,const int& newFlag, const Real& time);

    virtual vectorPtr_Type directionalVectorField (EMSolver<RegionMesh<LinearTetra>, EMMonodomainSolver<RegionMesh<LinearTetra> > >& solver,const boost::shared_ptr<FESpace<RegionMesh<LinearTetra>, MapEpetra >> dFeSpace, Vector3D& direction, const Real& disp, const Real& time);

    virtual const bool nodeOnPatchCurrent(const Vector3D& coord,const Real& time);

    virtual void modifyPatchBC(EMSolver<RegionMesh<LinearTetra>, EMMonodomainSolver<RegionMesh<LinearTetra> > >& solver, const Real& time, int& PatchFlag);
    
    virtual vector_Type displayDirectionalVectorField(EMSolver<RegionMesh<LinearTetra>, EMMonodomainSolver<RegionMesh<LinearTetra> > >& solver, const Real& time);
    
private:
    //vectorPtr_Type p2PositionVector;
    //LifeV::EssentialPatchBC::vector_Type p2PositionVector;
    //vector_Type p2PositionVector;
    vectorPtr_Type m_p2currentPositionVector;
    Real m_maxDisplacement;
    Vector3D normal_vector;
    Vector3D starting_point;
    Vector3D direction_to_axis; //2020.02.10 lg
    Vector3D axis_direction; //2020.02.10 lg
    Real distance_to_axis; //2020.02.10 lg
    Real maximum_angle; //2020.02.10 lg
    Real minimum_angle; //2020.02.10 lg
    Vector3D pointOnHeart;
    Real m_tduration;

    
    
};

REGISTER(EssentialPatchBC, EssentialPatchBCRotatingPlane);

}
#endif /* LIFEV_EM_EXAMPLES_EXAMPLE_EMHEART_ESSENTIALPATCHBCMOVINGPLANE_H_ */
