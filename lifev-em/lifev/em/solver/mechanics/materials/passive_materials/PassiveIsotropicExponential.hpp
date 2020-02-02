/*
 * EMPassiveNeoHookean.hpp
 *
 *  Created on: 29/apr/2014
 *      Author: srossi
 */

#ifndef PASSIVEISOTROPICEXPONENTIAL_HPP_
#define PASSIVEISOTROPICEXPONENTIAL_HPP_


#include <lifev/em/solver/mechanics/materials/EMPassiveMaterialType.hpp>

//#include <lifev/em/solver/mechanics/materials/EMMaterialFunctions.hpp>



//class EMPassiveNeoHookean : public virtual EMPassiveMaterial
//{
//public:
//  EMPassiveNeoHookean();
//  virtual ~EMPassiveNeoHookean() {}
//};

namespace LifeV
{

template<typename Mesh>
class PassiveIsotropicExponential : public virtual EMPassiveMaterialType<Mesh>
{
public:
    typedef EMMaterialType<Mesh> super;

    PassiveIsotropicExponential();
    virtual ~PassiveIsotropicExponential() {}
};

template<typename Mesh>
PassiveIsotropicExponential<Mesh>::PassiveIsotropicExponential() :
    super ("Passive Isotropic Exponential", 4)
{
    this -> M_materialFunctionList[0].reset (new MaterialFunctions::Volumetric<Mesh> (3500000.0)  );
    this -> M_materialFunctionList[1].reset (new MaterialFunctions::dVolumetric<Mesh> (3500000.0) );
    this -> M_materialFunctionList[2].reset (new MaterialFunctions::IsotropicExponential<Mesh>()  );
    this -> M_materialFunctionList[3].reset (new MaterialFunctions::dIsotropicExponential<Mesh>()  );
}


template <typename MeshType>
inline EMPassiveMaterialType<MeshType>* createPassiveIsotropicExponential()
{
    return new PassiveIsotropicExponential<MeshType>();
}
namespace
{
static bool registerEM_passiveIE = EMPassiveMaterialType<LifeV::RegionMesh<LinearTetra> >::EMPassiveMaterialFactory::instance().registerProduct ("PIE", &createPassiveIsotropicExponential<LifeV::RegionMesh<LinearTetra> > );
}

}//LifeV

#endif /* EMPASSIVENEOHOOKEAN_HPP_ */
