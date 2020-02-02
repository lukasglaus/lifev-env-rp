/*
 * MaterialsList.hpp
 *
 *  Created on: 30/apr/2014
 *      Author: srossi
 */

#ifndef MATERIALSLIST_HPP_
#define MATERIALSLIST_HPP_


#include <lifev/em/solver/mechanics/materials/passive_materials/PassiveHolzapfelOgden.hpp>
#include <lifev/em/solver/mechanics/materials/passive_materials/PassiveOrthotropicFung.hpp>
#include <lifev/em/solver/mechanics/materials/passive_materials/PassiveNeoHookean.hpp>
#include <lifev/em/solver/mechanics/materials/passive_materials/PassiveIsotropicExponential.hpp>
#include <lifev/em/solver/mechanics/materials/passive_materials/PassiveIsotropicExponentialWithShear.hpp>
#include <lifev/em/solver/mechanics/materials/passive_materials/PassiveIsotropicFung.hpp>
#include <lifev/em/solver/mechanics/materials/passive_materials/PassiveLinearizedNeoHookean.hpp>
#include <lifev/em/solver/mechanics/materials/passive_materials/PassiveVolumetric.hpp>
#include <lifev/em/solver/mechanics/materials/passive_materials/PassiveTransverselyIsotropicExponential.hpp>
#include <lifev/em/solver/mechanics/materials/passive_materials/PassiveTransverselyIsotropicFung.hpp>
#include <lifev/em/solver/mechanics/materials/passive_materials/PassiveMooneyRivlin.hpp>
#include <lifev/em/solver/mechanics/materials/passive_materials/PassiveMooneyRivlinCompact.hpp>

// ACTIVE STRESS MATERIAL
#include <lifev/em/solver/mechanics/materials/active_stress_materials/SimpleFibersActiveStress.hpp>
#include <lifev/em/solver/mechanics/materials/active_stress_materials/OrthotropicFibersActiveStress.hpp>

// ACTIVE STRAIN MATERIAL
#include <lifev/em/solver/mechanics/materials/active_strain_materials/ActiveNeoHookean.hpp>
#include <lifev/em/solver/mechanics/materials/active_strain_materials/ActiveIsotropicExponential.hpp>
#include <lifev/em/solver/mechanics/materials/active_strain_materials/ActiveHolzapfelOgden.hpp>
#include <lifev/em/solver/mechanics/materials/active_strain_materials/ActiveHolzapfelOgdenNIIsotropicSplitting.hpp>

#endif /* MATERIALSLIST_HPP_ */
