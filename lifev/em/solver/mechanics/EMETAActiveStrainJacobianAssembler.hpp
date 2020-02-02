/*
 * EMETAActiveStrainJaconbianAssembler.hpp
 *
 *  Created on: 29/apr/2014
 *      Author: srossi, Luca Barbarotta
 */

#ifndef EMETAACTIVESTRAINJACOBIANASSMEBLER_HPP_
#define EMETAACTIVESTRAINJACOBIANASSMEBLER_HPP_


#include <lifev/core/fem/FESpace.hpp>
#include <lifev/core/array/VectorEpetra.hpp>
#include <lifev/core/array/MatrixEpetra.hpp>

//ET include for assemblings
#include <lifev/eta/fem/ETFESpace.hpp>
#include <lifev/eta/expression/Integrate.hpp>

#include <lifev/em/solver/mechanics/EMMechanicalExpressions.hpp>
//#include <lifev/em/solver/mechanics/materials/EMMaterialFunctions.hpp>

#include <lifev/em/util/EMUtility.hpp>
#include <lifev/em/solver/EMETAFunctors.hpp>

//#include <boost/typeof/typeof.hpp>

namespace LifeV
{

    typedef VectorEpetra           vector_Type;
    typedef boost::shared_ptr<vector_Type>         vectorPtr_Type;

    typedef MatrixEpetra<Real>           matrix_Type;
    typedef boost::shared_ptr<matrix_Type>         matrixPtr_Type;


    namespace EMAssembler
    {
        namespace ActiveStrain
        {
            template <typename Mesh, typename FunctorPtr >
            void
            computeActiveStrainI1JacobianTerms (    const vector_Type& disp,
                                                    boost::shared_ptr<ETFESpace<Mesh, MapEpetra, 3, 3 > > dispETFESpace,
                                                    const vector_Type& fibers,
                                                    const vector_Type& sheets,
                                                    const vectorPtr_Type& gammaf,
                                                    const vectorPtr_Type& gammas,
                                                    const vectorPtr_Type& gamman,
                                                    boost::shared_ptr<ETFESpace<Mesh, MapEpetra, 3, 1 > >  activationETFESpace,
                                                    matrixPtr_Type           jacobianPtr,
                                                    FunctorPtr               W1,
                                                    Real orthotropicParameter = -666.)
            {
                //
                if(disp.comm().MyPID() == 0)
                    std::cout << "EMETA - Computing Isotropic Active Strain jacobian terms: ";

                using namespace ExpressionAssembly;

                auto I = _I;
                auto dF = grad(phi_j);
//                auto GradU = _Grad_u(dispETFESpace, disp, 0);
//                auto F = I + GradU;
                auto F = _F (dispETFESpace, disp, 0);

                auto f_0 = _v0 (dispETFESpace, fibers);
                auto s_0 = _v0 (dispETFESpace, sheets);

                boost::shared_ptr<orthonormalizeFibers> normalize0 (new orthonormalizeFibers);
                auto f0 = eval (normalize0, f_0);

                boost::shared_ptr<orthonormalizeFibers> normalize1 (new orthonormalizeFibers (1) );
                auto s0 = eval (normalize1, f0, s_0);

                boost::shared_ptr<CrossProduct> wedge (new CrossProduct);
                auto n0 = eval ( wedge, f0, s0);


                if(gammas && gamman)
                {
                    if(disp.comm().MyPID() == 0)
                        std::cout << " Anisotropic case ... \n";

                    auto gf = value (activationETFESpace, *gammaf);
                    auto gs = value (activationETFESpace, *gammas);
                    auto gn = value (activationETFESpace, *gamman);


                    auto FAinv = _FAinv(gf, gs, gn, f0, s0, n0);
                    auto FE =  F * FAinv;
                    auto dFE = _dFE(FAinv);

                    auto dP = eval (W1, FE ) * (_d2I1dF (dFE) ) * FAinv;

                    integrate ( elements ( dispETFESpace->mesh() ) ,
                                quadRuleTetra4pt,
                                dispETFESpace,
                                dispETFESpace,
                                dot ( dP , grad (phi_i) )
                        ) >> jacobianPtr;
                }
                else
                {
                    auto gf = value (activationETFESpace, *gammaf);

                    if(orthotropicParameter > 0 )
                    {
                        if(disp.comm().MyPID() == 0)
                            std::cout << " Orthotropic case ... \n";

                        auto k = value(orthotropicParameter);
                        auto FAinv = _FAinv(gf, k, f0, s0, n0);

                        auto FE =  F * FAinv;
                        auto dFE = _dFE(FAinv);
                        auto dP = eval (W1, FE ) * (_d2I1dF (dFE) ) * FAinv;

                        integrate ( elements ( dispETFESpace->mesh() ) ,
                                    quadRuleTetra4pt,
                                    dispETFESpace,
                                    dispETFESpace,
                                    dot ( dP , grad (phi_i) )
                            ) >> jacobianPtr;
                    }
                    else
                    {
                        if(disp.comm().MyPID() == 0)
                            std::cout << " Transversely isotropic case ... \n";

//    		auto FAinv = _FAinv(gf, f0, s0, n0);
                        auto FAinv = _FAinv(gf, f0);


                        auto FE =  F * FAinv;
                        auto dFE = _dFE(FAinv);
                        auto dP = eval (W1, FE ) * (_d2I1dF (dFE) ) * FAinv;

//    		auto dFE = dF * FAinv;
//    		auto dFET = FAinv * transpose(dF)
//
//    		auto mu = value(1000.);
//
//    		auto J = det(F);
//    		auto Jm23 = pow(J, -2./3.);
//    		auto FmT = minusT(F);
//    		auto FEmT = minusT(FE);
//    		auto dF = _dF;
//    		auto I1E = dot(FE, FE);
//
//    		auto varPE = mu * (FE + value(-1./3.) * I1E * FEmT);
//    		auto P = varPE * FAinv;
//
//    		auto dJm23dF = value(-2./3.) * Jm23 * dot( FmT, dF );
//
//    		auto dFEmTdFE = value(-1.) * FEmT * dFET * FEmT;
//    		auto dI1EdFE  = value(2.) * dot(FE, dFE);
//
//    		auto dP = dJm23dF * P + Jm23 * ( dFE + value(-1./3.) * ( dI1EdFE * FEmT + I1E * dFEmTdFE )  ) * FAinv;
//

                        integrate ( elements ( dispETFESpace->mesh() ) ,
                                    quadRuleTetra4pt,
                                    dispETFESpace,
                                    dispETFESpace,
                                    dot ( dP , grad (phi_i) )
                            ) >> jacobianPtr;
                    }
                }

            }
        
            template <typename Mesh, typename FunctorPtr >
            void
            computeActiveStrainI1JacobianTermsSecondDerivative (    const vector_Type& disp,
                                                                    boost::shared_ptr<ETFESpace<Mesh, MapEpetra, 3, 3 > > dispETFESpace,
                                                                    const vector_Type& fibers,
                                                                    const vector_Type& sheets,
                                                                    const vectorPtr_Type& gammaf,
                                                                    const vectorPtr_Type& gammas,
                                                                    const vectorPtr_Type& gamman,
                                                                    boost::shared_ptr<ETFESpace<Mesh, MapEpetra, 3, 1 > >  activationETFESpace,
                                                                    matrixPtr_Type           jacobianPtr,
                                                                    FunctorPtr               d2W,
                                                                    Real orthotropicParameter = -666.)
            {
                //
                if(disp.comm().MyPID() == 0)
                    std::cout << "EMETA - Computing Isotropic Active Strain jacobian second derivative terms: ";

                using namespace ExpressionAssembly;

                auto I = _I;
                auto dF = grad(phi_j);
//                auto GradU = _Grad_u(dispETFESpace, disp, 0);
//                auto F = I + GradU;
                auto F = _F (dispETFESpace, disp, 0);

                auto f_0 = _v0 (dispETFESpace, fibers);
                auto s_0 = _v0 (dispETFESpace, sheets);

                boost::shared_ptr<orthonormalizeFibers> normalize0 (new orthonormalizeFibers);
                auto f0 = eval (normalize0, f_0);

                boost::shared_ptr<orthonormalizeFibers> normalize1 (new orthonormalizeFibers (1) );
                auto s0 = eval (normalize1, f0, s_0);

                boost::shared_ptr<CrossProduct> wedge (new CrossProduct);
                auto n0 = eval ( wedge, f0, s0);


                if(gammas && gamman)
                {
                    if(disp.comm().MyPID() == 0)
                        std::cout << " Anisotropic case ... \n";

                    auto gf = value (activationETFESpace, *gammaf);
                    auto gs = value (activationETFESpace, *gammas);
                    auto gn = value (activationETFESpace, *gamman);


                    auto FAinv = _FAinv(gf, gs, gn, f0, s0, n0);
                    auto FE =  F * FAinv;
                    auto dFE = _dFE(FAinv);
                    auto dP = eval (d2W, FE ) * (_dI1dF (FE, dFE) ) * _dI1(FE) * FAinv;
		
                    integrate ( elements ( dispETFESpace->mesh() ) ,
                                quadRuleTetra4pt,
                                dispETFESpace,
                                dispETFESpace,
                                dot ( dP , grad (phi_i) )
                        ) >> jacobianPtr;
                }
                else
                {
                    auto gf = value (activationETFESpace, *gammaf);

                    if(orthotropicParameter > 0 )
                    {
                        if(disp.comm().MyPID() == 0)
                            std::cout << " Orthotropic case ... \n";

                        auto k = value(orthotropicParameter);
                        auto FAinv = _FAinv(gf, k, f0, s0, n0);

                        auto FE =  F * FAinv;
                        auto dFE = _dFE(FAinv);
                        auto dP = eval (d2W, FE ) * (_dI1dF (FE, dFE) ) * _dI1(FE) * FAinv;

                        integrate ( elements ( dispETFESpace->mesh() ) ,
                                    quadRuleTetra4pt,
                                    dispETFESpace,
                                    dispETFESpace,
                                    dot ( dP , grad (phi_i) )
                            ) >> jacobianPtr;
                    }
                    else
                    {
                        if(disp.comm().MyPID() == 0)
                            std::cout << " Transversely isotropic case ... \n";

//    		auto FAinv = _FAinv(gf, f0, s0, n0);
                        auto FAinv = _FAinv(gf, f0);


                        auto FE =  F * FAinv;
                        auto dFE = _dFE(FAinv);
                        auto dP = eval (d2W, FE ) * (_dI1dF (FE, dFE) ) * _dI1(FE) * FAinv;

//    		auto dFE = dF * FAinv;
//    		auto dFET = FAinv * transpose(dF)
//
//    		auto mu = value(1000.);
//
//    		auto J = det(F);
//    		auto Jm23 = pow(J, -2./3.);
//    		auto FmT = minusT(F);
//    		auto FEmT = minusT(FE);
//    		auto dF = _dF;
//    		auto I1E = dot(FE, FE);
//
//    		auto varPE = mu * (FE + value(-1./3.) * I1E * FEmT);
//    		auto P = varPE * FAinv;
//
//    		auto dJm23dF = value(-2./3.) * Jm23 * dot( FmT, dF );
//
//    		auto dFEmTdFE = value(-1.) * FEmT * dFET * FEmT;
//    		auto dI1EdFE  = value(2.) * dot(FE, dFE);
//
//    		auto dP = dJm23dF * P + Jm23 * ( dFE + value(-1./3.) * ( dI1EdFE * FEmT + I1E * dFEmTdFE )  ) * FAinv;
//

                        integrate ( elements ( dispETFESpace->mesh() ) ,
                                    quadRuleTetra4pt,
                                    dispETFESpace,
                                    dispETFESpace,
                                    dot ( dP , grad (phi_i) )
                            ) >> jacobianPtr;
                    }
                }

            }
        
            template <typename Mesh, typename FunctorPtr >
            void
            computeActiveStrainI4FibersJacobianTerms (    const vector_Type& disp,
                                                          boost::shared_ptr<ETFESpace<Mesh, MapEpetra, 3, 3 > > dispETFESpace,
                                                          const vector_Type& fibers,
                                                          const vector_Type& sheets,
                                                          const vectorPtr_Type& gammaf,
                                                          const vectorPtr_Type& gammas,
                                                          const vectorPtr_Type& gamman,
                                                          boost::shared_ptr<ETFESpace<Mesh, MapEpetra, 3, 1 > >  activationETFESpace,
                                                          matrixPtr_Type           jacobianPtr,
                                                          FunctorPtr               W4f,
                                                          Real orthotropicParameter = -666.)
            {
                //
                if(disp.comm().MyPID() == 0)
                    std::cout << "EMETA - Computing Anisotropic Fibers Active Strain jacobian terms: ";

                using namespace ExpressionAssembly;

                auto I = _I;
                auto dF = grad(phi_j);
                auto GradU = _Grad_u(dispETFESpace, disp, 0);
                auto F = I + GradU;

                auto f_0 = _v0 (dispETFESpace, fibers);
                auto s_0 = _v0 (dispETFESpace, sheets);

                boost::shared_ptr<orthonormalizeFibers> normalize0 (new orthonormalizeFibers);
                auto f0 = eval (normalize0, f_0);

                boost::shared_ptr<orthonormalizeFibers> normalize1 (new orthonormalizeFibers (1) );
                auto s0 = eval (normalize1, f0, s_0);

                boost::shared_ptr<CrossProduct> wedge (new CrossProduct);
                auto n0 = eval ( wedge, f0, s0);


                if(gammas && gamman)
                {
                    if(disp.comm().MyPID() == 0)
                        std::cout << " Anisotropic case ... \n";

                    auto gf = value (activationETFESpace, *gammaf);
                    auto gs = value (activationETFESpace, *gammas);
                    auto gn = value (activationETFESpace, *gamman);


                    auto FAinv = _FAinv(gf, gs, gn, f0, s0, n0);
                    auto FE =  F * FAinv;
                    auto dFE = _dFE(FAinv);

                    auto dP = eval (W4f, FE ) * (_d2I4dF (f0, dFE) ) * FAinv;

                    integrate ( elements ( dispETFESpace->mesh() ) ,
                                quadRuleTetra4pt,
                                dispETFESpace,
                                dispETFESpace,
                                dot ( dP , grad (phi_i) )
                        ) >> jacobianPtr;
                }
                else
                {
                    auto gf = value (activationETFESpace, *gammaf);

                    if(orthotropicParameter > 0 )
                    {
                        if(disp.comm().MyPID() == 0)
                            std::cout << " Orthotropic case ... \n";

                        auto k = value(orthotropicParameter);
                        auto FAinv = _FAinv(gf, k, f0, s0, n0);

                        auto FE =  F * FAinv;
                        auto dFE = _dFE(FAinv);
                        auto dP = eval (W4f, FE ) * (_d2I4dF (f0, dFE) ) * FAinv;

                        integrate ( elements ( dispETFESpace->mesh() ) ,
                                    quadRuleTetra4pt,
                                    dispETFESpace,
                                    dispETFESpace,
                                    dot ( dP , grad (phi_i) )
                            ) >> jacobianPtr;
                    }
                    else
                    {
                        if(disp.comm().MyPID() == 0)
                            std::cout << " Transversely isotropic case ... \n";

                        auto FAinv = _FAinv(gf, f0);

                        auto FE =  F * FAinv;
                        auto dFE = _dFE(FAinv);
                        auto dP = eval (W4f, FE ) * (_d2I4dF (f0, dFE) ) * FAinv;
                    
                        integrate ( elements ( dispETFESpace->mesh() ) ,
                                    quadRuleTetra4pt,
                                    dispETFESpace,
                                    dispETFESpace,
                                    dot ( dP , grad (phi_i) )
                            ) >> jacobianPtr;
                    }
                }

            }
        
            template <typename Mesh, typename FunctorPtr >
            void
            computeActiveStrainI4FibersJacobianTermsSecondDerivative (    const vector_Type& disp,
                                                                          boost::shared_ptr<ETFESpace<Mesh, MapEpetra, 3, 3 > > dispETFESpace,
                                                                          const vector_Type& fibers,
                                                                          const vector_Type& sheets,
                                                                          const vectorPtr_Type& gammaf,
                                                                          const vectorPtr_Type& gammas,
                                                                          const vectorPtr_Type& gamman,
                                                                          boost::shared_ptr<ETFESpace<Mesh, MapEpetra, 3, 1 > >  activationETFESpace,
                                                                          matrixPtr_Type           jacobianPtr,
                                                                          FunctorPtr               d2W4f,
                                                                          Real orthotropicParameter = -666.)
            {
                //
                if(disp.comm().MyPID() == 0)
                    std::cout << "EMETA - Computing Anisotropic Fibers Active Strain jacobian second derivative terms: ";

                using namespace ExpressionAssembly;

                auto I = _I;
                auto dF = grad(phi_j);
                auto GradU = _Grad_u(dispETFESpace, disp, 0);
                auto F = I + GradU;

                auto f_0 = _v0 (dispETFESpace, fibers);
                auto s_0 = _v0 (dispETFESpace, sheets);

                boost::shared_ptr<orthonormalizeFibers> normalize0 (new orthonormalizeFibers);
                auto f0 = eval (normalize0, f_0);

                boost::shared_ptr<orthonormalizeFibers> normalize1 (new orthonormalizeFibers (1) );
                auto s0 = eval (normalize1, f0, s_0);

                boost::shared_ptr<CrossProduct> wedge (new CrossProduct);
                auto n0 = eval ( wedge, f0, s0);

                if(gammas && gamman)
                {
                    if(disp.comm().MyPID() == 0)
                        std::cout << " Anisotropic case ... \n";

                    auto gf = value (activationETFESpace, *gammaf);
                    auto gs = value (activationETFESpace, *gammas);
                    auto gn = value (activationETFESpace, *gamman);

                    auto FAinv = _FAinv(gf, gs, gn, f0, s0, n0);
                    auto FE =  F * FAinv;
                    auto dFE = _dFE(FAinv);
                    auto dP = eval (d2W4f, FE ) * (_dI4dF (FE, f0, dFE) ) * _dI4(FE, f0) * FAinv;
		
                    integrate ( elements ( dispETFESpace->mesh() ) ,
                                quadRuleTetra4pt,
                                dispETFESpace,
                                dispETFESpace,
                                dot ( dP , grad (phi_i) )
                        ) >> jacobianPtr;
                }
                else
                {
                    auto gf = value (activationETFESpace, *gammaf);

                    if(orthotropicParameter > 0 )
                    {
                        if(disp.comm().MyPID() == 0)
                            std::cout << " Orthotropic case ... \n";

                        auto k = value(orthotropicParameter);
                        auto FAinv = _FAinv(gf, k, f0, s0, n0);

                        auto FE =  F * FAinv;
                        auto dFE = _dFE(FAinv);
                        auto dP = eval (d2W4f, FE ) * (_dI4dF (FE, f0, dFE) ) * _dI4(FE, f0) * FAinv;

                        integrate ( elements ( dispETFESpace->mesh() ) ,
                                    quadRuleTetra4pt,
                                    dispETFESpace,
                                    dispETFESpace,
                                    dot ( dP , grad (phi_i) )
                            ) >> jacobianPtr;
                    }
                    else
                    {
                        if(disp.comm().MyPID() == 0)
                            std::cout << " Transversely isotropic case ... \n";

                        auto FAinv = _FAinv(gf, f0);


                        auto FE =  F * FAinv;
                        auto dFE = _dFE(FAinv);
                        auto dP = eval (d2W4f, FE ) * (_dI4dF (FE, f0, dFE) ) * _dI4(FE, f0) * FAinv;

                        integrate ( elements ( dispETFESpace->mesh() ) ,
                                    quadRuleTetra4pt,
                                    dispETFESpace,
                                    dispETFESpace,
                                    dot ( dP , grad (phi_i) )
                            ) >> jacobianPtr;
                    }
                }

            }

            template <typename Mesh, typename FunctorPtr >
            void
            computeActiveStrainI4SheetsJacobianTerms (    const vector_Type& disp,
                                                          boost::shared_ptr<ETFESpace<Mesh, MapEpetra, 3, 3 > > dispETFESpace,
                                                          const vector_Type& fibers,
                                                          const vector_Type& sheets,
                                                          const vectorPtr_Type& gammaf,
                                                          const vectorPtr_Type& gammas,
                                                          const vectorPtr_Type& gamman,
                                                          boost::shared_ptr<ETFESpace<Mesh, MapEpetra, 3, 1 > >  activationETFESpace,
                                                          matrixPtr_Type           jacobianPtr,
                                                          FunctorPtr               W4s,
                                                          Real orthotropicParameter = -666.)
            {
                //
                if(disp.comm().MyPID() == 0)
                    std::cout << "EMETA - Computing Anisotropic Sheets Active Strain jacobian terms: ";
            
                using namespace ExpressionAssembly;

                auto I = _I;
                auto dF = grad(phi_j);
                auto GradU = _Grad_u(dispETFESpace, disp, 0);
                auto F = I + GradU;

                auto f_0 = _v0 (dispETFESpace, fibers);
                auto s_0 = _v0 (dispETFESpace, sheets);

                boost::shared_ptr<orthonormalizeFibers> normalize0 (new orthonormalizeFibers);
                auto f0 = eval (normalize0, f_0);

                boost::shared_ptr<orthonormalizeFibers> normalize1 (new orthonormalizeFibers (1) );
                auto s0 = eval (normalize1, f0, s_0);

                boost::shared_ptr<CrossProduct> wedge (new CrossProduct);
                auto n0 = eval ( wedge, f0, s0);


                if(gammas && gamman)
                {
                    if(disp.comm().MyPID() == 0)
                        std::cout << " Anisotropic case ... \n";

                    auto gf = value (activationETFESpace, *gammaf);
                    auto gs = value (activationETFESpace, *gammas);
                    auto gn = value (activationETFESpace, *gamman);


                    auto FAinv = _FAinv(gf, gs, gn, f0, s0, n0);
                    auto FE =  F * FAinv;
                    auto dFE = _dFE(FAinv);
                    auto dP = eval (W4s, FE ) * (_d2I4dF (s0, dFE) ) * FAinv;

                    integrate ( elements ( dispETFESpace->mesh() ) ,
                                quadRuleTetra4pt,
                                dispETFESpace,
                                dispETFESpace,
                                dot ( dP , grad (phi_i) )
                        ) >> jacobianPtr;
                }
                else
                {
                    auto gf = value (activationETFESpace, *gammaf);

                    if(orthotropicParameter > 0 )
                    {
                        if(disp.comm().MyPID() == 0)
                            std::cout << " Orthotropic case ... \n";

                        auto k = value(orthotropicParameter);
                        auto FAinv = _FAinv(gf, k, f0, s0, n0);

                        auto FE =  F * FAinv;
                        auto dFE = _dFE(FAinv);
                        auto dP = eval (W4s, FE ) * (_d2I4dF (s0, dFE) ) * FAinv;

                        integrate ( elements ( dispETFESpace->mesh() ) ,
                                    quadRuleTetra4pt,
                                    dispETFESpace,
                                    dispETFESpace,
                                    dot ( dP , grad (phi_i) )
                            ) >> jacobianPtr;
                    }
                    else
                    {
                        if(disp.comm().MyPID() == 0)
                            std::cout << " Transversely isotropic case ... \n";

                        auto FAinv = _FAinv(gf, f0);

                        auto FE =  F * FAinv;
                        auto dFE = _dFE(FAinv);
                        auto dP = eval (W4s, FE ) * (_d2I4dF (s0, dFE) ) * FAinv;
                    
                        integrate ( elements ( dispETFESpace->mesh() ) ,
                                    quadRuleTetra4pt,
                                    dispETFESpace,
                                    dispETFESpace,
                                    dot ( dP , grad (phi_i) )
                            ) >> jacobianPtr;
                    }
                }

            }
        
            template <typename Mesh, typename FunctorPtr >
            void
            computeActiveStrainI4SheetsJacobianTermsSecondDerivative (    const vector_Type& disp,
                                                                          boost::shared_ptr<ETFESpace<Mesh, MapEpetra, 3, 3 > > dispETFESpace,
                                                                          const vector_Type& fibers,
                                                                          const vector_Type& sheets,
                                                                          const vectorPtr_Type& gammaf,
                                                                          const vectorPtr_Type& gammas,
                                                                          const vectorPtr_Type& gamman,
                                                                          boost::shared_ptr<ETFESpace<Mesh, MapEpetra, 3, 1 > >  activationETFESpace,
                                                                          matrixPtr_Type           jacobianPtr,
                                                                          FunctorPtr               d2W4s,
                                                                          Real orthotropicParameter = -666.)
            {
                //
                if(disp.comm().MyPID() == 0)
                    std::cout << "EMETA - Computing Anisotropic Sheets Active Strain jacobian second derivative terms: ";

                using namespace ExpressionAssembly;

                auto I = _I;
                auto dF = grad(phi_j);
                auto GradU = _Grad_u(dispETFESpace, disp, 0);
                auto F = I + GradU;

                auto f_0 = _v0 (dispETFESpace, fibers);
                auto s_0 = _v0 (dispETFESpace, sheets);

                boost::shared_ptr<orthonormalizeFibers> normalize0 (new orthonormalizeFibers);
                auto f0 = eval (normalize0, f_0);

                boost::shared_ptr<orthonormalizeFibers> normalize1 (new orthonormalizeFibers (1) );
                auto s0 = eval (normalize1, f0, s_0);

                boost::shared_ptr<CrossProduct> wedge (new CrossProduct);
                auto n0 = eval ( wedge, f0, s0);

                if(gammas && gamman)
                {
                    if(disp.comm().MyPID() == 0)
                        std::cout << " Anisotropic case ... \n";

                    auto gf = value (activationETFESpace, *gammaf);
                    auto gs = value (activationETFESpace, *gammas);
                    auto gn = value (activationETFESpace, *gamman);

                    auto FAinv = _FAinv(gf, gs, gn, f0, s0, n0);
                    auto FE =  F * FAinv;
                    auto dFE = _dFE(FAinv);
                    auto dP = eval (d2W4s, FE ) * (_dI4dF (FE, s0, dFE) ) * _dI4(FE, s0) * FAinv;
		
                    integrate ( elements ( dispETFESpace->mesh() ) ,
                                quadRuleTetra4pt,
                                dispETFESpace,
                                dispETFESpace,
                                dot ( dP , grad (phi_i) )
                        ) >> jacobianPtr;
                }
                else
                {
                    auto gf = value (activationETFESpace, *gammaf);

                    if(orthotropicParameter > 0 )
                    {
                        if(disp.comm().MyPID() == 0)
                            std::cout << " Orthotropic case ... \n";

                        auto k = value(orthotropicParameter);
                        auto FAinv = _FAinv(gf, k, f0, s0, n0);

                        auto FE =  F * FAinv;
                        auto dFE = _dFE(FAinv);
                        auto dP = eval (d2W4s, FE ) * (_dI4dF (FE, s0, dFE) ) * _dI4(FE, s0) * FAinv;

                        integrate ( elements ( dispETFESpace->mesh() ) ,
                                    quadRuleTetra4pt,
                                    dispETFESpace,
                                    dispETFESpace,
                                    dot ( dP , grad (phi_i) )
                            ) >> jacobianPtr;
                    }
                    else
                    {
                        if(disp.comm().MyPID() == 0)
                            std::cout << " Transversely isotropic case ... \n";

                        auto FAinv = _FAinv(gf, f0);


                        auto FE =  F * FAinv;
                        auto dFE = _dFE(FAinv);
                        auto dP = eval (d2W4s, FE ) * (_dI4dF (FE, s0, dFE) ) * _dI4(FE, s0) * FAinv;

                        integrate ( elements ( dispETFESpace->mesh() ) ,
                                    quadRuleTetra4pt,
                                    dispETFESpace,
                                    dispETFESpace,
                                    dot ( dP , grad (phi_i) )
                            ) >> jacobianPtr;
                    }
                }

            }
        
            template <typename Mesh, typename FunctorPtr >
            void
            computeActiveStrainI8JacobianTerms (    const vector_Type& disp,
                                                    boost::shared_ptr<ETFESpace<Mesh, MapEpetra, 3, 3 > > dispETFESpace,
                                                    const vector_Type& fibers,
                                                    const vector_Type& sheets,
                                                    const vectorPtr_Type& gammaf,
                                                    const vectorPtr_Type& gammas,
                                                    const vectorPtr_Type& gamman,
                                                    boost::shared_ptr<ETFESpace<Mesh, MapEpetra, 3, 1 > >  activationETFESpace,
                                                    matrixPtr_Type           jacobianPtr,
                                                    FunctorPtr               W8fs,
                                                    Real orthotropicParameter = -666.)
            {
                //
                if(disp.comm().MyPID() == 0)
                    std::cout << "EMETA - Computing Anisotropic Shear Active Strain jacobian terms: ";

                using namespace ExpressionAssembly;

                auto I = _I;
                auto dF = grad(phi_j);
                auto GradU = _Grad_u(dispETFESpace, disp, 0);
                auto F = I + GradU;

                auto f_0 = _v0 (dispETFESpace, fibers);
                auto s_0 = _v0 (dispETFESpace, sheets);

                boost::shared_ptr<orthonormalizeFibers> normalize0 (new orthonormalizeFibers);
                auto f0 = eval (normalize0, f_0);

                boost::shared_ptr<orthonormalizeFibers> normalize1 (new orthonormalizeFibers (1) );
                auto s0 = eval (normalize1, f0, s_0);

                boost::shared_ptr<CrossProduct> wedge (new CrossProduct);
                auto n0 = eval ( wedge, f0, s0);


                if(gammas && gamman)
                {
                    if(disp.comm().MyPID() == 0)
                        std::cout << " Anisotropic case ... \n";

                    auto gf = value (activationETFESpace, *gammaf);
                    auto gs = value (activationETFESpace, *gammas);
                    auto gn = value (activationETFESpace, *gamman);


                    auto FAinv = _FAinv(gf, gs, gn, f0, s0, n0);
                    auto FE =  F * FAinv;
                    auto dFE = _dFE(FAinv);

                    auto dP = eval (W8fs, FE ) * (_d2I8dF (f0, s0, dFE) ) * FAinv;

                    integrate ( elements ( dispETFESpace->mesh() ) ,
                                quadRuleTetra4pt,
                                dispETFESpace,
                                dispETFESpace,
                                dot ( dP , grad (phi_i) )
                        ) >> jacobianPtr;
                }
                else
                {
                    auto gf = value (activationETFESpace, *gammaf);

                    if(orthotropicParameter > 0 )
                    {
                        if(disp.comm().MyPID() == 0)
                            std::cout << " Orthotropic case ... \n";

                        auto k = value(orthotropicParameter);
                        auto FAinv = _FAinv(gf, k, f0, s0, n0);

                        auto FE =  F * FAinv;
                        auto dFE = _dFE(FAinv);
                        auto dP = eval (W8fs, FE ) * (_d2I8dF (f0, s0, dFE) ) * FAinv;

                        integrate ( elements ( dispETFESpace->mesh() ) ,
                                    quadRuleTetra4pt,
                                    dispETFESpace,
                                    dispETFESpace,
                                    dot ( dP , grad (phi_i) )
                            ) >> jacobianPtr;
                    }
                    else
                    {
                        if(disp.comm().MyPID() == 0)
                            std::cout << " Transversely isotropic case ... \n";

                        auto FAinv = _FAinv(gf, f0);


                        auto FE =  F * FAinv;
                        auto dFE = _dFE(FAinv);
                        auto dP = eval (W8fs, FE ) * (_d2I8dF (f0, s0, dFE) ) * FAinv;

                        integrate ( elements ( dispETFESpace->mesh() ) ,
                                    quadRuleTetra4pt,
                                    dispETFESpace,
                                    dispETFESpace,
                                    dot ( dP , grad (phi_i) )
                            ) >> jacobianPtr;
                    }
                }

            }
        
            template <typename Mesh, typename FunctorPtr >
            void
            computeActiveStrainI8JacobianTermsSecondDerivative (    const vector_Type& disp,
                                                                    boost::shared_ptr<ETFESpace<Mesh, MapEpetra, 3, 3 > > dispETFESpace,
                                                                    const vector_Type& fibers,
                                                                    const vector_Type& sheets,
                                                                    const vectorPtr_Type& gammaf,
                                                                    const vectorPtr_Type& gammas,
                                                                    const vectorPtr_Type& gamman,
                                                                    boost::shared_ptr<ETFESpace<Mesh, MapEpetra, 3, 1 > >  activationETFESpace,
                                                                    matrixPtr_Type           jacobianPtr,
                                                                    FunctorPtr               d2W8fs,
                                                                    Real orthotropicParameter = -666.)
            {
                //
                if(disp.comm().MyPID() == 0)
                    std::cout << "EMETA - Computing Anisotropic Shear Active Strain jacobian second derivative terms: ";

                using namespace ExpressionAssembly;

                auto I = _I;
                auto dF = grad(phi_j);
                auto GradU = _Grad_u(dispETFESpace, disp, 0);
                auto F = I + GradU;

                auto f_0 = _v0 (dispETFESpace, fibers);
                auto s_0 = _v0 (dispETFESpace, sheets);

                boost::shared_ptr<orthonormalizeFibers> normalize0 (new orthonormalizeFibers);
                auto f0 = eval (normalize0, f_0);

                boost::shared_ptr<orthonormalizeFibers> normalize1 (new orthonormalizeFibers (1) );
                auto s0 = eval (normalize1, f0, s_0);

                boost::shared_ptr<CrossProduct> wedge (new CrossProduct);
                auto n0 = eval ( wedge, f0, s0);


                if(gammas && gamman)
                {
                    if(disp.comm().MyPID() == 0)
                        std::cout << " Anisotropic case ... \n";

                    auto gf = value (activationETFESpace, *gammaf);
                    auto gs = value (activationETFESpace, *gammas);
                    auto gn = value (activationETFESpace, *gamman);

                    auto FAinv = _FAinv(gf, gs, gn, f0, s0, n0);
                    auto FE =  F * FAinv;
                    auto dFE = _dFE(FAinv);
                    auto dP = eval (d2W8fs, FE ) * (_dI8dF (FE, f0, s0, dFE) ) * _dI8(FE, f0, s0) * FAinv;
		
                    integrate ( elements ( dispETFESpace->mesh() ) ,
                                quadRuleTetra4pt,
                                dispETFESpace,
                                dispETFESpace,
                                dot ( dP , grad (phi_i) )
                        ) >> jacobianPtr;
                }
                else
                {
                    auto gf = value (activationETFESpace, *gammaf);

                    if(orthotropicParameter > 0 )
                    {
                        if(disp.comm().MyPID() == 0)
                            std::cout << " Orthotropic case ... \n";

                        auto k = value(orthotropicParameter);
                        auto FAinv = _FAinv(gf, k, f0, s0, n0);

                        auto FE =  F * FAinv;
                        auto dFE = _dFE(FAinv);
                        auto dP = eval (d2W8fs, FE ) * (_dI8dF (FE, f0, s0, dFE) ) * _dI8(FE, f0, s0) * FAinv;

                        integrate ( elements ( dispETFESpace->mesh() ) ,
                                    quadRuleTetra4pt,
                                    dispETFESpace,
                                    dispETFESpace,
                                    dot ( dP , grad (phi_i) )
                            ) >> jacobianPtr;
                    }
                    else
                    {
                        if(disp.comm().MyPID() == 0)
                            std::cout << " Transversely isotropic case ... \n";

                        auto FAinv = _FAinv(gf, f0);


                        auto FE =  F * FAinv;
                        auto dFE = _dFE(FAinv);
                        auto dP = eval (d2W8fs, FE ) * (_dI8dF (FE, f0, s0, dFE) ) * _dI8(FE, f0, s0) * FAinv;

                        integrate ( elements ( dispETFESpace->mesh() ) ,
                                    quadRuleTetra4pt,
                                    dispETFESpace,
                                    dispETFESpace,
                                    dot ( dP , grad (phi_i) )
                            ) >> jacobianPtr;
                    }
                }

            }
        } // End ActiveStrain namespace

        namespace ActiveStrainNearlyIncompressible
        {
            template <typename Mesh, typename FunctorPtr >
            void
            computeActiveStrainI1barJacobianTerms (    const vector_Type& disp,
                                                    boost::shared_ptr<ETFESpace<Mesh, MapEpetra, 3, 3 > > dispETFESpace,
                                                    const vector_Type& fibers,
                                                    const vector_Type& sheets,
                                                    const vectorPtr_Type& gammaf,
                                                    const vectorPtr_Type& gammas,
                                                    const vectorPtr_Type& gamman,
                                                    boost::shared_ptr<ETFESpace<Mesh, MapEpetra, 3, 1 > >  activationETFESpace,
                                                    matrixPtr_Type           jacobianPtr,
                                                    FunctorPtr               W1,
                                                    Real orthotropicParameter = -666.)
            {
                //
                if(disp.comm().MyPID() == 0)
                    std::cout << "EMETA - Computing Isotropic Active Strain jacobian terms: ";

                using namespace ExpressionAssembly;

                auto I = _I;
                auto dF = grad(phi_j);
//                auto GradU = _Grad_u(dispETFESpace, disp, 0);
//                auto F = I + GradU;
                auto F = _F (dispETFESpace, disp, 0);

                auto f_0 = _v0 (dispETFESpace, fibers);
                auto s_0 = _v0 (dispETFESpace, sheets);

                boost::shared_ptr<orthonormalizeFibers> normalize0 (new orthonormalizeFibers);
                auto f0 = eval (normalize0, f_0);

                boost::shared_ptr<orthonormalizeFibers> normalize1 (new orthonormalizeFibers (1) );
                auto s0 = eval (normalize1, f0, s_0);

                boost::shared_ptr<CrossProduct> wedge (new CrossProduct);
                auto n0 = eval ( wedge, f0, s0);


                if(gammas && gamman)
                {
                    if(disp.comm().MyPID() == 0)
                        std::cout << " Anisotropic case ... \n";

                    auto gf = value (activationETFESpace, *gammaf);
                    auto gs = value (activationETFESpace, *gammas);
                    auto gn = value (activationETFESpace, *gamman);


                    auto FAinv = _FAinv(gf, gs, gn, f0, s0, n0);
                    auto FE =  F * FAinv;
                    auto dFE = _dFE(FAinv);

                    auto dP = eval (W1, FE ) * (_d2I1bardF (FE, dFE) ) * FAinv;

                    integrate ( elements ( dispETFESpace->mesh() ) ,
                                quadRuleTetra4pt,
                                dispETFESpace,
                                dispETFESpace,
                                dot ( dP , grad (phi_i) )
                        ) >> jacobianPtr;
                }
                else
                {
                    auto gf = value (activationETFESpace, *gammaf);

                    if(orthotropicParameter > 0 )
                    {
                        if(disp.comm().MyPID() == 0)
                            std::cout << " Orthotropic case ... \n";

                        auto k = value(orthotropicParameter);
                        auto FAinv = _FAinv(gf, k, f0, s0, n0);

                        auto FE =  F * FAinv;
                        auto dFE = _dFE(FAinv);
                        auto dP = eval (W1, FE ) * (_d2I1bardF (FE, dFE) ) * FAinv;

                        integrate ( elements ( dispETFESpace->mesh() ) ,
                                    quadRuleTetra4pt,
                                    dispETFESpace,
                                    dispETFESpace,
                                    dot ( dP , grad (phi_i) )
                            ) >> jacobianPtr;
                    }
                    else
                    {
                        if(disp.comm().MyPID() == 0)
                            std::cout << " Transversely isotropic case ... \n";

//    		auto FAinv = _FAinv(gf, f0, s0, n0);
                        auto FAinv = _FAinv(gf, f0);


                        auto FE =  F * FAinv;
                        auto dFE = _dFE(FAinv);
                        auto dP = eval (W1, FE ) * (_d2I1bardF (FE, dFE) ) * FAinv;

//    		auto dFE = dF * FAinv;
//    		auto dFET = FAinv * transpose(dF)
//
//    		auto mu = value(1000.);
//
//    		auto J = det(F);
//    		auto Jm23 = pow(J, -2./3.);
//    		auto FmT = minusT(F);
//    		auto FEmT = minusT(FE);
//    		auto dF = _dF;
//    		auto I1E = dot(FE, FE);
//
//    		auto varPE = mu * (FE + value(-1./3.) * I1E * FEmT);
//    		auto P = varPE * FAinv;
//
//    		auto dJm23dF = value(-2./3.) * Jm23 * dot( FmT, dF );
//
//    		auto dFEmTdFE = value(-1.) * FEmT * dFET * FEmT;
//    		auto dI1EdFE  = value(2.) * dot(FE, dFE);
//
//    		auto dP = dJm23dF * P + Jm23 * ( dFE + value(-1./3.) * ( dI1EdFE * FEmT + I1E * dFEmTdFE )  ) * FAinv;
//

                        integrate ( elements ( dispETFESpace->mesh() ) ,
                                    quadRuleTetra4pt,
                                    dispETFESpace,
                                    dispETFESpace,
                                    dot ( dP , grad (phi_i) )
                            ) >> jacobianPtr;
                    }
                }

            }
        
            template <typename Mesh, typename FunctorPtr >
            void
            computeActiveStrainI1barJacobianTermsSecondDerivative (    const vector_Type& disp,
                                                                    boost::shared_ptr<ETFESpace<Mesh, MapEpetra, 3, 3 > > dispETFESpace,
                                                                    const vector_Type& fibers,
                                                                    const vector_Type& sheets,
                                                                    const vectorPtr_Type& gammaf,
                                                                    const vectorPtr_Type& gammas,
                                                                    const vectorPtr_Type& gamman,
                                                                    boost::shared_ptr<ETFESpace<Mesh, MapEpetra, 3, 1 > >  activationETFESpace,
                                                                    matrixPtr_Type           jacobianPtr,
                                                                    FunctorPtr               d2W,
                                                                    Real orthotropicParameter = -666.)
            {
                //
                if(disp.comm().MyPID() == 0)
                    std::cout << "EMETA - Computing Isotropic Active Strain jacobian second derivative terms: ";

                using namespace ExpressionAssembly;

                auto I = _I;
                auto dF = grad(phi_j);
                auto GradU = _Grad_u(dispETFESpace, disp, 0);
                auto F = I + GradU;

                auto f_0 = _v0 (dispETFESpace, fibers);
                auto s_0 = _v0 (dispETFESpace, sheets);

                boost::shared_ptr<orthonormalizeFibers> normalize0 (new orthonormalizeFibers);
                auto f0 = eval (normalize0, f_0);

                boost::shared_ptr<orthonormalizeFibers> normalize1 (new orthonormalizeFibers (1) );
                auto s0 = eval (normalize1, f0, s_0);

                boost::shared_ptr<CrossProduct> wedge (new CrossProduct);
                auto n0 = eval ( wedge, f0, s0);


                if(gammas && gamman)
                {
                    if(disp.comm().MyPID() == 0)
                        std::cout << " Anisotropic case ... \n";

                    auto gf = value (activationETFESpace, *gammaf);
                    auto gs = value (activationETFESpace, *gammas);
                    auto gn = value (activationETFESpace, *gamman);


                    auto FAinv = _FAinv(gf, gs, gn, f0, s0, n0);
                    auto FE =  F * FAinv;
                    auto dFE = _dFE(FAinv);
                    auto dP = eval (d2W, FE ) * (_dI1bardF (FE, dFE) ) * _dI1bar(FE) * FAinv;
		
                    integrate ( elements ( dispETFESpace->mesh() ) ,
                                quadRuleTetra4pt,
                                dispETFESpace,
                                dispETFESpace,
                                dot ( dP , grad (phi_i) )
                        ) >> jacobianPtr;
                }
                else
                {
                    auto gf = value (activationETFESpace, *gammaf);

                    if(orthotropicParameter > 0 )
                    {
                        if(disp.comm().MyPID() == 0)
                            std::cout << " Orthotropic case ... \n";

                        auto k = value(orthotropicParameter);
                        auto FAinv = _FAinv(gf, k, f0, s0, n0);

                        auto FE =  F * FAinv;
                        auto dFE = _dFE(FAinv);
                        auto dP = eval (d2W, FE ) * (_dI1bardF (FE, dFE) ) * _dI1bar(FE) * FAinv;

                        integrate ( elements ( dispETFESpace->mesh() ) ,
                                    quadRuleTetra4pt,
                                    dispETFESpace,
                                    dispETFESpace,
                                    dot ( dP , grad (phi_i) )
                            ) >> jacobianPtr;
                    }
                    else
                    {
                        if(disp.comm().MyPID() == 0)
                            std::cout << " Transversely isotropic case ... \n";

//    		auto FAinv = _FAinv(gf, f0, s0, n0);
                        auto FAinv = _FAinv(gf, f0);


                        auto FE =  F * FAinv;
                        auto dFE = _dFE(FAinv);
                        auto dP = eval (d2W, FE ) * (_dI1bardF (FE, dFE) ) * _dI1bar(FE) * FAinv;

//    		auto dFE = dF * FAinv;
//    		auto dFET = FAinv * transpose(dF)
//
//    		auto mu = value(1000.);
//
//    		auto J = det(F);
//    		auto Jm23 = pow(J, -2./3.);
//    		auto FmT = minusT(F);
//    		auto FEmT = minusT(FE);
//    		auto dF = _dF;
//    		auto I1E = dot(FE, FE);
//
//    		auto varPE = mu * (FE + value(-1./3.) * I1E * FEmT);
//    		auto P = varPE * FAinv;
//
//    		auto dJm23dF = value(-2./3.) * Jm23 * dot( FmT, dF );
//
//    		auto dFEmTdFE = value(-1.) * FEmT * dFET * FEmT;
//    		auto dI1EdFE  = value(2.) * dot(FE, dFE);
//
//    		auto dP = dJm23dF * P + Jm23 * ( dFE + value(-1./3.) * ( dI1EdFE * FEmT + I1E * dFEmTdFE )  ) * FAinv;
//

                        integrate ( elements ( dispETFESpace->mesh() ) ,
                                    quadRuleTetra4pt,
                                    dispETFESpace,
                                    dispETFESpace,
                                    dot ( dP , grad (phi_i) )
                            ) >> jacobianPtr;
                    }
                }

            }
        
            template <typename Mesh, typename FunctorPtr >
            void
            computeActiveStrainI4barFibersJacobianTerms (    const vector_Type& disp,
                                                          boost::shared_ptr<ETFESpace<Mesh, MapEpetra, 3, 3 > > dispETFESpace,
                                                          const vector_Type& fibers,
                                                          const vector_Type& sheets,
                                                          const vectorPtr_Type& gammaf,
                                                          const vectorPtr_Type& gammas,
                                                          const vectorPtr_Type& gamman,
                                                          boost::shared_ptr<ETFESpace<Mesh, MapEpetra, 3, 1 > >  activationETFESpace,
                                                          matrixPtr_Type           jacobianPtr,
                                                          FunctorPtr               W4f,
                                                          Real orthotropicParameter = -666.)
            {
                //
                if(disp.comm().MyPID() == 0)
                    std::cout << "EMETA - Computing Anisotropic Fibers Active Strain jacobian terms: ";

                using namespace ExpressionAssembly;

                auto I = _I;
                auto dF = grad(phi_j);
                auto GradU = _Grad_u(dispETFESpace, disp, 0);
                auto F = I + GradU;

                auto f_0 = _v0 (dispETFESpace, fibers);
                auto s_0 = _v0 (dispETFESpace, sheets);

                boost::shared_ptr<orthonormalizeFibers> normalize0 (new orthonormalizeFibers);
                auto f0 = eval (normalize0, f_0);

                boost::shared_ptr<orthonormalizeFibers> normalize1 (new orthonormalizeFibers (1) );
                auto s0 = eval (normalize1, f0, s_0);

                boost::shared_ptr<CrossProduct> wedge (new CrossProduct);
                auto n0 = eval ( wedge, f0, s0);


                if(gammas && gamman)
                {
                    if(disp.comm().MyPID() == 0)
                        std::cout << " Anisotropic case ... \n";

                    auto gf = value (activationETFESpace, *gammaf);
                    auto gs = value (activationETFESpace, *gammas);
                    auto gn = value (activationETFESpace, *gamman);


                    auto FAinv = _FAinv(gf, gs, gn, f0, s0, n0);
                    auto FE =  F * FAinv;
                    auto dFE = _dFE(FAinv);

                    auto dP = eval (W4f, FE ) * (_d2I4bardF (FE, f0, dFE) ) * FAinv;

                    integrate ( elements ( dispETFESpace->mesh() ) ,
                                quadRuleTetra4pt,
                                dispETFESpace,
                                dispETFESpace,
                                dot ( dP , grad (phi_i) )
                        ) >> jacobianPtr;
                }
                else
                {
                    auto gf = value (activationETFESpace, *gammaf);

                    if(orthotropicParameter > 0 )
                    {
                        if(disp.comm().MyPID() == 0)
                            std::cout << " Orthotropic case ... \n";

                        auto k = value(orthotropicParameter);
                        auto FAinv = _FAinv(gf, k, f0, s0, n0);

                        auto FE =  F * FAinv;
                        auto dFE = _dFE(FAinv);
                        auto dP = eval (W4f, FE ) * (_d2I4bardF (FE, f0, dFE) ) * FAinv;

                        integrate ( elements ( dispETFESpace->mesh() ) ,
                                    quadRuleTetra4pt,
                                    dispETFESpace,
                                    dispETFESpace,
                                    dot ( dP , grad (phi_i) )
                            ) >> jacobianPtr;
                    }
                    else
                    {
                        if(disp.comm().MyPID() == 0)
                            std::cout << " Transversely isotropic case ... \n";

                        auto FAinv = _FAinv(gf, f0);

                        auto FE =  F * FAinv;
                        auto dFE = _dFE(FAinv);
                        auto dP = eval (W4f, FE ) * (_d2I4bardF (FE, f0, dFE) ) * FAinv;
                    
                        integrate ( elements ( dispETFESpace->mesh() ) ,
                                    quadRuleTetra4pt,
                                    dispETFESpace,
                                    dispETFESpace,
                                    dot ( dP , grad (phi_i) )
                            ) >> jacobianPtr;
                    }
                }

            }
        
            template <typename Mesh, typename FunctorPtr >
            void
            computeActiveStrainI4barFibersJacobianTermsSecondDerivative (    const vector_Type& disp,
                                                                          boost::shared_ptr<ETFESpace<Mesh, MapEpetra, 3, 3 > > dispETFESpace,
                                                                          const vector_Type& fibers,
                                                                          const vector_Type& sheets,
                                                                          const vectorPtr_Type& gammaf,
                                                                          const vectorPtr_Type& gammas,
                                                                          const vectorPtr_Type& gamman,
                                                                          boost::shared_ptr<ETFESpace<Mesh, MapEpetra, 3, 1 > >  activationETFESpace,
                                                                          matrixPtr_Type           jacobianPtr,
                                                                          FunctorPtr               d2W4f,
                                                                          Real orthotropicParameter = -666.)
            {
                //
                if(disp.comm().MyPID() == 0)
                    std::cout << "EMETA - Computing Anisotropic Fibers Active Strain jacobian second derivative terms: ";

                using namespace ExpressionAssembly;

                auto I = _I;
                auto dF = grad(phi_j);
                auto GradU = _Grad_u(dispETFESpace, disp, 0);
                auto F = I + GradU;

                auto f_0 = _v0 (dispETFESpace, fibers);
                auto s_0 = _v0 (dispETFESpace, sheets);

                boost::shared_ptr<orthonormalizeFibers> normalize0 (new orthonormalizeFibers);
                auto f0 = eval (normalize0, f_0);

                boost::shared_ptr<orthonormalizeFibers> normalize1 (new orthonormalizeFibers (1) );
                auto s0 = eval (normalize1, f0, s_0);

                boost::shared_ptr<CrossProduct> wedge (new CrossProduct);
                auto n0 = eval ( wedge, f0, s0);

                if(gammas && gamman)
                {
                    if(disp.comm().MyPID() == 0)
                        std::cout << " Anisotropic case ... \n";

                    auto gf = value (activationETFESpace, *gammaf);
                    auto gs = value (activationETFESpace, *gammas);
                    auto gn = value (activationETFESpace, *gamman);

                    auto FAinv = _FAinv(gf, gs, gn, f0, s0, n0);
                    auto FE =  F * FAinv;
                    auto dFE = _dFE(FAinv);
                    auto dP = eval (d2W4f, FE ) * (_dI4bardF (FE, f0, dFE) ) * _dI4bar(FE, f0) * FAinv;
		
                    integrate ( elements ( dispETFESpace->mesh() ) ,
                                quadRuleTetra4pt,
                                dispETFESpace,
                                dispETFESpace,
                                dot ( dP , grad (phi_i) )
                        ) >> jacobianPtr;
                }
                else
                {
                    auto gf = value (activationETFESpace, *gammaf);

                    if(orthotropicParameter > 0 )
                    {
                        if(disp.comm().MyPID() == 0)
                            std::cout << " Orthotropic case ... \n";

                        auto k = value(orthotropicParameter);
                        auto FAinv = _FAinv(gf, k, f0, s0, n0);

                        auto FE =  F * FAinv;
                        auto dFE = _dFE(FAinv);
                        auto dP = eval (d2W4f, FE ) * (_dI4bardF (FE, f0, dFE) ) * _dI4bar(FE, f0) * FAinv;

                        integrate ( elements ( dispETFESpace->mesh() ) ,
                                    quadRuleTetra4pt,
                                    dispETFESpace,
                                    dispETFESpace,
                                    dot ( dP , grad (phi_i) )
                            ) >> jacobianPtr;
                    }
                    else
                    {
                        if(disp.comm().MyPID() == 0)
                            std::cout << " Transversely isotropic case ... \n";

                        auto FAinv = _FAinv(gf, f0);


                        auto FE =  F * FAinv;
                        auto dFE = _dFE(FAinv);
                        auto dP = eval (d2W4f, FE ) * (_dI4bardF (FE, f0, dFE) ) * _dI4bar(FE, f0) * FAinv;

                        integrate ( elements ( dispETFESpace->mesh() ) ,
                                    quadRuleTetra4pt,
                                    dispETFESpace,
                                    dispETFESpace,
                                    dot ( dP , grad (phi_i) )
                            ) >> jacobianPtr;
                    }
                }

            }

            template <typename Mesh, typename FunctorPtr >
            void
            computeActiveStrainI4barSheetsJacobianTerms (    const vector_Type& disp,
                                                          boost::shared_ptr<ETFESpace<Mesh, MapEpetra, 3, 3 > > dispETFESpace,
                                                          const vector_Type& fibers,
                                                          const vector_Type& sheets,
                                                          const vectorPtr_Type& gammaf,
                                                          const vectorPtr_Type& gammas,
                                                          const vectorPtr_Type& gamman,
                                                          boost::shared_ptr<ETFESpace<Mesh, MapEpetra, 3, 1 > >  activationETFESpace,
                                                          matrixPtr_Type           jacobianPtr,
                                                          FunctorPtr               W4s,
                                                          Real orthotropicParameter = -666.)
            {
                //
                if(disp.comm().MyPID() == 0)
                    std::cout << "EMETA - Computing Anisotropic Sheets Active Strain jacobian terms: ";
            
                using namespace ExpressionAssembly;

                auto I = _I;
                auto dF = grad(phi_j);
                auto GradU = _Grad_u(dispETFESpace, disp, 0);
                auto F = I + GradU;

                auto f_0 = _v0 (dispETFESpace, fibers);
                auto s_0 = _v0 (dispETFESpace, sheets);

                boost::shared_ptr<orthonormalizeFibers> normalize0 (new orthonormalizeFibers);
                auto f0 = eval (normalize0, f_0);

                boost::shared_ptr<orthonormalizeFibers> normalize1 (new orthonormalizeFibers (1) );
                auto s0 = eval (normalize1, f0, s_0);

                boost::shared_ptr<CrossProduct> wedge (new CrossProduct);
                auto n0 = eval ( wedge, f0, s0);


                if(gammas && gamman)
                {
                    if(disp.comm().MyPID() == 0)
                        std::cout << " Anisotropic case ... \n";

                    auto gf = value (activationETFESpace, *gammaf);
                    auto gs = value (activationETFESpace, *gammas);
                    auto gn = value (activationETFESpace, *gamman);


                    auto FAinv = _FAinv(gf, gs, gn, f0, s0, n0);
                    auto FE =  F * FAinv;
                    auto dFE = _dFE(FAinv);
                    auto dP = eval (W4s, FE ) * (_d2I4bardF (FE, s0, dFE) ) * FAinv;

                    integrate ( elements ( dispETFESpace->mesh() ) ,
                                quadRuleTetra4pt,
                                dispETFESpace,
                                dispETFESpace,
                                dot ( dP , grad (phi_i) )
                        ) >> jacobianPtr;
                }
                else
                {
                    auto gf = value (activationETFESpace, *gammaf);

                    if(orthotropicParameter > 0 )
                    {
                        if(disp.comm().MyPID() == 0)
                            std::cout << " Orthotropic case ... \n";

                        auto k = value(orthotropicParameter);
                        auto FAinv = _FAinv(gf, k, f0, s0, n0);

                        auto FE =  F * FAinv;
                        auto dFE = _dFE(FAinv);
                        auto dP = eval (W4s, FE ) * (_d2I4bardF (FE, s0, dFE) ) * FAinv;

                        integrate ( elements ( dispETFESpace->mesh() ) ,
                                    quadRuleTetra4pt,
                                    dispETFESpace,
                                    dispETFESpace,
                                    dot ( dP , grad (phi_i) )
                            ) >> jacobianPtr;
                    }
                    else
                    {
                        if(disp.comm().MyPID() == 0)
                            std::cout << " Transversely isotropic case ... \n";

                        auto FAinv = _FAinv(gf, f0);

                        auto FE =  F * FAinv;
                        auto dFE = _dFE(FAinv);
                        auto dP = eval (W4s, FE ) * (_d2I4bardF (FE, s0, dFE) ) * FAinv;
                    
                        integrate ( elements ( dispETFESpace->mesh() ) ,
                                    quadRuleTetra4pt,
                                    dispETFESpace,
                                    dispETFESpace,
                                    dot ( dP , grad (phi_i) )
                            ) >> jacobianPtr;
                    }
                }

            }
        
            template <typename Mesh, typename FunctorPtr >
            void
            computeActiveStrainI4barSheetsJacobianTermsSecondDerivative (    const vector_Type& disp,
                                                                          boost::shared_ptr<ETFESpace<Mesh, MapEpetra, 3, 3 > > dispETFESpace,
                                                                          const vector_Type& fibers,
                                                                          const vector_Type& sheets,
                                                                          const vectorPtr_Type& gammaf,
                                                                          const vectorPtr_Type& gammas,
                                                                          const vectorPtr_Type& gamman,
                                                                          boost::shared_ptr<ETFESpace<Mesh, MapEpetra, 3, 1 > >  activationETFESpace,
                                                                          matrixPtr_Type           jacobianPtr,
                                                                          FunctorPtr               d2W4s,
                                                                          Real orthotropicParameter = -666.)
            {
                //
                if(disp.comm().MyPID() == 0)
                    std::cout << "EMETA - Computing Anisotropic Sheets Active Strain jacobian second derivative terms: ";

                using namespace ExpressionAssembly;

                auto I = _I;
                auto dF = grad(phi_j);
                auto GradU = _Grad_u(dispETFESpace, disp, 0);
                auto F = I + GradU;

                auto f_0 = _v0 (dispETFESpace, fibers);
                auto s_0 = _v0 (dispETFESpace, sheets);

                boost::shared_ptr<orthonormalizeFibers> normalize0 (new orthonormalizeFibers);
                auto f0 = eval (normalize0, f_0);

                boost::shared_ptr<orthonormalizeFibers> normalize1 (new orthonormalizeFibers (1) );
                auto s0 = eval (normalize1, f0, s_0);

                boost::shared_ptr<CrossProduct> wedge (new CrossProduct);
                auto n0 = eval ( wedge, f0, s0);

                if(gammas && gamman)
                {
                    if(disp.comm().MyPID() == 0)
                        std::cout << " Anisotropic case ... \n";

                    auto gf = value (activationETFESpace, *gammaf);
                    auto gs = value (activationETFESpace, *gammas);
                    auto gn = value (activationETFESpace, *gamman);

                    auto FAinv = _FAinv(gf, gs, gn, f0, s0, n0);
                    auto FE =  F * FAinv;
                    auto dFE = _dFE(FAinv);
                    auto dP = eval (d2W4s, FE ) * (_dI4bardF (FE, s0, dFE) ) * _dI4bar(FE, s0) * FAinv;
		
                    integrate ( elements ( dispETFESpace->mesh() ) ,
                                quadRuleTetra4pt,
                                dispETFESpace,
                                dispETFESpace,
                                dot ( dP , grad (phi_i) )
                        ) >> jacobianPtr;
                }
                else
                {
                    auto gf = value (activationETFESpace, *gammaf);

                    if(orthotropicParameter > 0 )
                    {
                        if(disp.comm().MyPID() == 0)
                            std::cout << " Orthotropic case ... \n";

                        auto k = value(orthotropicParameter);
                        auto FAinv = _FAinv(gf, k, f0, s0, n0);

                        auto FE =  F * FAinv;
                        auto dFE = _dFE(FAinv);
                        auto dP = eval (d2W4s, FE ) * (_dI4bardF (FE, s0, dFE) ) * _dI4bar(FE, s0) * FAinv;

                        integrate ( elements ( dispETFESpace->mesh() ) ,
                                    quadRuleTetra4pt,
                                    dispETFESpace,
                                    dispETFESpace,
                                    dot ( dP , grad (phi_i) )
                            ) >> jacobianPtr;
                    }
                    else
                    {
                        if(disp.comm().MyPID() == 0)
                            std::cout << " Transversely isotropic case ... \n";

                        auto FAinv = _FAinv(gf, f0);


                        auto FE =  F * FAinv;
                        auto dFE = _dFE(FAinv);
                        auto dP = eval (d2W4s, FE ) * (_dI4bardF (FE, s0, dFE) ) * _dI4bar(FE, s0) * FAinv;

                        integrate ( elements ( dispETFESpace->mesh() ) ,
                                    quadRuleTetra4pt,
                                    dispETFESpace,
                                    dispETFESpace,
                                    dot ( dP , grad (phi_i) )
                            ) >> jacobianPtr;
                    }
                }

            }
        
            template <typename Mesh, typename FunctorPtr >
            void
            computeActiveStrainI8barJacobianTerms (    const vector_Type& disp,
                                                    boost::shared_ptr<ETFESpace<Mesh, MapEpetra, 3, 3 > > dispETFESpace,
                                                    const vector_Type& fibers,
                                                    const vector_Type& sheets,
                                                    const vectorPtr_Type& gammaf,
                                                    const vectorPtr_Type& gammas,
                                                    const vectorPtr_Type& gamman,
                                                    boost::shared_ptr<ETFESpace<Mesh, MapEpetra, 3, 1 > >  activationETFESpace,
                                                    matrixPtr_Type           jacobianPtr,
                                                    FunctorPtr               W8fs,
                                                    Real orthotropicParameter = -666.)
            {
                //
                if(disp.comm().MyPID() == 0)
                    std::cout << "EMETA - Computing Anisotropic Shear Active Strain jacobian terms: ";

                using namespace ExpressionAssembly;

                auto I = _I;
                auto dF = grad(phi_j);
                auto GradU = _Grad_u(dispETFESpace, disp, 0);
                auto F = I + GradU;

                auto f_0 = _v0 (dispETFESpace, fibers);
                auto s_0 = _v0 (dispETFESpace, sheets);

                boost::shared_ptr<orthonormalizeFibers> normalize0 (new orthonormalizeFibers);
                auto f0 = eval (normalize0, f_0);

                boost::shared_ptr<orthonormalizeFibers> normalize1 (new orthonormalizeFibers (1) );
                auto s0 = eval (normalize1, f0, s_0);

                boost::shared_ptr<CrossProduct> wedge (new CrossProduct);
                auto n0 = eval ( wedge, f0, s0);


                if(gammas && gamman)
                {
                    if(disp.comm().MyPID() == 0)
                        std::cout << " Anisotropic case ... \n";

                    auto gf = value (activationETFESpace, *gammaf);
                    auto gs = value (activationETFESpace, *gammas);
                    auto gn = value (activationETFESpace, *gamman);


                    auto FAinv = _FAinv(gf, gs, gn, f0, s0, n0);
                    auto FE =  F * FAinv;
                    auto dFE = _dFE(FAinv);

                    auto dP = eval (W8fs, FE ) * (_d2I8bardF (FE, f0, s0, dFE) ) * FAinv;

                    integrate ( elements ( dispETFESpace->mesh() ) ,
                                quadRuleTetra4pt,
                                dispETFESpace,
                                dispETFESpace,
                                dot ( dP , grad (phi_i) )
                        ) >> jacobianPtr;
                }
                else
                {
                    auto gf = value (activationETFESpace, *gammaf);

                    if(orthotropicParameter > 0 )
                    {
                        if(disp.comm().MyPID() == 0)
                            std::cout << " Orthotropic case ... \n";

                        auto k = value(orthotropicParameter);
                        auto FAinv = _FAinv(gf, k, f0, s0, n0);

                        auto FE =  F * FAinv;
                        auto dFE = _dFE(FAinv);
                        auto dP = eval (W8fs, FE ) * (_d2I8bardF (FE, f0, s0, dFE) ) * FAinv;

                        integrate ( elements ( dispETFESpace->mesh() ) ,
                                    quadRuleTetra4pt,
                                    dispETFESpace,
                                    dispETFESpace,
                                    dot ( dP , grad (phi_i) )
                            ) >> jacobianPtr;
                    }
                    else
                    {
                        if(disp.comm().MyPID() == 0)
                            std::cout << " Transversely isotropic case ... \n";

                        auto FAinv = _FAinv(gf, f0);


                        auto FE =  F * FAinv;
                        auto dFE = _dFE(FAinv);
                        auto dP = eval (W8fs, FE ) * (_d2I8bardF (FE, f0, s0, dFE) ) * FAinv;

                        integrate ( elements ( dispETFESpace->mesh() ) ,
                                    quadRuleTetra4pt,
                                    dispETFESpace,
                                    dispETFESpace,
                                    dot ( dP , grad (phi_i) )
                            ) >> jacobianPtr;
                    }
                }

            }
        
            template <typename Mesh, typename FunctorPtr >
            void
            computeActiveStrainI8barJacobianTermsSecondDerivative (    const vector_Type& disp,
                                                                    boost::shared_ptr<ETFESpace<Mesh, MapEpetra, 3, 3 > > dispETFESpace,
                                                                    const vector_Type& fibers,
                                                                    const vector_Type& sheets,
                                                                    const vectorPtr_Type& gammaf,
                                                                    const vectorPtr_Type& gammas,
                                                                    const vectorPtr_Type& gamman,
                                                                    boost::shared_ptr<ETFESpace<Mesh, MapEpetra, 3, 1 > >  activationETFESpace,
                                                                    matrixPtr_Type           jacobianPtr,
                                                                    FunctorPtr               d2W8fs,
                                                                    Real orthotropicParameter = -666.)
            {
                //
                if(disp.comm().MyPID() == 0)
                    std::cout << "EMETA - Computing Anisotropic Shear Active Strain jacobian second derivative terms: ";

                using namespace ExpressionAssembly;

                auto I = _I;
                auto dF = grad(phi_j);
                auto GradU = _Grad_u(dispETFESpace, disp, 0);
                auto F = I + GradU;

                auto f_0 = _v0 (dispETFESpace, fibers);
                auto s_0 = _v0 (dispETFESpace, sheets);

                boost::shared_ptr<orthonormalizeFibers> normalize0 (new orthonormalizeFibers);
                auto f0 = eval (normalize0, f_0);

                boost::shared_ptr<orthonormalizeFibers> normalize1 (new orthonormalizeFibers (1) );
                auto s0 = eval (normalize1, f0, s_0);

                boost::shared_ptr<CrossProduct> wedge (new CrossProduct);
                auto n0 = eval ( wedge, f0, s0);


                if(gammas && gamman)
                {
                    if(disp.comm().MyPID() == 0)
                        std::cout << " Anisotropic case ... \n";

                    auto gf = value (activationETFESpace, *gammaf);
                    auto gs = value (activationETFESpace, *gammas);
                    auto gn = value (activationETFESpace, *gamman);

                    auto FAinv = _FAinv(gf, gs, gn, f0, s0, n0);
                    auto FE =  F * FAinv;
                    auto dFE = _dFE(FAinv);
                    auto dP = eval (d2W8fs, FE ) * (_dI8bardF (FE, f0, s0, dFE) ) * _dI8bar(FE, f0, s0) * FAinv;
		
                    integrate ( elements ( dispETFESpace->mesh() ) ,
                                quadRuleTetra4pt,
                                dispETFESpace,
                                dispETFESpace,
                                dot ( dP , grad (phi_i) )
                        ) >> jacobianPtr;
                }
                else
                {
                    auto gf = value (activationETFESpace, *gammaf);

                    if(orthotropicParameter > 0 )
                    {
                        if(disp.comm().MyPID() == 0)
                            std::cout << " Orthotropic case ... \n";

                        auto k = value(orthotropicParameter);
                        auto FAinv = _FAinv(gf, k, f0, s0, n0);

                        auto FE =  F * FAinv;
                        auto dFE = _dFE(FAinv);
                        auto dP = eval (d2W8fs, FE ) * (_dI8bardF (FE, f0, s0, dFE) ) * _dI8bar(FE, f0, s0) * FAinv;

                        integrate ( elements ( dispETFESpace->mesh() ) ,
                                    quadRuleTetra4pt,
                                    dispETFESpace,
                                    dispETFESpace,
                                    dot ( dP , grad (phi_i) )
                            ) >> jacobianPtr;
                    }
                    else
                    {
                        if(disp.comm().MyPID() == 0)
                            std::cout << " Transversely isotropic case ... \n";

                        auto FAinv = _FAinv(gf, f0);


                        auto FE =  F * FAinv;
                        auto dFE = _dFE(FAinv);
                        auto dP = eval (d2W8fs, FE ) * (_dI8bardF (FE, f0, s0, dFE) ) * _dI8bar(FE, f0, s0) * FAinv;

                        integrate ( elements ( dispETFESpace->mesh() ) ,
                                    quadRuleTetra4pt,
                                    dispETFESpace,
                                    dispETFESpace,
                                    dot ( dP , grad (phi_i) )
                            ) >> jacobianPtr;
                    }
                }

            }
        }
        
    }//EMAssembler

}//LifeV

#endif /* EMETAACTIVESTRAINJACOBIANASSMEBLER_HPP_ */
