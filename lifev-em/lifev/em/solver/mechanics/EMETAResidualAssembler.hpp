/*
 * EMETAAssmebler.hpp
 *
 *  Created on: 29/apr/2014
 *      Author: srossi
 */

#ifndef EMETARESIDUALASSMEBLER_HPP_
#define EMETARESIDUALASSMEBLER_HPP_

#include <lifev/core/fem/FESpace.hpp>
#include <lifev/core/array/VectorEpetra.hpp>
#include <lifev/core/array/MatrixEpetra.hpp>

//ET include for assemblings
#include <lifev/eta/fem/ETFESpace.hpp>
#include <lifev/eta/expression/Integrate.hpp>

#include <lifev/em/solver/mechanics/EMMechanicalExpressions.hpp>

#include <lifev/em/solver/EMETAFunctors.hpp>

namespace LifeV
{

typedef VectorEpetra           vector_Type;
typedef boost::shared_ptr<vector_Type>         vectorPtr_Type;

typedef MatrixEpetra<Real>           matrix_Type;
typedef boost::shared_ptr<matrix_Type>         matrixPtr_Type;

namespace EMAssembler
{


template <typename Mesh, typename FunctorPtr >
void
computeOrthotropicActiveStressResidualTerms ( const vector_Type& disp,
                                       boost::shared_ptr<ETFESpace<Mesh, MapEpetra, 3, 3 > >  dispETFESpace,
                                       const vector_Type& fibers,
                                       const vector_Type& sheets,
                                       const vector_Type& activation,
                                       boost::shared_ptr<ETFESpace<Mesh, MapEpetra, 3, 1 > >  activationETFESpace,
                                       vectorPtr_Type           residualVectorPtr,
                                       FunctorPtr               W,
                                       const std::string& field)
{
    using namespace ExpressionAssembly;

    auto I = _I;
    
    auto f_0 = _v0 (dispETFESpace, fibers);
    boost::shared_ptr<orthonormalizeFibers> normalize0 (new orthonormalizeFibers);
    auto f0 = eval (normalize0, f_0);

    auto F = _F (dispETFESpace, disp, 0);
    auto H = value (activationETFESpace, activation);
    auto Wa = eval (W, H);
    auto Wm = eval (W, I);
    
    if ( field == "Sheets" )
    {
        auto s_0 = _v0 (dispETFESpace, sheets);
        auto s_00 = s_0 - dot (f0, s_0) * f0;
        boost::shared_ptr<orthonormalizeFibers> normalize1 (new orthonormalizeFibers (1) );
        auto s0 = eval (normalize1, s_00);
        
        auto i = F * s0;
        auto ixi0 = outerProduct (i, f0);
        
        auto P = Wa * /*Wm*/ ixi0;
        integrate ( elements ( dispETFESpace->mesh() ) ,
                   quadRule(),
                   dispETFESpace,
                   dot ( P , grad (phi_i) )
                   ) >> residualVectorPtr;
    }
    else if ( field == "Normal" )
    {
        auto s_0 = _v0 (dispETFESpace, sheets);
        auto s_00 = s_0 - dot (f0, s_0) * f0;
        boost::shared_ptr<orthonormalizeFibers> normalize1 (new orthonormalizeFibers (1) );
        auto s0 = eval (normalize1, s_00);
        
        boost::shared_ptr<CrossProduct> wedge (new CrossProduct);
        auto n0 = eval ( wedge, f0, s0);
        
        auto i = F * n0;
        auto ixi0 = outerProduct (i, f0);
        
        auto P = Wa * /*Wm*/ ixi0;
        integrate ( elements ( dispETFESpace->mesh() ) ,
                   quadRule(),
                   dispETFESpace,
                   dot ( P , grad (phi_i) )
                   ) >> residualVectorPtr;
    }
    else
    {
        auto i = F * f0;
        auto ixi0 = outerProduct (i, f0);
        
        auto P = Wa * /*Wm*/ ixi0;
        integrate ( elements ( dispETFESpace->mesh() ) ,
                   quadRule(),
                   dispETFESpace,
                   dot ( P , grad (phi_i) )
                   ) >> residualVectorPtr;
    }
}

    
template <typename Mesh, typename FunctorPtr >
void
computeFiberActiveStressResidualTerms ( const vector_Type& disp,
                                        boost::shared_ptr<ETFESpace<Mesh, MapEpetra, 3, 3 > >  dispETFESpace,
                                        const vector_Type& fibers,
                                        const vector_Type& activation,
                                        boost::shared_ptr<ETFESpace<Mesh, MapEpetra, 3, 1 > >  activationETFESpace,
                                        vectorPtr_Type           residualVectorPtr,
                                        FunctorPtr               W)
{
    //
	//if(disp.comm().MyPID() == 0)
    //std::cout << "EMETA - Computing Fibers Active Stress residual terms ... \n";

    using namespace ExpressionAssembly;

    auto F = _F (dispETFESpace, disp, 0);
    auto I = _I;

    auto f0 = value (dispETFESpace, fibers);
    auto f = F * f0;
    auto fxf0 = outerProduct (f, f0);
    auto H = value (activationETFESpace, activation);
    auto Wa = eval (W, H);
    auto Wm = eval (W, I);
    auto P = Wa /* Wm */ *  fxf0;

    integrate ( elements ( dispETFESpace->mesh() ) ,
                quadRule(),
                dispETFESpace,
                dot ( Wa * outerProduct ( F * f0, f0), grad (phi_i) )
              ) >> residualVectorPtr;

}


template <typename Mesh, typename FunctorPtr >
void
computeModifiedFiberActiveStressResidualTerms ( const vector_Type& disp,
												boost::shared_ptr<ETFESpace<Mesh, MapEpetra, 3, 3 > >  dispETFESpace,
												const vector_Type& fibers,
												const vector_Type& activation,
												boost::shared_ptr<ETFESpace<Mesh, MapEpetra, 3, 1 > >  activationETFESpace,
												vectorPtr_Type           residualVectorPtr,
												FunctorPtr               W)
{
    //
	if(disp.comm().MyPID() == 0)
    std::cout << "EMETA - Computing Fibers Active Stress residual terms ... \n";

    using namespace ExpressionAssembly;

    auto F = _F (dispETFESpace, disp, 0);
    auto I = _I;

    auto f0 = value (dispETFESpace, fibers);
    auto f = F * f0;
    auto fxf0 = outerProduct (f, f0);
    auto H = value (activationETFESpace, activation);
    auto Wa = eval (W, H);
    auto Wm = eval (W, I);
//    auto P = Wa /* Wm */ *  fxf0;

    auto P = Wa /* Wm */ * _dI4bar (F, f0);

    integrate ( elements ( dispETFESpace->mesh() ) ,
                quadRule(),
                dispETFESpace,
                dot ( P, grad (phi_i) )
              ) >> residualVectorPtr;

}


template< typename Mesh, typename FunctorPtr >
void
computeI4ResidualTerms ( const vector_Type& disp,
                         boost::shared_ptr<ETFESpace<Mesh, MapEpetra, 3, 3 > > dispETFESpace,
                         const vector_Type& fibers,
                         vectorPtr_Type     residualVectorPtr,
                         FunctorPtr         W4)
{
    using namespace ExpressionAssembly;
    //
	//if(disp.comm().MyPID() == 0)
    //std::cout << "EMETA - Computing I4 f residual terms ... \n";

    auto f_0 = _v0 (dispETFESpace, fibers);
    boost::shared_ptr<orthonormalizeFibers> normalize0 (new orthonormalizeFibers);
    boost::shared_ptr<ShowValue> sv ( new ShowValue() );

    auto f0 = eval (normalize0, f_0);
	auto F = _F (dispETFESpace, disp, 0);

    auto P = eval (W4, _I4 ( F, f0 ) )
             * _dI4 (F, f0);
//
//    auto P = eval (W4, _I4bar ( dispETFESpace, disp, 0, f0 ) )
//             * _dI4bar (dispETFESpace, disp, 0, f0);

    integrate ( elements ( dispETFESpace->mesh() ) ,
                quadRule(),
                dispETFESpace,
                dot ( P , grad (phi_i) )
              ) >> residualVectorPtr;
}


template< typename Mesh, typename FunctorPtr >
void
computeI4ResidualTerms ( const vector_Type& disp,
                         boost::shared_ptr<ETFESpace<Mesh, MapEpetra, 3, 3 > > dispETFESpace,
                         const vector_Type& fibers,
                         const vector_Type& sheets,
                         vectorPtr_Type     residualVectorPtr,
                         FunctorPtr         W4)
{
    using namespace ExpressionAssembly;
    //
	//if(disp.comm().MyPID() == 0)
    //std::cout << "EMETA - Computing I4 s residual terms ... \n";

    auto f_0 = _v0 (dispETFESpace, fibers);
    auto s_0 = _v0 (dispETFESpace, sheets);

    boost::shared_ptr<orthonormalizeFibers> normalize0 (new orthonormalizeFibers);
    boost::shared_ptr<orthonormalizeFibers> normalize1 (new orthonormalizeFibers (1) );
    auto f0 = eval (normalize0, f_0);

    auto s_00 = s_0 - dot (f0, s_0) * f0;

    auto s0 = eval (normalize1, s_00);

	auto F = _F (dispETFESpace, disp, 0);

    auto P = eval (W4, _I4 ( F, s0 ) )
             * _dI4 ( F, s0);


    integrate ( elements ( dispETFESpace->mesh() ) ,
                quadRule(),
                dispETFESpace,
                dot ( P , grad (phi_i) )
              ) >> residualVectorPtr;
}


template< typename Mesh, typename FunctorPtr >
void
computeI4ResidualTermsFung ( const vector_Type& disp,
                             boost::shared_ptr<ETFESpace<Mesh, MapEpetra, 3, 3 > >  dispETFESpace,
                             const vector_Type& fibers,
                             const vector_Type& sheets,
                             vectorPtr_Type     residualVectorPtr,
                             FunctorPtr         W4,
                             Int Case = 0)
{
    using namespace ExpressionAssembly;
    //
	//if(disp.comm().MyPID() == 0)
    //std::cout << "EMETA - Computing I4 s residual terms (Fung) ... \n";

    auto f_0 = _v0 (dispETFESpace, fibers);
    auto s_0 = _v0 (dispETFESpace, sheets);

    boost::shared_ptr<orthonormalizeFibers> normalize0 (new orthonormalizeFibers);
    boost::shared_ptr<orthonormalizeFibers> normalize1 (new orthonormalizeFibers (1) );
    auto f0 = eval (normalize0, f_0);

    auto s_00 = s_0 - dot (f0, s_0) * f0;

    auto s0 = eval (normalize1, s_00);
	auto F = _F (dispETFESpace, disp, 0);

    if (Case == 0)
    {
        auto P = eval (W4, F, f0, s0 )
                 * _dI4 (F, f0);


        integrate ( elements ( dispETFESpace->mesh() ) ,
                    quadRule(),
                    dispETFESpace,
                    dot ( P , grad (phi_i) )
                  ) >> residualVectorPtr;
    }
    else
    {
        auto P = eval (W4, F, f0, s0 )
                 * _dI4 ( F, s0);


        integrate ( elements ( dispETFESpace->mesh() ) ,
                    quadRule(),
                    dispETFESpace,
                    dot ( P , grad (phi_i) )
                  ) >> residualVectorPtr;
    }

}

template< typename Mesh, typename FunctorPtr >
void
computeI8ResidualTerms ( const vector_Type& disp,
                         boost::shared_ptr<ETFESpace<Mesh, MapEpetra, 3, 3 > >  dispETFESpace,
                         const vector_Type& fibers,
                         const vector_Type& sheets,
                         vectorPtr_Type           residualVectorPtr,
                         FunctorPtr                  W8,
                         bool orthonormalize = true)
{
    using namespace ExpressionAssembly;

	auto F = _F (dispETFESpace, disp, 0);

    if (orthonormalize)
    {
        auto f_0 = _v0 (dispETFESpace, fibers);
        auto s_0 = _v0 (dispETFESpace, sheets);

        boost::shared_ptr<orthonormalizeFibers> normalize0 (new orthonormalizeFibers);
        auto f0 = eval (normalize0, f_0);

        boost::shared_ptr<orthonormalizeFibers> normalize1 (new orthonormalizeFibers (1) );
        auto s0 = eval (normalize1, f0, s_0);

        auto P = eval (W8, _I8 ( F, f0, s0 ) )
                 * _dI8 ( F, f0, s0);

    	//if(disp.comm().MyPID() == 0)
        //std::cout << "EMETA - Computing I8 residual terms while orthonormalizing ... \n";

        integrate ( elements ( dispETFESpace->mesh() ),
                    quadRule(),
                    dispETFESpace,
                    dot (  P, grad (phi_i) )
                  ) >> residualVectorPtr;
    }
    else
    {

        auto f_0 = _v0 (dispETFESpace, fibers);
        auto s_0 = _v0 (dispETFESpace, sheets);

        boost::shared_ptr<orthonormalizeFibers> normalize0 (new orthonormalizeFibers);
        boost::shared_ptr<orthonormalizeFibers> normalize1 (new orthonormalizeFibers (1) );

        auto f0 = eval (normalize0, f_0);
        auto s0 = eval (normalize1, s_0);

        auto P = eval (W8, _I8 ( F, f0, s0 ) )
                 * _dI8 (F, f0, s0);

    	//if(disp.comm().MyPID() == 0)
        //std::cout << "EMETA - Computing I8 residual terms ... \n";

        integrate ( elements ( dispETFESpace->mesh() ),
                    quadRule(),
                    dispETFESpace,
                    dot (  P, grad (phi_i) )
                  ) >> residualVectorPtr;
    }


}

template< typename Mesh, typename FunctorPtr >
void
computeI8ResidualTermsFung ( const vector_Type& disp,
                             boost::shared_ptr<ETFESpace<Mesh, MapEpetra, 3, 3 > >  dispETFESpace,
                             const vector_Type& fibers,
                             const vector_Type& sheets,
                             vectorPtr_Type           residualVectorPtr,
                             FunctorPtr                  W8,
                             Int Case = 0,
                             bool orthonormalize = true)
{
    using namespace ExpressionAssembly;
	auto F = _F (dispETFESpace, disp, 0);

    if (orthonormalize)
    {
        auto f_0 = _v0 (dispETFESpace, fibers);
        auto s_0 = _v0 (dispETFESpace, sheets);

        boost::shared_ptr<orthonormalizeFibers> normalize0 (new orthonormalizeFibers);
        boost::shared_ptr<orthonormalizeFibers> normalize1 (new orthonormalizeFibers (1) );
        auto f0 = eval (normalize0, f_0);

        auto s_00 = s_0 - dot (f0, s_0) * f0;

        auto s0 = eval (normalize1, s_00);

        if (Case == 0)
        {
            auto P = eval (W8, _I8 ( F, f0, s0 ) )
                     * _dI8 ( F, f0, s0);

        	if(disp.comm().MyPID() == 0)
            std::cout << "EMETA - Computing I8 fs residual terms Fung orthonormalize... \n";


            integrate ( elements ( dispETFESpace->mesh() ),
                        quadRule(),
                        dispETFESpace,
                        dot (  P, grad (phi_i) )
                      ) >> residualVectorPtr;
        }
        else
        {
            boost::shared_ptr<CrossProduct> wedge (new CrossProduct);
            auto n0 = eval ( wedge, f0, s0);

            if (Case == 1)
            {
                auto P = eval (W8, F, f0, s0 )
                         * _dI8 ( F, f0, n0);

            	if(disp.comm().MyPID() == 0)
                std::cout << "EMETA - Computing I8 fn residual terms Fung orthonormalize ... \n";


                integrate ( elements ( dispETFESpace->mesh() ),
                            quadRule(),
                            dispETFESpace,
                            dot (  P, grad (phi_i) )
                          ) >> residualVectorPtr;
            }
            else
            {
                auto P = eval (W8, F, f0, s0 )
                         * _dI8 (F, s0, n0);

            	if(disp.comm().MyPID() == 0)
                std::cout << "EMETA - Computing I8 sn residual terms Fung orthonormalize ... \n";


                integrate ( elements ( dispETFESpace->mesh() ),
                            quadRule(),
                            dispETFESpace,
                            dot (  P, grad (phi_i) )
                          ) >> residualVectorPtr;

            }
        }
    }
    else
    {

        auto f_0 = _v0 (dispETFESpace, fibers);
        auto s_0 = _v0 (dispETFESpace, sheets);

        boost::shared_ptr<orthonormalizeFibers> normalize0 (new orthonormalizeFibers);
        boost::shared_ptr<orthonormalizeFibers> normalize1 (new orthonormalizeFibers (1) );
        boost::shared_ptr<orthonormalizeFibers> normalize2 (new orthonormalizeFibers (2) );

        auto f0 = eval (normalize0, f_0);
        auto s0 = eval (normalize1, s_0);

        if (Case == 0)
        {
            auto P = eval (W8, _I8 (F, f0, s0 ) )
                     * _dI8 (F, f0, s0);

        	if(disp.comm().MyPID() == 0)
            std::cout << "EMETA - Computing I8 fs residual terms Fung ... \n";


            integrate ( elements ( dispETFESpace->mesh() ),
                        quadRule(),
                        dispETFESpace,
                        dot (  P, grad (phi_i) )
                      ) >> residualVectorPtr;
        }
        else
        {
            boost::shared_ptr<CrossProduct> wedge (new CrossProduct);
            auto n_0 = eval ( wedge, f0, s0);
            auto n0 = eval (normalize2, n_0);


            if (Case == 1)
            {
                auto P = eval (W8, F , f0, s0 )
                         * _dI8 (F, f0, n0);

            	if(disp.comm().MyPID() == 0)
                std::cout << "EMETA - Computing I8 fn residual terms Fung ... \n";


                integrate ( elements ( dispETFESpace->mesh() ),
                            quadRule(),
                            dispETFESpace,
                            dot (  P, grad (phi_i) )
                          ) >> residualVectorPtr;
            }
            else
            {
                auto P = eval (W8, F, f0, s0 )
                         * _dI8 (F, s0, n0);

            	if(disp.comm().MyPID() == 0)
            		std::cout << "EMETA - Computing I8 sn residual terms Fung ... \n";


                integrate ( elements ( dispETFESpace->mesh() ),
                            quadRule(),
                            dispETFESpace,
                            dot (  P, grad (phi_i) )
                          ) >> residualVectorPtr;

            }
        }

    }


}




}//EMAssembler

}//LifeV

#endif /* EMETAASSMEBLER_HPP_ */
