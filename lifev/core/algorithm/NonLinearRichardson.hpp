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
    @brief Preconditioned relaxed solver for non linear problems.

    @author
    @date 01-01-2004

    Given a functional this method tries to find a root near an
    initial guess sol. The two methods that the functional must
    have are evalRes and solveJac. If solveJac inverts the exact
    Jacobian of the functional, than this is equivalent to a Newton
    method. If convergence do not appens in the expected ratio,
    a NonLinearLineSearch algorithm is called (uadratic or parabolic).

    @contributor Simone Deparis <simone.deparis@epfl.ch>
    @maintainer Simone Deparis <simone.deparis@epfl.ch>


 */

#ifndef _NONLINEARRICHARDSON_HPP
#define _NONLINEARRICHARDSON_HPP

#include <algorithm> // for min and max
#include <lifev/core/algorithm/NonLinearLineSearch.hpp>
#include <lifev/core/array/VectorEpetra.hpp>

namespace LifeV
{
//! Preconditioned relaxed solver for non linear problems.
/*!
    Add more details about the constructor.
    NOTE: short description is automatically added before this part.
    @param sol            :  the solution
    @param maxit          :  input: maximum iterations, output: nb of iterations
    @param abstol, reltol :  the stoping criteria is abstol+reltol*norm(residual_0),
    @param eta_max        :  Maximum error tolerance for residual in linear solver.

    The linear solver terminates when the relative
    linear residual is smaller than eta*| f(sol) |.

    The value linear_rel_tol send for the relative tolerance
    to the linear solver is therefore eta. eta is determined
    by the modified Eisenstat-Walker formula if etamax > 0.

    @param NonLinearLineSearch     :  for now consider only the case NonLinearLineSearch=0
      (coded but not theoretically analysed)

    @param omega          :  default relaxation parameter to be passed to Aitken.
    if omega is negative, then its absolute value is
    taken as constant relaxation parameter

 */

template < class Fct >
Int NonLinearRichardson ( VectorEpetra& sol,
                          Fct&        functional,
                          Real        abstol,
                          Real        reltol,
                          UInt&       maxit,
                          Real        eta_max,
                          Int         NonLinearLineSearch,
                          UInt iter = UInt (0),
                          UInt        verboseLevel = 0,
                          std::ostream& output = std::cout,
                          const Real& time = 0
                        )
{
    /*x
        */

    /*
      max_increase_res: maximum number of successive increases in residual
      before failure is reported
    */

    //    const Int max_increase_res = 5;

    /*
      Parameters for the linear solver, gamma: Default value = 0.9
    */

    const Real gamma   = 0.9;

    //----------------------------------------------------------------------

    if ( sol.comm().MyPID() != 0 )
    {
        verboseLevel = 0;
    }

    VectorEpetra residual ( sol.map() );
    VectorEpetra step     ( sol.map() );

    step *= 0.;

    Real normResOld = 1;

    if ( verboseLevel > 0 )
    {
        std::cout << std::endl;
        std::cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
        std::cout << "      Non-Linear Richardson: starting          " << std::endl;
        std::cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
        std::cout << std::endl;
    }

    functional.evalResidual ( residual, sol, iter );

    Real normRes      = residual.normInf();
    Real stop_tol     = abstol + reltol * normRes;
    Real linearRelTol = std::fabs (eta_max);
    Real eta_old;
    Real eta_new;
    Real ratio;
    Real slope;
    Real linres;
    Real lambda;

    //

    Real solNormInf (sol.normInf() );
    Real stepNormInf;
    if ( verboseLevel > 1 )
    {
        output << std::scientific;
        output << "# time = ";
        output << time << "   " << "initial norm_res " <<  normRes
               << " stop tol = " << stop_tol
               << "initial norm_sol "
               << solNormInf << std::endl;
        output << "#iter      disp_norm       step_norm       residual_norm" << std::endl;
    }
    while ( normRes > stop_tol && iter < maxit )
    {
        if ( verboseLevel > 0 )
        {
            std::cout << std::endl;
            std::cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
            std::cout << "      Non-Linear Richardson: iteration  =      " << iter << std::endl
                      << "                             residual   =      " << normRes << std::endl
                      << "                             tolerance  =      " << stop_tol << std::endl;
            std::cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
            std::cout << std::endl;
        }

        if ( iter > (maxit-2) )
        {
            break;
            
            if ( verboseLevel > 0 )
            {
                std::cout << std::endl;
                std::cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
                std::cout << " Non-Linear Richardson: desired tolerance level not accomplished, go to next timestep" << std::endl;
                std::cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
                std::cout << std::endl;
            }
        }
        
        iter++;

        ratio      = normRes / normResOld;
        normResOld = normRes;
        normRes    = residual.normInf();

        residual *= -1;

        functional.solveJac (step, residual, linearRelTol); // J*step = -R

        solNormInf = sol.normInf();
        stepNormInf = step.normInf();
        if ( verboseLevel > 1 )
        {
            output   << std::setw (5) << iter
                     << std::setw (15) << solNormInf
                     << std::setw (15) << stepNormInf;
        }
        linres = linearRelTol;
	

	////////////////////////////////////////////////HERE I HAVE CHANGED
	// lambda = 1.;
	// /////////////////////////////////////////////
	//lambda = 0.1;
	lambda = 0.2; //2020.03.25: 0.2 -> 1
	//lambda = 0.2;//previously 0.4 lg

        slope  = normRes * normRes * ( linres * linres - 1 );

        Int status (EXIT_SUCCESS);
        switch ( NonLinearLineSearch )
        {
            case 0: // no NonLinearLineSearch
                //if ( iter > 20 ) lambda = 0.6;
                if (iter > 7 ) lambda = 1.0;
                if(iter > 30) lambda = 1.0;
                if(normRes < 2000 && iter > 5) lambda = 1.0;//>5 war >3 lg
                else if(iter > 17 && iter < 25 ) lambda = 0.5;//>17 war >10 lg
                else if(iter > 25) lambda = 1.0;
                else if(normRes < 20000 && iter > 20) lambda = 0.6;

                //if(iter > 25) lambda = 1.0;
                //if(normRes < 4000 && iter > 20) lambda = 1.0; //Hier aendern zu lambda = 0.6
                //if(normRes < 2000 && iter > 10) lambda = 1.0;
                //if(normRes < 4000 || iter > 40 ) lambda = 1.0;
                sol += lambda * step;
                functional.evalResidual ( residual, sol, iter);
                //                normRes = residual.NormInf();
                break;
            case 1:
                status = NonLinearLineSearchParabolic ( functional, residual, sol, step, normRes, lambda, iter, verboseLevel );
                break;
            case 2:  // recommended
                status = NonLinearLineSearchCubic ( functional, residual, sol, step, normRes, lambda, slope, iter, verboseLevel );
                break;
            default:
                std::cout << "Unknown NonLinearLineSearch \n";
                status = EXIT_FAILURE;
        }

        if (status == EXIT_FAILURE)
        {
            return status;
        }

        normRes = residual.normInf();

        if ( verboseLevel > 1 )
        {
            output << std::setw (15) << normRes << std::endl;
        }

        if ( eta_max > 0 )
        {
            eta_old = linearRelTol;
            eta_new = gamma * ratio * ratio;
            if ( gamma * eta_old * eta_old > .1 )
            {
                eta_new = std::max<Real> ( eta_new, gamma * eta_old * eta_old );
            }
            linearRelTol = std::min<Real> ( eta_new, eta_max );
            linearRelTol = std::min<Real> ( eta_max,
                                            std::max<Real> ( linearRelTol,
                                                             .5 * stop_tol / normRes ) );
            //if ( verboseLevel > 0 )
            //    std::cout << "    Newton: forcing term eta = " << linearRelTol << std::endl;
            //std::cout<<"\nVerbose = "<<verboseLevel<<std::endl;
        }

    }

    if ( normRes > stop_tol )
    {
        if ( verboseLevel > 0 )
        {
            std::cout << "!!! NonLinRichardson: convergence fails" << std::endl;
        }
        maxit = iter;
        return EXIT_FAILURE;
    }


    if ( verboseLevel > 0 )
    {
        std::cout << std::endl;
        std::cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
        std::cout << "      Non-Linear Richardson: convergence =     " << normRes << std::endl
                  << "                             iterations  =     " << iter << std::endl;
        std::cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
        std::cout << std::endl;
    }
    maxit = iter;

    return EXIT_SUCCESS;
}
}
#endif
