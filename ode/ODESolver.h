/*
 * File:   ODESolver.h
 * Author: oberhuber
 *
 * Created on February 25, 2016, 8:48 AM
 */

#ifndef ODESOLVER_H
#define	ODESOLVER_H

#include <iostream>
#include <cmath>
#include "ODESolution.h"

template< typename Problem,
          typename Integrator >
class ODESolver
{
   public:

      ODESolver( Problem& problem,
                 Integrator& integrator )
      : problem( problem ),
        integrator( integrator )
      {
         this->u = new double[ problem.getDegreesOfFreedom() ];
      }

      void setInitialCondition( const double* initialCondition )
      {
         for( int i = 0; i < problem.getDegreesOfFreedom(); i++ )
            u[ i ] = initialCondition[ i ];
      }

      bool solve( ODESolution& solution,
                  const double& initialTime,
                  const double& finalTime,
                  const double& timeStep )
      {
         const int timeStepsCount = std::ceil( std::max( 0.0, finalTime - initialTime ) / timeStep );
         solution.setup( problem.getDegreesOfFreedom(),
                         timeStepsCount );
         // na pocatku se do reseni vlozi pocatecni podminka
         for( int dof = 0; dof < this->problem.getDegreesOfFreedom(); dof++ )
            solution.setElement( 0, dof, u[ dof ] );

         // samotne reseni
         double time( initialTime );
         for( int k = 1; k <= timeStepsCount; k++ )
         {
            std::cout << "Solving time step " << k << " / " << timeStepsCount << " => " <<
               ( double ) k / ( double ) timeStepsCount * 100.0 << "% " << std::endl;

            // tady se resi samotny krok
            this->integrator.setTime( time );
            this->integrator.setStopTime( time + timeStep );
            if( ! this->integrator.solve( problem, u ) )
               return false;

            // krok hezky po slozkach ulozime
            for( int dof = 0; dof < this->problem.getDegreesOfFreedom(); dof++ )
               solution.setElement( k, dof, u[ dof ] );
            time += timeStep;
         }
         std::cout << "Done." << std::endl;
         return true;
      }

   protected:

      Problem& problem;

      Integrator& integrator;

      double* u;
};


#endif	/* ODESOLVER_H */
