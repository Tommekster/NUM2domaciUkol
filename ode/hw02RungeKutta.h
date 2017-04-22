/*
 * File:   hw02RungeKutta.h
 * Author: zikmuto2
 *
 * Created on April 21, 2017, 4:46 PM
 */

#ifndef RUNGEKUTTA_H
#define	RUNGEKUTTA_H

#include "IntegratorBase.h"

template< typename Problem >
class RungeKutta : public IntegratorBase
{
   public:

      RungeKutta( Problem& problem )
      {
         this->k1 = new double[ problem.getDegreesOfFreedom() ];
         this->k2 = new double[ problem.getDegreesOfFreedom() ];
         this->aux = new double[ problem.getDegreesOfFreedom() ];
      }

      bool solve( Problem& problem,
                  double* u )
      {
         const int dofs = problem.getDegreesOfFreedom();
         double tau = std::min( this->integrationTimeStep, this->stopTime - this->time ); // rekl bych, zo toto znacime 'h'
         long int iteration( 0 );
         while( this->time < this->stopTime )
         {
            /****
             * Compute k1
             */
            // k1 = f(t,u)
            problem.getRightHandSide( this->time/* t */, u /* u */, k1 /* result */ );

            /****
             * Compute k2
             */
            // aux = u + tau*k1
            for( int i = 0; i < dofs; i++ )
               aux[ i ] = u[ i ] + tau * k1[ i ] ;

            // k2 = f(t + tau/3, u+tau/3*k1)
            problem.getRightHandSide( this->time + tau /* t + tau */, aux, k2 );

            /***
             * vysledek: \delta u = tau * ( p1*k1(tau) + p2*k2(tau) )
             * p1 = p2 = 1/2
             */
             for( int i = 0; i < dofs; i++ )
                u[ i ] += tau / 5.0 * ( k1[ i ] + k2[ i ] );
             this->time += tau;
             iteration++;
             if( iteration > 100000 )
             {
                std::cerr << "The solver has reached the maximum number of iteratoins. " << std::endl;
                return false;
             }
            //if( this->adaptivity )
            //   tau *= 0.8 * pow( this->adaptivity / eps, 0.2 );
            tau = std::min( tau, this->stopTime - this->time );
            //std::cout << "ITER: " << iteration << " \t tau = " << tau << " \t time= " << time << "         \r " << std::flush;
         }
         //std::cout << std::endl;
         return true;
      }

      ~RungeKutta()
      {
         delete[] k1;
         delete[] k2;
         delete[] aux;
      }

   protected:

      double *k1, *k2, *aux;

};

#endif	/* RUNGEKUTTA_H */
