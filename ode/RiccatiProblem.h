/*
 * File:   RiccatiProblem.h
 * Author: oberhuber
 *
 * Created on February 25, 2016, 10:41 AM
 */

#ifndef RICCATIPROBLEM_H
#define RICCATIPROBLEM_H

#include<cmath>
#include<iostream>
#include<fstream>

class RiccatiProblem
{
   public:

      const int getDegreesOfFreedom() { return 1; }

      void getRightHandSide( const double& t /* parametr */, const double* _u /* pole dlouhe pocet st. volnosti */, double* fu )
      {
         const double& u = _u[ 0 ];
         fu[ 0 ] = pow( t, -4.0 ) * exp( t ) + u + 2.0 * exp( -t ) * u * u;
         //std::cout << " t = " << t << " " << fu[ 0 ] << " ";
      }

      double getExactSolution( const double& t ,
                               const double& c = 1.0 )
      {
         return exp( t ) * ( 1.0 / ( sqrt( 2.0 ) * t * t ) * tan( sqrt( 2.0 ) * ( c - 1.0 / t ) ) - 1.0 / ( 2.0 * t ) );
      }

      bool writeExactSolution( const char*  fileName,
                               const double& initialTime, // a resiciho intervalu
                               const double& finalTime, // b resiciho intervalu
                               const double& timeStep, // ukladame si to pro vykresleni po casovem kroku
                               const double& c = 1.0 )
      {
         std::fstream file;
         file.open( fileName, std::ios::out );
         double t = initialTime;
         while( t < finalTime )
         {
            file << t << " " << this->getExactSolution( t, c ) << std::endl;
            t = std::min( t + timeStep, finalTime );
         }
         return true;
      }
};


#endif /* RICCATIPROBLEM_H */
