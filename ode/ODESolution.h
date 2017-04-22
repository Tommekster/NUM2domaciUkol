/*
 * File:   ODESolution.h
 * Author: oberhuber
 *
 * Created on February 25, 2016, 8:37 AM
 */

#ifndef ODESOLUTION_H
#define	ODESOLUTION_H

#include<assert.h>
#include<fstream>

class ODESolution
{
   public:

      ODESolution()
         : data( 0 ), degreesOfFreedom( 0 ), timeStepsCount( 0 ) {};

      ODESolution( const int degreesOfFreedom,
                   const int timeStepsCount )
      {
         setup( degreesOfFreedom, timeStepsCount );
      }

      bool setup( const int degreesOfFreedom,
                  const int timeStepsCount )
      {
         this->degreesOfFreedom = degreesOfFreedom;
         this->timeStepsCount = timeStepsCount;
         data = new double[ degreesOfFreedom * ( timeStepsCount + 1 ) ];
         if( ! data )
            return false;
         return true;
      }

      void setElement( const int timeStep,
                       const int dof, // degreesOfFreedom
                       const double& value )
      {
         assert( timeStep <= this->timeStepsCount );
         assert( dof < this->degreesOfFreedom );
         this->data[ timeStep * this->degreesOfFreedom + dof ] = value;
      }

      double getElement( const int timeStep,
                         const int dof )
      {
         assert( timeStep <= this->timeStepsCount );
         assert( dof < this->degreesOfFreedom );
         // vektory jsou ukladany po slozkach za sebou [a0, a1, .. , a_n, b0, b1, ... ]
         return this->data[ timeStep * this->degreesOfFreedom + dof ];
      }

      bool write( const char* fileName,
                  const double& initialTime = 0.0,
                  const double& timeStep = 0.0,
                  const int dofToWrite = -1 )
      {
         assert( dofToWrite == -1 || ( dofToWrite >= 0 && dofToWrite < degreesOfFreedom ) );
         std::fstream file;
         file.open( fileName, std::ios::out );
         if( !file )
         {
            std::cerr << "Unable to open a file " << fileName << "." << std::endl;
            return false;
         }
         double time( initialTime );
         for( int k = 0; k <= timeStepsCount; k++ )
         {
            if( timeStep != 0.0 )
               file << time << " ";
            time += timeStep;
            if( dofToWrite == -1 )
               for( int dof = 0; dof < degreesOfFreedom; dof++ )
                  file << this->getElement( k, dof ) << " ";
            else
               file << this->getElement( k, dofToWrite );
            file << std::endl;
         }
         return true;
      }

      ~ODESolution()
      {
         delete[] data;
      }

   protected:

      double* data;

      int degreesOfFreedom;

      int timeStepsCount;
};

#endif	/* ODESOLUTION_H */
