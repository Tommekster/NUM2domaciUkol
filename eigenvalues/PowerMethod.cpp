#include "PowerMethod.h"

PowerMethod::PowerMethod( const Matrix& A )
:  A( A ), max_iterations( 100 ), convergence_residue( 1.0e-8 )
{
}

void PowerMethod::setMaxIterations( const int max_iterations )
{
   this->max_iterations = max_iterations;
}

void PowerMethod::setConvergenceResidue( const Real& convergence_residue )
{
   this->convergence_residue = convergence_residue;
}

bool PowerMethod::getLargestEigenvalue( Vector& x, Real& lambda, int verbose )
{
   int iteration( 0 );
   int norm_index( 0 );
   Real rho, residue( this->convergence_residue + 1.0 );
   Vector y, r;
   y.setSize( x.getSize() );
   r.setSize( x.getSize() );
   while( iteration < this->max_iterations )
   {
      A.vectorMultiplication( x, y );
      if( y[ norm_index ] = 0.0 )
      {
         norm_index++;
         if( norm_index == y.getSize() )
            norm_index = 0;
         continue;
      }
      rho = y[ norm_index ];
      for( int i = 0; i < x.getSize(); i++ )
         x[ i ] = y[ i ] / rho;
      if( iteration % 10 == 0 )
      {
         A.getEigenvalueResidue( x, rho, r );
         residue = r.maxNorm();
      }
      if( residue < this->convergence_residue )
       return true;
   }
   return false;
}

bool PowerMethod::getSmallestEigenvalue( Vector& x, Real& lambda, int verbose )
{
 
}

