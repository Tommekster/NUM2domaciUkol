/* 
 * File:   denseMatrix.cpp
 * Author: oberhuber
 * 
 * Created on September 28, 2016, 5:31 PM
 */

#include <string>
#include <sstream>
#include <iomanip>
#include <assert.h>
#include "DenseMatrix.h"

DenseMatrix::DenseMatrix()
: Matrix( 0 ,0 )
{
}
DenseMatrix::DenseMatrix( const int rows, const int columns )
: Matrix( rows, columns )
{
   this->elements.resize( rows * columns );
}
      
bool DenseMatrix::setDimensions( const int rows, const int columns )
{
   this->rows = rows;
   this->columns = columns;
   this->elements.resize( rows * columns );
   return true;
}

   
Real& DenseMatrix::operator()( const int row, const int column )
{
   return this->elements[ row * this->columns + column ];
}
      
const Real& DenseMatrix::operator()( const int row, const int column ) const
{
   return this->elements[ row * this->columns + column ];
}

void DenseMatrix::vectorMultiplication( const Vector& in_vector,
                                        Vector& out_vector ) const
{
   assert( in_vector.getSize() == this->columns );
   assert( out_vector.getSize() == this->rows );
   
   int idx( 0 );
   for( int row = 0; row < this->rows; row++ )
   {
      Real result( 0.0 );
      for( int column = 0; column < this->columns; column++ )
         result += this->elements[ idx++ ] * in_vector[ column ];
      out_vector[ row ]= result;
   }
}

void DenseMatrix::performJacobiIteration( const Vector& b,
                                          const Vector& x,
                                          Vector& aux ) const
{
   assert( x.getSize() == this->columns );
   assert( b.getSize() == this->columns );
   assert( this->columns == this->rows );
   
   int idx( 0 );
   
   for( int row = 0; row < this->rows; row++ )
   {
      Real sum( 0.0 ), a_ii( 0.0 );
      for( int column = 0; column < this->columns; column++ )
         if( column != row )
            sum += this->elements[ idx++ ] * x[ column ];
         else
            a_ii = this->elements[ idx++ ];
      aux[ row ]= 1.0 / a_ii * ( b[ row ] - sum );
   }
   
}

DenseMatrix& DenseMatrix::operator=( const DenseMatrix& m )
{
   this->setDimensions( m.getRows(), m.getColumns() );
   this->elements = m.elements;
}

DenseMatrix& DenseMatrix::operator-=( const DenseMatrix& m )
{
   for( int i = 0; i < this->elements.size(); i++ )
      this->elements[ i ] -= m.elements[ i ];
}
