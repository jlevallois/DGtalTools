/**
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 *  License, or  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 **/
/**
 * @file 3dCurvatureViewer.cpp
 * @ingroup surfaceTools
 * @author Jérémy Levallois (\c jeremy.levallois@liris.cnrs.fr )
 * Laboratoire d'InfoRmatique en Image et Systèmes d'information - LIRIS (CNRS, UMR 5205), INSA-Lyon, France
 * LAboratoire de MAthématiques - LAMA (CNRS, UMR 5127), Université de Savoie, France
 *
 * @date 2013/01/10
 *
 * Vol file viewer, with curvature (mean or Gaussian, see parameters) information on surface.
 * Blue color means lowest curvature
 * Yellow color means highest curvature
 * Red means the in-between
 *
 * Uses IntegralInvariantCurvatureEstimation
 * @see related article:
 *       Coeurjolly, D.; Lachaud, J.O; Levallois, J., (2013). Integral based Curvature
 *       Estimators in Digital Geometry. DGCI 2013. Retrieved from
 *       https://liris.cnrs.fr/publis/?id=5866
 *
 * This file is part of the DGtalTools.
 */

///////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/math/KMeans.h"
#include <cstring>

using namespace DGtal;

///////////////////////////////////////////////////////////////////////////////

template <typename Point>
struct Euclidean
{
  double distance(const Point &a, const Point &b) const
  {
    return (double)(b-a).norm();
  }
};


template< typename P >
void myKMeans( P histogram[10], unsigned int seeds[3] )
{
  unsigned int size_seeds = 3;
  unsigned int size_histogram = 10;

  for( unsigned int ii = 0; ii < size_seeds; ++ii )
    {
      histogram[ seeds[ ii ]].second = ii;
    }

  unsigned int tries = 0;
  while( tries < 10 )
    {
      for( unsigned int ii = 0; ii < size_histogram; ++ii )
        {
          unsigned int min_seed = size_seeds;
          unsigned int min = size_histogram * size_histogram;
          for( unsigned int jj = 0; jj < size_seeds; ++jj )
            {
              unsigned int distance = ( seeds[ jj ] - ii ) * ( seeds[ jj ] - ii );
              if( distance < min )
                {
                  min = distance;
                  min_seed = ii;
                }
            }
          histogram[ii].second = min_seed;
        }

      unsigned int current_seed = 0;
      unsigned int begin = 0;
      unsigned int cumul = 0;

      for( unsigned int ii = 0; ii < size_histogram; ++ii )
        {
          if( histogram[ii].second != current_seed )
            {
              double split = cumul / 2.0;
              unsigned int temp_cumul = 0;
              for( ; begin < ii; ++begin )
                {
                  temp_cumul += histogram[begin].first;
                  if( temp_cumul > split )
                    {
                      seeds[current_seed] = begin;
                    }
                }
              cumul = histogram[ii].first;
              ++current_seed;
            }
          else
            {
              cumul += histogram[ii].first;
            }
        }
      {
        double split = cumul / 2.0;
        unsigned int temp_cumul = 0;
        for( ; begin < size_histogram; ++begin )
          {
            temp_cumul += histogram[begin].first;
            if( temp_cumul > split )
              {
                seeds[current_seed] = begin;
              }
          }
      }

      ++tries;
    }
}

int main( int argc, char** argv )
{
  unsigned int size_histogram = 10;
  typedef std::pair< unsigned int, unsigned int > P;
  P histogram[ 10 ] = { P(4,0), P(2,0), P(0,0), P(8,0), P(9,0), P(5,0), P(5,0), P(3,0), P(4,0), P(1,0) };
  typedef Z2i::Point Point;
  std::vector< Point > p_histogram;
  {
    for( int ii = 0; ii < 4; ++ii )
      p_histogram.push_back( Point( 1, 0 ));

    for( int ii = 0; ii < 2; ++ii )
      p_histogram.push_back( Point( 2, 0 ));

    for( int ii = 0; ii < 0; ++ii )
      p_histogram.push_back( Point( 3, 0 ));

    for( int ii = 0; ii < 8; ++ii )
      p_histogram.push_back( Point( 4, 0 ));

    for( int ii = 0; ii < 9; ++ii )
      p_histogram.push_back( Point( 5, 0 ));

    for( int ii = 0; ii < 5; ++ii )
      p_histogram.push_back( Point( 6, 0 ));

    for( int ii = 0; ii < 5; ++ii )
      p_histogram.push_back( Point( 7, 0 ));

    for( int ii = 0; ii < 3; ++ii )
      p_histogram.push_back( Point( 8, 0 ));

    for( int ii = 0; ii < 4; ++ii )
      p_histogram.push_back( Point( 9, 0 ));

    for( int ii = 0; ii < 1; ++ii )
      p_histogram.push_back( Point( 10, 0 ));
  }

  const unsigned int size_p_histogram = p_histogram.size();
  std::vector< unsigned int > reg;
  std::vector< Point > seeds;
  unsigned int size_seeds = 3;

  KMeans( p_histogram, size_seeds, Euclidean<Point>(), reg, seeds);

  for( unsigned int ii = 0; ii < size_seeds ; ++ii )
    {
      std::cout << "seed " << seeds[ii] << std::endl;

      for( unsigned jj = 0; jj < size_p_histogram; ++jj )
        {
          if( reg[jj] == ii )
            {
              std::cout << "point " << p_histogram[jj] << std::endl;
            }
        }
    }
  /*unsigned int size_seeds = 3;
  unsigned int seeds[ 3 ] = { 0, 5, 9 };

  myKMeans< P >( histogram, seeds );

  for( unsigned int ii = 0; ii < size_seeds; ++ii )
  {
    std::cout << "center[" << ii << "] = " << seeds[ ii ] << " aka value = " << histogram[ seeds[ ii ]].first << std::endl;
  }*/
}

///////////////////////////////////////////////////////////////////////////////
