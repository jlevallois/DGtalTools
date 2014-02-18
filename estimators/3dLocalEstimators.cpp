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
 * @file 3dLocalEstimators.cpp
 * @ingroup Tools
 * @author Jeremy Levallois (\c jeremy.levallois@liris.cnrs.fr )
 * Laboratoire d'InfoRmatique en Image et Systèmes d'information - LIRIS (CNRS, UMR 5205), Universite de Lyon, France
 * LAboratoire de MAthematiques - LAMA (CNRS, UMR 5807), Universite de Savoie, France
 *
 * @date 2012/06/20
 *
 * DGtal 3D curvature shape comparator
 *
 * This file is part of the DGtalTools library.
 */

///////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <iomanip>
#include <string>

#include "DGtal/base/Common.h"
#include "DGtal/base/Clock.h"
#include "DGtal/helpers/StdDefs.h"

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

#include "DGtal/math/KMeans.h"

//shapes
#include "DGtal/shapes/implicit/ImplicitBall.h"
#include "DGtal/math/MPolynomial.h"
#include "DGtal/io/readers/MPolynomialReader.h"
#include "DGtal/shapes/implicit/ImplicitPolynomial3Shape.h"

//Digitizer
#include "DGtal/shapes/GaussDigitizer.h"
#include "DGtal/topology/LightImplicitDigitalSurface.h"
#include "DGtal/topology/DigitalSurface.h"
#include "DGtal/graph/DepthFirstVisitor.h"
#include "DGtal/geometry/volumes/KanungoNoise.h"
#include "DGtal/topology/CanonicSCellEmbedder.h"

#include "DGtal/math/Statistic.h"

// Segments
#include "DGtal/geometry/curves/SaturatedSegmentation.h"
#include "DGtal/geometry/curves/ArithmeticalDSSComputer.h"
#include "DGtal/topology/DigitalSurface2DSlice.h"
#include "DGtal/geometry/surfaces/ChordGenericNaivePlaneComputer.h"

//Estimators
#include "DGtal/kernel/BasicPointFunctors.h"
#include "DGtal/geometry/surfaces/FunctorOnCells.h"

#include "DGtal/geometry/curves/estimation/TrueLocalEstimatorOnPoints.h"
#include "DGtal/geometry/surfaces/estimation/IntegralInvariantMeanCurvatureEstimator.h"
#include "DGtal/geometry/surfaces/estimation/IntegralInvariantGaussianCurvatureEstimator.h"

#include "DGtal/geometry/surfaces/estimation/LocalEstimatorFromSurfelFunctorAdapter.h"
#include "DGtal/geometry/surfaces/estimation/estimationFunctors/MongeJetFittingMeanCurvatureEstimator.h"
#include "DGtal/geometry/surfaces/estimation/estimationFunctors/MongeJetFittingGaussianCurvatureEstimator.h"
#include "DGtal/geometry/surfaces/estimation/estimationFunctors/MongeJetFittingPrincipalCurvaturesEstimator.h"

using namespace DGtal;

struct OptionsIntegralInvariant
{
  double constante;
  std::string tests;
  unsigned int nbKernels;
  std::string modeSegments;
  std::string typeSegments;
};


template < typename Shape, typename KSpace, typename ConstIterator, typename OutputIterator >
void
estimateTrueMeanCurvatureQuantity( const ConstIterator & it_begin,
                                   const ConstIterator & it_end,
                                   OutputIterator & output,
                                   const KSpace & K,
                                   const double & h,
                                   Shape * aShape )
{
  typedef typename KSpace::Space::RealPoint RealPoint;
  typedef CanonicSCellEmbedder< KSpace > Embedder;

  Embedder embedder( K );
  RealPoint currentRealPoint;

  for ( ConstIterator it = it_begin; it != it_end; ++it )
  {
    currentRealPoint = embedder( *it_begin ) * h;
    output = aShape->meanCurvature( currentRealPoint );
    ++output;
  }
}

template < typename Shape, typename KSpace, typename ConstIterator, typename OutputIterator >
void
estimateTrueGaussianCurvatureQuantity( const ConstIterator & it_begin,
                                       const ConstIterator & it_end,
                                       OutputIterator & output,
                                       const KSpace & K,
                                       const double & h,
                                       Shape * aShape )
{
  typedef typename KSpace::Space::RealPoint RealPoint;
  typedef CanonicSCellEmbedder< KSpace > Embedder;

  Embedder embedder( K );
  RealPoint currentRealPoint;

  for ( ConstIterator it = it_begin; it != it_end; ++it )
  {
    currentRealPoint = embedder( *it_begin ) * h;
    output = aShape->gaussianCurvature( currentRealPoint );
    ++output;
  }
}

template < typename Shape, typename KSpace, typename ConstIterator, typename OutputIterator >
void
estimateTruePrincipalCurvaturesQuantity( const ConstIterator & it_begin,
                                         const ConstIterator & it_end,
                                         OutputIterator & output,
                                         const KSpace & K,
                                         const double & h,
                                         Shape * aShape )
{
  typedef typename KSpace::Space::RealPoint RealPoint;
  typedef CanonicSCellEmbedder< KSpace > Embedder;

  Embedder embedder( K );
  RealPoint currentRealPoint;

  for ( ConstIterator it = it_begin; it != it_end; ++it )
  {
    currentRealPoint = embedder( *it_begin ) * h;
    double k1, k2;
    aShape->principalCurvatures( currentRealPoint, k1, k2 );
    CurvatureInformations result;
    result.k1 = k1;
    result.k2 = k2;
    output = result;
    ++output;
  }
}

/**
 * return the position inside surfels, size of surfels else
 */
template< typename Surfel >
Dimension findSurfel( const std::vector< Surfel > & surfels,
                      const Surfel & s )
{
  const Dimension surfels_size = surfels.size();
  bool found = false;
  Dimension position = surfels_size;

  for( Dimension ii = 0; !found && ii < surfels_size; ++ii )
  {
    if( surfels[ii] == s )
    {
      found = true;
      position = ii;
    }
  }

  return position;
}


template <typename K3D, typename K2D>
class SCellToMyPoint
{
public:
  typedef K3D KSpace3D;
  typedef K2D KSpace2D;
  typedef typename Z2i::Point Output;
  typedef Output Value;
  typedef typename KSpace3D::SCell Input;
  typedef typename KSpace3D::Point Point;

private:
  /**
       * Aliasing pointer on the Khalimsky space.
      */
  const KSpace3D* myK;
  Dimension i;

public:

  /**
       * Default constructor.
      */
  SCellToMyPoint() : myK(NULL), i(0) { }
  /**
       * Constructor.
       * @param aK a Khalimsky space
      */
  SCellToMyPoint(const KSpace3D& aK, Dimension myDimension) : myK(&aK), i(myDimension) { }

  /**
     * Copy constructor.
     * @param other any SCellToPoint functor
     */
  SCellToMyPoint(const SCellToMyPoint& other)
    : myK(other.myK), i(other.i) { }

  /**
     * Assignment.
     *
     * @param other the object to copy.
     * @return a reference on 'this'.
     */
  SCellToMyPoint& operator= ( const SCellToMyPoint & other )
  {
    if (this != &other)
    {
      myK = other.myK;
      i = other.i;
    }
    return *this;
  }

  /**
     * Returns a point (with integer coordinates)
     * from a scell (with khalimsky coordinates)
     * @param aSCell a scell
     * @return the corresponding point.
     */
  Output operator()(const Input& aSCell) const
  {
    Input s = aSCell;
    Dimension k = myK->sOrthDir( s );
    Dimension j = (k+1)%3;
    if ( j == i ) j = (i+1)%3;
    Input next_linel = myK->sDirectIncident( s, j );
    Input base_pointel = myK->sIncident( next_linel, i, false );
    Point p = myK->sCoords( base_pointel );
    Output q( p[(i+1)%3], p[(i+2)%3] );
    return q;
  }

}; // end of class SCellToMyPoint

template< typename KSpace, typename Iterator >
void analyseAllLengthMS( std::vector< Statistic<double> > & statE,
                         Iterator itb,
                         Iterator ite )
{
  typedef typename KSpace::Space Space;
  typedef typename Space::Point Point;
  typedef typename Space::Vector Vector;
  typedef ArithmeticalDSSComputer< Iterator, int, 4 > SegmentComputer;
  typedef SaturatedSegmentation< SegmentComputer > Decomposition;
  typedef typename Decomposition::SegmentComputerIterator SegmentComputerIterator;
  typedef std::vector< SegmentComputer > VectorOfSegmentComputer;
  typedef std::map< Point, VectorOfSegmentComputer > Pmap;

  // Computes the tangential cover
  SegmentComputer algo;
  Iterator itbegin = itb;
  Iterator itend = ite;
  Decomposition theDecomposition( itbegin, itend, algo);

  Pmap map;
  //for( itbegin = itb; itbegin != itend; ++itbegin )
  do
  {
    map.insert( std::pair< Point, VectorOfSegmentComputer >( *itbegin, VectorOfSegmentComputer() ) );
    ++itbegin;
  } while( itbegin != itend );


  for ( SegmentComputerIterator scIt = theDecomposition.begin(), scItEnd = theDecomposition.end();
        scIt != scItEnd; ++scIt )
  {
    const SegmentComputer & sc = *scIt;
    for ( Iterator ptIt = sc.begin(), ptItEnd = sc.end(); ptIt != ptItEnd; ++ptIt )
    {
      typename Pmap::iterator mloc = map.find( *ptIt );
      if( mloc != map.end() )
      {
        mloc->second.push_back( sc );
      }
      else
      {
        trace.error() << "not found ?" << std::endl;
      }
    }
  }

  itbegin = itb;
  itend = ite;
  Dimension ii = 0;
  //for( itbegin = itb; itbegin != itend; ++itbegin )
  do
  {
    //statD[ii].clear();
    statE[ii].clear();
    typename Pmap::iterator mloc = map.find( *itbegin );
    ASSERT(( mloc != map.end() ));

    /////////////
    for( typename VectorOfSegmentComputer::iterator scIt = mloc->second.begin(), scItEnd = mloc->second.end(); scIt != scItEnd; ++scIt )
    {
      const SegmentComputer & sc = *scIt;
      /*int64_t l = 0;
          for ( Iterator ptIt = sc.begin(), ptItEnd = sc.end(); ptIt != ptItEnd; ++ptIt )
            ++l;
          statD[ii].addValue( (double) l );*/
      double v = (sc.back( ) - sc.front()).norm1() ; // sc.size();  //*( sc.end() - 1 ) - *( sc.begin() );
      statE[ii].addValue( v )  ; //v.norm( ) );
      //          std::cout << " v=" << v.norm() << std::endl;
    }

    /////////////

    ++ii;
    ++itbegin;
  } while( itbegin != itend );
}

template< typename ImplicitDigitalSurface, typename Surfel >
void computeDPSSegments( std::vector< Surfel > & surfels,
                         std::map< Surfel, std::pair< Statistic< double >, Statistic< double > > > & segments,
                         const Z3i::KSpace & K,
                         const ImplicitDigitalSurface & impDigitalSurface,
                         const unsigned int widthNum,
                         const double Rmax)
{
  //l_2 graph visitor
  typedef CanonicSCellEmbedder<Z3i::KSpace> VertexEmbedder;
  typedef VertexEmbedder::Value RealPoint;
  typedef RealPoint::Coordinate Scalar;
  typedef ExactPredicateLpSeparableMetric<Z3i::Space,2> Distance;
  typedef std::binder1st< Distance > DistanceToPoint;
  typedef Composer<VertexEmbedder, DistanceToPoint, Scalar> VertexFunctor;
  typedef DistanceBreadthFirstVisitor< ImplicitDigitalSurface, VertexFunctor, std::set< Z3i::KSpace::SCell > > Visitor;
  VertexEmbedder embedder( K );
  Distance distance;

  //DPS
  typedef DGtal::int64_t InternalInteger;
  typedef ChordGenericNaivePlaneComputer<Z3i::Space,Z3i::Point, InternalInteger> PlaneComputer;


  Statistic<double> myStat(true);
  Dimension axis;
  double currentDistance = 0.0;
  Z3i::Point p;
  PlaneComputer planeComputer;

  for ( typename std::vector<Surfel>::const_iterator it = surfels.begin(), itE= surfels.end(); it != itE; ++it )
  {
    Surfel v = *it;
    currentDistance = 0.0;

    DistanceToPoint distanceToPoint = std::bind1st( distance, embedder( v ) );
    VertexFunctor vfunctor( embedder, distanceToPoint );
    axis = K.sOrthDir( v );
    planeComputer.init( widthNum, 1 );
    Visitor visitor( impDigitalSurface, vfunctor, v );
    Surfel vv;
    while ( ! visitor.finished() )
    {
      typename Visitor::Node node = visitor.current();
      vv = node.first;
      axis = K.sOrthDir( vv );
      p = K.sCoords( K.sDirectIncident( vv, axis ) );
      bool isExtended = planeComputer.extend( p );
      if (( isExtended ) && (node.second < Rmax))
      {
        // surfel is in plane.
        visitor.expand();
        currentDistance = node.second;
      }
      else // surfel is not in plane and should not be used in the visit.
        visitor.ignore();
    }

    myStat.addValue(currentDistance);
    segments[ v ] = std::pair<Statistic<double>,Statistic<double> >(myStat, myStat);
    myStat.clear();
  }
}

template< typename ImplicitDigitalSurface, typename Surfel >
void computeSegments( std::vector< Surfel > & surfels,
                      std::map< Surfel, std::pair< Statistic< double >, Statistic< double > > > & segments,
                      const Z3i::KSpace & K,
                      const ImplicitDigitalSurface & impDigitalSurface )
{
  typedef typename ImplicitDigitalSurface::Tracker Tracker;
  typedef DigitalSurface2DSlice< Tracker > MySlice;
  typedef std::pair< Statistic< double >, Statistic< double > > PairOfStatistics;
  typedef std::map< Surfel, PairOfStatistics > SurfelMap;
  typedef std::map< Surfel, bool > MarqueMap;

  const Dimension surfels_size = surfels.size();

  for( Dimension dim = 0; dim < 3; ++dim )
  {
    MarqueMap marque;
    for( Dimension ii = 0; ii < surfels_size; ++ii )
    {
      marque[ surfels[ ii ] ] = false;
    }

    for( Dimension ii = 0; ii < surfels_size; ++ii )
    {
      Surfel currentSurfel = surfels[ ii ];
      if( marque[ surfels[ ii ] ] == true )
      {
        continue;
      }
      else if( K.sOrthDir( currentSurfel ) == dim )
      {
        continue;
      }
      else
      {
        Tracker ptrTracker( impDigitalSurface, currentSurfel );
        Dimension otherdim = 0;
        {
          bool dims[3] = { false, false, false };
          dims[dim] = true;
          dims[K.sOrthDir( currentSurfel )] = true;

          while( dims[otherdim] == true )
          {
            ++otherdim;
          }
        }

        //          unsigned int bitAll = 7;
        /*unsigned int bitSuppr = std::pow(2, dim) + std::pow(2, K.sOrthDir( currentSurfel ));
        unsigned int bitResu = 7 ^ bitSuppr;

        bitResu >> (sizeof(unsigned int) * CHAR_BIT - 1);

        unsigned int aaa = 0;
        unsigned int bbb = 1;
        unsigned int ccc = 2;
        unsigned int ddd = 4;

        aaa = aaa >> (sizeof(unsigned int) * CHAR_BIT - 1);
        bbb = bbb >> (sizeof(unsigned int) * CHAR_BIT - 1);
        ccc = ccc >> (sizeof(unsigned int) * CHAR_BIT - 1);
        ddd = ddd >> (sizeof(unsigned int) * CHAR_BIT - 1);

        std::cout << aaa << " " << bbb << " " << ccc << " " << ddd << std::endl;*/

        //          std::cout << "dim: " << dim << " --- orth: " << K.sOrthDir( currentSurfel ) << " --- result " << otherdim << std::endl;

        //dimSlice
        MySlice slice( &ptrTracker, otherdim );

        typedef typename MySlice::ConstIterator ConstIterator3D;
        typedef typename MySlice::ConstCirculator ConstCirculator3D;
        typedef SCellToMyPoint< typename MySlice::KSpace, Z2i::KSpace > Functor2D;
        typedef ConstIteratorAdapter< ConstIterator3D, Functor2D > ConstIterator3D2P;

        Functor2D pointFunctor(K, dim);

        ConstIterator3D2P pbegin( slice.begin(), pointFunctor );
        ConstIterator3D2P pend( slice.end(), pointFunctor );

        const Dimension size_slice = slice.size();
        std::vector< Statistic< double > > v_statMSEL( size_slice );
        for( Dimension ii = 0; ii < size_slice; ++ii )
        {
          v_statMSEL[ii] = Statistic<double>(true);
        }

        Circulator< ConstIterator3D2P > cbegin( pbegin, pbegin, pend );
        Circulator< ConstIterator3D2P > cend( cbegin );

        analyseAllLengthMS<Z2i::KSpace>( v_statMSEL, cbegin, cend );

        Dimension iii = 0;
        for( ConstIterator3D sit = slice.begin(), send = slice.end(); sit != send; ++sit )
        {
          Dimension surfel_pos = findSurfel( surfels, *sit );
          ASSERT( surfel_pos != surfels_size );
          if( marque[ surfels[surfel_pos] ] == false )
          {
            marque[ surfels[surfel_pos] ] = true;
            ASSERT(( marque.size() == surfels_size ));
            PairOfStatistics & otherpair = segments[ surfels[surfel_pos] ];
            ASSERT(( segments.size() == surfels_size ));
            if( otherpair.first.samples() == 0 )
            {
              otherpair.first = v_statMSEL[ iii ];
              //              segments[ (Surfel*)&(*sit) ] = otherpair;
            }
            else if ( otherpair.second.samples() == 0 )
            {
              otherpair.second = v_statMSEL[ iii ];
              //              segments[ (Surfel*)&(*sit) ] = otherpair;
            }
            else
            {
              FATAL_ERROR_MSG( false, "ALREADY FILLED" );
              trace.error() << "ALREADY FILLED" << std::endl;
            }
          }
          else
          {
            trace.error() << "WAHAAAAT?" << std::endl;
          }
          ++iii;
        }
      }
    }
  }

}

template< typename Surfel >
void computeGlobalSegment( std::vector< Surfel > & surfels,
                           std::map< Surfel, std::pair< Statistic< double >, Statistic< double > > > & segments,
                           Statistic< double > & allSegments,
                           const std::string & mode )
{
  const Dimension surfels_size = surfels.size();

  for( Dimension ii = 0; ii < surfels_size; ++ii )
  {
    Statistic< double > & stata = segments[ surfels[ii] ].first;
    Statistic< double > & statb = segments[ surfels[ii] ].second;

    if( mode == "max" )
    {
      if( stata.max() > statb.max() )
      {
        allSegments.addValues< typename Statistic<double>::Iterator >( stata.begin(), stata.end() );
      }
      else
      {
        allSegments.addValues< typename Statistic<double>::Iterator >( statb.begin(), statb.end() );
      }
    }
    else if( mode == "min" )
    {
      if( stata.max() < statb.max() )
      {
        allSegments.addValues< typename Statistic<double>::Iterator >( stata.begin(), stata.end() );
      }
      else
      {
        allSegments.addValues< typename Statistic<double>::Iterator >( statb.begin(), statb.end() );
      }
    }
    else if( mode == "mean" )
    {
      allSegments.addValues< typename Statistic<double>::Iterator >( stata.begin(), stata.end() );
      allSegments.addValues< typename Statistic<double>::Iterator >( statb.begin(), statb.end() );
    }
    else
    {
      trace.error() << "I dont understand " << mode << ". I need {min, mean, max} only." << std::endl;
    }
  }
}

void checkSizeRadius( double & re,
                      const double h,
                      const double globalRadius )
{
  double a = (1.0/3.0)*globalRadius;
  double b = 3.0*globalRadius;

  if( re < a )
  {
    re = a;
  }
  if( re > b )
  {
    re = b;
  }
}


struct Euclidean
{
  double distance(const double &a, const double &b) const
  {
    return (double)std::sqrt((b-a)*(b-a));
  }

  template< typename Point >
  double distance(const Point &a, const Point &b) const
  {
    return (double)(b-a).norm();
  }
};

void suggestedRadiusForIntegralInvariantEstimators( const std::vector< double > & radius,
                                                    std::vector< Dimension > & registration,
                                                    std::vector< double > & chosenRadius,
                                                    const Dimension nbRadius )
{
  KMeans<double, Euclidean>( radius, nbRadius, Euclidean(), registration, chosenRadius);
}

template< typename Surfel >
Dimension computeRadius( std::vector< Surfel > & surfels,
                         std::map< Surfel, std::pair< Statistic< double >, Statistic< double > > > & segments,
                         const double constante,
                         const double h,
                         std::vector< double > & radius,
                         const double globalRadius,
                         const std::string & prop,
                         const std::string & mode )
{
  const Dimension surfels_size = surfels.size();
  std::map< double, unsigned int > nbKernelRadius;

  ASSERT(( radius.size() == surfels_size ));
  ASSERT(( segments.size() == surfels_size ));

  for( Dimension ii = 0; ii < surfels_size; ++ii )
  {
    Statistic< double > & stata = segments[ surfels[ii] ].first;
    Statistic< double > & statb = segments[ surfels[ii] ].second;

    Statistic< double > stat(true);

    stata.terminate();
    statb.terminate();

    if( mode == "max" )
    {
      if( stata.max() > statb.max() )
      {
        stat = stata;
      }
      else
      {
        stat = statb;
      }
    }
    else if( mode == "min" )
    {
      if( stata.max() < statb.max() )
      {
        stat = stata;
      }
      else
      {
        stat = statb;
      }
    }
    else if( mode == "mean" )
    {
      stat.addValues< typename Statistic<double>::Iterator >( stata.begin(), stata.end() );
      stat.addValues< typename Statistic<double>::Iterator >( statb.begin(), statb.end() );
    }
    else
    {
      trace.error() << "I dont understand " << mode << ". I need {min, mean, max} only." << std::endl;
    }

    double result = -1.0;
    if( prop == "min" )
    {
      result = stat.min();
    }
    else if( prop == "mean" )
    {
      result = stat.mean();
    }
    else if( prop == "median" )
    {
      result = stat.median();
    }
    else if( prop == "max" )
    {
      result = stat.max();
    }

    ASSERT(( result != -1.0 ));

    double re = constante * result * result * h;

    checkSizeRadius( re, h, globalRadius );

    radius[ii] = re;

    if( nbKernelRadius.find(re) == nbKernelRadius.end() )
    {
      nbKernelRadius[ re ] = 1;
    }
    else
    {
      nbKernelRadius[ re ] += 1;
    }
  }

  return nbKernelRadius.size();
}

template< typename Estimator, typename Quantity, typename Functor, typename Surfel >
void computeCurvatureWithLocalSegments( const Z3i::KSpace & K,
                                        const Functor & functor,
                                        const std::vector< Surfel > & surfels,
                                        const std::vector< double > & radius,
                                        std::vector< Quantity > & curvatures )
{
  const Dimension surfels_size = surfels.size();

  ASSERT(( radius.size() == surfels_size ));
  ASSERT(( curvatures.size() == surfels_size ));

#ifdef WITH_OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
  for( Dimension ii = 0; ii < surfels_size; ++ii )
  {
    double re = radius[ ii ];
    re = ( re < 10.0 )? 10.0 : re;

    Estimator estimator ( K, functor );
    estimator.init( 1.0, radius[ ii ]);
    curvatures[ ii ] = estimator.eval( &surfels[ ii ] );
  }
}

template <typename Space, typename Shape>
bool
compareShapeEstimators( const std::string & filename,
                        const Shape * aShape,
                        const typename Space::RealPoint & border_min,
                        const typename Space::RealPoint & border_max,
                        const double & h,
                        //                        const double & radius_kernel,
                        //                        const double & alpha,
                        const std::string & options,
                        const std::string & properties,
                        const bool & lambda_optimized,
                        struct OptionsIntegralInvariant optionsII,
                        double noiseLevel = 0.0 )
{
  typedef typename Space::RealPoint RealPoint;
  typedef GaussDigitizer< Z3i::Space, Shape > DigitalShape;
  typedef Z3i::Domain Domain;
  typedef Z3i::KSpace KSpace;
  typedef typename KSpace::SCell SCell;
  typedef typename KSpace::Surfel Surfel;

  bool withNoise = ( noiseLevel <= 0.0 ) ? false : true;

  ASSERT (( noiseLevel < 1.0 ));
  // Digitizer
  DigitalShape dshape;
  dshape.attach( *aShape );
  dshape.init( border_min, border_max, h );
  Domain domain = dshape.getDomain();

  KSpace K;
  if ( ! K.init( domain.lowerBound(), domain.upperBound(), true ) )
  {
    std::cerr << "[3dLocalEstimators] error in creating KSpace." << std::endl;
    return false;
  }

  try
  {
    if ( withNoise )
    {
      /*
      typedef KanungoNoise< DigitalShape, Z3i::Domain > KanungoPredicate;
      typedef LightImplicitDigitalSurface< KSpace, KanungoPredicate > Boundary;
      typedef DigitalSurface< Boundary > MyDigitalSurface;
      typedef typename MyDigitalSurface::ConstIterator ConstIterator;

      typedef DepthFirstVisitor< MyDigitalSurface > Visitor;
      typedef GraphVisitorRange< Visitor > VisitorRange;
      typedef typename VisitorRange::ConstIterator VisitorConstIterator;

      typedef PointFunctorFromPointPredicateAndDomain< KanungoPredicate, Z3i::Domain, unsigned int > MyPointFunctor;
      typedef FunctorOnCells< MyPointFunctor, KSpace > MySpelFunctor;

      // Extracts shape boundary
      KanungoPredicate * noisifiedObject = new KanungoPredicate( *dshape, dshape->getDomain(), noiseLevel );
      SCell bel = Surfaces< KSpace >::findABel( K, *noisifiedObject, 10000 );
      Boundary * boundary = new Boundary( K, *noisifiedObject, SurfelAdjacency< KSpace::dimension >( true ), bel );
      MyDigitalSurface surf ( *boundary );

      double minsize = dshape->getUpperBound()[0] - dshape->getLowerBound()[0];
      unsigned int tries = 0;
      while( surf.size() < 2 * minsize || tries > 150 )
      {
        delete boundary;
        bel = Surfaces< KSpace >::findABel( K, *noisifiedObject, 10000 );
        boundary = new Boundary( K, *noisifiedObject, SurfelAdjacency< KSpace::dimension >( true ), bel );
        surf = MyDigitalSurface( *boundary );
        ++tries;
      }

      if( tries > 150 )
      {
        std::cerr << "Can't found a proper bel. So .... I ... just ... kill myself." << std::endl;
        return false;
      }

      VisitorRange * range;
      VisitorConstIterator ibegin;
      VisitorConstIterator iend;

      // Estimations
      Clock c;

      // True
      if( options.at( 0 ) != '0' )
      {
        // True Mean Curvature
        if( properties.at( 0 ) != '0' )
        {
          trace.beginBlock( "True mean curvature" );

          char full_filename[360];
          sprintf( full_filename, "%s%s", filename.c_str(), "_True_mean.dat" );
          std::ofstream file( full_filename );
          file << "# h = " << h << std::endl;
          file << "# True Mean Curvature estimation" << std::endl;

          std::ostream_iterator< double > out_it_true_mean( file, "\n" );

          range = new VisitorRange( new Visitor( surf, *surf.begin() ));
          ibegin = range->begin();
          iend = range->end();

          c.startClock();

          estimateTrueMeanCurvatureQuantity( ibegin,
                                             iend,
                                             out_it_true_mean,
                                             K,
                                             h,
                                             aShape );

          double TTrueMeanCurv = c.stopClock();
          file << "# time = " << TTrueMeanCurv << std::endl;

          file.close();
          delete range;

          trace.endBlock();
        }

        // True Gaussian Curvature
        if( properties.at( 1 ) != '0' )
        {
          trace.beginBlock( "True Gausian curvature" );

          char full_filename[360];
          sprintf( full_filename, "%s%s", filename.c_str(), "_True_gaussian.dat" );
          std::ofstream file( full_filename );
          file << "# h = " << h << std::endl;
          file << "# True Gaussian Curvature estimation" << std::endl;

          std::ostream_iterator< double > out_it_true_gaussian( file, "\n" );

          range = new VisitorRange( new Visitor( surf, *surf.begin() ));
          ibegin = range->begin();
          iend = range->end();

          c.startClock();

          estimateTrueGaussianCurvatureQuantity( ibegin,
                                                 iend,
                                                 out_it_true_gaussian,
                                                 K,
                                                 h,
                                                 aShape );

          double TTrueGaussianCurv = c.stopClock();
          file << "# time = " << TTrueGaussianCurv << std::endl;

          file.close();
          delete range;

          trace.endBlock();
        }

        // True Principal Curvatures
        if( properties.at( 2 ) != '0' )
        {
          trace.beginBlock( "True principal curvatures" );

          char full_filename[360];
          sprintf( full_filename, "%s%s", filename.c_str(), "_True_principal.dat" );
          std::ofstream file( full_filename );
          file << "# h = " << h << std::endl;
          file << "# True Gaussian Curvature estimation" << std::endl;

          std::ostream_iterator< CurvatureInformations > out_it_true_pc( file, "\n" );

          range = new VisitorRange( new Visitor( surf, *surf.begin() ));
          ibegin = range->begin();
          iend = range->end();

          c.startClock();

          estimateTruePrincipalCurvaturesQuantity( ibegin,
                                                   iend,
                                                   out_it_true_pc,
                                                   K,
                                                   h,
                                                   aShape );

          double TTruePrincCurv = c.stopClock();
          file << "# time = " << TTruePrincCurv << std::endl;

          file.close();

          delete range;

          trace.endBlock();
        }
      }

      double re_convolution_kernel = radius_kernel * std::pow( h, alpha ); // to obtains convergence results, re must follow the rule re=kh^(1/3)

      // II
      if( options.at( 1 ) != '0' )
      {
        MyPointFunctor * pointFunctor = new MyPointFunctor( noisifiedObject, dshape->getDomain(), 1, 0 );
        MySpelFunctor * functor = new MySpelFunctor( *pointFunctor, K );

        // Integral Invariant Mean Curvature
        if( properties.at( 0 ) != '0' )
        {
          trace.beginBlock( "II mean curvature" );

          IntegralInvariantMeanCurvatureEstimator< KSpace, MySpelFunctor > * IIMeanCurvatureEstimator = new IntegralInvariantMeanCurvatureEstimator< KSpace, MySpelFunctor >( K, *functor );

          c.startClock();
          IIMeanCurvatureEstimator->init ( h, re_convolution_kernel );

          char full_filename[360];
          sprintf( full_filename, "%s%s", filename.c_str(), "_II_mean.dat" );
          std::ofstream file( full_filename );
          file << "# h = " << h << std::endl;
          file << "# Mean Curvature estimation from the Integral Invariant" << std::endl;
          file << "# computed kernel radius = " << re_convolution_kernel << std::endl;

          range = new VisitorRange( new Visitor( surf, *surf.begin() ));
          ibegin = range->begin();
          iend = range->end();

          std::ostream_iterator< double > out_it_ii_mean( file, "\n" );
          IIMeanCurvatureEstimator->eval( ibegin, iend, out_it_ii_mean );

          double TIIMeanCurv = c.stopClock();
          file << "# time = " << TIIMeanCurv << std::endl;
          file.close();

          delete range;
          delete IIMeanCurvatureEstimator;

          trace.endBlock();
        }

        // Integral Invariant Gaussian Curvature
        if( properties.at( 1 ) != '0' )
        {
          trace.beginBlock( "II Gaussian curvature" );

          IntegralInvariantGaussianCurvatureEstimator< KSpace, MySpelFunctor > * IIGaussianCurvatureEstimator = new IntegralInvariantGaussianCurvatureEstimator< KSpace, MySpelFunctor >( K, *functor );

          c.startClock();
          IIGaussianCurvatureEstimator->init( h, re_convolution_kernel );

          char full_filename[360];
          sprintf( full_filename, "%s%s", filename.c_str(), "_II_gaussian.dat" );
          std::ofstream file( full_filename );
          file << "# h = " << h << std::endl;
          file << "# Gaussian Curvature estimation from the Integral Invariant" << std::endl;
          file << "# computed kernel radius = " << re_convolution_kernel << std::endl;

          std::ostream_iterator< double > out_it_ii_gaussian( file, "\n" );

          range = new VisitorRange( new Visitor( surf, *surf.begin() ));
          ibegin = range->begin();
          iend = range->end();

          IIGaussianCurvatureEstimator->eval ( ibegin, iend, out_it_ii_gaussian );

          double TIIGaussCurv = c.stopClock();
          file << "# time = " << TIIGaussCurv << std::endl;
          file.close();

          delete range;
          delete IIGaussianCurvatureEstimator;

          trace.endBlock();
        }

        // Integral Invariant Principal Curvatures
        if( properties.at( 2 ) != '0' )
        {
          trace.beginBlock( "II Principal curvatures" );

          IntegralInvariantGaussianCurvatureEstimator< KSpace, MySpelFunctor > * IIGaussianCurvatureEstimator = new IntegralInvariantGaussianCurvatureEstimator< KSpace, MySpelFunctor >( K, *functor );

          c.startClock();
          IIGaussianCurvatureEstimator->init( h, re_convolution_kernel );

          char full_filename[360];
          sprintf( full_filename, "%s%s", filename.c_str(), "_II_principal.dat" );
          std::ofstream file( full_filename );
          file << "# h = " << h << std::endl;
          file << "# Gaussian Curvature estimation from the Integral Invariant" << std::endl;
          file << "# computed kernel radius = " << re_convolution_kernel << std::endl;

          std::ostream_iterator< CurvatureInformations > out_it_ii_principal( file, "\n" );

          range = new VisitorRange( new Visitor( surf, *surf.begin() ));
          ibegin = range->begin();
          iend = range->end();

          IIGaussianCurvatureEstimator->evalPrincipalCurvatures ( ibegin, iend, out_it_ii_principal );

          double TIIGaussCurv = c.stopClock();
          file << "# time = " << TIIGaussCurv << std::endl;
          file.close();

          delete range;
          delete IIGaussianCurvatureEstimator;

          trace.endBlock();
        }
        //delete noisifiedObject;
        delete boundary;
        delete pointFunctor;
        delete functor;
      }

      // Monge
      if( options.at( 2 ) != '0' )
      {
        // Monge Mean Curvature
        if( properties.at( 0 ) != '0' )
        {
          trace.beginBlock( "Monge mean curvature" );

          typedef MongeJetFittingMeanCurvatureEstimator<Surfel, CanonicSCellEmbedder<KSpace> > FunctorMean;
          typedef ConstValueFunctor< double > ConvFunctor;
          typedef LocalEstimatorFromSurfelFunctorAdapter<typename MyDigitalSurface::DigitalSurfaceContainer, Z3i::L2Metric, FunctorMean, ConvFunctor> ReporterH;
          CanonicSCellEmbedder<KSpace> embedder( K );
          FunctorMean estimatorH( embedder, h );
          ConvFunctor convFunc(1.0);
          ReporterH reporterH( surf.container(), Z3i::l2Metric, estimatorH, convFunc);
          c.startClock();
          reporterH.init( h , re_convolution_kernel / h  );

          char full_filename[360];
          sprintf( full_filename, "%s%s", filename.c_str(), "_MongeJetFitting_mean.dat" );
          std::ofstream file( full_filename );
          file << "# h = " << h << std::endl;
          file << "# Mean Curvature estimation from CGAL Monge from and Jet Fitting" << std::endl;
          file << "# computed kernel radius = " << re_convolution_kernel << std::endl;
          std::ostream_iterator< double > out_it_monge_mean( file, "\n" );

          range = new VisitorRange( new Visitor( surf, *surf.begin() ));
          ibegin = range->begin();
          iend = range->end();
          //typename ReporterH::SurfelConstIterator aabegin = surf.begin();
          //typename ReporterH::SurfelConstIterator aaend = surf.end();
          reporterH.eval(ibegin, iend, out_it_monge_mean);
          double TMongeMeanCurv = c.stopClock();
          file << "# time = " << TMongeMeanCurv << std::endl;
          file.close();
          delete range;


          trace.endBlock();
        }

        // Monge Gaussian Curvature
        if( properties.at( 1 ) != '0' )
        {
          trace.beginBlock( "Monge Gaussian curvature" );

          typedef MongeJetFittingGaussianCurvatureEstimator<Surfel, CanonicSCellEmbedder<KSpace> > FunctorGaussian;
          typedef ConstValueFunctor< double > ConvFunctor;
          typedef LocalEstimatorFromSurfelFunctorAdapter<typename MyDigitalSurface::DigitalSurfaceContainer, Z3i::L2Metric, FunctorGaussian, ConvFunctor> ReporterK;
          CanonicSCellEmbedder<KSpace> embedder( K );
          FunctorGaussian estimatorK( embedder, h );
          ConvFunctor convFunc(1.0);
          ReporterK reporterK(surf.container(), Z3i::l2Metric, estimatorK, convFunc);
          c.startClock();
          reporterK.init( h , re_convolution_kernel / h  );

          //typename ReporterK::SurfelConstIterator aaabegin = surf.begin();
          //typename ReporterK::SurfelConstIterator aaaend = surf.end();

          range = new VisitorRange( new Visitor( surf, *surf.begin() ));
          ibegin = range->begin();
          iend = range->end();

          char full_filename[360];
          sprintf( full_filename, "%s%s", filename.c_str(), "_MongeJetFitting_gaussian.dat" );
          std::ofstream file( full_filename );
          file << "# h = " << h << std::endl;
          file << "# Gaussian Curvature estimation from CGAL Monge from and Jet Fitting" << std::endl;
          file << "# computed kernel radius = " << re_convolution_kernel << std::endl;
          std::ostream_iterator< double > out_it_monge_gaussian( file, "\n" );
          reporterK.eval(ibegin, iend , out_it_monge_gaussian);
          double TMongeGaussCurv = c.stopClock();
          file << "# time = " << TMongeGaussCurv << std::endl;
          file.close();
          delete range;


          trace.endBlock();
        }

        // Monge Principal Curvatures
        if( properties.at( 2 ) != '0' )
        {
          trace.beginBlock( "Monge Principal Curvature" );

          typedef MongeJetFittingPrincipalCurvaturesEstimator<Surfel, CanonicSCellEmbedder<KSpace> > FunctorPrincCurv;
          typedef ConstValueFunctor< double > ConvFunctor;
          typedef LocalEstimatorFromSurfelFunctorAdapter<typename MyDigitalSurface::DigitalSurfaceContainer, Z3i::L2Metric, FunctorPrincCurv, ConvFunctor> ReporterK;
          CanonicSCellEmbedder<KSpace> embedder( K );
          FunctorPrincCurv estimatorK( embedder, h );
          ConvFunctor convFunc(1.0);
          ReporterK reporterK(surf.container(), Z3i::l2Metric, estimatorK, convFunc);
          c.startClock();
          reporterK.init( h , re_convolution_kernel / h  );

          range = new VisitorRange( new Visitor( surf, *surf.begin() ));
          ibegin = range->begin();
          iend = range->end();

          char full_filename[360];
          sprintf( full_filename, "%s%s", filename.c_str(), "_MongeJetFitting_principal.dat" );
          std::ofstream file( full_filename );
          file << "# h = " << h << std::endl;
          file << "# Gaussian Curvature estimation from CGAL Monge from and Jet Fitting" << std::endl;
          file << "# computed kernel radius = " << re_convolution_kernel << std::endl;
          std::ostream_iterator< CurvatureInformations > out_it_monge_gaussian( file, "\n" );
          reporterK.eval(ibegin, iend , out_it_monge_gaussian);
          double TMongeGaussCurv = c.stopClock();
          file << "# time = " << TMongeGaussCurv << std::endl;
          file.close();
          delete range;


          trace.endBlock();
        }
      }
      */
    }
    else
    {
      typedef LightImplicitDigitalSurface< KSpace, DigitalShape > Boundary;
      typedef DigitalSurface< Boundary > MyDigitalSurface;
      typedef typename MyDigitalSurface::ConstIterator ConstIterator;

      typedef DepthFirstVisitor< MyDigitalSurface > Visitor;
      typedef GraphVisitorRange< Visitor > VisitorRange;
      typedef typename VisitorRange::ConstIterator VisitorConstIterator;

      typedef PointFunctorFromPointPredicateAndDomain< DigitalShape, Z3i::Domain, unsigned int > MyPointFunctor;
      typedef FunctorOnCells< MyPointFunctor, KSpace > MySpelFunctor;

      // Extracts shape boundary
      SCell bel = Surfaces<KSpace>::findABel ( K, dshape, 10000 );
      Boundary boundary( K, dshape, SurfelAdjacency< KSpace::dimension >( true ), bel );
      MyDigitalSurface surf ( boundary );

      VisitorRange * range;
      VisitorConstIterator ibegin;
      VisitorConstIterator iend;

      //      unsigned int cntIn = 0;
      //      typedef typename Domain::ConstIterator DomainConstIterator;
      //      DomainConstIterator d_ite = domain.end();
      //      for( DomainConstIterator it = domain.begin(); it != d_ite; ++it )
      //      {
      //        if( dshape(*it) )
      //        {
      //          ++cntIn;
      //        }
      //      }
      //
      //      std::cout << "boundary:" << surf.size() << std::endl;
      //      std::cout << "full:" << cntIn << std::endl;

      // Estimations
      Clock c;

      // True
      if( options.at( 0 ) != '0' )
      {
        // True Mean Curvature
        if( properties.at( 0 ) != '0' )
        {
          trace.beginBlock( "True mean curvature" );

          char full_filename[360];
          sprintf( full_filename, "%s%s", filename.c_str(), "_True_mean.dat" );
          std::ofstream file( full_filename );
          file << "# h = " << h << std::endl;
          file << "# True Mean Curvature estimation" << std::endl;

          std::ostream_iterator< double > out_it_true_mean( file, "\n" );

          range = new VisitorRange( new Visitor( surf, *surf.begin() ));
          ibegin = range->begin();
          iend = range->end();

          c.startClock();

          estimateTrueMeanCurvatureQuantity( ibegin,
                                             iend,
                                             out_it_true_mean,
                                             K,
                                             h,
                                             aShape );

          double TTrueMeanCurv = c.stopClock();
          file << "# time = " << TTrueMeanCurv << std::endl;

          file.close();
          //          delete range;

          trace.endBlock();
        }

        // True Gaussian Curvature
        if( properties.at( 1 ) != '0' )
        {
          trace.beginBlock( "True Gaussian curvature" );

          char full_filename[360];
          sprintf( full_filename, "%s%s", filename.c_str(), "_True_gaussian.dat" );
          std::ofstream file( full_filename );
          file << "# h = " << h << std::endl;
          file << "# True Gaussian Curvature estimation" << std::endl;

          std::ostream_iterator< double > out_it_true_gaussian( file, "\n" );

          range = new VisitorRange( new Visitor( surf, *surf.begin() ));
          ibegin = range->begin();
          iend = range->end();

          c.startClock();

          estimateTrueGaussianCurvatureQuantity( ibegin,
                                                 iend,
                                                 out_it_true_gaussian,
                                                 K,
                                                 h,
                                                 aShape );

          double TTrueGaussianCurv = c.stopClock();
          file << "# time = " << TTrueGaussianCurv << std::endl;

          file.close();

          //          delete range;

          trace.endBlock();
        }

        // True Principal Curvatures
        if( properties.at( 2 ) != '0' )
        {
          trace.beginBlock( "True principal curvatures" );

          char full_filename[360];
          sprintf( full_filename, "%s%s", filename.c_str(), "_True_principal.dat" );
          std::ofstream file( full_filename );
          file << "# h = " << h << std::endl;
          file << "# True Gaussian Curvature estimation" << std::endl;

          std::ostream_iterator< CurvatureInformations > out_it_true_pc( file, "\n" );

          range = new VisitorRange( new Visitor( surf, *surf.begin() ));
          ibegin = range->begin();
          iend = range->end();

          c.startClock();

          estimateTruePrincipalCurvaturesQuantity( ibegin,
                                                   iend,
                                                   out_it_true_pc,
                                                   K,
                                                   h,
                                                   aShape );

          double TTruePrincCurv = c.stopClock();
          file << "# time = " << TTruePrincCurv << std::endl;

          file.close();

          //          delete range;

          trace.endBlock();
        }
      }

      double re_convolution_kernel = 4.0;
      //      double re_convolution_kernel = radius_kernel * std::pow( h, alpha ); // to obtains convergence results, re must follow the rule re=kh^(1/3)

      // II
      if( options.at( 1 ) != '0' )
      {
        MyPointFunctor pointFunctor( dshape, domain, 1, 0 );
        MySpelFunctor functor( pointFunctor, K );

        double minRadiusAABB = ( border_max[0] - border_min[0] ) / 2.0; // this is a box.

        std::vector< Surfel > surfels;
        std::map< Surfel, std::pair< Statistic< double >, Statistic< double > > > segments;
        Statistic< double > allSegments( true );

        trace.beginBlock("Extracting all surfels...");

        {
          VisitorRange ranger( new Visitor( surf, *surf.begin() ) );
          VisitorConstIterator abegin = ranger.begin();
          VisitorConstIterator aend = ranger.end();

          while( abegin != aend )
          {
            surfels.push_back( *abegin );
            ++abegin;
          }
        }

        const Dimension size_surfels = surfels.size();

        trace.endBlock();

        trace.beginBlock("Analyse segments and Mapping segments <-> Surfels...");

        for( Dimension ii = 0; ii < size_surfels; ++ii )
        {
          segments[ surfels[ii] ] = std::pair< Statistic< double >, Statistic< double > >( Statistic< double >( true ), Statistic< double >( true ) );
        }

        if( optionsII.typeSegments == "normal")
        {
          computeSegments< Boundary, Surfel >( surfels, segments, K, boundary );
        }
        else if( optionsII.typeSegments == "DPS" )
        {
          computeDPSSegments< Boundary, Surfel >( surfels, segments, K, boundary, 1, 0.5*minRadiusAABB );
        }
        else
        {
          trace.error() << "Unknown option typeSegment " << optionsII.typeSegments << std::endl;
          trace.endBlock();
          return false;
        }
        ASSERT(( segments.size() == size_surfels ));

        trace.endBlock();

        // Global contour analysis
//        if( optionsII.tests.at(0) != '0' || optionsII.tests.at(1) != '0' )
        {
          computeGlobalSegment( surfels, segments, allSegments, optionsII.modeSegments );
          allSegments.terminate();
        }

        double constante = optionsII.constante;
        double g_mean = allSegments.mean();
        const double global_mean = constante * g_mean * g_mean * h;


        // Integral Invariant Mean Curvature
        if( properties.at( 0 ) != '0' )
        {
          /////// global max
          if( optionsII.tests.at(0) != '0' )
          {
            c.startClock();

            char full_filename[360];
            sprintf( full_filename, "%s%s%s%s", filename.c_str(), "_II_mean_curvature_global_max_segments_", optionsII.modeSegments.c_str() ,".dat" );
            std::ofstream file( full_filename );

            file << "# h = " << h << std::endl;
            file << "# Mean Curvature estimation from the Integral Invariant with global max segment analysis" << std::endl;

            double max = allSegments.max();
            double re = constante * max * max * h;
//            checkSizeRadius( re, h, minRadiusAABB );

            file << "# computed kernel radius = " << re << std::endl;

            trace.beginBlock("II mean curvature computation with global max segments...");

            typedef IntegralInvariantMeanCurvatureEstimator< KSpace, MySpelFunctor > Estimator;
            Estimator estimator( K, functor );
            estimator.init( h, re );

            VisitorRange ranger( new Visitor( surf, *surf.begin() ) );
            VisitorConstIterator abegin = ranger.begin();
            VisitorConstIterator aend = ranger.end();

            std::ostream_iterator< double > out( file, "\n" );

            estimator.eval( abegin, aend, out );

            trace.endBlock();
          }

          /////// global mean
          if( optionsII.tests.at(1) != '0' )
          {
            c.startClock();

            char full_filename[360];
            sprintf( full_filename, "%s%s%s%s", filename.c_str(), "_II_mean_curvature_global_mean_segments_", optionsII.modeSegments.c_str() ,".dat" );
            std::ofstream file( full_filename );

            file << "# h = " << h << std::endl;
            file << "# Mean Curvature estimation from the Integral Invariant with global mean segment analysis" << std::endl;


            double re = global_mean;
//            checkSizeRadius( re, h, minRadiusAABB );

            file << "# computed kernel radius = " << re << std::endl;

            trace.beginBlock("II mean curvature computation with global mean segments...");

            typedef IntegralInvariantMeanCurvatureEstimator< KSpace, MySpelFunctor > Estimator;
            Estimator estimator( K, functor );
            estimator.init( h, re );

            VisitorRange ranger( new Visitor( surf, *surf.begin() ) );
            VisitorConstIterator abegin = ranger.begin();
            VisitorConstIterator aend = ranger.end();

            std::ostream_iterator< double > out( file, "\n" );

            estimator.eval( abegin, aend, out );

            trace.endBlock();
          }

          /////// global median
          if( optionsII.tests.at(2) != '0' )
          {
            c.startClock();

            char full_filename[360];
            sprintf( full_filename, "%s%s%s%s", filename.c_str(), "_II_mean_curvature_global_median_segments_", optionsII.modeSegments.c_str() ,".dat" );
            std::ofstream file( full_filename );

            file << "# h = " << h << std::endl;
            file << "# Mean Curvature estimation from the Integral Invariant with global median segment analysis" << std::endl;

            double median = allSegments.median();
            double re = constante * median * median * h;
//            checkSizeRadius( re, h, minRadiusAABB );
            file << "# computed kernel radius = " << re << std::endl;

            trace.beginBlock("II mean curvature computation with global median segments...");

            typedef IntegralInvariantMeanCurvatureEstimator< KSpace, MySpelFunctor > Estimator;
            Estimator estimator( K, functor );
            estimator.init( h, re );

            VisitorRange ranger( new Visitor( surf, *surf.begin() ) );
            VisitorConstIterator abegin = ranger.begin();
            VisitorConstIterator aend = ranger.end();

            std::ostream_iterator< double > out( file, "\n" );

            estimator.eval( abegin, aend, out );

            trace.endBlock();
          }

          /////// global min
          if( optionsII.tests.at(3) != '0' )
          {
            c.startClock();

            char full_filename[360];
            sprintf( full_filename, "%s%s%s%s", filename.c_str(), "_II_mean_curvature_global_min_segments_", optionsII.modeSegments.c_str() ,".dat" );
            std::ofstream file( full_filename );

            file << "# h = " << h << std::endl;
            file << "# Mean Curvature estimation from the Integral Invariant with global min segment analysis" << std::endl;

            double min = allSegments.min();
            double re = constante * min * min * h;
//            checkSizeRadius( re, h, minRadiusAABB );
            file << "# computed kernel radius = " << re << std::endl;

            trace.beginBlock("II mean curvature computation with global min segments...");

            typedef IntegralInvariantMeanCurvatureEstimator< KSpace, MySpelFunctor > Estimator;
            Estimator estimator( K, functor );
            estimator.init( h, re );

            VisitorRange ranger( new Visitor( surf, *surf.begin() ) );
            VisitorConstIterator abegin = ranger.begin();
            VisitorConstIterator aend = ranger.end();

            std::ostream_iterator< double > out( file, "\n" );

            estimator.eval( abegin, aend, out );

            trace.endBlock();
          }

          /////// local max
          if( optionsII.tests.at(4) != '0' )
          {
            char full_filename[360];
            sprintf( full_filename, "%s%s%s%s%u%s", filename.c_str(), "_II_mean_curvature_local_max_segments_", optionsII.modeSegments.c_str(), "_with_", optionsII.nbKernels, "kernels.dat" );
            std::ofstream file( full_filename );

            file << "# h = " << h << std::endl;
            file << "# Mean Curvature estimation from the Integral Invariant with local max segment analysis" << std::endl;

            trace.beginBlock("Computation of radius...");

            std::vector< double > v_estimated_radius;
            v_estimated_radius.resize( size_surfels );

            Dimension nbKernelsRadius = computeRadius< Surfel >( surfels, segments, constante, h, v_estimated_radius, global_mean, "max", optionsII.modeSegments );
            if( nbKernelsRadius < optionsII.nbKernels )
            {
              optionsII.nbKernels = nbKernelsRadius;
            }

            ASSERT(( v_estimated_radius.size() == size_surfels ));
            trace.endBlock();

            typedef IntegralInvariantMeanCurvatureEstimator< Z3i::KSpace, MySpelFunctor > Estimator;

            trace.beginBlock("Sorting radius & pre-computing estimators...");

            std::vector< double > v_radius;
            std::vector< Dimension > v_registration;
            std::vector< Estimator* > v_estimators;

            if( optionsII.nbKernels > 0 )
            {
              v_estimators.resize( optionsII.nbKernels );
              suggestedRadiusForIntegralInvariantEstimators( v_estimated_radius, v_registration, v_radius, optionsII.nbKernels );
              ASSERT(( v_radius.size() == optionsII.nbKernels ));

#ifdef WITH_OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
              for( Dimension ii = 0; ii < optionsII.nbKernels; ++ii )
              {
                v_estimators[ii] = new Estimator( K, functor );
                v_estimators[ii]->init( h, v_radius[ii] );

                std::cout << "estimator #" << ii << " of radius " << v_radius[ii] << " initialized ..." << std::endl;
              }
            }

            trace.endBlock();


            trace.beginBlock("II mean curvature computation with local max segments...");


            typedef double Quantity;
            std::vector< Quantity > v_curvatures( size_surfels, Quantity(0) );

#ifdef WITH_OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
            for( Dimension ii = 0; ii < size_surfels; ++ii )
            {
              if( optionsII.nbKernels > 0 )
              {
                if( v_radius[ v_registration[ii]] == -42 )
                {
                  v_curvatures[ii] = -42;
                }
                else
                {
                  v_curvatures[ii] = v_estimators[ v_registration[ ii ]]->eval( surfels.begin() + ii );
                }
              }
              else
              {
                if( v_estimated_radius[ii] == -42 )
                {
                  v_curvatures[ii] = -42;
                }
                else
                {
                  Estimator estimator( K, functor );
                  estimator.init( h, v_estimated_radius[ii] );
                  v_curvatures[ii] = estimator.eval( surfels.begin() + ii );
                }
              }
            }

            trace.endBlock();

            trace.beginBlock("Exporting results...");

            for( unsigned int ii = 0; ii < size_surfels; ++ii )
            {
              file << v_curvatures[ii] << std::endl;
              //                            file2 << v_estimated_radius[ii] << std::endl;
              //                            file4 << v_radius[v_registration[ii]] << std::endl;
            }

            trace.endBlock();
          }

          /////// local mean
          if( optionsII.tests.at(5) != '0' )
          {
            char full_filename[360];
            sprintf( full_filename, "%s%s%s%s%u%s", filename.c_str(), "_II_mean_curvature_local_mean_segments_", optionsII.modeSegments.c_str(), "_with_", optionsII.nbKernels, "kernels.dat" );
            std::ofstream file( full_filename );

            file << "# h = " << h << std::endl;
            file << "# Mean Curvature estimation from the Integral Invariant with local mean segment analysis" << std::endl;

            trace.beginBlock("Computation of radius...");

            std::vector< double > v_estimated_radius;
            v_estimated_radius.resize( size_surfels );

            Dimension nbKernelsRadius = computeRadius< Surfel >( surfels, segments, constante, h, v_estimated_radius, global_mean, "mean", optionsII.modeSegments );
            if( nbKernelsRadius < optionsII.nbKernels )
            {
              optionsII.nbKernels = nbKernelsRadius;
            }

            ASSERT(( v_estimated_radius.size() == size_surfels ));
            trace.endBlock();

            typedef IntegralInvariantMeanCurvatureEstimator< Z3i::KSpace, MySpelFunctor > Estimator;

            trace.beginBlock("Sorting radius & pre-computing estimators...");

            std::vector< double > v_radius;
            std::vector< Dimension > v_registration;
            std::vector< Estimator* > v_estimators;

            if( optionsII.nbKernels > 0 )
            {
              v_estimators.resize( optionsII.nbKernels );
              suggestedRadiusForIntegralInvariantEstimators( v_estimated_radius, v_registration, v_radius, optionsII.nbKernels );
              ASSERT(( v_radius.size() == optionsII.nbKernels ));

#ifdef WITH_OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
              for( Dimension ii = 0; ii < optionsII.nbKernels; ++ii )
              {
                v_estimators[ii] = new Estimator( K, functor );
                v_estimators[ii]->init( h, v_radius[ii] );

                std::cout << "estimator #" << ii << " of radius " << v_radius[ii] << " initialized ..." << std::endl;
              }
            }

            trace.endBlock();


            trace.beginBlock("II mean curvature computation with local mean segments...");


            typedef double Quantity;
            std::vector< Quantity > v_curvatures( size_surfels, Quantity(0) );

#ifdef WITH_OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
            for( Dimension ii = 0; ii < size_surfels; ++ii )
            {
              if( optionsII.nbKernels > 0 )
              {
                if( v_radius[ v_registration[ii]] == -42 )
                {
                  v_curvatures[ii] = -42;
                }
                else
                {
                  v_curvatures[ii] = v_estimators[ v_registration[ ii ]]->eval( surfels.begin() + ii );
                }
              }
              else
              {
                if( v_estimated_radius[ii] == -42 )
                {
                  v_curvatures[ii] = -42;
                }
                else
                {
                  Estimator estimator( K, functor );
                  estimator.init( h, v_estimated_radius[ii] );
                  v_curvatures[ii] = estimator.eval( surfels.begin() + ii );
                }
              }
            }

            trace.endBlock();

            trace.beginBlock("Exporting results...");

            for( unsigned int ii = 0; ii < size_surfels; ++ii )
            {
              file << v_curvatures[ii] << std::endl;
              //                            file2 << v_estimated_radius[ii] << std::endl;
              //                            file4 << v_radius[v_registration[ii]] << std::endl;
            }

            trace.endBlock();
          }

          /////// local median
          if( optionsII.tests.at(6) != '0' )
          {
            char full_filename[360];
            sprintf( full_filename, "%s%s%s%s%u%s", filename.c_str(), "_II_mean_curvature_local_median_segments_", optionsII.modeSegments.c_str(), "_with_", optionsII.nbKernels, "kernels.dat" );
            std::ofstream file( full_filename );

            file << "# h = " << h << std::endl;
            file << "# Mean Curvature estimation from the Integral Invariant with local median segment analysis" << std::endl;

            trace.beginBlock("Computation of radius...");

            std::vector< double > v_estimated_radius;
            v_estimated_radius.resize( size_surfels );

            Dimension nbKernelsRadius = computeRadius< Surfel >( surfels, segments, constante, h, v_estimated_radius, global_mean, "median", optionsII.modeSegments );
            if( nbKernelsRadius < optionsII.nbKernels )
            {
              optionsII.nbKernels = nbKernelsRadius;
            }

            ASSERT(( v_estimated_radius.size() == size_surfels ));
            trace.endBlock();

            typedef IntegralInvariantMeanCurvatureEstimator< Z3i::KSpace, MySpelFunctor > Estimator;

            trace.beginBlock("Sorting radius & pre-computing estimators...");

            std::vector< double > v_radius;
            std::vector< Dimension > v_registration;
            std::vector< Estimator* > v_estimators;

            if( optionsII.nbKernels > 0 )
            {
              v_estimators.resize( optionsII.nbKernels );
              suggestedRadiusForIntegralInvariantEstimators( v_estimated_radius, v_registration, v_radius, optionsII.nbKernels );
              ASSERT(( v_radius.size() == optionsII.nbKernels ));

#ifdef WITH_OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
              for( Dimension ii = 0; ii < optionsII.nbKernels; ++ii )
              {
                v_estimators[ii] = new Estimator( K, functor );
                v_estimators[ii]->init( h, v_radius[ii] );

                std::cout << "estimator #" << ii << " of radius " << v_radius[ii] << " initialized ..." << std::endl;
              }
            }

            trace.endBlock();


            trace.beginBlock("II mean curvature computation with local median segments...");


            typedef double Quantity;
            std::vector< Quantity > v_curvatures( size_surfels, Quantity(0) );

#ifdef WITH_OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
            for( Dimension ii = 0; ii < size_surfels; ++ii )
            {
              if( optionsII.nbKernels > 0 )
              {
                if( v_radius[ v_registration[ii]] == -42 )
                {
                  v_curvatures[ii] = -42;
                }
                else
                {
                  v_curvatures[ii] = v_estimators[ v_registration[ ii ]]->eval( surfels.begin() + ii );
                }
              }
              else
              {
                if( v_estimated_radius[ii] == -42 )
                {
                  v_curvatures[ii] = -42;
                }
                else
                {
                  Estimator estimator( K, functor );
                  estimator.init( h, v_estimated_radius[ii] );
                  v_curvatures[ii] = estimator.eval( surfels.begin() + ii );
                }
              }
            }

            trace.endBlock();

            trace.beginBlock("Exporting results...");

            for( unsigned int ii = 0; ii < size_surfels; ++ii )
            {
              file << v_curvatures[ii] << std::endl;
              //                            file2 << v_estimated_radius[ii] << std::endl;
              //                            file4 << v_radius[v_registration[ii]] << std::endl;
            }

            trace.endBlock();
          }

          /////// local min
          if( optionsII.tests.at(7) != '0' )
          {
            char full_filename[360];
            sprintf( full_filename, "%s%s%s%s%u%s", filename.c_str(), "_II_mean_curvature_local_min_segments_", optionsII.modeSegments.c_str(), "_with_", optionsII.nbKernels, "kernels.dat" );
            std::ofstream file( full_filename );

            file << "# h = " << h << std::endl;
            file << "# Mean Curvature estimation from the Integral Invariant with local min segment analysis" << std::endl;

            trace.beginBlock("Computation of radius...");

            std::vector< double > v_estimated_radius;
            v_estimated_radius.resize( size_surfels );

            Dimension nbKernelsRadius = computeRadius< Surfel >( surfels, segments, constante, h, v_estimated_radius, global_mean, "min", optionsII.modeSegments );
            if( nbKernelsRadius < optionsII.nbKernels )
            {
              optionsII.nbKernels = nbKernelsRadius;
            }

            ASSERT(( v_estimated_radius.size() == size_surfels ));
            trace.endBlock();

            typedef IntegralInvariantMeanCurvatureEstimator< Z3i::KSpace, MySpelFunctor > Estimator;

            trace.beginBlock("Sorting radius & pre-computing estimators...");

            std::vector< double > v_radius;
            std::vector< Dimension > v_registration;
            std::vector< Estimator* > v_estimators;

            if( optionsII.nbKernels > 0 )
            {
              v_estimators.resize( optionsII.nbKernels );
              suggestedRadiusForIntegralInvariantEstimators( v_estimated_radius, v_registration, v_radius, optionsII.nbKernels );
              ASSERT(( v_radius.size() == optionsII.nbKernels ));

#ifdef WITH_OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
              for( Dimension ii = 0; ii < optionsII.nbKernels; ++ii )
              {
                v_estimators[ii] = new Estimator( K, functor );
                v_estimators[ii]->init( h, v_radius[ii] );

                std::cout << "estimator #" << ii << " of radius " << v_radius[ii] << " initialized ..." << std::endl;
              }
            }

            trace.endBlock();


            trace.beginBlock("II mean curvature computation with local min segments...");


            typedef double Quantity;
            std::vector< Quantity > v_curvatures( size_surfels, Quantity(0) );

#ifdef WITH_OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
            for( Dimension ii = 0; ii < size_surfels; ++ii )
            {
              if( optionsII.nbKernels > 0 )
              {
                if( v_radius[ v_registration[ii]] == -42 )
                {
                  v_curvatures[ii] = -42;
                }
                else
                {
                  v_curvatures[ii] = v_estimators[ v_registration[ ii ]]->eval( surfels.begin() + ii );
                }
              }
              else
              {
                if( v_estimated_radius[ii] == -42 )
                {
                  v_curvatures[ii] = -42;
                }
                else
                {
                  Estimator estimator( K, functor );
                  estimator.init( h, v_estimated_radius[ii] );
                  v_curvatures[ii] = estimator.eval( surfels.begin() + ii );
                }
              }
            }

            trace.endBlock();

            trace.beginBlock("Exporting results...");

            for( unsigned int ii = 0; ii < size_surfels; ++ii )
            {
              file << v_curvatures[ii] << std::endl;
              //                            file2 << v_estimated_radius[ii] << std::endl;
              //                            file4 << v_radius[v_registration[ii]] << std::endl;
            }

            trace.endBlock();
          }
        }

        // Integral Invariant Gaussian Curvature
        if( properties.at( 1 ) != '0' )
        {
          /////// global max
          if( optionsII.tests.at(0) != '0' )
          {
            c.startClock();

            char full_filename[360];
            sprintf( full_filename, "%s%s%s%s", filename.c_str(), "_II_gaussian_curvature_global_max_segments_", optionsII.modeSegments.c_str() ,".dat" );
            std::ofstream file( full_filename );

            file << "# h = " << h << std::endl;
            file << "# Gaussian Curvature estimation from the Integral Invariant with global max segment analysis" << std::endl;

            double max = allSegments.max();
            double re = constante * max * max * h;
//            checkSizeRadius( re, h, minRadiusAABB );

            file << "# computed kernel radius = " << re << std::endl;

            trace.beginBlock("II Gaussian curvature computation with global max segments...");

            typedef IntegralInvariantGaussianCurvatureEstimator< KSpace, MySpelFunctor > Estimator;
            Estimator estimator( K, functor );
            estimator.init( h, re );

            VisitorRange ranger( new Visitor( surf, *surf.begin() ) );
            VisitorConstIterator abegin = ranger.begin();
            VisitorConstIterator aend = ranger.end();

            std::ostream_iterator< double > out( file, "\n" );

            estimator.eval( abegin, aend, out );

            trace.endBlock();
          }

          /////// global mean
          if( optionsII.tests.at(1) != '0' )
          {
            c.startClock();

            char full_filename[360];
            sprintf( full_filename, "%s%s%s%s", filename.c_str(), "_II_gaussian_curvature_global_mean_segments_", optionsII.modeSegments.c_str() ,".dat" );
            std::ofstream file( full_filename );

            file << "# h = " << h << std::endl;
            file << "# Gaussian Curvature estimation from the Integral Invariant with global mean segment analysis" << std::endl;


            double mean = allSegments.mean();
            double re = constante * mean * mean * h;
//            checkSizeRadius( re, h, minRadiusAABB );

            file << "# computed kernel radius = " << re << std::endl;

            trace.beginBlock("II Gaussian curvature computation with global mean segments...");

            typedef IntegralInvariantGaussianCurvatureEstimator< KSpace, MySpelFunctor > Estimator;
            Estimator estimator( K, functor );
            estimator.init( h, re );

            VisitorRange ranger( new Visitor( surf, *surf.begin() ) );
            VisitorConstIterator abegin = ranger.begin();
            VisitorConstIterator aend = ranger.end();

            std::ostream_iterator< double > out( file, "\n" );

            estimator.eval( abegin, aend, out );

            trace.endBlock();
          }

          /////// global median
          if( optionsII.tests.at(2) != '0' )
          {
            c.startClock();

            char full_filename[360];
            sprintf( full_filename, "%s%s%s%s", filename.c_str(), "_II_gaussian_curvature_global_median_segments_", optionsII.modeSegments.c_str() ,".dat" );
            std::ofstream file( full_filename );

            file << "# h = " << h << std::endl;
            file << "# Gaussian Curvature estimation from the Integral Invariant with global median segment analysis" << std::endl;


            double median = allSegments.median();
            double re = constante * median * median * h;
//            checkSizeRadius( re, h, minRadiusAABB );
            file << "# computed kernel radius = " << re << std::endl;

            trace.beginBlock("II Gaussian curvature computation with global median segments...");

            typedef IntegralInvariantGaussianCurvatureEstimator< KSpace, MySpelFunctor > Estimator;
            Estimator estimator( K, functor );
            estimator.init( h, re );

            VisitorRange ranger( new Visitor( surf, *surf.begin() ) );
            VisitorConstIterator abegin = ranger.begin();
            VisitorConstIterator aend = ranger.end();

            std::ostream_iterator< double > out( file, "\n" );

            estimator.eval( abegin, aend, out );

            trace.endBlock();
          }

          /////// global min
          if( optionsII.tests.at(3) != '0' )
          {
            c.startClock();

            char full_filename[360];
            sprintf( full_filename, "%s%s%s%s", filename.c_str(), "_II_gaussian_curvature_global_min_segments_", optionsII.modeSegments.c_str() ,".dat" );
            std::ofstream file( full_filename );

            file << "# h = " << h << std::endl;
            file << "# Gaussian Curvature estimation from the Integral Invariant with global min segment analysis" << std::endl;


            double min = allSegments.min();
            double re = constante * min * min * h;
//            checkSizeRadius( re, h, minRadiusAABB );
            file << "# computed kernel radius = " << re << std::endl;

            trace.beginBlock("II Gaussian curvature computation with global min segments...");

            typedef IntegralInvariantGaussianCurvatureEstimator< KSpace, MySpelFunctor > Estimator;
            Estimator estimator( K, functor );
            estimator.init( h, re );

            VisitorRange ranger( new Visitor( surf, *surf.begin() ) );
            VisitorConstIterator abegin = ranger.begin();
            VisitorConstIterator aend = ranger.end();

            std::ostream_iterator< double > out( file, "\n" );

            estimator.eval( abegin, aend, out );

            trace.endBlock();
          }

          /////// local max
          if( optionsII.tests.at(4) != '0' )
          {
            char full_filename[360];
            sprintf( full_filename, "%s%s%s%s%u%s", filename.c_str(), "_II_mean_curvature_local_max_segments_", optionsII.modeSegments.c_str(), "_with_", optionsII.nbKernels, "kernels.dat" );
            std::ofstream file( full_filename );

            file << "# h = " << h << std::endl;
            file << "# Mean Curvature estimation from the Integral Invariant with local max segment analysis" << std::endl;

            trace.beginBlock("Computation of radius...");

            std::vector< double > v_estimated_radius;
            v_estimated_radius.resize( size_surfels );
            Dimension nbKernelsRadius = computeRadius< Surfel >( surfels, segments, constante, h, v_estimated_radius, global_mean, "max", optionsII.modeSegments );
            if( nbKernelsRadius < optionsII.nbKernels )
            {
              optionsII.nbKernels = nbKernelsRadius;
            }

            ASSERT(( v_estimated_radius.size() == size_surfels ));
            trace.endBlock();

            typedef IntegralInvariantMeanCurvatureEstimator< Z3i::KSpace, MySpelFunctor > Estimator;

            trace.beginBlock("Sorting radius & pre-computing estimators...");

            std::vector< double > v_radius;
            std::vector< Dimension > v_registration;
            std::vector< Estimator* > v_estimators;

            if( optionsII.nbKernels > 0 )
            {
              v_estimators.resize( optionsII.nbKernels );
              suggestedRadiusForIntegralInvariantEstimators( v_estimated_radius, v_registration, v_radius, optionsII.nbKernels );
              ASSERT(( v_radius.size() == optionsII.nbKernels ));

#ifdef WITH_OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
              for( Dimension ii = 0; ii < optionsII.nbKernels; ++ii )
              {
                v_estimators[ii] = new Estimator( K, functor );
                v_estimators[ii]->init( h, v_radius[ii] );

                std::cout << "estimator #" << ii << " of radius " << v_radius[ii] << " initialized ..." << std::endl;
              }
            }

            trace.endBlock();


            trace.beginBlock("II mean curvature computation with local max segments...");


            typedef double Quantity;
            std::vector< Quantity > v_curvatures( size_surfels, Quantity(0) );

#ifdef WITH_OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
            for( Dimension ii = 0; ii < size_surfels; ++ii )
            {
              if( optionsII.nbKernels > 0 )
              {
                v_curvatures[ii] = v_estimators[ v_registration[ ii ]]->eval( surfels.begin() + ii );
              }
              else
              {
                Estimator estimator( K, functor );
                estimator.init( h, v_estimated_radius[ii] );
                v_curvatures[ii] = estimator.eval( surfels.begin() + ii );
              }
            }

            trace.endBlock();

            trace.beginBlock("Exporting results...");

            for( unsigned int ii = 0; ii < size_surfels; ++ii )
            {
              file << v_curvatures[ii] << std::endl;
              //                            file2 << v_estimated_radius[ii] << std::endl;
              //                            file4 << v_radius[v_registration[ii]] << std::endl;
            }

            trace.endBlock();
          }

          /////// local mean
          if( optionsII.tests.at(5) != '0' )
          {
            char full_filename[360];
            sprintf( full_filename, "%s%s%s%s", filename.c_str(), "_II_gaussian_curvature_local_mean_segments_", optionsII.modeSegments.c_str() ,".dat" );
            std::ofstream file( full_filename );

            file << "# h = " << h << std::endl;
            file << "# Gaussian Curvature estimation from the Integral Invariant with local mean segment analysis" << std::endl;

            trace.beginBlock("Computation of radius...");

            std::vector< double > radius;
            radius.resize( size_surfels );
            computeRadius< Surfel >( surfels, segments, constante, h, radius, global_mean, "mean", optionsII.modeSegments );

            trace.endBlock();

            typedef IntegralInvariantGaussianCurvatureEstimator< Z3i::KSpace, MySpelFunctor > Estimator;

            ASSERT(( radius.size() == size_surfels ));

            trace.beginBlock("II Gaussian curvature computation with local mean segments...");


#ifdef WITH_OPENMP
            typedef double Quantity;
            std::vector< Quantity > curvatures( size_surfels, Quantity(0) );

#pragma omp parallel for schedule(dynamic)
            for( Dimension ii = 0; ii < size_surfels; ++ii )
            {
              Estimator estimator ( K, functor );
              estimator.init( h, radius[ ii ]);
              curvatures[ii] = estimator.eval( &surfels[ ii ] );
            }

            for( Dimension ii = 0; ii < size_surfels; ++ii )
            {
              file << curvatures[ii] << std::endl;
            }
#else
            for( Dimension ii = 0; ii < size_surfels; ++ii )
            {
              double re = radius[ ii ];

              Estimator estimator ( K, functor );
              estimator.init( h, radius[ ii ]);
              file << estimator.eval( &surfels[ ii ] ) << std::endl;
            }
#endif

            trace.endBlock();
          }

          /////// local median
          if( optionsII.tests.at(6) != '0' )
          {
            char full_filename[360];
            sprintf( full_filename, "%s%s%s%s", filename.c_str(), "_II_gaussian_curvature_local_median_segments_", optionsII.modeSegments.c_str() ,".dat" );
            std::ofstream file( full_filename );

            file << "# h = " << h << std::endl;
            file << "# Gaussian Curvature estimation from the Integral Invariant with local median segment analysis" << std::endl;

            trace.beginBlock("Computation of radius...");

            std::vector< double > radius;
            radius.resize( size_surfels );
            computeRadius< Surfel >( surfels, segments, constante, h, radius, global_mean, "median", optionsII.modeSegments );

            trace.endBlock();

            typedef IntegralInvariantGaussianCurvatureEstimator< Z3i::KSpace, MySpelFunctor > Estimator;

            ASSERT(( radius.size() == size_surfels ));

            trace.beginBlock("II Gaussian curvature computation with local median segments...");


#ifdef WITH_OPENMP
            typedef double Quantity;
            std::vector< Quantity > curvatures( size_surfels, Quantity(0) );

#pragma omp parallel for schedule(dynamic)
            for( Dimension ii = 0; ii < size_surfels; ++ii )
            {
              Estimator estimator ( K, functor );
              estimator.init( h, radius[ ii ]);
              curvatures[ii] = estimator.eval( &surfels[ ii ] );
            }

            for( Dimension ii = 0; ii < size_surfels; ++ii )
            {
              file << curvatures[ii] << std::endl;
            }
#else
            for( Dimension ii = 0; ii < size_surfels; ++ii )
            {
              double re = radius[ ii ];

              Estimator estimator ( K, functor );
              estimator.init( h, radius[ ii ]);
              file << estimator.eval( &surfels[ ii ] ) << std::endl;
            }
#endif

            trace.endBlock();
          }

          /////// local min
          if( optionsII.tests.at(7) != '0' )
          {
            char full_filename[360];
            sprintf( full_filename, "%s%s%s%s", filename.c_str(), "_II_gaussian_curvature_local_min_segments_", optionsII.modeSegments.c_str() ,".dat" );
            std::ofstream file( full_filename );

            file << "# h = " << h << std::endl;
            file << "# Gaussian Curvature estimation from the Integral Invariant with local min segment analysis" << std::endl;

            trace.beginBlock("Computation of radius...");

            std::vector< double > radius;
            radius.resize( size_surfels );
            computeRadius< Surfel >( surfels, segments, constante, h, radius, global_mean, "min", optionsII.modeSegments );

            trace.endBlock();

            typedef IntegralInvariantGaussianCurvatureEstimator< Z3i::KSpace, MySpelFunctor > Estimator;

            ASSERT(( radius.size() == size_surfels ));

            trace.beginBlock("II Gaussian curvature computation with local min segments...");


#ifdef WITH_OPENMP
            typedef double Quantity;
            std::vector< Quantity > curvatures( size_surfels, Quantity(0) );

#pragma omp parallel for schedule(dynamic)
            for( Dimension ii = 0; ii < size_surfels; ++ii )
            {
              Estimator estimator ( K, functor );
              estimator.init( h, radius[ ii ]);
              curvatures[ii] = estimator.eval( &surfels[ ii ] );
            }

            for( Dimension ii = 0; ii < size_surfels; ++ii )
            {
              file << curvatures[ii] << std::endl;
            }
#else
            for( Dimension ii = 0; ii < size_surfels; ++ii )
            {
              double re = radius[ ii ];

              Estimator estimator ( K, functor );
              estimator.init( h, radius[ ii ]);
              file << estimator.eval( &surfels[ ii ] ) << std::endl;
            }
#endif

            trace.endBlock();
          }
        }

        // Integral Invariant Principal Curvatures
        if( properties.at( 2 ) != '0' )
        {
          /////// global max
          if( optionsII.tests.at(0) != '0' )
          {
            c.startClock();

            char full_filename[360];
            sprintf( full_filename, "%s%s%s%s", filename.c_str(), "_II_principal_curvatures_global_max_segments_", optionsII.modeSegments.c_str() ,".dat" );
            std::ofstream file( full_filename );

            file << "# h = " << h << std::endl;
            file << "# Principal curvatures estimation from the Integral Invariant with global max segment analysis" << std::endl;

            double max = allSegments.max();
            double re = constante * max * max * h;
//            checkSizeRadius( re, h, minRadiusAABB );
            file << "# computed kernel radius = " << re << std::endl;

            trace.beginBlock("II Principal curvatures computation with global max segments...");

            typedef IntegralInvariantGaussianCurvatureEstimator< KSpace, MySpelFunctor > Estimator;
            Estimator estimator( K, functor );
            estimator.init( h, re );

            VisitorRange ranger( new Visitor( surf, *surf.begin() ) );
            VisitorConstIterator abegin = ranger.begin();
            VisitorConstIterator aend = ranger.end();

            std::ostream_iterator< CurvatureInformations > out( file, "\n" );

            estimator.evalPrincipalCurvatures( abegin, aend, out );

            trace.endBlock();
          }

          /////// global mean
          if( optionsII.tests.at(1) != '0' )
          {
            c.startClock();

            char full_filename[360];
            sprintf( full_filename, "%s%s%s%s", filename.c_str(), "_II_principal_curvatures_global_mean_segments_", optionsII.modeSegments.c_str() ,".dat" );
            std::ofstream file( full_filename );

            file << "# h = " << h << std::endl;
            file << "# Principal curvatures estimation from the Integral Invariant with global mean segment analysis" << std::endl;

            double mean = allSegments.mean();
            double re = constante * mean * mean * h;
//            checkSizeRadius( re, h, minRadiusAABB );
            file << "# computed kernel radius = " << re << std::endl;

            trace.beginBlock("II Principal curvatures computation with global mean segments...");

            typedef IntegralInvariantGaussianCurvatureEstimator< KSpace, MySpelFunctor > Estimator;
            Estimator estimator( K, functor );
            estimator.init( h, re );

            VisitorRange ranger( new Visitor( surf, *surf.begin() ) );
            VisitorConstIterator abegin = ranger.begin();
            VisitorConstIterator aend = ranger.end();

            std::ostream_iterator< CurvatureInformations > out( file, "\n" );

            estimator.evalPrincipalCurvatures( abegin, aend, out );

            trace.endBlock();
          }

          /////// global median
          if( optionsII.tests.at(2) != '0' )
          {
            c.startClock();

            char full_filename[360];
            sprintf( full_filename, "%s%s%s%s", filename.c_str(), "_II_principal_curvatures_global_median_segments_", optionsII.modeSegments.c_str() ,".dat" );
            std::ofstream file( full_filename );

            file << "# h = " << h << std::endl;
            file << "# Principal curvatures estimation from the Integral Invariant with global median segment analysis" << std::endl;

            double median = allSegments.median();
            double re = constante * median * median * h;
//            checkSizeRadius( re, h, minRadiusAABB );
            file << "# computed kernel radius = " << re << std::endl;

            trace.beginBlock("II Principal curvatures computation with global median segments...");

            typedef IntegralInvariantGaussianCurvatureEstimator< KSpace, MySpelFunctor > Estimator;
            Estimator estimator( K, functor );
            estimator.init( h, re );

            VisitorRange ranger( new Visitor( surf, *surf.begin() ) );
            VisitorConstIterator abegin = ranger.begin();
            VisitorConstIterator aend = ranger.end();

            std::ostream_iterator< CurvatureInformations > out( file, "\n" );

            estimator.evalPrincipalCurvatures( abegin, aend, out );

            trace.endBlock();
          }

          /////// global min
          if( optionsII.tests.at(3) != '0' )
          {
            c.startClock();

            char full_filename[360];
            sprintf( full_filename, "%s%s%s%s", filename.c_str(), "_II_principal_curvatures_global_min_segments_", optionsII.modeSegments.c_str() ,".dat" );
            std::ofstream file( full_filename );

            file << "# h = " << h << std::endl;
            file << "# Principal curvatures estimation from the Integral Invariant with global min segment analysis" << std::endl;

            double min = allSegments.min();
            double re = constante * min * min * h;
//            checkSizeRadius( re, h, minRadiusAABB );
            file << "# computed kernel radius = " << re << std::endl;

            trace.beginBlock("II Principal curvatures computation with global min segments...");

            typedef IntegralInvariantGaussianCurvatureEstimator< KSpace, MySpelFunctor > Estimator;
            Estimator estimator( K, functor );
            estimator.init( h, re );

            VisitorRange ranger( new Visitor( surf, *surf.begin() ) );
            VisitorConstIterator abegin = ranger.begin();
            VisitorConstIterator aend = ranger.end();

            std::ostream_iterator< CurvatureInformations > out( file, "\n" );

            estimator.evalPrincipalCurvatures( abegin, aend, out );

            trace.endBlock();
          }

          /////// local max
          if( optionsII.tests.at(4) != '0' )
          {
            char full_filename[360];
            sprintf( full_filename, "%s%s%s%s%u%s", filename.c_str(), "_II_principal_curvatures_local_max_segments_", optionsII.modeSegments.c_str(), "_with_", optionsII.nbKernels, "kernels.dat" );
            std::ofstream file( full_filename );

            file << "# h = " << h << std::endl;
            file << "# Principal curvatures estimation from the Integral Invariant with local max segment analysis" << std::endl;

            trace.beginBlock("Computation of radius...");

            std::vector< double > v_estimated_radius;
            v_estimated_radius.resize( size_surfels );

            Dimension nbKernelsRadius = computeRadius< Surfel >( surfels, segments, constante, h, v_estimated_radius, global_mean, "max", optionsII.modeSegments );
            if( nbKernelsRadius < optionsII.nbKernels )
            {
              optionsII.nbKernels = nbKernelsRadius;
            }

            ASSERT(( v_estimated_radius.size() == size_surfels ));

            trace.endBlock();

            typedef IntegralInvariantGaussianCurvatureEstimator< Z3i::KSpace, MySpelFunctor > Estimator;

            trace.beginBlock("Sorting radius & pre-computing estimators...");

            std::vector< double > v_radius;
            std::vector< Dimension > v_registration;
            std::vector< Estimator* > v_estimators;

            if( optionsII.nbKernels > 0 )
            {
              v_estimators.resize( optionsII.nbKernels );
              suggestedRadiusForIntegralInvariantEstimators( v_estimated_radius, v_registration, v_radius, optionsII.nbKernels );
              ASSERT(( v_radius.size() == optionsII.nbKernels ));

#ifdef WITH_OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
              for( Dimension ii = 0; ii < optionsII.nbKernels; ++ii )
              {
                v_estimators[ii] = new Estimator( K, functor );
                v_estimators[ii]->init( h, v_radius[ii] );

                std::cout << "estimator #" << ii << " of radius " << v_radius[ii] << " initialized ..." << std::endl;
              }
            }

            trace.endBlock();

            trace.beginBlock("II Principal curvatures computation with local max segments...");

            typedef CurvatureInformations Quantity;
            std::vector< Quantity > v_curvatures( size_surfels, Quantity() );

#ifdef WITH_OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
            for( Dimension ii = 0; ii < size_surfels; ++ii )
            {
              if( optionsII.nbKernels > 0 )
              {
                v_curvatures[ii] = v_estimators[ v_registration[ ii ]]->evalPrincipalCurvatures( surfels.begin() + ii );
              }
              else
              {

                  Estimator estimator( K, functor );
                  estimator.init( h, v_estimated_radius[ii] );
                  v_curvatures[ii] = estimator.evalPrincipalCurvatures( surfels.begin() + ii );

              }
            }

            trace.endBlock();

            trace.beginBlock("Exporting results...");

            for( unsigned int ii = 0; ii < size_surfels; ++ii )
            {
              file << v_curvatures[ii] << std::endl;
              //                            file2 << v_estimated_radius[ii] << std::endl;
              //                            file4 << v_radius[v_registration[ii]] << std::endl;
            }

            trace.endBlock();
          }

          /////// local mean
          if( optionsII.tests.at(5) != '0' )
          {
            char full_filename[360];
            sprintf( full_filename, "%s%s%s%s%u%s", filename.c_str(), "_II_principal_curvatures_local_mean_segments_", optionsII.modeSegments.c_str(), "_with_", optionsII.nbKernels, "kernels.dat" );
            std::ofstream file( full_filename );

            file << "# h = " << h << std::endl;
            file << "# Principal curvatures estimation from the Integral Invariant with local mean segment analysis" << std::endl;

            trace.beginBlock("Computation of radius...");

            std::vector< double > v_estimated_radius;
            v_estimated_radius.resize( size_surfels );

            Dimension nbKernelsRadius = computeRadius< Surfel >( surfels, segments, constante, h, v_estimated_radius, global_mean, "mean", optionsII.modeSegments );
            if( nbKernelsRadius < optionsII.nbKernels )
            {
              optionsII.nbKernels = nbKernelsRadius;
            }

            ASSERT(( v_estimated_radius.size() == size_surfels ));

            trace.endBlock();

            typedef IntegralInvariantGaussianCurvatureEstimator< Z3i::KSpace, MySpelFunctor > Estimator;

            trace.beginBlock("Sorting radius & pre-computing estimators...");

            std::vector< double > v_radius;
            std::vector< Dimension > v_registration;
            std::vector< Estimator* > v_estimators;

            if( optionsII.nbKernels > 0 )
            {
              v_estimators.resize( optionsII.nbKernels );
              suggestedRadiusForIntegralInvariantEstimators( v_estimated_radius, v_registration, v_radius, optionsII.nbKernels );
              ASSERT(( v_radius.size() == optionsII.nbKernels ));

#ifdef WITH_OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
              for( Dimension ii = 0; ii < optionsII.nbKernels; ++ii )
              {
                v_estimators[ii] = new Estimator( K, functor );
                v_estimators[ii]->init( h, v_radius[ii] );

                std::cout << "estimator #" << ii << " of radius " << v_radius[ii] << " initialized ..." << std::endl;
              }
            }

            trace.endBlock();

            trace.beginBlock("II Principal curvatures computation with local mean segments...");

            typedef CurvatureInformations Quantity;
            std::vector< Quantity > v_curvatures( size_surfels, Quantity() );

#ifdef WITH_OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
            for( Dimension ii = 0; ii < size_surfels; ++ii )
            {
              if( optionsII.nbKernels > 0 )
              {
                v_curvatures[ii] = v_estimators[ v_registration[ ii ]]->evalPrincipalCurvatures( surfels.begin() + ii );
              }
              else
              {

                  Estimator estimator( K, functor );
                  estimator.init( h, v_estimated_radius[ii] );
                  v_curvatures[ii] = estimator.evalPrincipalCurvatures( surfels.begin() + ii );

              }
            }

            trace.endBlock();

            trace.beginBlock("Exporting results...");

            for( unsigned int ii = 0; ii < size_surfels; ++ii )
            {
              file << v_curvatures[ii] << std::endl;
              //                            file2 << v_estimated_radius[ii] << std::endl;
              //                            file4 << v_radius[v_registration[ii]] << std::endl;
            }

            trace.endBlock();
          }

          /////// local median
          if( optionsII.tests.at(6) != '0' )
          {
            char full_filename[360];
            sprintf( full_filename, "%s%s%s%s%u%s", filename.c_str(), "_II_principal_curvatures_local_median_segments_", optionsII.modeSegments.c_str(), "_with_", optionsII.nbKernels, "kernels.dat" );
            std::ofstream file( full_filename );

            file << "# h = " << h << std::endl;
            file << "# Principal curvatures estimation from the Integral Invariant with local median segment analysis" << std::endl;

            trace.beginBlock("Computation of radius...");

            std::vector< double > v_estimated_radius;
            v_estimated_radius.resize( size_surfels );

            Dimension nbKernelsRadius = computeRadius< Surfel >( surfels, segments, constante, h, v_estimated_radius, global_mean, "median", optionsII.modeSegments );
            if( nbKernelsRadius < optionsII.nbKernels )
            {
              optionsII.nbKernels = nbKernelsRadius;
            }

            ASSERT(( v_estimated_radius.size() == size_surfels ));

            trace.endBlock();

            typedef IntegralInvariantGaussianCurvatureEstimator< Z3i::KSpace, MySpelFunctor > Estimator;

            trace.beginBlock("Sorting radius & pre-computing estimators...");

            std::vector< double > v_radius;
            std::vector< Dimension > v_registration;
            std::vector< Estimator* > v_estimators;

            if( optionsII.nbKernels > 0 )
            {
              v_estimators.resize( optionsII.nbKernels );
              suggestedRadiusForIntegralInvariantEstimators( v_estimated_radius, v_registration, v_radius, optionsII.nbKernels );
              ASSERT(( v_radius.size() == optionsII.nbKernels ));

#ifdef WITH_OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
              for( Dimension ii = 0; ii < optionsII.nbKernels; ++ii )
              {
                v_estimators[ii] = new Estimator( K, functor );
                v_estimators[ii]->init( h, v_radius[ii] );

                std::cout << "estimator #" << ii << " of radius " << v_radius[ii] << " initialized ..." << std::endl;
              }
            }

            trace.endBlock();

            trace.beginBlock("II Principal curvatures computation with local median segments...");

            typedef CurvatureInformations Quantity;
            std::vector< Quantity > v_curvatures( size_surfels, Quantity() );

#ifdef WITH_OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
            for( Dimension ii = 0; ii < size_surfels; ++ii )
            {
              if( optionsII.nbKernels > 0 )
              {
                v_curvatures[ii] = v_estimators[ v_registration[ ii ]]->evalPrincipalCurvatures( surfels.begin() + ii );

              }
              else
              {

                  Estimator estimator( K, functor );
                  estimator.init( h, v_estimated_radius[ii] );
                  v_curvatures[ii] = estimator.evalPrincipalCurvatures( surfels.begin() + ii );

              }
            }

            trace.endBlock();

            trace.beginBlock("Exporting results...");

            for( unsigned int ii = 0; ii < size_surfels; ++ii )
            {
              file << v_curvatures[ii] << std::endl;
              //                            file2 << v_estimated_radius[ii] << std::endl;
              //                            file4 << v_radius[v_registration[ii]] << std::endl;
            }

            trace.endBlock();
          }

          /////// local min
          if( optionsII.tests.at(7) != '0' )
          {
            char full_filename[360];
            sprintf( full_filename, "%s%s%s%s%u%s", filename.c_str(), "_II_principal_curvatures_local_min_segments_", optionsII.modeSegments.c_str(), "_with_", optionsII.nbKernels, "kernels.dat" );
            std::ofstream file( full_filename );

            file << "# h = " << h << std::endl;
            file << "# Principal curvatures estimation from the Integral Invariant with local min segment analysis" << std::endl;

            trace.beginBlock("Computation of radius...");

            std::vector< double > v_estimated_radius;
            v_estimated_radius.resize( size_surfels );

            Dimension nbKernelsRadius = computeRadius< Surfel >( surfels, segments, constante, h, v_estimated_radius, global_mean, "min", optionsII.modeSegments );
            if( nbKernelsRadius < optionsII.nbKernels )
            {
              optionsII.nbKernels = nbKernelsRadius;
            }

            ASSERT(( v_estimated_radius.size() == size_surfels ));

            trace.endBlock();

            typedef IntegralInvariantGaussianCurvatureEstimator< Z3i::KSpace, MySpelFunctor > Estimator;

            trace.beginBlock("Sorting radius & pre-computing estimators...");

            std::vector< double > v_radius;
            std::vector< Dimension > v_registration;
            std::vector< Estimator* > v_estimators;

            if( optionsII.nbKernels > 0 )
            {
              v_estimators.resize( optionsII.nbKernels );
              suggestedRadiusForIntegralInvariantEstimators( v_estimated_radius, v_registration, v_radius, optionsII.nbKernels );
              ASSERT(( v_radius.size() == optionsII.nbKernels ));

#ifdef WITH_OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
              for( Dimension ii = 0; ii < optionsII.nbKernels; ++ii )
              {
                v_estimators[ii] = new Estimator( K, functor );
                v_estimators[ii]->init( h, v_radius[ii] );

                std::cout << "estimator #" << ii << " of radius " << v_radius[ii] << " initialized ..." << std::endl;
              }
            }

            trace.endBlock();

            trace.beginBlock("II Principal curvatures computation with local min segments...");

            typedef CurvatureInformations Quantity;
            std::vector< Quantity > v_curvatures( size_surfels, Quantity() );

#ifdef WITH_OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
            for( Dimension ii = 0; ii < size_surfels; ++ii )
            {
              if( optionsII.nbKernels > 0 )
              {
                v_curvatures[ii] = v_estimators[ v_registration[ ii ]]->evalPrincipalCurvatures( surfels.begin() + ii );
              }
              else
              {

                  Estimator estimator( K, functor );
                  estimator.init( h, v_estimated_radius[ii] );
                  v_curvatures[ii] = estimator.evalPrincipalCurvatures( surfels.begin() + ii );

              }
            }

            trace.endBlock();

            trace.beginBlock("Exporting results...");

            for( unsigned int ii = 0; ii < size_surfels; ++ii )
            {
              file << v_curvatures[ii] << std::endl;
              //                            file2 << v_estimated_radius[ii] << std::endl;
              //                            file4 << v_radius[v_registration[ii]] << std::endl;
            }

            trace.endBlock();
          }

          /*trace.beginBlock( "II Principal curvatures" );

          typedef IntegralInvariantGaussianCurvatureEstimator< KSpace, MySpelFunctor > Estimator;
          Estimator estimator( K, functor );

          c.startClock();
          estimator.init( h, re_convolution_kernel );

          char full_filename[360];
          sprintf( full_filename, "%s%s", filename.c_str(), "_II_principal.dat" );
          std::ofstream file( full_filename );
          file << "# h = " << h << std::endl;
          file << "# Gaussian Curvature estimation from the Integral Invariant" << std::endl;
          file << "# computed kernel radius = " << re_convolution_kernel << std::endl;

          std::ostream_iterator< CurvatureInformations > out_it_ii_principal( file, "\n" );

          range = new VisitorRange( new Visitor( surf, *surf.begin() ));
          ibegin = range->begin();
          iend = range->end();

          estimator.evalPrincipalCurvatures ( ibegin, iend, out_it_ii_principal );

          double TIIGaussCurv = c.stopClock();
          file << "# time = " << TIIGaussCurv << std::endl;
          file.close();

          trace.endBlock();*/
        }

      }

      // Monge
      if( options.at( 2 ) != '0' )
      {
        // Monge Mean Curvature
        if( properties.at( 0 ) != '0' )
        {
          trace.beginBlock( "Monge mean curvature" );

          typedef MongeJetFittingMeanCurvatureEstimator<Surfel, CanonicSCellEmbedder<KSpace> > FunctorMean;
          typedef ConstValueFunctor< double > ConvFunctor;
          typedef LocalEstimatorFromSurfelFunctorAdapter<typename MyDigitalSurface::DigitalSurfaceContainer, Z3i::L2Metric, FunctorMean, ConvFunctor> ReporterH;
          CanonicSCellEmbedder<KSpace> embedder( K );
          FunctorMean estimatorH( embedder, h );
          ConvFunctor convFunc(1.0);
          ReporterH reporterH(surf.container(), Z3i::l2Metric, estimatorH, convFunc);
          c.startClock();
          reporterH.init( h , re_convolution_kernel / h  );

          char full_filename[360];
          sprintf( full_filename, "%s%s", filename.c_str(), "_MongeJetFitting_mean.dat" );
          std::ofstream file( full_filename );
          file << "# h = " << h << std::endl;
          file << "# Mean Curvature estimation from CGAL Monge from and Jet Fitting" << std::endl;
          file << "# computed kernel radius = " << re_convolution_kernel << std::endl;
          std::ostream_iterator< double > out_it_monge_mean( file, "\n" );

          range = new VisitorRange( new Visitor( surf, *surf.begin() ));
          ibegin = range->begin();
          iend = range->end();

          reporterH.eval(ibegin, iend , out_it_monge_mean);
          double TMongeMeanCurv = c.stopClock();
          file << "# time = " << TMongeMeanCurv << std::endl;
          file.close();
          //          delete range;


          trace.endBlock();
        }

        // Monge Gaussian Curvature
        if( properties.at( 1 ) != '0' )
        {
          trace.beginBlock( "Monge Gaussian curvature" );

          typedef MongeJetFittingGaussianCurvatureEstimator<Surfel, CanonicSCellEmbedder<KSpace> > FunctorGaussian;
          typedef ConstValueFunctor< double > ConvFunctor;
          typedef LocalEstimatorFromSurfelFunctorAdapter<typename MyDigitalSurface::DigitalSurfaceContainer, Z3i::L2Metric, FunctorGaussian, ConvFunctor> ReporterK;
          CanonicSCellEmbedder<KSpace> embedder( K );
          FunctorGaussian estimatorK( embedder, h );
          ConvFunctor convFunc(1.0);
          ReporterK reporterK(surf.container(), Z3i::l2Metric, estimatorK, convFunc);
          c.startClock();
          reporterK.init( h , re_convolution_kernel / h  );

          range = new VisitorRange( new Visitor( surf, *surf.begin() ));
          ibegin = range->begin();
          iend = range->end();

          char full_filename[360];
          sprintf( full_filename, "%s%s", filename.c_str(), "_MongeJetFitting_gaussian.dat" );
          std::ofstream file( full_filename );
          file << "# h = " << h << std::endl;
          file << "# Gaussian Curvature estimation from CGAL Monge from and Jet Fitting" << std::endl;
          file << "# computed kernel radius = " << re_convolution_kernel << std::endl;
          std::ostream_iterator< double > out_it_monge_gaussian( file, "\n" );
          reporterK.eval(ibegin, iend , out_it_monge_gaussian);
          double TMongeGaussCurv = c.stopClock();
          file << "# time = " << TMongeGaussCurv << std::endl;
          file.close();
          //          delete range;


          trace.endBlock();
        }

        // Monge Principal Curvatures
        if( properties.at( 2 ) != '0' )
        {
          trace.beginBlock( "Monge Principal Curvature" );

          typedef MongeJetFittingPrincipalCurvaturesEstimator<Surfel, CanonicSCellEmbedder<KSpace> > FunctorPrincCurv;
          typedef ConstValueFunctor< double > ConvFunctor;
          typedef LocalEstimatorFromSurfelFunctorAdapter<typename MyDigitalSurface::DigitalSurfaceContainer, Z3i::L2Metric, FunctorPrincCurv, ConvFunctor> ReporterK;
          CanonicSCellEmbedder<KSpace> embedder( K );
          FunctorPrincCurv estimatorK( embedder, h );
          ConvFunctor convFunc(1.0);
          ReporterK reporterK(surf.container(), Z3i::l2Metric, estimatorK, convFunc);
          c.startClock();
          reporterK.init( h , re_convolution_kernel / h  );

          range = new VisitorRange( new Visitor( surf, *surf.begin() ));
          ibegin = range->begin();
          iend = range->end();

          char full_filename[360];
          sprintf( full_filename, "%s%s", filename.c_str(), "_MongeJetFitting_principal.dat" );
          std::ofstream file( full_filename );
          file << "# h = " << h << std::endl;
          file << "# Gaussian Curvature estimation from CGAL Monge from and Jet Fitting" << std::endl;
          file << "# computed kernel radius = " << re_convolution_kernel << std::endl;
          std::ostream_iterator< CurvatureInformations > out_it_monge_gaussian( file, "\n" );
          reporterK.eval(ibegin, iend , out_it_monge_gaussian);
          double TMongeGaussCurv = c.stopClock();
          file << "# time = " << TMongeGaussCurv << std::endl;
          file.close();
          //          delete range;


          trace.endBlock();
        }

      }
    }
  }
  catch ( InputException e )
  {
    std::cerr << "[estimatorCurvatureComparator3D]"
              << " error."
              << e.what()
              << " "
              << filename << " "
              << border_min[0] << " "
              << border_max[0] << " "
              << h << " "
                 //              << radius_kernel << " "
              << lambda_optimized << " "
              << options << " "
                 //              << alpha << " "
              << std::endl;
    return false;
  }

  return true;
}


/**
 * Missing parameter error message.
 *
 * @param param
 */
void missingParam( std::string param )
{
  trace.error() << " Parameter: " << param << " is required.";
  trace.info() << std::endl;
  exit( 1 );
}

namespace po = boost::program_options;

int main( int argc, char** argv )
{


#ifndef WITH_CGAL
#error You need to have activated CGAL (WITH_CGAL) to include this file.
#endif
#ifndef WITH_EIGEN
#error You need to have activated EIGEN (WITH_EIGEN) to include this file.
#endif

  // parse command line ----------------------------------------------
  po::options_description general_opt("Allowed options are");
  general_opt.add_options()
      ("help,h", "display this message")
      ("shape,s", po::value< std::string >(), "Shape")
      ("output,o", po::value< std::string >(), "Output file")
      //      ("radius,r",  po::value< double >(), "Kernel radius for IntegralInvariant" )
      //      ("alpha",  po::value<double>()->default_value(1.0/3.0), "Alpha parameter for Integral Invariant computation" )
      ("h",  po::value< double >(), "Grid step" )
      ("minAABB,a",  po::value< double >()->default_value( -10.0 ), "Min value of the AABB bounding box (domain)" )
      ("maxAABB,A",  po::value< double >()->default_value( 10.0 ), "Max value of the AABB bounding box (domain)" )
      ("noise,n",  po::value<double>()->default_value(0.0), "Level of noise to perturb the shape" )
      ("lambda,l",  po::value< bool >()->default_value( false ), "Use the shape to get a better approximation of the surface (optional)" )
      ("properties",  po::value<std::string>()->default_value("110"), "the i-th property is disabled iff there is a 0 at position i" )
      ("testsII,i",  po::value<std::string>()->default_value("11111111"), "the i-th test for II estimator is disabled iff there is a 0 at position i" )
      ("nbKernels", po::value<unsigned int>()->default_value(0), "Nb of kernels to use. 0 by default (aka all the kernels computed")
      ("constante,k", po::value<double>()->default_value(0.1), "Constante")
      ("modeSegments,m",  po::value< std::string >()->default_value("max"), "min, mean, max" )
      ("typeSegments,t",  po::value< std::string >()->default_value("normal"), "normal, DPS" )
      ("estimators,e",  po::value< std::string >()->default_value("110"), "the i-th estimator is disabled iff there is a 0 at position i" );


  bool parseOK = true;
  po::variables_map vm;
  try
  {
    po::store( po::parse_command_line( argc, argv, general_opt ), vm );
  }
  catch( const std::exception & ex )
  {
    parseOK = false;
    trace.info() << "Error checking program options: " << ex.what() << std::endl;
  }
  po::notify( vm );
  if( !parseOK || vm.count("help") || argc <= 1 )
  {
    trace.info()<< "Compare local estimators on implicit shapes using DGtal library" <<std::endl
                << "Basic usage: "<<std::endl
                << "\t3dlocalEstimators --shape <shape> --h <h> --radius <radius> --estimators <binaryWord> --output <output>"<<std::endl
                << std::endl
                << "Below are the different available families of estimators: " << std::endl
                << "\t - Integral Invariant Mean" << std::endl
                << "\t - Integral Invariant Gaussian" << std::endl
                << "\t - Monge Jet Fitting Mean" << std::endl
                << "\t - Monge Jet Fitting Gaussian" << std::endl
                << std::endl
                << "The i-th family of estimators is enabled if the i-th character of the binary word is not 0. "
                << "The default binary word is '1100'. This means that the first family of estimators, "
                << "ie. Integral Invariant, is enabled, whereas the next ones are disabled. "
                << "Below are the different available properties: " << std::endl
                << "\t - Mean Curvature" << std::endl
                << "\t - Gaussian Curvature" << std::endl
                << "\t - k1/k2" << std::endl
                << "Below are the different available testII: " << std::endl
                << "\t - Global Max" << std::endl
                << "\t - Global Mean" << std::endl
                << "\t - Global Median" << std::endl
                << "\t - Global Min" << std::endl
                << "\t - Local Max" << std::endl
                << "\t - Local Mean" << std::endl
                << "\t - Local Median" << std::endl
                << "\t - Local Min" << std::endl
                << std::endl;
    return 0;
  }

  if (!(vm.count("output"))) missingParam("--output");
  if (!(vm.count("shape"))) missingParam("--shape");
  if (!(vm.count("h"))) missingParam("--h");

  std::string file_export = vm["output"].as< std::string >();
  int nb = 3;
  std::string options = vm["estimators"].as< std::string >();
  if (options.size() < nb)
  {
    trace.error() << " At least " << nb
                  << " characters are required "
                  << " with option --estimators.";
    trace.info() << std::endl;
    exit(1);
  }
  double h = vm["h"].as< double >();
  std::string poly_str = vm["shape"].as< std::string >();
  bool lambda_optimized = vm["lambda"].as< bool >();
  double noiseLevel = vm["noise"].as<double>();

  nb = 3; //number of available properties
  std::string properties = vm["properties"].as<std::string>();
  if (properties.size() < nb)
  {
    trace.error() << " At least " << nb
                  << " characters are required "
                  << " with option --properties.";
    trace.info() << std::endl;
    exit(1);
  }

  typedef Z3i::Space::RealPoint RealPoint;
  typedef Z3i::Space::RealPoint::Coordinate Ring;

  RealPoint border_min( vm["minAABB"].as< double >(), vm["minAABB"].as< double >(), vm["minAABB"].as< double >() );
  RealPoint border_max( vm["maxAABB"].as< double >(), vm["maxAABB"].as< double >(), vm["maxAABB"].as< double >() );

  /// Construction of the polynomial shape

  typedef MPolynomial< 3, Ring > Polynomial3;
  typedef MPolynomialReader<3, Ring> Polynomial3Reader;
  typedef ImplicitPolynomial3Shape<Z3i::Space> ImplicitShape;

  Polynomial3 poly;
  Polynomial3Reader reader;
  std::string::const_iterator iter = reader.read( poly, poly_str.begin(), poly_str.end() );
  if ( iter != poly_str.end() )
  {
    std::cerr << "ERROR: I read only <"
              << poly_str.substr( 0, iter - poly_str.begin() )
              << ">, and I built P=" << poly << std::endl;
    return 1;
  }

  ImplicitShape* shape = new ImplicitShape( poly );

  struct OptionsIntegralInvariant optII;
  optII.constante =  vm["constante"].as< double >();
  optII.nbKernels =  vm["nbKernels"].as< unsigned int >();
  nb = 8; //number of available tests for II
  optII.tests = vm["testsII"].as<std::string>();
  if( optII.tests.size() < nb )
  {
    trace.error() << " At least " << nb
                  << " characters are required "
                  << " with option --testsII.";
    trace.info() << std::endl;
    exit(1);
  }
  optII.modeSegments = vm["modeSegments"].as< std::string >();
  optII.typeSegments = vm["typeSegments"].as< std::string >();

  compareShapeEstimators< Z3i::Space, ImplicitShape > (
        file_export,
        shape,
        border_min, border_max,
        h,
        //        radius,
        //        alpha,
        options,
        properties,
        lambda_optimized,
        optII,
        noiseLevel );

  delete shape;
}
