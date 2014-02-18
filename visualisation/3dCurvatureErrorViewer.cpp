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
#include "DGtal/helpers/StdDefs.h"

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

#include "DGtal/math/KMeans.h"

//shapes
//#include "DGtal/shapes/implicit/ImplicitBall.h"
#include "DGtal/math/MPolynomial.h"
#include "DGtal/io/readers/MPolynomialReader.h"
#include "DGtal/shapes/implicit/ImplicitPolynomial3Shape.h"

//Digitizer
#include "DGtal/shapes/GaussDigitizer.h"
#include "DGtal/topology/LightImplicitDigitalSurface.h"
#include "DGtal/topology/DigitalSurface.h"
#include "DGtal/graph/DepthFirstVisitor.h"
#include "DGtal/graph/DistanceBreadthFirstVisitor.h"
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
//#include "DGtal/geometry/surfaces/estimation/IntegralInvariantGaussianCurvatureEstimator.h"

#if 0
#include "DGtal/geometry/surfaces/estimation/LocalEstimatorFromSurfelFunctorAdapter.h"
#include "DGtal/geometry/surfaces/estimation/estimationFunctors/MongeJetFittingMeanCurvatureEstimator.h"
#include "DGtal/geometry/surfaces/estimation/estimationFunctors/MongeJetFittingGaussianCurvatureEstimator.h"
#include "DGtal/geometry/surfaces/estimation/estimationFunctors/MongeJetFittingPrincipalCurvaturesEstimator.h"
#endif


// Drawing
//#include "DGtal/io/boards/Board3D.h"
#include "DGtal/io/viewers/Viewer3D.h"
#include "DGtal/io/colormaps/GradientColorMap.h"
#include <QtGui/QApplication>

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

/*template < typename Shape, typename KSpace, typename ConstIterator, typename OutputIterator >
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
*/
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

template<typename Point2, typename Point1>
struct Myfunctor
{
  typedef Point2 Output;
  Point2 operator()(const Point1&a) const
  {
    return Point2(a.myCoordinates[0]/2,a.myCoordinates[1]/2);
  }
};

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

/*Point2D plongement2D( Surfel s, Dimension i )
{
  KSpace & K = l'espace cellulaire
  Dimension k = K.sOrthDir( s );
  Dimension j = (k+1)%3;
  if ( j == i ) j = (i+1)%3;
  SCell next_linel = K.sDirectIncident( s, j );
  SCell base_pointel = K.sIncident( next_linel, i, false );
  Point3D p = K.sCoords( base_pointel ); // une des coords est inutile
  Point2D q( p[(i+1)%3], p[(i+2)%3] );
  return q;
}*/

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
      if( marque[ currentSurfel ] )
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
        typedef typename MySlice::Iterator Iterator3D;
//        typedef SCellProjector< KhalimskySpaceND<2,int> > Functor;
//        typedef CanonicSCellEmbedder< Functor::KSpace > Functor2;
        typedef SCellToMyPoint< typename MySlice::KSpace, Z2i::KSpace > Functor2;
//        typedef Myfunctor<typename Functor::KSpace::Point , typename Functor::KSpace::SCell> Functor2;


        /*typedef ConstIteratorAdapter< ConstCirculator3D, Functor, Functor::SCell > ConstIterator2D;
        typedef ConstIteratorAdapter< ConstIterator3D, Functor, Functor::SCell > ConstIterator2D2;
        typedef ConstIteratorAdapter< ConstIterator2D, Functor2, Functor2::Output > ConstIterator2DP;
        typedef ConstIteratorAdapter< ConstIterator2D2, Functor2, Functor2::Output > ConstIterator2DP2;*/

        typedef ConstIteratorAdapter< ConstIterator3D, Functor2 > ConstIterator3D2P;

        ConstIterator3D a = slice.begin();
        ConstIterator3D b = ++(slice.begin());
        Dimension dimm = 0;
        while( a->myCoordinates[dimm] != b->myCoordinates[dimm] )
        {
          ++dimm;
        }

        /*Functor projector;
        projector.initRemoveOneDim( dimm );
        ConstIterator2D xbegin( slice.c(), projector );
        ConstIterator2D xend( slice.c(), projector );

        ConstIterator2D2 xxbegin( slice.begin(), projector );
        ConstIterator2D2 xxend( slice.end(), projector );

        Functor::KSpace k2d;*/
//        Functor2 pointFunctor;
        Functor2 pointFunctor(K, dim);

        ConstIterator3D2P ppbegin( slice.begin(), pointFunctor );
        ConstIterator3D2P ppend( slice.end(), pointFunctor );

        /*ConstIterator2DP pbegin( xbegin, pointFunctor );
        ConstIterator2DP pend( xend, pointFunctor );

        ConstIterator2DP2 ppbegin( xxbegin, pointFunctor );
        ConstIterator2DP2 ppend( xxend, pointFunctor );*/


        const Dimension size_slice = slice.size();
        std::vector< Statistic< double > > v_statMSEL(size_slice);
//        std::vector< Statistic< double > > v_statMSEL2(size_slice);

        typedef Z2i::Point Point;
        std::vector< Point > pts;

//        if (K.sKCoords(*slice.begin())[0] == 1)
//        trace.warning()<< "Slice begin"<< K.sKCoords(*slice.begin())<<std::endl;
//        do
//        {
////          if (K.sKCoords(*slice.begin())[0] == 1)
////          {
////            trace.info() << *pbegin << " ";
////          }
//            pts.push_back( *ppbegin );
//          ++ppbegin;
//        } while( ppbegin != ppend );
//        if (K.sKCoords(*slice.begin())[0] == 1)
//          trace.info()<< "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"<<std::endl;

//        SCell bel = Surfaces< Functor::KSpace >::findABel( k2d, slice, 10000 );
//        SurfelAdjacency< Functor::KSpace::dimension > SAdj( true );
//        Surfaces< Functor::KSpace >::track2DBoundary( pts, k2d, SAdj, slice, bel );

//        typedef GridCurve<  Functor::KSpace >::PointsRange PointRange;
//        GridCurve<  Functor::KSpace > gridcurve;
//        gridcurve.initFromSCellsVector( pts );
//        PointRange pr = gridcurve.getPointsRange();



        Circulator< ConstIterator3D2P > cbegin( ppbegin, ppbegin, ppend );
        Circulator< ConstIterator3D2P > cend( cbegin );

        for( Dimension ii = 0; ii < size_slice; ++ii )
        {
          v_statMSEL[ii] = Statistic<double>(true);
        }


//        trace.error() << "size " << size_slice << std::endl;
//        if( size_slice <= 4 )
//        {
//          for( Dimension ii = 0; ii < size_slice; ++ii )
//          {
//            v_statMSEL[ii].addValue(0);
//          }
////          analyseAllLengthMS<Functor::KSpace>( v_statMSEL, ppbegin, ppend );
//        }
//        else
        {
          analyseAllLengthMS<Z2i::KSpace>( v_statMSEL, cbegin, cend );
//          Circulator< std::vector< Point >::iterator > cbegin2( pts.begin(), pts.begin(), pts.end() );
//          int shift = 0;//size_slice / 2;
//          for( int i = 0; i < shift; ++i )
//            ++cbegin2;
//          Circulator< std::vector< Point >::iterator > cend2( cbegin2 );
//          analyseAllLengthMS<Functor::KSpace>( v_statMSEL, cbegin2, cend2 );
//          for( Dimension ii = 0; ii < size_slice; ++ii )
//          {
//            v_statMSEL[ii].addValues( v_statMSEL2[ ( ii + shift ) % size_slice ].begin(), v_statMSEL2[ ( ii + shift ) % size_slice ].end());
//          }
        }

        //                for(Dimension ii = 0; ii < pr2size; ++ii )
        //                {
        //                  v_statMSEL[ii].terminate();
        //                }

        Dimension iii = 0;
//        ConstIterator3D sit = slice.begin();
//        ConstIterator3D send = slice.end();
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
            trace.error() << "WHHHHHAT" << std::endl;
          }
          ++iii;
//          ++sit;
        }// while( sit != send );
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
                      const double minRadiusAABB,
                      const unsigned int minValue = 5,
                      const double defaultValue = 10.0 )
{
  //    if( re/h <= minValue ) /// 	ridiculously small radius check
  //    {
  //  //    trace.error() << "re small " << re << std::endl;
  //  //    re = 5.0 * h;
  //      re = defaultValue;//*h;
  //    }
  if( re > ( 0.75 * minRadiusAABB ))
  {
    re = 0.75 * minRadiusAABB;
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
                         const double minRadiusAABB,
                         const std::string & prop,
                         const std::string & mode,
                         const unsigned int minValue = 5,
                         const double defaultValue = 10.0 )
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

    //if( mode == "max" )
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
   /* else if( mode == "min" )
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
    }*/

    double result = -1.0;
    /*if( prop == "min" )
    {
      result = stat.min();
    }
    else if( prop == "mean" )*/
    {
      result = stat.min();
//      if( result < minValue )
//      {
//        result = defaultValue;
//      }
    }
   /* else if( prop == "median" )
    {
      result = stat.median();
    }
    else if( prop == "max" )
    {
      result = stat.max();
    }*/

    ASSERT(( result != -1.0 ));

    double re = -1.0;
   /* if( result == defaultValue )
      re = defaultValue;
    else*/
      re = constante * result * result * h;

    //checkSizeRadius( re, h, minRadiusAABB );

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
compareShapeEstimators( int argc, char** argv,
                        const Shape * aShape,
                        const typename Space::RealPoint & border_min,
                        const typename Space::RealPoint & border_max,
                        const double & h,
                        struct OptionsIntegralInvariant optionsII )
{
  typedef typename Space::RealPoint RealPoint;
  typedef GaussDigitizer< Z3i::Space, Shape > DigitalShape;
  typedef Z3i::Domain Domain;
  typedef Z3i::KSpace KSpace;
  typedef typename KSpace::SCell SCell;
  typedef typename KSpace::Surfel Surfel;

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

      std::vector< double > v_curvatures_true;
      std::back_insert_iterator< std::vector< double > > v_out_it_true( v_curvatures_true );
      std::vector< double > v_curvatures_II;
      std::back_insert_iterator< std::vector< double > > v_out_it_II( v_curvatures_II );
      std::vector< Surfel > surfels;

      // Estimations
      Clock c;

      ///// True
      {
        // True Mean Curvature
        {
          trace.beginBlock( "True mean curvature" );

          range = new VisitorRange( new Visitor( surf, *surf.begin() ));
          ibegin = range->begin();
          iend = range->end();

          c.startClock();

          estimateTrueMeanCurvatureQuantity( ibegin,
                                             iend,
                                             v_out_it_true, //out_it_true_mean,
                                             K,
                                             h,
                                             aShape );


          trace.endBlock();
        }
      }

      ///// II
      {
        MyPointFunctor pointFunctor( dshape, domain, 1, 0 );
        MySpelFunctor functor( pointFunctor, K );

        double minRadiusAABB = ( border_max[0] - border_min[0] ) / 2.0; // this is a box.

        std::map< Surfel, std::pair< Statistic< double >, Statistic< double > > > segments;
//        Statistic< double > allSegments( true );

        trace.beginBlock("Extracting all surfels...");

        {
          VisitorRange ranger( new Visitor( surf, *surf.begin() ) );
          VisitorConstIterator abegin = ranger.begin();
          VisitorConstIterator aend = ranger.end();

          surfels.clear();
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
//        if( true ) //optionsII.tests.at(0) != '0' || optionsII.tests.at(1) != '0' )
//        {
//          computeGlobalSegment( surfels, segments, allSegments, "max" );
//          allSegments.terminate();

//          if( optionsII.tests.at(2) == '0' && optionsII.tests.at(3) == '0' )
//          {
//            //              surfels.clear();
//            //              segments.clear();
//          }
//        }

        // Integral Invariant Mean Curvature
        {
          /////// local mean
          {
            trace.beginBlock("Computation of radius...");

            std::vector< double > v_estimated_radius;
            v_estimated_radius.resize( size_surfels );
            //            double global_mean = allSegments.mean();
            //            double global_re = optionsII.constante * global_mean * global_mean * h;
            double min_re = 5;//optionsII.constante * 2 * 2 * h;
            trace.error() << "min re " << min_re << std::endl;
            Dimension nbKernelsRadius = computeRadius< Surfel >( surfels, segments, optionsII.constante, h, v_estimated_radius, minRadiusAABB, "mean", "max", min_re, -42 );
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
            v_curvatures_II.resize( size_surfels );

#ifdef WITH_OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
            for( Dimension ii = 0; ii < size_surfels; ++ii )
            {
              if( optionsII.nbKernels > 0 )
              {
                //                v_curvatures[ii] = v_estimators[ v_registration[ ii ]]->eval( surfels.begin() + ii );
//                if( v_radius[ v_registration[ii]] == -42 )
//                {
//                  v_curvatures_II[ii] = 0.25;//-42;
//                }
//                else
//                {
                  v_curvatures_II[ii] = v_estimators[ v_registration[ ii ]]->eval( surfels.begin() + ii );
//                  v_curvatures_II[ii] = v_radius[v_registration[ ii ]];
//                }
              }
              else
              {
//                if( v_estimated_radius[ii] == -42 )
//                {
//                  v_curvatures_II[ii] = -42;
//                }
//                else
//                {
                Estimator estimator( K, functor );
                estimator.init( h, v_estimated_radius[ii] );
//                                  v_curvatures[ii] = estimator.eval( surfels.begin() + ii );
                v_curvatures_II[ii] = estimator.eval( surfels.begin() + ii );
//                  v_curvatures_II[ii] = v_estimated_radius[ii];
//                }
              }
            }

            trace.endBlock();

            trace.beginBlock("Exporting results...");

            for( unsigned int ii = 0; ii < size_surfels; ++ii )
            {
              //              file << v_curvatures[ii] << std::endl;
              //                            file2 << v_estimated_radius[ii] << std::endl;
              //                            file4 << v_radius[v_registration[ii]] << std::endl;
            }

            trace.endBlock();
          }

        }
      }
      ///////////////////////
      std::vector< double > v_error( v_curvatures_II.size());
      double min = 9999999;
      double max = -500000;
      for( Dimension ii = 0; ii < v_curvatures_II.size(); ++ii )
      {
//        if(v_curvatures_II[ii] == -42)
//        {
//          v_error[ii] = -42;
//          continue;
//        }
//        v_error[ii] = v_curvatures_II[ii];
        v_error[ii] = std::abs ( v_curvatures_true[ii] - v_curvatures_II[ii] );

        if( v_error[ii] < min )
          min = v_error[ii];


        if( v_error[ii] > max )
          max = v_error[ii];
      }

      trace.info() << "min: " << min << " max: " << max << std::endl;

      QApplication application( argc, argv );
      typedef Viewer3D< Z3i::Space, Z3i::KSpace > Viewer;
      //typedef Board3D< Z3i::Space, Z3i::KSpace> Board3D;
      Viewer viewer( K );
      viewer.show();

      typedef GradientColorMap< double > Gradient;
      Gradient cmap_grad( min, max );
      cmap_grad.addColor( Color( 50, 50, 255 ) );
      cmap_grad.addColor( Color( 255, 0, 0 ) );
      cmap_grad.addColor( Color( 255, 255, 10 ) );

      viewer << DGtal::SetMode3D(bel.className(), "Basic");
      const Color  AXIS_COLOR_GREEN( 20, 200, 20, 255 );

      for ( unsigned int ii = 0; ii < surfels.size(); ++ii )
      {
     //  if( K.sKCoords(surfels[ii])[0] == 1)
        {

//        if( v_error[ii] == -42 )
//        {
//          viewer << CustomColors3D( AXIS_COLOR_GREEN, AXIS_COLOR_GREEN)
//                 << surfels[ii];
//        }
//        else
        {
          viewer << CustomColors3D( Color::Black, cmap_grad( v_error[ ii ] ))
                 << surfels[ii];
        }
      }
      }

      viewer << Viewer3D<>::updateDisplay;
      return application.exec();
      ///////////////////////
    }


  }
  catch ( InputException e )
  {
    std::cerr << "[estimatorCurvatureComparator3D]"
              << " error."
              << e.what()
              << " "
//              << filename << " "
              << border_min[0] << " "
              << border_max[0] << " "
              << h << " "
                 //              << radius_kernel << " "
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
//      ("output,o", po::value< std::string >(), "Output file")
      //      ("radius,r",  po::value< double >(), "Kernel radius for IntegralInvariant" )
      //      ("alpha",  po::value<double>()->default_value(1.0/3.0), "Alpha parameter for Integral Invariant computation" )
      ("h",  po::value< double >(), "Grid step" )
      ("minAABB,a",  po::value< double >()->default_value( -10.0 ), "Min value of the AABB bounding box (domain)" )
      ("maxAABB,A",  po::value< double >()->default_value( 10.0 ), "Max value of the AABB bounding box (domain)" )
      ("nbKernels", po::value<unsigned int>()->default_value(0), "Nb of kernels to use. 0 by default (aka all the kernels computed")
      ("constante,k", po::value<double>()->default_value(0.1), "Constante")
      ("typeSegments,t",  po::value< std::string >()->default_value("normal"), "normal, DPS" );


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

//  if (!(vm.count("output"))) missingParam("--output");
  if (!(vm.count("shape"))) missingParam("--shape");
  if (!(vm.count("h"))) missingParam("--h");

//  std::string file_export = vm["output"].as< std::string >();
  double h = vm["h"].as< double >();
  std::string poly_str = vm["shape"].as< std::string >();

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
  optII.typeSegments = vm["typeSegments"].as< std::string >();

  compareShapeEstimators< Z3i::Space, ImplicitShape > (
        argc, argv,
//        file_export,
        shape,
        border_min, border_max,
        h,
        optII );

  delete shape;
}
