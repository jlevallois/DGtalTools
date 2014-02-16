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
#include <cstring>

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

#include "DGtal/math/Statistic.h"

// Shape constructors
#include "DGtal/io/readers/VolReader.h"
#include "DGtal/images/ImageSelector.h"
#include "DGtal/images/imagesSetsUtils/SetFromImage.h"
#include "DGtal/images/imagesSetsUtils/SimpleThresholdForegroundPredicate.h"
#include "DGtal/topology/SurfelAdjacency.h"
#include "DGtal/topology/helpers/Surfaces.h"
#include "DGtal/topology/LightImplicitDigitalSurface.h"
#include "DGtal/images/ImageHelper.h"
#include "DGtal/topology/DigitalSurface.h"
#include "DGtal/graph/DepthFirstVisitor.h"
#include "DGtal/graph/GraphVisitorRange.h"

// Segments
#include "DGtal/geometry/curves/SaturatedSegmentation.h"

// Integral Invariant includes
#include "DGtal/geometry/surfaces/FunctorOnCells.h"
#include "DGtal/geometry/surfaces/estimation/IntegralInvariantMeanCurvatureEstimator.h"
#include "DGtal/geometry/surfaces/estimation/IntegralInvariantGaussianCurvatureEstimator.h"

// Drawing
#include "DGtal/io/boards/Board2D.h"
#include "DGtal/io/colormaps/GradientColorMap.h"
#include "DGtal/topology/DigitalSurface2DSlice.h"

#ifdef WITH_VISU3D_QGLVIEWER
#include <QtGui/QApplication>
#include "DGtal/io/viewers/Viewer3D.h"
#endif

using namespace std;
using namespace DGtal;

const Color  AXIS_COLOR_RED( 200, 20, 20, 255 );
const Color  AXIS_COLOR_GREEN( 20, 200, 20, 255 );
const Color  AXIS_COLOR_BLUE( 20, 20, 200, 255 );
const double AXIS_LINESIZE = 0.05;


///////////////////////////////////////////////////////////////////////////////

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
  typedef std::vector< SegmentComputerIterator > VectorOfSegmentComputerIterator;
  typedef std::map< Point, VectorOfSegmentComputerIterator > Pmap;

  // Computes the tangential cover
  SegmentComputer algo;
  Decomposition theDecomposition( itb, ite, algo);

  Pmap map;
  for( Iterator itc = itb; itc != ite; ++itc )
  {
    map.insert( std::pair< Point, VectorOfSegmentComputerIterator >( *itc, VectorOfSegmentComputerIterator() ) );
  }


  for ( SegmentComputerIterator scIt = theDecomposition.begin(), scItEnd = theDecomposition.end();
        scIt != scItEnd; ++scIt )
  {
    const SegmentComputer & sc = *scIt;
    for ( Iterator ptIt = sc.begin(), ptItEnd = sc.end(); ptIt != ptItEnd; ++ptIt )
    {
      typename Pmap::iterator mloc = map.find( *ptIt );
      if( mloc != map.end() )
      {
        mloc->second.push_back( scIt );
      }
    }
  }

  Dimension ii = 0;
  for( Iterator itc = itb; itc != ite; ++itc )
  {
    //statD[ii].clear();
    statE[ii].clear();
    typename Pmap::iterator mloc = map.find( *itc );
    ASSERT(( mloc != map.end() ));

    /////////////
    for( typename VectorOfSegmentComputerIterator::iterator scIt = mloc->second.begin(), scItEnd = mloc->second.end(); scIt != scItEnd; ++scIt )
    {
      const SegmentComputer & sc = *(*scIt);
      /*int64_t l = 0;
          for ( Iterator ptIt = sc.begin(), ptItEnd = sc.end(); ptIt != ptItEnd; ++ptIt )
            ++l;
          statD[ii].addValue( (double) l );*/
      Vector v = *( sc.end() - 1 ) - *( sc.begin() );
      statE[ii].addValue( v.norm() );
      //          std::cout << " v=" << v.norm() << std::endl;
    }
    /////////////

    ++ii;
  }
}

template< typename ImplicitDigitalSurface, typename Surfel >
void computeSegments( std::vector< Surfel > & surfels,
                      std::map< Surfel*, std::pair< Statistic< double >, Statistic< double > > > & segments,
                      const Z3i::KSpace & K,
                      const ImplicitDigitalSurface & impDigitalSurface )
{
  typedef typename ImplicitDigitalSurface::Tracker Tracker;
  typedef DigitalSurface2DSlice< Tracker > MySlice;
  typedef std::pair< Statistic< double >, Statistic< double > > PairOfStatistics;
  typedef std::map< Surfel*, PairOfStatistics > SurfelMap;
  typedef std::map< Surfel*, bool > MarqueMap;

  const Dimension surfels_size = surfels.size();

  for( Dimension dim = 0; dim < 3; ++dim )
  {
    MarqueMap marque;
    for( Dimension ii = 0; ii < surfels_size; ++ii )
    {
      marque[ &surfels[ ii ] ] = false;
    }

    for( Dimension ii = 0; ii < surfels_size; ++ii )
    {
      Surfel currentSurfel = surfels[ ii ];
      if( marque[ &surfels[ ii ] ] == true )
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
        typedef typename MySlice::Iterator Iterator3D;
        typedef SCellProjector< KhalimskySpaceND<2, int> > Functor;
        typedef SCellToPoint< Functor::KSpace > Functor2;
        typedef ConstIteratorAdapter< ConstIterator3D, Functor, Functor::SCell > ConstIterator2D;
        typedef ConstIteratorAdapter< ConstIterator2D, Functor2, Functor2::Output > ConstIterator2DP;

        Iterator3D a = slice.begin();
        Iterator3D b = ++(slice.begin());
        Dimension dimm = 0;
        while( a->myCoordinates[dimm] != b->myCoordinates[dimm] )
        {
          ++dimm;
        }

        Functor projector;
        projector.initRemoveOneDim( dimm );
        ConstIterator2D xbegin( slice.begin(), projector );
        ConstIterator2D xend( slice.end(), projector );

        Functor::KSpace k2d;
        Functor2 pointFunctor( k2d );

        ConstIterator2DP pbegin( xbegin, pointFunctor );
        ConstIterator2DP pend( xend, pointFunctor );

        const Dimension pr2size = slice.size();
        std::vector< Statistic< double > > v_statMSEL(pr2size);
        for( Dimension ii = 0; ii < pr2size; ++ii )
        {
          v_statMSEL[ii] = Statistic<double>(true);
        }

        analyseAllLengthMS<Functor::KSpace, ConstIterator2DP>( v_statMSEL, pbegin, pend );

        for(Dimension ii = 0; ii < pr2size; ++ii )
        {
          v_statMSEL[ii].terminate();
        }

        Dimension iii = 0;
        for( Iterator3D sit = slice.begin(), send = slice.end(); sit != send; ++sit )
        {
          Dimension surfel_pos = findSurfel( surfels, *sit );
          ASSERT( surfel_pos != surfels_size );
          if( marque[ &surfels[surfel_pos] ] == false )
          {
            marque[ &surfels[surfel_pos] ] = true;
            ASSERT(( marque.size() == surfels_size ));
            PairOfStatistics & otherpair = segments[ &surfels[surfel_pos] ];
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
          ++iii;
        }
      }
    }
  }

}

template< typename Surfel >
void computeRadius( std::vector< Surfel > & surfels,
                    std::map< Surfel*, std::pair< Statistic< double >, Statistic< double > > > & segments,
                    const double constante,
                    const double h,
                    std::vector< double > & radius )
{
  const Dimension surfels_size = surfels.size();

  ASSERT(( radius.size() == surfels_size ));
  ASSERT(( segments.size() == surfels_size ));

  for( Dimension ii = 0; ii < surfels_size; ++ii )
  {
    Statistic< double > stata = segments[ &surfels[ii] ].first;
    Statistic< double > statb = segments[ &surfels[ii] ].second;

    double median;
    if( stata.max() > statb.max() )
    {
      median = stata.mean();
    }
    else
    {
      median = statb.mean();
    }

    double re = constante * median * median * h;

//    if( re < 5.0 )
//    {
//      re = 5.0;
//    }

    radius[ii] = re;
  }
}

template< typename Quantity, typename Functor, typename Surfel >
void computeCurvature( const Z3i::KSpace & K,
                       const Functor & functor,
                       const std::vector< Surfel > & surfels,
                       const std::vector< double > & radius,
                       std::vector< Quantity > & curvatures )
{
  typedef IntegralInvariantMeanCurvatureEstimator< Z3i::KSpace, Functor > GaussEstimator;

  const Dimension surfels_size = surfels.size();

  ASSERT(( radius.size() == surfels_size ));
  ASSERT(( curvatures.size() == surfels_size ));

#ifdef WITH_OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
  for( Dimension ii = 0; ii < surfels_size; ++ii )
  {
    double re = radius[ ii ];
    //re = ( re < 10.0 )? 10.0 : re;

    GaussEstimator estimator ( K, functor );
    estimator.init( 1.0, radius[ ii ]);
    curvatures[ ii ] = estimator.eval( &surfels[ ii ] );
  }
}

template< typename Surfel, typename Quantity >
int viewCurvature( const std::vector< Surfel > & surfels,
                   const std::vector< Quantity > & curvatures,
                   const Z3i::KSpace & K,
                   int argc,
                   char** argv )
{
#ifdef WITH_VISU3D_QGLVIEWER
  const Dimension surfels_size = surfels.size();

  ASSERT(( curvatures.size() == surfels_size ));

  QApplication application( argc, argv );
  typedef Viewer3D< Z3i::Space, Z3i::KSpace > Viewer;
  Viewer viewer( K );
  viewer.show();
  viewer << DGtal::SetMode3D( surfels[0].className(), "Basic" );

  Quantity min = numeric_limits < Quantity >::max();
  Quantity max = numeric_limits < Quantity >::min();

  for( Dimension ii = 0; ii < surfels_size; ++ii )
  {
    Quantity curvature = curvatures[ii];
    if ( curvature < min )
    {
      min = curvature;
    }
    else if ( curvature > max )
    {
      max = curvature;
    }
  }

  typedef GradientColorMap< Quantity > Gradient;
  Gradient cmap_grad( min, max );
  cmap_grad.addColor( Color( 50, 50, 255 ) );
  cmap_grad.addColor( Color( 255, 0, 0 ) );
  cmap_grad.addColor( Color( 255, 255, 10 ) );

  for( Dimension ii = 0; ii < surfels_size; ++ii )
  {
    viewer << CustomColors3D( Color::Black, cmap_grad( curvatures[ii] ))
           << surfels[ii];
  }

  viewer << Viewer3D<>::updateDisplay;
  return application.exec();
#else
  trace.error() << "OpenGL isn't activate in CMake" << std::endl;
  return -1;
#endif
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
}

namespace po = boost::program_options;

int main2( int argc, char** argv )
{

  /*typedef Z3i::Point Point;
  std::vector< Point > v_test;
  v_test.push_back( Point( 1,1,1 ) );
  v_test.push_back( Point( 3,1,1 ) );
  v_test.push_back( Point( 3,3,1 ) );
  v_test.push_back( Point( 3,6,1 ) );
  v_test.push_back( Point( 3,1,1 ) );

  Projector< SpaceND<2, int> > proj;
  //proj.init( v_test.begin(), v_test.end() );
  proj.initRemoveOneDim( 1 );
  Z2i::Point result = proj( v_test[3] );
  std::cout << "Annnnd the resulting point iiiiiiiis : " << result << std::endl;

  return 0;

  typedef SCellTo2DSCell< std::vector< Point > > Embedder;
  Embedder TEST(v_test);
  typedef Embedder::ConstIterator ConstIterator;
  typedef Embedder::Iterator Iterator;
  ConstIterator itb = TEST.begin();
  ConstIterator ite = TEST.end();

  std::cout << "itrealb=" << *v_test.begin() << std::endl;
  std::cout << "itb=" << *itb << std::endl;
  std::cout << "ite=" << *ite << std::endl;

  ConstIterator itc = itb;
  for( ;
       itc != ite;
       ++itc )
  {
    std::cout << "itc="<< *itc << std::endl;
  }

  std::cout << "_itc=" << *itc << std::endl;
  ++itc;
  std::cout << "_itc=" << *itc << std::endl;

  return 0;*/

  // parse command line ----------------------------------------------
  po::options_description general_opt("Allowed options are");
  general_opt.add_options()
      ("help,h", "display this message")
      ("input-file,i", po::value< std::string >(), ".vol file");

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
  bool neededArgsGiven=true;
  if (!(vm.count("input-file"))){
    missingParam("--input-file");
    neededArgsGiven=false;
  }

  double h = 1.0;

  if(!neededArgsGiven || !parseOK || vm.count("help") || argc <= 1 )
  {
    trace.info()<< "Visualisation of 3d curvature from .vol file using curvature from Integral Invariant" <<std::endl
                << general_opt << "\n"
                << "Basic usage: "<<std::endl
                << "\t3dCurvatureViewer -i <file.vol> --radius <radius> --properties <\"mean\">"<<std::endl
                << std::endl
                << "Below are the different available properties: " << std::endl
                << "\t - \"mean\" for the mean curvature" << std::endl
                << "\t - \"gaussian\" for the Gaussian curvature" << std::endl
                << "\t - \"prindir1\" for the first principal curvature direction" << std::endl
                << "\t - \"prindir2\" for the second principal curvature direction" << std::endl
                << std::endl;
    return 0;
  }

  // Construction of the shape from vol file
  typedef Z3i::Space::RealPoint RealPoint;
  typedef Z3i::Point Point;
  typedef ImageSelector< Z3i::Domain, bool>::Type Image;
  typedef SimpleThresholdForegroundPredicate< Image > ImagePredicate;
  typedef Z3i::KSpace KSpace;
  typedef KSpace::SCell SCell;
  typedef KSpace::Cell Cell;
  typedef KSpace::Surfel KSurfel;
  typedef LightImplicitDigitalSurface< Z3i::KSpace, ImagePredicate > MyLightImplicitDigitalSurface;
  typedef DigitalSurface< MyLightImplicitDigitalSurface > MyDigitalSurface;

  std::string filename = vm["input-file"].as< std::string >();
  Image image = VolReader<Image>::importVol( filename );
  ImagePredicate predicate = ImagePredicate( image, 0 );

  Z3i::Domain domain = image.domain();

  Z3i::KSpace K;

  bool space_ok = K.init( domain.lowerBound(), domain.upperBound(), true );
  if (!space_ok)
  {
    trace.error() << "Error in the Khalimsky space construction."<<std::endl;
    return 2;
  }

  CanonicSCellEmbedder< KSpace > embedder( K );

  SurfelAdjacency< Z3i::KSpace::dimension > SAdj( true );
  KSurfel bel = Surfaces< Z3i::KSpace >::findABel( K, predicate, 100000 );
  MyLightImplicitDigitalSurface LightImplDigSurf( K, predicate, SAdj, bel );
  MyDigitalSurface digSurf( LightImplDigSurf );

  typedef DepthFirstVisitor<MyDigitalSurface> Visitor;
  typedef GraphVisitorRange< Visitor > VisitorRange;
  typedef VisitorRange::ConstIterator SurfelConstIterator;
  typedef ImageToConstantFunctor< Image, ImagePredicate > MyPointFunctor;
  typedef FunctorOnCells< MyPointFunctor, Z3i::KSpace > MyCellFunctor;

  MyPointFunctor pointFunctor( image, predicate, 1 );
  MyCellFunctor functor ( pointFunctor, K ); // Creation of a functor on Cells, returning true if the cell is inside the shape

#ifdef WITH_VISU3D_QGLVIEWER
  QApplication application( argc, argv );
  typedef Viewer3D<Z3i::Space, Z3i::KSpace> Viewer;
  Viewer viewer( K );
  viewer.show();
  //    viewer << SetMode3D(image.domain().className(), "BoundingBox") << image.domain();
#endif


  trace.beginBlock("curvature computation");
  {
    typedef MyLightImplicitDigitalSurface::Tracker Tracker;
    typedef DigitalSurface2DSlice< Tracker > MySlice;
    typedef Tracker::Surfel Surfel;

    VisitorRange range( new Visitor( digSurf, *digSurf.begin() ) );
    SurfelConstIterator abegin = range.begin();
    SurfelConstIterator aend = range.end();

    //    ImageContainerBySTLMap< Z3i::Domain, std::pair< Statistic< double >, Statistic< double > > > mapStats(  );
    std::map< Surfel, std::pair< Statistic< double >, Statistic< double > > > mapStat;

    Surfel oneSCell;
    Dimension ii = 0;
    while( abegin != aend )//for ( unsigned int i = 0; i < results.size(); ++i )
    {
      if( ii == 15 )
      {
        oneSCell = *abegin;
      }
      //viewer << CustomColors3D( Color::White, Color::White )
      //       << *abegin2;
      ++abegin;
      ++ii;
    }

    trace.error() << "oneSCell = " << oneSCell << std::endl;
    trace.error() << "oneSCell orth dir = " << K.sOrthDir( oneSCell ) << std::endl;

    //    Surfel pointel = K.sPointel( K.sKCoords( oneSCell ));

    Tracker* ptrTracker = new Tracker( LightImplDigSurf, oneSCell ); // some pointer on a tracker.
    MySlice slicex( ptrTracker, 0 ); // slice containing x-axis
    MySlice slicey( ptrTracker, 1 ); // slice containing y-axis
    MySlice slicez( ptrTracker, 2 ); // slice containing z-axis

    typedef MySlice::ConstIterator ConstIterator3D;
    typedef SCellProjector< KhalimskySpaceND<2, int> > Functor;
    typedef ConstIteratorAdapter< ConstIterator3D, Functor, Functor::SCell > ConstIterator2D;

    //////////////////// X
    ConstIterator3D a = slicex.begin();
    ConstIterator3D b = ++(slicex.begin());
    Dimension dim = 0;
    while( a->myCoordinates[dim] != b->myCoordinates[dim] )
    {
      ++dim;
    }

    Functor projectx;
    projectx.initRemoveOneDim( dim );
    ConstIterator2D xbegin( slicex.begin(), projectx );
    ConstIterator2D xend( slicex.end(), projectx );

    //////////////////// Y
    a = slicey.begin();
    b = ++(slicey.begin());
    dim = 0;
    while( a->myCoordinates[dim] != b->myCoordinates[dim] )
    {
      ++dim;
    }

    Functor projecty;
    projecty.initRemoveOneDim( dim );
    ConstIterator2D ybegin( slicey.begin(), projecty );
    ConstIterator2D yend( slicey.end(), projecty );

    //////////////////// Z
    a = slicez.begin();
    b = ++(slicez.begin());
    dim = 0;
    while( a->myCoordinates[dim] != b->myCoordinates[dim] )
    {
      ++dim;
    }

    Functor projectz;
    projectz.initRemoveOneDim( dim );
    ConstIterator2D zbegin( slicez.begin(), projectz );
    ConstIterator2D zend( slicez.end(), projectz );

    Board2D boardx;
    Board2D boardy;
    Board2D boardz;

#ifdef WITH_VISU3D_QGLVIEWER
    for( ConstIterator3D itb = slicex.begin(); itb != slicex.end(); ++itb )
    {
      viewer << CustomColors3D( Color::Green, Color::Green )
             << *itb;
    }
#endif
    for( ConstIterator2D itb = xbegin; itb != xend; ++itb )
    {
      boardx << *itb;
    }
#ifdef WITH_VISU3D_QGLVIEWER
    for( ConstIterator3D itb = slicey.begin(); itb != slicey.end(); ++itb )
    {
      viewer << CustomColors3D( Color::Red, Color::Red )
             << *itb;
    }
#endif
    for( ConstIterator2D itb = ybegin; itb != yend; ++itb )
    {
      boardy << *itb;
    }
#ifdef WITH_VISU3D_QGLVIEWER
    for( ConstIterator3D itb = slicez.begin(); itb != slicez.end(); ++itb )
    {
      viewer << CustomColors3D( Color::Blue, Color::Blue )
             << *itb;
    }
#endif
    for( ConstIterator2D itb = zbegin; itb != zend; ++itb )
    {
      boardz << *itb;
    }

    boardx.saveSVG ( "x.svg" );
    boardy.saveSVG ( "y.svg" );
    boardz.saveSVG ( "z.svg" );

    /*    viewer << CustomColors3D( Color::Red, Color::Red )
         << oneSCell*/;

  }

#ifdef WITH_VISU3D_QGLVIEWER
  viewer << Viewer3D<>::updateDisplay;
  return application.exec();
#else
  return 0;
#endif
}

int main( int argc, char** argv )
{
  // parse command line ----------------------------------------------
  po::options_description general_opt("Allowed options are");
  general_opt.add_options()
      ("help,h", "display this message")
      ("input-file,i", po::value< std::string >(), ".vol file");

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
  bool neededArgsGiven=true;
  if (!(vm.count("input-file"))){
    missingParam("--input-file");
    neededArgsGiven=false;
  }

  const double h = 1.0;
  const double constante = 0.1;

  if(!neededArgsGiven || !parseOK || vm.count("help") || argc <= 1 )
  {
    trace.info()<< "Visualisation of 3d curvature from .vol file using curvature from Integral Invariant" <<std::endl
                << general_opt << "\n"
                << "Basic usage: "<<std::endl
                << "\t3dCurvatureViewer -i <file.vol> --radius <radius> --properties <\"mean\">"<<std::endl
                << std::endl
                << "Below are the different available properties: " << std::endl
                << "\t - \"mean\" for the mean curvature" << std::endl
                << "\t - \"gaussian\" for the Gaussian curvature" << std::endl
                << "\t - \"prindir1\" for the first principal curvature direction" << std::endl
                << "\t - \"prindir2\" for the second principal curvature direction" << std::endl
                << std::endl;
    return 0;
  }

  // Construction of the shape from vol file
  typedef Z3i::Space::RealPoint RealPoint;
  typedef Z3i::Point Point;
  typedef ImageSelector< Z3i::Domain, bool>::Type Image;
  typedef SimpleThresholdForegroundPredicate< Image > ImagePredicate;
  typedef Z3i::KSpace KSpace;
  typedef KSpace::SCell SCell;
  typedef KSpace::Cell Cell;
  typedef KSpace::Surfel KSurfel;
  typedef LightImplicitDigitalSurface< Z3i::KSpace, ImagePredicate > MyLightImplicitDigitalSurface;
  typedef DigitalSurface< MyLightImplicitDigitalSurface > MyDigitalSurface;

  std::string filename = vm["input-file"].as< std::string >();
  Image image = VolReader<Image>::importVol( filename );
  ImagePredicate predicate = ImagePredicate( image, 0 );

  Z3i::Domain domain = image.domain();

  Z3i::KSpace K;

  bool space_ok = K.init( domain.lowerBound(), domain.upperBound(), true );
  if (!space_ok)
  {
    trace.error() << "Error in the Khalimsky space construction."<<std::endl;
    return 2;
  }

  CanonicSCellEmbedder< KSpace > embedder( K );

  SurfelAdjacency< Z3i::KSpace::dimension > SAdj( true );
  KSurfel bel = Surfaces< Z3i::KSpace >::findABel( K, predicate, 100000 );
  MyLightImplicitDigitalSurface LightImplDigSurf( K, predicate, SAdj, bel );
  MyDigitalSurface digSurf( LightImplDigSurf );

  typedef DepthFirstVisitor<MyDigitalSurface> Visitor;
  typedef GraphVisitorRange< Visitor > VisitorRange;
  typedef VisitorRange::ConstIterator SurfelConstIterator;
  typedef ImageToConstantFunctor< Image, ImagePredicate > MyPointFunctor;
  typedef FunctorOnCells< MyPointFunctor, Z3i::KSpace > MyCellFunctor;
  typedef KSpace::Surfel Surfel;

  MyPointFunctor pointFunctor( image, predicate, 1 );
  MyCellFunctor functor ( pointFunctor, K ); // Creation of a functor on Cells, returning true if the cell is inside the shape

  ///////////////////////////////////////////

  trace.beginBlock("Extracting all surfels...");

  std::vector< Surfel > surfels;
  {
    VisitorRange range( new Visitor( digSurf, *digSurf.begin() ) );
    SurfelConstIterator abegin = range.begin();
    SurfelConstIterator aend = range.end();

    while( abegin != aend )
    {
      surfels.push_back( *abegin );
      ++abegin;
    }
  }
  const Dimension size_surfels = surfels.size();

  trace.endBlock();

  ///////////////////////////////////////////

  trace.beginBlock("Analyse segments and Mapping segments <-> Surfels...");

  std::map< Surfel*, std::pair< Statistic< double >, Statistic< double > > > segments;
  for( Dimension ii = 0; ii < size_surfels; ++ii )
  {
    segments[ &surfels[ii] ] = std::pair< Statistic< double >, Statistic< double > >( Statistic< double >( true ), Statistic< double >( true ) );
  }
  computeSegments< MyLightImplicitDigitalSurface, Surfel >( surfels, segments, K, LightImplDigSurf );
  ASSERT(( segments.size() == size_surfels ));

  trace.endBlock();

  ///////////////////////////////////////////

  trace.beginBlock("Computation of radius...");

  std::vector< double > radius( size_surfels, 0.0 );
  computeRadius< Surfel >( surfels, segments, constante, h, radius );

  trace.endBlock();

  ///////////////////////////////////////////

  trace.beginBlock("Curvature computation...");

  typedef double Quantity;
  std::vector< Quantity > curvatures( size_surfels, Quantity(0) );
//  computeCurvature< double, MyCellFunctor, Surfel >( K, functor, surfels, radius, curvatures );

  trace.endBlock();


  /*{
    typedef MyLightImplicitDigitalSurface::Tracker Tracker;
    typedef DigitalSurface2DSlice< Tracker > MySlice;
    typedef Tracker::Surfel Surfel;

    VisitorRange range( new Visitor( digSurf, *digSurf.begin() ) );
    SurfelConstIterator abegin = range.begin();
    SurfelConstIterator aend = range.end();

    //    ImageContainerBySTLMap< Z3i::Domain, std::pair< Statistic< double >, Statistic< double > > > mapStats(  );
    typedef std::pair< Statistic< double >, Statistic< double > > PairOfStatistics;
    typedef std::map< Surfel, PairOfStatistics > SurfelMap;
    SurfelMap mapStat;

    while( abegin != aend )
    {
      mapStat.insert( std::pair< Surfel, PairOfStatistics >( *abegin, PairOfStatistics() ));
      ++abegin;
    }

    typedef std::map< Surfel, bool > MarqueMap;

    for( Dimension dim = 0; dim < 3; ++dim )
    {
      MarqueMap marque;

      for( SurfelMap::const_iterator mbegin = mapStat.begin(), mend = mapStat.end(); mbegin != mend; ++mbegin )
      {
        Surfel currentSurfel = mbegin->first;
        marque.insert( std::pair< Surfel, bool >( currentSurfel, false ));
      }

      for( SurfelMap::const_iterator mbegin = mapStat.begin(), mend = mapStat.end(); mbegin != mend; ++mbegin )
      {
        Surfel currentSurfel = mbegin->first;
        if( marque[ currentSurfel ] == true )
        {
          continue;
        }
        else if( K.sOrthDir( currentSurfel ) == dim )
        {
          continue;
        }
        else
        {
          Tracker ptrTracker( LightImplDigSurf, currentSurfel );
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
//          unsigned int bitSuppr = std::pow(2, dim) + std::pow(2, K.sOrthDir( currentSurfel ));
//          unsigned int bitResu = 7 ^ bitSuppr;

//          bitResu >> (sizeof(unsigned int) * CHAR_BIT - 1);

//          unsigned int aaa = 0;
//          unsigned int bbb = 1;
//          unsigned int ccc = 2;
//          unsigned int ddd = 4;

//          aaa = aaa >> (sizeof(unsigned int) * CHAR_BIT - 1);
//          bbb = bbb >> (sizeof(unsigned int) * CHAR_BIT - 1);
//          ccc = ccc >> (sizeof(unsigned int) * CHAR_BIT - 1);
//          ddd = ddd >> (sizeof(unsigned int) * CHAR_BIT - 1);

//          std::cout << aaa << " " << bbb << " " << ccc << " " << ddd << std::endl;

//          std::cout << "dim: " << dim << " --- orth: " << K.sOrthDir( currentSurfel ) << " --- result " << otherdim << std::endl;

          //dimSlice
          MySlice slice( &ptrTracker, otherdim );

          typedef MySlice::ConstIterator ConstIterator3D;
          typedef SCellProjector< KhalimskySpaceND<2, int> > Functor;
          typedef SCellToPoint< Functor::KSpace > Functor2;
          typedef ConstIteratorAdapter< ConstIterator3D, Functor, Functor::SCell > ConstIterator2D;
          typedef ConstIteratorAdapter< ConstIterator2D, Functor2, Functor2::Output > ConstIterator2DP;

          ConstIterator3D a = slice.begin();
          ConstIterator3D b = ++(slice.begin());
          Dimension dimm = 0;
          while( a->myCoordinates[dimm] != b->myCoordinates[dimm] )
          {
            ++dimm;
          }

          Functor projector;
          projector.initRemoveOneDim( dimm );
          ConstIterator2D xbegin( slice.begin(), projector );
          ConstIterator2D xend( slice.end(), projector );

          Functor::KSpace k2d;
          Functor2 pointFunctor( k2d );

          ConstIterator2DP pbegin( xbegin, pointFunctor );
          ConstIterator2DP pend( xend, pointFunctor );

          const Dimension pr2size = slice.size();
          std::vector< Statistic< double > > v_statMSEL(pr2size);
          for(Dimension ii = 0; ii < pr2size; ++ii )
          {
            v_statMSEL[ii] = Statistic<double>(true);
          }

          analyseAllLengthMS<Functor::KSpace, ConstIterator2DP>( v_statMSEL, pbegin, pend );

          for(Dimension ii = 0; ii < pr2size; ++ii )
          {
            v_statMSEL[ii].terminate();
          }

          Dimension iii = 0;
          for(ConstIterator3D sit = slice.begin(), send = slice.end(); sit != send; ++sit )
          {
            Surfel other = *sit;
            if( marque[ other ] == false )
            {
              marque[ other ] = true;
              PairOfStatistics otherpair = mapStat[ other ];
              if( otherpair.first.samples() == 0 )
              {
                otherpair.first = v_statMSEL[ iii ];
                mapStat[ other ] = otherpair;
              }
              else if ( otherpair.second.samples() == 0 )
              {
                otherpair.second = v_statMSEL[ iii ];
                mapStat[ other ] = otherpair;
              }
              else
              {
                trace.error() << "ALREADY FILLED" << std::endl;
              }
            }
            ++iii;
          }
        }
      }
    }

    VisitorRange range2( new Visitor( digSurf, *digSurf.begin() ) );
    abegin = range2.begin();

    typedef IntegralInvariantGaussianCurvatureEstimator< Z3i::KSpace, MyCellFunctor > MyIIGaussianEstimator;
    MyIIGaussianEstimator estimator ( K, functor );

    double k = 0.1;
    typedef double Quantity;
    typedef EigenValues3D< Quantity >::Matrix33 Matrix3x3;
    typedef EigenValues3D< Quantity >::Vector3 Vector3;
    typedef CurvatureInformations CurvInformation;
    std::vector< Quantity > results;
//    back_insert_iterator< std::vector< CurvInformation > > resultsIterator( results );

    while( abegin != aend )
    {
      Statistic<double> stata = mapStat[ *abegin ].first;
      Statistic<double> statb = mapStat[ *abegin ].second;

      double median;
      if( stata.max() > statb.max() )
      {
        median = stata.median();
      }
      else
      {
        median = statb.median();
      }
//      double mean = std::max( stata.mean(), statb.mean() );
      double re = (k * (median * median)) * 1.0;
      std::cout << re << std::endl;
      if( re < 10.0 )
      {
        re = 10.0;
//        trace.error() << stata.samples() << " ----- " << statb.samples() << " ----- " << stata.median() << " ----- " << statb.median() << std::endl;
      }
      if( lastre != re )
      {
        estimator.init ( h, re ); // Initialisation for a given Euclidean radius of the convolution kernel
      }

      Quantity result = estimator.eval ( abegin );
      results.push_back( result ); // Computation
//      trace.info() << "curvature: " << result << std::endl;
      ++abegin;
    }



//    Tracker* ptrTracker = new Tracker( LightImplDigSurf, oneSCell ); // some pointer on a tracker.
//    MySlice slicex( ptrTracker, 0 ); // slice containing x-axis
//    MySlice slicey( ptrTracker, 1 ); // slice containing y-axis
//    MySlice slicez( ptrTracker, 2 ); // slice containing z-axis

//    typedef MySlice::ConstIterator ConstIterator3D;
//    typedef SCellProjector< KhalimskySpaceND<2, int> > Functor;
//    typedef ConstIteratorAdapter< ConstIterator3D, Functor, Functor::SCell > ConstIterator2D;

//    //////////////////// X
//    ConstIterator3D a = slicex.begin();
//    ConstIterator3D b = ++(slicex.begin());
//    Dimension dim = 0;
//    while( a->myCoordinates[dim] != b->myCoordinates[dim] )
//    {
//      ++dim;
//    }

//    Functor projectx;
//    projectx.initRemoveOneDim( dim );
//    ConstIterator2D xbegin( slicex.begin(), projectx );
//    ConstIterator2D xend( slicex.end(), projectx );

//    //////////////////// Y
//    a = slicey.begin();
//    b = ++(slicey.begin());
//    dim = 0;
//    while( a->myCoordinates[dim] != b->myCoordinates[dim] )
//    {
//      ++dim;
//    }

//    Functor projecty;
//    projecty.initRemoveOneDim( dim );
//    ConstIterator2D ybegin( slicey.begin(), projecty );
//    ConstIterator2D yend( slicey.end(), projecty );

//    //////////////////// Z
//    a = slicez.begin();
//    b = ++(slicez.begin());
//    dim = 0;
//    while( a->myCoordinates[dim] != b->myCoordinates[dim] )
//    {
//      ++dim;
//    }

//    Functor projectz;
//    projectz.initRemoveOneDim( dim );
//    ConstIterator2D zbegin( slicez.begin(), projectz );
//    ConstIterator2D zend( slicez.end(), projectz );

//    Board2D boardx;
//    Board2D boardy;
//    Board2D boardz;

//#ifdef WITH_VISU3D_QGLVIEWER
//    for( ConstIterator3D itb = slicex.begin(); itb != slicex.end(); ++itb )
//    {
//      viewer << CustomColors3D( Color::Green, Color::Green )
//             << *itb;
//    }
//#endif
//    for( ConstIterator2D itb = xbegin; itb != xend; ++itb )
//    {
//      boardx << *itb;
//    }
//#ifdef WITH_VISU3D_QGLVIEWER
//    for( ConstIterator3D itb = slicey.begin(); itb != slicey.end(); ++itb )
//    {
//      viewer << CustomColors3D( Color::Red, Color::Red )
//             << *itb;
//    }
//#endif
//    for( ConstIterator2D itb = ybegin; itb != yend; ++itb )
//    {
//      boardy << *itb;
//    }
//#ifdef WITH_VISU3D_QGLVIEWER
//    for( ConstIterator3D itb = slicez.begin(); itb != slicez.end(); ++itb )
//    {
//      viewer << CustomColors3D( Color::Blue, Color::Blue )
//             << *itb;
//    }
//#endif
//    for( ConstIterator2D itb = zbegin; itb != zend; ++itb )
//    {
//      boardz << *itb;
//    }

//    boardx.saveSVG ( "x.svg" );
//    boardy.saveSVG ( "y.svg" );
//    boardz.saveSVG ( "z.svg" );


  }
  //trace.endBlock();*/

#ifdef WITH_VISU3D_QGLVIEWER
  ///////////////////////////////////////////

  trace.beginBlock("Viewing curvature...");

  return viewCurvature< Surfel, Quantity >( surfels, radius, K, argc, argv );

  trace.endBlock();
#else
  return 0;
#endif
}

///////////////////////////////////////////////////////////////////////////////
