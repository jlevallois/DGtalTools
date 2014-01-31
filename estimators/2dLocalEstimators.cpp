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
 * @file 2dLocalEstimators.cpp
 * @ingroup Tools
 * @author Jacques-Olivier Lachaud (\c jacques-olivier.lachaud@univ-savoie.fr )
 * Laboratory of Mathematics (CNRS, UMR 5807), University of Savoie, France
 * @author David Coeurjolly (\c david.coeurjolly@liris.cnrs.fr)
 * LIRIS (CNRS, UMR 5205),
 * @author Tristan Roussillon (\c tristan.roussillon@liris.cnrs.fr )
 * Laboratoire d'InfoRmatique en Image et Systèmes d'information - LIRIS (CNRS, UMR 5205), CNRS,
 * France
 * @author Jeremy Levallois (\c jeremy.levallois@liris.cnrs.fr )
 * Laboratoire d'InfoRmatique en Image et Systèmes d'information - LIRIS (CNRS, UMR 5205), Universite de Lyon, France
 * LAboratoire de MAthematiques - LAMA (CNRS, UMR 5807), Universite de Savoie, France
 *
 * @date 2011/07/04
 *
 * DGtal tangeant & curvature estimators comparator.
 * @WARNING IntegralInvariant curvature results are set in the reverse order in file. You need to reverse the order in order to compare with others.
 *
 * This file is part of the DGtal library.
 */

///////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <iomanip>
#include <vector>
#include <string>

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

#include "DGtal/base/Common.h"
#include "DGtal/base/Clock.h"

#include "DGtal/math/KMeans.h"

//shapes
#include "DGtal/shapes/ShapeFactory.h"
#include "DGtal/shapes/Shapes.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/topology/helpers/Surfaces.h"
#include "DGtal/topology/DigitalSurface.h"

//Digitizer
#include "DGtal/shapes/GaussDigitizer.h"
#include "DGtal/geometry/curves/GridCurve.h"
#include "DGtal/topology/LightImplicitDigitalSurface.h"
#include "DGtal/graph/DepthFirstVisitor.h"
#include "DGtal/graph/GraphVisitorRange.h"
#include "DGtal/geometry/volumes/KanungoNoise.h"

#include "DGtal/math/Statistic.h"
#include "DGtal/geometry/curves/ArithmeticalDSS.h"

//Estimators
#include "DGtal/geometry/curves/estimation/TrueLocalEstimatorOnPoints.h"
#include "DGtal/geometry/curves/estimation/TrueGlobalEstimatorOnPoints.h"
#include "DGtal/geometry/curves/estimation/ParametricShapeCurvatureFunctor.h"
#include "DGtal/geometry/curves/estimation/ParametricShapeTangentFunctor.h"
#include "DGtal/geometry/curves/estimation/ParametricShapeArcLengthFunctor.h"

#include "DGtal/geometry/curves/BinomialConvolver.h"
#include "DGtal/geometry/curves/estimation/MostCenteredMaximalSegmentEstimator.h"
#include "DGtal/geometry/curves/estimation/SegmentComputerEstimators.h"
#include "DGtal/geometry/curves/ArithmeticalDSSComputer.h"
#include "DGtal/geometry/curves/StabbingCircleComputer.h"

#include "DGtal/images/ImageHelper.h"
#include "DGtal/geometry/surfaces/FunctorOnCells.h"
#include "DGtal/geometry/surfaces/estimation/IntegralInvariantMeanCurvatureEstimator.h"

#include "DGtal/kernel/BasicPointFunctors.h"

using namespace DGtal;


/**
 * Global vectors to describe the available shapes and their
 * parameters.
 */
std::vector<std::string> shapes2D;
std::vector<std::string> shapesDesc;
std::vector<std::string> shapesParam1;
std::vector<std::string> shapesParam2;
std::vector<std::string> shapesParam3;
std::vector<std::string> shapesParam4;

template< typename RealPoint >
struct OptionsIntegralInvariant
{
  double alpha; // <! Alpha parameter for the convolution kernel. 1/3 by default
  double radius; // <! Radius of the convolution kernel.
  RealPoint center; // <! Center of the shape.
  bool lambda_optimized;
  double cste;
  std::string testsII;
};


/**
 * Create the static list of shapes.
 *
 */
void createList()
{
  shapes2D.push_back("ball");
  shapesDesc.push_back("Ball for the Euclidean metric.");
  shapesParam1.push_back("--radius [-R]");
  shapesParam2.push_back("");
  shapesParam3.push_back("");
  shapesParam4.push_back("");

  shapes2D.push_back("square");
  shapesDesc.push_back("square (no signature).");
  shapesParam1.push_back("--width [-w]");
  shapesParam2.push_back("");
  shapesParam3.push_back("");
  shapesParam4.push_back("");

  shapes2D.push_back("lpball");
  shapesDesc.push_back("Ball for the l_power metric (no signature).");
  shapesParam1.push_back("--radius [-R],");
  shapesParam2.push_back("--power [-p]");
  shapesParam3.push_back("");
  shapesParam4.push_back("");

  shapes2D.push_back("flower");
  shapesDesc.push_back("Flower with k petals with radius ranging from R+/-v.");
  shapesParam1.push_back("--radius [-R],");
  shapesParam2.push_back("--varsmallradius [-v],");
  shapesParam3.push_back("--k [-k],");
  shapesParam4.push_back("--phi");

  shapes2D.push_back("ngon");
  shapesDesc.push_back("Regular k-gon.");
  shapesParam1.push_back("--radius [-R],");
  shapesParam2.push_back("--k [-k],");
  shapesParam3.push_back("--phi");
  shapesParam4.push_back("");

  shapes2D.push_back("accflower");
  shapesDesc.push_back("Accelerated Flower with k petals.");
  shapesParam1.push_back("--radius [-R],");
  shapesParam2.push_back("--varsmallradius [-v],");
  shapesParam3.push_back("--k [-k],");
  shapesParam4.push_back("--phi");

  shapes2D.push_back("ellipse");
  shapesDesc.push_back("Ellipse.");
  shapesParam1.push_back("--axis1 [-A],");
  shapesParam2.push_back("--axis2 [-a],");
  shapesParam3.push_back("--phi");
  shapesParam4.push_back("");


}

/**
 * Display the shape list with parameters.
 *
 */
void displayList()
{
  trace.emphase()<<"2D Shapes:"<<std::endl;
  for(unsigned int i=0; i<shapes2D.size(); ++i)
    trace.info()<<"\t"<<shapes2D[i]<<"\t"
               <<shapesDesc[i]<<std::endl
              <<"\t\tRequired parameter(s): "
             << shapesParam1[i]<<" "
             << shapesParam2[i]<<" "
             << shapesParam3[i]<<" "
             << shapesParam4[i]<<std::endl;

}


/**
 * Check if a given shape is available. If not, we exit with an error.
 * If it is, we return the corresponding index in the global vectors.
 *
 * @param shapeName name of the shape to search.
 *
 * @return index of the shape in the shape vectors.
 */
unsigned int checkAndReturnIndex(const std::string &shapeName)
{
  unsigned int pos=0;

  while ((pos < shapes2D.size()) && (shapes2D[pos] != shapeName))
    pos++;

  if (pos == shapes2D.size())
  {
    trace.error() << "The specified shape has not found.";
    trace.info() << std::endl;
    exit(1);
  }

  return pos;
}

/**
 * Missing parameter error message.
 *
 * @param param
 */
void missingParam(std::string param)
{
  trace.error() <<" Parameter: "<<param<<" is required.";
  trace.info()<<std::endl;
  exit(1);
}

/**
 * Estimation error message.
 *
 * @param currentSize number of values returned by the estimator
 * @param expectedSize expected number of values
 */
void estimationError(int currentSize, int expectedSize)
{
  if (currentSize != expectedSize)
  {
    trace.error() << " error in the estimation"
                  << " (got " << currentSize << " values"
                  << " instead of " << expectedSize << ")";
    trace.info() << std::endl;
    exit(1);
  }

}

template <typename KSpace, typename Iterator>
void analyseAllLengthMS( std::vector< Statistic<double> > & statE,
                         Iterator itb, Iterator ite )
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

template <typename KSpace, typename Iterator>
void analyseLengthMS( /*Statistic<double> & statD,*/ Statistic<double> & statE,
                      Iterator itb, Iterator ite )
{
  typedef typename KSpace::Space Space;
  typedef typename Space::Point Point;
  typedef typename Space::Vector Vector;
  typedef ArithmeticalDSSComputer< Iterator, int, 4 > SegmentComputer;
  typedef SaturatedSegmentation< SegmentComputer > Decomposition;
  typedef typename Decomposition::SegmentComputerIterator SegmentComputerIterator;
  // Computes the tangential cover
  SegmentComputer algo;
  Decomposition theDecomposition( itb, ite, algo);
  //statD.clear();
  statE.clear();
  for ( SegmentComputerIterator scIt = theDecomposition.begin(), scItEnd = theDecomposition.end();
        scIt != scItEnd; ++scIt )
  {
    const SegmentComputer & sc = *scIt;
    /*int64_t l = 0;
      for ( Iterator ptIt = sc.begin(), ptItEnd = sc.end(); ptIt != ptItEnd; ++ptIt )
        ++l;
      statD.addValue( (double) l );*/
    Vector v = *( sc.end() - 1 ) - *( sc.begin() );
    statE.addValue( v.norm() );
  }
}

double kappaGridStep( double length )
{
  return 5.0*5.0*5.0/(length*length*length);
}

/**
 * Estimation. Merely call the init and eval methods of the
 * given estimator.
 *
 * @param estimator any local estimator
 * @param h the grid step
 * @param itb begin iterator
 * @param ite end iterator
 * @param ito output iterator on estimated quantities
 */
template <typename Estimator, typename ConstIterator, typename OutputIterator>
void
estimation( Estimator & estimator, double h,
            const ConstIterator& itb, const ConstIterator& ite, const OutputIterator& ito )
{
  Clock c;
  c.startClock();
  estimator.init( h, itb, ite );
  estimator.eval( itb, ite, ito );
  double time = c.stopClock();
  std::cout << "# Time: " << time << std::endl;
}


/**
 *
 * @return Euclidean radius for the convolver of Integral Invariant estimators
 */
template< typename ConstIteratorOnPoints, typename Point >
unsigned int suggestedSizeIntegralInvariant( const double h,
                                             const Point& center,
                                             const ConstIteratorOnPoints& itb,
                                             const ConstIteratorOnPoints& ite )
{
  typedef typename Point::Component TValue;

  ConstIteratorOnPoints it = itb;
  Point p( *it );
  Point distance = p - center;
  TValue minRadius = distance.norm();
  ++it;

  for ( ; it != ite; ++it )
  {
    p = *it;
    distance = p - center;
    if ( distance.norm() < minRadius )
    {
      minRadius = distance.norm();
    }
  }

  return minRadius * h;
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

/**
 * Estimation of tangents and curvature
 * from several different methods
 *
 * @param filename name of a file to save results ( will be postfix by name of estimators )
 * @param aShape shape
 * @param h grid step
 * @param optionsII options for Integral Invariants estimators
 * @param options estimators to use (BC, II, MDCA, ...)
 * @param properties properties of estimators (curvature and/or tangeant)
 * @param noiseLevel level to noised the shape. 0 <= noiseLevel < 1
 */
template <typename Space, typename Shape>
bool
computeLocalEstimations( const std::string & filename,
                         Shape * aShape,
                         const double & h,
                         struct OptionsIntegralInvariant< Z2i::RealPoint > optionsII,
                         const std::string & options,
                         const std::string & properties,
                         double noiseLevel = 0.0 )
{
  // Types
  typedef typename Space::Point Point;
  typedef typename Space::Vector Vector;
  typedef typename Space::RealPoint RealPoint;
  typedef typename Space::Integer Integer;
  typedef HyperRectDomain<Space> Domain;
  typedef KhalimskySpaceND<Space::dimension,Integer> KSpace;
  typedef typename KSpace::SCell SCell;
  typedef GaussDigitizer<Space,Shape> Digitizer;
  typedef KanungoNoise< Digitizer, Z2i::Domain > KanungoPredicate;

  bool withNoise = ( noiseLevel <= 0.0 ) ? false : true;
  /*if( withNoise )
        noiseLevel = std::pow(noiseLevel, h);*/

  ASSERT (( noiseLevel < 1.0 ));

  bool tangent = ( properties.at( 0 ) != '0' ) ? true : false;
  bool curvature = ( properties.at( 1 ) != '0' ) ? true : false;

  // Digitizer
  Digitizer* dig = new Digitizer();
  dig->attach( *aShape ); // attaches the shape.
  Vector vlow(-1,-1); Vector vup(1,1);
  dig->init( aShape->getLowerBound()+vlow, aShape->getUpperBound()+vup, h );
  Domain domain = dig->getDomain();

  //Noise

  Clock c;

  // Create cellular space
  KSpace K;
  bool ok = K.init( dig->getLowerBound(), dig->getUpperBound(), true );
  if ( ! ok )
  {
    std::cerr << "[2dLocalEstimators]"
              << " error in creating KSpace." << std::endl;
    return false;
  }
  try {

    // Extracts shape boundary
    SurfelAdjacency< KSpace::dimension > SAdj( true );
    SCell bel;
    std::vector< SCell > points;

    KanungoPredicate  *noisifiedObject;
    if ( withNoise )
    {
      noisifiedObject = new KanungoPredicate( *dig, domain, noiseLevel );
      bel = Surfaces< KSpace >::findABel( K, *noisifiedObject, 10000 );
      Surfaces< KSpace >::track2DBoundary( points, K, SAdj, *noisifiedObject, bel );

      double minsize = dig->getUpperBound()[0] - dig->getLowerBound()[0];
      while( points.size() < 2 * minsize )
      {
        points.clear();
        bel = Surfaces< KSpace >::findABel( K, *noisifiedObject, 10000 );
        Surfaces< KSpace >::track2DBoundary( points, K, SAdj, *noisifiedObject, bel );
      }
    }
    else
    {
      bel = Surfaces< KSpace >::findABel( K, *dig, 10000 );
      Surfaces< KSpace >::track2DBoundary( points, K, SAdj, *dig, bel );
    }

    // Create GridCurve
    GridCurve< KSpace > gridcurve;
    gridcurve.initFromSCellsVector( points );

    // Ranges
    typedef typename GridCurve< KSpace >::MidPointsRange MidPointsRange;
    MidPointsRange pointsRange = gridcurve.getMidPointsRange();
    typedef typename GridCurve< KSpace >::PointsRange PointsRange;
    PointsRange pointsRange2 = gridcurve.getPointsRange();

    // Estimations
    if (gridcurve.isClosed())
    {
      //////////////////////////
      //Statistic<double> statMSL( true );
      Statistic<double> statMSEL( true );
      analyseLengthMS<KSpace>( /*statMSL,*/ statMSEL, pointsRange2.begin(), pointsRange2.end() );
      //statMSL.terminate();
      /*Statistic<double> resolution( false );
        for ( Statistic<double>::ConstIterator
              it = statMSL.begin(), itE = statMSL.end(); it != itE; ++it )
          resolution.addValue( 1.0 / kappaGridStep( *it ) );
        Statistic<double> resolutionE( false );
        for ( Statistic<double>::ConstIterator
              it = statMSEL.begin(), itE = statMSEL.end(); it != itE; ++it )
          resolutionE.addValue( 1.0 / kappaGridStep( *it ) );*/

      trace.info() << "#Average " << statMSEL.mean() << std::endl;
      //////////////////////////

      //////////////////////////
      const Dimension pr2size = (pointsRange2.size());
      std::vector< Statistic< double > > v_statMSEL(pr2size);
      for(Dimension ii = 0; ii < pr2size; ++ii )
      {
        v_statMSEL[ii] = Statistic<double>(true);
      }

      trace.beginBlock("Analyse segments and Mapping segments <-> Surfels...");

      analyseAllLengthMS<KSpace>( v_statMSEL, pointsRange2.begin(), pointsRange2.end() );

      trace.endBlock();

      trace.info() << "done" << std::endl;
      //return 1;
      //////////////////////////

      if (options.at(0) != '0')
      {
        if( tangent )
        {
          char full_filename[360];
          sprintf( full_filename, "%s%s", filename.c_str(), "_True_tangeant.dat" );
          std::ofstream file( full_filename );

          file << "# h = " << h << std::endl;
          file << "# True tangents computation" << std::endl;
          file << "# range size = " << pointsRange.size() << std::endl;
          if ( withNoise )
          {
            file << "# noise level (init) = " << noiseLevel/h << std::endl;
            file << "# noise level (current) = " << noiseLevel << std::endl;
          }

          std::ostream_iterator< RealPoint > out_it( file, "\n" );

          typedef ParametricShapeTangentFunctor< Shape > TangentFunctor;
          typedef typename MidPointsRange::ConstCirculator C;
          TrueLocalEstimatorOnPoints< C, Shape, TangentFunctor >
              trueTangentEstimator;
          trueTangentEstimator.attach( aShape );
          estimation( trueTangentEstimator, h,
                      pointsRange.c(), pointsRange.c(),
                      out_it );

          file.close();

        }

        if( curvature )
        {
          char full_filename[360];
          sprintf( full_filename, "%s%s", filename.c_str(), "_True_curvature.dat" );
          std::ofstream file( full_filename );

          file << "# h = " << h << std::endl;
          file << "# True curvature computation" << std::endl;
          file << "# range size = " << pointsRange.size() << std::endl;
          if ( withNoise )
          {
            file << "# noise level (init) = " << noiseLevel/h << std::endl;
            file << "# noise level (current) = " << noiseLevel << std::endl;
          }

          std::ostream_iterator< double > out_it( file, "\n" );

          typedef ParametricShapeCurvatureFunctor< Shape > CurvatureFunctor;
          typedef typename MidPointsRange::ConstCirculator C;
          TrueLocalEstimatorOnPoints< C, Shape, CurvatureFunctor >
              trueCurvatureEstimator;
          trueCurvatureEstimator.attach( aShape );
          estimation( trueCurvatureEstimator, h,
                      pointsRange.c(), pointsRange.c(),
                      out_it );

          file.close();
        }
      }

      // Maximal Segments
      if (options.at(1) != '0')
      {
        if( tangent )
        {
          char full_filename[360];
          sprintf( full_filename, "%s%s", filename.c_str(), "_MDSS_tangeant.dat" );
          std::ofstream file( full_filename );

          file << "# h = " << h << std::endl;
          file << "# Most centered maximal DSS tangent estimation" << std::endl;
          file << "# range size = " << pointsRange.size() << std::endl;
          if ( withNoise )
          {
            file << "# noise level (init) = " << noiseLevel/h << std::endl;
            file << "# noise level (current) = " << noiseLevel << std::endl;
          }

          std::ostream_iterator< RealPoint > out_it( file, "\n" );

          typedef typename GridCurve< KSpace >::PointsRange PointsRange2;
          PointsRange2 pointsRange2 = gridcurve.getPointsRange();

          typedef typename PointsRange2::ConstCirculator C;
          typedef ArithmeticalDSSComputer< C, int, 4 > SegmentComputer;
          typedef TangentFromDSSEstimator<SegmentComputer> SCFunctor;
          SegmentComputer sc;
          SCFunctor f;


          MostCenteredMaximalSegmentEstimator<SegmentComputer,SCFunctor> MDSSTangentEstimator(sc, f);
          estimation( MDSSTangentEstimator, h,
                      pointsRange2.c(), pointsRange2.c(),
                      out_it );

          file.close();
        }
        if( curvature )
        {
          c.startClock();

          char full_filename[360];
          sprintf( full_filename, "%s%s", filename.c_str(), "_MDSSl_curvature.dat" );
          std::ofstream file( full_filename );

          file << "# h = " << h << std::endl;
          file << "# Most centered maximal DSS (length) curvature estimation" << std::endl;
          file << "# range size = " << pointsRange.size() << std::endl;
          if ( withNoise )
          {
            file << "# noise level (init) = " << noiseLevel/h << std::endl;
            file << "# noise level (current) = " << noiseLevel << std::endl;
          }

          std::ostream_iterator< double > out_it( file, "\n" );

          typedef typename GridCurve< KSpace >::PointsRange PointsRange2;
          PointsRange2 pointsRange2 = gridcurve.getPointsRange();

          typedef typename PointsRange2::ConstCirculator C;
          typedef ArithmeticalDSSComputer< C, int, 4 > SegmentComputer;
          typedef CurvatureFromDSSLengthEstimator<SegmentComputer> SCFunctor;
          SegmentComputer sc;
          SCFunctor f;
          MostCenteredMaximalSegmentEstimator<SegmentComputer,SCFunctor> MDSSCurvatureEstimator(sc, f);

          estimation( MDSSCurvatureEstimator, h,
                      pointsRange2.c(), pointsRange2.c(),
                      out_it );

          file.close();


          memset(&full_filename[0], 0, sizeof(full_filename));
          sprintf( full_filename, "%s%s", filename.c_str(), "_MDSSlw_curvature.dat" );
          file.open( full_filename );

          file << "# h = " << h << std::endl;
          file << "# Most centered maximal DSS (length & width) curvature estimation" << std::endl;
          file << "# range size = " << pointsRange.size() << std::endl;
          if ( withNoise )
          {
            file << "# noise level (init) = " << noiseLevel/h << std::endl;
            file << "# noise level (current) = " << noiseLevel << std::endl;
          }

          std::ostream_iterator< double > out_it2( file, "\n" );

          typedef CurvatureFromDSSEstimator<SegmentComputer> SCFunctor2;
          SegmentComputer sc2;
          SCFunctor2 f2;
          MostCenteredMaximalSegmentEstimator<SegmentComputer,SCFunctor2> MDSSCurvatureEstimator2(sc2, f2);
          estimation( MDSSCurvatureEstimator2, h,
                      pointsRange2.c(), pointsRange2.c(),
                      out_it2 );

          double time = c.stopClock();
          file << "# Time: " << time << std::endl;

          file.close();

        }
      }

      //Maximal circular arcs
      if (options.at(2) != '0')
      {
        if( tangent )
        {
          char full_filename[360];
          sprintf( full_filename, "%s%s", filename.c_str(), "_MDCA_tangent.dat" );
          std::ofstream file( full_filename );

          file << "# h = " << h << std::endl;
          file << "# Most centered maximal DCA tangents estimation" << std::endl;
          file << "# range size = " << pointsRange.size() << std::endl;
          if ( withNoise )
          {
            file << "# noise level (init) = " << noiseLevel/h << std::endl;
            file << "# noise level (current) = " << noiseLevel << std::endl;
          }

          std::ostream_iterator< RealPoint > out_it( file, "\n" );

          typedef typename GridCurve<KSpace>::IncidentPointsRange Range;
          typedef typename Range::ConstCirculator C;
          Range r = gridcurve.getIncidentPointsRange();
          typedef StabbingCircleComputer<C> SegmentComputer;
          typedef TangentFromDCAEstimator<SegmentComputer> SCFunctor;
          SegmentComputer sc;
          SCFunctor f;
          MostCenteredMaximalSegmentEstimator<SegmentComputer,SCFunctor> MDCATangentEstimator(sc, f);
          estimation( MDCATangentEstimator, h,
                      r.c(), r.c(),
                      out_it );

          file.close();
        }

        if( curvature )
        {
          c.startClock();

          char full_filename[360];
          sprintf( full_filename, "%s%s", filename.c_str(), "_MDCA_curvature.dat" );
          std::ofstream file( full_filename );

          file << "# h = " << h << std::endl;
          file << "# Most centered maximal DCA curvature estimation" << std::endl;
          file << "# range size = " << pointsRange.size() << std::endl;
          if ( withNoise )
          {
            file << "# noise level (init) = " << noiseLevel/h << std::endl;
            file << "# noise level (current) = " << noiseLevel << std::endl;
          }

          std::ostream_iterator< double > out_it( file, "\n" );

          typedef typename GridCurve<KSpace>::IncidentPointsRange Range;
          typedef typename Range::ConstCirculator C;
          Range r = gridcurve.getIncidentPointsRange();
          typedef StabbingCircleComputer<C> SegmentComputer;
          typedef CurvatureFromDCAEstimator<SegmentComputer, false> SCFunctor;
          SegmentComputer sc;
          SCFunctor f;
          MostCenteredMaximalSegmentEstimator<SegmentComputer,SCFunctor> MDCACurvatureEstimator(sc, f);
          estimation( MDCACurvatureEstimator, h,
                      r.c(), r.c(),
                      out_it );

          double time = c.stopClock();
          file << "# Time: " << time << std::endl;

          file.close();
        }
      }

      //Binomial convolver
      if (options.at(3) != '0')
      {
        if( tangent )
        {
          c.startClock();

          char full_filename[360];
          sprintf( full_filename, "%s%s", filename.c_str(), "_BC_tangeant.dat" );
          std::ofstream file( full_filename );

          file << "# h = " << h << std::endl;
          file << "# Tangents estimation from binomial convolution" << std::endl;
          file << "# range size = " << pointsRange.size() << std::endl;
          if ( withNoise )
          {
            file << "# noise level (init) = " << noiseLevel/h << std::endl;
            file << "# noise level (current) = " << noiseLevel << std::endl;
          }

          typedef typename MidPointsRange::ConstIterator I;
          typedef BinomialConvolver<I, double> MyBinomialConvolver;
          file << "# mask size = " <<
                  MyBinomialConvolver::suggestedSize( h, pointsRange.begin(), pointsRange.end() ) << std::endl;

          typedef TangentFromBinomialConvolverFunctor< MyBinomialConvolver, RealPoint >
              TangentBCFct;
          BinomialConvolverEstimator< MyBinomialConvolver, TangentBCFct> BCTangentEstimator;

          std::ostream_iterator< RealPoint > out_it( file, "\n" );

          BCTangentEstimator.init( h, pointsRange.begin(), pointsRange.end(), true );
          BCTangentEstimator.eval( pointsRange.begin(), pointsRange.end(), out_it );

          double time = c.stopClock();
          file << "# Time: " << time << std::endl;

          file.close();
        }

        if( curvature )
        {
          c.startClock();

          char full_filename[360];
          sprintf( full_filename, "%s%s", filename.c_str(), "_BC_curvature.dat" );
          std::ofstream file( full_filename );

          file << "# h = " << h << std::endl;
          file << "# Curvature estimation from binomial convolution" << std::endl;
          file << "# range size = " << pointsRange.size() << std::endl;
          if ( withNoise )
          {
            file << "# noise level (init) = " << noiseLevel/h << std::endl;
            file << "# noise level (current) = " << noiseLevel << std::endl;
          }

          typedef typename MidPointsRange::ConstIterator I;
          typedef BinomialConvolver<I, double> MyBinomialConvolver;
          file << "# mask size = " <<
                  MyBinomialConvolver::suggestedSize( h, pointsRange.begin(), pointsRange.end() ) << std::endl;

          std::ostream_iterator< double > out_it( file, "\n" );

          typedef CurvatureFromBinomialConvolverFunctor< MyBinomialConvolver, double >
              CurvatureBCFct;
          BinomialConvolverEstimator< MyBinomialConvolver, CurvatureBCFct> BCCurvatureEstimator;

          BCCurvatureEstimator.init( h, pointsRange.begin(), pointsRange.end(), true );
          BCCurvatureEstimator.eval( pointsRange.begin(), pointsRange.end(), out_it );

          double time = c.stopClock();
          file << "# Time: " << time << std::endl;

          file.close();
        }
      }

      /// <! @WARNING IntegralInvariant curvature results are set in the reverse order in file. You need to reverse the order in order to compare with others.
      //Integral Invariants
      if (options.at(4) != '0')
      {
        if( curvature )
        {
          /////// global mean
          if( optionsII.testsII.at(0) != '0' )
          {
            c.startClock();

            char full_filename[360];
            sprintf( full_filename, "%s%s", filename.c_str(), "_II_curvature_global_mean.dat" );
            std::ofstream file( full_filename );

            file << "# h = " << h << std::endl;
            file << "# Integral Invariant curvature estimation" << std::endl;
            file << "# range size = " << pointsRange.size() << std::endl;
            if ( withNoise )
            {
              file << "# noise level (init) = " << noiseLevel/h << std::endl;
              file << "# noise level (current) = " << noiseLevel << std::endl;
            }

            if( optionsII.radius <= 0.0 )
            {
              optionsII.radius = suggestedSizeIntegralInvariant( h, dig->round( optionsII.center ), pointsRange.begin(), pointsRange.end() );
              file << "# Estimated radius: " << optionsII.radius << std::endl;
            }

            double k = optionsII.cste;
            double mean = statMSEL.mean();
            double re_convolution_kernel = (k * (mean * mean)) * h;
            //          double re_convolution_kernel = optionsII.radius * std::pow( h, optionsII.alpha );
            file << "# full kernel (digital) size (with alpha = " << optionsII.alpha << ") = " <<
                    re_convolution_kernel / h << std::endl;

            std::ostream_iterator< double > out_it( file, "\n" );

            if ( withNoise )
            {
              typedef LightImplicitDigitalSurface< KSpace, KanungoPredicate > LightImplicitDigSurface;
              typedef DigitalSurface< LightImplicitDigSurface > DigSurface;

              LightImplicitDigSurface LightImplDigSurf( K, *noisifiedObject, SAdj, bel );
              DigSurface surf( LightImplDigSurf );

              typedef PointFunctorFromPointPredicateAndDomain< KanungoPredicate, Domain, unsigned int > KanungoFunctor;
              KanungoFunctor * noisifiedFunctor = new KanungoFunctor( noisifiedObject, domain, 1, 0 );

              typedef FunctorOnCells< KanungoFunctor, KSpace > CurvatureIIFct;
              CurvatureIIFct * functor = new CurvatureIIFct( *noisifiedFunctor, K );

              IntegralInvariantMeanCurvatureEstimator< KSpace, CurvatureIIFct> * IICurvatureEstimator = new IntegralInvariantMeanCurvatureEstimator< KSpace, CurvatureIIFct>( K, *functor );

              typedef DepthFirstVisitor< DigSurface > Visitor;
              typedef GraphVisitorRange< Visitor > VisitorRange;
              typedef typename VisitorRange::ConstIterator I;

              VisitorRange range( new Visitor( surf, *surf.begin() ) );
              I ibegin = range.begin();
              I iend = range.end();

              IICurvatureEstimator->init( h, re_convolution_kernel );
              IICurvatureEstimator->eval( points.begin(), points.end(), out_it );

              delete functor;
              delete noisifiedFunctor;
              delete IICurvatureEstimator;
            }
            else
            {
              typedef LightImplicitDigitalSurface< KSpace, Digitizer > LightImplicitDigSurface;
              typedef DigitalSurface< LightImplicitDigSurface > DigSurface;

              LightImplicitDigSurface LightImplDigSurf( K, *dig, SAdj, bel );
              DigSurface surf( LightImplDigSurf );

              typedef PointFunctorFromPointPredicateAndDomain< Digitizer, Domain, unsigned int > MyPointFunctor;
              MyPointFunctor * pointFunctor = new MyPointFunctor( dig, domain, 1, 0 );

              typedef FunctorOnCells< MyPointFunctor, KSpace > CurvatureIIFct;
              CurvatureIIFct * functor = new CurvatureIIFct( *pointFunctor, K );

              IntegralInvariantMeanCurvatureEstimator< KSpace, CurvatureIIFct> * IICurvatureEstimator = new IntegralInvariantMeanCurvatureEstimator< KSpace, CurvatureIIFct>( K, *functor );

              typedef DepthFirstVisitor< DigSurface > Visitor;
              typedef GraphVisitorRange< Visitor > VisitorRange;
              typedef typename VisitorRange::ConstIterator I;

              VisitorRange range( new Visitor( surf, *surf.begin() ) );
              I ibegin = range.begin();
              I iend = range.end();

              IICurvatureEstimator->init( h, re_convolution_kernel );
              IICurvatureEstimator->eval( ibegin, iend, out_it );

              delete functor;
              delete pointFunctor;
              delete IICurvatureEstimator;
            }

            double time = c.stopClock();
            file << "# Time: " << time << std::endl;

            file.close();
          }

          /////// global max
          if( optionsII.testsII.at(1) != '0' )
          {
            c.startClock();

            char full_filename[360];
            sprintf( full_filename, "%s%s", filename.c_str(), "_II_curvature_global_max.dat" );
            std::ofstream file( full_filename );

            file << "# h = " << h << std::endl;
            file << "# Integral Invariant curvature estimation" << std::endl;
            file << "# range size = " << pointsRange.size() << std::endl;
            if ( withNoise )
            {
              file << "# noise level (init) = " << noiseLevel/h << std::endl;
              file << "# noise level (current) = " << noiseLevel << std::endl;
            }

            if( optionsII.radius <= 0.0 )
            {
              optionsII.radius = suggestedSizeIntegralInvariant( h, dig->round( optionsII.center ), pointsRange.begin(), pointsRange.end() );
              file << "# Estimated radius: " << optionsII.radius << std::endl;
            }

            double k = optionsII.cste;
            double max = statMSEL.max();
            double re_convolution_kernel = (k * (max * max)) * h;
            //          double re_convolution_kernel = optionsII.radius * std::pow( h, optionsII.alpha );
            file << "# full kernel (digital) size (with alpha = " << optionsII.alpha << ") = " <<
                    re_convolution_kernel / h << std::endl;

            std::ostream_iterator< double > out_it( file, "\n" );

            if ( withNoise )
            {
              typedef LightImplicitDigitalSurface< KSpace, KanungoPredicate > LightImplicitDigSurface;
              typedef DigitalSurface< LightImplicitDigSurface > DigSurface;

              LightImplicitDigSurface LightImplDigSurf( K, *noisifiedObject, SAdj, bel );
              DigSurface surf( LightImplDigSurf );

              typedef PointFunctorFromPointPredicateAndDomain< KanungoPredicate, Domain, unsigned int > KanungoFunctor;
              KanungoFunctor * noisifiedFunctor = new KanungoFunctor( noisifiedObject, domain, 1, 0 );

              typedef FunctorOnCells< KanungoFunctor, KSpace > CurvatureIIFct;
              CurvatureIIFct * functor = new CurvatureIIFct( *noisifiedFunctor, K );

              IntegralInvariantMeanCurvatureEstimator< KSpace, CurvatureIIFct> * IICurvatureEstimator = new IntegralInvariantMeanCurvatureEstimator< KSpace, CurvatureIIFct>( K, *functor );

              typedef DepthFirstVisitor< DigSurface > Visitor;
              typedef GraphVisitorRange< Visitor > VisitorRange;
              typedef typename VisitorRange::ConstIterator I;

              VisitorRange range( new Visitor( surf, *surf.begin() ) );
              I ibegin = range.begin();
              I iend = range.end();

              IICurvatureEstimator->init( h, re_convolution_kernel );
              IICurvatureEstimator->eval( points.begin(), points.end(), out_it );

              delete functor;
              delete noisifiedFunctor;
              delete IICurvatureEstimator;
            }
            else
            {
              typedef LightImplicitDigitalSurface< KSpace, Digitizer > LightImplicitDigSurface;
              typedef DigitalSurface< LightImplicitDigSurface > DigSurface;

              LightImplicitDigSurface LightImplDigSurf( K, *dig, SAdj, bel );
              DigSurface surf( LightImplDigSurf );

              typedef PointFunctorFromPointPredicateAndDomain< Digitizer, Domain, unsigned int > MyPointFunctor;
              MyPointFunctor * pointFunctor = new MyPointFunctor( dig, domain, 1, 0 );

              typedef FunctorOnCells< MyPointFunctor, KSpace > CurvatureIIFct;
              CurvatureIIFct * functor = new CurvatureIIFct( *pointFunctor, K );

              IntegralInvariantMeanCurvatureEstimator< KSpace, CurvatureIIFct> * IICurvatureEstimator = new IntegralInvariantMeanCurvatureEstimator< KSpace, CurvatureIIFct>( K, *functor );

              typedef DepthFirstVisitor< DigSurface > Visitor;
              typedef GraphVisitorRange< Visitor > VisitorRange;
              typedef typename VisitorRange::ConstIterator I;

              VisitorRange range( new Visitor( surf, *surf.begin() ) );
              I ibegin = range.begin();
              I iend = range.end();

              IICurvatureEstimator->init( h, re_convolution_kernel );
              IICurvatureEstimator->eval( ibegin, iend, out_it );

              delete functor;
              delete pointFunctor;
              delete IICurvatureEstimator;
            }

            double time = c.stopClock();
            file << "# Time: " << time << std::endl;

            file.close();
          }

          /////// global median
          if( optionsII.testsII.at(2) != '0' )
          {
            c.startClock();

            char full_filename[360];
            sprintf( full_filename, "%s%s", filename.c_str(), "_II_curvature_global_median.dat" );
            std::ofstream file( full_filename );

            file << "# h = " << h << std::endl;
            file << "# Integral Invariant curvature estimation" << std::endl;
            file << "# range size = " << pointsRange.size() << std::endl;
            if ( withNoise )
            {
              file << "# noise level (init) = " << noiseLevel/h << std::endl;
              file << "# noise level (current) = " << noiseLevel << std::endl;
            }

            if( optionsII.radius <= 0.0 )
            {
              optionsII.radius = suggestedSizeIntegralInvariant( h, dig->round( optionsII.center ), pointsRange.begin(), pointsRange.end() );
              file << "# Estimated radius: " << optionsII.radius << std::endl;
            }

            double k = optionsII.cste;
            double median = statMSEL.median();
            double re_convolution_kernel = (k * (median * median)) * h;
            //          double re_convolution_kernel = optionsII.radius * std::pow( h, optionsII.alpha );
            file << "# full kernel (digital) size (with alpha = " << optionsII.alpha << ") = " <<
                    re_convolution_kernel / h << std::endl;

            std::ostream_iterator< double > out_it( file, "\n" );

            if ( withNoise )
            {
              typedef LightImplicitDigitalSurface< KSpace, KanungoPredicate > LightImplicitDigSurface;
              typedef DigitalSurface< LightImplicitDigSurface > DigSurface;

              LightImplicitDigSurface LightImplDigSurf( K, *noisifiedObject, SAdj, bel );
              DigSurface surf( LightImplDigSurf );

              typedef PointFunctorFromPointPredicateAndDomain< KanungoPredicate, Domain, unsigned int > KanungoFunctor;
              KanungoFunctor * noisifiedFunctor = new KanungoFunctor( noisifiedObject, domain, 1, 0 );

              typedef FunctorOnCells< KanungoFunctor, KSpace > CurvatureIIFct;
              CurvatureIIFct * functor = new CurvatureIIFct( *noisifiedFunctor, K );

              IntegralInvariantMeanCurvatureEstimator< KSpace, CurvatureIIFct> * IICurvatureEstimator = new IntegralInvariantMeanCurvatureEstimator< KSpace, CurvatureIIFct>( K, *functor );

              typedef DepthFirstVisitor< DigSurface > Visitor;
              typedef GraphVisitorRange< Visitor > VisitorRange;
              typedef typename VisitorRange::ConstIterator I;

              VisitorRange range( new Visitor( surf, *surf.begin() ) );
              I ibegin = range.begin();
              I iend = range.end();

              IICurvatureEstimator->init( h, re_convolution_kernel );
              IICurvatureEstimator->eval( points.begin(), points.end(), out_it );

              delete functor;
              delete noisifiedFunctor;
              delete IICurvatureEstimator;
            }
            else
            {
              typedef LightImplicitDigitalSurface< KSpace, Digitizer > LightImplicitDigSurface;
              typedef DigitalSurface< LightImplicitDigSurface > DigSurface;

              LightImplicitDigSurface LightImplDigSurf( K, *dig, SAdj, bel );
              DigSurface surf( LightImplDigSurf );

              typedef PointFunctorFromPointPredicateAndDomain< Digitizer, Domain, unsigned int > MyPointFunctor;
              MyPointFunctor * pointFunctor = new MyPointFunctor( dig, domain, 1, 0 );

              typedef FunctorOnCells< MyPointFunctor, KSpace > CurvatureIIFct;
              CurvatureIIFct * functor = new CurvatureIIFct( *pointFunctor, K );

              IntegralInvariantMeanCurvatureEstimator< KSpace, CurvatureIIFct> * IICurvatureEstimator = new IntegralInvariantMeanCurvatureEstimator< KSpace, CurvatureIIFct>( K, *functor );

              typedef DepthFirstVisitor< DigSurface > Visitor;
              typedef GraphVisitorRange< Visitor > VisitorRange;
              typedef typename VisitorRange::ConstIterator I;

              VisitorRange range( new Visitor( surf, *surf.begin() ) );
              I ibegin = range.begin();
              I iend = range.end();

              IICurvatureEstimator->init( h, re_convolution_kernel );
              IICurvatureEstimator->eval( ibegin, iend, out_it );

              delete functor;
              delete pointFunctor;
              delete IICurvatureEstimator;
            }

            double time = c.stopClock();
            file << "# Time: " << time << std::endl;

            file.close();
          }

          /////// local mean
          if( optionsII.testsII.at(3) != '0' )
          {
            c.startClock();

            char full_filename[360];
            sprintf( full_filename, "%s%s", filename.c_str(), "_II_curvature_local_mean.dat" );
            std::ofstream file( full_filename );

            file << "# h = " << h << std::endl;
            file << "# Integral Invariant curvature estimation" << std::endl;
            file << "# range size = " << pointsRange.size() << std::endl;
            if ( withNoise )
            {
              file << "# noise level (init) = " << noiseLevel/h << std::endl;
              file << "# noise level (current) = " << noiseLevel << std::endl;
            }

            if( optionsII.radius <= 0.0 )
            {
              optionsII.radius = suggestedSizeIntegralInvariant( h, dig->round( optionsII.center ), pointsRange.begin(), pointsRange.end() );
              file << "# Estimated radius: " << optionsII.radius << std::endl;
            }

            double k = optionsII.cste;

            ///////////
            char full_filename2[360];
            sprintf( full_filename2, "%s%s", filename.c_str(), "_II_curvature_local_mean_re.dat" );
            std::ofstream file2( full_filename2 );

            char full_filename3[360];
            sprintf( full_filename3, "%s%s", filename.c_str(), "_II_curvature_local_stats.dat" );
            std::ofstream file3( full_filename3 );

            char full_filename4[360];
            sprintf( full_filename4, "%s%s", filename.c_str(), "_II_curvature_local_mean_re_used.dat" );
            std::ofstream file4( full_filename2 );
            ///////////

            if ( withNoise )
            {
              trace.error() << "not yet" << std::endl;
            }
            else
            {
              typedef LightImplicitDigitalSurface< KSpace, Digitizer > LightImplicitDigSurface;
              typedef DigitalSurface< LightImplicitDigSurface > DigSurface;
              typedef PointFunctorFromPointPredicateAndDomain< Digitizer, Domain, unsigned int > MyPointFunctor;
              typedef FunctorOnCells< MyPointFunctor, KSpace > CurvatureIIFct;
              typedef IntegralInvariantMeanCurvatureEstimator< KSpace, CurvatureIIFct > Estimator;
              typedef DepthFirstVisitor< DigSurface > Visitor;
              typedef GraphVisitorRange< Visitor > VisitorRange;
              typedef typename VisitorRange::ConstIterator ConstIterator;

              LightImplicitDigSurface LightImplDigSurf( K, *dig, SAdj, bel );
              DigSurface surf( LightImplDigSurf );

              MyPointFunctor pointFunctor( dig, domain, 1, 0 );
              CurvatureIIFct functor( pointFunctor, K );

              VisitorRange range( new Visitor( surf, *surf.begin() ) );
              ConstIterator ibegin = range.begin();
              ConstIterator iend = range.end();

              trace.beginBlock("Extracting all surfels...");

              std::vector< SCell > contour;

              for( ; ibegin != iend; ++ibegin )
              {
                contour.push_back( *ibegin );
              }

              if( contour.size() != pr2size )
              {
                trace.error() << "ERROR! Not the same border size: " << contour.size() << " " << pr2size << std::endl;
                trace.endBlock();
                return 0;
              }

              trace.endBlock();

              std::vector< double > v_curvatures( contour.size() );
              std::vector< double > v_estimated_radius( contour.size() );

              trace.beginBlock("Computation of radius...");

              for( Dimension ii = 0; ii < pr2size; ++ii )
              {
                Dimension current_pos = pr2size - 1 - ii;
                double mean = v_statMSEL[ current_pos ].mean();
                v_estimated_radius[ii] = (k * (mean * mean)) * h;
                for( Statistic<double>::ConstIterator itb = v_statMSEL[current_pos].begin(), ite = v_statMSEL[current_pos].end(); itb != ite; ++itb )
                {
                  file3 << *itb << " ";
                }
                file3 << std::endl;
              }

              trace.endBlock();

              trace.beginBlock("Sorting radius & pre-computing estimators...");

              std::vector< double > v_radius;
              std::vector< Dimension > v_registration;

              const Dimension nbOfRadius = 10;

              suggestedRadiusForIntegralInvariantEstimators( v_estimated_radius, v_registration, v_radius, nbOfRadius );

              ASSERT(( v_radius.size() == nbOfRadius ));

              std::vector< Estimator* > v_estimators( nbOfRadius );

              for( Dimension ii = 0; ii < nbOfRadius; ++ii )
              {
                if(( v_radius[ii] / h ) < 10.0 ) /// 	ridiculously small radius check
                {
                  v_radius[ii] = 10.0 * h;
                }
                v_estimators[ii] = new Estimator( K, functor );
                v_estimators[ii]->init( h, v_radius[ii] );
              }

              trace.endBlock();

              ///WAIT

              trace.info() << "Curvature computation..." << std::endl;

#ifdef WITH_OPENMP
#pragma omp parallel for schedule(dynamic)
              for( unsigned int ii = 0; ii < pr2size; ++ii )
              {
                /*Estimator intinv( K, functor );
                intinv.init( h, re[ii] );
                Dimension position = 0;*/

                v_curvatures[ii] = v_estimators[ v_radius[ v_registration[ ii ]]]->eval( contour.begin() + ii );
              }
#else
              /*Estimator intinv( K, functor );
              double last_re = -1.0;
              for( unsigned int ii = 0; ii < pr2size; ++ii )
              {
                if( last_re != re[ii] )
                {
                  last_re = re[ii];
                  intinv.init( h, last_re );
                }
                v_curvatures[ii] = intinv.eval( contour.begin() + ii );*/

                v_curvatures[ii] = v_estimators[ v_radius[ v_registration[ ii ]]]->eval( contour.begin() + ii );
              }
#endif

              trace.endBlock();

//              trace.info() << "step C.1" << std::endl;
//              trace.endBlock();

              /// WAIT
              for( Dimension ii = 0; ii < nbOfRadius; ++ii )
              {
                delete v_estimators[ii];
              }

              trace.info() << "Exporting results..." << std::endl;

              for( unsigned int ii = 0; ii < pr2size; ++ii )
              {
                file << v_curvatures[ii] << std::endl;
                file2 << v_estimated_radius[ii] << std::endl;
                file4 << v_radius[ii] << std::endl;
              }

              trace.endBlock();

//              trace.info() << "step D." << std::endl;

//              trace.endBlock();

            }

            double time = c.stopClock();
            file << "# Time: " << time << std::endl;

            file.close();
            file2.close();
            file3.close();
            file4.close();
          }

          /////// local max
          if( optionsII.testsII.at(4) != '0' )
          {
            c.startClock();

            char full_filename[360];
            sprintf( full_filename, "%s%s", filename.c_str(), "_II_curvature_local_max.dat" );
            std::ofstream file( full_filename );

            file << "# h = " << h << std::endl;
            file << "# Integral Invariant curvature estimation" << std::endl;
            file << "# range size = " << pointsRange.size() << std::endl;
            if ( withNoise )
            {
              file << "# noise level (init) = " << noiseLevel/h << std::endl;
              file << "# noise level (current) = " << noiseLevel << std::endl;
            }

            if( optionsII.radius <= 0.0 )
            {
              optionsII.radius = suggestedSizeIntegralInvariant( h, dig->round( optionsII.center ), pointsRange.begin(), pointsRange.end() );
              file << "# Estimated radius: " << optionsII.radius << std::endl;
            }

            double k = optionsII.cste;

            ///////////
            char full_filename2[360];
            sprintf( full_filename2, "%s%s", filename.c_str(), "_II_curvature_local_max_re.dat" );
            std::ofstream file2( full_filename2 );
            ///////////

            if ( withNoise )
            {
              trace.error() << "not yet" << std::endl;
            }
            else
            {
              typedef LightImplicitDigitalSurface< KSpace, Digitizer > LightImplicitDigSurface;
              typedef DigitalSurface< LightImplicitDigSurface > DigSurface;
              typedef PointFunctorFromPointPredicateAndDomain< Digitizer, Domain, unsigned int > MyPointFunctor;
              typedef FunctorOnCells< MyPointFunctor, KSpace > CurvatureIIFct;
              typedef IntegralInvariantMeanCurvatureEstimator< KSpace, CurvatureIIFct > Estimator;
              typedef DepthFirstVisitor< DigSurface > Visitor;
              typedef GraphVisitorRange< Visitor > VisitorRange;
              typedef typename VisitorRange::ConstIterator ConstIterator;

              LightImplicitDigSurface LightImplDigSurf( K, *dig, SAdj, bel );
              DigSurface surf( LightImplDigSurf );

              MyPointFunctor pointFunctor( dig, domain, 1, 0 );
              CurvatureIIFct functor( pointFunctor, K );

              VisitorRange range( new Visitor( surf, *surf.begin() ) );
              ConstIterator ibegin = range.begin();
              ConstIterator iend = range.end();

              std::vector< SCell > contour;

              for( ; ibegin != iend; ++ibegin )
              {
                contour.push_back( *ibegin );
              }

              if( contour.size() != pr2size )
              {
                trace.error() << "ERROR! Not the same border size: " << contour.size() << " " << pr2size << std::endl;
                return 0;
              }

              std::vector< double > resultat( contour.size() );
              std::vector< double > re( contour.size() );

              trace.info() << "step A." << std::endl;

#ifdef WITH_OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
              for( unsigned int ii = 0; ii < pr2size; ++ii )
              {
                double max = v_statMSEL[ pr2size - 1 - ii ].max();
                re[ii] = (k * (max * max)) * h;
              }

              trace.info() << "step B.1" << std::endl;

              ///WAIT

              trace.info() << "step B.2" << std::endl;

#ifdef WITH_OPENMP
#pragma omp parallel for schedule(dynamic)

              for( unsigned int ii = 0; ii < pr2size; ++ii )
              {
                Estimator intinv( K, functor );
                intinv.init( h, re[ii] );
                resultat[ii] = intinv.eval( contour.begin() + ii );
              }
#else
              Estimator intinv( K, functor );
              double last_re = -1.0;
              for( unsigned int ii = 0; ii < pr2size; ++ii )
              {
                if( last_re != re[ii] )
                {
                  last_re = re[ii];
                  intinv.init( h, last_re );
                }
                resultat[ii] = intinv.eval( contour.begin() + ii );
              }
#endif

              trace.info() << "step C.1" << std::endl;

              /// WAIT

              trace.info() << "step C.2" << std::endl;

              for( unsigned int ii = 0; ii < pr2size; ++ii )
              {
                file << resultat[ii] << std::endl;
                file2 << re[ii] << std::endl;
              }

              trace.info() << "step D." << std::endl;

            }

            double time = c.stopClock();
            file << "# Time: " << time << std::endl;

            file.close();
            file2.close();
          }

          /////// local median
          if( optionsII.testsII.at(5) != '0' )
          {
            c.startClock();

            char full_filename[360];
            sprintf( full_filename, "%s%s", filename.c_str(), "_II_curvature_local_median.dat" );
            std::ofstream file( full_filename );

            file << "# h = " << h << std::endl;
            file << "# Integral Invariant curvature estimation" << std::endl;
            file << "# range size = " << pointsRange.size() << std::endl;
            if ( withNoise )
            {
              file << "# noise level (init) = " << noiseLevel/h << std::endl;
              file << "# noise level (current) = " << noiseLevel << std::endl;
            }

            if( optionsII.radius <= 0.0 )
            {
              optionsII.radius = suggestedSizeIntegralInvariant( h, dig->round( optionsII.center ), pointsRange.begin(), pointsRange.end() );
              file << "# Estimated radius: " << optionsII.radius << std::endl;
            }

            double k = optionsII.cste;

            ///////////
            char full_filename2[360];
            sprintf( full_filename2, "%s%s", filename.c_str(), "_II_curvature_local_median_re.dat" );
            std::ofstream file2( full_filename2 );
            ///////////

            if ( withNoise )
            {
              trace.error() << "not yet" << std::endl;
            }
            else
            {
              typedef LightImplicitDigitalSurface< KSpace, Digitizer > LightImplicitDigSurface;
              typedef DigitalSurface< LightImplicitDigSurface > DigSurface;
              typedef PointFunctorFromPointPredicateAndDomain< Digitizer, Domain, unsigned int > MyPointFunctor;
              typedef FunctorOnCells< MyPointFunctor, KSpace > CurvatureIIFct;
              typedef IntegralInvariantMeanCurvatureEstimator< KSpace, CurvatureIIFct > Estimator;
              typedef DepthFirstVisitor< DigSurface > Visitor;
              typedef GraphVisitorRange< Visitor > VisitorRange;
              typedef typename VisitorRange::ConstIterator ConstIterator;

              LightImplicitDigSurface LightImplDigSurf( K, *dig, SAdj, bel );
              DigSurface surf( LightImplDigSurf );

              MyPointFunctor pointFunctor( dig, domain, 1, 0 );
              CurvatureIIFct functor( pointFunctor, K );

              VisitorRange range( new Visitor( surf, *surf.begin() ) );
              ConstIterator ibegin = range.begin();
              ConstIterator iend = range.end();

              std::vector< SCell > contour;

              for( ; ibegin != iend; ++ibegin )
              {
                contour.push_back( *ibegin );
              }

              if( contour.size() != pr2size )
              {
                trace.error() << "ERROR! Not the same border size: " << contour.size() << " " << pr2size << std::endl;
                return 0;
              }

              std::vector< double > resultat( contour.size() );
              std::vector< double > re( contour.size() );

              trace.info() << "step A." << std::endl;

#ifdef WITH_OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
              for( unsigned int ii = 0; ii < pr2size; ++ii )
              {
                double median = v_statMSEL[ pr2size - 1 - ii ].median();
                re[ii] = (k * (median * median)) * h;
              }

              trace.info() << "step B.1" << std::endl;

              ///WAIT

              trace.info() << "step B.2" << std::endl;

#ifdef WITH_OPENMP
#pragma omp parallel for schedule(dynamic)

              for( unsigned int ii = 0; ii < pr2size; ++ii )
              {
                Estimator intinv( K, functor );
                intinv.init( h, re[ii] );
                resultat[ii] = intinv.eval( contour.begin() + ii );
              }
#else
              Estimator intinv( K, functor );
              double last_re = -1.0;
              for( unsigned int ii = 0; ii < pr2size; ++ii )
              {
                if( last_re != re[ii] )
                {
                  last_re = re[ii];
                  intinv.init( h, last_re );
                }
                resultat[ii] = intinv.eval( contour.begin() + ii );
              }
#endif

              trace.info() << "step C.1" << std::endl;

              /// WAIT

              trace.info() << "step C.2" << std::endl;

              for( unsigned int ii = 0; ii < pr2size; ++ii )
              {
                file << resultat[ii] << std::endl;
                file2 << re[ii] << std::endl;
              }

              trace.info() << "step D." << std::endl;

            }

            double time = c.stopClock();
            file << "# Time: " << time << std::endl;

            file.close();
            file2.close();
          }

        }
      }

      //delete noisifiedObject;
      delete dig;
    }
    else
    {
      //delete noisifiedObject;
      delete dig;
      std::cerr << "[computeLocalEstimations]"
                << " error: open digital curve found." << std::endl;
      return false;
    }
  }
  catch ( InputException e )
  {
    std::cerr << "[computeLocalEstimations]"
              << " error in finding a bel." << std::endl;
    return false;
  }
}


///////////////////////////////////////////////////////////////////////////////
namespace po = boost::program_options;

int main( int argc, char** argv )
{
  // parse command line ----------------------------------------------
  po::options_description general_opt("Allowed options are");
  general_opt.add_options()
      ("help,h", "display this message")
      ("list,l",  "List all available shapes")
      ("output,o", po::value<std::string>(), "Output")
      ("shape,s", po::value<std::string>(), "Shape name")
      ("radius,R",  po::value<double>(), "Radius of the shape" )
      //("kernelradius,K",  po::value<double>()->default_value(0.0), "Radius of the convolution kernel (Integral invariants estimators)" )
      ("alpha",  po::value<double>()->default_value(1.0/3.0), "Alpha parameter for Integral Invariant computation" )
      ("axis1,A",  po::value<double>(), "Half big axis of the shape (ellipse)" )
      ("axis2,a",  po::value<double>(), "Half small axis of the shape (ellipse)" )
      ("smallradius,r",  po::value<double>()->default_value(5), "Small radius of the shape" )
      ("varsmallradius,v",  po::value<double>()->default_value(5), "Variable small radius of the shape" )
      ("k,k",  po::value<unsigned int>()->default_value(3), "Number of branches or corners the shape" )
      ("phi",  po::value<double>()->default_value(0.0), "Phase of the shape (in radian)" )
      ("width,w",  po::value<double>()->default_value(10.0), "Width of the shape" )
      ("power,p",   po::value<double>()->default_value(2.0), "Power of the metric (double)" )
      ("center_x,x",   po::value<double>()->default_value(0.0), "x-coordinate of the shape center (double)" )
      ("center_y,y",   po::value<double>()->default_value(0.0), "y-coordinate of the shape center (double)" )
      ("gridstep,g",  po::value<double>()->default_value(1.0), "Grid step for the digitization" )
      ("cste",  po::value<double>()->default_value(0.1), "Constante for II estimator" )
      ("noise,n",  po::value<double>()->default_value(0.0), "Level of noise to perturb the shape" )
      ("properties",  po::value<std::string>()->default_value("11"), "the i-th property is disabled iff there is a 0 at position i" )
      ("estimators,e",  po::value<std::string>()->default_value("10000"), "the i-th estimator is disabled iff there is a 0 at position i" )
      ("testsII,t",  po::value<std::string>()->default_value("111111"), "the i-th test for II estimator is disabled iff there is a 0 at position i" )
      ("lambda,l",  po::value< bool >()->default_value( false ), "Use the shape to get a better approximation of the surface (optional)" );


  bool parseOK=true;
  po::variables_map vm;
  try{
    po::store(po::parse_command_line(argc, argv, general_opt), vm);
  }catch(const std::exception& ex){
    parseOK=false;
    trace.info()<< "Error checking program options: "<< ex.what()<< std::endl;
  }
  po::notify(vm);
  if(!parseOK || vm.count("help")||argc<=1)
  {
    trace.info()<< "Compare local estimators on implicit shapes using DGtal library" <<std::endl
                << "Basic usage: "<<std::endl
                << "\tlocalEstimators --output <output> --shape <shapeName> [required parameters] --estimators <binaryWord> --properties <binaryWord>"<<std::endl
                << std::endl
                << "Below are the different available families of estimators: " << std::endl
                << "\t - True estimators" << std::endl
                << "\t - Maximal DSS based estimators" << std::endl
                << "\t - Maximal DCA based estimators" << std::endl
                << "\t - Binomial convolver based estimators" << std::endl
                << "\t - Integral Invariants based estimators" << std::endl
                << std::endl
                << "The i-th family of estimators is enabled if the i-th character of the binary word is not 0. "
                << "The default binary word is '10000'. This means that the first family of estimators, "
                << "ie. true estimators, is enabled, whereas the next ones are disabled. "
                << std::endl
                << "Below are the different available properties: " << std::endl
                << "\t - Tangeant" << std::endl
                << "\t - Curvature" << std::endl
                << std::endl
                << "Below are the different tests available for Integral Invariant estimators: " << std::endl
                << "\t - Global mean maximal segment based radius" << std::endl
                << "\t - Global max maximal segment based radius" << std::endl
                << "\t - Global median maximal segment based radius" << std::endl
                << "\t - Local mean maximal segment based radius" << std::endl
                << "\t - Local max maximal segment based radius" << std::endl
                << "\t - Local median maximal segment based radius" << std::endl
                << std::endl
                << general_opt << std::endl;
    return 0;
  }

  //List creation
  createList();

  if (vm.count("list"))
  {
    displayList();
    return 0;
  }

  //Parse options
  if (!(vm.count("shape"))) missingParam("--shape");
  if (!(vm.count("output"))) missingParam("--output");

  std::string shapeName = vm["shape"].as<std::string>();
  std::string filename = vm["output"].as<std::string>();

  int nb = 5; //number of available methods
  std::string options = vm["estimators"].as< std::string >();
  if (options.size() < nb)
  {
    trace.error() << " At least " << nb
                  << " characters are required "
                  << " with option --estimators.";
    trace.info() << std::endl;
    exit(1);
  }

  nb = 2; //number of available properties
  std::string properties = vm["properties"].as<std::string>();
  if (properties.size() < nb)
  {
    trace.error() << " At least " << nb
                  << " characters are required "
                  << " with option --properties.";
    trace.info() << std::endl;
    exit(1);
  }


  //We check that the shape is known
  unsigned int id = checkAndReturnIndex(shapeName);

  // standard types
  typedef Z2i::Space Space;
  typedef Space::Point Point;
  typedef Space::RealPoint RealPoint;

  RealPoint center( vm["center_x"].as<double>(),
      vm["center_y"].as<double>() );
  double h = vm["gridstep"].as<double>();

  struct OptionsIntegralInvariant< RealPoint > optII;
  //optII.radius = vm["kernelradius"].as<double>();
  optII.alpha = vm["alpha"].as<double>();
  optII.lambda_optimized = vm["lambda"].as< bool >();
  optII.center = center;
  optII.cste = vm["cste"].as<double>();

  nb = 6; //number of available tests for II
  optII.testsII = vm["testsII"].as<std::string>();
  if (optII.testsII.size() < nb)
  {
    trace.error() << " At least " << nb
                  << " characters are required "
                  << " with option --testsII.";
    trace.info() << std::endl;
    exit(1);
  }

  double noiseLevel = vm["noise"].as<double>();

  if (id ==0)
  {
    if (!(vm.count("radius"))) missingParam("--radius");
    //if (!(vm.count("kernelradius"))) missingParam("--kernelradius");
    double radius = vm["radius"].as<double>();

    Ball2D<Space> * ball = new Ball2D<Space>( center, radius);
    computeLocalEstimations<Space>( filename, ball, h, optII, options, properties, noiseLevel );
    delete ball;
  }
  else if (id ==1)
  {
    if (!(vm.count("width"))) missingParam("--width");
    double width = vm["width"].as<double>();

    ImplicitHyperCube<Space> object(Z2i::Point(0,0), width/2);
    trace.error()<< "Not available.";
    trace.info()<<std::endl;
  }
  else if (id ==2)
  {
    if (!(vm.count("power"))) missingParam("--power");
    if (!(vm.count("radius"))) missingParam("--radius");
    double radius = vm["radius"].as<double>();
    double power = vm["power"].as<double>();

    ImplicitRoundedHyperCube<Space> ball( Z2i::Point(0,0), radius, power );
    trace.error()<< "Not available.";
    trace.info()<<std::endl;
  }
  else if (id ==3)
  {
    if (!(vm.count("varsmallradius"))) missingParam("--varsmallradius");
    if (!(vm.count("radius"))) missingParam("--radius");
    if (!(vm.count("k"))) missingParam("--k");
    if (!(vm.count("phi"))) missingParam("--phi");
    //if (!(vm.count("kernelradius"))) missingParam("--kernelradius");
    double radius = vm["radius"].as<double>();
    double varsmallradius = vm["varsmallradius"].as<double>();
    unsigned int k = vm["k"].as<unsigned int>();
    double phi = vm["phi"].as<double>();

    Flower2D<Space> * flower = new Flower2D<Space>( center, radius, varsmallradius, k, phi );
    computeLocalEstimations<Space>( filename, flower, h, optII, options, properties, noiseLevel );
    delete flower;
  }
  else if (id ==4)
  {
    if (!(vm.count("radius"))) missingParam("--radius");
    if (!(vm.count("k"))) missingParam("--k");
    if (!(vm.count("phi"))) missingParam("--phi");
    //if (!(vm.count("kernelradius"))) missingParam("--kernelradius");
    double radius = vm["radius"].as<double>();
    unsigned int k = vm["k"].as<unsigned int>();
    double phi = vm["phi"].as<double>();

    NGon2D<Space> * object = new NGon2D<Space>( center, radius, k, phi );
    computeLocalEstimations<Space>( filename, object, h, optII, options, properties, noiseLevel );
    delete object;
  }
  else if (id ==5)
  {
    if (!(vm.count("varsmallradius"))) missingParam("--varsmallradius");
    if (!(vm.count("radius"))) missingParam("--radius");
    if (!(vm.count("k"))) missingParam("--k");
    if (!(vm.count("phi"))) missingParam("--phi");
    //if (!(vm.count("kernelradius"))) missingParam("--kernelradius");
    double radius = vm["radius"].as<double>();
    double varsmallradius = vm["varsmallradius"].as<double>();
    unsigned int k = vm["k"].as<unsigned int>();
    double phi = vm["phi"].as<double>();

    AccFlower2D<Space> * accflower = new AccFlower2D<Space>( center, radius, varsmallradius, k, phi );
    computeLocalEstimations<Space>( filename, accflower, h, optII, options, properties, noiseLevel );
    delete accflower;
  }
  else if (id ==6)
  {
    if (!(vm.count("axis1"))) missingParam("--axis1");
    if (!(vm.count("axis2"))) missingParam("--axis2");
    if (!(vm.count("phi"))) missingParam("--phi");
    //if (!(vm.count("kernelradius"))) missingParam("--kernelradius");
    double a1 = vm["axis1"].as<double>();
    double a2 = vm["axis2"].as<double>();
    double phi = vm["phi"].as<double>();

    Ellipse2D<Space> * ellipse = new Ellipse2D<Space>( center, a1, a2, phi );
    computeLocalEstimations<Space>( filename, ellipse, h, optII, options, properties, noiseLevel );
    delete ellipse;
  }
}
