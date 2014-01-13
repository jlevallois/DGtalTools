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

// Integral Invariant includes
#include "DGtal/geometry/surfaces/FunctorOnCells.h"
#include "DGtal/geometry/surfaces/estimation/IntegralInvariantMeanCurvatureEstimator.h"
#include "DGtal/geometry/surfaces/estimation/IntegralInvariantGaussianCurvatureEstimator.h"

// Drawing
#include "DGtal/io/viewers/Viewer3D.h"
#include "DGtal/io/colormaps/GradientColorMap.h"
#include "DGtal/topology/DigitalSurface2DSlice.h"
#include <QtGui/QApplication>

using namespace std;
using namespace DGtal;

const Color  AXIS_COLOR_RED( 200, 20, 20, 255 );
const Color  AXIS_COLOR_GREEN( 20, 200, 20, 255 );
const Color  AXIS_COLOR_BLUE( 20, 20, 200, 255 );
const double AXIS_LINESIZE = 0.05;


///////////////////////////////////////////////////////////////////////////////

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

int main( int argc, char** argv )
{
  std::vector< double > v_test;
  v_test.push_back(1.5);
  v_test.push_back(2.5);
  v_test.push_back(3.5);
  v_test.push_back(4.5);
  v_test.push_back(5.5);

  typedef SCellTo2DSCell< std::vector< double > > Embedder;
  Embedder TEST(v_test);
  typedef Embedder::ConstIterator ConstIterator;
  typedef Embedder::Iterator Iterator;
  ConstIterator itb = TEST.begin();
  ConstIterator ite = TEST.end();

  std::cout << "itb=" << *v_test.begin() << std::endl;
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

  return 0;

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
  typedef KSpace::Surfel Surfel;
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
  Surfel bel = Surfaces< Z3i::KSpace >::findABel( K, predicate, 100000 );
  MyLightImplicitDigitalSurface LightImplDigSurf( K, predicate, SAdj, bel );
  MyDigitalSurface digSurf( LightImplDigSurf );

  typedef DepthFirstVisitor<MyDigitalSurface> Visitor;
  typedef GraphVisitorRange< Visitor > VisitorRange;
  typedef VisitorRange::ConstIterator SurfelConstIterator;
  VisitorRange range( new Visitor( digSurf, *digSurf.begin() ) );
  SurfelConstIterator abegin = range.begin();
  SurfelConstIterator aend = range.end();

  typedef ImageToConstantFunctor< Image, ImagePredicate > MyPointFunctor;
  MyPointFunctor pointFunctor( image, predicate, 1 );

  // Integral Invariant stuff

  typedef FunctorOnCells< MyPointFunctor, Z3i::KSpace > MyCellFunctor;
  MyCellFunctor functor ( pointFunctor, K ); // Creation of a functor on Cells, returning true if the cell is inside the shape

  QApplication application( argc, argv );
  typedef Viewer3D<Z3i::Space, Z3i::KSpace> Viewer;
  Viewer viewer( K );
  viewer.show();
  //    viewer << SetMode3D(image.domain().className(), "BoundingBox") << image.domain();

  VisitorRange range2( new Visitor( digSurf, *digSurf.begin() ) );
  SurfelConstIterator abegin2 = range2.begin();

  trace.beginBlock("curvature computation");
  {
    /*typedef double Quantity;
    std::vector< Quantity > results;
    back_insert_iterator< std::vector< Quantity > > resultsIterator( results ); // output iterator for results of Integral Invariante curvature computation

    typedef IntegralInvariantMeanCurvatureEstimator< Z3i::KSpace, MyCellFunctor > MyIIMeanEstimator;

    MyIIMeanEstimator estimator ( K, functor );
    estimator.init( h, re_convolution_kernel ); // Initialisation for a given Euclidean radius of the convolution kernel
    estimator.eval ( abegin, aend, resultsIterator ); // Computation
*/
    // Drawing results
    /*Quantity min = numeric_limits < Quantity >::max();
    Quantity max = numeric_limits < Quantity >::min();
    for ( unsigned int i = 0; i < results.size(); ++i )
    {
      if ( results[ i ] < min )
      {
        min = results[ i ];
      }
      else if ( results[ i ] > max )
      {
        max = results[ i ];
      }
    }

    typedef GradientColorMap< Quantity > Gradient;
    Gradient cmap_grad( min, max );
    cmap_grad.addColor( Color( 50, 50, 255 ) );
    cmap_grad.addColor( Color( 255, 0, 0 ) );
    cmap_grad.addColor( Color( 255, 255, 10 ) );*/

    typedef typename MyLightImplicitDigitalSurface::Tracker Tracker;
    typedef DigitalSurface2DSlice< Tracker > MySlice;

    typename Tracker::Surfel oneSCell;
    Dimension ii = 0;
    while( abegin2 != aend )//for ( unsigned int i = 0; i < results.size(); ++i )
    {
      if( ii == 15 )
      {
        oneSCell = *abegin2;
      }
      //viewer << CustomColors3D( Color::White, Color::White )
      //       << *abegin2;
      ++abegin2;
      ++ii;
    }

//    Surfel pointel = K.sPointel( K.sKCoords( oneSCell ));

    Tracker* ptrTracker = new Tracker( LightImplDigSurf, oneSCell ); // some pointer on a tracker.
    MySlice slicex( ptrTracker, 0 ); // slice containing x-axis
    MySlice slicey( ptrTracker, 1 ); // slice containing y-axis
    MySlice slicez( ptrTracker, 2 ); // slice containing z-axis
    for( typename MySlice::ConstIterator itb = slicex.begin(); itb != slicex.end(); ++itb )
    {
      viewer << CustomColors3D( Color::Green, Color::Green )
             << *itb;
    }
    for( typename MySlice::ConstIterator itb = slicey.begin(); itb != slicey.end(); ++itb )
    {
      viewer << CustomColors3D( Color::Red, Color::Red )
             << *itb;
    }
    for( typename MySlice::ConstIterator itb = slicez.begin(); itb != slicez.end(); ++itb )
    {
      viewer << CustomColors3D( Color::Blue, Color::Blue )
             << *itb;
    }

/*    viewer << CustomColors3D( Color::Red, Color::Red )
         << oneSCell*/;

  }

  viewer << Viewer3D<>::updateDisplay;
  return application.exec();
}

///////////////////////////////////////////////////////////////////////////////
