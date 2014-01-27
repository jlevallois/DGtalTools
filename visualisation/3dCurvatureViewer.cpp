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
#include "DGtal/io/readers/RawReader.h"
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
#include "DGtal/io/boards/Board3D.h"
#include "DGtal/io/viewers/Viewer3D.h"
#include "DGtal/io/colormaps/GradientColorMap.h"
#include <QtGui/QApplication>

using namespace std;
using namespace DGtal;

const Color  AXIS_COLOR_RED( 200, 20, 20, 255 );
const Color  AXIS_COLOR_GREEN( 20, 200, 20, 255 );
const Color  AXIS_COLOR_BLUE( 20, 20, 200, 255 );
const Color  AXIS_COLOR_BLACK( 200, 200, 200, 255 );
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

enum FileExtension
{
  VOL,
  RAW
};

template< typename Image >
Image FileReader( std::string filename, FileExtension extension )
{
  switch( extension )
  {
  case FileExtension::VOL:
    return VolReader<Image>::importVol( filename );
    break;

  case FileExtension::RAW:
    return RawReader< Image >::importRaw8( filename, Z3i::Domain::Vector( 512, 512, 512) );
    break;
  }

}

template< typename KSpace, typename Domain, typename Surfel >//= KSpace::Surfel >
bool CheckIfIsInBorder( const KSpace & k, const Domain & d, const Surfel & s, const unsigned int radius = 1 )
{
  typedef typename KSpace::Space::RealPoint MidPoint;
  typedef typename Domain::Point Point;
  typedef CanonicSCellEmbedder< KSpace > Embedder;

  Embedder embedder( k );

  MidPoint mp = embedder.embed( s );
  Point min = d.lowerBound() + Point( radius, radius, radius );
  Point max = d.upperBound() - Point( radius, radius, radius );
  const Dimension dimension = 3;
  for( Dimension ii = 0; ii < dimension; ++ii )
  {
    if( mp[ ii ] < min[ ii ] ||  mp[ ii ] > max [ ii ])
    {
      return true;
    }
  }
  return false;
}

int main( int argc, char** argv )
{
  // parse command line ----------------------------------------------
  po::options_description general_opt("Allowed options are");
  general_opt.add_options()
    ("help,h", "display this message")
    ("input-file,i", po::value< std::string >(), ".vol file")
    ("radius,r",  po::value< double >(), "Kernel radius for IntegralInvariant" )
    ("extension,e", po::value< std::string >()->default_value("vol"), "feiofjreio" )
    ("properties,p", po::value< std::string >()->default_value("mean"), "type of output : mean, gaussian, prindir1 or prindir2 (default mean)");

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
  if (!(vm.count("radius"))){ 
    missingParam("--radius");
    neededArgsGiven=false;
  }  
  double h = 1.0;

  FileExtension extension;
  std::string ext = vm["extension"].as< std::string >();
  if( ext == "raw" )
  {
    extension = FileExtension::RAW;
  }
  else if(ext == "vol" )
  {
    extension = FileExtension::VOL;
  }

  bool wrongMode = false;
  std::string mode = vm["properties"].as< std::string >();
  if (( mode.compare("gaussian") != 0 ) && ( mode.compare("mean") != 0 ) && ( mode.compare("prindir1") != 0 ) && ( mode.compare("prindir2") != 0 ))
  {
    wrongMode = true;
  }

  if(!neededArgsGiven ||  wrongMode || !parseOK || vm.count("help") || argc <= 1 )
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
  double re_convolution_kernel = vm["radius"].as< double >();
  
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
  Image image = FileReader< Image >( filename, extension );
  {
    Z3i::DigitalSet set_temp( image.domain() );
    SetFromImage< Z3i::DigitalSet >::append< Image >( set_temp, image, 0, 255 );

    image = ImageFromSet< Image >::create< Z3i::DigitalSet >( set_temp, 255, true, set_temp.begin(), set_temp.end() );
  }
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
  /*typedef Viewer3D<Z3i::Space, Z3i::KSpace> Viewer;
  Viewer viewer( K );
  viewer.show();*/
  typedef Board3D< Z3i::Space, Z3i::KSpace> Board3D;
  Board3D viewer( K );

  //    viewer << SetMode3D(image.domain().className(), "BoundingBox") << image.domain();

  VisitorRange range2( new Visitor( digSurf, *digSurf.begin() ) );
  SurfelConstIterator abegin2 = range2.begin();

  trace.beginBlock("curvature computation");
  if( ( mode.compare("gaussian") == 0 ) || ( mode.compare("mean") == 0 ) )
  {
    typedef double Quantity;
    std::vector< Quantity > results;
    back_insert_iterator< std::vector< Quantity > > resultsIterator( results ); // output iterator for results of Integral Invariante curvature computation

    if ( ( mode.compare("mean") == 0 ) )
    {
      typedef IntegralInvariantMeanCurvatureEstimator< Z3i::KSpace, MyCellFunctor > MyIIMeanEstimator;

      MyIIMeanEstimator estimator ( K, functor );
      estimator.init( h, re_convolution_kernel ); // Initialisation for a given Euclidean radius of the convolution kernel
      estimator.eval ( abegin, aend, resultsIterator ); // Computation
    }
    else if ( ( mode.compare("gaussian") == 0 ) )
    {
      typedef IntegralInvariantGaussianCurvatureEstimator< Z3i::KSpace, MyCellFunctor > MyIIGaussianEstimator;
      typedef CurvatureInformations CurvInformation;
      std::vector< CurvInformation > results2;
      back_insert_iterator< std::vector< CurvInformation > > resultsIterator2( results2 );


      MyIIGaussianEstimator estimator ( K, functor );
      estimator.init( h, re_convolution_kernel ); // Initialisation for a given Euclidean radius of the convolution kernel
      estimator.evalPrincipalCurvatures ( abegin, aend, resultsIterator2 ); // Computation

      for(int i = 0; i < results2.size(); ++i)
      {
        results.push_back(results2[i].k1 * results2[i].k2);
      }
    }

    // Drawing results
    VisitorRange range3( new Visitor( digSurf, *digSurf.begin() ) );
    SurfelConstIterator abegin3 = range3.begin();
    Quantity min = numeric_limits < Quantity >::max();
    Quantity max = numeric_limits < Quantity >::min();
    for ( unsigned int ii = 0; ii < results.size(); ++ii )
    {
      Surfel current = *abegin3;
      if( CheckIfIsInBorder< KSpace, Z3i::Domain, Surfel >( K, domain, current, 1 ) )
      {
        results[ ii ] = numeric_limits < Quantity >::min();
      }
      else
      {
        if ( results[ ii ] < min )
        {
          min = results[ ii ];
        }
        else if ( results[ ii ] > max )
        {
          max = results[ ii ];
        }
      }
      ++abegin3;
    }

    typedef GradientColorMap< Quantity > Gradient;
    Gradient cmap_grad( min, max );
    cmap_grad.addColor( Color( 50, 50, 255 ) );
    cmap_grad.addColor( Color( 255, 0, 0 ) );
    cmap_grad.addColor( Color( 255, 255, 10 ) );

    Cell currentCell;

    for ( unsigned int ii = 0; ii < results.size(); ++ii )
    {
      currentCell = Cell( K.sKCoords(*abegin2) );
      if( results[ ii ] != numeric_limits < Quantity >::min() )
      {
        viewer << CustomColors3D( Color::Black, cmap_grad( results[ ii ] ))
               << currentCell;
      }
      else
      {
        viewer << CustomColors3D( Color( 255, 255, 255, 120 ), Color( 255, 255, 255, 120 ) )
               << currentCell;
      }
      ++abegin2;
    }
  }
  else
  {
    typedef double Quantity;
    typedef EigenValues3D< Quantity >::Matrix33 Matrix3x3;
    typedef EigenValues3D< Quantity >::Vector3 Vector3;
    typedef CurvatureInformations CurvInformation;

    std::vector< CurvInformation > results;
    back_insert_iterator< std::vector< CurvInformation > > resultsIterator( results ); // output iterator for results of Integral Invariante curvature computation

    typedef IntegralInvariantGaussianCurvatureEstimator< Z3i::KSpace, MyCellFunctor > MyIIGaussianEstimator;

    MyIIGaussianEstimator estimator ( K, functor );
    estimator.init ( h, re_convolution_kernel ); // Initialisation for a given Euclidean radius of the convolution kernel
    estimator.evalPrincipalCurvatures ( abegin, aend, resultsIterator ); // Computation

    trace.endBlock();

    trace.beginBlock("viewer");

    // Drawing results
    typedef  Matrix3x3::RowVector RowVector;
    typedef  Matrix3x3::ColumnVector ColumnVector;

    for ( unsigned int i = 0; i < results.size(); ++i )
    {
      CurvInformation current = results[ i ];
      DGtal::Dimension kDim = K.sOrthDir( *abegin2 );
      SCell outer = K.sIndirectIncident( *abegin2, kDim);
      if ( predicate(embedder(outer)) )
      {
        outer = K.sDirectIncident( *abegin2, kDim);
      }

      Cell unsignedSurfel = K.uCell( K.sKCoords(*abegin2) );
      viewer << CustomColors3D( DGtal::Color(255,255,255,255),
                                DGtal::Color(255,255,255,255))
             << unsignedSurfel;


      //ColumnVector normal = current.vectors.column(0).getNormalized(); // don't show the normal
      ColumnVector curv1 = current.vectors.column(1).getNormalized();
      ColumnVector curv2 = current.vectors.column(2).getNormalized();

      double eps = 0.01;
      RealPoint center = embedder( outer );// + eps*embedder( *abegin2 );

      //            viewer.addLine ( center[0] - 0.5 * normal[ 0],
      //                             center[1] - 0.5 * normal[1],
      //                             center[2] - 0.5* normal[2],
      //                             center[0] +  0.5 * normal[0],
      //                             center[1] +  0.5 * normal[1],
      //                             center[2] +  0.5 * normal[2],
      //                             DGtal::Color ( 0,0,0 ), 5.0 ); // don't show the normal


      if( ( mode.compare("prindir1") == 0 ) )
      {
        viewer.setLineColor(AXIS_COLOR_BLUE);
        viewer.addLine (
              RealPoint(
                center[0] -  0.5 * curv1[0],
            center[1] -  0.5 * curv1[1],
            center[2] -  0.5 * curv1[2]
            ),
            RealPoint(
              center[0] +  0.5 * curv1[0],
            center[1] +  0.5 * curv1[1],
            center[2] +  0.5 * curv1[2]
            ),
            AXIS_LINESIZE );
      }
      else
      {
        viewer.setLineColor(AXIS_COLOR_RED);
        viewer.addLine (
              RealPoint(
                center[0] -  0.5 * curv2[0],
            center[1] -  0.5 * curv2[1],
            center[2] -  0.5 * curv2[2]
            ),
            RealPoint(
              center[0] +  0.5 * curv2[0],
            center[1] +  0.5 * curv2[1],
            center[2] +  0.5 * curv2[2]
            ),
            AXIS_LINESIZE );
      }

      ++abegin2;
    }
    trace.endBlock();
  }

  /*viewer << Viewer3D<>::updateDisplay;*/
  viewer.saveOBJ("/media/kha/7aa7fb31-a0b0-492d-8070-a05b66c1771a/snow.obj");
  return application.exec();
}

///////////////////////////////////////////////////////////////////////////////
