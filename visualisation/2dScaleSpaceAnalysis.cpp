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
 * @file 2dScaleSpaceAnalysis.cpp
 * @ingroup Tools
 * @author Jeremy Levallois (\c jeremy.levallois@liris.cnrs.fr )
 * Laboratoire d'InfoRmatique en Image et Systèmes d'information - LIRIS (CNRS, UMR 5205), INSA-Lyon, France
 * LAboratoire de MAthématiques - LAMA (CNRS, UMR 5127), Université de Savoie, France
 *
 * @date 2012/10/18
 *
 * IntegralInvariant curvature estimator for a given kernel radius.
 *
 * This file is part of the DGtalTools.
 */

///////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <sstream>
#include <limits>

#include "DGtal/base/Common.h"
#include "DGtal/base/Clock.h"

//shapes
#include "DGtal/shapes/ShapeFactory.h"
#include "DGtal/shapes/Shapes.h"
#include "DGtal/kernel/BasicPointFunctors.h"

//Digitizer
#include "DGtal/topology/LightImplicitDigitalSurface.h"
#include "DGtal/topology/DigitalSurface.h"
#include "DGtal/graph/DepthFirstVisitor.h"
#include "DGtal/graph/GraphVisitorRange.h"


//Estimator
#include "DGtal/geometry/surfaces/estimation/IntegralInvariantMeanCurvatureEstimator.h"
#include "DGtal/geometry/surfaces/FunctorOnCells.h"


//Image & Viewer3D
#include "DGtal/images/ImageSelector.h"
#include "DGtal/io/colormaps/GradientColorMap.h"
#include "DGtal/io/writers/PPMWriter.h"

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/parsers.hpp>

#define uint unsigned int

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



template< typename MyShape >
int Compute( const MyShape & shape,
             const std::vector< double > & re_array,
             const double h,
             const std::string & options,
             const std::string name )
{
  typedef ImageSelector< Z2i::Domain, float >::Type Image;

  Image * image;

  typedef Z2i::KSpace::Surfel Surfel;
  typedef GaussDigitizer< Z2i::Space, MyShape > Digitizer;
  typedef typename Digitizer::RealPoint RealPoint;

  Digitizer dig;
  dig.attach( shape );
  dig.init( shape.getLowerBound(), shape.getUpperBound(), h );

  Z2i::Domain domain = dig.getDomain();
  Z2i::KSpace K;
  bool space_ok = K.init( domain.lowerBound(), domain.upperBound(), true );
  if ( !space_ok )
  {
    trace.error() << "Error in the Khalimsky space construction."<<std::endl;
    return 0;
  }

  typedef LightImplicitDigitalSurface< Z2i::KSpace, Digitizer > LightDigitalSurface;
  typedef DigitalSurface< LightDigitalSurface > MyDigitalSurface;

  SurfelAdjacency< Z2i::KSpace::dimension > SAdj( true );
  Surfel bel = Surfaces< Z2i::KSpace >::findABel( K, dig, 100000 );
  LightDigitalSurface ldigsurf( K, dig, SAdj, bel );
  MyDigitalSurface surf( ldigsurf );

  typedef PointFunctorFromPointPredicateAndDomain< Digitizer, Z2i::Domain, unsigned int > MyPointFunctor;
  typedef FunctorOnCells< MyPointFunctor, Z2i::KSpace > MyCellFunctor;
  typedef IntegralInvariantMeanCurvatureEstimator< Z2i::KSpace, MyCellFunctor > Estimator;
  typedef typename Estimator::Quantity Quantity;
  typedef DepthFirstVisitor< MyDigitalSurface > Visitor;
  typedef GraphVisitorRange< Visitor > VisitorRange;
  typedef typename VisitorRange::ConstIterator SurfelConstIterator;
  typedef std::vector< Quantity > vQuantity;
  typedef typename vQuantity::const_iterator const_interatorQuantity;

  std::vector< Quantity > data;
  uint data_size = 0;
  uint surf_size = 0;

  const uint c_re_size = re_array.size ();

  Quantity min = std::numeric_limits< Quantity >::max();
  Quantity max = std::numeric_limits< Quantity >::min();
  MyPointFunctor pointFunctor( dig, domain, 1, 0 );
  MyCellFunctor cfunctor( pointFunctor, K );

  /// count surface
  {
    VisitorRange range( new Visitor( surf, *surf.begin() ) );
    SurfelConstIterator abegin = range.begin();
    SurfelConstIterator aend = range.end();

    while( abegin != aend )
    {
      ++surf_size;
      ++abegin;
    }

    Z2i::Point aa( 0, 0 );
    Z2i::Point bb( surf_size - 1, c_re_size - 1 );
    Z2i::Domain domain_image( aa, bb );
    image = new Image( domain_image );

    data_size = image->size();// surf_size * re_size;
    data.resize( data_size );
  }

  const uint c_data_size = data_size;
  const uint c_surf_size = surf_size;

  uint step_re;
  Visitor * vis;

  #pragma omp parallel shared(image,data,min,max,K,cfunctor,surf)
  {
    #pragma omp for schedule(dynamic,1) private(step_re,vis) nowait
    for ( step_re = 0; step_re < c_re_size; ++step_re )
    {
      Estimator estimator( K, cfunctor );

      Clock c;
      c.startClock();

//      Z2i::SCell center( Z2i::Point( 0, 0 ), true);
//      Quantity tcenter = cfunctor( center );
//      trace.warning() << tcenter << std::endl;

      estimator.init( h, re_array[ step_re ] );

//      results = vQuantity(surf_size, step_re); /// http://media.giphy.com/media/H8zg5nlGWWQ7K/giphy.gif
      vQuantity results;
      std::back_insert_iterator< vQuantity > resultIterator( results );
      vis = ::new Visitor( surf, *surf.begin() );
//      trace.info() << vis << std::endl;
      VisitorRange range( vis );
//      trace.info() << *range.begin() << std::endl;
//      trace.info() << &cfunctor << std::endl;
      SurfelConstIterator abegin = range.begin();
      SurfelConstIterator aend = range.end();
//      trace.info() << &aend << std::endl;
//      estimator.eval( abegin, aend, resultIterator );

      uint cc = 0;
      if( step_re % 2 == 0 )
      {
      while ( abegin != aend )
      {
        ++cc;
        trace.info() << step_re << " " << cc << std::endl;
        ++abegin;
      }
      }

      Quantity current_result = Quantity(0);

//      if( options.at( 0 ) != '0' ) // export as image
//      {
//        if ( step_re == 0 ) // compute image domain
//        {
//          Z2i::Point aa( 0, 0 );
//          Z2i::Point bb( rsize - 1, re_size - 1 );
//          Z2i::Domain domain_image( aa, bb );
//          image = new Image( domain_image );
//        }

        for ( uint ii = 0; ii < c_surf_size; ++ii )
        {
          current_result = results[ ii ];
          if ( current_result < min )
          {
            min = current_result;
          }
          else if ( current_result > max )
          {
            max = current_result;
          }

          uint current_pos = step_re * c_surf_size + ii;
          if( current_pos < c_data_size )
          {
            (*image)[ step_re * c_surf_size + ii ] = current_result;
          }
          else
          {
            (*image)[current_pos] = 0;
            trace.error() << "current_pos: " << current_pos << " c_data_size:" << c_data_size << std::endl;
          }
        }

//        typedef GradientColorMap< Quantity > Gradient;
//        Gradient cmap_grad( min, max );
//        cmap_grad.addColor( Color( 50, 50, 255 ) );
//        cmap_grad.addColor( Color( 255, 0, 0 ) );
//        cmap_grad.addColor( Color( 255, 255, 10 ) );

//        std::stringstream sstm;
//        sstm << name << ".ppm";
//        std::string exportname = sstm.str();
//        PPMWriter<Image, Gradient >::exportPPM( exportname, *image, cmap_grad ); // Erase previous image with old + new iteration values (the computation on full iteration take long time, so we can see result before ending)
//      }

      /*
      if( options.at( 1 ) != '0' ) // export as data file
      {
//        if ( step_re == 0 )
//        {
//          data_size = rsize * re_size;
//          data.resize( data_size );
//        }

        for ( ii = 0; ii < rsize; ++ii )
        {
          current_result = results[ ii ];
          data[ ii * re_size + step_re ] = current_result;
        }

        std::stringstream sstm;
        sstm << name << ".dat";
        std::string exportname = sstm.str();
        std::ofstream outf( exportname, std::ios::trunc );
        outf.flush();
        if( !outf )
        {
          trace.error() << "IO error with file " << exportname << std::endl;
          //return 0;
        }

        uint id = 0;
        for ( ii = 0; ii < data_size; ++ii )
        {
          if ( step_re == 0 )
          {
            outf << id;
            ++id;
          }

          outf << " ";

          if ( ii % re_size > step_re )
          {
            outf << "NA";
          }
          else
          {
            outf << data[ ii ];
          }
          if ( ii % re_size == re_size - 1 )
          {
            outf << '\n';
          }
        }
        outf.close();
      }
      */

//      if(((( step_re + 1 ) / re_size ) * 100 ) % 5 )
//      {
//        trace.progressBar( step_re + 1, re_size );
//      }
      double time = c.stopClock();
      trace.info() << "-- Done in " << time << " msec. --" << std::endl;
    }
  }

  typedef GradientColorMap< Quantity > Gradient;
  Gradient cmap_grad( min, max );
  trace.info() << "Min: " << min << " | Max: " << max << std::endl;
  cmap_grad.addColor( Color( 50, 50, 255 ) );
  cmap_grad.addColor( Color( 255, 0, 0 ) );
  cmap_grad.addColor( Color( 255, 255, 10 ) );

  std::stringstream sstm;
  sstm << name << ".ppm";
  std::string exportname = sstm.str();
  PPMWriter<Image, Gradient >::exportPPM( exportname, *image, cmap_grad ); // Erase previous image with old + new iteration values (the computation on full iteration take long time, so we can see result before ending)


  delete image;

  return 1;//image;
}

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

  shapes2D.push_back("cublipse");
  shapesDesc.push_back("Cublipse.");
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
 * Missing parameter error message.
 *
 * @param param
 */
void missingParam(std::string param)
{
  trace.error() << " Parameter: " << param << " is required.";
  trace.info() << std::endl;
  exit(1);
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


namespace po = boost::program_options;

int main ( int argc, char** argv )
{
  po::options_description general_opt("Allowed options are");
  general_opt.add_options()
      ("help,h", "display this message")
      ("list,l",  "List all available shapes")
      ("shape,s", po::value<std::string>(), "Shape name")
      ("radius,R",  po::value<double>(), "Radius of the shape" )
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

      ("re_min,re_min",  po::value<double>()->default_value(1.0), "min Euclidean radius of the convolution kernel" )
      ("re_max,re_max",  po::value<double>()->default_value(7.0), "max Euclidean radius of the convolution kernel" )
      ("re_step,re_step",  po::value<double>()->default_value(0.01), "step of the increase of the Euclidean radius of the convolution kernel" )
      ("gridstep_min,g_min",  po::value<double>()->default_value(0.05), "min Grid step for the digitization" )
      ("gridstep_max,g_max",  po::value<double>()->default_value(0.05), "max Grid step for the digitization" )
      ("gridstep_step,g_step",  po::value<double>()->default_value(0.01), "step of the increase of the Grid step for the digitization" )
      ("export,e",  po::value<std::string>()->default_value("11"), "the i-th estimator is disabled iff there is a 0 at position i" )
      ("filename,f",  po::value<std::string>(), "filename for the export" );

  bool parseOK = true;
  po::variables_map vm;
  try
  {
    po::store(po::parse_command_line(argc, argv, general_opt), vm);
  }
  catch( const std::exception& ex )
  {
    parseOK = false;
    trace.info() << "Error checking program options: " << ex.what() << std::endl;
  }

  po::notify(vm);
  if( !parseOK || vm.count("help") || argc <= 1)
  {
    trace.info()<< "Compare local estimators on implicit shapes using DGtal library" <<std::endl
                << "Basic usage: "<<std::endl
                << "\t" << argv[ 0 ] << " --shape <shapeName> [required parameters] --re <min_re> <max_re> <step_re> --gridstep <min_h> <max_h> <step_h> --export <binaryWord> <filename>"<<std::endl
                << std::endl
                << "Below are the different available exports: " << std::endl
                << "\t - As an image" << std::endl
                << "\t - As data" << std::endl
                << std::endl
                << "The i-th export type is used if the i-th character of the binary word is not 0. "
                << "The default binary word is '11'. This means that it export as an image and data file. "
                << std::endl
                << general_opt << std::endl;
    return 0;
  }

  //List creation

  createList();

  if( vm.count( "list" ))
  {
    displayList();
    return 0;
  }

  //Parse options
  if( !( vm.count( "shape" )))
  {
    missingParam( "--shape" );
  }
  std::string shapeName = vm[ "shape" ].as< std::string >();

  double re_init = vm[ "re_min" ].as< double >();
  double re_end  = vm[ "re_max" ].as< double >();
  double re_step = vm[ "re_step" ].as< double >();
  if ( re_init <= 0.0 || re_end <= 0.0 || re_step <= 0.0 )
  {
    trace.error() << " You need to specify positive values with option --re.";
    trace.info() << std::endl;
    exit(1);
  }
  std::vector< double > re_array;
  for ( ; re_init <= re_end; re_init += re_step )
  {
    re_array.push_back ( re_init );
  }

  double rh_init = vm[ "gridstep_min" ].as< double >();
  double rh_end  = vm[ "gridstep_max" ].as< double >();
  double rh_step = vm[ "gridstep_step" ].as< double >();
  if ( rh_init <= 0.0 || rh_end <= 0.0 || rh_step <= 0.0 )
  {
    trace.error() << " You need to specify positive values with option --gridstep.";
    trace.info() << std::endl;
    exit(1);
  }
  std::vector< double > rh_array;
  for ( ; rh_init >= rh_end; rh_init -= rh_step )
  {
    rh_array.push_back ( rh_init );
  }

  std::string options = vm[ "export" ].as< std::string >();
  std::string filename  = vm[ "filename" ].as< std::string >();
  if( options.size() < 2 )
  {
    trace.error() << " At least 2 characters are required with option --export.";
    trace.info() << std::endl;
    exit(1);
  }

  //We check that the shape is known
  unsigned int id = checkAndReturnIndex(shapeName);

  // standard types
  typedef Z2i::Space::Point Point;
  typedef Z2i::Space::RealPoint RealPoint;

  RealPoint center( vm[ "center_x" ].as< double >(),
      vm[ "center_y" ].as< double >());

  if( id == 0 )
  {
    if( !( vm.count( "radius" )))
    {
      missingParam( "--radius" );
    }
    double radius = vm[ "radius" ].as< double >();
    Ball2D< Z2i::Space > ball( Z2i::Point( 0, 0 ), radius );
    for ( uint i_h = 0; i_h < rh_array.size (); ++i_h )
    {
      std::cout << "computation for h=" << rh_array[ i_h ] << std::endl;
      std::stringstream sstm;
      sstm << "ScaleSpaceMean2D_" << filename << "_h_" << rh_array[ i_h ] << ".pgm";
      std::string exportname = sstm.str();

      Compute( ball, re_array, rh_array[ i_h ], options, exportname.c_str() );
    }
  }
  else if( id == 1 )
  {
    if( !( vm.count( "width" )))
    {
      missingParam( "--width" );
    }
    double width = vm[ "width" ].as< double >();
    ImplicitHyperCube< Z2i::Space > object( Z2i::Point( 0, 0 ), width / 2.0 );
    trace.error() << "Not available.";
    trace.info() << std::endl;
  }
  else if( id == 2 )
  {
    if( !( vm.count( "power" )))
    {
      missingParam( "--power" );
    }
    if( !( vm.count( "radius" )))
    {
      missingParam( "--radius" );
    }
    double radius = vm[ "radius" ].as< double >();
    double power = vm[ "power" ].as< double >();
    ImplicitRoundedHyperCube< Z2i::Space > ball( Z2i::Point( 0, 0 ), radius, power );
    trace.error() << "Not available.";
    trace.info() << std::endl;
  }
  else if( id == 3 )
  {
    if( !( vm.count( "varsmallradius" )))
    {
      missingParam( "--varsmallradius" );
    }
    if( !( vm.count( "radius" )))
    {
      missingParam( "--radius" );
    }
    if( !( vm.count( "k" )))
    {
      missingParam( "--k" );
    }
    if( !( vm.count( "phi" )))
    {
      missingParam( "--phi" );
    }
    double radius = vm[ "radius" ].as< double >();
    double varsmallradius = vm[ "varsmallradius" ].as< double >();
    unsigned int k = vm[ "k" ].as< unsigned int >();
    double phi = vm[ "phi" ].as< double >();
    Flower2D< Z2i::Space > flower( center, radius, varsmallradius, k, phi );
    for ( uint i_h = 0; i_h < rh_array.size (); ++i_h )
    {
      std::cout << "computation for h=" << rh_array[ i_h ] << std::endl;
      std::stringstream sstm;
      sstm << "ScaleSpaceMean2D_" << filename << "_h_" << rh_array[ i_h ] << ".pgm";
      std::string exportname = sstm.str();

      Compute( flower, re_array, rh_array[ i_h ], options, exportname.c_str() );
    }
  }
  else if( id == 4 )
  {
    if( !( vm.count( "radius" )))
    {
      missingParam( "--radius" );
    }
    if( !( vm.count( "k" )))
    {
      missingParam( "--k" );
    }
    if( !( vm.count( "phi" )))
    {
      missingParam( "--phi" );
    }
    double radius = vm[ "radius" ].as< double >();
    unsigned int k = vm[ "k" ].as< unsigned int >();
    double phi = vm[ "phi" ].as< double >();
    NGon2D< Z2i::Space > object( center, radius, k, phi );
    for ( uint i_h = 0; i_h < rh_array.size (); ++i_h )
    {
      std::cout << "computation for h=" << rh_array[ i_h ] << std::endl;
      std::stringstream sstm;
      sstm << "ScaleSpaceMean2D_" << filename << "_h_" << rh_array[ i_h ] << ".pgm";
      std::string exportname = sstm.str();

      Compute( object, re_array, rh_array[ i_h ], options, exportname.c_str() );
    }
  }
  else if( id == 5 )
  {
    if( !( vm.count( "varsmallradius" )))
    {
      missingParam( "--varsmallradius" );
    }
    if( !( vm.count( "radius" )))
    {
      missingParam( "--radius" );
    }
    if( !( vm.count( "k" )))
    {
      missingParam( "--k" );
    }
    if( !( vm.count( "phi" )))
    {
      missingParam( "--phi" );
    }
    double radius = vm[ "radius" ].as< double >();
    double varsmallradius = vm[ "varsmallradius" ].as< double >();
    unsigned int k = vm[ "k" ].as< unsigned int >();
    double phi = vm[ "phi" ].as< double >();
    AccFlower2D< Z2i::Space > accflower( center, radius, varsmallradius, k, phi );
    for ( uint i_h = 0; i_h < rh_array.size (); ++i_h )
    {
      std::cout << "computation for h=" << rh_array[ i_h ] << std::endl;
      std::stringstream sstm;
      sstm << "ScaleSpaceMean2D_" << filename << "_h_" << rh_array[ i_h ] << ".pgm";
      std::string exportname = sstm.str();

      Compute( accflower, re_array, rh_array[ i_h ], options, exportname.c_str() );
    }
  }
  else if( id == 6 )
  {
    if( !( vm.count( "axis1" )))
    {
      missingParam( "--axis1" );
    }
    if( !( vm.count( "axis2" )))
    {
      missingParam( "--axis2" );
    }
    if( !( vm.count( "phi" )))
    {
      missingParam( "--phi" );
    }
    double a1 = vm[ "axis1" ].as< double >();
    double a2 = vm[ "axis2" ].as< double >();
    double phi = vm[ "phi" ].as< double >();
    Ellipse2D< Z2i::Space > ellipse( center, a1, a2, phi );
    for ( uint i_h = 0; i_h < rh_array.size (); ++i_h )
    {
      std::cout << "computation for h=" << rh_array[ i_h ] << std::endl;
      std::stringstream sstm;
      sstm << "ScaleSpaceMean2D_" << filename << "_h_" << rh_array[ i_h ];
      std::string exportname = sstm.str();

      Compute( ellipse, re_array, rh_array[ i_h ], options, exportname.c_str() );
    }
  }
  else if( id == 7 )
  {
    if( !( vm.count( "axis1" )))
    {
      missingParam( "--axis1" );
    }
    if( !( vm.count( "axis2" )))
    {
      missingParam( "--axis2" );
    }
    if( !( vm.count( "phi" )))
    {
      missingParam( "--phi" );
    }
    double a1 = vm[ "axis1" ].as< double >();
    double a2 = vm[ "axis2" ].as< double >();
    double phi = vm[ "phi" ].as< double >();
    Cublipse2D< Z2i::Space > cublipse( center, a1, a2, phi );
    for ( uint i_h = 0; i_h < rh_array.size (); ++i_h )
    {
      std::cout << "computation for h=" << rh_array[ i_h ] << std::endl;
      std::stringstream sstm;
      sstm << "ScaleSpaceMean2D_" << filename << "_h_" << rh_array[ i_h ];
      std::string exportname = sstm.str();

      Compute( cublipse, re_array, rh_array[ i_h ], options, exportname.c_str() );
    }
  }

  return 42;
}

