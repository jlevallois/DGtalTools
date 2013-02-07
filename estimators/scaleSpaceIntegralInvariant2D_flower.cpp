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
 * @file scaleSpaceIntegralInvariant2D.cpp
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

#include "DGtal/base/Common.h"
#include "DGtal/base/Clock.h"

//shapes
#include "DGtal/shapes/parametric/Flower2D.h"

//Digitizer
#include "DGtal/topology/LightImplicitDigitalSurface.h"
#include "DGtal/topology/DigitalSurface.h"
#include "DGtal/topology/DepthFirstVisitor.h"


//Estimator
#include "DGtal/geometry/surfaces/estimation/IntegralInvariantMeanCurvatureEstimator.h"
#include "DGtal/geometry/surfaces/FunctorOnCells.h"


//Image & Viewer3D
#include "DGtal/images/ImageSelector.h"
#include "DGtal/io/colormaps/GradientColorMap.h"
#include "DGtal/io/writers/PPMWriter.h"

#define uint unsigned int

using namespace DGtal;

typedef Flower2D< Z2i::Space > MyShape;
typedef GaussDigitizer< Z2i::Space, MyShape > MyGaussDigitizer;
typedef FunctorOnCells< MyGaussDigitizer, Z2i::KSpace > MyFunctor;
typedef IntegralInvariantMeanCurvatureEstimator< Z2i::KSpace, MyFunctor > MyMeanCurvatureEstimator;
typedef typename MyMeanCurvatureEstimator::Quantity Quantity;
typedef ImageSelector< Z2i::Domain, float >::Type Image;

Image* Compute( MyShape & shape, std::vector< double > & re_array,
             double h,
             std::string name )
{
    Image * image;
    typedef Z2i::KSpace::Surfel Surfel;

    typedef typename MyGaussDigitizer::RealPoint RealPoint;
    MyGaussDigitizer gaussDigShape;
    gaussDigShape.attach( shape );
    gaussDigShape.init( shape.getLowerBound(), shape.getUpperBound(), h );

    Z2i::Domain domain = gaussDigShape.getDomain();
    Z2i::KSpace kSpace;
    bool space_ok = kSpace.init( domain.lowerBound(), domain.upperBound(), true );
    if ( !space_ok )
    {
        trace.error() << "Error in the Khalimsky space construction."<<std::endl;
        return image;
    }

    typedef LightImplicitDigitalSurface< Z2i::KSpace, MyGaussDigitizer > MyLightImplicitDigitalSurface;
    typedef DigitalSurface< MyLightImplicitDigitalSurface > MyDigitalSurface;

    SurfelAdjacency< Z2i::KSpace::dimension > SAdj( true );
    Surfel bel = Surfaces< Z2i::KSpace >::findABel( kSpace, gaussDigShape, 100000 );
    MyLightImplicitDigitalSurface lightImplDigSurf( kSpace, gaussDigShape, SAdj, bel );
    MyDigitalSurface digSurfShape( lightImplDigSurf );


    typedef DepthFirstVisitor< MyDigitalSurface > Visitor;
    typedef typename Visitor::VertexConstIterator SurfelConstIterator;

    typedef std::vector< Quantity > vQuantity;
    typedef typename vQuantity::const_iterator const_interatorQuantity;

    uint re_size = re_array.size ();

    Quantity min = numeric_limits< Quantity >::max();
    Quantity max = numeric_limits <Quantity >::min();

    for ( uint step_re = 0; step_re < re_size; ++step_re )
    {
        MyFunctor functor ( gaussDigShape, kSpace, domain, true );
        MyMeanCurvatureEstimator meanCurvatureEstimator ( kSpace, functor );
        //double k = re_array[ step_re ] / ( std::pow ( h, ( 4.0/3.0 )));
        //std::cout << "computation for re=" << re_array[ step_re ] << " k=" << k << std::endl;

        Clock c;
        c.startClock();

        meanCurvatureEstimator.init( h, re_array[ step_re ] );


        vQuantity results;
        back_insert_iterator< vQuantity > resultIterator( results );
        Visitor *depth = new Visitor( digSurfShape, *digSurfShape.begin() );
        SurfelConstIterator abegin = SurfelConstIterator( depth );
        SurfelConstIterator aend = SurfelConstIterator( 0 );
        meanCurvatureEstimator.eval( abegin, aend, resultIterator );

        uint rsize = results.size();

        if ( step_re == 0 ) // compute image domain
        {
            Z2i::Point aa( 0, 0 );
            Z2i::Point bb( rsize - 1, re_size - 1 );
            Z2i::Domain domain_image( aa, bb );
            image = new Image( domain_image );
        }

        Quantity current_result = Quantity(0);
        for ( uint i = 0; i < rsize; ++i )
        {
            current_result = results[ i ];
            if ( current_result < min )
            {
                min = current_result;
            }
            else if ( current_result > max )
            {
                max = current_result;
            }

            image->operator []( step_re * rsize + i ) = current_result;
        }

        typedef GradientColorMap< Quantity > Gradient;
        Gradient cmap_grad( min, max );
        cmap_grad.addColor( Color( 50, 50, 255 ) );
        cmap_grad.addColor( Color( 255, 0, 0 ) );
        cmap_grad.addColor( Color( 255, 255, 10 ) );

        PPMWriter<Image, Gradient >::exportPPM( name, *image, cmap_grad ); // Erase previous image with old + new iteration values (the computation on full iteration take long time, so we can see result before ending)

        trace.progressBar( step_re + 1, re_size );
        double time = c.stopClock();
        std::cout << "-- Done in " << time << " msec. --" << std::endl;
    }

    return image;
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


void usage( int argc, char** argv )
{
    std::cerr << "Usage: " << argv[ 0 ] << " <min_re> <max_re> <step_re> <min_h> <max_h> <step_h> <filename>" << std::endl;
    std::cerr << "\t - <min_re> is the min Euclidean radius of the convolution kernel computed." << std::endl;
    std::cerr << "\t - <max_re> is the max Euclidean radius of the convolution kernel computed." << std::endl;
    std::cerr << "\t - <step_re> is the step to increase the Euclidean radius of the convolution kernel from min to max." << std::endl;
    std::cerr << "\t - <min_h> is the min grid step computed." << std::endl;
    std::cerr << "\t - <max_h> is the max grid step computed." << std::endl;
    std::cerr << "\t - <step_h> is the step to increase the grid step from min to max." << std::endl;
    std::cerr << "\t - <\"mean\" || \"gaussian\"> compute mean or Gaussian curvature.";
    std::cerr << "\t - <filename> name of the pgm image created with the scale space analysis.";
    std::cerr << "\t - Example: " << argv[ 0 ] << "1.0 7.0 0.01 0.05 0.05 0.01 flower" << std::endl;
}

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

            ("re,r",  po::value<std::string>()->default_value("1.0 7.0 0.01"), "Grid step for the digitization" )
            ("gridstep,g",  po::value<std::string>()->default_value("0.05 0.05 0.01"), "Grid step for the digitization" )
            ("export,e",  po::value<std::string>()->default_value("11 filename"), "the i-th estimator is disabled iff there is a 0 at position i" );

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
    if( !parseOK || vm.count("help") || argc<=1)
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


    if ( argc < 8 )
    {
        usage( argc, argv );
        return 0;
    }

    double max_radius_shape = 20.00217;
    double min_radius_shape = 10.00217;
    MyShape shape( 0, 0, max_radius_shape, min_radius_shape, 4, 3.0 ); //Flower

    std::vector< double > rh_array;
    double rh_init = atof( argv[ 4 ] ); // 0.1
    double rh_end  = atof( argv[ 5 ] ); // 0.05
    double rh_step = atof( argv[ 6 ] ); // 0.01
    for ( ; rh_init >= rh_end; rh_init -= rh_step )
    {
        rh_array.push_back ( rh_init );
    };

    std::vector< double > re_array;
    double re_init = atof( argv[ 1 ] ); // 0.02
    double re_end  = atof( argv[ 2 ] ); // 0.06
    double re_step = atof( argv[ 3 ] ); // 0.0005
    for ( ; re_init <= re_end; re_init += re_step )
    {
        re_array.push_back ( re_init );
    }

    int nb = 2; //number of available methods
    std::string options = vm["export"].as<std::string>();
    if (options.size() < 2)
      {
        trace.error() << " At least " << nb
              << " characters are required "
              << " with option --export.";
        trace.info() << std::endl;
        exit(1);
      }

    for ( uint i_h = 0; i_h < rh_array.size (); ++i_h )
    {
        std::cout << "computation for h=" << rh_array[ i_h ] << std::endl;

        std::stringstream sstm;
        sstm << "ScaleSpaceMean2D_" << argv[ 7 ] << "_h_" << rh_array[ i_h ] << ".pgm";
        std::string exportname = sstm.str();
        Compute( shape, re_array, rh_array[ i_h ], exportname.c_str() );
    }

    return 42;
}

