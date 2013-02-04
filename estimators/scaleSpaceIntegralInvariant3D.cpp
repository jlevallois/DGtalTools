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
 * @file scaleSpaceIntegralInvariant3D.cpp
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
#include "DGtal/math/MPolynomial.h"
#include "DGtal/io/readers/MPolynomialReader.h"
#include "DGtal/shapes/implicit/ImplicitPolynomial3Shape.h"

//Digitizer
#include "DGtal/topology/LightImplicitDigitalSurface.h"
#include "DGtal/topology/DigitalSurface.h"
#include "DGtal/topology/DepthFirstVisitor.h"


//Estimator
#include "DGtal/geometry/surfaces/estimation/IntegralInvariantMeanCurvatureEstimator.h"
#include "DGtal/geometry/surfaces/estimation/IntegralInvariantGaussianCurvatureEstimator.h"
#include "DGtal/geometry/surfaces/FunctorOnCells.h"


//Image & Viewer3D
#include "DGtal/images/ImageSelector.h"
#include "DGtal/io/colormaps/GradientColorMap.h"
#include "DGtal/io/writers/PPMWriter.h"

#define uint unsigned int

using namespace DGtal;

typedef ImplicitPolynomial3Shape< Z3i::Space > ImplicitShape;
typedef ImageSelector< Z2i::Domain, float >::Type Image;

Image* Compute( ImplicitShape & shape, std::vector< double > & re_array,
             double * border_min, double * border_max, double h,
             std::string mode, std::string name )
{
    typedef Z3i::KSpace::Surfel Surfel;
    typedef double Quantity;
    typedef GaussDigitizer< Z3i::Space, ImplicitShape > MyGaussDigitizer;
    typedef FunctorOnCells< MyGaussDigitizer, Z3i::KSpace > MyFunctor;
    typedef typename MyGaussDigitizer::RealPoint RealPoint;
    typedef LightImplicitDigitalSurface< Z3i::KSpace, MyGaussDigitizer > MyLightImplicitDigitalSurface;
    typedef DigitalSurface< MyLightImplicitDigitalSurface > MyDigitalSurface;
    typedef DepthFirstVisitor< MyDigitalSurface > Visitor;
    typedef typename Visitor::VertexConstIterator SurfelConstIterator;
    typedef GradientColorMap< Quantity > Gradient;

    Image * image;

    MyGaussDigitizer gaussDigShape;
    gaussDigShape.attach( shape );
    gaussDigShape.init( RealPoint( border_min ), RealPoint( border_max ), h );

    Z3i::Domain domain = gaussDigShape.getDomain();
    Z3i::KSpace kSpace;
    bool space_ok = kSpace.init( domain.lowerBound(), domain.upperBound(), true );
    if ( !space_ok )
    {
        trace.error() << "Error in the Khalimsky space construction."<<std::endl;
        return image;
    }

    SurfelAdjacency< Z3i::KSpace::dimension > SAdj( true );
    Surfel bel = Surfaces< Z3i::KSpace >::findABel( kSpace, gaussDigShape, 100000 );
    MyLightImplicitDigitalSurface lightImplDigSurf( kSpace, gaussDigShape, SAdj, bel );
    MyDigitalSurface digSurfShape( lightImplDigSurf );

    uint re_size = re_array.size ();

    Quantity min = numeric_limits< Quantity >::max();
    Quantity max = numeric_limits <Quantity >::min();

    for ( uint step_re = 0; step_re < re_size; ++step_re )
    {
        MyFunctor functor ( gaussDigShape, kSpace, domain, true );

        Visitor *depth = new Visitor( digSurfShape, *digSurfShape.begin() );
        SurfelConstIterator abegin = SurfelConstIterator( depth );
        SurfelConstIterator aend = SurfelConstIterator( 0 );

        std::vector< Quantity > results;
        back_insert_iterator< std::vector< Quantity > > resultsIterator( results );

        Clock c;
        if ( mode == "mean" )
        {
            typedef IntegralInvariantMeanCurvatureEstimator< Z3i::KSpace, MyFunctor > MyIIMeanEstimator;

            MyIIMeanEstimator estimator ( kSpace, functor );

            c.startClock();

            estimator.init( h, re_array[ step_re ] );
            estimator.eval ( abegin, aend, resultsIterator );
        }
        else
        {
            typedef IntegralInvariantGaussianCurvatureEstimator< Z3i::KSpace, MyFunctor > MyIIGaussianEstimator;

            MyIIGaussianEstimator estimator ( kSpace, functor );

            c.startClock();

            estimator.init( h, re_array[ step_re ] );
            estimator.eval ( abegin, aend, resultsIterator );
        }

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

void usage( int argc, char** argv )
{
    std::cerr << "Usage: " << argv[ 0 ] << " <Polynomial> <Px> <Py> <Pz> <Qx> <Qy> <Qz> <min_re> <max_re> <step_re> <min_h> <max_h> <step_h> <\"mean\" || \"gaussian\"> <filename>" << std::endl;
    std::cerr << "\t - displays the boundary of a shape defined implicitly by a 3-polynomial <Polynomial>." << std::endl;
    std::cerr << "\t - P and Q defines the bounding box." << std::endl;
    std::cerr << "\t - <min_re> is the min Euclidean radius of the convolution kernel computed." << std::endl;
    std::cerr << "\t - <max_re> is the max Euclidean radius of the convolution kernel computed." << std::endl;
    std::cerr << "\t - <step_re> is the step to increase the Euclidean radius of the convolution kernel from min to max." << std::endl;
    std::cerr << "\t - <min_h> is the min grid step computed." << std::endl;
    std::cerr << "\t - <max_h> is the max grid step computed." << std::endl;
    std::cerr << "\t - <step_h> is the step to increase the grid step from min to max." << std::endl;
    std::cerr << "\t - <\"mean\" || \"gaussian\"> compute mean or Gaussian curvature.";
    std::cerr << "\t - <filename> name of the pgm image created with the scale space analysis.";
    std::cerr << "\t - Example: " << argv[ 0 ] << "\"x^4 + y^4 + z^4 - 0.25\" -1 -1 -1 1 1 1 0.05 0.4 0.005 0.01 0.01 0.01 mean roundedcube" << std::endl;
}

int main ( int argc, char** argv )
{
    if ( argc < 16 )
    {
        usage( argc, argv );
        return 0;
    }

    std::string mode = argv[ 14 ];
    if ( mode != "gaussian" && mode != "mean" )
    {
        usage( argc, argv );
        return 0;
    }

    double border_min[ 3 ];
    double border_max[ 3 ];
    for ( uint i = 0; i < 3; ++i )
    {
        border_min[ i ] = atof( argv[ 2+i ] );
        border_max[ i ] = atof( argv[ 5+i ] );
    }

    typedef Z3i::Space::RealPoint::Coordinate Ring;
    typedef MPolynomial< 3, Ring > Polynomial3;
    typedef MPolynomialReader< 3, Ring > Polynomial3Reader;

    Polynomial3 poly;
    Polynomial3Reader reader;
    std::string poly_str = argv[ 1 ];
    std::string::const_iterator iter = reader.read( poly, poly_str.begin(), poly_str.end() );
    if ( iter != poly_str.end() )
    {
        std::cerr << "ERROR: I read only <"
                  << poly_str.substr( 0, iter - poly_str.begin() )
                  << ">, and I built P=" << poly << std::endl;
        return 1;
    }
    ImplicitShape shape( poly );

    std::vector< double > rh_array;
    double rh_init = atof( argv[ 11 ] ); // 0.1
    double rh_end  = atof( argv[ 12 ] ); // 0.05
    double rh_step = atof( argv[ 13 ] ); // 0.01
    for ( ; rh_init >= rh_end; rh_init -= rh_step )
    {
        rh_array.push_back ( rh_init );
    };

    std::vector< double > re_array;
    double re_init = atof( argv[ 8 ] ); // 0.02
    double re_end  = atof( argv[ 9 ] ); // 0.06
    double re_step = atof( argv[ 10 ] ); // 0.0005
    for ( ; re_init <= re_end; re_init += re_step )
    {
        re_array.push_back ( re_init );
    }

    for ( uint i_h = 0; i_h < rh_array.size (); ++i_h )
    {
        std::cout << "computation for h=" << rh_array[ i_h ] << std::endl;

        std::stringstream sstm;
        sstm << "ScaleSpaceMean3D_" << argv[ 15 ] << "_h_" << rh_array[ i_h ] << ".pgm";
        std::string exportname = sstm.str();
        Compute( shape, re_array, border_min, border_max, rh_array[ i_h ], mode, exportname.c_str() );
    }

    return 42;
}

