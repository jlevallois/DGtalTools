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
#include <vector>
#include <string>

#include "DGtal/base/Common.h"
#include "DGtal/base/Clock.h"

//shapes
#include "DGtal/shapes/ShapeFactory.h"
#include "DGtal/shapes/Shapes.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/topology/helpers/Surfaces.h"

#include "DGtal/topology/SetOfSurfels.h"
#include "DGtal/topology/DigitalSurface.h"
#include "DGtal/kernel/sets/SetPredicate.h"

//Digitizer
#include "DGtal/shapes/GaussDigitizer.h"
#include "DGtal/geometry/curves/GridCurve.h"
#include "DGtal/graph/DepthFirstVisitor.h"
#include "DGtal/graph/GraphVisitorRange.h"
#include "DGtal/geometry/volumes/KanungoNoise.h"


//Estimators
#include "DGtal/geometry/curves/estimation/TrueLocalEstimatorOnPoints.h"

#include "DGtal/geometry/surfaces/estimation/IntegralInvariantMeanCurvatureEstimator_0memory.h"
#include "DGtal/geometry/surfaces/estimation/IntegralInvariantGaussianCurvatureEstimator_0memory.h"
#include "DGtal/geometry/surfaces/FunctorOnCells.h"

#include "DGtal/topology/LightImplicitDigitalSurface.h"
#include "DGtal/images/imagesSetsUtils/SimpleThresholdForegroundPredicate.h"
#include "DGtal/math/MPolynomial.h"
#include "DGtal/io/readers/MPolynomialReader.h"
#include "DGtal/shapes/implicit/ImplicitPolynomial3Shape.h"


//Vol Export
#include "DGtal/io/colormaps/GradientColorMap.h"
#include "DGtal/io/colormaps/GrayscaleColorMap.h"
#include "DGtal/io/writers/VolWriter.h"
#include "DGtal/images/ImageHelper.h"

#include "DGtal/kernel/BasicPointFunctors.h"

using namespace DGtal;

template < typename Shape, typename KSpace, typename ConstIterator, typename OutputIterator >
void
estimateTrueMeanCurvatureQuantity( const ConstIterator & it_begin,
                                   const ConstIterator & it_end,
                                   OutputIterator & output,
                                   const KSpace & K,
                                   const double & h,
                                   Shape * aShape )
{
    typedef typename KSpace::SCell Spel;
    typedef typename KSpace::Point Point;
    typedef typename KSpace::Space::RealPoint RealPoint;

    Spel currentInnerSpel;
    RealPoint currentInnerRealPoint;

    for ( ConstIterator it = it_begin; it != it_end; ++it )
    {
        currentInnerSpel = K.sDirectIncident( *it_begin, K.sOrthDir( *it_begin ) ); /// Spel on the border, but inside the shape
        currentInnerRealPoint = ((RealPoint)(K.sCoords( currentInnerSpel ))) * h;
        output = aShape->meanCurvature( currentInnerRealPoint );
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
    typedef typename KSpace::SCell Spel;
    typedef typename KSpace::Point Point;
    typedef typename KSpace::Space::RealPoint RealPoint;

    Spel currentInnerSpel;
    RealPoint currentInnerRealPoint;

    for ( ConstIterator it = it_begin; it != it_end; ++it )
    {
        currentInnerSpel = K.sDirectIncident( *it_begin, K.sOrthDir( *it_begin ) ); /// Spel on the border, but inside the shape
        currentInnerRealPoint = ((RealPoint)(K.sCoords( currentInnerSpel ))) * h;
        output = aShape->gaussianCurvature( currentInnerRealPoint );
        ++output;
    }
}

template <typename Space, typename Shape>
bool
compareShapeEstimators( const std::string & filename,
                        Shape * aShape,
                        const double border_min[],
                        const double border_max[],
                        const double & h,
                        const double & radius_kernel,
                        const bool & lambda_optimized )
{
    typedef typename Space::RealPoint RealPoint;
    typedef GaussDigitizer< Z3i::Space, Shape > DigitalShape;
    typedef Z3i::KSpace KSpace;
    typedef LightImplicitDigitalSurface< KSpace, DigitalShape > Boundary;
    typedef Z3i::Domain Domain;
    typedef typename KSpace::Surfel Surfel;
    typedef PointFunctorFromPointPredicateAndDomain< DigitalShape, Domain, unsigned int > MyPointFunctor;
    typedef FunctorOnCells< MyPointFunctor, KSpace > MyCellFunctor;

    DigitalShape* dshape = new DigitalShape();
    dshape->attach( *aShape );
    dshape->init( RealPoint( -10.0, -10.0, -10.0 ), RealPoint( 10.0, 10.0, 10.0 ), h );

    std::cout << "h=" << h << std::endl;

    KSpace K;
    bool ok = K.init( dshape->getLowerBound(), dshape->getUpperBound(), true );
    if ( ! ok )
    {
        std::cerr << "[3dLocalEstimators_0memory]" << " error in creating KSpace." << std::endl;
        return false;
    }

    try
    {
        // Extracts shape boundary
        Surfel bel = Surfaces<KSpace>::findABel ( K, *dshape, 10000 );
        Boundary boundary ( K, *dshape, SurfelAdjacency< KSpace::dimension >( true ), bel );

        // Estimations

        // True Mean Curvature
        {
            std::stringstream ss;
            ss << filename << "_True_mean" << ".dat";
            std::cout << ss.str().c_str() << std::endl;
            std::ofstream file( ss.str().c_str() );
            file.flags( std::ios_base::unitbuf );
            file << "# h = " << h << std::endl;
            file << "# True Mean Curvature estimation" << std::endl;

            std::ostream_iterator< double > out_it( file, "\n" );

            estimateTrueMeanCurvatureQuantity( boundary.begin(),
                                               boundary.end(),
                                               out_it,
                                               K,
                                               h,
                                               aShape );
            file.close();
        }

        // True Gaussian Curvature
        {
            std::stringstream ss;
            ss << filename << "_True_gaussian" << ".dat";
            std::ofstream file( ss.str().c_str() );
            file.flags( std::ios_base::unitbuf );
            file << "# h = " << h << std::endl;
            file << "# True Gaussian Curvature estimation" << std::endl;

            std::ostream_iterator< double > out_it( file, "\n" );

            estimateTrueGaussianCurvatureQuantity( boundary.begin(),
                                                   boundary.end(),
                                                   out_it,
                                                   K,
                                                   h,
                                                   aShape );
            file.close();
        }

        MyPointFunctor pointFunctor( dshape, dshape->getDomain(), 1, 0 );
        MyCellFunctor functor ( pointFunctor, K );
        Clock c;

        double re_convolution_kernel = radius_kernel * std::pow( h, 1.0/3.0 ); // to obtains convergence results, re must follow the rule re=kh^(1/3)

        // Integral Invariant Mean Curvature
        {
            IntegralInvariantMeanCurvatureEstimator_0memory< KSpace, MyCellFunctor > IIMeanCurvatureEstimator ( K, functor );

            c.startClock();
            IIMeanCurvatureEstimator.init ( h, 3.0 );//re_convolution_kernel );

            std::stringstream ss;
            ss << filename << "_II_mean" << ".dat";
            std::ofstream file( ss.str().c_str() );
            file.flags( std::ios_base::unitbuf );
            file << "# h = " << h << std::endl;
            file << "# Mean Curvature estimation from the Integral Invariant" << std::endl;
            file << "# computed kernel radius = " << re_convolution_kernel << std::endl;

            std::ostream_iterator< double > out_it( file, "\n" );
            /*if( !lambda_optimized )
            {
                IIMeanCurvatureEstimator.eval( boundary.begin(), boundary.end(), out_it );
            }
            else*/
            {
                IIMeanCurvatureEstimator.eval( boundary.begin(), boundary.end(), out_it, *aShape );
            }
            double TIIMeanCurv = c.stopClock();
            file << "# time = " << TIIMeanCurv << std::endl;
            file.close();
        }

        // Integral Invariant Gaussian Curvature

        /*{
            IntegralInvariantGaussianCurvatureEstimator_0memory< KSpace, MyCellFunctor > IIGaussianCurvatureEstimator ( K, functor );

            c.startClock();
            IIGaussianCurvatureEstimator.init( h, 3.0 );

            std::stringstream ss;
            ss << filename << "_II_gaussian" << ".dat";
            std::ofstream file( ss.str().c_str() );
            file.flags( std::ios_base::unitbuf );
            file << "# h = " << h << std::endl;
            file << "# Gaussian Curvature estimation from the Integral Invariant" << std::endl;
            file << "# computed kernel radius = " << re_convolution_kernel << std::endl;

            std::ostream_iterator< double > out_it( file, "\n" );


            if( !lambda_optimized )
            {
                IIGaussianCurvatureEstimator.eval ( boundary.begin(), boundary.end(), out_it );
            }
            else
            {
                IIGaussianCurvatureEstimator.eval ( boundary.begin(), boundary.end(), out_it, *aShape );
            }
            double TIIGaussCurv = c.stopClock();
            file << "# time = " << TIIGaussCurv << std::endl;
            file.close();
        }*/
    }
    catch ( InputException e )
    {
        std::cerr << "[estimatorCurvatureComparator3D]"
                  << " error."
                  << e.what() << std::endl;
        return false;
    }
    return true;
}

void usage( int /*argc*/, char** argv )
{
    std::cerr << "Usage: " << argv[ 0 ] << " <filename_export> <Polynomial> <Px> <Py> <Pz> <Qx> <Qy> <Qz> <step> <radius> <use optimised lambda (0 or 1)> (<path_to_save>) (<name_of_vol_file>)" << std::endl;
    std::cerr << "\t - displays the boundary of a shape defined implicitly by a 3-polynomial <Polynomial>." << std::endl;
    std::cerr << "\t - P and Q defines the bounding box." << std::endl;
    std::cerr << "\t - step is the grid step." << std::endl;
    std::cerr << "\t - radius is the kernel support k radius." << std::endl;
    std::cerr << "\t - path is optional. It's the path where you want to generate a .vol file of the polynomial shape (for external computation). If no path was set, we don't export as a .vol file." << std::endl;
    std::cerr << "\t - name is optional. It's the name of your .vol file you want to generate (for external computation). If no name was set, we don't export as a .vol file." << std::endl;
}

int main( int argc, char** argv )
{
    if ( argc < 12 )
    {
        usage( argc, argv );
        return 1;
    }
    std::string file_export = argv[ 1 ];

    double border_min[ 3 ];
    double border_max[ 3 ];
    for ( unsigned int i = 0; i < 3; ++i )
    {
        border_min[ i ] = atof( argv[ 3+i ] );
        border_max[ i ] = atof( argv[ 6+i ] );
    }
    double h = atof( argv[ 9 ] );
    double radius = atof( argv[ 10 ] );

    bool lambda_optimized = atoi( argv[ 11 ] );
    ASSERT(( lambda_optimized == 0 || lambda_optimized == 1 ));

    typedef Z3i::Space::RealPoint RealPoint;
    typedef Z3i::Space::RealPoint::Coordinate Ring;

    /// Construction of the polynomial shape

    std::string poly_str = "x^2+y^2+z^2-25";//argv[ 2 ];
    {
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

        compareShapeEstimators< Z3i::Space, ImplicitShape > (
                    file_export,
                    shape,
                    border_min, border_max,
                    h,
                    radius,
                    lambda_optimized );

        delete shape;
        shape = NULL;
    }


    //    bool export_vol = false;
    //    std::string pathToSaveVolFile = "";
    //    std::string nameVolFile = "";
    //    if( argc >= 12 )
    //    {
    //        export_vol = true;
    //        pathToSaveVolFile = argv[ 10 ];
    //        nameVolFile = argv[ 11 ];
    //    }

    //    ImplicitShape shape( poly );

    //    typedef ImplicitBall< Z3i::Space > ImpBall;
    //    ImpBall ball( Z3i::RealPoint(0,0,0), 5.0 );

    //    /// Computation of 3D curvature estimators (mean & Gaussian) on the implicit shape
    //    compareShapeEstimators< Z3i::Space, ImpBall >
    //            (
    //                poly_str,
    //                ball,//shape,
    //                border_min, border_max,
    //                h,
    //                radius,
    //                lambda_optimized );/*,
    //                export_vol, pathToSaveVolFile, nameVolFile
    //                );*/
}
