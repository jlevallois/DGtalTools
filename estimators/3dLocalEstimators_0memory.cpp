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
#include "DGtal/helpers/StdDefs.h"

//shapes
#include "DGtal/shapes/implicit/ImplicitBall.h"
#include "DGtal/math/MPolynomial.h"
#include "DGtal/io/readers/MPolynomialReader.h"
#include "DGtal/shapes/implicit/ImplicitPolynomial3Shape.h"

//Digitizer
#include "DGtal/shapes/GaussDigitizer.h"
//#include "DGtal/geometry/volumes/KanungoNoise.h"
#include "DGtal/topology/LightImplicitDigitalSurface.h"
#include "DGtal/topology/DigitalSurface.h"



//Estimators
#include "DGtal/kernel/BasicPointFunctors.h"
#include "DGtal/geometry/surfaces/FunctorOnCells.h"

#include "DGtal/geometry/curves/estimation/TrueLocalEstimatorOnPoints.h"
#include "DGtal/geometry/surfaces/estimation/IntegralInvariantMeanCurvatureEstimator_0memory.h"
#include "DGtal/geometry/surfaces/estimation/IntegralInvariantGaussianCurvatureEstimator_0memory.h"

#include "DGtal/geometry/surfaces/estimation/LocalEstimatorFromSurfelFunctorAdapter.h"
#include "DGtal/geometry/surfaces/estimation/estimationFunctors/MongeJetFittingGaussianCurvatureEstimator.h"
#include "DGtal/geometry/surfaces/estimation/estimationFunctors/MongeJetFittingMeanCurvatureEstimator.h"

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
                        const Shape * aShape,
                        const double border_min[],
                        const double border_max[],
                        const double & h,
                        const double & radius_kernel,
                        const bool & lambda_optimized )
{
    ////////
    typedef Z3i::Space::RealPoint::Coordinate Ring;
    typedef MPolynomial< 3, Ring > Polynomial3;
    typedef MPolynomialReader<3, Ring> Polynomial3Reader;
    typedef ImplicitPolynomial3Shape<Z3i::Space> ImplicitShape;
    ////////

    typedef typename Space::RealPoint RealPoint;
    typedef GaussDigitizer< Z3i::Space, Shape > DigitalShape;
    typedef Z3i::KSpace KSpace;
    typedef LightImplicitDigitalSurface< KSpace, DigitalShape > Boundary;
    typedef DigitalSurface< Boundary > MyDigitalSurface;
    typedef typename MyDigitalSurface::ConstIterator ConstIterator;

    typedef Z3i::Domain Domain;
    typedef typename KSpace::Surfel Surfel;
    typedef PointFunctorFromPointPredicateAndDomain< DigitalShape, Domain, unsigned int > MyPointFunctor;
    typedef FunctorOnCells< MyPointFunctor, KSpace > MyCellFunctor;

    DigitalShape* dshape = new DigitalShape();
    dshape->attach( *aShape );
    dshape->init( RealPoint( border_min[ 0 ], border_min[ 1 ], border_min[ 2 ] ), RealPoint( border_max[ 0 ], border_max[ 1 ], border_max[ 2 ] ), h );

    KSpace K;
    if ( ! K.init( dshape->getLowerBound(), dshape->getUpperBound(), true ) )
    {
        std::cerr << "[3dLocalEstimators_0memory]" << " error in creating KSpace." << std::endl;
        return false;
    }

    try
    {
        // Extracts shape boundary
        Surfel bel = Surfaces<KSpace>::findABel ( K, *dshape, 10000 );
        Boundary boundary ( K, *dshape, SurfelAdjacency< KSpace::dimension >( true ), bel );
        MyDigitalSurface surf ( boundary );

        // Estimations

        // True Mean Curvature

        std::stringstream ss;
        ss << filename << "_True_mean" << ".dat";
        std::ofstream file( ss.str().c_str() );
        file.flags( std::ios_base::unitbuf );
        file << "# h = " << h << std::endl;
        file << "# True Mean Curvature estimation" << std::endl;

        std::ostream_iterator< double > out_it_true_mean( file, "\n" );

        estimateTrueMeanCurvatureQuantity( boundary.begin(),
                                           boundary.end(),
                                           out_it_true_mean,
                                           K,
                                           h,
                                           aShape );
        file.close();


        // True Gaussian Curvature

        ss.str( std::string() );
        ss << filename << "_True_gaussian" << ".dat";
        file.open( ss.str().c_str() );
        file.flags( std::ios_base::unitbuf );
        file << "# h = " << h << std::endl;
        file << "# True Gaussian Curvature estimation" << std::endl;

        std::ostream_iterator< double > out_it_true_gaussian( file, "\n" );

        estimateTrueGaussianCurvatureQuantity( boundary.begin(),
                                               boundary.end(),
                                               out_it_true_gaussian,
                                               K,
                                               h,
                                               aShape );
        file.close();



        ss.str( std::string() );
        ss << filename << "_II_mean" << ".dat";
        file.open( ss.str().c_str() );
        file.flags( std::ios_base::unitbuf );
        file << "# h = " << h << std::endl;
        file << "# Mean Curvature estimation from the Integral Invariant" << std::endl;
        double re_convolution_kernel = radius_kernel * std::pow( h, 1.0/3.0 ); // to obtains convergence results, re must follow the rule re=kh^(1/3)
        file << "# computed kernel radius = " << re_convolution_kernel << std::endl;

        MyPointFunctor pointFunctor( dshape, dshape->getDomain(), 1, 0 );
        MyCellFunctor functor ( pointFunctor, K );
        Clock c;


        // Integral Invariant Mean Curvature
        IntegralInvariantMeanCurvatureEstimator_0memory< KSpace, MyCellFunctor > * IIMeanCurvatureEstimator = new IntegralInvariantMeanCurvatureEstimator_0memory< KSpace, MyCellFunctor >( K, functor );
        IIMeanCurvatureEstimator->init ( h, re_convolution_kernel );

        c.startClock();


        std::ostream_iterator< double > out_it_ii_mean( file, "\n" );
        if( !lambda_optimized )
        {
            IIMeanCurvatureEstimator->eval( boundary.begin(), boundary.end(), out_it_ii_mean );
        }
        else
        {
            IIMeanCurvatureEstimator->eval( boundary.begin(), boundary.end(), out_it_ii_mean, *aShape );
        }

        double TIIMeanCurv = c.stopClock();
        file << "# time = " << TIIMeanCurv << std::endl;
        file.close();

        delete IIMeanCurvatureEstimator;

        // Integral Invariant Gaussian Curvature
        IntegralInvariantGaussianCurvatureEstimator_0memory< KSpace, MyCellFunctor > * IIGaussianCurvatureEstimator = new IntegralInvariantGaussianCurvatureEstimator_0memory< KSpace, MyCellFunctor >( K, functor );

        c.startClock();
        IIGaussianCurvatureEstimator->init( h, re_convolution_kernel );

        ss.str( std::string() );
        ss << filename << "_II_gaussian" << ".dat";
        file.open( ss.str().c_str() );
        file.flags( std::ios_base::unitbuf );
        file << "# h = " << h << std::endl;
        file << "# Gaussian Curvature estimation from the Integral Invariant" << std::endl;
        file << "# computed kernel radius = " << re_convolution_kernel << std::endl;

        std::ostream_iterator< double > out_it_ii_gaussian( file, "\n" );


        if( !lambda_optimized )
        {
            IIGaussianCurvatureEstimator->eval ( boundary.begin(), boundary.end(), out_it_ii_gaussian );
        }
        else
        {
            IIGaussianCurvatureEstimator->eval ( boundary.begin(), boundary.end(), out_it_ii_gaussian, *aShape );
        }
        double TIIGaussCurv = c.stopClock();
        file << "# time = " << TIIGaussCurv << std::endl;
        file.close();

        delete IIGaussianCurvatureEstimator;


        // Monge Gaussian Curvature
        typedef MongeJetFittingGaussianCurvatureEstimator<Surfel, CanonicSCellEmbedder<KSpace> > FunctorGaussian;
        typedef LocalEstimatorFromSurfelFunctorAdapter<typename MyDigitalSurface::DigitalSurfaceContainer, Z3i::L2Metric, FunctorGaussian> ReporterK;
        FunctorGaussian estimatorK( (CanonicSCellEmbedder<KSpace>( K )), h );
        ReporterK reporterK(surf.container(), Z3i::l2Metric, estimatorK);
        c.startClock();
        reporterK.init( h , re_convolution_kernel / h  );

        typename ReporterK::SurfelConstIterator aaabegin = surf.container().begin();
        typename ReporterK::SurfelConstIterator aaaend = surf.container().end();

        ss.str( std::string() );
        ss << filename << "_MongeJetFitting_gaussian" << ".dat";
        file.open( ss.str().c_str() );
        file.flags( std::ios_base::unitbuf );
        file << "# h = " << h << std::endl;
        file << "# Gaussian Curvature estimation from CGAL Monge from and Jet Fitting" << std::endl;
        file << "# computed kernel radius = " << re_convolution_kernel << std::endl;
        std::ostream_iterator< double > out_it_monge_gaussian( file, "\n" );
        reporterK.eval(aaabegin, aaaend , out_it_monge_gaussian);
        double TMongeGaussCurv = c.stopClock();
        file << "# time = " << TMongeGaussCurv << std::endl;
        file.close();


        // Monge Mean Curvature
        typedef MongeJetFittingMeanCurvatureEstimator<Surfel, CanonicSCellEmbedder<KSpace> > FunctorMean;
        typedef LocalEstimatorFromSurfelFunctorAdapter<typename MyDigitalSurface::DigitalSurfaceContainer, Z3i::L2Metric, FunctorMean> ReporterH;
        FunctorMean estimatorH( (CanonicSCellEmbedder<KSpace>( K )), h );
        ReporterH reporterH(surf.container(), Z3i::l2Metric, estimatorH);
        c.startClock();
        reporterH.init( h , re_convolution_kernel / h  );

        ss.str( std::string() );
        ss << filename << "_MongeJetFitting_mean" << ".dat";
        file.open( ss.str().c_str() );
        file.flags( std::ios_base::unitbuf );
        file << "# h = " << h << std::endl;
        file << "# Mean Curvature estimation from CGAL Monge from and Jet Fitting" << std::endl;
        file << "# computed kernel radius = " << re_convolution_kernel << std::endl;
        std::ostream_iterator< double > out_it_monge_mean( file, "\n" );

        typename ReporterK::SurfelConstIterator aabegin = surf.container().begin();
        typename ReporterK::SurfelConstIterator aaend = surf.container().end();
        reporterK.eval(aabegin, aaend , out_it_monge_mean);
        double TMongeMeanCurv = c.stopClock();
        file << "# time = " << TMongeMeanCurv << std::endl;
        file.close();


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
    std::cerr << "Usage: " << argv[ 0 ] << " <filename_export> <Polynomial> <Px> <Py> <Pz> <Qx> <Qy> <Qz> <step> <radius> <use optimised lambda (0 or 1)>" << std::endl;
    std::cerr << "\t - displays the boundary of a shape defined implicitly by a 3-polynomial <Polynomial>." << std::endl;
    std::cerr << "\t - P and Q defines the bounding box." << std::endl;
    std::cerr << "\t - step is the grid step." << std::endl;
    std::cerr << "\t - radius is the kernel support k radius." << std::endl;
    std::cerr << "\t - path is optional. It's the path where you want to generate a .vol file of the polynomial shape (for external computation). If no path was set, we don't export as a .vol file." << std::endl;
    std::cerr << "\t - name is optional. It's the name of your .vol file you want to generate (for external computation). If no name was set, we don't export as a .vol file." << std::endl;
}

int main( int argc, char** argv )
{


#ifndef WITH_CGAL
#error You need to have activated CGAL (WITH_CGAL) to include this file.
#endif
#ifndef WITH_EIGEN
#error You need to have activated EIGEN (WITH_EIGEN) to include this file.
#endif

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

    std::string poly_str = argv[ 2 ];

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
}
