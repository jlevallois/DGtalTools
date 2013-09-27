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
#include "DGtal/base/Clock.h"
#include "DGtal/helpers/StdDefs.h"

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

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
#include "DGtal/graph/DepthFirstVisitor.h"
#include "DGtal/geometry/volumes/KanungoNoise.h"



//Estimators
#include "DGtal/kernel/BasicPointFunctors.h"
#include "DGtal/geometry/surfaces/FunctorOnCells.h"

#include "DGtal/geometry/curves/estimation/TrueLocalEstimatorOnPoints.h"
#include "DGtal/geometry/surfaces/estimation/IntegralInvariantMeanCurvatureEstimator.h"
#include "DGtal/geometry/surfaces/estimation/IntegralInvariantGaussianCurvatureEstimator.h"
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
    typedef typename KSpace::Space::RealPoint RealPoint;
    typedef SCellToMidPoint< KSpace > Embedder;

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
    typedef SCellToMidPoint< KSpace > Embedder;

    Embedder embedder( K );
    RealPoint currentRealPoint;

    for ( ConstIterator it = it_begin; it != it_end; ++it )
    {
        currentRealPoint = embedder( *it_begin ) * h;
        output = aShape->gaussianCurvature( currentRealPoint );
        ++output;
    }
}

/*template< typename TShape, typename TPredicat >
class ExtractDigitalSurface
{
    typedef TShape Shape;
    typedef TPredicat Predicat;
public:
    ExtractDigitalSurface( Shape * aShape )
    {

    }
};*/

template <typename Space, typename Shape>
bool
compareShapeEstimators( const std::string & filename,
                        const Shape * aShape,
                        const typename Space::RealPoint & border_min,
                        const typename Space::RealPoint & border_max,
                        const double & h,
                        const double & radius_kernel,
                        const double & alpha,
                        const std::string & options,
                        const bool & lambda_optimized,
                        double noiseLevel = 0.0 )
{
    typedef typename Space::RealPoint RealPoint;
    typedef GaussDigitizer< Z3i::Space, Shape > DigitalShape;
    typedef Z3i::KSpace KSpace;
    typedef typename KSpace::SCell SCell;
    typedef typename KSpace::Surfel Surfel;

    ASSERT (( noiseLevel < 1.0 ));
    bool withNoise = ( noiseLevel <= 0.0 ) ? false : true;
    if( withNoise )
        noiseLevel *= h;

    // Digitizer
    DigitalShape* dshape = new DigitalShape();
    dshape->attach( *aShape );
    dshape->init( border_min, border_max, h );

    KSpace K;
    if ( ! K.init( dshape->getLowerBound(), dshape->getUpperBound(), true ) )
    {
        std::cerr << "[3dLocalEstimators_0memory] error in creating KSpace." << std::endl;
        return false;
    }

    try
    {
        if ( withNoise )
        {
            typedef KanungoNoise< DigitalShape, Z3i::Domain > KanungoPredicate;
            typedef LightImplicitDigitalSurface< KSpace, KanungoPredicate > Boundary;
            typedef DigitalSurface< Boundary > MyDigitalSurface;
            typedef typename MyDigitalSurface::ConstIterator ConstIterator;

            typedef DepthFirstVisitor< MyDigitalSurface > Visitor;
            typedef GraphVisitorRange< Visitor > VisitorRange;
            typedef typename VisitorRange::ConstIterator VisitorConstIterator;

            typedef PointFunctorFromPointPredicateAndDomain< KanungoPredicate, Z3i::Domain, unsigned int > MyPointFunctor;
            typedef FunctorOnCells< MyPointFunctor, KSpace > MySpelFunctor;

            // Extracts shape boundary
            KanungoPredicate * noisifiedObject = new KanungoPredicate( *dshape, dshape->getDomain(), noiseLevel );
            SCell bel = Surfaces< KSpace >::findABel( K, *noisifiedObject, 10000 );
            Boundary * boundary = new Boundary( K, *noisifiedObject, SurfelAdjacency< KSpace::dimension >( true ), bel );
            MyDigitalSurface surf ( *boundary );

            double minsize = dshape->getUpperBound()[0] - dshape->getLowerBound()[0];
            unsigned int tries = 0;
            while( surf.size() < 2 * minsize || tries > 150 )
            {
                delete boundary;
                bel = Surfaces< KSpace >::findABel( K, *noisifiedObject, 10000 );
                boundary = new Boundary( K, *noisifiedObject, SurfelAdjacency< KSpace::dimension >( true ), bel );
                surf = MyDigitalSurface( *boundary );
                ++tries;
            }

            if( tries > 150 )
            {
                std::cerr << "Can't found a proper bel. So .... I ... just ... kill myself." << std::endl;
                return false;
            }

            VisitorRange * range;
            VisitorConstIterator ibegin;
            VisitorConstIterator iend;

            // Estimations

            // True Mean Curvature

            if( options.at( 0 ) != '0' )
            {
                char full_filename[360];
                sprintf( full_filename, "%s%s", filename.c_str(), "_True_mean.dat" );
                std::ofstream file( full_filename );
                file << "# h = " << h << std::endl;
                file << "# True Mean Curvature estimation" << std::endl;

                std::ostream_iterator< double > out_it_true_mean( file, "\n" );

                range = new VisitorRange( new Visitor( surf, *surf.begin() ));
                ibegin = range->begin();
                iend = range->end();
                estimateTrueMeanCurvatureQuantity( ibegin,
                                                   iend,
                                                   out_it_true_mean,
                                                   K,
                                                   h,
                                                   aShape );
                file.close();
                delete range;
            }

            // True Gaussian Curvature

            if( options.at( 1 ) != '0' )
            {
                char full_filename[360];
                sprintf( full_filename, "%s%s", filename.c_str(), "_True_gaussian.dat" );
                std::ofstream file( full_filename );
                file << "# h = " << h << std::endl;
                file << "# True Gaussian Curvature estimation" << std::endl;

                std::ostream_iterator< double > out_it_true_gaussian( file, "\n" );

                range = new VisitorRange( new Visitor( surf, *surf.begin() ));
                ibegin = range->begin();
                iend = range->end();
                estimateTrueGaussianCurvatureQuantity( ibegin,
                                                       iend,
                                                       out_it_true_gaussian,
                                                       K,
                                                       h,
                                                       aShape );
                file.close();
                delete range;
            }


            Clock c;
            double re_convolution_kernel = radius_kernel * std::pow( h, 1.0/3.0 ); // to obtains convergence results, re must follow the rule re=kh^(1/3)

            if( options.at( 2 ) != '0' || options.at( 3 ) != '0' )
            {
                MyPointFunctor * pointFunctor = new MyPointFunctor( noisifiedObject, dshape->getDomain(), 1, 0 );
                MySpelFunctor * functor = new MySpelFunctor( *pointFunctor, K );

                if( options.at( 2 ) != '0' )
                {
                    // Integral Invariant Mean Curvature
                    IntegralInvariantMeanCurvatureEstimator< KSpace, MySpelFunctor > * IIMeanCurvatureEstimator = new IntegralInvariantMeanCurvatureEstimator< KSpace, MySpelFunctor >( K, *functor );

                    c.startClock();
                    IIMeanCurvatureEstimator->init ( h, re_convolution_kernel );

                    char full_filename[360];
                    sprintf( full_filename, "%s%s", filename.c_str(), "_II_mean.dat" );
                    std::ofstream file( full_filename );
                    file << "# h = " << h << std::endl;
                    file << "# Mean Curvature estimation from the Integral Invariant" << std::endl;
                    file << "# computed kernel radius = " << re_convolution_kernel << std::endl;

                    range = new VisitorRange( new Visitor( surf, *surf.begin() ));
                    ibegin = range->begin();
                    iend = range->end();

                    std::ostream_iterator< double > out_it_ii_mean( file, "\n" );
                    if( !lambda_optimized )
                    {
                        IIMeanCurvatureEstimator->eval( ibegin, iend, out_it_ii_mean );
                    }
                    else
                    {
                        IIMeanCurvatureEstimator->eval( ibegin, iend, out_it_ii_mean, *aShape );
                    }
                    double TIIMeanCurv = c.stopClock();
                    file << "# time = " << TIIMeanCurv << std::endl;
                    file.close();

                    delete range;
                    delete IIMeanCurvatureEstimator;
                }

                if( options.at( 3 ) != '0' )
                {
                    // Integral Invariant Gaussian Curvature
                    IntegralInvariantGaussianCurvatureEstimator< KSpace, MySpelFunctor > * IIGaussianCurvatureEstimator = new IntegralInvariantGaussianCurvatureEstimator< KSpace, MySpelFunctor >( K, *functor );

                    c.startClock();
                    IIGaussianCurvatureEstimator->init( h, re_convolution_kernel );

                    char full_filename[360];
                    sprintf( full_filename, "%s%s", filename.c_str(), "_II_gaussian.dat" );
                    std::ofstream file( full_filename );
                    file << "# h = " << h << std::endl;
                    file << "# Gaussian Curvature estimation from the Integral Invariant" << std::endl;
                    file << "# computed kernel radius = " << re_convolution_kernel << std::endl;

                    std::ostream_iterator< double > out_it_ii_gaussian( file, "\n" );

                    range = new VisitorRange( new Visitor( surf, *surf.begin() ));
                    ibegin = range->begin();
                    iend = range->end();

                    if( !lambda_optimized )
                    {
                        IIGaussianCurvatureEstimator->eval ( ibegin, iend, out_it_ii_gaussian );
                    }
                    else
                    {
                        IIGaussianCurvatureEstimator->eval ( ibegin, iend, out_it_ii_gaussian, *aShape );
                    }
                    double TIIGaussCurv = c.stopClock();
                    file << "# time = " << TIIGaussCurv << std::endl;
                    file.close();

                    delete range;
                    delete IIGaussianCurvatureEstimator;
                }

                //delete noisifiedObject;
                delete boundary;
                delete pointFunctor;
                delete functor;
            }

            if( options.at( 5 ) != '0' )
            {
                // Monge Gaussian Curvature
                typedef MongeJetFittingGaussianCurvatureEstimator<Surfel, CanonicSCellEmbedder<KSpace> > FunctorGaussian;
                typedef LocalEstimatorFromSurfelFunctorAdapter<typename MyDigitalSurface::DigitalSurfaceContainer, Z3i::L2Metric, FunctorGaussian> ReporterK;
                CanonicSCellEmbedder<KSpace> * embedder = new CanonicSCellEmbedder<KSpace>( K );
                FunctorGaussian estimatorK( *embedder, h );
                ReporterK reporterK(surf.container(), Z3i::l2Metric, estimatorK);
                c.startClock();
                reporterK.init( h , re_convolution_kernel / h  );

                range = new VisitorRange( new Visitor( surf, *surf.begin() ));
                ibegin = range->begin();
                iend = range->end();

                char full_filename[360];
                sprintf( full_filename, "%s%s", filename.c_str(), "_MongeJetFitting_gaussian.dat" );
                std::ofstream file( full_filename );
                file << "# h = " << h << std::endl;
                file << "# Gaussian Curvature estimation from CGAL Monge from and Jet Fitting" << std::endl;
                file << "# computed kernel radius = " << re_convolution_kernel << std::endl;
                std::ostream_iterator< double > out_it_monge_gaussian( file, "\n" );
                reporterK.eval(ibegin, iend , out_it_monge_gaussian);
                double TMongeGaussCurv = c.stopClock();
                file << "# time = " << TMongeGaussCurv << std::endl;
                file.close();
                delete embedder;
            }

            if( options.at( 4 ) != '0' )
            {
                // Monge Mean Curvature
                typedef MongeJetFittingMeanCurvatureEstimator<Surfel, CanonicSCellEmbedder<KSpace> > FunctorMean;
                typedef LocalEstimatorFromSurfelFunctorAdapter<typename MyDigitalSurface::DigitalSurfaceContainer, Z3i::L2Metric, FunctorMean> ReporterH;
                CanonicSCellEmbedder<KSpace> * embedder = new CanonicSCellEmbedder<KSpace>( K );
                FunctorMean estimatorH( *embedder, h );
                ReporterH reporterH(surf.container(), Z3i::l2Metric, estimatorH);
                c.startClock();
                reporterH.init( h , re_convolution_kernel / h  );

                char full_filename[360];
                sprintf( full_filename, "%s%s", filename.c_str(), "_MongeJetFitting_mean.dat" );
                std::ofstream file( full_filename );
                file << "# h = " << h << std::endl;
                file << "# Mean Curvature estimation from CGAL Monge from and Jet Fitting" << std::endl;
                file << "# computed kernel radius = " << re_convolution_kernel << std::endl;
                std::ostream_iterator< double > out_it_monge_mean( file, "\n" );

                range = new VisitorRange( new Visitor( surf, *surf.begin() ));
                ibegin = range->begin();
                iend = range->end();
                reporterH.eval(ibegin, iend , out_it_monge_mean);
                double TMongeMeanCurv = c.stopClock();
                file << "# time = " << TMongeMeanCurv << std::endl;
                file.close();
                delete embedder;
            }
        }
        else
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
            SCell bel = Surfaces<KSpace>::findABel ( K, *dshape, 10000 );
            Boundary boundary( K, *dshape, SurfelAdjacency< KSpace::dimension >( true ), bel );
            MyDigitalSurface surf ( boundary );

            VisitorRange * range;
            VisitorConstIterator ibegin;
            VisitorConstIterator iend;

            // Estimations
            Clock c;
            // True Mean Curvature

            if( options.at( 0 ) != '0' )
            {
                char full_filename[360];
                sprintf( full_filename, "%s%s", filename.c_str(), "_True_mean.dat" );
                std::ofstream file( full_filename );
                file << "# h = " << h << std::endl;
                file << "# True Mean Curvature estimation" << std::endl;

                std::ostream_iterator< double > out_it_true_mean( file, "\n" );

                range = new VisitorRange( new Visitor( surf, *surf.begin() ));
                ibegin = range->begin();
                iend = range->end();

                c.startClock();

                estimateTrueMeanCurvatureQuantity( ibegin,
                                                   iend,
                                                   out_it_true_mean,
                                                   K,
                                                   h,
                                                   aShape );

                double TTrueMeanCurv = c.stopClock();
                file << "# time = " << TTrueMeanCurv << std::endl;

                file.close();
                delete range;
            }

            // True Gaussian Curvature

            if( options.at( 1 ) != '0' )
            {
                char full_filename[360];
                sprintf( full_filename, "%s%s", filename.c_str(), "_True_gaussian.dat" );
                std::ofstream file( full_filename );
                file << "# h = " << h << std::endl;
                file << "# True Gaussian Curvature estimation" << std::endl;

                std::ostream_iterator< double > out_it_true_gaussian( file, "\n" );

                range = new VisitorRange( new Visitor( surf, *surf.begin() ));
                ibegin = range->begin();
                iend = range->end();

                c.startClock();

                estimateTrueGaussianCurvatureQuantity( ibegin,
                                                       iend,
                                                       out_it_true_gaussian,
                                                       K,
                                                       h,
                                                       aShape );

                double TTrueGaussianCurv = c.stopClock();
                file << "# time = " << TTrueGaussianCurv << std::endl;

                file.close();

                delete range;
            }


            double re_convolution_kernel = radius_kernel * std::pow( h, 1.0/3.0 ); // to obtains convergence results, re must follow the rule re=kh^(1/3)

            if( options.at( 2 ) != '0' || options.at( 3 ) != '0' )
            {
                MyPointFunctor * pointFunctor = new MyPointFunctor( dshape, dshape->getDomain(), 1, 0 );
                MySpelFunctor * functor = new MySpelFunctor( *pointFunctor, K );

                if( options.at( 2 ) != '0' )
                {
                    // Integral Invariant Mean Curvature
                    IntegralInvariantMeanCurvatureEstimator< KSpace, MySpelFunctor > * IIMeanCurvatureEstimator = new IntegralInvariantMeanCurvatureEstimator< KSpace, MySpelFunctor >( K, *functor );

                    c.startClock();
                    IIMeanCurvatureEstimator->init ( h, re_convolution_kernel );

                    char full_filename[360];
                    sprintf( full_filename, "%s%s", filename.c_str(), "_II_mean.dat" );
                    std::ofstream file( full_filename );
                    file << "# h = " << h << std::endl;
                    file << "# Mean Curvature estimation from the Integral Invariant" << std::endl;
                    file << "# computed kernel radius = " << re_convolution_kernel << std::endl;

                    range = new VisitorRange( new Visitor( surf, *surf.begin() ));
                    ibegin = range->begin();
                    iend = range->end();

                    std::ostream_iterator< double > out_it_ii_mean( file, "\n" );
                    if( !lambda_optimized )
                    {
                        IIMeanCurvatureEstimator->eval( ibegin, iend, out_it_ii_mean );
                    }
                    else
                    {
                        IIMeanCurvatureEstimator->eval( ibegin, iend, out_it_ii_mean, *aShape );
                    }

                    double TIIMeanCurv = c.stopClock();
                    file << "# time = " << TIIMeanCurv << std::endl;
                    file.close();

                    delete range;
                    delete IIMeanCurvatureEstimator;
                }

                if( options.at( 3 ) != '0' )
                {
                    // Integral Invariant Gaussian Curvature
                    IntegralInvariantGaussianCurvatureEstimator< KSpace, MySpelFunctor > * IIGaussianCurvatureEstimator = new IntegralInvariantGaussianCurvatureEstimator< KSpace, MySpelFunctor >( K, *functor );

                    c.startClock();
                    IIGaussianCurvatureEstimator->init( h, re_convolution_kernel );

                    char full_filename[360];
                    sprintf( full_filename, "%s%s", filename.c_str(), "_II_gaussian.dat" );
                    std::ofstream file( full_filename );
                    file << "# h = " << h << std::endl;
                    file << "# Gaussian Curvature estimation from the Integral Invariant" << std::endl;
                    file << "# computed kernel radius = " << re_convolution_kernel << std::endl;

                    std::ostream_iterator< double > out_it_ii_gaussian( file, "\n" );

                    range = new VisitorRange( new Visitor( surf, *surf.begin() ));
                    ibegin = range->begin();
                    iend = range->end();

                    if( !lambda_optimized )
                    {
                        IIGaussianCurvatureEstimator->eval ( ibegin, iend, out_it_ii_gaussian );
                    }
                    else
                    {
                        IIGaussianCurvatureEstimator->eval ( ibegin, iend, out_it_ii_gaussian, *aShape );
                    }
                    double TIIGaussCurv = c.stopClock();
                    file << "# time = " << TIIGaussCurv << std::endl;
                    file.close();

                    delete range;
                    delete IIGaussianCurvatureEstimator;
                }

                delete pointFunctor;
                delete functor;
            }

            if( options.at( 5 ) != '0' )
            {
                // Monge Gaussian Curvature
                typedef MongeJetFittingGaussianCurvatureEstimator<Surfel, CanonicSCellEmbedder<KSpace> > FunctorGaussian;
                typedef LocalEstimatorFromSurfelFunctorAdapter<typename MyDigitalSurface::DigitalSurfaceContainer, Z3i::L2Metric, FunctorGaussian> ReporterK;
                CanonicSCellEmbedder<KSpace> * embedder = new CanonicSCellEmbedder<KSpace>( K );
                FunctorGaussian estimatorK( *embedder, h );
                ReporterK reporterK(surf.container(), Z3i::l2Metric, estimatorK);
                c.startClock();
                reporterK.init( h , re_convolution_kernel / h  );

                range = new VisitorRange( new Visitor( surf, *surf.begin() ));
                ibegin = range->begin();
                iend = range->end();

                char full_filename[360];
                sprintf( full_filename, "%s%s", filename.c_str(), "_MongeJetFitting_gaussian.dat" );
                std::ofstream file( full_filename );
                file << "# h = " << h << std::endl;
                file << "# Gaussian Curvature estimation from CGAL Monge from and Jet Fitting" << std::endl;
                file << "# computed kernel radius = " << re_convolution_kernel << std::endl;
                std::ostream_iterator< double > out_it_monge_gaussian( file, "\n" );
                reporterK.eval(ibegin, iend , out_it_monge_gaussian);
                double TMongeGaussCurv = c.stopClock();
                file << "# time = " << TMongeGaussCurv << std::endl;
                file.close();
                delete embedder;
            }

            if( options.at( 4 ) != '0' )
            {
                // Monge Mean Curvature
                typedef MongeJetFittingMeanCurvatureEstimator<Surfel, CanonicSCellEmbedder<KSpace> > FunctorMean;
                typedef LocalEstimatorFromSurfelFunctorAdapter<typename MyDigitalSurface::DigitalSurfaceContainer, Z3i::L2Metric, FunctorMean> ReporterH;
                CanonicSCellEmbedder<KSpace> * embedder = new CanonicSCellEmbedder<KSpace>( K );
                FunctorMean estimatorH( *embedder, h );
                ReporterH reporterH(surf.container(), Z3i::l2Metric, estimatorH);
                c.startClock();
                reporterH.init( h , re_convolution_kernel / h  );

                char full_filename[360];
                sprintf( full_filename, "%s%s", filename.c_str(), "_MongeJetFitting_mean.dat" );
                std::ofstream file( full_filename );
                file << "# h = " << h << std::endl;
                file << "# Mean Curvature estimation from CGAL Monge from and Jet Fitting" << std::endl;
                file << "# computed kernel radius = " << re_convolution_kernel << std::endl;
                std::ostream_iterator< double > out_it_monge_mean( file, "\n" );

                range = new VisitorRange( new Visitor( surf, *surf.begin() ));
                ibegin = range->begin();
                iend = range->end();
                reporterH.eval(ibegin, iend , out_it_monge_mean);
                double TMongeMeanCurv = c.stopClock();
                file << "# time = " << TMongeMeanCurv << std::endl;
                file.close();
                delete embedder;
            }
        }
    }
    catch ( InputException e )
    {
        std::cerr << "[estimatorCurvatureComparator3D]"
                  << " error."
                  << e.what()
                  << " "
                  << filename << " "
                  << border_min[0] << " "
                  << border_max[0] << " "
                  << h << " "
                  << radius_kernel << " "
                  << lambda_optimized << " "
                  << options << " "
                  << alpha << " "
                  << std::endl;
        return false;
    }

    delete dshape;

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
            ("output,o", po::value< std::string >(), "Output file")
            ("radius,r",  po::value< double >(), "Kernel radius for IntegralInvariant" )
            ("alpha",  po::value<double>()->default_value(1.0/3.0), "Alpha parameter for Integral Invariant computation" )
            ("h",  po::value< double >(), "Grid step" )
            ("minAABB,a",  po::value< double >()->default_value( -10.0 ), "Min value of the AABB bounding box (domain)" )
            ("maxAABB,A",  po::value< double >()->default_value( 10.0 ), "Max value of the AABB bounding box (domain)" )
            ("noise,n",  po::value<double>()->default_value(0.0), "Level of noise to perturb the shape" )
            ("lambda,l",  po::value< bool >()->default_value( false ), "Use the shape to get a better approximation of the surface (optional)" )
            ("estimators,e",  po::value< std::string >()->default_value("111100"), "the i-th estimator is disabled iff there is a 0 at position i" );


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
                    << "\t3dlocalEstimators_0memory --shape <shape> --h <h> --radius <radius> --estimators <binaryWord> --output <output>"<<std::endl
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
                    << std::endl;
        return 0;
    }

    if (!(vm.count("output"))) missingParam("--output");
    if (!(vm.count("shape"))) missingParam("--shape");
    if (!(vm.count("h"))) missingParam("--h");
    if (!(vm.count("radius"))) missingParam("--radius");


    std::string file_export = vm["output"].as< std::string >();
    std::string options = vm["estimators"].as< std::string >();
    double h = vm["h"].as< double >();
    double radius = vm["radius"].as< double >();
    double alpha = vm["alpha"].as< double >();
    std::string poly_str = vm["shape"].as< std::string >();
    bool lambda_optimized = vm["lambda"].as< bool >();
    double noiseLevel = vm["noise"].as<double>();

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

    compareShapeEstimators< Z3i::Space, ImplicitShape > (
                file_export,
                shape,
                border_min, border_max,
                h,
                radius,
                alpha,
                options,
                lambda_optimized,
                noiseLevel );

    delete shape;
}
