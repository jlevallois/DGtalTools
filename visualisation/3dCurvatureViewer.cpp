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

// Shape constructors
#include "DGtal/io/readers/VolReader.h"
#include "DGtal/images/ImageSelector.h"
#include "DGtal/images/imagesSetsUtils/SetFromImage.h"
#include "DGtal/topology/SurfelAdjacency.h"
#include "DGtal/topology/helpers/Surfaces.h"
#include "DGtal/topology/LightImplicitDigitalSurface.h"
#include "DGtal/topology/DigitalSurface.h"
#include "DGtal/topology/DepthFirstVisitor.h"

// Integral Invariant includes
#include "DGtal/geometry/surfaces/FunctorOnCells.h"
#include "DGtal/geometry/surfaces/estimation/IntegralInvariantMeanCurvatureEstimator.h"
#include "DGtal/geometry/surfaces/estimation/IntegralInvariantGaussianCurvatureEstimator.h"

// Drawing
#include "DGtal/io/viewers/Viewer3D.h"
#include "DGtal/io/colormaps/GradientColorMap.h"
#include <QtGui/QApplication>

using namespace std;
using namespace DGtal;


///////////////////////////////////////////////////////////////////////////////

void usage( int argc, char** argv )
{
    trace.error() << "Usage: " << argv[ 0 ]
                           << " <fileName.vol> <re> <\"mean\" || \"gaussian\" || \"princurv\">";
    trace.error() << "\t - <filename.vol> file you want to show the curvature information.";
    trace.error() << "\t - <re> Euclidean radius of the kernel.";
    trace.error() << "\t - <\"mean\" || \"gaussian\"> show mean or Gaussian curvature on shape.";
    trace.error() << "Example : "<< argv[ 0 ] << " Al.150.vol 7.39247665 \"mean\"";
}

int main( int argc, char** argv )
{
    if ( argc != 4 )
    {
        usage( argc, argv );
        return 0;
    }

    double h = 1.0;
    double re_convolution_kernel = atof(argv[2]);

    std::string mode = argv[ 3 ];
    if ( mode != "gaussian" && mode != "mean" && mode != "princurv" )
    {
        usage( argc, argv );
        return 0;
    }

    // Construction of the shape from vol file
    typedef ImageSelector< Z3i::Domain, bool>::Type Image;
    typedef Z3i::KSpace::Surfel Surfel;
    typedef LightImplicitDigitalSurface< Z3i::KSpace, Image > MyLightImplicitDigitalSurface;
    typedef DigitalSurface< MyLightImplicitDigitalSurface > MyDigitalSurface;

    std::string filename = argv[1];
    Image image = VolReader<Image>::importVol( filename );

    Z3i::Domain domain = image.domain();

    Z3i::KSpace KSpaceShape;

    bool space_ok = KSpaceShape.init( domain.lowerBound(), domain.upperBound(), true );
    if (!space_ok)
    {
      trace.error() << "Error in the Khalimsky space construction."<<std::endl;
      return 2;
    }

    SurfelAdjacency< Z3i::KSpace::dimension > SAdj( true );
    Surfel bel = Surfaces< Z3i::KSpace >::findABel( KSpaceShape, image, 100000 );
    MyLightImplicitDigitalSurface LightImplDigSurf( KSpaceShape, image, SAdj, bel );
    MyDigitalSurface digSurf( LightImplDigSurf );

    typedef DepthFirstVisitor<MyDigitalSurface> Visitor;
    typedef Visitor::VertexConstIterator SurfelConstIterator;
    Visitor *depth = new Visitor (digSurf, *digSurf.begin());
    SurfelConstIterator abegin = SurfelConstIterator(depth);
    SurfelConstIterator aend = SurfelConstIterator(0);


    // Integral Invariant stuff
    typedef FunctorOnCells< Image, Z3i::KSpace > MyFunctor;


    MyFunctor functor ( image, KSpaceShape, domain ); // Creation of a functor on Cells, returning true if the cell is inside the shape

    QApplication application( argc, argv );
    Viewer3D viewer;
    viewer.show();
//    viewer << SetMode3D(image.domain().className(), "BoundingBox") << image.domain();

    Visitor *depth2 = new Visitor (digSurf, *digSurf.begin());
    SurfelConstIterator abegin2 = SurfelConstIterator(depth2);

    if( mode == "mean" || mode == "gaussian" )
    {
        typedef double Quantity;
        std::vector< Quantity > results;
        back_insert_iterator< std::vector< Quantity > > resultsIterator( results ); // output iterator for results of Integral Invariante curvature computation

        if ( mode == "mean" )
        {
            typedef IntegralInvariantMeanCurvatureEstimator< Z3i::KSpace, MyFunctor > MyIIMeanEstimator;

            MyIIMeanEstimator estimator ( KSpaceShape, functor );
            estimator.init( h, re_convolution_kernel ); // Initialisation for a given Euclidean radius of the convolution kernel
            estimator.eval ( abegin, aend, resultsIterator ); // Computation
        }
        else if ( mode == "gaussian" )
        {
            typedef IntegralInvariantGaussianCurvatureEstimator< Z3i::KSpace, MyFunctor > MyIIGaussianEstimator;

            MyIIGaussianEstimator estimator ( KSpaceShape, functor );
            estimator.init( h, re_convolution_kernel ); // Initialisation for a given Euclidean radius of the convolution kernel
            estimator.eval ( abegin, aend, resultsIterator ); // Computation
        }

        // Drawing results
        Quantity min = numeric_limits < Quantity >::max();
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
        cmap_grad.addColor( Color( 255, 255, 10 ) );

        for ( unsigned int i = 0; i < results.size(); ++i )
        {
            viewer << CustomColors3D( Color::Black, cmap_grad( results[ i ] ))
                   << *abegin2;
            ++abegin2;
        }
    }
    else
    {
        typedef double Quantity;
        typedef SimpleMatrix< double, 3, 3 > Matrix3x3;
        typedef EigenValues3D::Vector3 Vector3;
        typedef CurvatureInformation< Quantity, Matrix3x3, Vector3 > CurvInformation;

        std::vector< CurvInformation > results;
        back_insert_iterator< std::vector< CurvInformation > > resultsIterator( results ); // output iterator for results of Integral Invariante curvature computation

        typedef IntegralInvariantGaussianCurvatureEstimator< Z3i::KSpace, MyFunctor > MyIIGaussianEstimator;

        MyIIGaussianEstimator estimator ( KSpaceShape, functor );
        estimator.init( h, re_convolution_kernel ); // Initialisation for a given Euclidean radius of the convolution kernel
        estimator.evalComplete ( abegin, aend, resultsIterator ); // Computation

        // Drawing results
        SCellToMidPoint< Z3i::KSpace > midpoint( KSpaceShape );
        typedef typename Matrix3x3::RowVector RowVector;
        for ( unsigned int i = 0; i < results.size(); ++i )
        {
            CurvInformation current = results[ i ];
            Z3i::Point center = midpoint( *abegin2 );

            viewer << CustomColors3D( DGtal::Color(255,255,255,255),
                                      DGtal::Color(255,255,255,255))
                   << *abegin2;

//            RowVector normal = current.eigenVectors.row(0).getNormalized(); // don't show the normal
            RowVector curv1 = current.eigenVectors.row(1).getNormalized();
            RowVector curv2 = current.eigenVectors.row(2).getNormalized();


//            viewer.addLine ( center[0],
//                             center[1],
//                             center[2],
//                             center[0] +  normal[0],
//                             center[1] +  normal[1],
//                             center[2] +  normal[2],
//                             DGtal::Color ( 20,200,200 ), 5.0 ); // don't show the normal


            viewer.addLine ( center[0] -  0.5 * curv1[0],
                             center[1] -  0.5 * curv1[1],
                             center[2] -  0.5 * curv1[2],
                             center[0] +  0.5 * curv1[0],
                             center[1] +  0.5 * curv1[1],
                             center[2] +  0.5 * curv1[2],
                             DGtal::Color ( 20,200,200 ), 5.0 );

            viewer.addLine ( center[0] -  0.5 * curv2[0],
                             center[1] -  0.5 * curv2[1],
                             center[2] -  0.5 * curv2[2],
                             center[0] +  0.5 * curv2[0],
                             center[1] +  0.5 * curv2[1],
                             center[2] +  0.5 * curv2[2],
                             DGtal::Color ( 200,20,20 ), 5.0 );

            ++abegin2;
        }
    }

    viewer << Viewer3D::updateDisplay;
    return application.exec();
}

///////////////////////////////////////////////////////////////////////////////
