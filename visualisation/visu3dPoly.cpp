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
#include "DGtal/topology/SurfelAdjacency.h"
#include "DGtal/topology/LightImplicitDigitalSurface.h"
#include "DGtal/topology/DigitalSurface.h"
#include "DGtal/graph/DepthFirstVisitor.h"
#include "DGtal/graph/GraphVisitorRange.h"
#include "DGtal/shapes/GaussDigitizer.h"
#include "DGtal/shapes/implicit/ImplicitPolynomial3Shape.h"
#include "DGtal/math/MPolynomial.h"
#include "DGtal/io/readers/MPolynomialReader.h"

// Drawing
#include "DGtal/io/viewers/Viewer3D.h"
#include "DGtal/io/colormaps/GradientColorMap.h"
#include <QtGui/QApplication>

using namespace std;
using namespace DGtal;


///////////////////////////////////////////////////////////////////////////////

void usage( int argc, char** argv )
{
    trace.info() << "Usage: " << argv[ 0 ]
          << " <polynomial form> <step>"<< std::endl;
    trace.info() << "\t - <polynomial form> - try \"0.03x^4 + 0.03y^4 + 0.03z^4 - 2x^2 - 2y^2 - 2z^2 + 8\" for example."<< std::endl;
    trace.info() << "\t - <step> grid step."<< std::endl;
    trace.info() << "Example : "<< argv[ 0 ] << " \"0.03x^4 + 0.03y^4 + 0.03z^4 - 2x^2 - 2y^2 - 2z^2 + 8\" 0.2 "<< std::endl;
}

int main( int argc, char** argv )
{
    if( argc < 3 )
    {
        usage( argc, argv );
        return -1;
    }

    std::string poly_str = argv[1];//"1000x^2y^2z^2 + 3x^2 + 3y^2 + z^2 - 1";
    double h = atof( argv[2] );

    trace.beginBlock("Euclidean construction");

    typedef Z3i::Space::RealPoint RealPoint;
    typedef Z3i::Space::RealPoint::Coordinate Ring;
    typedef MPolynomial< 3, Ring > Polynomial3;
    typedef MPolynomialReader<3, Ring> Polynomial3Reader;
    typedef ImplicitPolynomial3Shape<Z3i::Space> Shape;
    typedef Z3i::Space Space;

    RealPoint border_min( -10, -10, -10 );
    RealPoint border_max( 10, 10, 10 );

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

    Shape* aShape = new Shape( poly );

    trace.endBlock();

    trace.beginBlock("Digitization");

    typedef typename Space::RealPoint RealPoint;
    typedef GaussDigitizer< Z3i::Space, Shape > DigitalShape;
    typedef Z3i::KSpace KSpace;
    typedef typename KSpace::SCell SCell;
    typedef typename KSpace::Cell Cell;
    typedef typename KSpace::Surfel Surfel;

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

    typedef LightImplicitDigitalSurface< KSpace, DigitalShape > Boundary;
    typedef DigitalSurface< Boundary > MyDigitalSurface;
    typedef typename MyDigitalSurface::ConstIterator ConstIterator;

    typedef DepthFirstVisitor< MyDigitalSurface > Visitor;
    typedef GraphVisitorRange< Visitor > VisitorRange;
    typedef typename VisitorRange::ConstIterator VisitorConstIterator;

    // Extracts shape boundary
    SCell bel = Surfaces< KSpace >::findABel( K, *dshape, 10000 );
    Boundary * boundary = new Boundary( K, *dshape, SurfelAdjacency< KSpace::dimension >( true ), bel );
    MyDigitalSurface surf ( *boundary );

    trace.endBlock();

    VisitorRange * range;
    VisitorConstIterator ibegin;
    VisitorConstIterator iend;

    range = new VisitorRange( new Visitor( surf, *surf.begin() ));
    ibegin = range->begin();
    iend = range->end();

    trace.beginBlock("viewer");

    QApplication application( argc, argv );
    Viewer3D viewer;

    viewer.show();
    //viewer << SetMode3D( dshape->getDomain().className(), "BoundingBox" ) << dshape->getDomain();

    for( ; ibegin != iend; ++ibegin )
    {
        // viewer << CustomColors3D( Color::Black, cmap_grad( b ));
        Cell unsignedSurfel = K.uCell( K.sKCoords(*ibegin) );
        viewer << unsignedSurfel;
    }

    delete range;

    viewer << Viewer3D::updateDisplay;

    trace.endBlock();
    return application.exec();
}

///////////////////////////////////////////////////////////////////////////////
