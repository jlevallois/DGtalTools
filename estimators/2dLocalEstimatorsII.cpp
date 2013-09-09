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
 * @file 2dLocalEstimators.cpp
 * @ingroup Tools
 * @author Jacques-Olivier Lachaud (\c jacques-olivier.lachaud@univ-savoie.fr )
 * Laboratory of Mathematics (CNRS, UMR 5807), University of Savoie, France
 * @author David Coeurjolly (\c david.coeurjolly@liris.cnrs.fr)
 * LIRIS (CNRS, UMR 5205),
 * @author Tristan Roussillon (\c tristan.roussillon@liris.cnrs.fr )
 * Laboratoire d'InfoRmatique en Image et Systèmes d'information - LIRIS (CNRS, UMR 5205), CNRS,
 * France
 * @author Jeremy Levallois (\c jeremy.levallois@liris.cnrs.fr )
 * Laboratoire d'InfoRmatique en Image et Systèmes d'information - LIRIS (CNRS, UMR 5205), Universite de Lyon, France
 * LAboratoire de MAthematiques - LAMA (CNRS, UMR 5807), Universite de Savoie, France
 *
 * @date 2011/07/04
 *
 * DGtal shape generator
 *
 * This file is part of the DGtal library.
 */

///////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <iomanip>
#include <vector>
#include <string>

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

#include "DGtal/base/Common.h"
#include "DGtal/base/Clock.h"

//shapes
#include "DGtal/shapes/ShapeFactory.h"
#include "DGtal/shapes/Shapes.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/topology/helpers/Surfaces.h"
#include "DGtal/topology/DigitalSurface.h"

//Digitizer
#include "DGtal/shapes/GaussDigitizer.h"
#include "DGtal/geometry/curves/GridCurve.h"
#include "DGtal/topology/LightImplicitDigitalSurface.h"
#include "DGtal/graph/DepthFirstVisitor.h"
#include "DGtal/graph/GraphVisitorRange.h"
#include "DGtal/geometry/volumes/KanungoNoise.h"


//Estimators
#include "DGtal/geometry/curves/estimation/TrueLocalEstimatorOnPoints.h"
#include "DGtal/geometry/curves/estimation/TrueGlobalEstimatorOnPoints.h"
#include "DGtal/geometry/curves/estimation/ParametricShapeCurvatureFunctor.h"
#include "DGtal/geometry/curves/estimation/ParametricShapeTangentFunctor.h"
#include "DGtal/geometry/curves/estimation/ParametricShapeArcLengthFunctor.h"

#include "DGtal/geometry/curves/BinomialConvolver.h"
#include "DGtal/geometry/curves/estimation/MostCenteredMaximalSegmentEstimator.h"
#include "DGtal/geometry/curves/estimation/SegmentComputerEstimators.h"
#include "DGtal/geometry/curves/ArithmeticalDSS.h"
#include "DGtal/geometry/curves/GeometricalDCA.h"

#include "DGtal/images/ImageHelper.h"
#include "DGtal/geometry/surfaces/FunctorOnCells.h"
#include "DGtal/geometry/surfaces/estimation/IntegralInvariantMeanCurvatureEstimator.h"

#include "DGtal/kernel/BasicPointFunctors.h"

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

template< typename RealPoint >
struct OptionsIntegralInvariant
{
    std::vector< double > alpha; // <! Alpha parameter for the convolution kernel. 1/3 by default
    double radius; // <! Radius of the convolution kernel.
    RealPoint center; // <! Center of the shape.
};


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

/** 
 * Missing parameter error message.
 *
 * @param param
 */
void missingParam(std::string param)
{
    trace.error() <<" Parameter: "<<param<<" is required.";
    trace.info()<<std::endl;
    exit(1);
}

/** 
 * Estimation error message.
 *
 * @param currentSize number of values returned by the estimator
 * @param expectedSize expected number of values
 */
void estimationError(int currentSize, int expectedSize)
{
    if (currentSize != expectedSize)
    {
        trace.error() << " error in the estimation"
                      << " (got " << currentSize << " values"
                      << " instead of " << expectedSize << ")";
        trace.info() << std::endl;
        exit(1);
    }

}

/** 
 * Estimation. Merely call the init and eval methods of the
 * given estimator.
 *
 * @param estimator any local estimator
 * @param h the grid step
 * @param itb begin iterator
 * @param ite end iterator
 * @param ito output iterator on estimated quantities
 */
template <typename Estimator, typename ConstIterator, typename OutputIterator>
void
estimation( Estimator & estimator, double h,  
            const ConstIterator& itb, const ConstIterator& ite, const OutputIterator& ito )
{
    Clock c;
    c.startClock();
    estimator.init( h, itb, ite );
    estimator.eval( itb, ite, ito );
    double time = c.stopClock();
    std::cout << "# Time: " << time << std::endl;
}


/**
 *
 * @return Euclidean radius for the convolver of Integral Invariant estimators
 */
template< typename ConstIteratorOnPoints, typename Point >
unsigned int suggestedSizeIntegralInvariant( const double h,
                                             const Point& center,
                                             const ConstIteratorOnPoints& itb,
                                             const ConstIteratorOnPoints& ite )
{
    typedef typename Point::Component TValue;

    ConstIteratorOnPoints it = itb;
    Point p( *it );
    Point distance = p - center;
    TValue minRadius = distance.norm();
    ++it;

    for ( ; it != ite; ++it )
    {
        p = *it;
        distance = p - center;
        if ( distance.norm() < minRadius )
        {
            minRadius = distance.norm();
        }
    }

    return minRadius * h;
}


/**
 * Estimation of tangents and curvature
 * from several different methods
 *
 * @param name shape name
 * @param aShape shape
 * @param h grid step
 * @param noiseLevel level to noised the shape. 0 <= noiseLevel < 1
 * @param radiusKernel optional - Euclidean radius of the convolution kernel for Integral Invariants estimators. 0.0 by default (aka no noise).
 * @param alphaII optional - Specification of alpha param for the radius of the convolution kernel for Integral Invariants estimators. 1/3 by default.
 *
 */
template< typename Space, typename Shape >
bool
computeLocalEstimations( Shape & aShape,
                         double h,
                         struct OptionsIntegralInvariant< Z2i::RealPoint > optionsII )
{
    /// Shape creation

    typedef Z2i::Domain Domain;
    typedef Z2i::KSpace KSpace;
    typedef typename KSpace::SCell SCell;
    typedef typename Space::Point Point;
    typedef typename Space::Vector Vector;
    typedef typename Space::RealPoint RealPoint;
    typedef typename Space::Integer Integer;
    typedef GaussDigitizer<Space,Shape> Digitizer;

    Digitizer* dig = new Digitizer();
    dig->attach( aShape ); // attaches the shape.
    Vector vlow(-1,-1); Vector vup(1,1);
    dig->init( aShape.getLowerBound()+vlow, aShape.getUpperBound()+vup, h );
    Domain domain = dig->getDomain();

    KSpace K;
    bool ok = K.init( dig->getLowerBound(), dig->getUpperBound(), true );
    if ( ! ok )
    {
        std::cerr << "[2dLocalEstimatorsII]"
                  << " error in creating KSpace." << std::endl;
        return false;
    }

    SurfelAdjacency< KSpace::dimension > SAdj( true );
    SCell bel;
    try
    {
        bel = Surfaces< KSpace >::findABel( K, *dig, 10000 );
    }
    catch( ... )
    {
        exit( -42 );
    }

    typedef LightImplicitDigitalSurface< KSpace, Digitizer > LightImplicitDigSurface;
    typedef DigitalSurface< LightImplicitDigSurface > DigSurface;
    typedef DepthFirstVisitor< DigSurface > Visitor;
    typedef GraphVisitorRange< Visitor > VisitorRange;
    typedef typename VisitorRange::ConstIterator ConstIterator;
    typedef PointFunctorFromPointPredicateAndDomain< Digitizer, Domain, unsigned int > MyPointFunctor;
    typedef FunctorOnCells< MyPointFunctor, KSpace > CurvatureIIFct;
    typedef IntegralInvariantMeanCurvatureEstimator< KSpace, CurvatureIIFct > IntegralInvariantCurvatureEstimator;
    typedef typename IntegralInvariantCurvatureEstimator::Quantity Quantity;

    LightImplicitDigSurface LightImplDigSurf( K, *dig, SAdj, bel );
    DigSurface digSurf( LightImplDigSurf );
    MyPointFunctor pointFunctor( dig, domain, 1, 0 );
    CurvatureIIFct functor ( pointFunctor, K );

    std::cout << "# h = " << h << std::endl;
    std::cout << "# domain size = " << domain.size() << std::endl;

    std::vector< double > trueCurvatures;
    std::vector< std::vector< Quantity > > estimatedCurvatures;

    /// True values
    {
        std::cout << "# ---------- True curvature ----------" << std::endl;

        std::vector< Point > points;
        Surfaces< KSpace >::track2DBoundaryPoints( points, K, SAdj, *dig, bel );
        GridCurve< KSpace > gridcurve;
        gridcurve.initFromVector( points );

        typedef typename GridCurve< KSpace >::PointsRange PointsRange;
        PointsRange pointsRange = gridcurve.getPointsRange();

        typedef ParametricShapeCurvatureFunctor< Shape > CurvatureFunctor;
        typedef typename PointsRange::ConstCirculator C;
        TrueLocalEstimatorOnPoints< C, Shape, CurvatureFunctor > trueCurvatureEstimator;

        Clock c;
        c.startClock();

        trueCurvatureEstimator.attach( &aShape );
        trueCurvatureEstimator.init( h, pointsRange.c(), pointsRange.c() );
        trueCurvatureEstimator.eval( pointsRange.c(), pointsRange.c(), std::back_inserter( trueCurvatures ));

        double time = c.stopClock();
        std::cout << "# Time (s): " << time / 1000 << std::endl;

        estimationError( trueCurvatures.size(), pointsRange.size() );

        if( optionsII.radius <= 0.0 )
        {
            VisitorRange range( new Visitor( digSurf, *digSurf.begin() ) );
            optionsII.radius = suggestedSizeIntegralInvariant( h, dig->round( optionsII.center ), pointsRange.begin(), pointsRange.end() );
        }
    }

    /// Integral Inv.

    for( int i = 0; i < optionsII.alpha.size(); ++i )
    {
        estimatedCurvatures.push_back( std::vector< Quantity >() );

        double alpha = optionsII.alpha[ i ];
        double re_convolution_kernel = optionsII.radius * std::pow( h, alpha );

        std::cout << "# ---------- Curvature estimation from integral invariants with alpha=" << alpha  << "----------"<< std::endl;
        std::cout << "# Alpha: " << alpha << std::endl;
        std::cout << "# Euclidean radius: " << optionsII.radius << std::endl;
        std::cout << "# Radius: " << re_convolution_kernel << std::endl;

        VisitorRange range( new Visitor( digSurf, *digSurf.begin() ) );
        IntegralInvariantCurvatureEstimator IICurvatureEstimator( K, functor );

        Clock c;
        c.startClock();

        IICurvatureEstimator.init( h, re_convolution_kernel );

        std::back_insert_iterator< std::vector< Quantity > > bckInserterIterator( estimatedCurvatures[ i ] );
        IICurvatureEstimator.eval( range.begin(), range.end(), bckInserterIterator );

        double time = c.stopClock();
        std::cout << "# Time (s): " << time / 1000 << std::endl;

        estimationError( estimatedCurvatures[ i ].size(), trueCurvatures.size() );
    }

    /// Output

    std::cout << "# id k";
    for( int i = 0; i < optionsII.alpha.size(); ++i )
    {
        std::cout << " alpha" << optionsII.alpha[ i ] << "k";
    }
    std::cout << std::endl;

    unsigned int prsize = trueCurvatures.size() - 1;
    for( int i = 0; i <= prsize; ++i )
    {
        std::cout << i << " "
                  << trueCurvatures[ prsize - i ] << " ";
        for( int j = 0; j < optionsII.alpha.size(); ++j )
        {
            std::cout << estimatedCurvatures[ j ][ i ] << " ";
        }

        std::cout << std::endl;
    }

    delete dig;

    return true;
}


///////////////////////////////////////////////////////////////////////////////
namespace po = boost::program_options;

int main( int argc, char** argv )
{
    // parse command line ----------------------------------------------
    po::options_description general_opt("Allowed options are");
    general_opt.add_options()
            ("help,h", "display this message")
            ("list,l",  "List all available shapes")
            ("shape,s", po::value<std::string>(), "Shape name")
            ("radius,R",  po::value<double>(), "Radius of the shape" )
            ("kernelradius,K",  po::value<double>()->default_value(0.0), "Radius of the convolution kernel (Integral invariants estimators)" )
            ("alpha", po::value<std::vector<double> >()->multitoken(), "List of Alpha parameter for Integral Invariant computation.")
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
            ("gridstep,g",  po::value<double>()->default_value(1.0), "Grid step for the digitization" )
            ("noise,n",  po::value<double>()->default_value(0.0), "Level of noise to perturb the shape" )
            ("properties",  po::value<std::string>()->default_value("11"), "the i-th property is disabled iff there is a 0 at position i" )
            ("estimators,e",  po::value<std::string>()->default_value("1000"), "the i-th estimator is disabled iff there is a 0 at position i" );
    

    bool parseOK=true;
    po::variables_map vm;
    try{
        po::store(po::parse_command_line(argc, argv, general_opt), vm);
    }catch(const std::exception& ex){
        parseOK=false;
        trace.info()<< "Error checking program options: "<< ex.what()<< std::endl;
    }
    po::notify(vm);
    if(!parseOK || vm.count("help")||argc<=1)
    {
        trace.info()<< "Compare local estimators on implicit shapes using DGtal library" <<std::endl
                    << "Basic usage: "<<std::endl
                    << "\tlocalEstimators --shape <shapeName> [required parameters] --estimators <binaryWord> --properties <binaryWord>"<<std::endl
                    << std::endl
                    << "Below are the different available families of estimators: " << std::endl
                    << "\t - Maximal DSS based estimators" << std::endl
                    << "\t - Maximal DCA based estimators" << std::endl
                    << "\t - Binomial convolver based estimators" << std::endl
                    << "\t - Integral Invariants based estimators" << std::endl
                    << std::endl
                    << "The i-th family of estimators is enabled if the i-th character of the binary word is not 0. "
                    << "The default binary word is '100'. This means that the first family of estimators, "
                    << "ie. maximal DSS based estimators, is enabled, whereas the next ones are disabled. "
                    << std::endl
                    << "Below are the different available properties: " << std::endl
                    << "\t - Tangeant" << std::endl
                    << "\t - Curvature" << std::endl
                    << std::endl
                    << general_opt << std::endl;
        return 0;
    }

    //List creation
    createList();

    if (vm.count("list"))
    {
        displayList();
        return 0;
    }

    //Parse options
    if (!(vm.count("shape"))) missingParam("--shape");
    std::string shapeName = vm["shape"].as<std::string>();
    
    int nb = 4; //number of available methods
    std::string options = vm["estimators"].as<std::string>();
    if (options.size() < nb)
    {
        trace.error() << " At least " << nb
                      << " characters are required "
                      << " with option --estimators.";
        trace.info() << std::endl;
        exit(1);
    }

    nb = 2; //number of available properties
    std::string properties = vm["properties"].as<std::string>();
    if (properties.size() < nb)
    {
        trace.error() << " At least " << nb
                      << " characters are required "
                      << " with option --properties.";
        trace.info() << std::endl;
        exit(1);
    }
    
    //We check that the shape is known
    unsigned int id = checkAndReturnIndex(shapeName);

    // standard types
    typedef Z2i::Space Space;
    typedef Space::Point Point;
    typedef Space::RealPoint RealPoint;

    RealPoint center( vm["center_x"].as<double>(),
            vm["center_y"].as<double>() );
    double h = vm["gridstep"].as<double>();

    struct OptionsIntegralInvariant< RealPoint > optII;
    optII.radius = vm["kernelradius"].as<double>();
    optII.alpha = vm["alpha"].as<std::vector<double> >();
    optII.center = center;

    double noiseLevel = vm["noise"].as<double>();

    if (id ==0)
    {
        if (!(vm.count("radius"))) missingParam("--radius");
        //if (!(vm.count("kernelradius"))) missingParam("--kernelradius");
        double radius = vm["radius"].as<double>();

        Ball2D<Space> ball(Z2i::Point(0,0), radius);
        computeLocalEstimations<Space>( ball, h, optII );
    }
    else if (id ==1)
    {
        if (!(vm.count("width"))) missingParam("--width");
        double width = vm["width"].as<double>();

        ImplicitHyperCube<Space> object(Z2i::Point(0,0), width/2);
        trace.error()<< "Not available.";
        trace.info()<<std::endl;
    }
    else if (id ==2)
    {
        if (!(vm.count("power"))) missingParam("--power");
        if (!(vm.count("radius"))) missingParam("--radius");
        double radius = vm["radius"].as<double>();
        double power = vm["power"].as<double>();

        ImplicitRoundedHyperCube<Space> ball( Z2i::Point(0,0), radius, power );
        trace.error()<< "Not available.";
        trace.info()<<std::endl;
    }
    else if (id ==3)
    {
        if (!(vm.count("varsmallradius"))) missingParam("--varsmallradius");
        if (!(vm.count("radius"))) missingParam("--radius");
        if (!(vm.count("k"))) missingParam("--k");
        if (!(vm.count("phi"))) missingParam("--phi");
        //if (!(vm.count("kernelradius"))) missingParam("--kernelradius");
        double radius = vm["radius"].as<double>();
        double varsmallradius = vm["varsmallradius"].as<double>();
        unsigned int k = vm["k"].as<unsigned int>();
        double phi = vm["phi"].as<double>();

        Flower2D<Space> flower( center, radius, varsmallradius, k, phi );
        computeLocalEstimations<Space>( flower, h, optII );
    }
    else if (id ==4)
    {
        if (!(vm.count("radius"))) missingParam("--radius");
        if (!(vm.count("k"))) missingParam("--k");
        if (!(vm.count("phi"))) missingParam("--phi");
        //if (!(vm.count("kernelradius"))) missingParam("--kernelradius");
        double radius = vm["radius"].as<double>();
        unsigned int k = vm["k"].as<unsigned int>();
        double phi = vm["phi"].as<double>();

        NGon2D<Space> object( center, radius, k, phi );
        computeLocalEstimations<Space>( object, h, optII );
    }
    else if (id ==5)
    {
        if (!(vm.count("varsmallradius"))) missingParam("--varsmallradius");
        if (!(vm.count("radius"))) missingParam("--radius");
        if (!(vm.count("k"))) missingParam("--k");
        if (!(vm.count("phi"))) missingParam("--phi");
        //if (!(vm.count("kernelradius"))) missingParam("--kernelradius");
        double radius = vm["radius"].as<double>();
        double varsmallradius = vm["varsmallradius"].as<double>();
        unsigned int k = vm["k"].as<unsigned int>();
        double phi = vm["phi"].as<double>();

        AccFlower2D<Space> accflower( center, radius, varsmallradius, k, phi );
        computeLocalEstimations<Space>( accflower, h, optII );
    }
    else if (id ==6)
    {
        if (!(vm.count("axis1"))) missingParam("--axis1");
        if (!(vm.count("axis2"))) missingParam("--axis2");
        if (!(vm.count("phi"))) missingParam("--phi");
        //if (!(vm.count("kernelradius"))) missingParam("--kernelradius");
        double a1 = vm["axis1"].as<double>();
        double a2 = vm["axis2"].as<double>();
        double phi = vm["phi"].as<double>();

        Ellipse2D<Space> ellipse( center, a1, a2, phi );
        computeLocalEstimations<Space>( ellipse, h, optII );
    }
}
