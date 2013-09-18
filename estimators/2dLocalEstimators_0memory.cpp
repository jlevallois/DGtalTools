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
#include "DGtal/geometry/surfaces/estimation/IntegralInvariantMeanCurvatureEstimator_0memory.h"

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
 * Estimation. Merely call the init and eval methods of the
 * given estimator.
 *
 * @param estimator any local estimator
 * @param h the grid step
 * @param itb begin iterator
 * @param ite end iterator
 * @param ito output iterator on estimated quantities
 */
template < typename Shape, typename KSpace, typename ConstIterator, typename OutputIterator >
void
estimationTrueCurvature( const ConstIterator & it_begin,
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
    Point currentInnerRealPoint;

    for ( ConstIterator it = it_begin; it != it_end; ++it )
    {
        currentInnerSpel = K.sDirectIncident( *it_begin, K.sOrthDir( *it_begin ) ); /// Spel on the border, but inside the shape
        currentInnerRealPoint = ((RealPoint)(K.sCoords( currentInnerSpel ))) * h;
        output = aShape->curvature( aShape->parameter( currentInnerRealPoint ));
        ++output;
    }
}

/**
 * Estimation of tangents and curvature
 * from several different methods
 *
 * @param aShape shape
 * @param h grid step
 * @param noiseLevel level to noised the shape. 0 <= noiseLevel < 1
 * @param radiusKernel optional - Euclidean radius of the convolution kernel for Integral Invariants estimators. 0.0 by default (aka no noise).
 * @param alphaII optional - Specification of alpha param for the radius of the convolution kernel for Integral Invariants estimators. 1/3 by default.
 *
 */
template <typename Space, typename Shape>
bool
computeLocalEstimations( const std::string & filename,
                         Shape * aShape,
                         const double & h,
                         const double & radius_kernel,
                         const double & alpha,
                         const bool & lambda_optimized )
{
    // Types
    typedef typename Space::Point Point;
    typedef typename Space::Vector Vector;
    typedef typename Space::RealPoint RealPoint;
    typedef typename Space::Integer Integer;
    typedef Z2i::KSpace KSpace;
    typedef typename KSpace::SCell SCell;
    typedef GaussDigitizer< Space, Shape > DigitalShape;
    typedef LightImplicitDigitalSurface< KSpace, DigitalShape > Boundary;
    typedef DigitalSurface< Boundary > MyDigitalSurface;

    typedef PointFunctorFromPointPredicateAndDomain< DigitalShape, Z2i::Domain, unsigned int > MyPointFunctor;
    typedef FunctorOnCells< MyPointFunctor, KSpace > MySpelFunctor;

    // Digitizer
    DigitalShape* dshape = new DigitalShape();
    dshape->attach( *aShape );
    Vector vlow(-1, -1); Vector vup(1, 1);
    dshape->init( aShape->getLowerBound() + vlow, aShape->getUpperBound() + vup, h );

    KSpace K;
    if ( ! K.init( dshape->getLowerBound(), dshape->getUpperBound(), true ) )
    {
        std::cerr << "[2dLocalEstimators_0memory] error in creating KSpace." << std::endl;
        return false;
    }
    try
    {

        // Extracts shape boundary
        SCell bel = Surfaces< KSpace >::findABel( K, *dshape, 10000 );
        Boundary boundary( K, *dshape, SurfelAdjacency< KSpace::dimension >( true ), bel );
        MyDigitalSurface surf ( boundary );

        // Estimations

        // True
        char full_filename[360];
        sprintf( full_filename, "%s%s", filename.c_str(), "_True_curvature.dat" );
        std::ofstream file( full_filename );
        file.flags( std::ios_base::unitbuf );
        file << "# h = " << h << std::endl;
        file << "# True Curvature estimation" << std::endl;

        std::ostream_iterator< double > out_it_true( file, "\n" );

        estimationTrueCurvature( surf.begin(),
                                 surf.end(),
                                 out_it_true,
                                 K,
                                 h,
                                 aShape );

        file.close();

        //Integral Invariants

        Clock c;
        double re_convolution_kernel = radius_kernel * std::pow( h, alpha );

        MyPointFunctor * pointFunctor = new MyPointFunctor( dshape, dshape->getDomain(), 1, 0 );
        MySpelFunctor * functor = new MySpelFunctor( *pointFunctor, K );

        IntegralInvariantMeanCurvatureEstimator_0memory< KSpace, MySpelFunctor > * IIMeanCurvatureEstimator = new IntegralInvariantMeanCurvatureEstimator_0memory< KSpace, MySpelFunctor >( K, *functor );

        c.startClock();
        IIMeanCurvatureEstimator->init ( h, re_convolution_kernel );


        memset(&full_filename[0], 0, sizeof(full_filename));
        sprintf( full_filename, "%s%s", filename.c_str(), "_II_curvature.dat" );
        file.open( full_filename );
        file.flags( std::ios_base::unitbuf );
        file << "# h = " << h << std::endl;
        file << "# Curvature estimation from the Integral Invariant" << std::endl;
        file << "# computed kernel radius = " << re_convolution_kernel << std::endl;

        std::ostream_iterator< double > out_it_ii( file, "\n" );
        if( !lambda_optimized )
        {
            IIMeanCurvatureEstimator->eval( surf.begin(), surf.end(), out_it_ii );
        }
        else
        {
            IIMeanCurvatureEstimator->eval( surf.begin(), surf.end(), out_it_ii, *aShape );
        }

        double TIICurv = c.stopClock();
        file << "# time = " << TIICurv << std::endl;
        file.close();

        delete pointFunctor;
        delete functor;
        delete IIMeanCurvatureEstimator;

        return true;
    }
    catch ( InputException e )
    {
        std::cerr << "[computeLocalEstimations]"
                  << " error in finding a bel." << std::endl;
        return false;
    }
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
            ("output,o", po::value< std::string >(), "Output")
            ("shape,s", po::value< std::string >(), "Shape name")
            ("radius,R",  po::value<double>(), "Radius of the shape" )
            ("kernelradius,K",  po::value<double>()->default_value(0.0), "Radius of the convolution kernel (Integral invariants estimators)" )
            ("alpha",  po::value<double>()->default_value(1.0/3.0), "Alpha parameter for Integral Invariant computation" )
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
            ("lambda,l",  po::value< bool >()->default_value( false ), "Use the shape to get a better approximation of the surface (optional)" );


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
                    << "\tlocalEstimators --shape <shapeName> [required parameters] --kernelradius <radius> --gridstep <h> --output <output>"<<std::endl
                    << std::endl;
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
    if (!(vm.count("output"))) missingParam("--output");

    std::string shapeName = vm["shape"].as< std::string >();
    std::string filename = vm["output"].as< std::string >();

    //We check that the shape is known
    unsigned int id = checkAndReturnIndex(shapeName);

    // standard types
    typedef Z2i::Space Space;
    typedef Space::Point Point;
    typedef Space::RealPoint RealPoint;

    RealPoint center( vm["center_x"].as< double >(),
            vm["center_y"].as< double >() );
    double h = vm["gridstep"].as< double >();

    double radius = vm["kernelradius"].as< double >();
    double alpha = vm["alpha"].as< double >();
    bool lambda = vm["lambda"].as< bool >();

    if (id ==0)
    {
        if (!(vm.count("radius"))) missingParam("--radius");
        //if (!(vm.count("kernelradius"))) missingParam("--kernelradius");
        double radius = vm["radius"].as<double>();

        Ball2D<Space> * ball = new Ball2D<Space>( center, radius );
        computeLocalEstimations<Space>( filename, ball, h, radius, alpha, lambda );
        delete ball;
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

        Flower2D<Space> * flower = new Flower2D<Space>( center, radius, varsmallradius, k, phi );
        computeLocalEstimations<Space>( filename, flower, h, radius, alpha, lambda );
        delete flower;
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

        NGon2D<Space> * object = new NGon2D<Space>( center, radius, k, phi );
        computeLocalEstimations<Space>( filename, object, h, radius, alpha, lambda );
        delete object;
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

        AccFlower2D<Space> * accflower = new AccFlower2D<Space>( center, radius, varsmallradius, k, phi );
        computeLocalEstimations<Space>( filename, accflower, h, radius, alpha, lambda );
        delete accflower;
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

        Ellipse2D<Space> * ellipse = new Ellipse2D<Space>( center, a1, a2, phi );
        computeLocalEstimations<Space>( filename, ellipse, h, radius, alpha, lambda );
        delete ellipse;
    }
}
