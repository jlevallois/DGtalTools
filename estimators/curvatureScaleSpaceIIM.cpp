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
 * @file curvatureScaleSpace.cpp
 * @ingroup estimators
 *
 * @author Bertrand Kerautret (\c kerautre@loria.fr )
 * LORIA (CNRS, UMR 7503), University of Nancy, France
 *
 *
 * @date 2014/04/02
 *
 * Output the image of curvature scale space
 * 
 * This file is part of the DGtal library.
 */

///////////////////////////////////////////////////////////////////////////////
#include <iostream>

#include "DGtal/base/Common.h"

#include "DGtal/shapes/ShapeFactory.h"
#include "DGtal/shapes/Shapes.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/topology/helpers/Surfaces.h"

//Grid curve
#include "DGtal/geometry/curves/FreemanChain.h"
#include "DGtal/geometry/curves/GridCurve.h"

//Estimators
#include "DGtal/geometry/surfaces/estimation/IntegralInvariantMeanCurvatureEstimator.h"
#include "DGtal/geometry/surfaces/FunctorOnCells.h"

// Generation of resulting image
#include "DGtal/images/ImageContainerBySTLVector.h"
#include "DGtal/io/writers/GenericWriter.h"
#include "DGtal/io/colormaps/GradientColorMap.h"
#include "DGtal/io/writers/PPMWriter.h"
#include "DGtal/io/colormaps/HueShadeColorMap.h"

// images
#include "DGtal/images/ImageHelper.h"
#include "DGtal/images/imagesSetsUtils/SimpleThresholdForegroundPredicate.h"

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

#include "DGtal/io/boards/Board2D.h"
#include "DGtal/io/colormaps/GradientColorMap.h"
#include "DGtal/topology/LightImplicitDigitalSurface.h"
#include "DGtal/topology/DigitalSurface.h"
#include "DGtal/graph/DepthFirstVisitor.h"
#include "DGtal/graph/GraphVisitorRange.h"

#include <vector>
#include <string>
#include <iomanip>


 using namespace DGtal;


 DGtal::ImageContainerBySTLVector<Z2i::Domain, unsigned char> 
 getImageFromFC(double scaleMax, const FreemanChain<Z2i::Integer> &fc, 
   std::vector<Z2i::KSpace::SCell> &boundarySCell)
 {
  int minx, miny, maxx, maxy;
  fc.computeBoundingBox(minx, miny, maxx, maxy);
  minx-=scaleMax;
  miny-=scaleMax;
  maxx+=scaleMax;
  maxy+=scaleMax;
  Z2i::KSpace aKSpace;
  aKSpace.init(Z2i::Point(minx, miny), Z2i::Point(maxx, maxy), false);
  
  std::set<Z2i::KSpace::SCell> bSCell;
  FreemanChain<Z2i::Integer>::getContourSCell(aKSpace, fc, bSCell, true); 
  std::set<Z2i::KSpace::Cell> interiorCell;
  Surfaces<Z2i::KSpace>::uComputeInterior(aKSpace, bSCell, interiorCell, false);  
  ImageContainerBySTLVector<Z2i::Domain, unsigned char> imageResult (Z2i::Domain(Z2i::Point(minx, miny), Z2i::Point(maxx, maxy))); 
  for(std::set<Z2i::KSpace::Cell>::const_iterator it = interiorCell.begin(); it!=interiorCell.end(); ++it){
    imageResult.setValue(aKSpace.uCoords(*it), 1);
  }

  ///////////////
  boundarySCell.clear();

  typedef Z2i::KSpace::Surfel Surfel;
  typedef PointPredicateFromPointFunctor< ImageContainerBySTLVector<Z2i::Domain, unsigned char>, Z2i::Domain, unsigned char > PP;
  typedef LightImplicitDigitalSurface< Z2i::KSpace, PP > LightImplicitDigSurface;
  typedef DigitalSurface< LightImplicitDigSurface > MyDigitalSurface;
  Z2i::KSpace KSpaceShape;
  KSpaceShape.init( imageResult.domain().lowerBound(), imageResult.domain().upperBound(), true );
  SurfelAdjacency<Z2i::KSpace::dimension> SAdj( true );
  PP pp( imageResult, 1 );
  Surfel bel = Surfaces<Z2i::KSpace>::findABel( KSpaceShape, pp, 100000 );
  LightImplicitDigSurface LightImplDigSurf( KSpaceShape, pp, SAdj, bel );
  MyDigitalSurface digSurf( LightImplDigSurf );

  typedef DepthFirstVisitor< MyDigitalSurface > Visitor;
  typedef GraphVisitorRange< Visitor > VisitorRange;
  typedef VisitorRange::ConstIterator SurfelConstIterator;

  VisitorRange range( new Visitor( digSurf, *digSurf.begin() ) );
  SurfelConstIterator abegin = range.begin();
  SurfelConstIterator aend = range.end();

  while (abegin != aend )
  {
    boundarySCell.push_back( *abegin);
    ++abegin;
  }
  ///////////////
  return imageResult; 
}


void
computeCurvatureIIMC(double h, double scale, ImageContainerBySTLVector<Z2i::Domain, unsigned char>  &image,
 const std::vector<Z2i::KSpace::SCell> &boundarySCell, std::vector<double> &resCurvature)
{
  typedef DGtal::KhalimskySpaceND< 2, int > KSpace;
  typedef DGtal::ImageContainerBySTLVector<Z2i::Domain, unsigned char> Image2D ; 
  typedef KhalimskySpaceND<2, int>::Cell Cell;
  typedef KhalimskySpaceND<2, int>::SCell SCell;
  /*typedef SimpleThresholdForegroundPredicate< Image2D > ImageThreshPredicate;
  typedef ImageToConstantFunctor< Image2D, ImageThreshPredicate > MyPointFunctor;*/
  typedef FunctorOnCells< Image2D, KSpace > MyCellFunctor;

  KSpace aKSpace;
  aKSpace.init(image.domain().lowerBound(), image.domain().upperBound(), false);
  
  
  /*ImageThreshPredicate predicate = ImageThreshPredicate( image, 0 );
  MyPointFunctor pointFunctor( image, predicate, 1 ); */   
  MyCellFunctor functor ( image, aKSpace );  
  typedef IntegralInvariantMeanCurvatureEstimator< Z2i::KSpace, MyCellFunctor > MyIIMeanEstimator;
  MyIIMeanEstimator estimator ( aKSpace, functor );

  estimator.init(h , scale ); // Initialisation for a given Euclidean radius of the convolution kernel
  std::back_insert_iterator< std::vector< double > > resultsIterator( resCurvature );
  estimator.eval ( boundarySCell.begin(), boundarySCell.end(), resultsIterator ); // Computation

}

struct MyDrawStyleCustomColor : public DrawableWithBoard2D
{
  Color myPenColor;
  Color myFillColor;
  MyDrawStyleCustomColor( const Color & penColor,
        const Color & fillColor )
    : myPenColor( penColor ), myFillColor( fillColor )
  {}

  virtual void setStyle( Board2D & aboard) const
  {
    aboard.setFillColor( myFillColor);
    aboard.setPenColor( myPenColor );
  }
};


///////////////////////////////////////////////////////////////////////////////
namespace po = boost::program_options;

int main( int argc, char** argv )
{

  // parse command line ----------------------------------------------
  po::options_description general_opt("Allowed options are: ");
  general_opt.add_options()
  ("help,h", "display this message")
  ("FreemanChain,f", po::value<std::string>(), "FreemanChain file name")
  ("gridStep,s", po::value<double>()->default_value(1.0), "Grid step size")
  ("radiusInit", po::value<double>()->default_value(1.0), "initial radius")
  ("radiusIncrement", po::value<double>()->default_value(1.0), "radius increment")
  ("radiusFinal", po::value<double>()->default_value(1.0), "final radius") 
  ("output,o ", po::value<std::string>(), "set the output name ")
  ("curvatureCutOff,c", po::value<double>()->default_value(10.0), "set the curvature limits to better display");
  

  typedef DGtal::ImageContainerBySTLVector<Z2i::Domain, unsigned char> Image2D ; 
  typedef ImageContainerBySTLVector<Z2i::Domain, double > ImageCurvature;


  bool parseOK=true;
  po::variables_map vm;
  try{
    po::store(po::parse_command_line(argc, argv, general_opt), vm);  
  }catch(const std::exception& ex){
    parseOK=false;
    trace.info()<< "Error checking program options: "<< ex.what()<< std::endl;
  }
  po::notify(vm);    
  if(!parseOK || vm.count("help")||argc<=1 || (!(vm.count("FreemanChain"))) )
  {
    trace.info()<< "Generate the Curvature Scale Sapce image using a binomial convolver based estimator." <<std::endl
    << "The x axis is associated to the contour point and the y axis to the scale. The color represent the curvature values included between the cutoff values (set to 10 by default)."
    <<std::endl << "Basic usage: "<<std::endl
    << "\t curvatureScaleSpaceIMM -f ${DGtal}/examples/samples/contourS.fc  --radiusInit 2 --radiusIncrement  0.5 --radiusFinal 100 -o cssResu.ppm "<<std::endl
    << general_opt << "\n";
    return 0;
  }
  double radius_initial = vm["radiusInit"].as<double>();
  double radius_increment = vm["radiusIncrement"].as<double>();
  double radius_final = vm["radiusFinal"].as<double>();
  double curvatureCutOff = vm["curvatureCutOff"].as<double>();
  double gridStep = vm["gridStep"].as<double>();
  
  if(vm.count("FreemanChain")){
    std::string fileName = vm["FreemanChain"].as<std::string>();
    
    std::vector< DGtal::FreemanChain<Z2i::Integer>  > vectFcs =  PointListReader< Z2i::Point >:: getFreemanChainsFromFile<Z2i::Integer> (fileName);     
    bool isClosed = vectFcs.at(0).isClosed(); 
    
    // Preparing resulting image:
    unsigned int height =  (int)((radius_final-radius_initial)/radius_increment);
    // We add one point since the 
    Z2i::Domain domain (Z2i::Point(0,0), Z2i::Point(vectFcs.at(0).size()+(isClosed? 0: 1), height-1));
    ImageCurvature cssImage(domain);    
    HueShadeColorMap<double>  gradCurvature (-curvatureCutOff, curvatureCutOff);
    
    trace.progressBar(0, height);
    double rad= radius_initial;
    std::vector<Z2i::KSpace::SCell> boundarySCell;
    Image2D image = getImageFromFC( radius_final/gridStep, vectFcs.at(0), boundarySCell );

////////////////////////////
    typename std::vector<Z2i::KSpace::SCell>::const_iterator dit = boundarySCell.begin();
    typename std::vector<Z2i::KSpace::SCell>::const_iterator ditEnd = boundarySCell.end();
    int nb = 0;
    for(; dit != ditEnd; ++dit)
    {
      ++nb;
    }

    typedef GradientColorMap< int > Gradient;
    Gradient cmap_grad( 0, nb );
    cmap_grad.addColor( Color( 50, 50, 255 ) );
    cmap_grad.addColor( Color( 255, 0, 0 ) );
    cmap_grad.addColor( Color( 255, 255, 10 ) );
    /*Color red( 255, 0, 0 );
    Color dred( 192, 0, 0 );*/

    dit = boundarySCell.begin();
    Board2D board;
    board << SetMode( (*dit).className(), "Paving" );
    std::string specificStyle = (*dit).className() + "/Paving";
    for(int i = 0; dit != ditEnd; ++dit)
    {
      board << CustomStyle( specificStyle, new CustomColors( Color::Black, cmap_grad( i )))//CustomStyle( (*dit).className(), new MyDrawStyleCustomColor( red, dred ) )
            << *dit;
      ++i;
    }
    board.saveEPS( "export.eps");
//////////////////////////////

  for(double l= 0; l < height; l++ ){

    trace.progressBar(l, height);
    std::vector<double> curvaturesIIM;

    computeCurvatureIIMC(gridStep, rad, image, boundarySCell,  curvaturesIIM);      
      // Output
    unsigned int j = 0;
    for ( std::vector<double>::const_iterator it = curvaturesIIM.begin(), it_end = curvaturesIIM.end();
      it != it_end; ++it, ++j ) {
      double c = (*it);//*gridStep;
    c = c<-curvatureCutOff? -curvatureCutOff: c;
    c = c>curvatureCutOff? curvatureCutOff: c;
    cssImage.setValue(Z2i::Point(j, l), c); 
  }      
  rad=rad+radius_increment;
}

trace.progressBar(height, height);
trace.info() <<std::endl;

DGtal::GenericWriter<ImageCurvature, 2, double, HueShadeColorMap<double> >::exportFile(vm["output"].as<std::string>(), cssImage, gradCurvature );       
}

return 0;
}

