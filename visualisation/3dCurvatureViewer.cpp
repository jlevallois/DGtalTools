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
#include <cstring>

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

#include "DGtal/math/Statistic.h"
#include "DGtal/geometry/curves/SaturatedSegmentation.h"
#include "DGtal/geometry/curves/ArithmeticalDSSComputer.h"
#include "DGtal/topology/DigitalSurface2DSlice.h"
#include "DGtal/topology/DigitalSurface.h"

// Shape constructors
#include "DGtal/io/readers/VolReader.h"
#include "DGtal/io/readers/RawReader.h"
#include "DGtal/images/ImageSelector.h"
#include "DGtal/images/imagesSetsUtils/SetFromImage.h"
#include "DGtal/images/imagesSetsUtils/SimpleThresholdForegroundPredicate.h"
#include "DGtal/topology/SurfelAdjacency.h"
#include "DGtal/topology/helpers/Surfaces.h"
#include "DGtal/topology/LightImplicitDigitalSurface.h"
#include "DGtal/images/ImageHelper.h"
#include "DGtal/topology/DigitalSurface.h"
#include "DGtal/graph/DepthFirstVisitor.h"
#include "DGtal/graph/GraphVisitorRange.h"

// Integral Invariant includes
#include "DGtal/geometry/surfaces/FunctorOnCells.h"
#include "DGtal/geometry/surfaces/estimation/IntegralInvariantMeanCurvatureEstimator.h"
#include "DGtal/geometry/surfaces/estimation/IntegralInvariantGaussianCurvatureEstimator.h"

// Drawing
#include "DGtal/io/boards/Board3D.h"
#include "DGtal/io/viewers/Viewer3D.h"
#include "DGtal/io/colormaps/GradientColorMap.h"
#include <QtGui/QApplication>

using namespace std;
using namespace DGtal;

const Color  AXIS_COLOR_RED( 200, 20, 20, 255 );
const Color  AXIS_COLOR_GREEN( 20, 200, 20, 255 );
const Color  AXIS_COLOR_BLUE( 20, 20, 200, 255 );
const Color  AXIS_COLOR_BLACK( 200, 200, 200, 255 );
const double AXIS_LINESIZE = 0.15;


///////////////////////////////////////////////////////////////////////////////

/**
 * Missing parameter error message.
 *
 * @param param
 */
void missingParam( std::string param )
{
  trace.error() << " Parameter: " << param << " is required.";
  trace.info() << std::endl;
}

namespace po = boost::program_options;

enum FileExtension
{
  VOL,
  RAW
};

template< typename Image >
Image FileReader( std::string filename, FileExtension extension )
{
  switch( extension )
  {
  case FileExtension::VOL:
    return VolReader<Image>::importVol( filename );
    break;

  case FileExtension::RAW:
    return RawReader< Image >::importRaw8( filename, Z3i::Domain::Vector( 512, 512, 512) );
    break;
  }

}

template< typename KSpace, typename Iterator >
void analyseAllLengthMS( std::vector< Statistic<double> > & statE,
                         Iterator itb,
                         Iterator ite )
{
  typedef typename KSpace::Space Space;
  typedef typename Space::Point Point;
  typedef typename Space::Vector Vector;
  typedef ArithmeticalDSSComputer< Iterator, int, 4 > SegmentComputer;
  typedef SaturatedSegmentation< SegmentComputer > Decomposition;
  typedef typename Decomposition::SegmentComputerIterator SegmentComputerIterator;
  typedef std::vector< SegmentComputer > VectorOfSegmentComputer;
  typedef std::map< Point, VectorOfSegmentComputer > Pmap;

  // Computes the tangential cover
  SegmentComputer algo;
  Iterator itbegin = itb;
  Iterator itend = ite;
  Decomposition theDecomposition( itbegin, itend, algo);

  Pmap map;
  //for( itbegin = itb; itbegin != itend; ++itbegin )
  do
  {
    map.insert( std::pair< Point, VectorOfSegmentComputer >( *itbegin, VectorOfSegmentComputer() ) );
    ++itbegin;
  } while( itbegin != itend );


  for ( SegmentComputerIterator scIt = theDecomposition.begin(), scItEnd = theDecomposition.end();
        scIt != scItEnd; ++scIt )
  {
    const SegmentComputer & sc = *scIt;
    for ( Iterator ptIt = sc.begin(), ptItEnd = sc.end(); ptIt != ptItEnd; ++ptIt )
    {
      typename Pmap::iterator mloc = map.find( *ptIt );
      if( mloc != map.end() )
      {
        mloc->second.push_back( sc );
      }
      else
      {
        trace.error() << "not found ?" << std::endl;
      }
    }
  }

  itbegin = itb;
  itend = ite;
  Dimension ii = 0;
  //for( itbegin = itb; itbegin != itend; ++itbegin )
  do
  {
    //statD[ii].clear();
    statE[ii].clear();
    typename Pmap::iterator mloc = map.find( *itbegin );
    ASSERT(( mloc != map.end() ));

    /////////////
    for( typename VectorOfSegmentComputer::iterator scIt = mloc->second.begin(), scItEnd = mloc->second.end(); scIt != scItEnd; ++scIt )
    {
      const SegmentComputer & sc = *scIt;
      /*int64_t l = 0;
          for ( Iterator ptIt = sc.begin(), ptItEnd = sc.end(); ptIt != ptItEnd; ++ptIt )
            ++l;
          statD[ii].addValue( (double) l );*/
      double v = (sc.back( ) - sc.front()).norm1() ; // sc.size();  //*( sc.end() - 1 ) - *( sc.begin() );
      statE[ii].addValue( v )  ; //v.norm( ) );
      //          std::cout << " v=" << v.norm() << std::endl;
    }

    /////////////

    ++ii;
    ++itbegin;
  } while( itbegin != itend );
}


template< typename Surfel >
Dimension findSurfel( const std::vector< Surfel > & surfels,
                      const Surfel & s )
{
  const Dimension surfels_size = surfels.size();
  bool found = false;
  Dimension position = surfels_size;

  for( Dimension ii = 0; !found && ii < surfels_size; ++ii )
  {
    if( surfels[ii] == s )
    {
      found = true;
      position = ii;
    }
  }

  return position;
}



template <typename K3D, typename K2D>
class SCellToMyPoint
{
public:
  typedef K3D KSpace3D;
  typedef K2D KSpace2D;
  typedef typename Z2i::Point Output;
  typedef Output Value;
  typedef typename KSpace3D::SCell Input;
  typedef typename KSpace3D::Point Point;

private:
  /**
       * Aliasing pointer on the Khalimsky space.
      */
  const KSpace3D* myK;
  Dimension i;

public:

  /**
       * Default constructor.
      */
  SCellToMyPoint() : myK(NULL), i(0) { }
  /**
       * Constructor.
       * @param aK a Khalimsky space
      */
  SCellToMyPoint(const KSpace3D& aK, Dimension myDimension) : myK(&aK), i(myDimension) { }

  /**
     * Copy constructor.
     * @param other any SCellToPoint functor
     */
  SCellToMyPoint(const SCellToMyPoint& other)
    : myK(other.myK), i(other.i) { }

  /**
     * Assignment.
     *
     * @param other the object to copy.
     * @return a reference on 'this'.
     */
  SCellToMyPoint& operator= ( const SCellToMyPoint & other )
  {
    if (this != &other)
    {
      myK = other.myK;
      i = other.i;
    }
    return *this;
  }

  /**
     * Returns a point (with integer coordinates)
     * from a scell (with khalimsky coordinates)
     * @param aSCell a scell
     * @return the corresponding point.
     */
  Output operator()(const Input& aSCell) const
  {
    Input s = aSCell;
    Dimension k = myK->sOrthDir( s );
    Dimension j = (k+1)%3;
    if ( j == i ) j = (i+1)%3;
    Input next_linel = myK->sDirectIncident( s, j );
    Input base_pointel = myK->sIncident( next_linel, i, false );
    Point p = myK->sCoords( base_pointel );
    Output q( p[(i+1)%3], p[(i+2)%3] );
    return q;
  }

}; // end of class SCellToMyPoint


template< typename ImplicitDigitalSurface, typename Surfel >
void computeSegments( std::vector< Surfel > & surfels,
                      std::map< Surfel, std::pair< Statistic< double >, Statistic< double > > > & segments,
                      const Z3i::KSpace & K,
                      const ImplicitDigitalSurface & impDigitalSurface )
{
  typedef typename ImplicitDigitalSurface::Tracker Tracker;
  typedef DigitalSurface2DSlice< Tracker > MySlice;
  typedef std::pair< Statistic< double >, Statistic< double > > PairOfStatistics;
  typedef std::map< Surfel, PairOfStatistics > SurfelMap;
  typedef std::map< Surfel, bool > MarqueMap;

  const Dimension surfels_size = surfels.size();

  for( Dimension dim = 0; dim < 3; ++dim )
  {
    MarqueMap marque;
    for( Dimension ii = 0; ii < surfels_size; ++ii )
    {
      marque[ surfels[ ii ] ] = false;
    }

    for( Dimension ii = 0; ii < surfels_size; ++ii )
    {
      Surfel currentSurfel = surfels[ ii ];
      if( marque[ currentSurfel ] )
      {
        continue;
      }
      else if( K.sOrthDir( currentSurfel ) == dim )
      {
        continue;
      }
      else
      {
        Tracker ptrTracker( impDigitalSurface, currentSurfel );
        Dimension otherdim = 0;
        {
          bool dims[3] = { false, false, false };
          dims[dim] = true;
          dims[K.sOrthDir( currentSurfel )] = true;

          while( dims[otherdim] == true )
          {
            ++otherdim;
          }
        }

        //          unsigned int bitAll = 7;
        /*unsigned int bitSuppr = std::pow(2, dim) + std::pow(2, K.sOrthDir( currentSurfel ));
        unsigned int bitResu = 7 ^ bitSuppr;

        bitResu >> (sizeof(unsigned int) * CHAR_BIT - 1);

        unsigned int aaa = 0;
        unsigned int bbb = 1;
        unsigned int ccc = 2;
        unsigned int ddd = 4;

        aaa = aaa >> (sizeof(unsigned int) * CHAR_BIT - 1);
        bbb = bbb >> (sizeof(unsigned int) * CHAR_BIT - 1);
        ccc = ccc >> (sizeof(unsigned int) * CHAR_BIT - 1);
        ddd = ddd >> (sizeof(unsigned int) * CHAR_BIT - 1);

        std::cout << aaa << " " << bbb << " " << ccc << " " << ddd << std::endl;*/

        //          std::cout << "dim: " << dim << " --- orth: " << K.sOrthDir( currentSurfel ) << " --- result " << otherdim << std::endl;

        //dimSlice
        MySlice slice( &ptrTracker, otherdim );


        typedef typename MySlice::ConstIterator ConstIterator3D;
        typedef typename MySlice::ConstCirculator ConstCirculator3D;
        typedef typename MySlice::Iterator Iterator3D;
//        typedef SCellProjector< KhalimskySpaceND<2,int> > Functor;
//        typedef CanonicSCellEmbedder< Functor::KSpace > Functor2;
        typedef SCellToMyPoint< typename MySlice::KSpace, Z2i::KSpace > Functor2;
//        typedef Myfunctor<typename Functor::KSpace::Point , typename Functor::KSpace::SCell> Functor2;


        /*typedef ConstIteratorAdapter< ConstCirculator3D, Functor, Functor::SCell > ConstIterator2D;
        typedef ConstIteratorAdapter< ConstIterator3D, Functor, Functor::SCell > ConstIterator2D2;
        typedef ConstIteratorAdapter< ConstIterator2D, Functor2, Functor2::Output > ConstIterator2DP;
        typedef ConstIteratorAdapter< ConstIterator2D2, Functor2, Functor2::Output > ConstIterator2DP2;*/

        typedef ConstIteratorAdapter< ConstIterator3D, Functor2 > ConstIterator3D2P;

        ConstIterator3D a = slice.begin();
        ConstIterator3D b = ++(slice.begin());
        Dimension dimm = 0;
        while( a->myCoordinates[dimm] != b->myCoordinates[dimm] )
        {
          ++dimm;
        }

        /*Functor projector;
        projector.initRemoveOneDim( dimm );
        ConstIterator2D xbegin( slice.c(), projector );
        ConstIterator2D xend( slice.c(), projector );

        ConstIterator2D2 xxbegin( slice.begin(), projector );
        ConstIterator2D2 xxend( slice.end(), projector );

        Functor::KSpace k2d;*/
//        Functor2 pointFunctor;
        Functor2 pointFunctor(K, dim);

        ConstIterator3D2P ppbegin( slice.begin(), pointFunctor );
        ConstIterator3D2P ppend( slice.end(), pointFunctor );

        /*ConstIterator2DP pbegin( xbegin, pointFunctor );
        ConstIterator2DP pend( xend, pointFunctor );

        ConstIterator2DP2 ppbegin( xxbegin, pointFunctor );
        ConstIterator2DP2 ppend( xxend, pointFunctor );*/


        const Dimension size_slice = slice.size();
        std::vector< Statistic< double > > v_statMSEL(size_slice);
//        std::vector< Statistic< double > > v_statMSEL2(size_slice);

        typedef Z2i::Point Point;
        std::vector< Point > pts;

//        if (K.sKCoords(*slice.begin())[0] == 1)
//        trace.warning()<< "Slice begin"<< K.sKCoords(*slice.begin())<<std::endl;
//        do
//        {
////          if (K.sKCoords(*slice.begin())[0] == 1)
////          {
////            trace.info() << *pbegin << " ";
////          }
//            pts.push_back( *ppbegin );
//          ++ppbegin;
//        } while( ppbegin != ppend );
//        if (K.sKCoords(*slice.begin())[0] == 1)
//          trace.info()<< "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"<<std::endl;

//        SCell bel = Surfaces< Functor::KSpace >::findABel( k2d, slice, 10000 );
//        SurfelAdjacency< Functor::KSpace::dimension > SAdj( true );
//        Surfaces< Functor::KSpace >::track2DBoundary( pts, k2d, SAdj, slice, bel );

//        typedef GridCurve<  Functor::KSpace >::PointsRange PointRange;
//        GridCurve<  Functor::KSpace > gridcurve;
//        gridcurve.initFromSCellsVector( pts );
//        PointRange pr = gridcurve.getPointsRange();



        Circulator< ConstIterator3D2P > cbegin( ppbegin, ppbegin, ppend );
        Circulator< ConstIterator3D2P > cend( cbegin );

        for( Dimension ii = 0; ii < size_slice; ++ii )
        {
          v_statMSEL[ii] = Statistic<double>(true);
        }


//        trace.error() << "size " << size_slice << std::endl;
//        if( size_slice <= 4 )
//        {
//          for( Dimension ii = 0; ii < size_slice; ++ii )
//          {
//            v_statMSEL[ii].addValue(0);
//          }
////          analyseAllLengthMS<Functor::KSpace>( v_statMSEL, ppbegin, ppend );
//        }
//        else
        {
          analyseAllLengthMS<Z2i::KSpace>( v_statMSEL, cbegin, cend );
//          Circulator< std::vector< Point >::iterator > cbegin2( pts.begin(), pts.begin(), pts.end() );
//          int shift = 0;//size_slice / 2;
//          for( int i = 0; i < shift; ++i )
//            ++cbegin2;
//          Circulator< std::vector< Point >::iterator > cend2( cbegin2 );
//          analyseAllLengthMS<Functor::KSpace>( v_statMSEL, cbegin2, cend2 );
//          for( Dimension ii = 0; ii < size_slice; ++ii )
//          {
//            v_statMSEL[ii].addValues( v_statMSEL2[ ( ii + shift ) % size_slice ].begin(), v_statMSEL2[ ( ii + shift ) % size_slice ].end());
//          }
        }

        //                for(Dimension ii = 0; ii < pr2size; ++ii )
        //                {
        //                  v_statMSEL[ii].terminate();
        //                }

        Dimension iii = 0;
//        ConstIterator3D sit = slice.begin();
//        ConstIterator3D send = slice.end();
        for( ConstIterator3D sit = slice.begin(), send = slice.end(); sit != send; ++sit )
        {
          Dimension surfel_pos = findSurfel( surfels, *sit );
          ASSERT( surfel_pos != surfels_size );
          if( marque[ surfels[surfel_pos] ] == false )
          {
            marque[ surfels[surfel_pos] ] = true;
            ASSERT(( marque.size() == surfels_size ));
            PairOfStatistics & otherpair = segments[ surfels[surfel_pos] ];
            ASSERT(( segments.size() == surfels_size ));
            if( otherpair.first.samples() == 0 )
            {
              otherpair.first = v_statMSEL[ iii ];
              //              segments[ (Surfel*)&(*sit) ] = otherpair;
            }
            else if ( otherpair.second.samples() == 0 )
            {
              otherpair.second = v_statMSEL[ iii ];
              //              segments[ (Surfel*)&(*sit) ] = otherpair;
            }
            else
            {
              FATAL_ERROR_MSG( false, "ALREADY FILLED" );
              trace.error() << "ALREADY FILLED" << std::endl;
            }
          }
          else
          {
            trace.error() << "WHHHHHAT" << std::endl;
          }
          ++iii;
//          ++sit;
        }// while( sit != send );
      }
    }
  }
}

void checkSizeRadius( double & re,
                      const double h,
                      const double minRadiusAABB,
                      const unsigned int minValue = 5,
                      const double defaultValue = 10.0 )
{
      if( re <= minValue ) /// 	ridiculously small radius check
      {
    //    trace.error() << "re small " << re << std::endl;
    //    re = 5.0 * h;
        re = defaultValue;//*h;
      }
  if( re > ( 0.75 * minRadiusAABB ))
  {
    re = 0.75 * minRadiusAABB;
  }
}


template< typename Surfel >
Dimension computeRadius( std::vector< Surfel > & surfels,
                         std::map< Surfel, std::pair< Statistic< double >, Statistic< double > > > & segments,
                         const double constante,
                         const double h,
                         std::vector< double > & radius,
                         const double minRadiusAABB,
                         const std::string & prop,
                         const std::string & mode,
                         const unsigned int minValue = 5,
                         const double defaultValue = 10.0 )
{
  const Dimension surfels_size = surfels.size();
  std::map< double, unsigned int > nbKernelRadius;

  ASSERT(( radius.size() == surfels_size ));
  ASSERT(( segments.size() == surfels_size ));

  for( Dimension ii = 0; ii < surfels_size; ++ii )
  {
    Statistic< double > & stata = segments[ surfels[ii] ].first;
    Statistic< double > & statb = segments[ surfels[ii] ].second;

    Statistic< double > stat(true);

    stata.terminate();
    statb.terminate();

    //if( mode == "max" )
    {
      if( stata.max() > statb.max() )
      {
        stat = stata;
      }
      else
      {
        stat = statb;
      }
    }
   /* else if( mode == "min" )
    {
      if( stata.max() < statb.max() )
      {
        stat = stata;
      }
      else
      {
        stat = statb;
      }
    }
    else if( mode == "mean" )
    {
      stat.addValues< typename Statistic<double>::Iterator >( stata.begin(), stata.end() );
      stat.addValues< typename Statistic<double>::Iterator >( statb.begin(), statb.end() );
    }
    else
    {
      trace.error() << "I dont understand " << mode << ". I need {min, mean, max} only." << std::endl;
    }*/

    double result = -1.0;
    /*if( prop == "min" )
    {
      result = stat.min();
    }
    else if( prop == "mean" )*/
    {
      result = stat.mean();
//      if( result < minValue )
//      {
//        result = defaultValue;
//      }
    }
   /* else if( prop == "median" )
    {
      result = stat.median();
    }
    else if( prop == "max" )
    {
      result = stat.max();
    }*/

    ASSERT(( result != -1.0 ));

    double re = -1.0;
   /* if( result == defaultValue )
      re = defaultValue;
    else*/
      re = constante * result * result * h;

    checkSizeRadius( re, h, minRadiusAABB );

    radius[ii] = re;

    if( nbKernelRadius.find(re) == nbKernelRadius.end() )
    {
      nbKernelRadius[ re ] = 1;
    }
    else
    {
      nbKernelRadius[ re ] += 1;
    }
  }

  return nbKernelRadius.size();
}


template< typename KSpace, typename Domain, typename Surfel >//= KSpace::Surfel >
bool CheckIfIsInBorder( const KSpace & k, const Domain & d, const Surfel & s, const unsigned int radius = 1 )
{
  typedef typename KSpace::Space::RealPoint MidPoint;
  typedef typename Domain::Point Point;
  typedef CanonicSCellEmbedder< KSpace > Embedder;

  Embedder embedder( k );

  MidPoint mp = embedder.embed( s );
  Point min = d.lowerBound() + Point( radius, radius, radius );
  Point max = d.upperBound() - Point( radius, radius, radius );
  const Dimension dimension = 3;
  for( Dimension ii = 0; ii < dimension; ++ii )
  {
    if( mp[ ii ] < min[ ii ] ||  mp[ ii ] > max [ ii ])
    {
      return true;
    }
  }
  return false;
}

int main( int argc, char** argv )
{
  // parse command line ----------------------------------------------
  po::options_description general_opt("Allowed options are");
  general_opt.add_options()
    ("help,h", "display this message")
    ("input-file,i", po::value< std::string >(), ".vol file")
    ("radius,r",  po::value< double >(), "Kernel radius for IntegralInvariant" )
    ("extension,e", po::value< std::string >()->default_value("vol"), "feiofjreio" )
    ("properties,p", po::value< std::string >()->default_value("mean"), "type of output : mean, gaussian, prindir1 or prindir2 (default mean)");

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
  bool neededArgsGiven=true;
  if (!(vm.count("input-file"))){ 
    missingParam("--input-file");
    neededArgsGiven=false;
  }
  if (!(vm.count("radius"))){ 
    missingParam("--radius");
    neededArgsGiven=false;
  }  
  double h = 1.0;

  FileExtension extension;
  std::string ext = vm["extension"].as< std::string >();
  if( ext == "raw" )
  {
    extension = FileExtension::RAW;
  }
  else if(ext == "vol" )
  {
    extension = FileExtension::VOL;
  }

  bool wrongMode = false;
  std::string mode = vm["properties"].as< std::string >();
  if (( mode.compare("gaussian") != 0 ) && ( mode.compare("mean") != 0 ) && ( mode.compare("prindir1") != 0 ) && ( mode.compare("prindir2") != 0 ))
  {
    wrongMode = true;
  }

  if(!neededArgsGiven ||  wrongMode || !parseOK || vm.count("help") || argc <= 1 )
  {
    trace.info()<< "Visualisation of 3d curvature from .vol file using curvature from Integral Invariant" <<std::endl
                << general_opt << "\n"
                << "Basic usage: "<<std::endl
                << "\t3dCurvatureViewer -i <file.vol> --radius <radius> --properties <\"mean\">"<<std::endl
                << std::endl
                << "Below are the different available properties: " << std::endl
                << "\t - \"mean\" for the mean curvature" << std::endl
                << "\t - \"gaussian\" for the Gaussian curvature" << std::endl
                << "\t - \"prindir1\" for the first principal curvature direction" << std::endl
                << "\t - \"prindir2\" for the second principal curvature direction" << std::endl
                << std::endl;
    return 0;
  }
  double re_convolution_kernel = vm["radius"].as< double >();
  
  // Construction of the shape from vol file
  typedef Z3i::Space::RealPoint RealPoint;
  typedef Z3i::Point Point;
  typedef ImageSelector< Z3i::Domain, bool>::Type Image;
  typedef SimpleThresholdForegroundPredicate< Image > ImagePredicate;
  typedef Z3i::KSpace KSpace;
  typedef KSpace::SCell SCell;
  typedef KSpace::Cell Cell;
  typedef KSpace::Surfel Surfel;
  typedef LightImplicitDigitalSurface< Z3i::KSpace, ImagePredicate > MyLightImplicitDigitalSurface;
  typedef DigitalSurface< MyLightImplicitDigitalSurface > MyDigitalSurface;

  std::string filename = vm["input-file"].as< std::string >();
  Image image = FileReader< Image >( filename, extension );
  {
    Z3i::DigitalSet set_temp( image.domain() );
    SetFromImage< Z3i::DigitalSet >::append< Image >( set_temp, image, 0, 255 );

    image = ImageFromSet< Image >::create< Z3i::DigitalSet >( set_temp, 255, true, set_temp.begin(), set_temp.end() );
  }
  ImagePredicate predicate = ImagePredicate( image, 0 );

  Z3i::Domain domain = image.domain();

  Z3i::KSpace K;

  bool space_ok = K.init( domain.lowerBound(), domain.upperBound(), true );
  if (!space_ok)
  {
    trace.error() << "Error in the Khalimsky space construction."<<std::endl;
    return 2;
  }

  CanonicSCellEmbedder< KSpace > embedder( K );

  SurfelAdjacency< Z3i::KSpace::dimension > SAdj( true );
  Surfel bel = Surfaces< Z3i::KSpace >::findABel( K, predicate, 100000 );
  MyLightImplicitDigitalSurface LightImplDigSurf( K, predicate, SAdj, bel );
  MyDigitalSurface digSurf( LightImplDigSurf );

  typedef DepthFirstVisitor<MyDigitalSurface> Visitor;
  typedef GraphVisitorRange< Visitor > VisitorRange;
  typedef VisitorRange::ConstIterator SurfelConstIterator;
  VisitorRange range( new Visitor( digSurf, *digSurf.begin() ) );
  SurfelConstIterator abegin = range.begin();
  SurfelConstIterator aend = range.end();

  typedef ImageToConstantFunctor< Image, ImagePredicate > MyPointFunctor;
  MyPointFunctor pointFunctor( image, predicate, 1 );

  // Integral Invariant stuff

  typedef FunctorOnCells< MyPointFunctor, Z3i::KSpace > MyCellFunctor;
  MyCellFunctor functor ( pointFunctor, K ); // Creation of a functor on Cells, returning true if the cell is inside the shape

  QApplication application( argc, argv );
  typedef Viewer3D<Z3i::Space, Z3i::KSpace> Viewer;
  //typedef Board3D< Z3i::Space, Z3i::KSpace> Board3D;
  Viewer viewer( K );
  viewer.show();

  viewer << DGtal::SetMode3D(bel.className(), "Basic");

  //    viewer << SetMode3D(image.domain().className(), "BoundingBox") << image.domain();

  VisitorRange range2( new Visitor( digSurf, *digSurf.begin() ) );
  SurfelConstIterator abegin2 = range2.begin();

  trace.beginBlock("Extracting all surfels...");

  std::vector< Surfel > surfels;
  {
    VisitorRange ranger( new Visitor( digSurf, *digSurf.begin() ) );
    SurfelConstIterator abeginer = ranger.begin();
    SurfelConstIterator aender = ranger.end();

    surfels.clear();
    while( abeginer != aender )
    {
      surfels.push_back( *abeginer );
      ++abeginer;
    }
  }

  const Dimension size_surfels = surfels.size();

  trace.endBlock();

  trace.beginBlock("Analyse segments and Mapping segments <-> Surfels...");

  std::map< Surfel, std::pair< Statistic< double >, Statistic< double > > > segments;
  for( Dimension ii = 0; ii < size_surfels; ++ii )
  {
    segments[ surfels[ii] ] = std::pair< Statistic< double >, Statistic< double > >( Statistic< double >( true ), Statistic< double >( true ) );
  }

  computeSegments< MyLightImplicitDigitalSurface, Surfel >( surfels, segments, K, LightImplDigSurf );

  ASSERT(( segments.size() == size_surfels ));

  trace.endBlock();

  trace.beginBlock("Computation of radius...");

  std::vector< double > v_estimated_radius;
  v_estimated_radius.resize( size_surfels );
  //            double global_mean = allSegments.mean();
  //            double global_re = optionsII.constante * global_mean * global_mean * h;
  double min_re = 5;//optionsII.constante * 2 * 2 * h;
//  trace.error() << "min re " << min_re << std::endl;
  Dimension nbKernelsRadius = computeRadius< Surfel >( surfels, segments, 0.1, h, v_estimated_radius, 10, "mean", "max", min_re, -42 );
//  if( nbKernelsRadius < optionsII.nbKernels )
//  {
//    optionsII.nbKernels = nbKernelsRadius;
//  }

  trace.endBlock();

  trace.beginBlock("curvature computation");
  if( ( mode.compare("gaussian") == 0 ) || ( mode.compare("mean") == 0 ) )
  {
    typedef double Quantity;
    std::vector< Quantity > results(size_surfels);
//    back_insert_iterator< std::vector< Quantity > > resultsIterator( results ); // output iterator for results of Integral Invariante curvature computation

    if ( ( mode.compare("mean") == 0 ) )
    {
      typedef IntegralInvariantMeanCurvatureEstimator< Z3i::KSpace, MyCellFunctor > MyIIMeanEstimator;

//      MyIIMeanEstimator estimator ( K, functor );
//      estimator.init( h, re_convolution_kernel ); // Initialisation for a given Euclidean radius of the convolution kernel
//      estimator.eval ( abegin, aend, resultsIterator ); // Computation

#ifdef WITH_OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
      for( Dimension ii = 0; ii < size_surfels; ++ii )
      {
        MyIIMeanEstimator estimator( K, functor );
        estimator.init( h, v_estimated_radius[ii] );
        results[ii] = estimator.eval( surfels.begin() + ii );
      }
    }
    else if ( ( mode.compare("gaussian") == 0 ) )
    {
      typedef IntegralInvariantGaussianCurvatureEstimator< Z3i::KSpace, MyCellFunctor > MyIIGaussianEstimator;
      typedef CurvatureInformations CurvInformation;
      std::vector< CurvInformation > results2;
      back_insert_iterator< std::vector< CurvInformation > > resultsIterator2( results2 );


      MyIIGaussianEstimator estimator ( K, functor );
      estimator.init( h, re_convolution_kernel ); // Initialisation for a given Euclidean radius of the convolution kernel
      estimator.evalPrincipalCurvatures ( abegin, aend, resultsIterator2 ); // Computation

      for(int i = 0; i < results2.size(); ++i)
      {
        results.push_back(results2[i].k1 * results2[i].k2);
      }
    }

    // Drawing results
    VisitorRange range3( new Visitor( digSurf, *digSurf.begin() ) );
    SurfelConstIterator abegin3 = range3.begin();
    Quantity min = numeric_limits < Quantity >::max();
    Quantity max = numeric_limits < Quantity >::min();
    for ( unsigned int ii = 0; ii < results.size(); ++ii )
    {
      Surfel current = *abegin3;
      if( CheckIfIsInBorder< KSpace, Z3i::Domain, Surfel >( K, domain, current, 1 ) )
      {
        results[ ii ] = numeric_limits < Quantity >::min();
      }
      else
      {
        if ( results[ ii ] < min )
        {
          min = results[ ii ];
        }
        else if ( results[ ii ] > max )
        {
          max = results[ ii ];
        }
      }
      ++abegin3;
    }

    typedef GradientColorMap< Quantity > Gradient;
    Gradient cmap_grad( min, max );
    cmap_grad.addColor( Color( 50, 50, 255 ) );
    cmap_grad.addColor( Color( 255, 0, 0 ) );
    cmap_grad.addColor( Color( 255, 255, 10 ) );

    Cell currentCell;

    for ( unsigned int ii = 0; ii < results.size(); ++ii )
    {
      currentCell = Cell( K.sKCoords(*abegin2) );
      if( results[ ii ] != numeric_limits < Quantity >::min() )
      {
        viewer << CustomColors3D( Color::Black, cmap_grad( results[ ii ] ))
               << currentCell;
      }
      else
      {
        viewer << CustomColors3D( Color( 255, 255, 255, 120 ), Color( 255, 255, 255, 120 ) )
               << currentCell;
      }
      ++abegin2;
    }
  }
  else
  {
    typedef double Quantity;
    typedef EigenValues3D< Quantity >::Matrix33 Matrix3x3;
    typedef EigenValues3D< Quantity >::Vector3 Vector3;
    typedef CurvatureInformations CurvInformation;

    std::vector< CurvInformation > results(size_surfels);
    back_insert_iterator< std::vector< CurvInformation > > resultsIterator( results ); // output iterator for results of Integral Invariante curvature computation

    typedef IntegralInvariantGaussianCurvatureEstimator< Z3i::KSpace, MyCellFunctor > MyIIGaussianEstimator;

#ifdef WITH_OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
      for( Dimension ii = 0; ii < size_surfels; ++ii )
      {
        MyIIGaussianEstimator estimator( K, functor );
        estimator.init( h, v_estimated_radius[ii] );
        results[ii] = estimator.evalPrincipalCurvatures( surfels.begin() + ii );
      }
//    MyIIGaussianEstimator estimator ( K, functor );
//    estimator.init ( h, re_convolution_kernel ); // Initialisation for a given Euclidean radius of the convolution kernel
//    estimator.evalPrincipalCurvatures ( abegin, aend, resultsIterator ); // Computation

    trace.endBlock();

    trace.beginBlock("viewer");

    // Drawing results
    typedef  Matrix3x3::RowVector RowVector;
    typedef  Matrix3x3::ColumnVector ColumnVector;

    for ( unsigned int i = 0; i < results.size(); ++i )
    {
      CurvInformation current = results[ i ];
      DGtal::Dimension kDim = K.sOrthDir( *abegin2 );
      SCell outer = K.sIndirectIncident( *abegin2, kDim);
      if ( predicate(embedder(outer)) )
      {
        outer = K.sDirectIncident( *abegin2, kDim);
      }

      Cell unsignedSurfel = K.uCell( K.sKCoords(*abegin2) );

      if( !CheckIfIsInBorder< KSpace, Z3i::Domain, Surfel >( K, domain, *abegin2, 1 ) )
      {
        viewer << CustomColors3D( DGtal::Color(255,255,255,255),
                                  DGtal::Color(255,255,255,255))
               << unsignedSurfel;

        //ColumnVector normal = current.vectors.column(0).getNormalized(); // don't show the normal
        ColumnVector curv1 = current.vectors.column(1).getNormalized();
        ColumnVector curv2 = current.vectors.column(2).getNormalized();

        double eps = 0.01;
        RealPoint center = embedder( outer );// + eps*embedder( *abegin2 );

        //            viewer.addLine ( center[0] - 0.5 * normal[ 0],
        //                             center[1] - 0.5 * normal[1],
        //                             center[2] - 0.5* normal[2],
        //                             center[0] +  0.5 * normal[0],
        //                             center[1] +  0.5 * normal[1],
        //                             center[2] +  0.5 * normal[2],
        //                             DGtal::Color ( 0,0,0 ), 5.0 ); // don't show the normal


        if( ( mode.compare("prindir1") == 0 ) )
        {
          viewer.setLineColor(AXIS_COLOR_BLUE);
          viewer.addLine (
                RealPoint(
                  center[0] -  0.5 * curv1[0],
                  center[1] -  0.5 * curv1[1],
                  center[2] -  0.5 * curv1[2]
                ),
                RealPoint(
                  center[0] +  0.5 * curv1[0],
                  center[1] +  0.5 * curv1[1],
                  center[2] +  0.5 * curv1[2]
                ),
                AXIS_LINESIZE );
        }
        else
        {
          viewer.setLineColor(AXIS_COLOR_RED);
          viewer.addLine (
                RealPoint(
                  center[0] -  0.5 * curv2[0],
                  center[1] -  0.5 * curv2[1],
                  center[2] -  0.5 * curv2[2]
                ),
                RealPoint(
                  center[0] +  0.5 * curv2[0],
                  center[1] +  0.5 * curv2[1],
                  center[2] +  0.5 * curv2[2]
                ),
                AXIS_LINESIZE );
        }
      }
      else
      {
        viewer << CustomColors3D( DGtal::Color(170,170,170,255),
                                  DGtal::Color(170,170,170,255))
               << unsignedSurfel;
      }

      ++abegin2;
    }
    trace.endBlock();
  }

  viewer << Viewer3D<>::updateDisplay;
//  viewer.saveOBJ("/Volumes/Jeremy/snow.obj");
  return application.exec();
}

///////////////////////////////////////////////////////////////////////////////
