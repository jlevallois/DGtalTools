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
 * @file statisticsEstimators.cpp
 * @ingroup Tools
 * @author Jeremy Levallois (\c jeremy.levallois@liris.cnrs.fr)
 * Laboratoire d'InfoRmatique en Image et Systèmes d'information - LIRIS (CNRS, UMR 5205), CNRS,
 * France
 *
 * @date 2012/06/12
 *
 * Compute statistics between two estimators
 *
 * This file is part of the DGtal tools.
 */

///////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <math.h>
#include <limits>
#include <boost/foreach.hpp>
#include <boost/tokenizer.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

#include "DGtal/base/Common.h"

using namespace DGtal;

/*
int LoadingStringFromFile ( std::string & parseDataFile, std::string & parseParamFile, const std::string filename )
{
    std::ifstream myfile;
    std::stringstream ss_data;
    std::stringstream ss_param;

    myfile.open ( filename.c_str() );
    if ( myfile.is_open() )
    {
        std::string currentLine;
        while ( myfile.good() )
        {
            getline ( myfile, currentLine );
            if (currentLine.size() > 0 && currentLine[0] == '#')
            {
                ss_param << currentLine << '\n';
            }
            else
            {
                ss_data << currentLine << '\n';
            }
        }
        myfile.close();
        parseDataFile = ss_data.str();
        parseParamFile = ss_param.str();
    }
    else
    {
        //std::cout << "estimatorStatistics@LoadingStringFromFile error : Can't open file " << filename << std::endl;
        return 0; //silent exit
    }
    return 1;
}
*/
bool LoadingStringFromFile_0memory( std::ifstream & file, std::string & value )
{
    if( file.good() )
    {
        std::getline( file, value );
        return true;
    }
    return false;
}

int ComputeStatistics_0memory ( const std::string & inputdata1,
                                const std::string & inputdata2,
                                const unsigned int & idColumnData1,
                                const unsigned int & idColumnData2,
                                const bool & isMongeMean,
                                std::ofstream & output )
{
    std::ifstream file1( inputdata1.c_str() );
    std::ifstream file2( inputdata2.c_str() );

    double absd1d2;
    double L1 = 0.0;
    double L2 = 0.0;
    double Linf = 0.0;

    std::string s1, s2;
    double v1, v2;
    double h = - std::numeric_limits<double>::max();

    unsigned int nb_elements = 0;
    while( LoadingStringFromFile_0memory( file1, s1 ) && LoadingStringFromFile_0memory( file2, s2 ))
    {
        if ( s1[ 0 ] == '#' )
        {
            int p = s1.find( "# h = " );
            if ( p != std::string::npos )
            {
                h = atof((s1.erase( p, 5 )).c_str());
                continue;
            }
        }
        if ( s1 == "NA" || s1 == "-nan" || s1 == "-inf" || s1 == "inf" || s1 == "" || s1 == " " )
            continue;
        if ( s2 == "NA" || s2 == "-nan" || s2 == "-inf" || s2 == "inf" || s2 == "" || s2 == " " )
            continue;

        v1 = atof( s1.c_str() );
        v2 = atof( s2.c_str() );

        if( isMongeMean && !(( v1 >= 0.0 ) ^ ( v2 >= 0.0 ))) // hack for Monge. Can be reversed.
        {
            v2 = -v2;
        }

        absd1d2 = std::abs ( v1 - v2 );
        if ( Linf < absd1d2 )
        {
            Linf = absd1d2;
        }
        L1 += absd1d2;
        L2 += absd1d2 * absd1d2;

        ++nb_elements;
    }

    double meanL1 = L1 / (double)nb_elements;
    double meanL2 = ( sqrt ( L2 )) / (double)nb_elements;

    output << h << " "
           << meanL1 << " "
           << meanL2 << " "
           << Linf
           << std::endl;

    return 1; //I Guess (dc)
}
/*
int ComputeStatisticsFromString ( const unsigned int idColumnData1, const unsigned int idColumnData2, const std::string & inputdata, const std::string & inputparam )
{

    std::vector<double> data1, data2;

    boost::char_separator<char> sep_lines("\n");
    boost::char_separator<char> sep_column(" ");

    boost::tokenizer< boost::char_separator<char> > tokens_lines(inputdata, sep_lines);
    BOOST_FOREACH (const std::string & line, tokens_lines)
    {
        boost::tokenizer< boost::char_separator<char> > tokens(line, sep_column);
        boost::tokenizer< boost::char_separator<char> >::iterator beg = tokens.begin();
        boost::tokenizer< boost::char_separator<char> >::iterator current = beg;

        boost::tokenizer< boost::char_separator<char> >::iterator itdata1;
        boost::tokenizer< boost::char_separator<char> >::iterator itdata2;

        for ( int i = 0; i < idColumnData1; ++i )
            ++current;
        itdata1 = current;
        current = beg;
        for ( int i = 0; i < idColumnData2; ++i )
            ++current;
        itdata2 = current;

        std::string cstring = ( (std::string)(*itdata1) );
        std::string cstring2 = ( (std::string)(*itdata2) );

        if ( cstring == "NA" || cstring == "-nan" || cstring == "-inf" || cstring == "inf" || cstring == "" || cstring2 == "NA" || cstring2 == "-nan" || cstring2 == "-inf" || cstring2 == "inf" || cstring2 == "" )
            continue;

        data1.push_back ( atof ( cstring.c_str() ) );
        data2.push_back ( atof ( cstring2.c_str() ) );

    }

    unsigned int sizeVector1 = data1.size();
    unsigned int sizeVector2 = data2.size();
    if ( sizeVector1 == 0 || sizeVector2 == 0 )
    {
        return 0;
    }

    double h = 0.0;
    double radius = 0.0;

    boost::tokenizer< boost::char_separator<char> > tokens_param(inputparam, sep_lines);
    BOOST_FOREACH (const std::string & cline, tokens_param)
    {
        std::string line = cline;
        if ( h != 0.0 )//&& radius != 0.0 )
            break;

        std::string::size_type pos;

        pos = line.find("# h = ");
        if ( pos != std::string::npos )
        {
            h = atof ( (line.erase ( pos, 5 )).c_str() );
            continue;
        }

        //    pos = line.find("# computed kernel radius = ");
        //    if ( pos != std::string::npos)
        //    {
        //      radius = atof ( (line.erase ( pos, 26 )).c_str() );
        //      continue;
        //    }
    }

    if (h == 0.0 )
    {
        std::cout << "estimatorStatistics@ComputeStatisticsFromString error : h param can't found. h = " << h << std::endl;
        return 0;
    }

    if ( sizeVector1 != sizeVector2 )
    {
        std::cout << "estimatorStatistics@ComputeStatisticsFromString error : data1 & data 2 haven't the same size." << sizeVector1 << " " << sizeVector2 << std::endl;
        return 0;
    }

    double absd1d2;

    double L1 = 0.0;
    double L2 = 0.0;
    double Linf = 0.0;

    for ( int index = 0; index < sizeVector1; ++index )
    {
        absd1d2 = std::abs ( (double)( data1[index] - data2[index] ));
        if ( Linf < absd1d2 )
            Linf = absd1d2;
        if( Linf > 10 )
            std::cout << "here " << data1[index] << " " << data2[index] << std::endl;
        L1 += absd1d2;
        L2 += absd1d2 * absd1d2;
    }

    std::cout
            << "# h | "
            << "L1 Mean Error | "
            << "L2 Mean Error | "
            << "Loo Mean Error"
            << std::endl;

    double meanL1 = L1 / (double)sizeVector1;
    double meanL2 = ( sqrt ( L2 )) / (double)sizeVector1;
    std::cout
            << h << " "
            << radius << " "
            << meanL1 << " "
            << meanL2 << " "
            << Linf
            << std::endl;

    return 1;
}
*/

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
    po::options_description general_opt("Allowed options are");
    general_opt.add_options()
            ("help,h", "display this message")
            ("file1,f", po::value< std::string >(), "File 1")
            ("file2,F", po::value< std::string >(), "File 2")
            ("column1,c",  po::value< unsigned int >(), "Column of file 1" )
            ("column2,C",  po::value< unsigned int >(), "Column of file 2" )
            ("output,o", po::value< std::string >(), "Output file")
            ("monge,m",  po::value< bool >()->default_value( false ), "Is from Monge mean computation (optional)" );


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
        trace.info()<< "Compute satistics (L1, L2, Loo) from results of two extimators" <<std::endl
                    << "Basic usage: "<<std::endl
                    << "\tstatisticsEstimators_0memory --file1 <file1> --column1 <column1> --file2 <file2> --column2 <column2> --output <output>"<<std::endl
                    << std::endl;

        return 0;
    }


    if (!(vm.count("file1"))) missingParam("--file1");
    if (!(vm.count("file2"))) missingParam("--file2");
    if (!(vm.count("column1"))) missingParam("--column1");
    if (!(vm.count("column2"))) missingParam("--column2");
    if (!(vm.count("output"))) missingParam("--output");

    std::string filename1 = vm["file1"].as< std::string >();
    std::string filename2 = vm["file2"].as< std::string >();
    unsigned int column1 = vm["column1"].as< unsigned int >();
    unsigned int column2 = vm["column2"].as< unsigned int >();
    std::string output_filename = vm["output"].as< std::string >();
    bool isMongeMean = vm["monge"].as< bool >();

    std::ofstream file( output_filename.c_str(), std::ofstream::out | std::ofstream::app );
    file.flags( std::ios_base::unitbuf );

    if( file.eof() )
    {
        file << "# h | "
             << "L1 Mean Error | "
             << "L2 Mean Error | "
             << "Loo Mean Error"
             << std::endl;
    }

    if ( ComputeStatistics_0memory( filename1, filename2, column1, column2, isMongeMean, file ) == 0 )
    {
        file.close();
        return -1;
    }

    file.close();
    return 1;

    /*if ( LoadingStringFromFile ( inputData, inputParam, filename ) == 0 )
    return -2;
  if ( ComputeStatisticsFromString ( column1, column2, inputData, inputParam ) == 0 )
    return -3;*/
}
