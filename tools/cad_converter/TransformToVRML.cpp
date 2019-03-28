#include "TransformToVRML.hpp"

#include <iostream>
#include "IGES2Vrml.hpp"
#include "STEP2Vrml.hpp"
#include "boost/algorithm/string.hpp"
#include <vector>
#include <boost/unordered_map.hpp>

using namespace boost;

TransformToVRML:: TransformToVRML
(
  const std::string & inputFileName,
  const std::string & outputFileName )
:
  _inputFileName ( inputFileName ),
  _outputFileName ( outputFileName ),
  _inputFileType ( UNDEFINED )
{}

bool TransformToVRML:: parseFileType ()
{
  if ( (_inputFileName.length() < 1) || (_outputFileName.length() < 1) ) {
    return false;
  }

  std::string suffix; // file type:  .igs .iges .stp .step

  suffix = _inputFileName.substr ( _inputFileName.length() - 4 );
  boost::to_lower( suffix );
  if ( suffix == ".igs") {
    _inputFileType = IGES;
  }
  else if ( suffix == ".stp") {
    _inputFileType = STEP;
  }
  else if ( suffix == ".wrl" ) {
    _inputFileType = VRML;
  }
  else {
    suffix = _inputFileName.substr ( _inputFileName.length() - 5 );
    boost::to_lower( suffix );
    if ( suffix == ".iges") {
      _inputFileType = IGES;
    }
    else if ( suffix == ".step") {
      _inputFileType = STEP;
    }
  }
  return _inputFileType != UNDEFINED;
}

bool TransformToVRML:: doTransform ()
{
  if ( _inputFileType == IGES ) {
    transformIGES ();
  }
  else if ( _inputFileType == STEP ) {
    transformSTEP ();
  }
  else if ( _inputFileType == VRML ) {
    // do nothing
  }
  else {
    return false;
  }
  eliminateRedundancy ();
  return true;
}

bool TransformToVRML:: transformIGES ()
{
  IGES2Vrml convertIGES;
  if ( ! convertIGES.readFile(_inputFileName) ) {
    return false;
  }
  convertIGES.loadReport ();
  return convertIGES.transfer ( _outputFileName );
}

bool TransformToVRML:: transformSTEP ()
{
  STEP2Vrml convertSTEP;
  if ( convertSTEP.readFile(_inputFileName) ) {
    return false;
  }
  convertSTEP.loadReport ();
  return convertSTEP.transfer ( _outputFileName );
}

void TransformToVRML:: eliminateRedundancy ()
{
  std::fstream fp;
  char buf[1024] = {0};
  std::vector<std::string> stringCoordinate3;
  std::vector<std::string> vrmlText;

  fp.open( _outputFileName.c_str(), std::ios_base::in );

  while ( fp.getline(buf, sizeof(buf)) ) {
    std::string str = buf;
    vrmlText.push_back ( str );
    if ( str.find("point [") != std::string::npos ) {
      break;
    }
  }

  while ( fp.getline ( buf, sizeof ( buf ) ) ) {
    std::string str = buf;
    // not last line
    if ( str.find(']') == std::string::npos ) {
      str.erase ( str.length()-1 );
      stringCoordinate3.push_back ( str );
    }
    // last line
    else {
      str.erase( str.length()-1 );
      str.erase( str.length()-1 );
      stringCoordinate3.push_back( str );
      break;
    }
  }

  // count redundancy
  typedef boost::unordered_map<std::string,int> Map;
  Map oneMap;
  std::vector<int> newOrder;
  std::vector<int> deleteList;

  std::cout << "before elemination: " << stringCoordinate3.size() << " points "
            << '\n';

  size_t j = 0;
  for ( size_t i=0; i < stringCoordinate3.size (); i++ ) {
    Map::iterator iter;
    iter = oneMap.find ( stringCoordinate3 [ i ] );
    if ( iter == oneMap.end() ) {
      oneMap.insert ( Map::value_type(stringCoordinate3[i], j) );
      newOrder.push_back ( j++ );
    }
    else {
      newOrder.push_back ( iter->second );
      deleteList.push_back ( i );
    }
  }

  std::cout << "after elemination: " << oneMap.size() << " points \n";

  // eliminate redundancy
  for ( size_t i=deleteList.size ()-1; i >= 0; i-- ) {
    stringCoordinate3.erase ( stringCoordinate3.begin() + deleteList[i] );
  }

  while ( fp.getline(buf, sizeof(buf)) ) {
    std::string str = buf;
    if ( str.find("IndexedFaceSet") != std::string::npos ) {
      fp.getline ( buf, sizeof(buf) );
      break;
    }
  }

  int s0, s1, s2;
  std::vector<int> vectorFaceIndices;
  char kommer;
  std::string str;

  do {
    fp >> s0 >> kommer >> s1 >> kommer >> s2 ;

    vectorFaceIndices.push_back( s0 );
    vectorFaceIndices.push_back( s1 );
    vectorFaceIndices.push_back( s2 );

    fp.getline ( buf, sizeof ( buf ) );
    str = buf;

  } while ( str.compare(",-1,") == 0  );


  for ( size_t i=0; i < vectorFaceIndices.size(); i++ ) {
    vectorFaceIndices[i] = newOrder[vectorFaceIndices[i]];
  }
  fp.close ();

  fp.open ( _outputFileName.c_str(), std::ios_base::out );

  for ( size_t i=0; i < vrmlText.size(); i++ ) {
    fp << vrmlText [ i ] << '\n';
  }

  for ( size_t i=0; i < stringCoordinate3.size()-1; i++ ) {
    fp << stringCoordinate3[i] << ",\n";
  }

  fp << stringCoordinate3[stringCoordinate3.size()-1] << " ]\n";

  fp << "}\n";
  fp << "ShapeHints {\n";
  fp << "}\n";
  fp << "IndexedFaceSet {\n";
  fp << "coordIndex [\n";

  for ( size_t i=0; i < vectorFaceIndices.size() - 1; i += 3 ) {
    fp << vectorFaceIndices [ i + 0 ] << ","
    << vectorFaceIndices [ i + 1 ] << ","
    << vectorFaceIndices [ i + 2 ] << ","
    << "-1,\n";
  }

  fp << vectorFaceIndices [ vectorFaceIndices.size() - 3 ] << ","
  << vectorFaceIndices [ vectorFaceIndices.size() - 2 ] << ","
  << vectorFaceIndices [ vectorFaceIndices.size() - 1 ] << ","
  << "-1\n";

  fp << "]\n";
  fp << "}\n";
  fp << "}\n";
  fp << "}\n";
  fp << "}\n";
  fp << "}\n";

  return;
}











