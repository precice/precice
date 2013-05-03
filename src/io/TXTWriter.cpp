// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#include "TXTWriter.hpp"
#include "utils/Globals.hpp"
#include "tarch/la/DynamicMatrix.h"
#include "tarch/la/DynamicVector.h"

namespace precice {
namespace io {

tarch::logging::Log TXTWriter:: _log ( "precice::io::TXTWriter" );

TXTWriter:: TXTWriter
(
  const std::string& filename )
:
  _file()
{
  _file.open(filename.c_str());
  if (not _file){
    preciceError("TXTWriter()", "Could not open file \"" << filename
                 << "\" for txt writing!");
  }
  _file.setf ( std::ios::showpoint );
  _file.setf ( std::ios::fixed );
  _file << std::setprecision(16);
}

TXTWriter:: ~TXTWriter()
{
  if (_file){
    _file.close ();
  }
}

//void TXTWriter:: write
//(
//  const tarch::la::DynamicMatrix<double>& matrix )
//  //const std::string&                      filename )
//{
//  //std::ofstream outputStream;
//  //outputStream.open ( _filename.c_str() );
//  //if ( not outputStream ) {
//  //  preciceError ( "TXTWriter()", "Could not open file \"" << filename
//  //                 << "\" for writing matrix!" );
//  //}
//  //_file.setf ( std::ios::showpoint );
//  //_file.setf ( std::ios::fixed );
//  //_file << std::setprecision(16);
//  for ( int i=0; i < matrix.rows(); i++ ){
//    for ( int j=0; j < matrix.cols(); j++ ){
//      _file << matrix(i,j) << " ";
//    }
//    _file << std::endl;
//  }
//  //_file.close ();
//}

//void TXTWriter:: write
//(
//  const tarch::la::DynamicVector<double>& vector )
//  //const std::string&                      filename )
//{
//  //std::ofstream _file;
//  //_file.open ( filename.c_str() );
//  //if ( not _file ) {
//  //  preciceError ( "TXTWriter()", "Could not open file \"" << filename
//  //                 << "\" for writing vector!" );
//  //}
//  //_file.setf ( std::ios::showpoint );
//  //_file.setf ( std::ios::fixed );
//  //_file << std::setprecision(16);
//  for ( int i=0; i < vector.size(); i++ ){
//    _file << vector[i] << " ";
//  }
//  _file << std::endl;
//  //_file.close ();
//}

}} // namespace precice, io
