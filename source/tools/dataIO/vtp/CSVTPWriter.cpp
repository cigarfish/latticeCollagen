////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  File Name:  CSHVTPWriter.cpp                                              //
//                                                                            //
//     Author:  Johannes Neitsch <johannes@neitsch.de>                        //
//    Created:  2014-11-27 13:34:28                                           //
//                                                                            //
//  This file is part of the CellSys7 code.                                   //
//                                                                            //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig,                         //
//                   and INRIA, Paris-Rocquencourt.                           //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include <string>
#include <iostream>
#include <fstream>

#include <cstdio>

#include "CSVTPWriter.h"

#include "../../model/Model/ModelCellsSpherical/ModelCellsSpherical.h"

CSVTPWriter::CSVTPWriter( const std::string & fileName , int limit)
    : mOutputFileName( fileName ),
      mLimit( limit ),
      mOutputRadius( false ),
      mOutputQuiescence( false ),
      mOutputPressure( false ),
      mOutput_sphericalCells_1D_cut( false ),
      mOutput_sphericalCells_1D_cutLength( 0 ),
      mCurrent( 0 )
{

}

CSVTPWriter::~CSVTPWriter()
{
  //  mpHDF5File->close();
 //   delete mpHDF5File;
}

void CSVTPWriter::reset(){
  this->mCurrent = 0.;
}

void CSVTPWriter::setOutputRadius( bool outputRadius ){
  this->mOutputRadius = outputRadius;
}
void CSVTPWriter::setOutputQuiescence( bool outputQuiescence ){
  this->mOutputQuiescence = outputQuiescence;
}
void CSVTPWriter::setOutputConcentrations( bool outputConcentrations ){
  this->mOutputConcentrations = outputConcentrations;
}
void CSVTPWriter::setOutputGradients( bool outputGradients ){
  this->mOutputGradients = outputGradients;
}
void CSVTPWriter::setOutputPressure( bool outputPressure ){
  this->mOutputPressure = outputPressure;
}
void CSVTPWriter::setOutputPolarVector( bool outputPolarVector ){
  this->mOutputPolarVector = outputPolarVector;
}

void CSVTPWriter::setOutput_sphericalCells_1D_cut( double cutLength ){

  if( cutLength < 0 ){
    this->mOutput_sphericalCells_1D_cut = false;
  }
  else{
    this->mOutput_sphericalCells_1D_cut = true;
    this->mOutput_sphericalCells_1D_cutLength = cutLength;

  }

}

void CSVTPWriter::exec( ModelCellsSpherical *model )
{

  std::ofstream vtp;
  std::ostringstream outStream;
  outStream<< setfill('0') << setw(7) << mCurrent;
  std::string vtpOutputFileName=mOutputFileName+outStream.str()+".vtp";

  vtp.open(vtpOutputFileName.c_str(),std::ios::out|std::ios::app);

  vtp << "<?xml version=\"1.0\" encoding=\"utf-8\"?>" << std::endl;
  vtp << "<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"BigEndian\">" << std::endl;
  vtp << "<PolyData GhostLevel=\"0\">" << std::endl;

  vtp << "<Piece NumberOfPoints=\""
          << model->cells.size()
          << "\" NumberOfVerts=\""
          << 0
          << "\" NumberOfLines=\""
          << 0
          << "\" NumberOfStrips=\"0\" NumberOfPolys=\""
          << 0
          << "\">" << std::endl;

  vtp << "<Points>" << std::endl;
  vtp << "<DataArray Name=\"xyz\" NumberOfComponents=\"3\" type=\"Float64\" format=\"ascii\">" << std::endl;

    
  if( this->mOutput_sphericalCells_1D_cut == true && model->dimension == 1 && model->cells.size() != 0 ){

    double x_left  = model->cells[0]->position.x;
    double x_right = model->cells[0]->position.x;

    for( unsigned int c = 1 ; c < model->cells.size() ; c++ ){
      if( model->cells[c]->position.x < x_left )
        x_left = model->cells[c]->position.x;
      if( model->cells[c]->position.x > x_right )
        x_right = model->cells[c]->position.x;
    }

    double x_middle = 0.5 * ( x_left + x_right );



    for( unsigned int i = 0 ; i < model->cells.size() ; i++ ){

      double x_tmp = model->cells[i]->position.x - x_middle;

      int shift_rows = x_tmp/this->mOutput_sphericalCells_1D_cutLength;

      if( shift_rows >= 0 ){
        if( x_tmp > 0 )
          shift_rows += 1;
        vtp << x_tmp - (shift_rows+0.5)*this->mOutput_sphericalCells_1D_cutLength << " ";
      }
      else
        vtp << x_tmp - (shift_rows+0.5)*this->mOutput_sphericalCells_1D_cutLength << " ";
      vtp << shift_rows*model->defaultInitialCellRadius * 2.3 << " "
          << model->cells[i]->position.z << std::endl;
    }
  }else{
    for( unsigned int i = 0 ; i < model->cells.size() ; i++ )
    vtp << model->cells[i]->position.x << " " <<
           model->cells[i]->position.y << " "<<
           model->cells[i]->position.z << std::endl;
  }

  vtp << "</DataArray>" << std::endl;
  vtp << "</Points>"    << std::endl;
  vtp << "<PointData>"  << std::endl;

  //cell index in cells container
    vtp << "<DataArray Name=\"Cell index\" NumberOfComponents=\"1\" type=\"Int32\" format=\"ascii\">" << std::endl;
    for( unsigned int i = 0 ; i < model->cells.size() ; i++ )
      vtp << i << std::endl;
    vtp << "</DataArray>" << std::endl;

//radius
  if( this->mOutputRadius ){
    vtp << "<DataArray Name=\"Radius\" NumberOfComponents=\"1\" type=\"Float64\" format=\"ascii\">" << std::endl;
    for( unsigned int i = 0 ; i < model->cells.size() ; i++ )
      vtp << model->cells[i]->mRadius << std::endl;
    vtp << "</DataArray>" << std::endl;
  }

//proliferation status
  if( this->mOutputQuiescence){
    vtp << "<DataArray Name=\"Quiescence\" NumberOfComponents=\"1\" type=\"Int32\" format=\"ascii\">" << std::endl;
    for( unsigned int i = 0 ; i < model->cells.size() ; i++ )
      vtp << model->cells[i]->getState( Cell::StateQuiescent ) << std::endl;
    vtp << "</DataArray>" << std::endl;
  }

//pressure of cells
  if( this->mOutputPressure){

    double correct = model->biolink->energy_scale/(model->biolink->length_scale*model->biolink->length_scale*model->biolink->length_scale);

    vtp << "<DataArray Name=\"Pressure\" NumberOfComponents=\"1\" type=\"Float64\" format=\"ascii\">" << std::endl;
    for( unsigned int i = 0 ; i < model->cells.size() ; i++ )
      vtp << model->cells[i]->lastPressure * correct << std::endl;
    vtp << "</DataArray>" << std::endl;


  }


  if ( this->mOutputPolarVector )
  {
      vtp << "<DataArray Name=\"PolarVector\" NumberOfComponents=\"3\" type=\"Float64\" format=\"ascii\">" << std::endl;
      for( auto cell: model->cells )
      {
          if ( cell->mType == ModelElement::TypeCellSphericalPolar )
              vtp << static_cast<CellSphericalPolar *>(cell)->mPolarDirection.x << " "
                  << static_cast<CellSphericalPolar *>(cell)->mPolarDirection.y << " "
                  << static_cast<CellSphericalPolar *>(cell)->mPolarDirection.z << " "
                  << std::endl;
          else {
              vtp << "0 0 0" << std::endl;
          }
      }
      vtp << "</DataArray>" << std::endl;
  }


  vtp << "</PointData>" << std::endl;
  vtp << "</Piece>"     << std::endl;
  vtp << "</PolyData>"  << std::endl;
  vtp << "</VTKFile>";

  vtp.close();

  mCurrent++;

}
