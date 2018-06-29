////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  File Name:  CSHVTPWriter.h                                                //
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


#ifndef CS_VTP_WRITER_H
#define CS_VTP_WRITER_H

#include "../../model/Model/ModelCellsSpherical/ModelCellsSpherical.h"
#include <string>
#include<iostream>
#include<fstream>

class CSModel;
class ModelCellsSpherical;

class CSVTPWriter
{
public:
  CSVTPWriter( const std::string & outputFileName , int limit=10000000);
  ~CSVTPWriter();

  void reset();

  void setOutputRadius( bool outputRadius );
  void setOutputQuiescence( bool outputQuiescence );
  void setOutputConcentrations( bool outputConcentrations );
  void setOutputGradients( bool outputGradients );
  void setOutputPressure( bool outputPressure );
  void setOutputPolarVector( bool outputPolarVector );

  void setOutput_sphericalCells_1D_cut( double cutLength );

  void exec( ModelCellsSpherical *model );

private:

    std::string mOutputFileName;
    int mLimit;
    int mCurrent;

    bool mOutputRadius;
    bool mOutputQuiescence;
    bool mOutputConcentrations;
    bool mOutputGradients;
    bool mOutputPressure;
    bool mOutputPolarVector;

    bool   mOutput_sphericalCells_1D_cut;
    double mOutput_sphericalCells_1D_cutLength;

};

#endif // CS_VTP_WRITER_H
