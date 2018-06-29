///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  vtkPOVExporterWireFrame.h                                            //
//                                                                                   //
//     Author:  Johannes Neitsch <johannes.neitsch@uni-leipzig.de>                   //
//    Created:  2013-08-30                                                           //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

#ifndef __vtkPOVExporterwireframe_h
#define __vtkPOVExporterwireframe_h

#include "vtkPOVExporter.h"

class vtkRenderer;
class vtkActor;
class vtkCamera;
class vtkLight;
class vtkPolyData;
class vtkProperty;
class vtkTexture;
class vtkPOVInternalsWireFrame;

class vtkPOVExporterWireFrame : public vtkPOVExporter
{
public:
    static vtkPOVExporterWireFrame *New();
    vtkTypeMacro(vtkPOVExporterWireFrame, vtkPOVExporter);
    void PrintSelf(ostream& os, vtkIndent indent);

    //Description:
    //The filename to save into. 
    vtkSetStringMacro(FileName);
    vtkGetStringMacro(FileName);

protected:
    vtkPOVExporterWireFrame();
    ~vtkPOVExporterWireFrame();
    
    void WriteData();
    virtual void WriteHeader(vtkRenderer *renderer);
    void WriteCamera(vtkCamera *camera);
    void WriteLight(vtkLight *light);
    void WriteProperty(vtkProperty *property);
    void WritePolygons(vtkPolyData *polydata, bool scalar_visible);
    void WriteTriangleStrips(vtkPolyData *strip, bool scalar_visible);

    virtual void WriteActor(vtkActor *actor);

    char *FileName;
    FILE *FilePtr;

private:
    vtkPOVExporterWireFrame(const vtkPOVExporterWireFrame&);  // Not implemented.
    void operator=(const vtkPOVExporterWireFrame&);  // Not implemented.    

    vtkPOVInternalsWireFrame *Internals;
};

#endif