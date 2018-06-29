/*
 * ImageViewer2D.h
 *
 *  Created on: Mar 20, 2012
 *      Author: friebel
 */

#ifndef IMAGEVIEWER2D_H_
#define IMAGEVIEWER2D_H_

#include "vtkImageViewer2.h"

class SeedPicker;
class InteractorObserver;
class vtkTextActor;
class vtkDataSet;

class ImageViewer2D : public vtkImageViewer2
{
public:
  static ImageViewer2D *New();
  vtkTypeRevisionMacro(ImageViewer2D, vtkImageViewer2);
  void PrintSelf(ostream& os, vtkIndent indent);

  void AddDataSet(vtkDataSet*, bool=false);

  virtual void SetupInteractor(vtkRenderWindowInteractor*);

  SeedPicker* GetSeedPicker()
  {
    return seedpicker;
  };

  void SetPickable(bool);
  bool GetPickable();

  vtkSetMacro(Loud, int);
  vtkBooleanMacro(Loud, int);

  vtkSetMacro(NoImage, int);
  vtkBooleanMacro(NoImage, int);

  void SetMode(bool m) { mode = m; };

  void ResetWindowLevel();
  void ResetInformation();

protected:
  ImageViewer2D();
  virtual ~ImageViewer2D();

  SeedPicker *seedpicker;
  InteractorObserver* observer;
  vtkTextActor *TextActor;
  vtkTextActor *StatusTextActor;

  int NoImage;
  int *dim;
  double *range;
  int Loud;
  bool mode;
  char *infoText;
  char *statusInfoText;

private:
  ImageViewer2D(const ImageViewer2D&);  // Not implemented.
  void operator=(const ImageViewer2D&);  // Not implemented.
};

#endif /* IMAGEVIEWER2D_H_ */
