
#include <QGraphicsSceneMoveEvent>

#include "QCS2DImageItem.h"

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageToVTKImageFilter.h"

#include "vtkRenderer.h"
#include "vtkImageMapper.h"
#include "vtkActor2D.h"
#include "vtkGenericOpenGLRenderWindow.h"
#include "vtkRendererCollection.h"

#include "../images/CSImage.h"
#include "../images/filters/imageFilters/exampleITKFilter.h"

#include <iostream>
#include <exception>


QCS2DImageItem::QCS2DImageItem( QGLContext * glCtxt, vtkImageData *data, QGraphicsItem * parent ) :
  QVTKGraphicsItem(glCtxt,parent),
  mWidth(0),mHeight(0)
{
  mpImage = new CSImage();

  if (data)
    {
      mpImage->setImageData(data);

      mWidth  = data->GetDimensions()[0];
      mHeight = data->GetDimensions()[1];

      setMinimumSize(mWidth,mHeight);

      vtkImageMapper * imgMapper = vtkImageMapper::New();
      imgMapper->SetInput(data);
 
      imgMapper->SetColorLevel(138.5);
      imgMapper->SetColorWindow(255);

      vtkActor2D * imgActor = vtkActor2D::New();
      imgActor->SetMapper(imgMapper);

      vtkRenderer * imgRenderer =  mWin->GetRenderers()->GetFirstRenderer();
      if ( ! imgRenderer )
       	{
	  imgRenderer = vtkRenderer::New();
      	  mWin->AddRenderer(imgRenderer);
      	}
      imgRenderer->AddActor(imgActor);

      mWin->AddRenderer(imgRenderer);

      mWin->Render();
    }
}


/*!
  \brief Load image data file.

  \param fileName The name of the file to load

  This function will load the image file in form of an itkImage.
*/
void
QCS2DImageItem::loadFromFile( const char * fileName )
{
  mpImage->loadFromFile(fileName);
}


/*!
  \brief Load image data file.

  \param fileName The name of the file to load

  This function will load the image file in form of an itkImage.
*/
void
QCS2DImageItem::loadFromFile( std::string fileName )
{ loadFromFile( fileName.c_str() ); }


void
QCS2DImageItem::display()
{
  mpImage->update();

  vtkImageData * img = mpImage->getImageData();

  mWidth  = img->GetDimensions()[0];
  mHeight = img->GetDimensions()[1];

  setMinimumSize(mWidth,mHeight);

  vtkImageMapper * imgMapper = vtkImageMapper::New();
  imgMapper->SetInput(img);

  imgMapper->SetColorLevel(138.5);
  imgMapper->SetColorWindow(235);

  vtkActor2D * imgActor = vtkActor2D::New();
  imgActor->SetMapper(imgMapper);

  vtkRenderer * imgRenderer =  mWin->GetRenderers()->GetFirstRenderer();
  if ( ! imgRenderer )
    {
      imgRenderer = vtkRenderer::New();
      mWin->AddRenderer(imgRenderer);
    }
  imgRenderer->AddActor(imgActor);

  mWin->AddRenderer(imgRenderer);

  mWin->Render();
}


/*!
  \brief The bounding rectangle of the item.

  Returns the rectangle with the upper left corner at (0,0).
*/
QRectF
QCS2DImageItem::boundingRect() const
{
  return QRectF( 0, 0, mWidth, mHeight );
}

