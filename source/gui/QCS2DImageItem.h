#ifndef QCS_2D_IMAGEITEM_H
#define QCS_2D_IMAGEITEM_H

#include <string>

#include <QGraphicsItem>
#include "QVTKGraphicsItem.h"

#include "itkRGBPixel.h"
#include "itkImage.h"
#include "itkImageToVTKImageFilter.h"

class QGLContext;
class CSImage;

typedef itk::Image< itk::RGBPixel<unsigned char>, 2 >   itkImageType;
typedef itk::ImageToVTKImageFilter<itkImageType>      itkToVTKFilter;


class QCS2DImageItem : public QVTKGraphicsItem
{
 public:

  QCS2DImageItem( QGLContext * glCtxt, vtkImageData *data=0, QGraphicsItem * parent=0 );

  void loadFromFile( const char * fileName );
  void loadFromFile( std::string fileName );

  QRectF boundingRect() const;

  void display();
  
 private:

  int mWidth, mHeight;
  CSImage * mpImage;
  /* vtkImageData * mpVTKImageData; */
  /* itkImageType::Pointer mpITKImage; */
  /* itkToVTKFilter::Pointer mpConnector; */
};

#endif // QCS_2D_IMAGEITEM_H
