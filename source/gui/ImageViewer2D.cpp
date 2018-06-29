/*
 * ImageViewer2D.cpp
 *
 *  Created on: Mar 20, 2012
 *      Author: friebel
 */

#include "ImageViewer2D.h"

#include "SeedPicker.h"

#include <vtkObjectFactory.h>

#include <vtkCommand.h>
#include <vtkCamera.h>
#include <vtkImageActor.h>
#include <vtkImageData.h>
#include <vtkImageMapToWindowLevelColors.h>
#include <vtkInteractorStyleImage.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkTextActor.h>
#include <vtkTextProperty.h>

#include <vtkDataSet.h>
#include <vtkDataSetMapper.h>
#include <vtkActor.h>
#include <vtkProperty.h>

#include <vtkSmartPointer.h>


using namespace std;


class InteractorObserver: public vtkCommand
{
public:
    static InteractorObserver* New()
    {
        return new InteractorObserver;
    };

    void Setup(ImageViewer2D *g, SeedPicker *sp)
    {
        this->viewer = g;
        this->seedpicker = sp;
        this->Picking = 0;
        this->Pickable = 1;
    };

    void SetPickable(bool p)
    {
        this->Pickable = p;
    };

    bool GetPickable()
    {
        return this->Pickable;
    };

    void Execute(vtkObject* caller, unsigned long eventId, void* callData)
    {
        vtkRenderWindowInteractor *iren = static_cast<vtkRenderWindowInteractor*> (caller);

        char code = iren->GetKeyCode();
        string keysym = iren->GetKeySym();
        int slice = this->viewer->GetSlice();

        if(keysym == "less" || keysym == "Down" || code == '<') {
            if(slice > viewer->GetSliceMin())
                slice--;
        }
        else if(keysym == "greater" || keysym == "Up" || code == '>') {
            if(slice < viewer->GetSliceMax())
                slice++;
        }
        else if(code == 's') {
            Picking = !Picking;
            if(Pickable && !Picking)
                seedpicker->DisableAll();
        }
//        else if(code == 'c') {
//            this->viewer->ResetWindowLevel();
//        }
//        else if(code == 'p') {
//            this->seedpicker->SaveSeedPoints("pts.txt");
//        }
        else if(code == 'd') {
            this->seedpicker->ClearSeedPoints();
        }
        else {
            return;
        }

        this->viewer->SetSlice(slice);
        this->viewer->SetMode(!Picking);
        this->viewer->ResetInformation();
        this->viewer->Render();

        if(Pickable && Picking)
            seedpicker->EnableSlice(slice);
    };

protected:
    InteractorObserver()
{
        viewer = 0;
};

    virtual ~InteractorObserver() {};

private:
    ImageViewer2D *viewer;
    SeedPicker *seedpicker;
    bool Picking, Pickable;
};



vtkCxxRevisionMacro(ImageViewer2D, "$Revision: 0.0 $")

vtkStandardNewMacro(ImageViewer2D)


ImageViewer2D::ImageViewer2D()
{
    this->Loud = 0;
    this->NoImage = 0;
    this->seedpicker = 0;
    this->observer = 0;
    this->TextActor = 0;
    this->StatusTextActor = 0;
    this->infoText = new char[200];
    this->statusInfoText = new char[200];
    this->mode = true;
}


ImageViewer2D::~ImageViewer2D()
{
    if(Loud)
        cout << "ImageViewer::~ImageViewer" << endl;
    if(this->seedpicker) {
        delete this->seedpicker;
        this->seedpicker = 0;
    }
    if(this->observer) {
        this->observer->Delete();
        this->observer = 0;
    }
    if(this->TextActor) {
        this->TextActor->Delete();
        this->TextActor = 0;
    }
    if(this->StatusTextActor) {
        this->StatusTextActor->Delete();
        this->StatusTextActor = 0;
    }
    if(infoText) {
        delete[] this->infoText;
        this->infoText = 0;
    }
    if(statusInfoText) {
        delete[] this->statusInfoText;
        this->statusInfoText = 0;
    }
}


void ImageViewer2D::ResetInformation()
{
    if(Loud)
        cout << "ImageViewer::ResetInformation" << endl;

    vtkImageData *image = this->GetInput();
    if(NoImage || !image)
        return;

    int *dim = image->GetDimensions();
    double *range = image->GetScalarRange();
    double *spacing = image->GetSpacing();
    double *origin = image->GetOrigin();

    sprintf(this->infoText,
            "Range: (%g, %g) \nDimensions: %d %d %d \nOrigin: %g %g %g \nSpacing: %g %g %g \nSlice: %d",
            range[0], range[1], dim[0], dim[1], dim[2], origin[0], origin[1], origin[2], spacing[0], spacing[1],
            spacing[2], this->GetSlice());

    if(mode) {
        sprintf(this->statusInfoText, "Please hit 's' to enter seed point selection mode.\nUse the arrow up and down keys to navigate the stack in z-direction.");
        this->StatusTextActor->GetTextProperty()->SetColor(1, 0, 0);
    }
    else {
        sprintf(this->statusInfoText, "Please place one seed point in the cavity of each vein.\nUse 'd' to delete the placed seed points.\nPress 's' after selection and confirm your selection by hitting the button in the right panel.");
        this->StatusTextActor->GetTextProperty()->SetColor(0, 1, 0);
    }

    this->TextActor->SetInput(this->infoText);
    this->StatusTextActor->SetInput(this->statusInfoText);
}


void ImageViewer2D::ResetWindowLevel()
{
    if(Loud)
        cout << "ImageViewer::ResetWindowLevel" << endl;

    this->GetInput()->Update(); // ?

    vtkImageData *image = this->GetInput();
    this->range = image->GetScalarRange();
    this->SetColorLevel(0.5 * (range[0] + range[1]));
    this->SetColorWindow(range[1] - range[0]);
}


void ImageViewer2D::SetupInteractor(vtkRenderWindowInteractor *rwi)
{
    if(Loud)
        cout << "ImageViewer::SetupInteractor" << endl;

    vtkImageData *image = this->GetInput();
    if(!NoImage && image) {
        this->Superclass::SetupInteractor(rwi);

        this->Render();
        this->GetRenderer()->ResetCamera();
        this->GetRenderer()->GetActiveCamera()->SetClippingRange(0.1, 100000);

        this->TextActor = vtkTextActor::New();
        this->TextActor->SetDisplayPosition(10, 10);
        this->TextActor->GetTextProperty()->SetJustificationToLeft();
        this->TextActor->GetTextProperty()->SetOpacity(0.7);
        this->TextActor->GetTextProperty()->SetColor(1, 1, 1);

        this->StatusTextActor = vtkTextActor::New();
        this->StatusTextActor->SetDisplayPosition(200, 10);
        this->StatusTextActor->GetTextProperty()->SetJustificationToLeft();
        this->StatusTextActor->GetTextProperty()->SetOpacity(0.7);
        this->StatusTextActor->GetTextProperty()->SetColor(1, 0, 0);

        this->ResetInformation();
        this->GetRenderer()->AddActor(this->TextActor);
        this->GetRenderer()->AddActor(this->StatusTextActor);

        this->ImageActor->InterpolateOff();

        this->ResetWindowLevel();

        this->seedpicker = new SeedPicker(rwi, image->GetDimensions(), image->GetSpacing(), image->GetOrigin());
        this->observer = InteractorObserver::New();
        observer->Setup(this, seedpicker);
        rwi->AddObserver(vtkCommand::KeyPressEvent, observer);

    }
    else if(this->ImageActor) { // no image to show
        this->ImageActor->VisibilityOff();
    }
}


void ImageViewer2D::AddDataSet(vtkDataSet *ds, bool noimage)
{
    if(Loud)
        cout << "ImageViewer::AddDataSet" << endl;

    NoImage = noimage;

    vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
    vtkSmartPointer<vtkDataSetMapper> mapper = vtkSmartPointer<vtkDataSetMapper>::New();

    actor->GetProperty()->EdgeVisibilityOn();
    actor->GetProperty()->SetLineWidth(1);
    actor->GetProperty()->SetEdgeColor(0,0,1);
    mapper->SetInput(ds);
    actor->SetMapper(mapper);
    this->GetRenderer()->AddActor(actor);

    if(NoImage && this->ImageActor) {
        vtkSmartPointer<vtkImageData> image = vtkSmartPointer<vtkImageData>::New();

        this->SetInput(image);
        this->ImageActor->VisibilityOff();
    }
}


void ImageViewer2D::SetPickable(bool p)
{
    this->observer->SetPickable(p);
}


bool ImageViewer2D::GetPickable()
{
    return this->observer->GetPickable();
}


void ImageViewer2D::PrintSelf(ostream& os, vtkIndent indent) {
    this->Superclass::PrintSelf(os, indent);
}
