/*
 * SeedPicker.h
 *
 *  Created on: Mar 20, 2012
 *      Author: friebel
 */

#ifndef SEEDPICKER_H_
#define SEEDPICKER_H_

#include <itkIndex.h>

#include <vtkSeedWidget.h>
#include <vtkSeedRepresentation.h>
#include <vtkPointHandleRepresentation2D.h>
#include <vtkImageActor.h>
#include <vtkRenderWindowInteractor.h>

#include <vector>

class SeedPicker
{
public:

    SeedPicker(vtkRenderWindowInteractor* iren, int* dim, double* spacing, double* origin);
    ~SeedPicker(void);
    void DisableAll();
    void EnableSlice(int slice);
    void GetSeedPoints(std::vector<itk::Index<3> >& indices);
    int  SaveSeedPoints(const char*);
    void ClearSeedPoints();
    void SetInteractor(vtkRenderWindowInteractor*);
    void SetYAxisFlipCorrection(bool y) { m_yAxisCorr = y; };

private:

    struct SeedWidget
    {
        vtkSeedWidget* seedWidget;
        vtkSeedRepresentation* seedRepresentation;
        vtkPointHandleRepresentation2D* handleRepresentation;
    };

    vtkRenderWindowInteractor* m_iren;
    int                        m_dim [3];
    double                     m_spacing [3];
    double                     m_origin [3];
    bool                       m_yAxisCorr;

    std::vector< SeedWidget >  m_seedWidgets;
};

#endif /* SEEDPICKER_H_ */
