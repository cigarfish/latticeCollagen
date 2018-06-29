/*
 * SeedPicker.cpp
 *
 *  Created on: Mar 20, 2012
 *      Author: friebel
 */

#include "SeedPicker.h"

using namespace std;



SeedPicker::SeedPicker(vtkRenderWindowInteractor* iren, int* dim, double* spacing, double* origin)
{
    m_iren = iren;
    for(int i=0; i<3; i++) {
        m_dim[i] = dim[i];
        m_spacing[i] = spacing[i];
        m_origin[i] = origin[i];
    }

    m_seedWidgets.resize(m_dim[2]);
    for(int i=0; i<m_dim[2]; i++) {
        m_seedWidgets[i].seedWidget = vtkSeedWidget::New();
        m_seedWidgets[i].seedRepresentation = vtkSeedRepresentation::New();
        m_seedWidgets[i].handleRepresentation = vtkPointHandleRepresentation2D::New();
        m_seedWidgets[i].seedRepresentation->SetHandleRepresentation(m_seedWidgets[i].handleRepresentation);
        m_seedWidgets[i].seedWidget->SetInteractor(m_iren);
        m_seedWidgets[i].seedWidget->SetRepresentation(m_seedWidgets[i].seedRepresentation);
        m_seedWidgets[i].seedWidget->Off();
    }

    m_yAxisCorr = false;
}


SeedPicker::~SeedPicker(void)
{
    for(int i=0; i<m_dim[2]; i++) {
        m_seedWidgets[i].seedWidget->Delete();
        m_seedWidgets[i].seedRepresentation->Delete();
        m_seedWidgets[i].handleRepresentation->Delete();
    }
}


void SeedPicker::SetInteractor(vtkRenderWindowInteractor *iren)
{
    m_iren = iren;
    for(int i=0; i<m_dim[2]; i++)
        m_seedWidgets[i].seedWidget->SetInteractor(m_iren);
}


void SeedPicker::GetSeedPoints(std::vector<itk::Index<3> >& indices)
{
    indices.clear();
    itk::Index<3> idx;

    double pos [3];
    for(idx[2]=0; idx[2]<m_dim[2]; idx[2]++) {
        for(int j=0; j<m_seedWidgets[idx[2]].seedRepresentation->GetNumberOfSeeds(); j++) {
            m_seedWidgets[idx[2]].seedRepresentation->GetSeedWorldPosition(j, pos);
            idx[0] = (int)( (pos[0] - m_origin[0]) / m_spacing[0]);
            if(!m_yAxisCorr)    idx[1] = (int)( (pos[1] - m_origin[1]) / m_spacing[1]);
            else                idx[1] = (int)(( m_dim[1] - (pos[1] - m_origin[1])) / m_spacing[1]);
            indices.push_back(idx);
        }
    }
}


int SeedPicker::SaveSeedPoints(const char*fname)
{
    std::vector<itk::Index<3> > indices;
    this->GetSeedPoints(indices);

    ofstream ofs(fname);
    if(!ofs.is_open()) {
        cerr << "Warning: cannot open file " + string(fname) << endl;
        return 0;
    }
    for(unsigned int i=0; i<indices.size(); i++)
        ofs << indices[i][0] << " " << indices[i][1] << " " << indices[i][2] << "\n";
    ofs.close();

    cout << "Seed points written to " << fname << endl;

    return 1;
}


void SeedPicker::ClearSeedPoints()
{
    // currently only can clear them by deleting the see widgets
    for(int i=0; i<m_dim[2]; i++) {
        m_seedWidgets[i].seedWidget->Delete();
        m_seedWidgets[i].seedRepresentation->Delete();
        m_seedWidgets[i].handleRepresentation->Delete();
    }

    m_seedWidgets.resize(m_dim[2]);
    for(int i=0; i<m_dim[2]; i++) {
        m_seedWidgets[i].seedWidget = vtkSeedWidget::New();
        m_seedWidgets[i].seedRepresentation = vtkSeedRepresentation::New();
        m_seedWidgets[i].handleRepresentation = vtkPointHandleRepresentation2D::New();
        m_seedWidgets[i].seedRepresentation->SetHandleRepresentation(m_seedWidgets[i].handleRepresentation);
        m_seedWidgets[i].seedWidget->SetInteractor(m_iren);
        m_seedWidgets[i].seedWidget->SetRepresentation(m_seedWidgets[i].seedRepresentation);
        m_seedWidgets[i].seedWidget->Off();
    }
}


void SeedPicker::DisableAll()
{
    for(int i=0; i<m_dim[2]; i++) {
        m_seedWidgets[i].seedWidget->SetEnabled(0);
    }
}


void SeedPicker::EnableSlice(int slice)
{
    for(int i=0; i<m_dim[2]; i++) {
        if(i!=slice)
            m_seedWidgets[i].seedWidget->SetEnabled(0);
        else
            m_seedWidgets[i].seedWidget->SetEnabled(1);
    }
}
