////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  File Name:  monolayerDemo.h                                               //
//                                                                            //
//     Author:  Tim Johann <tim.johann@uni-leipzig.de>                        //
//    Created:  2015-02-22 23:23:58                                           //
//                                                                            //
//  This file is part of the CellSys7 code.                                   //
//                                                                            //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig,                         //
//                   and INRIA, Paris-Rocquencourt.                           //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////


#ifndef MONOLAYER_DEMO_H
#define MONOLAYER_DEMO_H

#include <iostream>

#include "../tabMonolayer/monolayer.h"

#include "../../model/Model/ModelCellsSpherical/ModelCellsSpherical.h"
#include "../../tools/parameters/CSParameterContext.h"


class monolayerDemo : public monolayer
{
public:
    monolayerDemo(QWidget *parent=0)
        : monolayer(parent)
        {
            groupBoxSimulation->setEnabled(false);
            groupBoxParameters->setEnabled(false);
            groupBoxScenario->setEnabled(false);
            checkBoxEnableObservation->setChecked(false);
            checkBoxEnableObservation->setEnabled(false);
            checkBoxEnablePOVOutput->setEnabled(false);
            labelOutputObservablesEvery->setEnabled(false);
            doubleSpinBoxObserveEvery->setEnabled(false);
            labelOutputObservablesEveryUnit->setEnabled(false);
            observeCellPopulationSnapshotButton->setEnabled(false);
            buttonWritePovray->setEnabled(false);
            groupBoxVisualizationParameter->setEnabled(false);
            groupBoxPovRayOptions->setEnabled(false);
            saveMXFButton->setEnabled(false);

            ModelCellsSpherical * model = new ModelCellsSpherical();

            model->SetName("Spherical Cells");

            model->mScenario = ModelCellsSpherical::ScenarioSingleCell;

            CSParameterContext * parms = model->GetParameters();

            ModelCellsSpherical::DefaultParameters(parms);

            CSParameter * contactModelParm = parms->findParameter("Contact Model");
            contactModelParm->setData(ModelCellsSpherical::ContactModelJKR);

            CSParameter * parmCellDiameter = parms->findParameter("Cell Diameter");
            if ( parmCellDiameter )
                parmCellDiameter->setData( 1.4e-5 );

            CSParameter * parmDumbBell = parms->findParameter("Use Dumbbell Division");
            if (parmDumbBell)
                parmDumbBell->setData(true);

            CSParameter * parmPolarity = parms->findParameter("Use Cell Polarity");
            if ( parmPolarity )
                parmPolarity->setData(false);

            model->Reset();
            monolayer::initFromModel(model);
        };

};

#endif // MONOLAYER_DEMO_H
