///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  GraphEdgePythonSVMClassifier.tpp                                     //
//                                                                                   //
//     Author:  Adrian Friebel <adrian.friebel@uni-leipzig.de>                       //
//    Created:  2014-08-08                                                           //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

#include "GraphEdgePythonSVMClassifier.h"

#include <vtkDataSetAttributes.h>
#include <vtkDoubleArray.h>
#include <vtkEdgeListIterator.h>
#include <vtkFloatArray.h>
#include <vtkIntArray.h>
#include <vtkOutEdgeIterator.h>
#include <vtkUnsignedLongArray.h>



template< unsigned int VImageDimension > GraphEdgePythonSVMClassifier< VImageDimension >::GraphEdgePythonSVMClassifier()
    : GraphEdgeBaseClassifier< VImageDimension >()
{
    // TODO Auto-generated constructor stub
}


template< unsigned int VImageDimension > GraphEdgePythonSVMClassifier< VImageDimension >::~GraphEdgePythonSVMClassifier()
{
//#ifdef CS_BUILD_PYTHON
//    Py_XDECREF(mpTrainDBHandlerInstance);
//    Py_XDECREF(mpSVMHandlingInstance);
//    Py_Finalize();
//#endif
}


template< unsigned int VImageDimension > void GraphEdgePythonSVMClassifier< VImageDimension >::Init()
{
    std::cout << "GraphEdgePythonSVMClassifier::Init here" << std::endl;
    GraphEdgeBaseClassifier<VImageDimension>::Init();
    std::cout << "GraphEdgeBaseClassifier::Init finished" << std::endl;
    SetupClassifier();
}


template< unsigned int VImageDimension > void GraphEdgePythonSVMClassifier< VImageDimension >::SetupClassifier()
{
    std::cout << "GraphEdgePythonSVMClassifier::SetupClassifier started" << std::endl;
#ifdef CS_BUILD_PYTHON
    PyObject *pName, *pSysPath, *pPath, *pModule, *pDict, *pArglist, *pTrainDBHandlerClass, *pSVMHandlingClass, *pValue1, *pValue2, *pAttributeName1, *pAttributeName2;

    pSysPath = PySys_GetObject("path");
    if(pSysPath == NULL)
        std::cout << "pSysPath still NULL" << std::endl;

    pPath = PyString_FromString(".");
    if(PyList_Append(pSysPath, pPath) < 0)
        std::cout << "PYList_Append < 0" << std::endl;

    pPath = PyString_FromString("./source/images/filters/classifiers/pythonScripts");
    if(PyList_Append(pSysPath, pPath) < 0)
        std::cout << "PYList_Append < 0" << std::endl;

    pName = PyString_FromString("testReadDataSetupSVMScript");      //name of file

    pModule = PyImport_Import(pName);
    if(pModule==NULL)
        std::cout << "pModule still NULL" << std::endl;

    pDict = PyModule_GetDict(pModule);

    pTrainDBHandlerClass = PyDict_GetItemString(pDict, "TrainingDataHandling");
    pArglist = Py_BuildValue("(s)", this->mTrainingDatabaseFilename.c_str());
    if(PyCallable_Check(pTrainDBHandlerClass))
        mpTrainDBHandlerInstance = PyObject_CallObject(pTrainDBHandlerClass, pArglist);
    else
        std::cout << "TrainingDataHandling can not be called" << std::endl;


    unsigned int sizeFeatureVector = 2*this->mNumActiveUnaryFeatures + this->mNumActiveBinaryFeatures;

    PyObject* list = PyList_New(sizeFeatureVector);

    int count=0;
    for(int j=0; j<LabelMapGraphType::NumberUnaryFeatures; ++j) {
        if(this->mLabelMapGraph->GetUnaryFeatureState(static_cast<UnaryFeatureType>(j)) == LabelMapGraphType::IsActive) {
            std::string f1 = LabelMapGraphType::UnaryFeatureName[j]+".1";
            std::string f2 = LabelMapGraphType::UnaryFeatureName[j]+".2";
            PyList_SetItem( list, count, PyString_FromString(f1.c_str()) );
            PyList_SetItem( list, count+1, PyString_FromString(f2.c_str()) );
            count+=2;
        }
    }
    for(int j=0; j<LabelMapGraphType::NumberBinaryFeatures; ++j) {
        if(this->mLabelMapGraph->GetBinaryFeatureState(static_cast<BinaryFeatureType>(j)) == LabelMapGraphType::IsActive) {
            std::string f = LabelMapGraphType::BinaryFeatureName[j];
            PyList_SetItem( list, count, PyString_FromString(f.c_str()) );
            count++;
        }
    }

    PyObject_CallMethod(mpTrainDBHandlerInstance, "loadDataset3", "(O)", list);
    PyObject_CallMethod(mpTrainDBHandlerInstance, "printSummary", "(s)", "py.script: before preprocessing: ");
    PyObject_CallMethod(mpTrainDBHandlerInstance, "preprocessData", NULL);
    PyObject_CallMethod(mpTrainDBHandlerInstance, "printSummary", "(s)", "py.script: after preprocessing: ");


    pAttributeName1 = PyString_FromString("trainingData");
    pValue1 = PyObject_GetAttr(mpTrainDBHandlerInstance, pAttributeName1);
    pAttributeName2 = PyString_FromString("scaler");
    pValue2 = PyObject_GetAttr(mpTrainDBHandlerInstance, pAttributeName2);
    pArglist = Py_BuildValue("(OO)", pValue1, pValue2);


    pSVMHandlingClass = PyDict_GetItemString(pDict, "SVMHandling");
    if(PyCallable_Check(pSVMHandlingClass))
        mpSVMHandlingInstance = PyObject_CallObject(pSVMHandlingClass, pArglist);
    else
        std::cout << "SVMHandling can not be called" << std::endl;

    PyObject_CallMethod(mpSVMHandlingInstance, "fitSVMToData", NULL);
#endif
    std::cout << "GraphEdgePythonSVMClassifier::SetupClassifier finished" << std::endl;
}


template< unsigned int VImageDimension > std::vector< std::pair<vtkIdType,vtkIdType> > GraphEdgePythonSVMClassifier< VImageDimension >::EvaluateAllLabelObjectPairs()
{
    vtkSmartPointer<vtkUnsignedLongArray> pedigreeIdArray = vtkUnsignedLongArray::SafeDownCast(this->mGraph->GetVertexData()->GetPedigreeIds());
    std::cout << "pedigreeIdArray->GetSize() = " << pedigreeIdArray->GetSize() << std::endl;

    this->mLabelObjectPairsOfTargetClass.clear();

#ifdef CS_BUILD_PYTHON
    PyObject *pList, *pValue;

    unsigned int sizeFeatureVector = 2*this->mNumActiveUnaryFeatureComponents + this->mNumActiveBinaryFeatures;

    std::cout << "GraphEdgePythonSVMClassifier::EvaluateAllLabelObjectPairs here" << std::endl;

    vtkSmartPointer<vtkEdgeListIterator> it = vtkSmartPointer<vtkEdgeListIterator>::New();
    this->mGraph->GetEdges(it);
    while (it->HasNext()) {
        vtkEdgeType e = it->Next();

        PyObject* list = PyList_New(sizeFeatureVector);

        int count=0;
        for(int j=0; j<LabelMapGraphType::NumberUnaryFeatures; ++j) {
            if(this->mLabelMapGraph->GetUnaryFeatureState(static_cast<UnaryFeatureType>(j)) == LabelMapGraphType::IsActive) {
                for(unsigned int k=0; k<LabelMapGraphType::UnaryFeatureComponents[j]; k++) {
                    std::stringstream featureName;
                    featureName << LabelMapGraphType::UnaryFeatureName[j];
                    if(LabelMapGraphType::UnaryFeatureComponents[j] > 1)
                        featureName << k;

                    PyList_SetItem( list, count, Py_BuildValue("f", vtkFloatArray::SafeDownCast(this->mGraph->GetVertexData()->GetArray(featureName.str().c_str()))->GetValue(e.Source) ));
                    PyList_SetItem( list, count+1, Py_BuildValue("f", vtkFloatArray::SafeDownCast(this->mGraph->GetVertexData()->GetArray(featureName.str().c_str()))->GetValue(e.Target) ));
                    count+=2;
                }
            }
        }
        for(int j=0; j<LabelMapGraphType::NumberBinaryFeatures; ++j) {
            if(this->mLabelMapGraph->GetBinaryFeatureState(static_cast<BinaryFeatureType>(j)) == LabelMapGraphType::IsActive) {
                PyList_SetItem( list, count, Py_BuildValue("f", vtkFloatArray::SafeDownCast(this->mGraph->GetEdgeData()->GetArray(LabelMapGraphType::BinaryFeatureName[j].c_str()))->GetValue(e.Id) ));
                count++;
            }
        }

        pValue = PyObject_CallMethod(mpSVMHandlingInstance, "predictFromMeasurement", "(O)", list);
        ClassType predClass = (ClassType)PyInt_AsLong(pValue);

        if(predClass == this->mTargetClass)
            this->mLabelObjectPairsOfTargetClass.push_back(std::pair<vtkIdType,vtkIdType>(pedigreeIdArray->GetValue(e.Source), pedigreeIdArray->GetValue(e.Target)));
    }
    std::cout << "GraphEdgePythonSVMClassifier::EvaluateAllLabelObjectPairs finished" << std::endl;
#endif

    return this->mLabelObjectPairsOfTargetClass;
}


template< unsigned int VImageDimension > std::vector< std::pair<vtkIdType,vtkIdType> > GraphEdgePythonSVMClassifier< VImageDimension >::EvaluateModifiedLabelObjectPairs(std::set<vtkIdType> &modifiedInLastIteration)
{
    vtkSmartPointer<vtkUnsignedLongArray> pedigreeIdArray = vtkUnsignedLongArray::SafeDownCast(this->mGraph->GetVertexData()->GetPedigreeIds());
    std::cout << "pedigreeIdArray->GetSize() = " << pedigreeIdArray->GetSize() << std::endl;

    std::set<vtkIdType> edgesAlreadyVisited;
    this->mLabelObjectPairsOfTargetClass.clear();

#ifdef CS_BUILD_PYTHON
    PyObject *pList, *pValue;

    unsigned int sizeFeatureVector = 2*this->mNumActiveUnaryFeatureComponents + this->mNumActiveBinaryFeatures;

    std::cout << "GraphEdgePythonSVMClassifier::EvaluateModifiedLabelObjectPairs here" << std::endl;

    for(std::set<vtkIdType>::iterator it=modifiedInLastIteration.begin(); it!=modifiedInLastIteration.end(); ++it) {
        vtkIdType source = this->mGraph->FindVertex(*it);

        if(source != -1) {
            vtkSmartPointer<vtkOutEdgeIterator> outIt = vtkSmartPointer<vtkOutEdgeIterator>::New();
            this->mGraph->GetOutEdges(source, outIt);

            while (outIt->HasNext()) {
                vtkOutEdgeType e = outIt->Next();

                if(edgesAlreadyVisited.count(e.Id)==0) {
                    vtkIdType target = e.Target;

                    PyObject* list = PyList_New(sizeFeatureVector);

                    int count=0;
                    for(int j=0; j<LabelMapGraphType::NumberUnaryFeatures; ++j) {
                        if(this->mLabelMapGraph->GetUnaryFeatureState(static_cast<UnaryFeatureType>(j)) == LabelMapGraphType::IsActive) {
                            for(unsigned int k=0; k<LabelMapGraphType::UnaryFeatureComponents[j]; k++) {
                                std::stringstream featureName;
                                featureName << LabelMapGraphType::UnaryFeatureName[j];
                                if(LabelMapGraphType::UnaryFeatureComponents[j] > 1)
                                    featureName << k;

                                PyList_SetItem( list, count, Py_BuildValue("f", vtkFloatArray::SafeDownCast(this->mGraph->GetVertexData()->GetArray(featureName.str().c_str()))->GetValue(source) ));
                                PyList_SetItem( list, count+1, Py_BuildValue("f", vtkFloatArray::SafeDownCast(this->mGraph->GetVertexData()->GetArray(featureName.str().c_str()))->GetValue(target) ));
                                count+=2;
                            }
                        }
                    }
                    for(int j=0; j<LabelMapGraphType::NumberBinaryFeatures; ++j) {
                        if(this->mLabelMapGraph->GetBinaryFeatureState(static_cast<BinaryFeatureType>(j)) == LabelMapGraphType::IsActive) {
                            PyList_SetItem( list, count, Py_BuildValue("f", vtkFloatArray::SafeDownCast(this->mGraph->GetEdgeData()->GetArray(LabelMapGraphType::BinaryFeatureName[j].c_str()))->GetValue(e.Id) ));
                            count++;
                        }
                    }

                    pValue = PyObject_CallMethod(mpSVMHandlingInstance, "predictFromMeasurement", "(O)", list);
                    ClassType predClass = (ClassType)PyInt_AsLong(pValue);

                    if(predClass == this->mTargetClass)
                        this->mLabelObjectPairsOfTargetClass.push_back(std::pair<vtkIdType,vtkIdType>(pedigreeIdArray->GetValue(source), pedigreeIdArray->GetValue(target)));

                    edgesAlreadyVisited.insert(e.Id);
                }
            }
        }
    }
#endif

    return this->mLabelObjectPairsOfTargetClass;
}
