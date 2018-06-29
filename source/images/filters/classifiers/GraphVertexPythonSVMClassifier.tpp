///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  GraphVertexPythonSVMClassifier.tpp                                   //
//                                                                                   //
//     Author:  Adrian Friebel <adrian.friebel@uni-leipzig.de>                       //
//    Created:  2015-03-23                                                           //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

#include "GraphVertexPythonSVMClassifier.h"

#include <vtkDataSetAttributes.h>
#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>
#include <vtkIntArray.h>
#include <vtkOutEdgeIterator.h>
#include <vtkUnsignedLongArray.h>
#include <vtkVertexListIterator.h>



template< unsigned int VImageDimension > GraphVertexPythonSVMClassifier< VImageDimension >::GraphVertexPythonSVMClassifier()
            : GraphVertexBaseClassifier< VImageDimension >()
{
    // TODO Auto-generated constructor stub

}


template< unsigned int VImageDimension > GraphVertexPythonSVMClassifier< VImageDimension >::~GraphVertexPythonSVMClassifier()
{
//#ifdef CS_BUILD_PYTHON
//    Py_XDECREF(mpTrainDBHandlerInstance);
//    Py_XDECREF(mpSVMHandlingInstance);
//    Py_Finalize();
//#endif
}


template< unsigned int VImageDimension > void GraphVertexPythonSVMClassifier< VImageDimension >::Init()
{
    std::cout << "GraphVertexPythonSVMClassifier::Init here" << std::endl;
    GraphVertexBaseClassifier<VImageDimension>::Init();
    std::cout << "GraphVertexBaseClassifier::Init finished" << std::endl;
    SetupClassifier();
}


template< unsigned int VImageDimension > void GraphVertexPythonSVMClassifier< VImageDimension >::SetupClassifier()
{
    std::cout << "GraphVertexPythonSVMClassifier::SetupClassifier started" << std::endl;
#ifdef CS_BUILD_PYTHON
    PyObject *pName, *pSysPath, *pPath, *pModule, *pDict, *pArglist, *pTrainDBHandlerClass, *pSVMHandlingClass, *pValue1, *pValue2, *pAttributeName1, *pAttributeName2, *pList;

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


    unsigned int sizeFeatureVector = this->mNumActiveUnaryFeatures;

    PyObject* list = PyList_New(sizeFeatureVector);

    int count=0;
    for(int j=0; j<LabelMapGraphType::NumberUnaryFeatures; ++j) {
        if(this->mLabelMapGraph->GetUnaryFeatureState(static_cast<UnaryFeatureType>(j)) == LabelMapGraphType::IsActive) {
            std::string f = LabelMapGraphType::UnaryFeatureName[j];
            PyList_SetItem( list, count, PyString_FromString(f.c_str()) );
            count++;
        }
    }

    PyObject_CallMethod(mpTrainDBHandlerInstance, "loadUnaryTrainingFeatures", "(O)", list);
    PyObject_CallMethod(mpTrainDBHandlerInstance, "printSummary", "(s)", "py.script: before preprocessing: ");
    PyObject_CallMethod(mpTrainDBHandlerInstance, "preprocessData", NULL);
    PyObject_CallMethod(mpTrainDBHandlerInstance, "printSummary", "(s)", "py.script: after preprocessing: ");

    pList = PyObject_CallMethod(mpTrainDBHandlerInstance, "getNamesOfTargetClasses", NULL);
    std::map<vtkIdType, float> probabilities;
    for(unsigned int i=0; i<PyList_Size(pList); i++)
        this->mClassProbabilites.insert( std::pair<std::string, std::map<vtkIdType, float> >("ClassProbability_" + std::string(PyString_AsString(PyList_GetItem(pList, i))), probabilities) );

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
    std::cout << "GraphVertexPythonSVMClassifier::SetupClassifier finished" << std::endl;
}


template< unsigned int VImageDimension > std::set<vtkIdType> GraphVertexPythonSVMClassifier< VImageDimension >::EvaluateAllLabelObjects()
{
    vtkSmartPointer<vtkUnsignedLongArray> pedigreeIdArray = vtkUnsignedLongArray::SafeDownCast(this->mGraph->GetVertexData()->GetPedigreeIds());
    std::cout << "pedigreeIdArray->GetSize() = " << pedigreeIdArray->GetSize() << std::endl;

    this->mLabelObjectsOfTargetClass.clear();

#ifdef CS_BUILD_PYTHON
    PyObject *pList, *pValue;

    unsigned int sizeFeatureVector = this->mNumActiveUnaryFeatureComponents;

    std::cout << "GraphVertexPythonSVMClassifier::EvaluateAllLabelObjects here" << std::endl;

    vtkSmartPointer<vtkVertexListIterator> it = vtkSmartPointer<vtkVertexListIterator>::New();
    this->mGraph->GetVertices(it);
    while (it->HasNext()) {
        vtkIdType v = it->Next();

        PyObject* list = PyList_New(sizeFeatureVector);

        int count=0;
        for(int j=0; j<LabelMapGraphType::NumberUnaryFeatures; ++j) {
            if(this->mLabelMapGraph->GetUnaryFeatureState(static_cast<UnaryFeatureType>(j)) == LabelMapGraphType::IsActive) {
                for(int k=0; k<LabelMapGraphType::UnaryFeatureComponents[j]; k++) {
                    std::stringstream featureName;
                    featureName << LabelMapGraphType::UnaryFeatureName[j];
                    if(LabelMapGraphType::UnaryFeatureComponents[j] > 1)
                        featureName << k;

                    PyList_SetItem( list, count, Py_BuildValue("f", vtkFloatArray::SafeDownCast(this->mGraph->GetVertexData()->GetArray(featureName.str().c_str()))->GetValue(v) ));
                    count++;
                }
            }
        }
        pList = PyObject_CallMethod(mpSVMHandlingInstance, "predictProbaFromMeasurement", "(O)", list);
        std::map<std::string, std::map<vtkIdType, float> >::iterator it = this->mClassProbabilites.begin();
        for(unsigned int i=0; i<PyList_Size(pList); i++) {
            if(it!=this->mClassProbabilites.end()) {
                float prob = PyFloat_AsDouble(PyList_GetItem(pList, i));
                it->second.insert( std::pair<vtkIdType, float>(v, prob) );
                ++it;
            }
        }

        pValue = PyObject_CallMethod(mpSVMHandlingInstance, "predictFromMeasurement", "(O)", list);
        ClassType predClass = (ClassType)PyInt_AsLong(pValue);

        if(predClass == this->mTargetClass)
            this->mLabelObjectsOfTargetClass.insert(pedigreeIdArray->GetValue(v));
    }
    std::cout << "GraphVertexPythonSVMClassifier::EvaluateAllLabelObjects finished" << std::endl;
#endif

    return this->mLabelObjectsOfTargetClass;
}


template< unsigned int VImageDimension > std::set<vtkIdType> GraphVertexPythonSVMClassifier< VImageDimension >::EvaluateModifiedLabelObjects(std::set<vtkIdType> &modifiedInLastIteration)
{
//    vtkSmartPointer<vtkUnsignedLongArray> pedigreeIdArray = vtkUnsignedLongArray::SafeDownCast(this->mGraph->GetVertexData()->GetPedigreeIds());
//    std::cout << "pedigreeIdArray->GetSize() = " << pedigreeIdArray->GetSize() << std::endl;
//
//    std::set<vtkIdType> edgesAlreadyVisited;
//    this->mLabelObjectsOfTargetClass.clear();
//
//#ifdef CS_BUILD_PYTHON
//    PyObject *pList, *pValue;
//
//    unsigned int sizeFeatureVector = this->mNumActiveUnaryFeatureComponents;
//
//    std::cout << "GraphVertexPythonSVMClassifier::EvaluateModifiedLabelObjects here" << std::endl;
//
//    for(std::set<vtkIdType>::iterator it=modifiedInLastIteration.begin(); it!=modifiedInLastIteration.end(); ++it) {
//        vtkIdType source = this->mGraph->FindVertex(*it);
//
//        if(source != -1) {
//            vtkSmartPointer<vtkOutEdgeIterator> outIt = vtkSmartPointer<vtkOutEdgeIterator>::New();
//            this->mGraph->GetOutEdges(source, outIt);
//
//            while (outIt->HasNext()) {
//                vtkOutEdgeType e = outIt->Next();
//
//                if(edgesAlreadyVisited.count(e.Id)==0) {
//                    vtkIdType target = e.Target;
//
//                    PyObject* list = PyList_New(sizeFeatureVector);
//
//                    int count=0;
//                    for(int j=0; j<LabelMapGraphType::NumberUnaryFeatures; ++j) {
//                        if(this->mLabelMapGraph->GetUnaryFeatureState(static_cast<UnaryFeatureType>(j)) == LabelMapGraphType::IsActive) {
//                            for(unsigned int k=0; k<LabelMapGraphType::UnaryFeatureComponents[j]; k++) {
//                                std::stringstream featureName;
//                                featureName << LabelMapGraphType::UnaryFeatureName[j];
//                                if(LabelMapGraphType::UnaryFeatureComponents[j] > 1)
//                                    featureName << k;
//
//                                PyList_SetItem( list, count, Py_BuildValue("f", vtkFloatArray::SafeDownCast(this->mGraph->GetVertexData()->GetArray(featureName.str().c_str()))->GetValue(source) ));
//                                PyList_SetItem( list, count+1, Py_BuildValue("f", vtkFloatArray::SafeDownCast(this->mGraph->GetVertexData()->GetArray(featureName.str().c_str()))->GetValue(target) ));
//                                count+=2;
//                            }
//                        }
//                    }
//
//                    pValue = PyObject_CallMethod(mpSVMHandlingInstance, "predictFromMeasurement", "(O)", list);
//                    ClassType predClass = (ClassType)PyInt_AsLong(pValue);
//
//                    if(predClass == this->mTargetClass)
//                        this->mLabelObjectsOfTargetClass.push_back(std::pair<vtkIdType,vtkIdType>(pedigreeIdArray->GetValue(source), pedigreeIdArray->GetValue(target)));
//
//                    edgesAlreadyVisited.insert(e.Id);
//                }
//            }
//        }
//    }
//#endif
//
//    return this->mLabelObjectsOfTargetClass;
}
