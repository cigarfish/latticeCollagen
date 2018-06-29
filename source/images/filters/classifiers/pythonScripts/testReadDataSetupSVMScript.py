'''testReadDataSetupSVMScript.py - Python source designed to test basic functionality'''

import os
import csv
import shutil
import warnings
from os import environ
from os.path import dirname
from os.path import join
from os.path import exists
from os.path import expanduser
from os.path import isdir
from os import listdir
from os import makedirs

import numpy as np

import sklearn
from sklearn import preprocessing
from sklearn import svm


class Bunch(dict):
    """Container object for datasets: dictionary-like object that
       exposes its keys as attributes."""

    def __init__(self, **kwargs):
        dict.__init__(self, kwargs)
        self.__dict__ = self


class TrainingDataHandling:
    
    def __init__(self, pathToDB):
        print "py.script: init TrainingDataHandling"
        self.pathToTrainingDB = pathToDB 
        self.trainingData = Bunch()
        self.scaler = sklearn.preprocessing.StandardScaler()
    
    #DEPRECATED: first version of training dataset load method
    def loadDataset(self):
        print "py.script: loadDataset from path = {}".format(self.pathToTrainingDB)
        with open(self.pathToTrainingDB) as csv_file:
            data_file = csv.reader(csv_file)
            temp = next(data_file)
            n_samples = int(temp[1])
            n_features = int(temp[2])
            n_classes = int(temp[3])
            feature_names = np.array(temp[4:(4+n_features)])
            target_names = np.array(temp[(4+n_features):((4+n_features+n_classes))])
            data = np.empty((n_samples, n_features))
            target = np.empty((n_samples,), dtype=np.int)

            for i, ir in enumerate(data_file):                            
                data[i] = np.asarray(ir[0:n_features], dtype=np.float)
                target[i] = np.asarray(ir[n_features], dtype=np.int)

            self.trainingData = Bunch(data=data, target=target, target_names=target_names, feature_names=feature_names)
            
    #DEPRECATED: load method abel of per-default single-component feature training dataset handling
    def loadDataset2(self, selectedFeatures):
        print "py.script: loadDataset2 from path = {}".format(self.pathToTrainingDB)
        selectedFeatures = [feature.strip() for feature in selectedFeatures]
        
        with open(self.pathToTrainingDB) as csv_file:
            data_file = csv.reader(csv_file)
            temp = next(data_file)
            n_samples = int(temp[1])
            n_features = int(temp[2])
            n_classes = int(temp[3])
            feature_names = np.array(temp[4:(4+n_features)])
            feature_names = [feature.strip() for feature in feature_names]
            
            featureFilter = np.empty(n_features, dtype=bool)
            for i, ir in enumerate(feature_names):
                featureFilter[i] = self.find_helper(selectedFeatures, ir)
                        
            target_names = np.array(temp[(4+n_features):((4+n_features+n_classes))])
            data = np.empty((n_samples, n_features))
            target = np.empty((n_samples,), dtype=np.int)

            for i, ir in enumerate(data_file):                            
                data[i] = np.asarray(ir[0:n_features], dtype=np.float)
                target[i] = np.asarray(ir[n_features], dtype=np.int)
                
            selectedData = data[:,featureFilter]
            selectedDataNP = np.array(selectedData)
            selectedFeaturesNP = np.array(selectedFeatures)
            self.trainingData = Bunch(data=selectedDataNP, target=target, target_names=target_names, feature_names=selectedFeaturesNP)
            
    #DEPRECATED: load method able of multi-component feature training dataset handling
    def loadDataset3(self, selectedFeatures):
        print "py.script: loadDataset3 path = {}".format(self.pathToTrainingDB)
        print "py.script: loadDataset3 selectedFeatures = {}".format(selectedFeatures)
        selectedFeatures = [feature.strip() for feature in selectedFeatures]
        
        with open(self.pathToTrainingDB) as csv_file:
            data_file = csv.reader(csv_file)
            temp = next(data_file)
            n_samples = int(temp[1])
            n_unaryFeatures = int(temp[2])
            n_unaryFeatureComponents = int(temp[3])
            n_binaryFeatures = int(temp[4])
            n_classes = int(temp[5])
            feature_names = np.array(temp[6 : (6+2*(n_unaryFeatures+n_binaryFeatures)) : 2])
            feature_names = [feature.strip() for feature in feature_names]
            feature_components = np.array(temp[7 : (7+2*(n_unaryFeatures+n_binaryFeatures)) : 2], dtype=np.int)
            
            print "feature_names = {}".format(feature_names)
            print "feature_components = {}".format(feature_components)

            feature_component_names = np.empty(0)
            for(i, ir) in enumerate(feature_names):
                for j in range(feature_components[i]):
                    feature_component_names = np.append(feature_component_names, ir)
            
            featureFilter = np.empty(n_unaryFeatureComponents+n_binaryFeatures, dtype=bool)
            for i, ir in enumerate(feature_component_names):
                featureFilter[i] = self.find_helper(selectedFeatures, ir)
                
            target_names = np.array(temp[(6+2*(n_unaryFeatures+n_binaryFeatures)):((6+2*(n_unaryFeatures+n_binaryFeatures)+n_classes))])
            
            print "target_names = {}".format(target_names)
            
            data = np.empty((n_samples, n_unaryFeatureComponents+n_binaryFeatures))
            target = np.empty((n_samples,), dtype=np.int)

            for i, ir in enumerate(data_file):                            
                data[i] = np.asarray(ir[0:(n_unaryFeatureComponents+n_binaryFeatures)], dtype=np.float)
                target[i] = np.asarray(ir[(n_unaryFeatureComponents+n_binaryFeatures)], dtype=np.int)
                
            selectedData = data[:,featureFilter]
            selectedDataNP = np.array(selectedData)
            selectedFeaturesNP = np.array(selectedFeatures)
            self.trainingData = Bunch(data=selectedDataNP, target=target, target_names=target_names, feature_names=selectedFeaturesNP)
            
    def loadUnaryTrainingFeatures(self, selectedFeatures):
        print "py.script: loadDataset3 path = {}".format(self.pathToTrainingDB)
        print "py.script: loadDataset3 selectedFeatures = {}".format(selectedFeatures)
        selectedFeatures = [feature.strip() for feature in selectedFeatures]
        
        with open(self.pathToTrainingDB) as csv_file:
            data_file = csv.reader(csv_file)
            temp = next(data_file)
            n_unarySamples = int(temp[2])
            n_binarySamples = int(temp[4])
            n_unaryFeatures = int(temp[6])
            n_unaryFeatureComponents = int(temp[8])
            n_binaryFeatures = int(temp[10])
            n_classes = int(temp[12])
            
            unaryFeatureNamesStart = 14
            unaryFeatureNamesEnd = unaryFeatureNamesStart+2*n_unaryFeatures
            unaryFeatureCompStart = unaryFeatureNamesStart+1
            unaryFeatureCompEnd = unaryFeatureCompStart+2*n_unaryFeatures
            binaryFeatureNamesStart = unaryFeatureNamesEnd+1
            binaryFeatureNamesEnd = binaryFeatureNamesStart + 2*n_binaryFeatures
            binaryFeatureCompStart = binaryFeatureNamesStart+1
            binaryFeatureCompEnd = binaryFeatureCompStart + 2*n_binaryFeatures
            classStart = binaryFeatureNamesEnd+1
            classEnd = classStart + n_classes
            
            
            unaryFeatureNames = np.array(temp[unaryFeatureNamesStart : unaryFeatureNamesEnd : 2])
            unaryFeatureNames = [feature.strip() for feature in unaryFeatureNames]
            unaryFeatureComponents = np.array(temp[unaryFeatureCompStart : unaryFeatureCompEnd : 2], dtype=np.int)
            binaryFeatureNames = np.array(temp[binaryFeatureNamesStart : binaryFeatureNamesEnd : 2])
            binaryFeatureNames = [feature.strip() for feature in binaryFeatureNames]
            binaryFeatureComponents = np.array(temp[binaryFeatureCompStart : binaryFeatureCompEnd : 2], dtype=np.int)
            
            print "unaryFeatureNames = {}".format(unaryFeatureNames)
            print "unaryFeatureComponents = {}".format(unaryFeatureComponents)
            print "binaryFeatureNames = {}".format(binaryFeatureNames)
            print "binaryFeatureComponents = {}".format(binaryFeatureComponents)

            feature_component_names = np.empty(0)
            for(i, ir) in enumerate(unaryFeatureNames):
                for j in range(unaryFeatureComponents[i]):
                    feature_component_names = np.append(feature_component_names, ir)
            
            featureFilter = np.empty(n_unaryFeatureComponents, dtype=bool)
            for i, ir in enumerate(feature_component_names):
                featureFilter[i] = self.find_helper(selectedFeatures, ir)
                
            target_names = np.array(temp[classStart : classEnd])
            for i, ir in enumerate(target_names):
                target_names[i] = ir.strip()
            print "target_names = {}".format(target_names)
            
            data = np.empty((n_unarySamples, n_unaryFeatureComponents))
            target = np.empty((n_unarySamples,), dtype=np.int)

            for i, ir in enumerate(data_file):                   
                if(i<n_unarySamples):
                    data[i] = np.asarray(ir[1:(1+n_unaryFeatureComponents)], dtype=np.float)
                    target[i] = np.asarray(ir[(1+n_unaryFeatureComponents)], dtype=np.int)
                
            selectedData = data[:,featureFilter]
            selectedDataNP = np.array(selectedData)
            selectedFeaturesNP = np.array(selectedFeatures)
            self.trainingData = Bunch(data=selectedDataNP, target=target, target_names=target_names, feature_names=selectedFeaturesNP)

    def find_helper(self, list, element):
        for i, x in enumerate(list):
            if x==element:
                return True
        return False
    
    def preprocessData(self):
        print "py.script: preprocessData"
#        self.trainingData.data = sklearn.preprocessing.scale(self.trainingData.data)
        self.scaler = sklearn.preprocessing.StandardScaler().fit(self.trainingData.data)
        self.trainingData.data = self.scaler.transform(self.trainingData.data)
        
    def getNamesOfTargetClasses(self):
        return self.trainingData.target_names[np.unique( self.trainingData.target )].tolist()
        
    def printSummary(self, prefix):
        print prefix + " data {}".format(self.trainingData.data)
        print prefix + " targets {}".format(self.trainingData.target)
        print prefix + " target names {}".format(self.trainingData.target_names)
        print prefix + " target names used {}".format(self.getNamesOfTargetClasses())
        print prefix + " feature names {}".format(self.trainingData.feature_names)
        

class SVMHandling:
    
    def __init__(self, bunch, scaler):
        self.trainingData = Bunch()
        self.trainingData = bunch
        self.scaler = sklearn.preprocessing.StandardScaler()
        self.scaler = scaler
        self.svm = sklearn.svm.SVC(gamma=0.001, C=100., probability=True)
        print "py.script: SVMHandling here: {}".format(self.trainingData)
        
    def fitSVMToData(self):
        self.svm.fit(self.trainingData.data[:-1], self.trainingData.target[:-1])
        print "py.script: SVM fitted"
        
    def predictFromMeasurement(self, measurement):
        #print "py.script: measurement = "
        #print measurement
        npm = np.array(measurement)
        npm = self.scaler.transform(npm)
        #print "py.script: measurement scaled = "
        #print npm
        return self.svm.predict(npm)
    
    def predictProbaFromMeasurement(self, measurement):
        #print "py.script: measurement = "
        #print measurement
        npm = np.array(measurement)
        npm = self.scaler.transform(npm)
        #print "py.script: measurement scaled = "
        #print npm
        return self.svm.predict_proba(npm)[0].tolist()
    
    def printArray(self, a):
        print "py.script: Contents of a :"
        print a
        npm = np.array(a)
        print "py.script: As array:"
        print npm
        

                 