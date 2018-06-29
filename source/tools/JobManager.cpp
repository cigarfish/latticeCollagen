///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  JobManager.cpp                                                       //
//                                                                                   //
//     Author:  Adrian Friebel <adrian.friebel@uni-leipzig.de>                       //
//    Created:  2012-12-20                                                           //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


#include "JobManager.h"

#include <fstream>
#include <string>
#include <sstream>

#include "../images/filters/analysisFilters/AnalyzeCellsFilter.h"
#include "../images/filters/analysisFilters/AnalyzeLobuleFilter.h"
#include "../images/filters/analysisFilters/AnalyzeNucleiFilter.h"
#include "../../images/filters/graphFilters/AnalyzeBileNetworkFilter.h"
#include "../images/pipelines/CropDataSet.h"
#include "../images/pipelines/EstimateCellShape.h"
#include "../images/pipelines/EstimateLobuleShape.h"
#include "../images/pipelines/ExtractGraph.h"
#include "../images/pipelines/SegmentBileNetworkOnDPPIV20x.h"
#include "../images/pipelines/SegmentBileNetworkOnDPPIV60x.h"
#include "../images/pipelines/SegmentCellMembraneOnBCat60x.h"
#include "../images/pipelines/SegmentNecroticRegionOnDM.h"
#include "../images/pipelines/SegmentNucleiOnDAPI20x.h"
#include "../images/pipelines/SegmentNucleiOnDAPI60x.h"
#include "../images/pipelines/SegmentNucleiOnDAPIWithHough.h"
#include "../images/pipelines/SegmentAndClassifyStellateCells.h"
#include "../images/pipelines/SegmentSinusoidalNetworkOnTwoChannels20x.h"
#include "../images/pipelines/SegmentSinusoidalNetworkOnTwoChannels60x.h"
#include "../images/pipelines/ObjectBasedSegmentation.h"
#include "../images/pipelines/Skeletonization.h"
#include "../images/pipelines/SuperpixelPreparation.h"
#include "../images/pipelines/PreprocessDataSet.h"
#include "../images/tools/FileFormatConverter.h"

#include "../tools/input/FilenameParser.h"

#include "parameters/CSParameterContext.h"
#include "parameters/CSParameterContextTemporary.h"
#include "parameters/CSParameterTemporary.h"

#include <QtCore>

#include <QTableView>
#include "../../gui/tabImageProcessing/JobManagerModel.h"


std::string JobManager::JobTypeString[] = {
    "CLAHEJob",
    "BackgroundEliminationJob",
    "HConvexImageFilterJob",
    "CropDataSetJob",
    "SegmentNecroticRegionJob",
    "SegmentNuclei20xJob",
    "SegmentNuclei60xJob",
    "SegmentNucleiWithHoughJob",
    "EstimateCellShapeJob",
    "EstimateLobuleShapeJob",
    "AnalyzeCellsJob",
    "AnalyzeLobuleJob",
    "AnalyzeNucleiJob",
    "SegmentAndClassifyStellateCellsJob",
    "SegmentSinusoidsBile20xJob",
    "SegmentSinusoidsBile60xJob",
    "SuperpixelPreparationJob",
    "ObjectBasedSegmentationJob",
    "SkeletonizationJob",
    "ExtractAndAnalyzeGraphJob",
    "SegmentCellMembraneJob",
    "FileFormatConversionJob",
    ""
};


JobManager::JobManager(void)
    : qtModel( NULL )
{
#ifdef __BUILD_MAC__
    jobFilePath = QDir::home().absolutePath().toStdString();
    jobFilePath += "/Library/Application Support/TiQuant/JobManager/";
    QDir("/").mkpath( QString(jobFilePath.c_str()) );
#else
    jobFilePath = "./JobManager/";
    QDir().mkpath( QString(jobFilePath.c_str()) );
#endif
    jobQueue = new JobQueueType( "JobManager" );
}


JobManager::~JobManager(void)
{
    if ( qtModel )
        delete qtModel;
}


void JobManager::AddJob(JobType type, CSParameterContext * parameters)
{
    int id;
    std::stringstream idstr;

	if(!jobQueue->hasChildren()) id = 0;
	else id = atoi( jobQueue->getChildren().back()->name().c_str() ) + 1;

    idstr << id;

    std::string name = idstr.str();
    Job * newJob = new Job( name );

	newJob->addParameter("type", CSParameter::Int, new int(type), "");

    CSParameterContextTemporary * payLoad
        = new CSParameterContextTemporary( parameters );

    newJob->addContext( payLoad );

    emit aboutToChange();

	jobQueue->addContext(newJob);

    emit changed();

	RewriteToDoFile();
}


void JobManager::AddJobs(JobType type, CSParameterContext * parameters)
{
    std::vector<CSParameter *> params, paramsDatasetID;
    std::vector<std::string> fileNames;
    std::vector<std::string> datasetIDName;
    std::vector< std::vector<std::string> > concretePaths;
    bool setDatasetID = false;

    params = parameters->findParametersByDatatype(CSParameter::FileName);
    paramsDatasetID = parameters->findParametersBySubstring("Dataset ID");

    for(unsigned int i=0; i<params.size(); i++)
        std::cout << "parameters[" << i << "] = " << params[i]->name() << std::endl;

    bool isAllRegExp = true;
    bool isNoneRegExp = true;
    std::vector<bool> isRegExp;
    for(unsigned int i=0; i<params.size(); i++) {
        fileNames.push_back( *(std::string*)(parameters->findParameter(params[i]->name(), 0)->dataPointer()) );
        std::size_t found = fileNames[i].find("*");
        if(found!=std::string::npos)
            isRegExp.push_back(true);
        else
            isRegExp.push_back(false);
        isAllRegExp = isAllRegExp && isRegExp[i];
        isNoneRegExp = isNoneRegExp && !isRegExp[i];
    }
    std::string datasetID;
    if(paramsDatasetID.size()>0) {
        datasetID = *(std::string*)(parameters->findParameter(paramsDatasetID[0]->name(), 0)->dataPointer());
        std::size_t found = datasetID.find("*");
        if(found!=std::string::npos)
            setDatasetID = true;
    }

    QFileInfo infoFile;
    infoFile.setFile(QString::fromStdString(fileNames[0]));

    QString path = infoFile.path() + QString("/");
    QString filename = infoFile.baseName();
    QString ext = QString(".") + infoFile.suffix();

//    if(!isAllRegExp && !isNoneRegExp)
//        throw std::string("Please do not mix regExp and explicit filenames.");
    if(isNoneRegExp)
        this->AddJob(type, parameters);
    if(isAllRegExp) {
        if(!path.contains("*"))                                                         //infoFile.path() is a complete path
            ;
        if(path.contains("*") && !path.endsWith("*/"))
            ;//throw std::string("Please do not mix regExp and explicit filenames.");      //not handled e.g. /path/*/to/bla
        if(path.contains("*") && path.endsWith("*/")) {
            QString pathTemp = infoFile.path();
            pathTemp.chop(1);
            QDir rootDir,dir;
            rootDir.setPath(pathTemp);
            rootDir.setFilter(QDir::AllDirs | QDir::NoDotAndDotDot | QDir::NoSymLinks);

            QStringList dirList = rootDir.entryList();
            for(int i=0; i<dirList.size(); ++i) {
                std::vector<std::string> cPaths;
                bool allPathsOK = true;
                QString newPath = QString("%1/%2").arg(rootDir.absolutePath()).arg(dirList.at(i)) + QString("/");
                std::string datasetID = dirList.at(i).toStdString();
                std::cout << "Found dir: " << newPath.toStdString() << std::endl;

                dir.setPath(newPath);
                for(unsigned int k=0; k<params.size(); k++) {
                    QFileInfo infoFileTemp;
                    infoFileTemp.setFile(QString::fromStdString(fileNames[k]));
                    dir.setNameFilters(QStringList(infoFileTemp.baseName() + QString(".") + infoFileTemp.suffix()));
                    dir.setFilter(QDir::Files | QDir::NoDotAndDotDot | QDir::NoSymLinks);

                    std::cout << "Scanning: " << dir.path().toStdString() << std::endl;

                    QStringList fileList = dir.entryList();
                    for(int j=0; j<fileList.count(); j++) {
                        std::cout << "Found file: " << fileList[j].toStdString() << std::endl;

                        QFileInfo infoFileTemp2;
                        infoFileTemp2.setFile(dir.path() + QString("/") + fileList[j]);
                        if(infoFileTemp2.exists())
                            cPaths.push_back(dir.path().toStdString() + "/" + fileList[j].toStdString());
                        else
                            allPathsOK = false;
                    }
                    if(fileList.count()==0)
                        allPathsOK = false;
                }
                if(allPathsOK) {
                    concretePaths.push_back(cPaths);
                    datasetIDName.push_back(datasetID);
                }
            }
        }

        for(unsigned int i=0; i<concretePaths.size(); i++) {
            if(setDatasetID)
                parameters->findParameter(paramsDatasetID[0]->name(), 0)->setData(datasetIDName[i]);
            for(unsigned j=0; j<params.size(); j++)
                parameters->findParameter(params[j]->name(), 0)->setData(concretePaths[i][j]);
            this->AddJob(type, parameters);
        }
    }
}


void JobManager::DeleteJob( Job * job )
{
    int id = atoi( job->name().c_str() );

    if ( job )
    {
        emit aboutToChange();

        std::cout << "Found job with Id " << id << std::endl;
        jobQueue->removeContext( job );
        std::cout << "Deleting job with Id " << id << std::endl;

        delete job;

        emit changed();

        RewriteToDoFile();
    }
    else
        std::cout << "There is no job with id " << id << std::endl;
}


int JobManager::GetNumberOfQueuedJobs() const
{ return jobQueue->getChildren().size(); }


JobManager::Job * JobManager::GetJob( int id ) const
{
    std::stringstream idstr;
    idstr << id;
    return (Job *) jobQueue->findContext( idstr.str() );
}


void JobManager::InitWithJobQueueFile(std::string path)
{
    std::string jobFileName;
    std::string jobFileExt;
    bool f1 = FilenameParser::ParseFilename(path, jobFilePath, jobFileName, jobFileExt);


    QFile infile( QString((jobFilePath + jobFileName + jobFileExt).c_str()) );

    infile.open( QFile::ReadOnly | QFile::Text );

    QXmlStreamReader xmlReader( &infile );

    xmlReader.readNextStartElement();

    // ToDo:  Test if we really have opened a jobTodoFile

    CSParameterContextTemporary * ctxt =
        CSParameterContextTemporary::createFromXML( &xmlReader );

    infile.close();

    if ( ! ctxt )
        std::cout << "Got empty context\n";

    emit aboutToChange();

    std::vector<CSParameterContext *> jobs = ctxt->getChildren();
    std::vector<CSParameterContext *>::iterator jobIter = jobs.begin();
    while ( jobIter != jobs.end() )
    {
        ctxt->removeContext( *jobIter );
        jobQueue->addContext( *jobIter );
        ++jobIter;
    }

    emit changed();

    ctxt->getChildren().clear();
    delete ctxt;

}


void JobManager::Start()
{
    unsigned int ImageDimension = 3;

    try{
        while(jobQueue->hasChildren())
        {
            std::cout << "Oh no, there are still jobs to process!" << std::endl;

            Job * j = (Job *) jobQueue->getChildren()[0];
            JobType type = *( (JobType *) j->findParameter("type")->dataPointer() );

            CSParameter *dataDimParam = j->findParameter("Data dimensionality", 0);
            if(dataDimParam != NULL) {
                std::string dataDimMode = ( (CSParameterChoice*)(dataDimParam->dataPointer()) )->currentString();
                if( dataDimMode.compare("2D")==0 )
                    ImageDimension = 2;
                else if( dataDimMode.compare("3D")==0 )
                    ImageDimension = 3;
                else
                    throw std::string("Unknown data dimensionality.");
            }

            switch(type)
            {
            case CLAHEJob:
                std::cout << "Start PreprocessDataSet Pipeline!" << std::endl;
                startPreprocessDataSet( j, 0, ImageDimension );
                std::cout << "Finished PreprocessDataSet Pipeline!" << std::endl;
                break;
            case BackgroundEliminationJob:
                std::cout << "Start PreprocessDataSet Pipeline!" << std::endl;
                startPreprocessDataSet( j, 1, ImageDimension );
                std::cout << "Finished PreprocessDataSet Pipeline!" << std::endl;
                break;
            case HConvexImageFilterJob:
                std::cout << "Start PreprocessDataSet Pipeline!" << std::endl;
                startPreprocessDataSet( j, 2, ImageDimension );
                std::cout << "Finished PreprocessDataSet Pipeline!" << std::endl;
                break;
            case CropDataSetJob:
                std::cout << "Start CropDataSet Pipeline!" << std::endl;
                startCropDataSet( j, ImageDimension );
                std::cout << "Finished CropDataSet Pipeline!" << std::endl;
                break;
            case SegmentNecroticRegionJob:
                std::cout << "Start SegmentNecroticRegion Pipeline!" << std::endl;
                startSegmentNecroticRegion( j );
                std::cout << "Finished SegmentNecroticRegion Pipeline!" << std::endl;
                break;
            case SegmentNuclei20xJob:
                std::cout << "Start Nuclei20xSegmentation Pipeline!" << std::endl;
                startSegmentNuclei20x( j );
                std::cout << "Finished Nuclei20xSegmentation Pipeline!" << std::endl;
                break;
            case SegmentNuclei60xJob:
                std::cout << "Start Nuclei60xSegmentation Pipeline!" << std::endl;
                startSegmentNuclei60x( j, ImageDimension );
                std::cout << "Finished Nuclei60xSegmentation Pipeline!" << std::endl;
                break;
            case SegmentNucleiWithHoughJob:
                std::cout << "Start NucleiWithHoughSegmentation Pipeline!" << std::endl;
                startSegmentNucleiWithHough( j );
                std::cout << "Finished NucleiWithHoughSegmentation Pipeline!" << std::endl;
                break;
            case EstimateCellShapeJob:
                std::cout << "Start EstimateCellShape Pipeline!" << std::endl;
                startEstimateCellShape( j );
                std::cout << "Finished EstimateCellShape Pipeline!" << std::endl;
                break;
            case EstimateLobuleShapeJob:
                std::cout << "Start EstimateLobuleShape Pipeline!" << std::endl;
                startEstimateLobuleShape( j, ImageDimension );
                std::cout << "Finished EstimateLobuleShape Pipeline!" << std::endl;
                break;
            case AnalyzeCellsJob:
                std::cout << "Start AnalyzeCells Pipeline!" << std::endl;
                startAnalyzeCells( j );
                std::cout << "Finished AnalyzeCells Pipeline!" << std::endl;
                break;
            case AnalyzeLobuleJob:
                std::cout << "Start AnalyzeLobule Pipeline!" << std::endl;
                startAnalyzeLobule( j );
                std::cout << "Finished AnalyzeLobule Pipeline!" << std::endl;
                break;
            case AnalyzeNucleiJob:
                std::cout << "Start AnalyzeNuclei Pipeline!" << std::endl;
                startAnalyzeNuclei( j, ImageDimension );
                std::cout << "Finished AnalyzeNuclei Pipeline!" << std::endl;
                break;
            case SegmentAndClassifyStellateCellsJob:
                std::cout << "Start StellateCellSegmentation Pipeline!" << std::endl;
                startSegmentAndClassifyStellateCells( j );
                std::cout << "Finished StellateCellSegmentation Pipeline!" << std::endl;
                break;
            case SegmentSinusoidsBile20xJob:
                std::cout << "Start SinusoidBile20xSegmentation Pipeline!" << std::endl;
                startSinusoidBileSegmentation20x( j );
                std::cout << "Finished SinusoidBile20xSegmentation Pipeline!" << std::endl;
                break;
            case SegmentSinusoidsBile60xJob:
                std::cout << "Start SinusoidBile60xSegmentation Pipeline!" << std::endl;
                startSinusoidBileSegmentation60x( j );
                std::cout << "Finished SinusoidBile60xSegmentation Pipeline!" << std::endl;
                break;
            case SuperpixelPreparationJob:
                std::cout << "Start SuperpixelPreparation Pipeline!" << std::endl;
                startSuperpixelPreparation( j, ImageDimension );
                std::cout << "Finished SuperpixelPreparation Pipeline!" << std::endl;
                break;
            case ObjectBasedSegmentationJob:
                std::cout << "Start ObjectBasedSegmentation Pipeline!" << std::endl;
                startObjectBasedSegmentation( j, ImageDimension );
                std::cout << "Finished ObjectBasedSegmentation Pipeline!" << std::endl;
                break;
            case SkeletonizationJob:
                std::cout << "Start Skeletonization Pipeline!" << std::endl;
                startSkeletonization( j );
                std::cout << "Finished Skeletonization Pipeline!" << std::endl;
                break;
            case ExtractAndAnalyzeGraphJob:
                std::cout << "Start ExtractAndAnalyzeGraph Pipeline!" << std::endl;
                startExtractAndAnalyzeGraph( j );
                std::cout << "Finished ExtractAndAnalyzeGraph Pipeline!" << std::endl;
                break;
            case SegmentCellMembraneJob:
                std::cout << "Start SegmentCellMembrane Pipeline!" << std::endl;
                startSegmentCellMembrane( j );
                std::cout << "Finished SegmentCellMembrane Pipeline!" << std::endl;
                break;
            case FileFormatConversionJob:
                std::cout << "Start File Format Conversion Pipeline!" << std::endl;
                startFileFormatConversionJob( j );
                std::cout << "Finished File Format Conversion Pipeline!" << std::endl;
                break;

            default:
                std::cout << "Undefined job type!" << std::endl;
                break;
            }

            WriteJobToJobsFinishedFile( j );
            DeleteJob( j );
            //		RewriteToDoFile();
        }
    }
    catch( std::string& err )
    {
        std::stringstream s;
        s << "Sorry. An error occured.\n" << err << std::endl;
        QString messErr = QString::fromStdString(s.str());

        emit aborted( messErr );
    }
    catch( itk::ExceptionObject & err )
    {
        std::stringstream s;
        s << "Sorry. An error occured.\n" << err << std::endl;
        QString messErr = QString::fromStdString(s.str());

        emit aborted( messErr );
    }
    catch( std::exception& err )
    {
        std::stringstream s;
        s << "Sorry. An error occured.\n" << err.what() << std::endl;
        QString messErr = QString::fromStdString(s.str());

        emit aborted( messErr );
    }

	std::cout << "Oh joy, I finished all jobs!" << std::endl;
    emit finished();
}


void JobManager::SetupGUI( QTableView * tableView )
{
    qtModel = new JobManagerModel( this );
    tableView->setModel( qtModel );
    connect( qtModel, SIGNAL( modelReset() ),
             tableView, SLOT( resizeRowsToContents() ) );
    connect( qtModel, SIGNAL( modelReset() ),
             tableView, SLOT( resizeColumnsToContents() ) );
}


void JobManager::RewriteToDoFile() const
{
    QFile outfile( QString((jobFilePath + "_temp_jobsToDo.txt").c_str()) );

    outfile.open( QFile::WriteOnly | QFile::Truncate | QFile::Text );

    QXmlStreamWriter xmlStream( &outfile);
    xmlStream.setAutoFormatting(true);
    jobQueue->writeXML( &xmlStream );

    outfile.close();

	std::remove((jobFilePath + "jobsToDo.txt").c_str());
	std::rename((jobFilePath + "_temp_jobsToDo.txt").c_str(), (jobFilePath + "jobsToDo.txt").c_str());
}


void JobManager::WriteJobToJobsFinishedFile(Job * job)
{
	std::fstream file;

	file.open((jobFilePath + "jobsFinished.txt").c_str(), std::fstream::out | std::fstream::app);

    CSParameter * typeParm = job->findParameter("type");
    std::string type = JobTypeString[ *(int*)typeParm->dataPointer() ];

	file.width(10);
	file << job->name();
	file.width(30);
	file << type;
	file.width(100);

    std::stack<CSParameterContext *> contextStack;
    contextStack.push( (CSParameterContext *) job );

    while ( !contextStack.empty() )
    {
        CSParameterContext * currentContext = contextStack.top();
        contextStack.pop();

        std::vector<CSParameter *> parms = currentContext->getParameters();
        std::vector<CSParameter *>::iterator parmsIt = parms.begin();
        while ( parmsIt != parms.end() )
        {
            if ( (*parmsIt)->name() == "type" )
            {
                ++parmsIt;
                continue;
            }

            file << (*parmsIt)->name() << "=";
            file << (*parmsIt)->dataString() << ";";
            ++parmsIt;
        }

        std::vector<CSParameterContext *> ctxts = currentContext->getChildren();
        std::vector<CSParameterContext *>::iterator ctxt = ctxts.begin();
        while ( ctxt != ctxts.end() )
        {
            contextStack.push( *ctxt );
            ++ctxt;
        }
    }

    file << "\n";

	file.close();
}


void JobManager::startPreprocessDataSet(Job * job, unsigned int algorithm, unsigned int dataDim)
{
    if(dataDim == 2) {                                          //template arguments have to be const :/
        PreprocessDataSet<2>* preprocessDataSet = new PreprocessDataSet<2>();

        preprocessDataSet->SetAlgorithm(PreprocessDataSet<2>::AlgorithmType(algorithm));
        preprocessDataSet->SetParameterContext(job);
        preprocessDataSet->Update();

        delete preprocessDataSet;
    }
    else if(dataDim == 3) {
        PreprocessDataSet<3>* preprocessDataSet = new PreprocessDataSet<3>();

        preprocessDataSet->SetAlgorithm(PreprocessDataSet<3>::AlgorithmType(algorithm));
        preprocessDataSet->SetParameterContext(job);
        preprocessDataSet->Update();

        delete preprocessDataSet;
    }
}

void JobManager::startCropDataSet(Job * job, unsigned int dataDim)
{
    if(dataDim == 2) {                                          //template arguments have to be const :/
        CropDataSet<2>* cropDataSet = new CropDataSet<2>();

        cropDataSet->SetParameterContext(job);
        cropDataSet->Update();

        delete cropDataSet;
    }
    else if(dataDim == 3) {
        CropDataSet<3>* cropDataSet = new CropDataSet<3>();

        cropDataSet->SetParameterContext(job);
        cropDataSet->Update();

        delete cropDataSet;
    }
}

void JobManager::startSegmentNecroticRegion(Job * job)
{
    SegmentNecroticRegionOnDM* extractNecroticRegion = new SegmentNecroticRegionOnDM();

    extractNecroticRegion->SetParameterContext(job);
    extractNecroticRegion->Update();
}

void JobManager::startSegmentNuclei20x(Job * job)
{
    SegmentNucleiOnDAPI20x* segmentNucleiOnDAPI20x = new SegmentNucleiOnDAPI20x();

    segmentNucleiOnDAPI20x->SetParameterContext(job);
    segmentNucleiOnDAPI20x->Update();
}

void JobManager::startSegmentNuclei60x(Job * job, unsigned int dataDim)
{
    if(dataDim == 2) {                                          //template arguments have to be const :/
        SegmentNucleiOnDAPI60x<2>* segmentNucleiOnDAPI60x = new SegmentNucleiOnDAPI60x<2>();

        segmentNucleiOnDAPI60x->SetParameterContext(job);
        segmentNucleiOnDAPI60x->Update();

        delete segmentNucleiOnDAPI60x;
    }
    else if(dataDim == 3) {
        SegmentNucleiOnDAPI60x<3>* segmentNucleiOnDAPI60x = new SegmentNucleiOnDAPI60x<3>();

        segmentNucleiOnDAPI60x->SetParameterContext(job);
        segmentNucleiOnDAPI60x->Update();

        delete segmentNucleiOnDAPI60x;
    }
}

void JobManager::startSegmentNucleiWithHough(Job * job)
{
    SegmentNucleiOnDAPIWithHough* segmentNucleiOnDAPIWithHough = new SegmentNucleiOnDAPIWithHough();

    segmentNucleiOnDAPIWithHough->SetParameterContext(job);
    segmentNucleiOnDAPIWithHough->Update();
}

void JobManager::startAnalyzeCells(Job * job)
{
    AnalyzeCellsFilter* anaylzeCells = new AnalyzeCellsFilter();

    anaylzeCells->SetParameterContext(job);
    anaylzeCells->Update();
}

void JobManager::startAnalyzeLobule(Job * job)
{
    AnalyzeLobuleFilter* analyzeLobule = new AnalyzeLobuleFilter();

    analyzeLobule->SetParameterContext(job);
    analyzeLobule->Update();
}

void JobManager::startAnalyzeNuclei(Job * job, unsigned int dataDim)
{
    if(dataDim == 2) {                                   //template arguments have to be const :/
        AnalyzeNucleiFilter<2>* analyzeNuclei = new AnalyzeNucleiFilter<2>();

        analyzeNuclei->SetParameterContext(job);
        analyzeNuclei->Update();

        delete analyzeNuclei;
    }
    else if(dataDim == 3) {                                   //template arguments have to be const :/
        AnalyzeNucleiFilter<3>* analyzeNuclei = new AnalyzeNucleiFilter<3>();

        analyzeNuclei->SetParameterContext(job);
        analyzeNuclei->Update();

        delete analyzeNuclei;
    }
}


void JobManager::startEstimateCellShape(Job * job)
{
    EstimateCellShape* estimateCellShape = new EstimateCellShape();

    estimateCellShape->SetParameterContext(job);
    estimateCellShape->Update();

    delete estimateCellShape;
}

void JobManager::startEstimateLobuleShape(Job * job, unsigned int dataDim)
{
    if(dataDim == 2) {                                   //template arguments have to be const :/
        EstimateLobuleShape<2>* estimateLobuleShape = new EstimateLobuleShape<2>();

        estimateLobuleShape->SetParameterContext(job);
        estimateLobuleShape->Update();

        delete estimateLobuleShape;
    }
    else if(dataDim == 3) {
        EstimateLobuleShape<3>* estimateLobuleShape = new EstimateLobuleShape<3>();

        estimateLobuleShape->SetParameterContext(job);
        estimateLobuleShape->Update();

        delete estimateLobuleShape;
    }
}

void JobManager::startSegmentAndClassifyStellateCells(Job * job)
{
    SegmentAndClassifyStellateCells* segmentStellateCells = new SegmentAndClassifyStellateCells();

    segmentStellateCells->SetParameterContext(job);
    segmentStellateCells->Update();
}

void JobManager::startSinusoidBileSegmentation20x(Job * job)
{
    SegmentSinusoidalNetworkOnTwoChannels20x* extractSinusNetwork = new SegmentSinusoidalNetworkOnTwoChannels20x();
    SegmentBileNetworkOnDPPIV20x* extractBileNetwork = new SegmentBileNetworkOnDPPIV20x();

    extractSinusNetwork->SetParameterContext(job);
    extractSinusNetwork->Update();

    extractBileNetwork->SetParameterContext(job);
    extractBileNetwork->Update();
}

void JobManager::startSinusoidBileSegmentation60x(Job * job)
{
    SegmentSinusoidalNetworkOnTwoChannels60x* extractSinusNetwork = new SegmentSinusoidalNetworkOnTwoChannels60x();
    SegmentBileNetworkOnDPPIV60x* extractBileNetwork = new SegmentBileNetworkOnDPPIV60x();

    extractSinusNetwork->SetParameterContext(job);
    extractSinusNetwork->Update();

    extractBileNetwork->SetParameterContext(job);
    extractBileNetwork->Update();
}

void JobManager::startSuperpixelPreparation(Job * job, unsigned int dataDim)
{
    if(dataDim == 2) {                                   //template arguments have to be const :/
        SuperpixelPreparation<2>* superpixelPreparationPipeline = new SuperpixelPreparation<2>();

        superpixelPreparationPipeline->SetParameterContext(job);
        superpixelPreparationPipeline->Update();

        delete superpixelPreparationPipeline;
    }
    else {
        SuperpixelPreparation<3>* superpixelPreparationPipeline = new SuperpixelPreparation<3>();

        superpixelPreparationPipeline->SetParameterContext(job);
        superpixelPreparationPipeline->Update();

        delete superpixelPreparationPipeline;
    }
}

void JobManager::startObjectBasedSegmentation(Job * job, unsigned int dataDim)
{
    if(dataDim == 2) {                                   //template arguments have to be const :/
        ObjectBasedSegmentation<2>* objectBasedSegPipeline = new ObjectBasedSegmentation<2>();

        objectBasedSegPipeline->SetParameterContext(job);
        objectBasedSegPipeline->Update();

        delete objectBasedSegPipeline;
    }
    else {
        ObjectBasedSegmentation<3>* objectBasedSegPipeline = new ObjectBasedSegmentation<3>();

        objectBasedSegPipeline->SetParameterContext(job);
        objectBasedSegPipeline->Update();

        delete objectBasedSegPipeline;
    }
}

void JobManager::startSkeletonization(Job * job)
{
    SkeletonizationPipeline* skeletonizeDataSet = new SkeletonizationPipeline();

    skeletonizeDataSet->SetParameterContext(job);
    skeletonizeDataSet->Update();
}

void JobManager::startExtractAndAnalyzeGraph(Job * job)
{
    ExtractGraph *pipeline = new ExtractGraph();
    pipeline->SetParameterContext(job);
    pipeline->Update();

    std::string analysisMode = ( (CSParameterChoice*)(job->findParameter("Analysis mode", 0)->dataPointer()) )->currentString();
    bool withAnalysis = (analysisMode.compare("With Analysis")==0);

    AnalyzeBileNetworkFilter *analyzeBile = new AnalyzeBileNetworkFilter();
    if(withAnalysis) {
        analyzeBile->SetParameterContext(job);
        analyzeBile->SetGraphsToAnalyze(pipeline->GetFinalGraphs());
        analyzeBile->SaveGraphsWithArrays(true);
        analyzeBile->SetGraphFilename(pipeline->GetNameOfFirstGraph());
        analyzeBile->Update();
    }
}

void JobManager::startSegmentCellMembrane(Job * job)
{
    SegmentCellMembraneOnBCat60x* segmentCellMembrane = new SegmentCellMembraneOnBCat60x();

    segmentCellMembrane->SetParameterContext(job);
    segmentCellMembrane->Update();
}

void JobManager::startFileFormatConversionJob(Job * job)
{
    FileFormatConverter* pFileFormatConverter = new FileFormatConverter();

    pFileFormatConverter->SetParameterContext(job);
    pFileFormatConverter->Update();
}
