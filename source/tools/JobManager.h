///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  JobManager.h                                                         //
//                                                                                   //
//     Author:  Adrian Friebel <adrian.friebel@uni-leipzig.de>                       //
//    Created:  2012-12-20                                                           //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

#include <string>
#include <vector>

#include <QtCore>

class QTableView;
class CSParameterContext;
class CSParameterContextTemporary;
class JobManagerModel;


class JobManager : public QObject
{
    Q_OBJECT

    typedef CSParameterContextTemporary Job;
	typedef CSParameterContext JobQueueType;

public:
    enum JobType 
    {
        CLAHEJob,
        BackgroundEliminationJob,
        HConvexImageFilterJob,
        CropDataSetJob,
        SegmentNecroticRegionJob,
        SegmentNuclei20xJob,
        SegmentNuclei60xJob,
        SegmentNucleiWithHoughJob,
        EstimateCellShapeJob,
        EstimateLobuleShapeJob,
        AnalyzeCellsJob,
        AnalyzeLobuleJob,
        AnalyzeNucleiJob,
        SegmentAndClassifyStellateCellsJob,
        SegmentSinusoidsBile20xJob,
        SegmentSinusoidsBile60xJob,
        SuperpixelPreparationJob,
        ObjectBasedSegmentationJob,
        SkeletonizationJob,
        ExtractAndAnalyzeGraphJob,
        SegmentCellMembraneJob,
        FileFormatConversionJob
    };

    // update the following (in the .cpp) when adding JobTypes!
    static std::string JobTypeString[];


public:
	JobManager(void);
	~JobManager(void);

	void AddJob(JobType type, CSParameterContext * parameters);
	void AddJobs(JobType type, CSParameterContext * parameters);
	void DeleteJob( Job * job );

    JobQueueType *Queue() const { return jobQueue; };
	int GetNumberOfQueuedJobs() const;
	Job * GetJob(int id) const;

	void InitWithJobQueueFile(std::string path);

    void SetupGUI( QTableView * );

    std::string GetJobFilePath() const {return jobFilePath;};

public slots:
    void Start();


signals:
    void aboutToChange();
    void changed();
    void aborted(QString reason);
    void finished();

private:
	void RewriteToDoFile() const;
    void WriteJobToJobsFinishedFile(Job * job);

	void startPreprocessDataSet(Job * job, unsigned int algorithm, unsigned int dataDim);
	void startCropDataSet(Job * job, unsigned int dataDim);
	void startSegmentNecroticRegion(Job * job);
	void startSegmentNuclei20x(Job * job);
	void startSegmentNuclei60x(Job * job, unsigned int dataDim);
	void startSegmentNucleiWithHough(Job * job);
	void startAnalyzeCells(Job * job);
	void startAnalyzeLobule(Job * job);
	void startAnalyzeNuclei(Job * job, unsigned int dataDim);
	void startEstimateCellShape(Job * job);
	void startEstimateLobuleShape(Job * job, unsigned int dataDim);
	void startSegmentAndClassifyStellateCells(Job * job);
	void startSinusoidBileSegmentation20x(Job * job);
	void startSinusoidBileSegmentation60x(Job * job);
	void startSuperpixelPreparation(Job * job,  unsigned int dataDim);
	void startObjectBasedSegmentation(Job * job,  unsigned int dataDim);
	void startSkeletonization(Job * job);
	void startExtractAndAnalyzeGraph(Job * job);
	void startSegmentCellMembrane(Job * job);
	void startFileFormatConversionJob(Job * job);

	std::string jobFilePath;

	JobQueueType * jobQueue;
    JobManagerModel * qtModel;
};

