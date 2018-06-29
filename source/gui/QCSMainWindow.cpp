
#include "QCSMainWindow.h"

#include <QtGui>

#include "QCSGLDisplay.h"
#include "QCSCentralWidget.h"

#include "../model/Model/CSModel.h"

#include "../tools/dataIO/CSDataWriter.h"
#include "../tools/dataIO/CSDataReader.h"

#include "tabMonolayer/monolayer.h"

#if QT_VERSION >= 0x050000
# include <QtWidgets/QtWidgets>
#else
# include <QtCore>
#endif

#ifndef __BUILD_MAC__
# include <qapplication.h>
#endif

#include <sstream>
#include <iostream>


/*
  \brief The constructor

  setting the central widget and constructing the menus in the menuBar
*/
QCSMainWindow::QCSMainWindow(QWidget * parent, Qt::WindowFlags flags)
    : QMainWindow(parent,flags),
      mDataDirName("")
{
    resize(1000, 600);

    setCentralWidget(new QCSCentralWidget(this));

    // Adding menus to the menu bar

#ifndef CS_TI_QUANT_ONLY
    // the File menu:
    QMenu * fileMenu = menuBar()->addMenu("File");

    QAction * actionFileOpen = new QAction(tr("&Open.."),this);
    actionFileOpen->setShortcuts(QKeySequence::Open);
    connect(actionFileOpen,SIGNAL(triggered()),this,SLOT(fileOpenDialog()));
    fileMenu->addAction(actionFileOpen);

    QAction * actionFileSave = new QAction("&Save",this);
    actionFileSave->setShortcuts(QKeySequence::Save);
    connect(actionFileSave,SIGNAL(triggered()),this,SLOT(saveToFile()));
    fileMenu->addAction(actionFileSave);

    QAction * actionFileSaveAs = new QAction("&Save as..",this);
    actionFileSaveAs->setShortcuts(QKeySequence::SaveAs);
    connect(actionFileSaveAs,SIGNAL(triggered()),this,SLOT(saveToOtherFile()));
    fileMenu->addAction(actionFileSaveAs);

    // Under Mac OS X a quit function is added automatically
#ifndef __BUILD_MAC__
    fileMenu->addSeparator();

    QAction * actionQuit = new QAction("&Quit",this);
    actionQuit->setShortcuts(QKeySequence::Quit);
    connect(actionQuit, SIGNAL(triggered()), qApp, SLOT(quit()));
    fileMenu->addAction(actionQuit);
#endif //__BUILD_MAC__

    #pragma region View menu

    // View menu
    // The view menu should be used to:
    // Switch between what is shown in the 2D or 3D tabs (e.g. in 2D it could be image processing or model cells)
    // Reset view and assume some predefined views
    QMenu *viewMenu = menuBar()->addMenu("View");

    // Preliminary: Put this in submenus
    QAction * actionViewMonolayer2D = new QAction("2D: Monolayer",this);
    connect(actionViewMonolayer2D,SIGNAL(triggered()),this,SLOT(messageNotImplementedYet()));
    viewMenu->addAction(actionViewMonolayer2D);

    QAction * actionViewImage2D = new QAction("2D: Image",this);
    connect(actionViewImage2D,SIGNAL(triggered()),this,SLOT(messageNotImplementedYet()));
    viewMenu->addAction(actionViewImage2D);

    #pragma endregion

    // the Tools menu
    QMenu *toolsMenu = menuBar()->addMenu("Tools");

    QAction * actionTools3DCut = new QAction("3D &Cut",this);
    QList<QKeySequence> keyList3DCut;
    keyList3DCut.append(QKeySequence(tr("Ctrl+c")));
    actionTools3DCut->setShortcuts( keyList3DCut );
    connect(actionTools3DCut,SIGNAL(triggered()),this,SLOT(messageNotImplementedYet()));
    toolsMenu->addAction(actionTools3DCut);

    // dummy entries after a seperator
    toolsMenu->addSeparator();

    QAction * actionDummy1 = new QAction("Dummy Entry 1",this);
    connect(actionDummy1,SIGNAL(triggered()),this,SLOT(messageNotImplementedYet()));
    toolsMenu->addAction(actionDummy1);

    QAction * actionDummy2 = new QAction("Dummy Entry 2",this);
    connect(actionDummy2,SIGNAL(triggered()),this,SLOT(messageNotImplementedYet()));
    toolsMenu->addAction(actionDummy2);


    // help!!
    QMenu * helpMenu = menuBar()->addMenu("Help");

    QAction *actionHelp = new QAction("CellSys &Help",this);
    actionHelp->setShortcuts(QKeySequence::HelpContents);
    connect(actionHelp,SIGNAL(triggered()),this,SLOT(messageNotImplementedYet()));
    helpMenu->addAction(actionHelp);
#endif // CS_TI_QUANT_ONLY
}


QCSMainWindow::~QCSMainWindow()
{}


/*!
  \brief Set the active display
  \param display new active display
*/
void
QCSMainWindow::setDisplay(QCSGLDisplay *display)
{
    if (!display)
        return;

    // first, allow only one display
    if (mpActiveDisplay != NULL)
        delete mpActiveDisplay;

    mpActiveDisplay = display;
}


/*!
  \brief Get the active display
 */
QCSGLDisplay *
QCSMainWindow::getDisplay()
{
    return mpActiveDisplay;
}


/*!
  \brief Pop up a dialog to choose a file to open

  used as a slot
*/
void
QCSMainWindow::fileOpenDialog()
{
    QFileDialog openDialog(this);
#ifndef __BUILD_MAC__
    openDialog.setFileMode( QFileDialog::Directory );
#endif
    openDialog.setNameFilter( tr("CellSys data bundle (*.csx *.csys)") );
    openDialog.setDirectory( QDir::currentPath() );
    openDialog.setAcceptMode( QFileDialog::AcceptOpen );

    QString dataDirName;

    if ( openDialog.exec() )
        dataDirName = openDialog.selectedFiles().at(0);
    else
        return;

    if ( dataDirName.isEmpty() )
        return;

    //messageNotImplementedYet();
    QFile inputFile( dataDirName );
    QFileInfo finfo(inputFile);
    if ( !finfo.isDir() )
    {
        // maybe call the old xml file parser from here before issuing an error.
        QMessageBox::warning( this, tr("Open CellSys data bundle"),
                              tr("Cannot open file %1:\n%2")
                              .arg( dataDirName )
                              .arg("Not a CellSys Data bundle.") );
        return;
    }

    mDataDirName = dataDirName;
    std::string dataDirNameStdString = mDataDirName.toStdString();

    CSDataReader dataReader( dataDirNameStdString );
    std::stringstream errors;
    std::stringstream warnings;
    std::vector<CSModel *> newModels = dataReader.exec( errors, warnings );

    if (errors.str().size())
    {
        QMessageBox warning;
        warning.setText(QString(errors.str().c_str()));

        warning.exec();
    }

    std::cerr << errors.str() << std::endl;

    monolayer * tabMonolayer = dynamic_cast<monolayer *>( findChild<monolayer *>("tabMonolayer") );
    std::vector<CSModel *>::iterator modelIt;
    std::cout << "Got the following Models:\n";
    for (modelIt = newModels.begin(); modelIt != newModels.end(); ++modelIt)
    {
        std::cout << (*modelIt)->name /* << " of type " << (*modelIt)->xmlType*/ << std::endl;
        if ( tabMonolayer )
            tabMonolayer->initFromModel(*modelIt);
    }
}


/*!
  \brief Save to the currently chosen file.

  If there is no open file, open the same action as 'save as', i.e. call saveToOtherFile().
 */
void
QCSMainWindow::saveToFile()
{
    if ( mDataDirName.size() )
    {
        if ( !mSaved )
        {
            QMessageBox forewarn;
            forewarn.setStandardButtons( QMessageBox::Cancel | QMessageBox::Save );
            forewarn.setDefaultButton( QMessageBox::Cancel );
            forewarn.setText(tr("You have not saved the data yet.\n\nAre you sure you want to overwrite the existing data?"));
            if ( forewarn.exec() == QMessageBox::Cancel )
                return;
        }

        // QMessageBox warn;
        // warn.setText(tr("Save what??\nThe World?  Our Souls?\nHaha, you really thought we can save your work?"));
        // warn.exec();
        mSaved = true;
    }
    else
        saveToOtherFile();
}


/*!
  \brief Pop up a dialog to choose a file name or file to save under
*/
void
QCSMainWindow::saveToOtherFile()
{
    QString directoryName;

    QFileDialog saveDialog(this);
    saveDialog.setFileMode( QFileDialog::AnyFile );
    saveDialog.setNameFilter( tr("CellSys data bundle (*.csx *.csys)") );
    saveDialog.setDirectory( QDir::currentPath() );
    saveDialog.setAcceptMode( QFileDialog::AcceptSave );

    if ( saveDialog.exec() )
        directoryName = saveDialog.selectedFiles().at(0);
    else
        return;

    if (directoryName.size())
    {
        QFile testFile( directoryName );

        // has it got the right suffix?  add lmao by default if not!
        if ( !directoryName.endsWith(".csx") && !directoryName.endsWith(".csys") )
            directoryName.append(".csx");

        QFileInfo finfo(testFile);

        // was a directory chosen?
        if ( testFile.exists() )
            if ( !finfo.isDir() )
            {
                QMessageBox::warning(this, tr("Saving Data"),
                                     tr("Cannot create data directory %1:\n%2.")
                                     .arg(mDataDirName)
                                     .arg("A file of the same name exists!"));
                return;
            }

        QDir upperDirectory = finfo.absoluteDir();

        std::cout << "creating path " << directoryName.toStdString() << std::endl;

        if ( ! upperDirectory.mkpath( directoryName ) )
          {
            QMessageBox::warning(this, tr("Saving Data"),
                                 tr("Cannot create data directory %1:\n%2.")
                                 .arg(mDataDirName)
                                 .arg("Unknown error."));
            return;
          }

        mDataDirName = directoryName;

        std::string directoryNameStdString = mDataDirName.toStdString();

        CSDataWriter dataWriter( directoryNameStdString );
        dataWriter.exec();

        mSaved = true;
    }
}


/*!
  \brief Display a message for not-yet-implemented features
*/
void
QCSMainWindow::messageNotImplementedYet()
{
    QMessageBox msg;
    msg.setText(tr("Sorry!\n\nThis functionality is not implemented yet."));
    msg.exec();
}


void
QCSMainWindow::closeEvent( QCloseEvent * event )
{
    QMessageBox message;

#ifndef CS_TI_QUANT_ONLY
    message.setText("You are about to quit CellSys.");
#else
    message.setText("You are about to quit TiQuant.");
#endif
    message.setStandardButtons( QMessageBox::Ok | QMessageBox::Cancel );

    int close = message.exec();

    if ( close == QMessageBox::Cancel )
    {
        event->ignore();
        return;
    }

    qApp->quit();
}
