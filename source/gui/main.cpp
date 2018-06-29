
#define _CS_MAIN_

#include <qapplication.h>

#ifdef __BUILD_MAC__
# include <GLUT/glut.h>
#else
# include <GL/glut.h>
#endif

#include <iostream>

#include "QCSMainWindow.h"
#include "../tools/batchJob/CSBatchJob.h"

#include "../Core.h"

#if !defined ( __BUILD_WINDOWS__ ) && ! defined ( __BUILD_XCODE__ )
// build information: about commit, branch, and modifications.
#include "../tools/dataIO/buildInformation.h"
#endif

int main(int argc, char **argv)
{
    core = new Core();
    core->init();
#ifdef CS_BUILD_PYTHON
    core->initPythonInterpreter();
#endif

#if !defined( __BUILD_WINDOWS__ ) && !defined( __BUILD_XCODE__ )
    buildInformation(argv[0]);
#endif

    for ( int i=1; i<argc; ++i )
    {
        std::string argument = argv[i];
        if ( argument == "--batch" )
        {
            ++i;
            if ( i >= argc )
            {
                std::cout << "\n\nCellSys Error:  \"--batch\" needs a file from which to read the job description.\n\n";
                delete core;
                return 1;
            }
            // read next Argument -> xml model description
            argument = argv[i];
            CSBatchJob job( argument );

            if ( !job.isValid() )
                return 2;

            return job.run();
            // return 0;
        }
    }

    // Init and Start QT
    QApplication mainApp(argc, argv);

#if ! defined( __BUILD_MAC__ )
    glutInit(&argc,argv);
#endif

    QMainWindow * mainDisplay = new QCSMainWindow();

    mainDisplay->show();

    return mainApp.exec();
}
