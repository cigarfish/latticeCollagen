

#include "Core.h"
#include "tools/parameters/CSParameterContext.h"
#include "tools/Tools.h"


#ifdef CS_BUILD_PYTHON
    PyThreadState * Core::pMainThreadState = NULL;
#endif

Core::Core()
{
 // Tools (Random, OutputText, ...)
    tools = new Tools();
}

void Core::init()
{
    //  global parameter tree
    // initialize the root node before anything else!
    parameters = new CSParameterContext( "Root" );
// Init 2D monolayer model
//    models["ModelCellsSpherical"] = new ModelCellsSpherical();
}

bool Core::initPythonInterpreter()
{
    // Initialize the Python Interpreter
#ifdef CS_BUILD_PYTHON
    if(pMainThreadState == NULL) {
        Py_Initialize();
        PyEval_InitThreads();

        pMainThreadState = PyThreadState_Get();
        PyEval_ReleaseLock();

        return true;
    }
#endif
    return false;
}
