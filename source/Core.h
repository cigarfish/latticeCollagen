
#ifndef CORE_H
#define CORE_H

#include <iomanip>
#include <map>

#ifdef CS_BUILD_PYTHON
#include <Python.h>
#endif

// due to some circular include snafu, we have to do this:
class CSModel;
class Tools;
class CSParameterContext;


//! Collection of all backend classes
class Core
{
public:
    Core();

	void init();
	bool initPythonInterpreter();

	// Tools (Random, OutputText, ...)
    Tools * tools;

    // Model with spherical cells (Strict monolayer, Piling up monolayer and tumor spheroids)
    std::map<std::string, CSModel *> models;

    // The global tree of parameters and contexts:
    CSParameterContext * parameters;

#ifdef CS_BUILD_PYTHON
    static PyThreadState * pMainThreadState;
#endif
};


#ifdef _CS_MAIN_
Core *core;
#else
extern Core *core;
#endif // _CS_MAIN_

#endif
