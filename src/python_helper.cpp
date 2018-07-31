
#include <string>
#include <vector>
#include <iostream>

#include <Python.h>
#include "python_helper.h"

PythonInterpreter::PythonInterpreter(int argc, char *argv[]) {
    this->program = Py_DecodeLocale(argv[0], NULL);
    if (this->program == NULL) {
        fprintf(stderr, "Fatal error: cannot decode argv[0]\n");
        exit(1);
    }
    Py_SetProgramName(this->program);  /* optional but recommended */
    Py_Initialize();
    this->main = PyImport_AddModule("__main__");
}

PythonInterpreter::PythonInterpreter() {
    Py_Initialize();
    this->main = PyImport_AddModule("__main__");
}

PythonInterpreter::~PythonInterpreter() {
    Py_Finalize();
    if (this->program != NULL)
        PyMem_RawFree(this->program);
}

int PythonInterpreter::run_string(std::string command) {
    return PyRun_SimpleString(command.c_str());
}

PyObject* PythonInterpreter::run_string_get_obj(std::string command){
    PyObject *pGlobal, *pLocal, *result;
    pGlobal = PyModule_GetDict(this->main);
    pLocal = PyDict_New();
    result = PyRun_String(command.c_str(), Py_eval_input,pGlobal , PyImport_GetModuleDict() );
    return result;
}

int PythonInterpreter::run_strings( std::vector< std::string > statements){
    int success;
    for (int i=0; i< statements.size(); i++)    {
        success = this->run_string(statements[i]);
        if (success != 0)
            break;
    }
    return success;
}


