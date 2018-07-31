#ifndef PYTHON_HELPER_H
#define PYTHON_HELPER_H

#include <string>
#include <vector>

#include <Python.h>


class PythonInterpreter {

public:
    PyObject* main;
    wchar_t *program = NULL;
    PythonInterpreter(int argc, char *argv[]);
    PythonInterpreter();
    ~PythonInterpreter();
    bool import_module(std::string module);
    int run_string(std::string command);
    PyObject* run_string_get_obj(std::string command);
    int run_strings( std::vector< std::string > statements);
    int run_file(FILE *fp, const char *filename);
};

#endif // PYTHON_HELPER_H

