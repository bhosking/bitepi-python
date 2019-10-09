#define PY_SSIZE_T_CLEAN
#include <Python.h>

#include <vector>

#include "BitEpi.cpp"

static PyObject * main_wrapper(PyObject * self, PyObject * args)
{
    int result;
    char *delim;
    char *argchars;
    char *token;
    char *saveptr;
    PyObject *ret;
    std::vector<char *> argvect;

    // parse arguments
    if (!PyArg_ParseTuple(args, "ss", &delim, &argchars)) {
        return NULL;
    }
    int j;
    for (j=1; ; j++, argchars = NULL)
    {
        token = strtok_r(argchars, delim, &saveptr);
        if (token == NULL)
        {
            break;
        }
        argvect.push_back(token);
    }

    // run the actual function
    result = main(int(argvect.size()), &argvect[0]);

    // build the resulting string into a Python object.
    ret = PyLong_FromLong(result);

    return ret;
}

static PyMethodDef BitEpiMethods[] = {
    { "bitepi", main_wrapper, METH_VARARGS, "Analyse with BitEpi" },
    { NULL, NULL, 0, NULL }
};

static struct PyModuleDef BitEpiModule = {
    PyModuleDef_HEAD_INIT,
    "bitepimodule",
    NULL,
    -1,
    BitEpiMethods
};

PyMODINIT_FUNC PyInit_bitepimodule(void)
{
    return PyModule_Create(&BitEpiModule);
}
