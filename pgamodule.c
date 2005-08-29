#include <Python.h>
#include <pgapack.h>
#include <stdio.h>

static PyObject *context = NULL;

/*
 * Need a hash table of mapping ctx to PGA objects. Look up the
 * appropriate object and call its PGA_evaluate
 */
static double evaluate (PGAContext *ctx, int p, int pop)
{
    double retval;
    PyObject *PGA_ctx, *o, *result;

    PGA_ctx = PyCObject_FromVoidPtr  (ctx, NULL);
    o       = PyObject_GetItem       (context, PGA_ctx);
    result  = PyNumber_Float
        (PyObject_CallMethod (o, "evaluate", "ii", p, pop));
    /* FIXME: decrement refcount? */
    retval  = PyArg_Parse (result, "d");
    return retval;
}

/* name, type, length, maximize */
static PyObject *PGA_init (PyObject *self, PyObject *args)
{
    int argc = 0, maximize = 0, length = 0;
    char *name;
    PyObject *type, *max = NULL, *PGA_ctx;
    char *argv [] = {NULL, NULL};
    PGAContext *ctx;

    if (!PyArg_ParseTuple (args, "sOi|O", &name, &type, &length, &max))
        return NULL;
    if (max)
    {
        maximize = PyObject_IsTrue (max);
    }
    argv [0] = name;
    
    ctx = PGACreate
        ( &argc
        , argv
        , PGA_DATATYPE_BINARY
        , length
        , maximize ? PGA_MAXIMIZE : PGA_MINIMIZE
        );
    PGASetUp (ctx);
    PGAPrintContextVariable (ctx, stdout);
    PGA_ctx = PyCObject_FromVoidPtr ((void  *) ctx, NULL);
    PyObject_SetItem       (context, PGA_ctx, self);
    PyObject_SetAttrString (self, "context", PGA_ctx);
    return Py_BuildValue("");
}

static PyObject *PGA_run (PyObject *self, PyObject *args)
{
    PyObject   *PGA_ctx = PyObject_GetAttrString (self, "context");
    PGAContext *ctx     = (PGAContext *) PyCObject_AsVoidPtr (PGA_ctx);

    if (!PyArg_ParseTuple(args, ""))
        return NULL;
    PGARun (ctx, evaluate);
    return Py_BuildValue("");
}

static PyObject *PGA_len (PyObject *self, PyObject *args)
{
    if (!PyArg_ParseTuple(args, ""))
        return NULL;
    return Py_BuildValue("i", 123);
}

static PyObject *PGA_evaluate (PyObject *self, PyObject *args)
{
    int p, pop;
    if (!PyArg_ParseTuple(args, "ii", &p, &pop))
        return NULL;
    return Py_BuildValue ("i", 0);
}

static PyObject *PGA_get_allele (PyObject *self, PyObject *args)
{
    int p, pop, i, allele;
    PyObject   *PGA_ctx = PyObject_GetAttrString (self, "context");
    PGAContext *ctx     = (PGAContext *) PyCObject_AsVoidPtr (PGA_ctx);

    if (!PyArg_ParseTuple(args, "iii", &p, &pop, &i))
        return NULL;
    allele = PGAGetBinaryAllele (ctx, p, pop, i);
    return Py_BuildValue ("i", allele);
}


/*
__del__
    PGADestroy (ctx)
*/

static PyMethodDef PGA_Methods [] =
{ {"__init__",   PGA_init,       METH_VARARGS, "Init object"}
, {"__len__",    PGA_len,        METH_VARARGS, "Return length of gene"}
, {"evaluate",   PGA_evaluate,   METH_VARARGS, "Evaluate"}
, {"get_allele", PGA_get_allele, METH_VARARGS, "Get allele"}
, {"run",        PGA_run,        METH_VARARGS, "Run optimization"}
, {NULL, NULL, 0, NULL}
};

PyMODINIT_FUNC initpga (void)
{
    (void) Py_InitModule("pga", PGA_Methods);
    context = Py_BuildValue("{}");
}
