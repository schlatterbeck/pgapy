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

    PGA_ctx = Py_BuildValue ("i", (int) ctx);
    o       = PyObject_GetItem       (context, PGA_ctx);
    result  = PyObject_CallMethod (o, "evaluate", "ii", p, pop);
    result  = PyNumber_Float (result);
    /* FIXME: decrement refcount? */
    PyArg_Parse (result, "d", &retval);
    return retval;
}

/* name, type, length, maximize */
static PyObject *PGA_init (PyObject *self0, PyObject *args)
{
    int argc = 0, maximize = 0, length = 0;
    char *name;
    PyObject *self, *type, *max = NULL, *PGA_ctx;
    char *argv [] = {NULL, NULL};
    PGAContext *ctx;

    if (!PyArg_ParseTuple (args, "OsOi|O", &self, &name, &type, &length, &max))
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
    printf ("ctx: %08X; self: %08X\n", (int)ctx, (int)self);
    fflush (stdout);
    PGASetUp (ctx);
    PGA_ctx = Py_BuildValue ("i", (int) ctx);
    printf ("init: PGA_ctx: %08X\n", (int)PGA_ctx);
    fflush (stdout);
    PyObject_SetItem       (context, PGA_ctx, self);
    PyObject_SetAttrString (self, "context", PGA_ctx);
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *PGA_run (PyObject *self0, PyObject *args)
{
    PyObject *self;
    PyObject   *PGA_ctx;
    PGAContext *ctx;

    if (!PyArg_ParseTuple(args, "O", &self))
        return NULL;
    PGA_ctx = PyObject_GetAttrString (self, "context");
    printf ("PGA_run, PGA_ctx: %08X\n", (int)PGA_ctx);
    PyArg_Parse (PGA_ctx, "i", &ctx);
    printf ("PGA_run, ctx: %08X\n", (int)ctx);
    fflush (stdout);
    PGARun (ctx, evaluate);
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *PGA_len (PyObject *self0, PyObject *args)
{
    PyObject *self;
    PyObject   *PGA_ctx;
    PGAContext *ctx;

    if (!PyArg_ParseTuple(args, "O", &self))
        return NULL;
    PGA_ctx = PyObject_GetAttrString (self, "context");
    PyArg_Parse (PGA_ctx, "i", &ctx);
    return Py_BuildValue ("i", PGAGetStringLength (ctx));
}

static PyObject *PGA_evaluate (PyObject *self0, PyObject *args)
{
    PyObject *self;
    int p, pop;

    if (!PyArg_ParseTuple(args, "Oii", &self, &p, &pop))
        return NULL;
    /* FIXME: should raise NotImplementedError */
    printf ("Ooops\n");
    fflush (stdout);
    return Py_BuildValue ("i", 4711);
}

static PyObject *PGA_get_allele (PyObject *self0, PyObject *args)
{
    PyObject *self;
    PyObject   *PGA_ctx;
    PGAContext *ctx;
    int p, pop, i, allele;

    if (!PyArg_ParseTuple(args, "Oiii", &self, &p, &pop, &i))
        return NULL;
    PGA_ctx = PyObject_GetAttrString (self, "context");
    PyArg_Parse (PGA_ctx, "i", &ctx);
    allele  = PGAGetBinaryAllele (ctx, p, pop, i);
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

static PyMethodDef Module_Methods[] = { {NULL, NULL, 0, NULL} };

PyMODINIT_FUNC initpga (void)
{
    PyMethodDef *def;
    PyObject *module      = Py_InitModule       ("pga", Module_Methods);
    PyObject *module_Dict = PyModule_GetDict    (module);
    PyObject *class_Dict  = PyDict_New          ();
    PyObject *class_Name  = PyString_FromString ("PGA");
    PyObject *pga_Class   = PyClass_New         (NULL, class_Dict, class_Name);
    PyDict_SetItemString (module_Dict, "PGA", pga_Class);
    context               = Py_BuildValue       ("{}");

    Py_DECREF(class_Dict);
    Py_DECREF(class_Name);
    Py_DECREF(pga_Class);
    
    /* add methods to class */
    for (def = PGA_Methods; def->ml_name != NULL; def++) {
        PyObject *func   = PyCFunction_New (def,  NULL);
        PyObject *method = PyMethod_New    (func, NULL, pga_Class);
        PyDict_SetItemString (class_Dict, def->ml_name, method);
        Py_DECREF (func);
        Py_DECREF (method);
    }
}
