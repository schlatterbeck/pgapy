#include <Python.h>
#include <pgapack.h>
#include <stdio.h>
#undef NDEBUG
#include <assert.h>

static PyObject *context        = NULL;
static int       error_occurred = 0;

# if 0
static void prc (PyObject *o, char *str)
{
    fprintf (stderr, "%s: %08X", str, (int)o);
    if (o)
    {
        fprintf (stderr, "%s: Refcount: %d\n", str, o->ob_refcnt);
    }
    else
    {
        fprintf (stderr, "\n");
    }
    fflush  (stderr);
}
# endif

/*
 * Need a hash table of mapping ctx to PGA objects. Look up the
 * appropriate object and call its PGA_evaluate
 */
static double evaluate (PGAContext *ctx, int p, int pop)
{
    double retval;
    PyObject *PGA_ctx, *o, *res1, *res2;

    if (error_occurred)
    {
        return (double) 0;
    }
    PGA_ctx = Py_BuildValue       ("i", (int) ctx);
    assert (PGA_ctx);
    o       = PyObject_GetItem    (context, PGA_ctx);
    assert (o);
    res1    = PyObject_CallMethod (o, "evaluate", "ii", p, pop);
    if (!res1)
    {
        error_occurred = 1;
        return (double) 0;
    }
    res2    = PyNumber_Float      (res1);
    if (!res2)
    {
        error_occurred = 1;
        return (double) 0;
    }
    PyArg_Parse (res2, "d", &retval);
    Py_DECREF (o);
    Py_DECREF (PGA_ctx);
    Py_DECREF (res1);
    Py_DECREF (res2);
    return retval;
}

static int check_stop (PGAContext *ctx)
{
    if (error_occurred)
    {
        return PGA_TRUE;
    }
    return PGACheckStoppingConditions (ctx);
}

static PyObject *PGA_init (PyObject *self0, PyObject *args, PyObject *kw)
{
    int argc = 0, max = 0, length = 0, pop_size = 0, pga_type = 0;
    PyObject *PGA_ctx;
    PyObject *self = NULL, *type = NULL, *maximize = NULL;
    char *argv [] = {NULL, NULL};
    PGAContext *ctx;
    static char *kwlist[] =
        {"self", "type", "length", "maximize", "pop_size", NULL};

    if  (!PyArg_ParseTupleAndKeywords 
            ( args
            , kw
            , "OOi|Oi"
            , kwlist
            , &self
            , &type
            , &length
            , &maximize
            , &pop_size
            )
        )
    {
        return NULL;
    }

    if (PyObject_IsSubclass      (type, (PyObject *)&PyBool_Type))
    {
        pga_type = PGA_DATATYPE_BINARY;
    }
    else if (PyObject_IsSubclass (type, (PyObject *)&PyInt_Type))
    {
        pga_type = PGA_DATATYPE_INTEGER;
    }
    else if (PyObject_IsSubclass (type, (PyObject *)&PyFloat_Type))
    {
        pga_type = PGA_DATATYPE_REAL;
    }
    else if (PyObject_IsSubclass (type, (PyObject *)&PyString_Type))
    {
        pga_type = PGA_DATATYPE_CHARACTER;
    }
    else
    {
        /* FIXME: Implement PGA_DATATYPE_USER */
        PyErr_SetString \
            ( PyExc_NotImplementedError
            , "Gene type must currently be one of [bool, int, real, string]"
            );
        return NULL;
    }

    if (maximize)
    {
        max = PyObject_IsTrue (maximize);
    }
    if (length <= 0)
    {
        PyErr_SetString \
            ( PyExc_ValueError
            , "Gene length must be at least 1"
            );
        return NULL;
    }
    /*
     * Get class name (FIXME -- needed if debug options should be
     * supported)
     */
    argv [0] = "huhu";
    
    ctx = PGACreate
        ( &argc
        , argv
        , pga_type
        , length
        , max ? PGA_MAXIMIZE : PGA_MINIMIZE
        );
    PGASetUserFunction (ctx, PGA_USERFUNCTION_STOPCOND, (void *)check_stop);
    if (pop_size)
    {
        if (pop_size <= 1)
        {
            PyErr_SetString \
                ( PyExc_ValueError
                , "Population size must be at least 2"
                );
            return NULL;
        }
        PGASetPopSize (ctx, pop_size);
    }
    PGASetUp (ctx);
    PGA_ctx = Py_BuildValue ("i", (int) ctx);
    PyObject_SetItem       (context, PGA_ctx, self);
    PyObject_SetAttrString (self, "context", PGA_ctx);
    Py_DECREF (PGA_ctx);
    Py_INCREF (Py_None);
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
    PyArg_Parse (PGA_ctx, "i", &ctx);
    Py_DECREF   (PGA_ctx);
    PGARun      (ctx, evaluate);
    if (error_occurred)
    {
        return NULL;
    }
    Py_INCREF   (Py_None);
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
    PyArg_Parse          (PGA_ctx, "i", &ctx);
    Py_DECREF            (PGA_ctx);
    return Py_BuildValue ("i", PGAGetStringLength (ctx));
}

static PyObject *PGA_evaluate (PyObject *self0, PyObject *args)
{
    PyObject *self;
    int p, pop;

    if (!PyArg_ParseTuple(args, "Oii", &self, &p, &pop))
        return NULL;
    PyErr_SetString \
        ( PyExc_NotImplementedError
        , "You must define \"evaluate\" in a derived class"
        );
    return NULL;
}

static PyObject *PGA_get_allele (PyObject *self0, PyObject *args)
{
    PyObject *self;
    PyObject   *PGA_ctx;
    PGAContext *ctx;
    int p, pop, i;

    if (!PyArg_ParseTuple(args, "Oiii", &self, &p, &pop, &i))
        return NULL;
    PGA_ctx = PyObject_GetAttrString (self, "context");
    PyArg_Parse                      (PGA_ctx, "i", &ctx);
    Py_DECREF                        (PGA_ctx);
    PGA_ctx = NULL;

    if (pop != PGA_OLDPOP && pop != PGA_NEWPOP)
    {
        char x [50];
        sprintf (x, "%d: invalid population", pop);
        PyErr_SetString (PyExc_ValueError, x);
        return NULL;
    }
    if (p < 0 || p >= PGAGetStringLength (ctx))
    {
        char x [50];
        sprintf (x, "%d: invalid index", p);
        PyErr_SetString (PyExc_ValueError, x);
        return NULL;
    }
    switch (PGAGetDataType (ctx)) {
    case PGA_DATATYPE_BINARY:
    {
        int allele = PGAGetBinaryAllele (ctx, p, pop, i);
        return Py_BuildValue ("i", allele);
        break;
    }
    case PGA_DATATYPE_CHARACTER:
    {
        char allele = PGAGetCharacterAllele (ctx, p, pop, i);
        return Py_BuildValue ("c", allele);
        break;
    }
    case PGA_DATATYPE_INTEGER:
    {
        int allele = PGAGetIntegerAllele (ctx, p, pop, i);
        return Py_BuildValue ("i", allele);
        break;
    }
    case PGA_DATATYPE_REAL:
    {
        double allele = PGAGetRealAllele (ctx, p, pop, i);
        return Py_BuildValue ("d", allele);
        break;
    }
    default:
        assert (0);
    }
    Py_INCREF   (Py_None);
    return Py_None;
}


/*
__del__
    PGADestroy (ctx)
*/

static PyMethodDef PGA_Methods [] =
{ { "__init__",   (PyCFunction)PGA_init, METH_VARARGS | METH_KEYWORDS
  , "Init object"
  }
, { "__len__",    PGA_len,               METH_VARARGS
  , "Return length of gene"
  }
, { "evaluate",   PGA_evaluate,          METH_VARARGS
  , "Evaluate"
  }
, { "get_allele", PGA_get_allele,        METH_VARARGS
  , "Get allele"
  }
, { "run",        PGA_run,               METH_VARARGS
  , "Run optimization"
  }
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
    context               = Py_BuildValue       ("{}");
    PyDict_SetItemString (module_Dict, "PGA",     pga_Class);
    PyDict_SetItemString (module_Dict, "context", context);

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
