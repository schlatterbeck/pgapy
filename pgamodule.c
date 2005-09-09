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

static PGAContext *get_context (PyObject *self)
{
    PyObject   *PGA_ctx = PyObject_GetAttrString (self, "context");
    PGAContext *ctx;
    if (!PGA_ctx)
        return NULL;
    if (!PyArg_Parse (PGA_ctx, "i", &ctx))
    {
        Py_DECREF   (PGA_ctx);
        return NULL;
    }
    Py_DECREF   (PGA_ctx);
    return ctx;
}

static PyObject *get_self (PGAContext *ctx)
{
    PyObject *self, *PGA_ctx = Py_BuildValue    ("i", (int) ctx);

    if (!PGA_ctx)
        return NULL;
    self = PyObject_GetItem (context, PGA_ctx);
    Py_DECREF (PGA_ctx);
    return self;
}

#define ERR_CHECK(x,r) do {      \
    if (!(x)) {                  \
        error_occurred = 1;    \
        return (r);     \
    }                          \
} while (0)

/*
 * Need a hash table of mapping ctx to PGA objects. Look up the
 * appropriate object and call its PGA_evaluate
 */
static double evaluate (PGAContext *ctx, int p, int pop)
{
    double retval;
    PyObject *self, *res1, *res2;
    int r;

    ERR_CHECK (!error_occurred, (double)0);
    self    = get_self (ctx);
    ERR_CHECK (self, (double)0);
    res1    = PyObject_CallMethod (self, "evaluate", "ii", p, pop);
    ERR_CHECK (res1, (double)0);
    res2    = PyNumber_Float      (res1);
    ERR_CHECK (res2, (double)0);
    r = PyArg_Parse (res2, "d", &retval);
    ERR_CHECK (r, (double)0);
    Py_DECREF (self);
    Py_DECREF (res1);
    Py_DECREF (res2);
    return retval;
}

static int check_stop (PGAContext *ctx)
{
    PyObject *self;
    ERR_CHECK (!error_occurred, PGA_TRUE);
    self = get_self (ctx);
    ERR_CHECK (self, PGA_TRUE);
    if (PyObject_HasAttrString (self, "stop_cond"))
    {
        PyObject *r = PyObject_CallMethod (self, "stop_cond", NULL);
        int retval, rr;
        ERR_CHECK (r, PGA_TRUE);
        rr = PyArg_Parse (r, "i", &retval);
        ERR_CHECK (rr, PGA_TRUE);
        return !!retval;
    }
    return PGACheckStoppingConditions (ctx);
}

/*
 * Used if the calling object has a mutation method.
 */
static int mutation (PGAContext *ctx, int p, int pop, double mr)
{
    PyObject *self, *r;
    int retval, rr;
    ERR_CHECK (!error_occurred, 0);
    self = get_self (ctx);
    ERR_CHECK (self, PGA_TRUE);
    r    = PyObject_CallMethod (self, "mutation", "iid", p, pop, mr);
    ERR_CHECK (r, PGA_TRUE);
    rr = PyArg_Parse (r, "i", &retval);
    ERR_CHECK (rr, PGA_TRUE);
    return retval;
}

static PyObject *PGA_init (PyObject *self0, PyObject *args, PyObject *kw)
{
    int argc = 0, max = 0, length = 0, pop_size = 0, pga_type = 0;
    int random_seed = 0, max_GA_iter = 0, max_no_change = 0;
    int num_replace = -1, pop_replace_type = -1;
    int max_similarity = 0;
    double mutation_prob = 0;
    PyObject *PGA_ctx;
    PyObject *self = NULL, *type = NULL, *maximize = NULL, *init = NULL;
    PyObject *init_percent = NULL, *stopping_rule_types = NULL;
    char *argv [] = {NULL, NULL};
    PGAContext *ctx;
    static char *kwlist[] =
        { "self"
        , "type"
        , "length"
        , "maximize"
        , "pop_size"
        , "init"
        , "init_percent"
        , "random_seed"
        , "max_GA_iter"
        , "max_no_change"
        , "max_similarity"
        , "mutation_prob"
        , "stopping_rule_types"
        , "num_replace"
        , "pop_replace_type"
        , NULL
        };

    if  (!PyArg_ParseTupleAndKeywords 
            ( args
            , kw
            , "OOi|OiOOiiiidOii"
            , kwlist
            , &self
            , &type
            , &length
            , &maximize
            , &pop_size
            , &init
            , &init_percent
            , &random_seed
            , &max_GA_iter
            , &max_no_change
            , &max_similarity
            , &mutation_prob
            , &stopping_rule_types
            , &num_replace
            , &pop_replace_type
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
    if (PyObject_HasAttrString (self, "mutation"))
    {
        PGASetUserFunction (ctx, PGA_USERFUNCTION_MUTATION, (void *)mutation);
    }
    if (random_seed)
    {
        PGASetRandomSeed   (ctx, random_seed);
    }
    if (mutation_prob)
    {
        PGASetMutationProb (ctx, mutation_prob);
    }
    if (num_replace >= 0)
    {
        PGASetNumReplaceValue (ctx, num_replace);
    }
    if (pop_replace_type >= 0)
    {
        if (  pop_replace_type != PGA_POPREPL_BEST
           && pop_replace_type != PGA_POPREPL_RANDOM_NOREP
           && pop_replace_type != PGA_POPREPL_RANDOM_REP
           )
        {
            PyErr_SetString (PyExc_ValueError, "invalid pop_replace_type");
            return NULL;
        }
        PGASetPopReplaceType (ctx, pop_replace_type);
    }

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
    if (stopping_rule_types)
    {
        int i, len = PySequence_Length (stopping_rule_types);
        if (len < 0)
        {
            return NULL;
        }
        for (i = 0; i < len; i++)
        {
            PyObject *x = PySequence_GetItem (stopping_rule_types, i);
            int val;
            if (!x)
                return NULL;
            if (!PyArg_Parse (x, "i", &val))
                return NULL;
            if (  val != PGA_STOP_MAXITER
               && val != PGA_STOP_NOCHANGE
               && val != PGA_STOP_TOOSIMILAR
               )
            {
                PyErr_SetString 
                    (PyExc_ValueError, "Invalid stoppping_rule_type");
                return NULL;
            }
            PGASetStoppingRuleType (ctx, val);
        }
    }
    if (init || init_percent)
    {
        int datatype = PGAGetDataType (ctx);
        int is_real  = (datatype == PGA_DATATYPE_REAL);
        int i, len;
        void *i_low, *i_high;
        PyObject *initvals = (init ? init : init_percent);

        if (datatype != PGA_DATATYPE_INTEGER && datatype != PGA_DATATYPE_REAL)
        {
            PyErr_SetString (PyExc_ValueError, "Init only for int/real");
            return NULL;
        }
        if (init && init_percent)
        {
            PyErr_SetString (PyExc_ValueError, "Only one of init/init_percent");
            return NULL;
        }
        if (init_percent && !is_real)
        {
            PyErr_SetString (PyExc_ValueError, "init_percent only for float");
            return NULL;
        }
        len = PySequence_Length (initvals);
        if (len < 0)
        {
            return NULL;
        }
        if (len != PGAGetStringLength (ctx))
        {
            PyErr_SetString (PyExc_ValueError, "Init length != string length");
            return NULL;
        }
        i_low  = malloc (len * (is_real ? sizeof (double) : sizeof (int)));
        i_high = malloc (len * (is_real ? sizeof (double) : sizeof (int)));
        assert (i_low && i_high);
        for (i = 0; i < len; i++)
        {
            PyObject *x = PySequence_GetItem (initvals, i);
            PyObject *low, *high;
            if (!x)
                return NULL;
            low  = PySequence_GetItem (x, 0);
            if (!low)
                return NULL;
            high = PySequence_GetItem (x, 1);
            if (!high)
                return NULL;
            if (is_real)
            {
                PyObject *l, *h;
                double hi;
                l = PyNumber_Float (low);
                if (!l)
                    return NULL;
                h = PyNumber_Float (high);
                if (!h)
                    return NULL;
                if (  !PyArg_Parse (h, "d", ((double *)i_high) + i)
                   || !PyArg_Parse (l, "d", ((double *)i_low)  + i)
                   )
                    return NULL;
                hi = ((double *)i_high) [i];
                if (init_percent && (hi <= 0 || hi > 1))
                {
                    PyErr_SetString 
                        (PyExc_ValueError, "Percentage must be 0 < p <= 1");
                    return NULL;
                }
            }
            else
            {
                PyObject *l, *h;
                l = PyNumber_Int (low);
                if (!l)
                    return NULL;
                h = PyNumber_Int (high);
                if (!h)
                    return NULL;
                if (  !PyArg_Parse (h, "i", ((int *)i_high) + i)
                   || !PyArg_Parse (l, "i", ((int *)i_low)  + i)
                   )
                    return NULL;
            }
        }
        if (is_real)
        {
            if (init)
            {
                PGASetRealInitRange    (ctx, i_low, i_high);
            }
            else
            {
                PGASetRealInitPercent  (ctx, i_low, i_high);
            }
        }
        else
        {
            PGASetIntegerInitRange (ctx, i_low, i_high);
        }
        free (i_low);
        free (i_high);
    }
    PGASetUp (ctx);
    PGA_ctx = Py_BuildValue ("i", (int) ctx);
    PyObject_SetItem       (context, PGA_ctx, self);
    PyObject_SetAttrString (self, "context", PGA_ctx);
    Py_DECREF (PGA_ctx);
    Py_INCREF (Py_None);
    return Py_None;
}

static PyObject *PGA_check_stopping_conditions (PyObject *self0, PyObject *args)
{
    PyObject *self;
    PGAContext *ctx;

    if (!PyArg_ParseTuple(args, "O", &self))
        return NULL;
    if (!(ctx = get_context (self)))
        return NULL;
    return Py_BuildValue ("i", PGACheckStoppingConditions (ctx));
}

static PyObject *PGA_del (PyObject *self0, PyObject *args)
{
    PyObject *self;
    PyObject   *PGA_ctx;
    PGAContext *ctx;

    if (!PyArg_ParseTuple(args, "O", &self))
        return NULL;
    Py_INCREF (Py_None);
    PGA_ctx = PyObject_GetAttrString (self, "context");
    /*
    fprintf (stderr, "After PGA_ctx in PGA_del: %08X\n", (int) PGA_ctx);
    fflush  (stderr);
    */
    if (!PGA_ctx)
        return Py_None;
    if (!PyArg_Parse (PGA_ctx, "i", &ctx))
        return Py_None;
    PyObject_DelItem     (context, PGA_ctx);
    Py_DECREF            (PGA_ctx);
    PGADestroy           (ctx);
    return Py_None;
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

static PyObject *PGA_len (PyObject *self0, PyObject *args)
{
    PyObject *self;
    PGAContext *ctx;

    if (!PyArg_ParseTuple(args, "O", &self))
        return NULL;
    if (!(ctx = get_context (self)))
        return NULL;
    return Py_BuildValue ("i", PGAGetStringLength (ctx));
}

static int check_allele (PGAContext *ctx, int p, int pop, int i)
{
    if (pop != PGA_OLDPOP && pop != PGA_NEWPOP)
    {
        char x [50];
        sprintf (x, "%d: invalid population", pop);
        PyErr_SetString (PyExc_ValueError, x);
        return 0;
    }
    if ((p < 0 || p >= PGAGetPopSize (ctx)) && p != PGA_TEMP1 && p != PGA_TEMP2)
    {
        char x [50];
        sprintf (x, "%d: invalid population index", p);
        PyErr_SetString (PyExc_ValueError, x);
        return 0;
    }
    if (i < 0 || i >= PGAGetStringLength (ctx))
    {
        char x [50];
        sprintf (x, "%d: allele index out of range", i);
        PyErr_SetString (PyExc_ValueError, x);
        return 0;
    }
    return 1;
}

/*
 * Get and Set methods.
 */
static PyObject *PGA_get_allele (PyObject *self0, PyObject *args)
{
    PyObject *self;
    PGAContext *ctx;
    int p, pop, i;

    if (!PyArg_ParseTuple(args, "Oiii", &self, &p, &pop, &i))
        return NULL;
    if (!(ctx = get_context (self)))
        return NULL;
    if (!check_allele (ctx, p, pop, i))
        return NULL;

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

static PyObject *PGA_set_allele (PyObject *self0, PyObject *args)
{
    PyObject *self, *val;
    PGAContext *ctx;
    int p, pop, i;

    if (!PyArg_ParseTuple(args, "OiiiO", &self, &p, &pop, &i, &val))
        return NULL;
    if (!(ctx = get_context (self)))
        return NULL;
    if (!check_allele (ctx, p, pop, i))
        return NULL;

    switch (PGAGetDataType (ctx)) {
    case PGA_DATATYPE_BINARY:
    {
        int allele;
        if (!PyArg_Parse (val, "i", &allele))
            return NULL;
        PGASetBinaryAllele (ctx, p, pop, i, allele);
        break;
    }
    case PGA_DATATYPE_CHARACTER:
    {
        char allele;
        if (!PyArg_Parse (val, "c", &allele))
            return NULL;
        PGASetCharacterAllele (ctx, p, pop, i, allele);
        break;
    }
    case PGA_DATATYPE_INTEGER:
    {
        int allele;
        if (!PyArg_Parse (val, "i", &allele))
            return NULL;
        PGASetIntegerAllele (ctx, p, pop, i, allele);
        break;
    }
    case PGA_DATATYPE_REAL:
    {
        double allele;
        if (!PyArg_Parse (val, "d", &allele))
            return NULL;
        PGASetRealAllele (ctx, p, pop, i, allele);
        break;
    }
    default:
        assert (0);
    }
    Py_INCREF   (Py_None);
    return Py_None;
}

static PyObject *PGA_get_best_index (PyObject *self0, PyObject *args)
{
    PyObject *self;
    PGAContext *ctx;
    int pop;

    if (!PyArg_ParseTuple(args, "Oi", &self, &pop))
        return NULL;
    if (!(ctx = get_context (self)))
        return NULL;
    return Py_BuildValue ("i", PGAGetBestIndex (ctx, pop));
}

static PyObject *PGA_set_random_seed (PyObject *self0, PyObject *args)
{
    PyObject *self;
    PGAContext *ctx;
    int seed;

    if (!PyArg_ParseTuple(args, "Oi", &self, &seed))
        return NULL;
    if (!(ctx = get_context (self)))
        return NULL;
    PGASetRandomSeed (ctx, seed);
    Py_INCREF   (Py_None);
    return Py_None;
}

static PyObject *PGA_random_01 (PyObject *self0, PyObject *args)
{
    PyObject   *self;
    PGAContext *ctx;

    if (!PyArg_ParseTuple(args, "O", &self))
        return NULL;
    if (!(ctx = get_context (self)))
        return NULL;
    return Py_BuildValue ("d", PGARandom01 (ctx, 0));
}

static int check_probability (double probability)
{
    if (probability < 0 || probability > 1)
    {
        PyErr_SetString 
            (PyExc_ValueError, "Probability must be 0 <= p <= 1");
        return 0;
    }
    return 1;
}

static PyObject *PGA_random_flip (PyObject *self0, PyObject *args)
{
    PyObject   *self;
    PGAContext *ctx;
    double     probability;

    if (!PyArg_ParseTuple(args, "Od", &self, &probability))
        return NULL;
    if (!check_probability (probability))
        return NULL;
    if (!(ctx = get_context (self)))
        return NULL;
    return Py_BuildValue ("i", PGARandomFlip (ctx, probability));
}

#define check_interval(l,r) do { \
    if (l > r)                                                             \
    {                                                                      \
        PyErr_SetString                                                    \
            (PyExc_ValueError, "interval_left must be <= interval_right"); \
        return NULL;                                                       \
    }                                                                      \
} while (0)

static PyObject *PGA_random_interval (PyObject *self0, PyObject *args)
{
    PyObject   *self;
    PGAContext *ctx;
    int        l, r;

    if (!PyArg_ParseTuple(args, "Oii", &self, &l, &r))
        return NULL;
    check_interval (l, r);
    if (!(ctx = get_context (self)))
        return NULL;
    return Py_BuildValue ("i", PGARandomInterval (ctx, l, r));
}

static PyObject *PGA_random_uniform (PyObject *self0, PyObject *args)
{
    PyObject   *self;
    PGAContext *ctx;
    double     l, r;

    if (!PyArg_ParseTuple(args, "Odd", &self, &l, &r))
        return NULL;
    check_interval (l, r);
    if (!(ctx = get_context (self)))
        return NULL;
    return Py_BuildValue ("d", PGARandomUniform (ctx, l, r));
}

static PyObject *PGA_random_gaussian (PyObject *self0, PyObject *args)
{
    PyObject   *self;
    PGAContext *ctx;
    double     l, r;

    if (!PyArg_ParseTuple(args, "Odd", &self, &l, &r))
        return NULL;
    if (!(ctx = get_context (self)))
        return NULL;
    return Py_BuildValue ("d", PGARandomGaussian (ctx, l, r));
}

static PyObject *PGA_run (PyObject *self0, PyObject *args)
{
    PyObject *self;
    PGAContext *ctx;

    if (!PyArg_ParseTuple(args, "O", &self))
        return NULL;
    if (!(ctx = get_context (self)))
        return NULL;
    PGARun      (ctx, evaluate);
    if (error_occurred)
    {
        return NULL;
    }
    Py_INCREF   (Py_None);
    return Py_None;
}

static PyMethodDef PGA_Methods [] =
{ { "__del__",                   PGA_del,                       METH_VARARGS
  , "Delete object and PGA data structures"
  }
, { "__init__",     (PyCFunction)PGA_init,       METH_VARARGS | METH_KEYWORDS
  , "Init object"
  }
, { "__len__",                   PGA_len,                       METH_VARARGS
  , "Return length of gene"
  }
, { "check_stopping_conditions", PGA_check_stopping_conditions, METH_VARARGS
  , "Return original stop condition check"
  }
, { "evaluate",                  PGA_evaluate,                  METH_VARARGS
  , "Evaluate"
  }
, { "get_allele",                PGA_get_allele,                METH_VARARGS
  , "Get allele"
  }
, { "get_best_index",            PGA_get_best_index,            METH_VARARGS
  , "Get best index in population pop"
  }
, { "random01",                  PGA_random_01,                 METH_VARARGS
  , "Random float 0 <= f <= 1"
  }
, { "random_flip",               PGA_random_flip,               METH_VARARGS
  , "Random int 0/1 with probability p"
  }
, { "random_interval",           PGA_random_interval,           METH_VARARGS
  , "Random int [l,r]"
  }
, { "random_uniform",            PGA_random_uniform,            METH_VARARGS
  , "Random float [l,r]"
  }
, { "random_gaussian",           PGA_random_gaussian,           METH_VARARGS
  , "Random value from gaussian distribution with mean, std_deviation"
  }
, { "run",                       PGA_run,                       METH_VARARGS
  , "Run optimization"
  }
, { "set_allele",                PGA_set_allele,                METH_VARARGS
  , "Set allele"
  }
, { "set_random_seed",           PGA_set_random_seed,           METH_VARARGS
  , "Set random seed to integer value"
  }
, {NULL, NULL, 0, NULL}
};

typedef struct
{
    char *cd_name;
    int   cd_value;
} constdef_t;

static constdef_t constdef [] =
    { {"PGA_NEWPOP",               PGA_NEWPOP                }
    , {"PGA_OLDPOP",               PGA_OLDPOP                }
    , {"PGA_POPREPL_BEST",         PGA_POPREPL_BEST          }
    , {"PGA_POPREPL_RANDOM_REP",   PGA_POPREPL_RANDOM_REP    }
    , {"PGA_POPREPL_RANDOM_NOREP", PGA_POPREPL_RANDOM_NOREP  }
    , {"PGA_STOP_NOCHANGE",        PGA_STOP_NOCHANGE         }
    , {"PGA_STOP_MAXITER",         PGA_STOP_MAXITER          }
    , {"PGA_STOP_TOOSIMILAR",      PGA_STOP_TOOSIMILAR       }
    , {NULL,                       0                         }
    };

static PyMethodDef Module_Methods[] = { {NULL, NULL, 0, NULL} };

PyMODINIT_FUNC initpga (void)
{
    PyMethodDef *def;
    constdef_t  *cd;
    PyObject *module      = Py_InitModule       ("pga", Module_Methods);
    PyObject *module_Dict = PyModule_GetDict    (module);
    PyObject *class_Dict  = PyDict_New          ();
    PyObject *class_Name  = PyString_FromString ("PGA");
    PyObject *pga_Class   = PyClass_New         (NULL, class_Dict, class_Name);
    context               = Py_BuildValue       ("{}");
    PyDict_SetItemString (module_Dict, "PGA",     pga_Class);
    PyDict_SetItemString (module_Dict, "context", context);

    for (cd = constdef; cd->cd_name; cd++) {
        PyObject *constant = Py_BuildValue    ("i", cd->cd_value);
        PyDict_SetItemString (module_Dict, cd->cd_name, constant);
        Py_DECREF (constant);
    }

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
