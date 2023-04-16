/* Copyright (C) 2005-23 Dr. Ralf Schlatterbeck Open Source Consulting.
 * Reichergasse 131, A-3411 Weidling.
 * Web: http://www.runtux.com Email: office@runtux.com
 * ****************************************************************************
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
 * IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 * PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
 * TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * ****************************************************************************
 */

#include <Python.h>
#include <pgapack.h>
#include <stdio.h>
#undef NDEBUG
#include <stddef.h>
#include <assert.h>
#include <Version.h>

#define IS_PY3 (PY_VERSION_HEX >= 0x3000000)

/* Conditional compilation for python2 vs python3 */
# if IS_PY3
# define PyInt_Type_Compat PyLong_Type
# define PyString_Check_Compat PyUnicode_Check
# define PyString_FromString_Compat(x) PyUnicode_FromString((x))
# define PyFileObject PyObject
# define FAIL NULL
# else /* Python 2 */
# define FAIL
# define PyString_Check_Compat PyString_Check
# define PyString_FromString_Compat(x) PyString_FromString((x))
# define PyInt_Type_Compat PyInt_Type
# endif /* Python 2 */

/* Necessary forward declarations */
static char **parse_argv (PyObject *argv, int *argcp);

/*
 * Global data
 */
/* This is a dictionary for retrieving Python PGA objects by PGA ctx */
static PyObject *contexts       = NULL;

/* Error handling macros */
#define SET_ERR(ctx) (*((int *)ctx->ga.CustomData) = 1)
#define HAS_ERR(ctx) (*((int *)(ctx)->ga.CustomData))
#define ERR_CHECK(ctx,x,r) do {             \
    if (!(x)) {                             \
        SET_ERR(ctx);                       \
        return r;                           \
    }                                       \
} while (0)
#define ERR_CHECK_ERRNO(ctx,x,r) do {                            \
    if (!(x)) {                                                  \
        PyErr_SetString (PyExc_RuntimeError, strerror (errno));  \
        SET_ERR(ctx);                                            \
        return r;                                                \
    }                                                            \
} while (0)

#define ERR_CHECK_OCCURRED(ctx,r) do {                                 \
    if (HAS_ERR(ctx)) {                                                \
        if (!PyErr_Occurred ()) {                                      \
            PyErr_SetString                                            \
                (PyExc_RuntimeError, "PGA object is in status error"); \
        }                                                              \
        return r;                                                      \
    }                                                                  \
} while (0)

#define ERR_CHECK_X(ctx,x) do {             \
    if (!(x)) {                             \
        SET_ERR(ctx);                       \
        goto errout;                        \
    }                                       \
} while (0)
#define ERR_CHECK_X_OCCURRED(ctx) do {      \
    if (HAS_ERR(ctx)) {                     \
        goto errout;                        \
    }                                       \
} while (0)

#define ERR_CHECK_RET(x) do {  \
    if (!(x)) {                \
        goto errout;           \
    }                          \
} while (0)

#define ERR_CHECK_VALUE_ERROR(x,e) do {          \
    if (!(x)) {                                  \
        PyErr_SetString (PyExc_ValueError, (e)); \
        goto errout;                             \
    }                                            \
} while (0)

#define CHECK_VALUE(cond, errmsg) \
    CHECK_VALUE_EXCEPTION ((cond), (errmsg), PyExc_ValueError, INIT_FAIL)
#define CHECK_VALUE_EXCEPTION(cond, errmsg, exc, ret) do {  \
    if (!(cond)) {                                          \
        PyErr_SetString ((exc), (errmsg));                  \
        return (ret);                                       \
    }                                                       \
} while (0)
#define CHECK_VALUE_AND_FREE(cond, errmsg, exc, tofree) do { \
    if (!(cond)) {                                           \
        free ((tofree));                                     \
        PyErr_SetString ((exc), (errmsg));                   \
        return INIT_FAIL;                                    \
    }                                                        \
} while (0)

#define ERR_RET(cond,ret) do {    \
    if (!(cond)) {                \
        return (ret);             \
    }                             \
} while (0)

#define ERR_DECREF_RET(cond, dec, ret) do {  \
    if (!(cond)) {                           \
        Py_CLEAR (dec);                      \
        return (ret);                        \
    }                                        \
} while (0)
#define ERR_DECREF_2_RET(cond, dec1, dec2, ret) do {  \
    if (!(cond)) {                                    \
        Py_CLEAR (dec1);                              \
        Py_CLEAR (dec2);                              \
        return (ret);                                 \
    }                                                 \
} while (0)
#define ERR_DECREF_3_RET(cond, dec1, dec2, dec3, ret) do {  \
    if (!(cond)) {                                          \
        Py_CLEAR (dec1);                                    \
        Py_CLEAR (dec2);                                    \
        Py_CLEAR (dec3);                                    \
        return (ret);                                       \
    }                                                       \
} while (0)

/*********************************
 * Convenience functions in module
 *********************************/
/* Visual C disable warning about unused variable "self" */
#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable:4100)
#endif

static PyObject *das_dennis (PyObject *self, PyObject *args, PyObject *kw)
#ifdef _MSC_VER
#pragma warning(pop)
#endif
{
    int i, j;
    int dim, npart;
    double scale = 1;
    double *dir = NULL;
    PyObject *direction = NULL;
    PyObject *tuple = NULL;
    PyObject *inner_tuple = NULL;
    PyObject *ele = NULL;
    void *result = NULL;
    PyObject *res = NULL, *res2 = NULL;
    int npoints;
    static char *kwlist [] =
        { "dimension"
        , "npartitions"
        , "scale"
        , "direction"
        , NULL
        };
    if (!PyArg_ParseTupleAndKeywords
        (args, kw, "ii|dO", kwlist, &dim, &npart, &scale, &direction)
       )
    {
        return NULL;
    }
    ERR_CHECK_VALUE_ERROR (dim >= 1, "Dimension must be >= 1");
    ERR_CHECK_VALUE_ERROR (npart >= 1, "npartitions must be >= 1");
    if (direction != NULL) {
        Py_ssize_t length;
        ERR_CHECK_VALUE_ERROR \
            ( PySequence_Check (direction)
            , "Expected sequence for direction parameter"
            );
        length = PySequence_Length (direction);
        ERR_CHECK_VALUE_ERROR \
            ( length == dim
            , "Direction must be sequence with length=dimension"
            );
        if ((dir = malloc (sizeof (double) * dim)) == NULL) {
            return PyErr_NoMemory ();
        }
        for (i=0; i<length; i++) {
            res = PySequence_GetItem (direction, i);
            ERR_CHECK_RET (res != NULL);
            res2 = PyNumber_Float (res);
            ERR_CHECK_RET (res2 != NULL);
            Py_CLEAR (res);
            dir [i] = PyFloat_AsDouble (res2);
            ERR_CHECK_RET (!PyErr_Occurred ());
            Py_CLEAR (res2);
        }
    }
    npoints = LIN_binom (dim + npart - 1, npart);
    LIN_dasdennis (dim, npart, &result, 0, scale, dir);
    tuple = PyTuple_New (npoints);
    ERR_CHECK_RET (tuple != NULL);
    PyTuple_SetItem (tuple, 0, ele);
    for (i=0; i<npoints; i++) {
        int r;
        double (*points)[3] = (double (*)[3]) result;
        inner_tuple = PyTuple_New (dim);
        ERR_CHECK_RET (inner_tuple != NULL);
        r = PyTuple_SetItem (tuple, i, inner_tuple);
        assert (r == 0);
        for (j=0; j<dim; j++) {
            ele = Py_BuildValue ("d", points [i][j]);
            ERR_CHECK_RET (ele != NULL);
            r = PyTuple_SetItem (inner_tuple, j, ele);
            assert (r == 0);
        }
    }
    free (result);
    return tuple;
errout:
    if (dir != NULL) {
        free (dir);
    }
    Py_CLEAR (tuple);
    Py_CLEAR (res);
    Py_CLEAR (res2);
    return NULL;
}

/***********************
 * Wrapped MPI functions
 ***********************/

/* These are needed whenever we want to perform *several* PGA
 * optimization runs: In that case we need to do MPI_Init first (and
 * finally an MPI_finit) because MPI_Init may be called only once.
 * These are *module* methods (not PGA object methods)
 */

static PyObject *PGA_MPI_Init (PyObject *self, PyObject *args, PyObject *kw)
{
    int argc;
    int ret = 0;
    char **c_argv = NULL;
    PyObject *pyargv;
    static char *kwlist [] = {"argv", NULL};

    if (!PyArg_ParseTupleAndKeywords (args, kw, "O", kwlist, &pyargv)) {
        return NULL;
    }
    c_argv = parse_argv (pyargv, &argc);
    Py_DECREF (pyargv);
    ret = MPI_Init (&argc, &c_argv);
    return Py_BuildValue ("i", ret);
}

static PyObject *PGA_MPI_Abort (PyObject *self, PyObject *args, PyObject *kw)
{
    int errorcode = 0;
    int ret = 0;
    static char *kwlist [] =
        { "errorcode"
        , NULL
        };
    if (!PyArg_ParseTupleAndKeywords (args, kw, "i", kwlist, &errorcode)) {
        return NULL;
    }
    ret = MPI_Abort (MPI_COMM_WORLD, errorcode);
    return Py_BuildValue ("i", ret);
}

/* We don't care if someone calls this with arguments */
static PyObject *PGA_MPI_Finalize (PyObject *self, PyObject *args, PyObject *kw)
{
    int ret = MPI_Finalize ();
    return Py_BuildValue ("i", ret);
}

static PyMethodDef Module_Methods [] =
{ { "das_dennis", (PyCFunction)das_dennis, METH_VARARGS | METH_KEYWORDS
  , "Return Das/Dennis points"
  }
, { "MPI_Abort", (PyCFunction)PGA_MPI_Abort, METH_VARARGS | METH_KEYWORDS
  , "Abort MPI"
  }
, { "MPI_Finalize", (PyCFunction)PGA_MPI_Finalize, METH_VARARGS | METH_KEYWORDS
  , "Finalize MPI"
  }
, { "MPI_Init", (PyCFunction)PGA_MPI_Init, METH_VARARGS | METH_KEYWORDS
  , "Initialize MPI"
  }
, { NULL } /* EMPTY VALUE AS END-MARKER */
};


/***********************
 * Constants (in Module)
 ***********************/

typedef struct
{
    char *cd_name;
    int   cd_value;
} constdef_t;

/* These need to be kept sorted */
static constdef_t constdef [] =
    { {"PGA_CINIT_LOWER",           PGA_CINIT_LOWER           }
    , {"PGA_CINIT_MIXED",           PGA_CINIT_MIXED           }
    , {"PGA_CINIT_UPPER",           PGA_CINIT_UPPER           }
    , {"PGA_CROSSOVER_EDGE",        PGA_CROSSOVER_EDGE        }
    , {"PGA_CROSSOVER_ONEPT",       PGA_CROSSOVER_ONEPT       }
    , {"PGA_CROSSOVER_SBX",         PGA_CROSSOVER_SBX         }
    , {"PGA_CROSSOVER_TWOPT",       PGA_CROSSOVER_TWOPT       }
    , {"PGA_CROSSOVER_UNIFORM",     PGA_CROSSOVER_UNIFORM     }
    , {"PGA_DE_CROSSOVER_BIN",      PGA_DE_CROSSOVER_BIN      }
    , {"PGA_DE_CROSSOVER_EXP",      PGA_DE_CROSSOVER_EXP      }
    , {"PGA_DE_VARIANT_BEST",       PGA_DE_VARIANT_BEST       }
    , {"PGA_DE_VARIANT_EITHER_OR",  PGA_DE_VARIANT_EITHER_OR  }
    , {"PGA_DE_VARIANT_RAND",       PGA_DE_VARIANT_RAND       }
    , {"PGA_FITNESSMIN_CMAX",       PGA_FITNESSMIN_CMAX       }
    , {"PGA_FITNESSMIN_RECIPROCAL", PGA_FITNESSMIN_RECIPROCAL }
    , {"PGA_FITNESS_NORMAL",        PGA_FITNESS_NORMAL        }
    , {"PGA_FITNESS_RANKING",       PGA_FITNESS_RANKING       }
    , {"PGA_FITNESS_RAW",           PGA_FITNESS_RAW           }
    , {"PGA_MIX_MUTATE_AND_CROSS",  PGA_MIX_MUTATE_AND_CROSS  }
    , {"PGA_MIX_MUTATE_ONLY",       PGA_MIX_MUTATE_ONLY       }
    , {"PGA_MIX_MUTATE_OR_CROSS",   PGA_MIX_MUTATE_OR_CROSS   }
    , {"PGA_MIX_TRADITIONAL",       PGA_MIX_TRADITIONAL       }
    , {"PGA_MUTATION_CONSTANT",     PGA_MUTATION_CONSTANT     }
    , {"PGA_MUTATION_DE",           PGA_MUTATION_DE           }
    , {"PGA_MUTATION_GAUSSIAN",     PGA_MUTATION_GAUSSIAN     }
    , {"PGA_MUTATION_PERMUTE",      PGA_MUTATION_PERMUTE      }
    , {"PGA_MUTATION_POLY",         PGA_MUTATION_POLY         }
    , {"PGA_MUTATION_RANGE",        PGA_MUTATION_RANGE        }
    , {"PGA_MUTATION_UNIFORM",      PGA_MUTATION_UNIFORM      }
    , {"PGA_NEWPOP",                PGA_NEWPOP                }
    , {"PGA_OLDPOP",                PGA_OLDPOP                }
    , {"PGA_POPREPL_BEST",          PGA_POPREPL_BEST          }
    , {"PGA_POPREPL_NSGA_II",       PGA_POPREPL_NSGA_II       }
    , {"PGA_POPREPL_NSGA_III",      PGA_POPREPL_NSGA_III      }
    , {"PGA_POPREPL_PAIRWISE_BEST", PGA_POPREPL_PAIRWISE_BEST }
    , {"PGA_POPREPL_RANDOM_NOREP",  PGA_POPREPL_RANDOM_NOREP  }
    , {"PGA_POPREPL_RANDOM_REP",    PGA_POPREPL_RANDOM_REP    }
    , {"PGA_POPREPL_RTR",           PGA_POPREPL_RTR           }
    , {"PGA_REPORT_AVERAGE",        PGA_REPORT_AVERAGE        }
    , {"PGA_REPORT_HAMMING",        PGA_REPORT_HAMMING        }
    , {"PGA_REPORT_OFFLINE",        PGA_REPORT_OFFLINE        }
    , {"PGA_REPORT_ONLINE",         PGA_REPORT_ONLINE         }
    , {"PGA_REPORT_STRING",         PGA_REPORT_STRING         }
    , {"PGA_REPORT_WORST",          PGA_REPORT_WORST          }
    , {"PGA_SELECT_LINEAR",         PGA_SELECT_LINEAR         }
    , {"PGA_SELECT_PROPORTIONAL",   PGA_SELECT_PROPORTIONAL   }
    , {"PGA_SELECT_PTOURNAMENT",    PGA_SELECT_PTOURNAMENT    }
    , {"PGA_SELECT_SUS",            PGA_SELECT_SUS            }
    , {"PGA_SELECT_TOURNAMENT",     PGA_SELECT_TOURNAMENT     }
    , {"PGA_SELECT_TRUNCATION",     PGA_SELECT_TRUNCATION     }
    , {"PGA_STOP_MAXITER",          PGA_STOP_MAXITER          }
    , {"PGA_STOP_NOCHANGE",         PGA_STOP_NOCHANGE         }
    , {"PGA_STOP_TOOSIMILAR",       PGA_STOP_TOOSIMILAR       }
    , {NULL,                        0                         }
    };

int compare_constdef (const void *v1, const void *v2)
{
    const constdef_t *p1 = v1, *p2 = v2;
    return strcmp (p1->cd_name, p2->cd_name);
}

/***********************
 * Exported Symbols
 ***********************/

/*
 * Compute __all__, all exported symbols
 * Return 0 on success, -1 on error
 */
static int all_symbols (PyObject *module)
{
    PyObject *tuple = NULL;
    PyObject *ele = NULL;
    int nfun = sizeof (Module_Methods) / sizeof (PyMethodDef) - 1;
    int nobj = 2; /* the PGA class and VERSION */
    int ncon = sizeof (constdef) / sizeof (constdef_t) - 1;
    int idx  = 2;
    int r    = -1;
    constdef_t  *cd;
    PyMethodDef *funp = NULL;

    tuple = PyTuple_New (nfun + nobj + ncon);
    ERR_RET (tuple != NULL, -1);

    ele = Py_BuildValue ("s", "PGA");
    ERR_CHECK_RET (ele != NULL);
    r = PyTuple_SetItem (tuple, 0, ele);
    ERR_CHECK_RET (r == 0);
    ele = Py_BuildValue ("s", "VERSION");
    ERR_CHECK_RET (ele != NULL);
    r = PyTuple_SetItem (tuple, 1, ele);
    ERR_CHECK_RET (r == 0);

    /* Add constants */
    for (cd = constdef; cd->cd_name; cd++) {
        ele = Py_BuildValue ("s", cd->cd_name);
        ERR_CHECK_RET (ele != NULL);
        r = PyTuple_SetItem (tuple, idx++, ele);
        ERR_CHECK_RET (r == 0);
    }
    for (funp = Module_Methods; funp->ml_name; funp++) {
        ele = Py_BuildValue ("s", funp->ml_name);
        ERR_CHECK_RET (ele != NULL);
        r = PyTuple_SetItem (tuple, idx++, ele);
        ERR_CHECK_RET (r == 0);
    }
    ele = NULL;
    r = PyModule_AddObject (module, "__all__", tuple);
    ERR_CHECK_RET (r == 0);
    return 0;
errout:
    Py_CLEAR (ele);
    Py_CLEAR (tuple);
    return -1;
}

/******************************************************
 * Helpers for retriving ctx from self or self from ctx
 ******************************************************/

/*
 * Retrieve the PGApack ctx from a PGA object
 */
static PGAContext *get_context (PyObject *self)
{
    PyObject   *PGA_ctx = NULL;
    long long llctx;

    PGA_ctx = PyObject_GetAttrString (self, "context");
    if (!PGA_ctx) {
        return NULL;
    }
    ERR_DECREF_RET (PyArg_Parse (PGA_ctx, "L", &llctx), PGA_ctx, NULL);
    Py_DECREF (PGA_ctx);
    /* If an error occurred */
    if (*((int *)((PGAContext *)llctx)->ga.CustomData)) {
        return NULL;
    }
    /* Visual C disable warning about size */
    #ifdef _MSC_VER
    #pragma warning(push)
    #pragma warning(disable:4305)
    #endif
    return (PGAContext *)llctx;
    #ifdef _MSC_VER
    #pragma warning(pop)
    #endif
}

/*
 * Retrieve the FILE *fp from the PGA object
 */
static FILE *get_fp (PyObject *self, PyObject *file)
{
    PyObject *pyfp = PyObject_GetAttrString (self, "_fp");
    PyObject *fno = NULL;
    long long lfp;
    int fd_from_self = -1;
    int fd_from_file = -1;
    if (!pyfp) {
        return NULL;
    }
    ERR_DECREF_RET (PyArg_Parse (pyfp, "L", &lfp), pyfp, NULL);
    Py_DECREF (pyfp);
    fd_from_self = fileno ((FILE *)lfp);
    /* Now extract fd from file */
    fno = PyObject_CallMethod (file, "fileno", "");
    /* If an error happens here (should not) we return the fp anyway */
    if (fno != NULL) {
        if (PyArg_Parse (fno, "i", &fd_from_file)) {
            ERR_DECREF_RET (fd_from_self == fd_from_file, fno, NULL);
        }
        Py_DECREF (fno);
    }
    /* Visual C disable warning about size */
    #ifdef _MSC_VER
    #pragma warning(push)
    #pragma warning(disable:4305)
    #endif
    return (FILE *)lfp;
    #ifdef _MSC_VER
    #pragma warning(pop)
    #endif
}

static PyObject *module_os = NULL;

/* Return a FILE created from dup'd fd from file, caller needs to close */
static FILE *get_fp_from_file (PyObject *file)
{
    FILE *fp = NULL;
    int fd = -1;
    PyObject *fno = PyObject_CallMethod (file, "fileno", "");
    PyObject *fdp = NULL;
    if (fno == NULL) {
        return NULL;
    }
    if (module_os == NULL) {
        module_os = PyImport_ImportModule ("os");
        if (!module_os) {
            return NULL;
        }
    }
    fdp = PyObject_CallMethod (module_os, "dup", "O", fno);
    Py_DECREF (fno);
    if (fdp == NULL) {
        return NULL;
    }
    ERR_DECREF_RET (PyArg_Parse (fdp, "i", &fd), fdp, NULL);
    Py_DECREF (fdp);
    fp = fdopen (fd, "w");
    CHECK_VALUE_EXCEPTION \
        (fp != NULL, strerror (errno), PyExc_RuntimeError, NULL);
    return fp;
}

/*
 * Retrieve the PGA object via the PGApack ctx
 */
static PyObject *get_self (PGAContext *ctx)
{
    PyObject *self, *PGA_ctx = Py_BuildValue ("L", (long long) ctx);

    if (!PGA_ctx) {
        return NULL;
    }
    self = PyObject_GetItem (contexts, PGA_ctx);
    Py_DECREF (PGA_ctx);
    return self;
}

/**************************************************
 * PGApack callback functions
 * These get a PGApack ctx and retrieve the object.
 * Then they call into an object method
 **************************************************/

/*
 * Used if the calling object has an endofgen method.
 * Otherwise use built-in default for the datatype.
 * Called after each generation.
 * Do something useful, display the population on a graphics output,
 * let the user adjust the population, etc.
 */
static void endofgen (PGAContext *ctx)
{
    PyObject *self = NULL, *r = NULL;
    ERR_CHECK_X_OCCURRED (ctx);
    self = get_self (ctx);
    ERR_CHECK_X (ctx, self);
    r    = PyObject_CallMethod (self, "endofgen", "");
    ERR_CHECK_X (ctx, r);
errout:
    Py_CLEAR (r);
    Py_CLEAR (self);
    return;
}

/*
 * Need a hash table of mapping ctx to PGA objects. Look up the
 * appropriate object and call its PGA_evaluate
 */
static double evaluate (PGAContext *ctx, int p, int pop, double *aux)
{
    double retval = 0.0;
    PyObject *self = NULL, *res1 = NULL, *res2 = NULL, *res3 = NULL;
    int r;
    Py_ssize_t length, i;

    ERR_CHECK_X_OCCURRED (ctx);
    self    = get_self (ctx);
    ERR_CHECK_X (ctx, self);
    res1    = PyObject_CallMethod (self, "evaluate", "ii", p, pop);
    ERR_CHECK_X (ctx, res1);
    if (PySequence_Check (res1)) {
        length = PySequence_Length (res1);
        if (length != ctx->ga.NumAuxEval + 1) {
            char x [60];
            sprintf
                ( x, "Invalid length %zd of evaluations, expect %d"
                , length, ctx->ga.NumAuxEval + 1
                );
            PyErr_SetString (PyExc_ValueError, x);
            SET_ERR (ctx);
            goto errout;
        }
        res2 = PySequence_GetItem(res1, 0);
        ERR_CHECK_X (ctx, res2);
        res3 = PyNumber_Float (res2);
        Py_CLEAR (res2);
        ERR_CHECK_X (ctx, res3);
        r = PyArg_Parse (res3, "d", &retval);
        ERR_CHECK_X (ctx, r);
        Py_CLEAR (res3);
        for (i=1; i<length; i++) {
            res2 = PySequence_GetItem (res1, i);
            ERR_CHECK_X (ctx, res2);
            res3 = PyNumber_Float (res2);
            Py_CLEAR (res2);
            ERR_CHECK_X (ctx, res3);
            r = PyArg_Parse (res3, "d", aux + (i - 1));
            ERR_CHECK_X (ctx, r);
            Py_CLEAR (res3);
        }
    } else {
        if (ctx->ga.NumAuxEval) {
            char x [50];
            sprintf (x, "Expected %d evaluations", ctx->ga.NumAuxEval + 1);
            PyErr_SetString (PyExc_ValueError, x);
            SET_ERR (ctx);
            goto errout;
        }
        res2 = PyNumber_Float (res1);
        ERR_CHECK_X (ctx, res2);
        r = PyArg_Parse (res2, "d", &retval);
        ERR_CHECK_X (ctx, r);
    }
errout:
    Py_CLEAR (self);
    Py_CLEAR (res1);
    Py_CLEAR (res2);
    Py_CLEAR (res3);
    return retval;
}

/*
 * Hash function for duplicate checking
 * This is used if we have user-defined datatypes
 * The user datatype must be hashable (e.g. define a __hash__ method)
 */
static PGAHash build_hash (PGAContext *ctx, int p, int pop)
{
    Py_hash_t hash = 0;
    PyObject *self = NULL;
    PGAIndividual *ind = NULL;
    ERR_CHECK_X_OCCURRED (ctx);
    self = get_self (ctx);
    ERR_CHECK_X (ctx, self);
    if (PyObject_HasAttrString (self, "hash")) {
        int rr;
        PyObject *r = PyObject_CallMethod (self, "hash", "ii", p, pop);
        ERR_CHECK_X (ctx, r);
        rr = PyArg_Parse (r, "L", &hash);
        ERR_CHECK_X (ctx, rr);
        hash = ((hash >> 32) ^ hash) & 0xFFFFFFFF;
        return hash;
    }
    /* No custom hash method */

    /* We checked before installing the user function if either a custom
     * method is defined or we have user datatypes, so this assert
     * should not trigger here.
     */
    assert (ctx->ga.datatype == PGA_DATATYPE_USER);
    ind = PGAGetIndividual (ctx, p, pop);
    if (ind->chrom == NULL) {
        PyErr_SetString (PyExc_ValueError, "This gene is not set");
        goto errout;
    }
    hash = PyObject_Hash ((PyObject *)ind->chrom);
    ERR_CHECK_X (ctx, hash != -1);
    Py_CLEAR (self);
    hash = ((hash >> 32) ^ hash) & 0xFFFFFFFF;
    return (PGAHash)hash;
errout:
    Py_CLEAR (self);
    return 0;
}

/*
 * Used if the calling object has a check_duplicate method.
 * Otherwise use built-in default for the datatype.
 * Perform duplicate checking, compare strings p1 and p2 to see if they
 * are different. If they are, return non-zero, else return 0.
 */
static int check_duplicate (PGAContext *ctx, int p1, int pop1, int p2, int pop2)
{
    PyObject *self = NULL, *r = NULL;
    int rr, retval = 0;
    ERR_CHECK_X_OCCURRED (ctx);
    self = get_self (ctx);
    ERR_CHECK_X (ctx, self);
    r    = PyObject_CallMethod
        (self, "check_duplicate", "iiii", p1, pop1, p2, pop2);
    ERR_CHECK_X (ctx, r);
    rr = PyArg_Parse (r, "i", &retval);
    ERR_CHECK_X (ctx, rr);
errout:
    Py_CLEAR (r);
    Py_CLEAR (self);
    return !!retval;
}

/*
 * Check stopping criteria, this is always active.
 * User can set a stop_cond method to add stopping criteria.
 * We perform the magic here that if during one of the callback
 * functions (calling into python, e.g. evaluate) an error occurs, we
 * check the error flag here and stop. This way we can return an error
 * in one of the callback functions to python and raise the appropriate
 * exception there.
 */
static int check_stop (PGAContext *ctx)
{
    PyObject *self = NULL;
    ERR_CHECK_OCCURRED (ctx, PGA_TRUE);
    self = get_self (ctx);
    ERR_CHECK (ctx, self, PGA_TRUE);
    if (PyObject_HasAttrString (self, "stop_cond")) {
        int retval = PGA_TRUE, rr;
        PyObject *r = PyObject_CallMethod (self, "stop_cond", NULL);
        ERR_CHECK_X (ctx, r);
        rr = PyArg_Parse (r, "i", &retval);
        ERR_CHECK_X (ctx, rr);
    errout:
        Py_CLEAR (r);
        Py_CLEAR (self);
        return !!retval;
    }
    Py_CLEAR (self);
    return PGACheckStoppingConditions (ctx);
}

/*
 * Used only for user defined data type.
 * Copy the python object and update refcounts.
 * We explicitly do *NOT* do a deepcopy!
 * This is because the python implementation will in many cases allocate
 * a new object anyway, so copying the to-be-overwritten object once
 * doesn't make sense from a performance perspective.
 */
static void copystring (PGAContext *ctx, int p1, int pop1, int p2, int pop2)
{
    PGAIndividual *src = PGAGetIndividual (ctx, p1, pop1);
    PGAIndividual *dst = PGAGetIndividual (ctx, p2, pop2);
    ERR_CHECK_X_OCCURRED (ctx);
    ERR_CHECK_X (ctx, src->chrom);
    if (src == dst) {
        return;
    }
    if (dst->chrom != NULL) {
        Py_DECREF (dst->chrom);
        dst->chrom = NULL;
    }
    assert (src->chrom != NULL);
    dst->chrom = src->chrom;
    Py_INCREF (src->chrom);
errout:
    return;
}

/*
 * Used if the calling object has an initstring method.
 * Otherwise use built-in default for the datatype.
 * Note: In python this is also used when createstring requests
 * initialization of a string.
 */
static void initstring (PGAContext *ctx, int p, int pop)
{
    PyObject *self = NULL, *r = NULL;
    ERR_CHECK_X_OCCURRED (ctx);
    self = get_self (ctx);
    ERR_CHECK_X (ctx, self);
    r    = PyObject_CallMethod (self, "initstring", "ii", p, pop);
    ERR_CHECK_X (ctx, r);
errout:
    Py_CLEAR (r);
    Py_CLEAR (self);
    return;
}

/*
 * Used only for user defined data type.
 * If the initflag is zero, we do no initialization and assert that the
 * chromosome is a NULL pointer. Otherwise we call self.initstring.
 */
static void createstring (PGAContext *ctx, int p, int pop, int initflag)
{
    PGAIndividual *ind = PGAGetIndividual (ctx, p, pop);
    if (!initflag) {
        assert (ind->chrom == NULL);
        return;
    }
    initstring (ctx, p, pop);
}

/*
 * Used if the calling object has a crossover method.
 * Otherwise use built-in default for the datatype.
 * Perform crossover from p1 and p2 into c1 and c2
 */
static void crossover
    (PGAContext *ctx, int p1, int p2, int p_pop, int c1, int c2, int c_pop)
{
    PyObject *self = NULL, *r = NULL;
    ERR_CHECK_X_OCCURRED (ctx);
    self = get_self (ctx);
    ERR_CHECK_X (ctx, self);
    r    = PyObject_CallMethod
        (self, "crossover", "iiiiii", p1, p2, p_pop, c1, c2, c_pop);
    ERR_CHECK_X (ctx, r);
errout:
    Py_CLEAR (r);
    Py_CLEAR (self);
    return;
}

/*
 * Used if the calling object has a mutation method.
 * Otherwise use built-in default for the datatype.
 * Insert code to mutate Data.
 * Remember to count the number of mutations that happen, and return
 * that value!
 */
static int mutation (PGAContext *ctx, int p, int pop, double mr)
{
    PyObject *self = NULL, *r = NULL;
    int retval = 0, rr;
    ERR_CHECK_X_OCCURRED (ctx);
    self = get_self (ctx);
    ERR_CHECK_X (ctx, self);
    r    = PyObject_CallMethod (self, "mutation", "iid", p, pop, mr);
    ERR_CHECK_X (ctx, r);
    rr = PyArg_Parse (r, "i", &retval);
    ERR_CHECK_X (ctx, rr);
errout:
    Py_CLEAR (r);
    Py_CLEAR (self);
    return retval;
}

/*
 * Used if the calling object has a gene_distance method.
 * Otherwise use built-in default for the datatype.
 * Insert code to compute genetic difference of two individuals.
 */
static double gene_distance
    (PGAContext *ctx, int p1, int pop1, int p2, int pop2)
{
    PyObject *self = NULL, *r = NULL;
    int rr;
    double retval = 0.0;
    ERR_CHECK_X_OCCURRED (ctx);
    self = get_self (ctx);
    ERR_CHECK_X (ctx, self);
    r = PyObject_CallMethod
        (self, "gene_distance", "iiii", p1, pop1, p2, pop2);
    ERR_CHECK_X (ctx, r);
    rr = PyArg_Parse (r, "d", &retval);
    ERR_CHECK_X (ctx, rr);
errout:
    Py_CLEAR (r);
    Py_CLEAR (self);
    return retval;
}

/*
 * Used if the calling object has a pre_eval method.
 */
static void pre_eval (PGAContext *ctx, int pop)
{
    PyObject *self = NULL, *r = NULL;
    ERR_CHECK_X_OCCURRED (ctx);
    self = get_self (ctx);
    ERR_CHECK_X (ctx, self);
    r = PyObject_CallMethod (self, "pre_eval", "i", pop);
    ERR_CHECK_X (ctx, r);
errout:
    Py_CLEAR (r);
    Py_CLEAR (self);
    return;
}

static PyObject *get_file_from_fp (PGAContext *ctx, PyObject *self, FILE *fp)
{
    int retval = -1;
    PyObject *file = NULL;
    PyObject *pyfp = NULL;
    int fd = fileno (fp);
    
    /* Should never happen unles fp is not a valid stream */
    ERR_CHECK_ERRNO (ctx, fd >= 0, NULL);
#if IS_PY3
    file = PyFile_FromFd (fd, "", "w", -1, "utf-8", NULL, NULL, 0);
#else
    {
        FILE *fp2 = fdopen (fd, "w");
        ERR_CHECK_ERRNO (ctx, fp2 != NULL, NULL);
        file = PyFile_FromFile (fp2, "<PGA_file>", "w", NULL);
    }
#endif
    ERR_CHECK (ctx, file, NULL);
    pyfp = Py_BuildValue ("L", (long long) fp);
    ERR_CHECK (ctx, pyfp, NULL);
    retval = PyObject_SetAttrString (self, "_fp", pyfp);
    ERR_CHECK (ctx, retval >= 0, NULL);
    return file;
}

/*
 * Low-level gene print function
 * Note: We do a hack here and store both, the original FILE *fp and the
 * generated python file object into member variables _fp and _file,
 * respectively. These are reused in the print_gene function (and the
 * called print_string python method). This avoids garbage collection
 * and/or memory leak issues with files and FILE * variables. There is
 * no way to temporarily create a FILE * and avoiding the unterlying
 * file descriptor being closed in C. The previous version used dup2 to
 * duplicate the file descriptor but this is non-portable. The issue is
 * further complicated by the fact that Python3 (as opposed to Python2)
 * uses file descriptors not FILE* for generating a file object.
 */
static void print_gene (PGAContext *ctx, FILE *fp, int p, int pop)
{
    PyObject *self = NULL, *file = NULL, *r = NULL;
    int retval = -1;
    ERR_CHECK_X_OCCURRED (ctx);
    self    = get_self (ctx);
    ERR_CHECK_X (ctx, self);
    fflush (fp);

    file = get_file_from_fp (ctx, self, fp);
    ERR_CHECK_X (ctx, file);
    r = PyObject_CallMethod (self, "print_string", "Oii", file, p, pop);
    ERR_CHECK_X (ctx, r);
    Py_CLEAR (r);
    /* Flush file */
    r = PyObject_CallMethod (file, "flush", "");
    ERR_CHECK_X (ctx, r);
    Py_CLEAR (r);
    /* File is set to not close unterlying fd but close anyway */
    r = PyObject_CallMethod (file, "close", "");
    ERR_CHECK_X (ctx, r);
    Py_CLEAR (r);
    /* Now remove self._fp */
    retval = PyObject_DelAttrString (self, "_fp");
    ERR_CHECK_X (ctx, retval == 0);
errout:
    Py_CLEAR (file);
    Py_CLEAR (r);
    Py_CLEAR (self);
}

/******************
 * Serialization
 ******************/

static void *serialize_object  = NULL;
static void *serialize_inner   = NULL;
static PyObject *pickle_module = NULL;

/*
 * Used only for user defined data type.
 * This implementation relies on a serialization to be immediately used
 * and freed afterwards.
 */
static size_t serialize (PGAContext *ctx, int p, int pop, void **ser)
{
    PGAIndividual *ind = PGAGetIndividual (ctx, p, pop);
    Py_ssize_t serial_size = 0;
    ERR_CHECK_X_OCCURRED (ctx);
    ERR_CHECK_X (ctx, ind->chrom != NULL);
    ERR_CHECK_X (ctx, serialize_object == NULL);
    ERR_CHECK_X (ctx, serialize_inner  == NULL);
    if (pickle_module == NULL) {
        pickle_module = PyImport_ImportModule ("pickle");
    }
    ERR_CHECK_X (ctx, pickle_module);
    serialize_object = PyObject_CallMethod
        (pickle_module, "dumps", "O", ind->chrom);
    ERR_CHECK_X (ctx, serialize_object);
    serial_size = PyByteArray_Size (serialize_object);
    ERR_CHECK_X (ctx, serial_size >= 0);
    serialize_inner = PyBytes_AsString (serialize_object);
    ERR_CHECK_X (ctx, serialize_inner);
    goto out;
errout:
    Py_CLEAR (serialize_object);
    serialize_object = serialize_inner = NULL;
out:
    *ser = serialize_inner;
    return (size_t)serial_size;
}

/*
 * The pointer is immediately freed after being used, so we free the
 * corresponding python object here.
 */
static void serialize_free (void *p)
{
    PyObject *o = serialize_object;
    assert (o != NULL);
    assert (p == serialize_inner);
    Py_DECREF (o);
    serialize_object = NULL;
    serialize_inner  = NULL;
}

/*
 * Free the gene
 */
static void chrom_free (PGAIndividual *ind)
{
    if (ind->chrom != NULL) {
        Py_DECREF (ind->chrom);
        ind->chrom = NULL;
    }
}

static void deserialize
    (PGAContext *ctx, int p, int pop, const void *serial, size_t size)
{
    PGAIndividual *ind = PGAGetIndividual (ctx, p, pop);
    PyObject *serialized = NULL;
    PyObject *obj = NULL;
    ERR_CHECK_X_OCCURRED (ctx);
    serialized = PyBytes_FromStringAndSize (serial, size);
    ERR_CHECK_X (ctx, serialized);
    if (pickle_module == NULL) {
        pickle_module = PyImport_ImportModule ("pickle");
    }
    ERR_CHECK_X (ctx, pickle_module);
    obj = PyObject_CallMethod (pickle_module, "loads", "O", serialized);
    ERR_CHECK_X (ctx, obj);
    if (ind->chrom != NULL) {
        Py_DECREF (ind->chrom);
        ind->chrom = NULL;
    }
    /* No need to increment refcount */
    ind->chrom = obj;
errout:
    Py_CLEAR (serialized);
}

/******************
 * Helper functions
 ******************/

static int check_allele (PGAContext *ctx, int p, int pop, int i)
{
    ERR_CHECK_OCCURRED (ctx, 1);
    if (pop != PGA_OLDPOP && pop != PGA_NEWPOP) {
        char x [50];
        sprintf (x, "%d: invalid population", pop);
        PyErr_SetString (PyExc_ValueError, x);
        return 0;
    }
    if ((p < 0 || p >= ctx->ga.PopSize) && p != PGA_TEMP1 && p != PGA_TEMP2)
    {
        char x [50];
        sprintf (x, "%d: invalid population index", p);
        PyErr_SetString (PyExc_ValueError, x);
        return 0;
    }
    if (i < 0 || i >= ctx->ga.StringLen) {
        char x [50];
        sprintf (x, "%d: allele index out of range", i);
        PyErr_SetString (PyExc_ValueError, x);
        return 0;
    }
    return 1;
}

static int check_probability (double probability)
{
    CHECK_VALUE_EXCEPTION
        ( 0 <= probability && probability <= 1
        , "Probability must be 0 <= p <= 1"
        , PyExc_ValueError
        , 0
        );
    return 1;
}

#define check_interval(l,r) do { \
    if (l > r) {                                                           \
        PyErr_SetString                                                    \
            (PyExc_ValueError, "interval_left must be <= interval_right"); \
        return NULL;                                                       \
    }                                                                      \
} while (0)

/*
 * Helper method for calling a PGA C function from a Python sequence.
 * We get a callback method that calls a PGA function with a context as
 * the first parameter.
 * The sequence contains integers. We get a pointer into the constants
 * array and a length (number of valid constants for this PGA callback)
 * Constants are checked against the constants in the array. Optionally
 * an additional checkfun is called for each object.
 */
static int init_sequence
    ( PyObject   *sequence
    , PGAContext *ctx
    , void      (*fun)(PGAContext *, int)
    , char       *firstconst
    , int         constlen
    , char       *name
    , int       (*checkfun)(PGAContext *, int, char *)
    )
{
    constdef_t *constcheck, key;

    key.cd_name = firstconst;
    constcheck = bsearch
            ( &key
            , constdef
            , sizeof (constdef) / sizeof (constdef_t) - 1
            , sizeof (constdef_t)
            , compare_constdef
            );
    assert (constcheck);
    if (sequence) {
        Py_ssize_t i, len = PySequence_Length (sequence);
        if (len < 0) {
            return 0;
        }
        for (i = 0; i < len; i++) {
            PyObject *x = PySequence_GetItem (sequence, i);
            int val, ci;
            if (!x) {
                return 0;
            }
            if (!PyArg_Parse (x, "i", &val)) {
                return 0;
            }
            for (ci = 0; ci < constlen; ci++) {
                if (val == constcheck [ci].cd_value) {
                    break;
                }
            }
            CHECK_VALUE_EXCEPTION
                ( ci != constlen
                , name
                , PyExc_ValueError
                , 0
                );
            if (checkfun && !checkfun (ctx, val, name)) {
                return 0;
            }
            fun (ctx, val);
        }
    }
    return 1;
}

static int check_hamming (PGAContext *ctx, int val, char *name)
{
    if  (  val == PGA_REPORT_HAMMING
        && PGAGetDataType (ctx) != PGA_DATATYPE_BINARY
        )
    {
        char errmsg [1024];
        sprintf
            (errmsg, "%s: no hamming reporting for non-binary datatype", name);
        PyErr_SetString
            (PyExc_ValueError, errmsg);
        return 0;
    }
    return 1;
}

Py_ssize_t parse_points (int dim, PyObject *points, void **result)
{
    int i, j;
    Py_ssize_t length;
    PyObject *point = NULL;
    PyObject *res = NULL, *res2 = NULL;
    double *refpoints = NULL;

    CHECK_VALUE_EXCEPTION
        ( PySequence_Check (points)
        , "Expected sequence for points parameter"
        , PyExc_ValueError
        , 0
        );
    length = PySequence_Length (points);
    CHECK_VALUE_EXCEPTION
        ( length >= 1
        , "Must at least specify one point"
        , PyExc_ValueError
        , 0
        );
    refpoints = malloc (sizeof (double) * length * dim);
    if (refpoints == NULL) {
        PyErr_NoMemory ();
        goto errout;
    }
    for (i=0; i<length; i++) {
        Py_ssize_t l;

        point = PySequence_GetItem (points, i);
        ERR_CHECK_RET (point);
        l = PySequence_Length (point);
        if (l != dim) {
            char x [100];
            sprintf (x, "Points must have dimension %d", dim);
            PyErr_SetString (PyExc_ValueError, x);
            goto errout;
        }
        for (j=0; j<dim; j++) {
            res = PySequence_GetItem (point, j);
            ERR_CHECK_RET (res);
            res2 = PyNumber_Float (res);
            ERR_CHECK_RET (res2);
            Py_CLEAR (res);
            refpoints [dim * i + j] = PyFloat_AsDouble (res2);
            if (PyErr_Occurred ()) {
                goto errout;
            }
            Py_CLEAR (res2);
        }
        Py_CLEAR (point);
    }
    *result = refpoints;
    return length;
errout:
    Py_CLEAR (point);
    Py_CLEAR (res);
    Py_CLEAR (res2);
    if (refpoints != NULL) {
        free (refpoints);
    }
    return 0;
}

/* Parse an array of strings into a C argv vector, return argc in argcp */
static char **parse_argv (PyObject *argv, int *argcp)
{
    int i;
    char **c_argv;
    Py_ssize_t argc_py = PySequence_Length (argv);
    /* Avoid warnings on size of int vs Py_ssize_t */
    if (argc_py < 0) {
        return NULL;
    }
    assert (argc_py < INT_MAX);
    *argcp = (int)argc_py;
    c_argv = malloc ((*argcp + 1) * sizeof (char *));
    c_argv [*argcp] = NULL;
    for (i = 0; i < *argcp; i++) {
        Py_ssize_t len = 0;
        PyObject *b = NULL;
        PyObject *s = PySequence_GetItem (argv, i);
        if (!s) {
            return NULL;
        }
        b = PyUnicode_AsEncodedString (s, "utf-8", "strict");
        Py_DECREF (s);
        if (!b) {
            return NULL;
        }
        len = PyBytes_Size (b);
        c_argv [i] = malloc (len + 1);
        strcpy (c_argv [i], PyBytes_AsString (b));
        Py_DECREF (b);
    }
    return c_argv;
}

/*************
 * Constructor
 *************/

#define INIT_FAIL -1

/* Registered as atexit function and also called by destructor */
static void exitfunc (void)
{
    int mpi_finalized = 0;
    MPI_Finalized (&mpi_finalized);
    if (!mpi_finalized) {
        MPI_Finalize ();
    }
}

static int PGA_init (PyObject *self, PyObject *args, PyObject *kw)
{
    int argc = 0, max = 0, length = 0, pop_size = 0, pga_type = 0;
    int num_eval = 1;
    int random_seed = 0, max_GA_iter = 0, max_no_change = 0;
    int num_replace = -1, pop_replace_type = -1;
    int max_similarity = 0, crossover_type = -1, select_type = -1;
    int print_frequency = 0;
    int restart_frequency = 0;
    int mutation_type = 0;
    int fitness_type = 0;
    int fitness_min_type = 0;
    double tournament_size = 0;
    int rtr_window_size = 0;
    int tournament_with_replacement = -1;
    int randomize_select = -1;
    int num_constraint = -1;
    int sum_constraints = -1;
    double truncation_proportion = 0.0;
    double mutation_value = 0.0;
    double mutation_poly_eta = -1;
    double mutation_poly_value = -1;
    int DE_variant = 0;
    int DE_num_diffs = 0;
    int DE_crossover_type = 0;
    double DE_scale_factor = -1;
    double DE_aux_factor = -1;
    double DE_crossover_prob = -1;
    double DE_jitter = -1;
    double DE_dither = -1;
    double DE_probability_EO = -1;
    PyObject *DE_dither_per_individual = NULL;
    double mutation_prob = -1;
    double crossover_prob = -1;
    PyObject *crossover_bounded = NULL;
    PyObject *crossover_bounce_back = NULL;
    PyObject *crossover_SBX_once_per_string = NULL;
    double uniform_crossover_prob = -1.0;
    double crossover_SBX_eta = -1.0;
    double p_tournament_prob = -1.0;
    double max_fitness_rank = -1.0;
    double fitness_cmax = -1.0;
    PyObject *PGA_ctx = NULL;
    PyObject *type = NULL, *maximize = NULL, *init = NULL;
    PyObject *init_percent = NULL, *stopping_rule_types = NULL;
    PyObject *integer_init_permute = NULL;
    PyObject *print_options = NULL;
    PyObject *no_duplicates = NULL;
    PyObject *restart = NULL;
    PyObject *mutation_bounded = NULL;
    PyObject *mutation_bounce_back = NULL;
    PyObject *mutation_and_crossover = NULL;
    PyObject *mutation_or_crossover = NULL;
    PyObject *mutation_only = NULL;
    int mixing_type = 0;
    PyObject *refpoints = NULL;
    PyObject *refdirs = NULL;
    int refdir_partitions = 0;
    double refdir_scale = 0.05;
    PyObject *argv = NULL;
    char **c_argv = NULL;
    PGAContext *ctx;
    PyObject *fixed_edges = NULL;
    int fixed_edges_symmetric = PGA_TRUE;
    int epsilon_generation = 0;
    double epsilon_exponent = -1;
    int epsilon_theta = -1;
    int multi_obj_precision = -1;
    int retval = -1;
    int mpi_initialized = 0;
    int char_init_type = -1;
    PyObject *Py_MPI_Initialized = NULL;
    PyObject *output_file = NULL;
    static char *kwlist[] =
        { "type"
        , "length"
        , "maximize"
        , "pop_size"
        , "num_eval"
        , "init"
        , "init_percent"
        , "integer_init_permute"
        , "random_seed"
        , "max_GA_iter"
        , "max_no_change"
        , "max_similarity"
        , "mutation_prob"
        , "stopping_rule_types"
        , "num_replace"
        , "pop_replace_type"
        , "print_options"
        , "no_duplicates"
        , "DE_dither_per_individual"
        , "crossover_type"
        , "select_type"
        , "crossover_prob"
        , "crossover_bounded"
        , "crossover_bounce_back"
        , "crossover_SBX_once_per_string"
        , "crossover_SBX_eta"
        , "uniform_crossover_prob"
        , "print_frequency"
        , "restart"
        , "restart_frequency"
        , "mutation_bounded"
        , "mutation_bounce_back"
        , "mutation_type"
        , "mutation_value"
        , "mutation_poly_eta"
        , "mutation_poly_value"
        , "DE_variant"
        , "DE_num_diffs"
        , "DE_crossover_type"
        , "DE_scale_factor"
        , "DE_aux_factor"
        , "DE_crossover_prob"
        , "DE_jitter"
        , "DE_dither"
        , "DE_probability_EO"
        , "mutation_and_crossover"
        , "mutation_or_crossover"
        , "mutation_only"
        , "mixing_type"
        , "p_tournament_prob"
        , "fitness_type"
        , "max_fitness_rank"
        , "fitness_cmax"
        , "fitness_min_type"
        , "tournament_size"
        , "rtr_window_size"
        , "tournament_with_replacement"
        , "truncation_proportion"
        , "randomize_select"
        , "num_constraint"
        , "sum_constraints"
        , "reference_points"
        , "reference_directions"
        , "refdir_partitions"
        , "refdir_scale"
        , "argv"
        , "fixed_edges"
        , "fixed_edges_symmetric"
        , "epsilon_generation"
        , "epsilon_exponent"
        , "epsilon_theta"
        , "multi_obj_precision"
        , "output_file"
        , "char_init_type"
        , NULL
        };

    if  (!PyArg_ParseTupleAndKeywords
            ( args
            , kw
            , "Oi|OiiOOOiiiidOiiOOOiidOOOddiOiOOiddd"
              "iiiddddddOOOididdidiidiiiOOidOOiidiiOi"
            , kwlist
            , &type
            , &length
            , &maximize
            , &pop_size
            , &num_eval
            , &init
            , &init_percent
            , &integer_init_permute
            , &random_seed
            , &max_GA_iter
            , &max_no_change
            , &max_similarity
            , &mutation_prob
            , &stopping_rule_types
            , &num_replace
            , &pop_replace_type
            , &print_options
            , &no_duplicates
            , &DE_dither_per_individual
            , &crossover_type
            , &select_type
            , &crossover_prob
            , &crossover_bounded
            , &crossover_bounce_back
            , &crossover_SBX_once_per_string
            , &crossover_SBX_eta
            , &uniform_crossover_prob
            , &print_frequency
            , &restart
            , &restart_frequency
            , &mutation_bounded
            , &mutation_bounce_back
            , &mutation_type
            , &mutation_value
            , &mutation_poly_eta
            , &mutation_poly_value
            , &DE_variant
            , &DE_num_diffs
            , &DE_crossover_type
            , &DE_scale_factor
            , &DE_aux_factor
            , &DE_crossover_prob
            , &DE_jitter
            , &DE_dither
            , &DE_probability_EO
            , &mutation_and_crossover
            , &mutation_or_crossover
            , &mutation_only
            , &mixing_type
            , &p_tournament_prob
            , &fitness_type
            , &max_fitness_rank
            , &fitness_cmax
            , &fitness_min_type
            , &tournament_size
            , &rtr_window_size
            , &tournament_with_replacement
            , &truncation_proportion
            , &randomize_select
            , &num_constraint
            , &sum_constraints
            , &refpoints
            , &refdirs
            , &refdir_partitions
            , &refdir_scale
            , &argv
            , &fixed_edges
            , &fixed_edges_symmetric
            , &epsilon_generation
            , &epsilon_exponent
            , &epsilon_theta
            , &multi_obj_precision
            , &output_file
            , &char_init_type
            )
        )
    {
        return INIT_FAIL;
    }

    if (PyObject_IsSubclass (type, (PyObject *)&PyBool_Type)) {
        pga_type = PGA_DATATYPE_BINARY;
    }
    else if (  PyObject_IsSubclass (type, (PyObject *)&PyInt_Type_Compat)
            || PyObject_IsSubclass (type, (PyObject *)&PyLong_Type)
            )
    {
        pga_type = PGA_DATATYPE_INTEGER;
    } else if (PyObject_IsSubclass (type, (PyObject *)&PyFloat_Type)) {
        pga_type = PGA_DATATYPE_REAL;
    } else if (PyObject_IsSubclass (type, (PyObject *)&PyBytes_Type)) {
        pga_type = PGA_DATATYPE_CHARACTER;
    } else {
        pga_type = PGA_DATATYPE_USER;
    }

    if (maximize) {
        max = PyObject_IsTrue (maximize);
    }
    CHECK_VALUE (length > 1, "Gene length must be at least 2");
    /* If user didn't specify argv we get sys.argv */
    if (!argv) {
        PyObject *sys = PyImport_ImportModule ("sys");
        if (!sys) {
            return INIT_FAIL;
        }
        argv = PyObject_GetAttrString (sys, "argv");
        Py_DECREF (sys);
        if (!argv) {
            return INIT_FAIL;
        }
    }
    c_argv = parse_argv (argv, &argc);
    if (c_argv == NULL) {
        return INIT_FAIL;
    }

    /* Only initialize MPI if not already done. Install atexit handler
     * only if MPI is not yet initialized.
     */
    MPI_Initialized (&mpi_initialized);
    if (!mpi_initialized) {
        MPI_Init (&argc, &c_argv);
        retval = atexit (exitfunc);
        CHECK_VALUE_EXCEPTION
            ( retval == 0
            , "Cannot register exit function"
            , PyExc_RuntimeError
            , INIT_FAIL
            );
    }
    Py_MPI_Initialized = Py_BuildValue ("i", mpi_initialized);
    if (PyObject_SetAttrString
           (self, "mpi_initialized", Py_MPI_Initialized) < 0
       )
    {
        return INIT_FAIL;
    }

    ctx = PGACreate
        ( &argc
        , c_argv
        , pga_type
        , length
        , max ? PGA_MAXIMIZE : PGA_MINIMIZE
        );
    /* handle context pointer */
    PGA_ctx = Py_BuildValue ("L", (long long) ctx);
    if (PGA_ctx == NULL) {
        return INIT_FAIL;
    }
    assert (contexts);
    assert (self);
    PyObject_SetItem (contexts, PGA_ctx, self);

    ERR_DECREF_RET
        ( PyObject_SetAttrString (self, "context", PGA_ctx) >= 0
        , PGA_ctx, INIT_FAIL
        );
    Py_CLEAR (PGA_ctx);
    /*
     * Allocate data structure for error indicator, for now this is just
     * an integer flag. This is used for terminating the search and
     * raising an error outside
     */
    assert (ctx->ga.CustomData == NULL);
    ctx->ga.CustomData = malloc (sizeof (int));
    if (ctx->ga.CustomData == NULL) {
        PyErr_NoMemory ();
        return INIT_FAIL;
    }
    *((int *)ctx->ga.CustomData) = 0;

    /* If using userdefined datatypes we also set the user functions
     * because PGAPack requires these and for many use-cases they are
     * not used. The effect when such a function is needed (and called)
     * is that a traceback is raised that the method is non-existing.
     * This is considered better than PGAPack terminating with an error.
     * Some of the functions are *only* defined internally and do not
     * call into python methods.
     */
    if (  PyObject_HasAttrString (self, "check_duplicate")
       || ctx->ga.datatype == PGA_DATATYPE_USER
       )
    {
        PGASetUserFunction
            (ctx, PGA_USERFUNCTION_DUPLICATE, (void *)check_duplicate);
    }
    if (  PyObject_HasAttrString (self, "hash")
       || ctx->ga.datatype == PGA_DATATYPE_USER
       )
    {
        PGASetUserFunction (ctx, PGA_USERFUNCTION_HASH, (void *)build_hash);
    }
    PGASetUserFunction (ctx, PGA_USERFUNCTION_STOPCOND, (void *)check_stop);
    if (ctx->ga.datatype == PGA_DATATYPE_USER) {
        PGASetUserFunction
            (ctx, PGA_USERFUNCTION_COPYSTRING, (void *)copystring);
    }
    if (ctx->ga.datatype == PGA_DATATYPE_USER) {
        PGASetUserFunction
            (ctx, PGA_USERFUNCTION_CREATESTRING, (void *)createstring);
    }
    if (  PyObject_HasAttrString (self, "crossover")
       || ctx->ga.datatype == PGA_DATATYPE_USER
       )
    {
        PGASetUserFunction (ctx, PGA_USERFUNCTION_CROSSOVER, (void *)crossover);
    }
    if (ctx->ga.datatype == PGA_DATATYPE_USER) {
        PGASetUserFunction
            (ctx, PGA_USERFUNCTION_DESERIALIZE, (void *)deserialize);
    }
    if (PyObject_HasAttrString (self, "endofgen")) {
        PGASetUserFunction (ctx, PGA_USERFUNCTION_ENDOFGEN, (void *)endofgen);
    }
    if (  PyObject_HasAttrString (self, "gene_distance")
       || ctx->ga.datatype == PGA_DATATYPE_USER
       )
    {
        PGASetUserFunction
            (ctx, PGA_USERFUNCTION_GEN_DISTANCE, (void *)gene_distance);
    }
    if (  PyObject_HasAttrString (self, "initstring")
       || ctx->ga.datatype == PGA_DATATYPE_USER
       )
    {
        PGASetUserFunction
            (ctx, PGA_USERFUNCTION_INITSTRING, (void *)initstring);
    }
    if (  PyObject_HasAttrString (self, "mutation")
       || ctx->ga.datatype == PGA_DATATYPE_USER
       )
    {
        PGASetUserFunction (ctx, PGA_USERFUNCTION_MUTATION, (void *)mutation);
    }
    if (PyObject_HasAttrString (self, "pre_eval")) {
        PGASetUserFunction
            (ctx, PGA_USERFUNCTION_PRE_EVAL, (void *)pre_eval);
    }
    PGASetUserFunction (ctx, PGA_USERFUNCTION_PRINTSTRING, (void *)print_gene);
    if (ctx->ga.datatype == PGA_DATATYPE_USER) {
        PGASetUserFunction
            (ctx, PGA_USERFUNCTION_SERIALIZE, (void *)serialize);
    }
    if (ctx->ga.datatype == PGA_DATATYPE_USER) {
        PGASetUserFunction
            (ctx, PGA_USERFUNCTION_SERIALIZE_FREE, (void *)serialize_free);
    }
    if (ctx->ga.datatype == PGA_DATATYPE_USER) {
        PGASetUserFunction
            (ctx, PGA_USERFUNCTION_CHROM_FREE, (void *)chrom_free);
    }

    if (crossover_prob >= 0) {
        PGASetCrossoverProb (ctx, crossover_prob);
    }
    if (crossover_bounded && PyObject_IsTrue (crossover_bounded)) {
        PGASetCrossoverBoundedFlag (ctx, PGA_TRUE);
    }
    if (crossover_bounce_back && PyObject_IsTrue (crossover_bounce_back)) {
        PGASetCrossoverBounceBackFlag (ctx, PGA_TRUE);
    }
    if (  crossover_SBX_once_per_string
       && PyObject_IsTrue (crossover_SBX_once_per_string)
       )
    {
        PGASetCrossoverSBXOncePerString (ctx, PGA_TRUE);
    }
    if (crossover_SBX_eta >= 0) {
        PGASetCrossoverSBXEta (ctx, crossover_SBX_eta);
    }
    if (crossover_type >= 0) {
        CHECK_VALUE
            ( (  crossover_type == PGA_CROSSOVER_ONEPT
              || crossover_type == PGA_CROSSOVER_TWOPT
              || crossover_type == PGA_CROSSOVER_UNIFORM
              || crossover_type == PGA_CROSSOVER_EDGE
              )
            , "invalid crossover_type"
            );
        PGASetCrossoverType (ctx, crossover_type);
    }
    if (max_GA_iter) {
        CHECK_VALUE (max_GA_iter >= 2, "Iteration count must be at least 2");
        PGASetMaxGAIterValue (ctx, max_GA_iter);
    }
    if (max_no_change) {
        PGASetMaxNoChangeValue (ctx, max_no_change);
    }
    if (max_similarity) {
        PGASetMaxSimilarityValue (ctx, max_similarity);
    }
    if (mutation_prob >= 0) {
        PGASetMutationProb (ctx, mutation_prob);
    }
    if (DE_num_diffs > 0) {
        CHECK_VALUE (DE_num_diffs <= 2, "invalid DE_num_diffs");
        PGASetDENumDiffs (ctx, DE_num_diffs);
    }
    if  (  DE_dither_per_individual
        && PyObject_IsTrue (DE_dither_per_individual)
        )
    {
        PGASetDEDitherPerIndividual (ctx, PGA_TRUE);
    }
    if (DE_crossover_type > 0) {
        CHECK_VALUE
            ( ( DE_crossover_type == PGA_DE_CROSSOVER_BIN
              || DE_crossover_type == PGA_DE_CROSSOVER_EXP
              )
            , "invalid DE crossover_type"
            );
        PGASetDECrossoverType (ctx, DE_crossover_type);
    }
    if (DE_scale_factor > 0) {
        CHECK_VALUE (DE_scale_factor <= 2, "invalid DE_scale_factor");
        PGASetDEScaleFactor (ctx, DE_scale_factor);
    }
    if (DE_aux_factor > 0) {
        CHECK_VALUE (DE_aux_factor <= 2, "invalid DE_aux_factor");
        PGASetDEAuxFactor (ctx, DE_aux_factor);
    }
    if (DE_crossover_prob > 0) {
        CHECK_VALUE (DE_crossover_prob <= 1, "invalid DE_crossover_prob");
        PGASetDECrossoverProb (ctx, DE_crossover_prob);
    }
    if (DE_jitter > 0) {
        CHECK_VALUE (DE_jitter <= 1, "invalid DE_jitter");
        PGASetDEJitter (ctx, DE_jitter);
    }
    if (DE_dither > 0) {
        CHECK_VALUE (DE_dither <= 1, "invalid DE_dither");
        PGASetDEDither (ctx, DE_dither);
    }
    if (DE_probability_EO > 0) {
        CHECK_VALUE (DE_probability_EO <= 1, "invalid DE_probability_EO");
        PGASetDEProbabilityEO (ctx, DE_probability_EO);
    }
    if (no_duplicates && PyObject_IsTrue (no_duplicates)) {
        PGASetNoDuplicatesFlag (ctx, PGA_TRUE);
    }
    if (num_replace >= 0) {
        PGASetNumReplaceValue (ctx, num_replace);
    }
    if (pop_replace_type >= 0) {
        CHECK_VALUE
            ( (  pop_replace_type == PGA_POPREPL_BEST
              || pop_replace_type == PGA_POPREPL_RANDOM_NOREP
              || pop_replace_type == PGA_POPREPL_RANDOM_REP
              || pop_replace_type == PGA_POPREPL_RTR
              || pop_replace_type == PGA_POPREPL_PAIRWISE_BEST
              || pop_replace_type == PGA_POPREPL_NSGA_II
              || pop_replace_type == PGA_POPREPL_NSGA_III
              )
            , "invalid pop_replace_type"
            );
        PGASetPopReplaceType (ctx, pop_replace_type);
    }
    if (pop_size) {
        CHECK_VALUE (pop_size > 1, "Population size must be at least 2");
        CHECK_VALUE (!(pop_size & 1), "Population size must be an even number");
        PGASetPopSize (ctx, pop_size);
    }
    if (num_eval != 1) {
        CHECK_VALUE (num_eval >= 1, "Number of evaluations must be at least 1");
        PGASetNumAuxEval (ctx, num_eval - 1);
    }
    if (print_frequency) {
        PGASetPrintFrequencyValue (ctx, print_frequency);
    }
    if (random_seed) {
        PGASetRandomSeed   (ctx, random_seed);
    }
    if (select_type >= 0) {
        CHECK_VALUE
            ( (  select_type == PGA_SELECT_PROPORTIONAL
              || select_type == PGA_SELECT_SUS
              || select_type == PGA_SELECT_TOURNAMENT
              || select_type == PGA_SELECT_PTOURNAMENT
              || select_type == PGA_SELECT_TRUNCATION
              || select_type == PGA_SELECT_LINEAR
              )
            , "invalid select_type"
            );
        PGASetSelectType (ctx, select_type);
    }
    if (uniform_crossover_prob >= 0) {
        PGASetUniformCrossoverProb (ctx, uniform_crossover_prob);
    }
    if (p_tournament_prob >= 0) {
        PGASetPTournamentProb (ctx, p_tournament_prob);
    }
    if (max_fitness_rank >= 0) {
        PGASetMaxFitnessRank (ctx, max_fitness_rank);
    }
    if (fitness_cmax >= 0) {
        PGASetFitnessCmaxValue (ctx, fitness_cmax);
    }

    if  (!init_sequence
            ( stopping_rule_types
            , ctx
            , PGASetStoppingRuleType
            , "PGA_STOP_MAXITER"
            , 3
            , "stopping_rule_types"
            , NULL
            )
        )
    {
        return INIT_FAIL;
    }
    if  (!init_sequence
            ( print_options
            , ctx
            , PGASetPrintOptions
            , "PGA_REPORT_AVERAGE"
            , 6
            , "print_options"
            , check_hamming
            )
        )
    {
        return INIT_FAIL;
    }
    if (integer_init_permute) {
        int i;
        Py_ssize_t len = PySequence_Length (integer_init_permute);
        int permute_lh [2];
        if (len < 0) {
            return INIT_FAIL;
        }
        CHECK_VALUE (len == 2, "Need lower, upper for integer_init_permute");
        for (i = 0; i < 2; i++) {
            PyObject *x = PySequence_GetItem (integer_init_permute, i);
            PyObject *l = NULL;
            if (!x) {
                return INIT_FAIL;
            }
            l = PyNumber_Long (x);
            Py_CLEAR (x);
            if (!l) {
                return INIT_FAIL;
            }
            ERR_DECREF_RET (PyArg_Parse (l, "i", permute_lh + i), l, INIT_FAIL);
            Py_CLEAR (l);
        }
        CHECK_VALUE
            ( permute_lh [1] - permute_lh [0] == length - 1
            , "integer_init_permute: upper - lower must be gene length"
            );
        PGASetIntegerInitPermute (ctx, permute_lh [0], permute_lh [1]);
    }
    if (init || init_percent) {
        int datatype = PGAGetDataType (ctx);
        int is_real  = (datatype == PGA_DATATYPE_REAL);
        int i;
        Py_ssize_t len;
        void *i_low, *i_high;
        PyObject *initvals = (init ? init : init_percent);

        CHECK_VALUE
            ( (  datatype == PGA_DATATYPE_INTEGER
              || datatype == PGA_DATATYPE_REAL
              )
            , "Init only for int/real"
            );
        CHECK_VALUE (!(init && init_percent), "Only one of init/init_percent");
        CHECK_VALUE
            (!(init_percent && !is_real), "init_percent only for float");
        len = PySequence_Length (initvals);
        if (len < 0) {
            return INIT_FAIL;
        }
        CHECK_VALUE
            (len == PGAGetStringLength (ctx), "Init length != string length");
        i_low  = malloc (len * (is_real ? sizeof (double) : sizeof (int)));
        i_high = malloc (len * (is_real ? sizeof (double) : sizeof (int)));
        if (i_low == NULL || i_high == NULL) {
            PyErr_NoMemory ();
            return INIT_FAIL;
        }
        for (i = 0; i < len; i++) {
            PyObject *x = PySequence_GetItem (initvals, i);
            PyObject *low = NULL, *high = NULL;
            if (!x) {
                return INIT_FAIL;
            }
            low  = PySequence_GetItem (x, 0);
            ERR_DECREF_RET (low, x, INIT_FAIL);
            high = PySequence_GetItem (x, 1);
            ERR_DECREF_2_RET (high, x, low, INIT_FAIL);
            Py_CLEAR (x);
            if (is_real) {
                PyObject *l = NULL, *h = NULL;
                double hi;
                l = PyNumber_Float (low);
                ERR_DECREF_2_RET (l, low, high, INIT_FAIL);
                h = PyNumber_Float (high);
                ERR_DECREF_3_RET (h, low, high, l, INIT_FAIL);
                Py_CLEAR (low);
                Py_CLEAR (high);
                ERR_DECREF_2_RET
                    ( (  PyArg_Parse (h, "d", ((double *)i_high) + i)
                      && PyArg_Parse (l, "d", ((double *)i_low)  + i)
                      )
                    , l, h, INIT_FAIL
                    );
                hi = ((double *)i_high) [i];
                if (init_percent) {
                    CHECK_VALUE
                        (0 <= hi && hi <= 1, "Percentage must be 0 < p <= 1");
                }
            } else {
                PyObject *l = NULL, *h = NULL;
                l = PyNumber_Long (low);
                ERR_DECREF_2_RET (l, low, high, INIT_FAIL);
                h = PyNumber_Long (high);
                ERR_DECREF_3_RET (h, low, high, l, INIT_FAIL);
                Py_CLEAR (low);
                Py_CLEAR (high);
                ERR_DECREF_2_RET
                    ( (  PyArg_Parse (h, "i", ((int *)i_high) + i)
                      && PyArg_Parse (l, "i", ((int *)i_low)  + i)
                      )
                    , l, h, INIT_FAIL
                    );
            }
        }
        if (is_real) {
            if (init) {
                PGASetRealInitRange    (ctx, i_low, i_high);
            } else {
                PGASetRealInitPercent  (ctx, i_low, i_high);
            }
        } else {
            PGASetIntegerInitRange (ctx, i_low, i_high);
        }
        free (i_low);
        free (i_high);
    }
    if (restart && PyObject_IsTrue (restart)) {
        PGASetRestartFlag (ctx, PGA_TRUE);
    }
    if (restart_frequency) {
        PGASetRestartFrequencyValue (ctx, restart_frequency);
    }
    if (mutation_bounded && PyObject_IsTrue (mutation_bounded)) {
        PGASetMutationBoundedFlag (ctx, PGA_TRUE);
    }
    if (mutation_bounce_back && PyObject_IsTrue (mutation_bounce_back)) {
        PGASetMutationBounceBackFlag (ctx, PGA_TRUE);
    }
    if (mutation_and_crossover && PyObject_IsTrue (mutation_and_crossover)) {
        PGASetMutationAndCrossoverFlag (ctx, PGA_TRUE);
    }
    if (mutation_or_crossover && PyObject_IsTrue (mutation_or_crossover)) {
        PGASetMutationOrCrossoverFlag (ctx, PGA_TRUE);
    }
    if (mutation_only && PyObject_IsTrue (mutation_only)) {
        PGASetMutationOnlyFlag (ctx, PGA_TRUE);
    }
    if (mixing_type != 0) {
        CHECK_VALUE
            ( (  mixing_type == PGA_MIX_TRADITIONAL
              || mixing_type == PGA_MIX_MUTATE_OR_CROSS
              || mixing_type == PGA_MIX_MUTATE_AND_CROSS
              || mixing_type == PGA_MIX_MUTATE_ONLY
              )
            , "invalid mixing_type"
            );
        PGASetMixingType (ctx, mixing_type);
    }
    if (mutation_type) {
        int data_type = PGAGetDataType (ctx);
        CHECK_VALUE
            ( (  data_type == PGA_DATATYPE_REAL
              || data_type == PGA_DATATYPE_INTEGER
              )
            , "mutation_type only for int and float"
            );
        if (PGAGetDataType (ctx) == PGA_DATATYPE_REAL) {
            CHECK_VALUE
                ( (  mutation_type == PGA_MUTATION_CONSTANT
                  || mutation_type == PGA_MUTATION_GAUSSIAN
                  || mutation_type == PGA_MUTATION_RANGE
                  || mutation_type == PGA_MUTATION_UNIFORM
                  || mutation_type == PGA_MUTATION_DE
                  )
                , "invalid mutation_type for float"
                );
        }
        if (PGAGetDataType (ctx) == PGA_DATATYPE_INTEGER) {
            CHECK_VALUE
                ( (  mutation_type == PGA_MUTATION_CONSTANT
                  || mutation_type == PGA_MUTATION_PERMUTE
                  || mutation_type == PGA_MUTATION_RANGE
                  )
                , "invalid mutation_type for int"
                );
        }
        PGASetMutationType (ctx, mutation_type);
    }
    if (mutation_value) {
        if (PGAGetDataType (ctx) == PGA_DATATYPE_REAL) {
            PGASetMutationRealValue (ctx, mutation_value);
        } else {
            PGASetMutationIntegerValue (ctx, (int)mutation_value);
        }
    }
    if (mutation_poly_eta >= 0) {
        PGASetMutationPolyEta (ctx, mutation_poly_eta);
    }
    if (mutation_poly_value >= 0) {
        PGASetMutationPolyValue (ctx, mutation_poly_value);
    }
    if (DE_variant) {
        CHECK_VALUE
            ( (  DE_variant == PGA_DE_VARIANT_RAND
              || DE_variant == PGA_DE_VARIANT_BEST
              || DE_variant == PGA_DE_VARIANT_EITHER_OR
              )
            , "invalid DE_variant"
            );
        PGASetDEVariant (ctx, DE_variant);
    }
    if (fitness_type) {
        CHECK_VALUE
            ( (  fitness_type == PGA_FITNESS_RANKING
              || fitness_type == PGA_FITNESS_RAW
              || fitness_type == PGA_FITNESS_NORMAL
              )
            , "invalid fitness_type"
            );
        PGASetFitnessType (ctx, fitness_type);
    }
    if (fitness_min_type) {
        CHECK_VALUE
            ( (  fitness_min_type == PGA_FITNESSMIN_CMAX
              || fitness_min_type == PGA_FITNESSMIN_RECIPROCAL
              )
            , "invalid fitness_min_type"
            );
        PGASetFitnessMinType (ctx, fitness_min_type);
    }
    if (tournament_size) {
        PGASetTournamentSize (ctx, tournament_size);
    }
    if (rtr_window_size) {
        PGASetRTRWindowSize (ctx, rtr_window_size);
    }
    if (tournament_with_replacement >= 0) {
        int v = tournament_with_replacement ? PGA_TRUE : PGA_FALSE;
        PGASetTournamentWithReplacement (ctx, v);
    }
    if (truncation_proportion) {
        PGASetTruncationProportion (ctx, truncation_proportion);
    }
    if (randomize_select >= 0) {
        int v = randomize_select ? PGA_TRUE : PGA_FALSE;
        PGASetRandomizeSelect (ctx, v);
    }
    if (num_constraint >= 0) {
        PGASetNumConstraint (ctx, num_constraint);
    }
    if (sum_constraints >= 0) {
        PGASetSumConstraintsFlag (ctx, sum_constraints);
    }
    if (refpoints != NULL) {
        Py_ssize_t npoints = 0;
        void *ref = NULL;
        npoints = parse_points (num_eval - num_constraint, refpoints, &ref);
        if (npoints == 0) {
            return INIT_FAIL;
        }
        PGASetReferencePoints (ctx, npoints, ref);
    }
    if (refdirs != NULL) {
        Py_ssize_t ndirs = 0;
        void *ref = NULL;
        CHECK_VALUE
            ( 0 < refdir_scale && refdir_scale <= 1
            , "Need 0 < refdir_scale <= 1"
            );
        CHECK_VALUE (refdir_partitions > 0, "Need refdir_partitions > 0");
        ndirs = parse_points (num_eval - num_constraint, refdirs, &ref);
        if (ndirs == 0) {
            return INIT_FAIL;
        }
        PGASetReferenceDirections
            (ctx, ndirs, ref, refdir_partitions, refdir_scale);
    }
    if (fixed_edges) {
        Py_ssize_t i, j, len = PySequence_Length (fixed_edges);
        PGAInteger (*fix)[2];
        CHECK_VALUE
            ( len > 0
            , "Fixed edges must be a sequence of sequences of length 2"
            );
        fix = malloc (len * 2 * sizeof (PGAInteger));
        if (fix == NULL) {
            PyErr_NoMemory ();
            return INIT_FAIL;
        }
        for (i=0; i<len; i++) {
            PyObject *s = PySequence_GetItem (fixed_edges, i);
            CHECK_VALUE_AND_FREE
                ( PySequence_Length (s) == 2
                , "Fixed edges must be a sequence of sequences of length 2"
                , PyExc_ValueError
                , fix
                );
            for (j=0; j<2; j++) {
                PyObject *v = PySequence_GetItem (s, j);
                int ret = PyArg_Parse (v, "l", fix [i] + j);
                CHECK_VALUE (ret != 0, "Failed to parse fixed edge value");
                Py_DECREF (v);
            }
            Py_DECREF (s);
        }
        PGAIntegerSetFixedEdges (ctx, len, fix, fixed_edges_symmetric);
        free (fix);
    }
    if (epsilon_generation) {
        PGASetEpsilonGeneration (ctx, epsilon_generation);
    }
    if (epsilon_exponent > 0) {
        PGASetEpsilonExponent (ctx, epsilon_exponent);
    }
    if (epsilon_theta > 0) {
        PGASetEpsilonTheta (ctx, epsilon_theta);
    }
    if (multi_obj_precision > 0) {
        PGASetMultiObjPrecision (ctx, multi_obj_precision);
    }
    if (output_file != NULL) {
        PyObject *b = NULL;
        char *filename = NULL;

        b = PyUnicode_AsEncodedString (output_file, "utf-8", "strict");
        if (!b) {
            return INIT_FAIL;
        }
        filename = PyBytes_AsString (b);
        if (filename == NULL) {
            return INIT_FAIL;
        }
        PGASetOutputFile (ctx, filename);
        Py_DECREF (b);
    }
    if (char_init_type >= 0) {
        CHECK_VALUE
            ( (  char_init_type == PGA_CINIT_UPPER
              || char_init_type == PGA_CINIT_LOWER
              || char_init_type == PGA_CINIT_MIXED
              )
            , "invalid char_init_type"
            );
        PGASetCharacterInitType (ctx, char_init_type);
    }

    PGASetUp (ctx);

    return 0;
}

/***********************
 * Other special methods
 ***********************/

/* __len__ */
Py_ssize_t PGA_len (PyObject *self)
{
    PGAContext *ctx;

    if (!(ctx = get_context (self))) {
        return -1;
    }
    return PGAGetStringLength (ctx);
}


/****************
 * Member Methods
 ****************/

static PyObject *PGA_check_stopping_conditions (PyObject *self, PyObject *args)
{
    PGAContext *ctx;

    if (!(ctx = get_context (self))) {
        return NULL;
    }
    return Py_BuildValue ("i", PGACheckStoppingConditions (ctx));
}

static PyObject *PGA_encode_int_as_binary (PyObject *self, PyObject *args)
{
    PGAContext *ctx = NULL;
    int p, pop, frm, to, val;

    if (!PyArg_ParseTuple (args, "iiiii", &p, &pop, &frm, &to, &val)) {
        return NULL;
    }
    if (!(ctx = get_context (self))) {
        return NULL;
    }
    if (!check_allele (ctx, p, pop, frm)) {
        return NULL;
    }
    if (!check_allele (ctx, p, pop, to)) {
        return NULL;
    }
    CHECK_VALUE_EXCEPTION
        ( PGAGetDataType (ctx) == PGA_DATATYPE_BINARY
        , "Only valid for binary allele"
        , PyExc_ValueError
        , NULL
        );
    PGAEncodeIntegerAsBinary (ctx, p, pop, frm, to, val);
    Py_INCREF   (Py_None);
    return Py_None;
}

static PyObject *PGA_encode_int_as_gray_code (PyObject *self, PyObject *args)
{
    PGAContext *ctx = NULL;
    int p, pop, frm, to, val;

    if (!PyArg_ParseTuple (args, "iiiii", &p, &pop, &frm, &to, &val)) {
        return NULL;
    }
    if (!(ctx = get_context (self))) {
        return NULL;
    }
    if (!check_allele (ctx, p, pop, frm)) {
        return NULL;
    }
    if (!check_allele (ctx, p, pop, to)) {
        return NULL;
    }
    CHECK_VALUE_EXCEPTION
        ( PGAGetDataType (ctx) == PGA_DATATYPE_BINARY
        , "Only valid for binary allele"
        , PyExc_ValueError
        , NULL
        );
    PGAEncodeIntegerAsGrayCode (ctx, p, pop, frm, to, val);
    Py_INCREF   (Py_None);
    return Py_None;
}

static PyObject *PGA_encode_real_as_binary (PyObject *self, PyObject *args)
{
    PGAContext *ctx = NULL;
    int p, pop, frm, to;
    double l, u, val;

    if (!PyArg_ParseTuple (args, "iiiiddd", &p, &pop, &frm, &to, &l, &u, &val))
    {
        return NULL;
    }
    if (!(ctx = get_context (self))) {
        return NULL;
    }
    if (!check_allele (ctx, p, pop, frm)) {
        return NULL;
    }
    if (!check_allele (ctx, p, pop, to)) {
        return NULL;
    }
    CHECK_VALUE_EXCEPTION
        ( PGAGetDataType (ctx) == PGA_DATATYPE_BINARY
        , "Only valid for binary allele"
        , PyExc_ValueError
        , NULL
        );
    PGAEncodeRealAsBinary (ctx, p, pop, frm, to, l, u, val);
    Py_INCREF   (Py_None);
    return Py_None;
}

static PyObject *PGA_encode_real_as_gray_code (PyObject *self, PyObject *args)
{
    PGAContext *ctx = NULL;
    int p, pop, frm, to;
    double l, u, val;

    if (!PyArg_ParseTuple (args, "iiiiddd", &p, &pop, &frm, &to, &l, &u, &val))
    {
        return NULL;
    }
    if (!(ctx = get_context (self))) {
        return NULL;
    }
    if (!check_allele (ctx, p, pop, frm)) {
        return NULL;
    }
    if (!check_allele (ctx, p, pop, to)) {
        return NULL;
    }
    CHECK_VALUE_EXCEPTION
        ( PGAGetDataType (ctx) == PGA_DATATYPE_BINARY
        , "Only valid for binary allele"
        , PyExc_ValueError
        , NULL
        );
    PGAEncodeRealAsGrayCode (ctx, p, pop, frm, to, l, u, val);
    Py_INCREF   (Py_None);
    return Py_None;
}

static PyObject *PGA_euclidian_distance (PyObject *self, PyObject *args)
{
    int p1, pop1, p2, pop2;
    double dist = 0.0;
    PGAContext *ctx;

    if (!PyArg_ParseTuple (args, "iiii", &p1, &pop1, &p2, &pop2)) {
        return NULL;
    }
    if (!(ctx = get_context (self))) {
        return NULL;
    }
    if (PGAGetDataType (ctx) == PGA_DATATYPE_INTEGER) {
        dist = PGAIntegerEuclidianDistance (ctx, p1, pop1, p2, pop2);
    } else if (PGAGetDataType (ctx) == PGA_DATATYPE_REAL) {
        dist = PGARealEuclidianDistance (ctx, p1, pop1, p2, pop2);
    } else {
        PyErr_SetString \
            ( PyExc_NotImplementedError
            , "No euclidian distance method for this data type"
            );
        return NULL;
    }
    return Py_BuildValue ("d", dist);
}

static PyObject *PGA_evaluate (PyObject *self, PyObject *args)
{
    int p, pop;

    if (!PyArg_ParseTuple (args, "ii", &p, &pop)) {
        return NULL;
    }
    PyErr_SetString \
        ( PyExc_NotImplementedError
        , "You must define \"evaluate\" in a derived class"
        );
    return NULL;
}

/* (Re)compute fitness */
static PyObject *PGA_fitness (PyObject *self, PyObject *args)
{
    PGAContext *ctx;
    int pop;

    if (!PyArg_ParseTuple (args, "i", &pop)) {
        return NULL;
    }
    if (!(ctx = get_context (self))) {
        return NULL;
    }
    PGAFitness (ctx, pop);
    Py_INCREF   (Py_None);
    return Py_None;
}

/* Used to set the gene for user defined datatypes.
 * If no user defined datatypes are in use, raise a ValueError.
 * This gets p, pop, userdata, the last being a PyObject (an instance of
 * the user defined datatype)
 */
static PyObject *PGA_set_gene (PyObject *self, PyObject *args)
{
    PGAContext *ctx;
    int p, pop;
    PyObject *gene = NULL;
    PGAIndividual *ind = NULL;
    if (!PyArg_ParseTuple (args, "iiO", &p, &pop, &gene)) {
        return NULL;
    }
    if (!(ctx = get_context (self))) {
        return NULL;
    }
    CHECK_VALUE_EXCEPTION
        ( ctx->ga.datatype == PGA_DATATYPE_USER
        , "This method is used only for user-defined datatypes"
        , PyExc_ValueError
        , NULL
        );
    ind = PGAGetIndividual (ctx, p, pop);
    /* Decrement refcount of stored user data object if existing */
    if (ind->chrom != NULL) {
        Py_DECREF (ind->chrom);
        ind->chrom = NULL;
    }
    ind->chrom = gene;
    Py_INCREF (gene);
    Py_INCREF (Py_None);
    return Py_None;
}

/* Used to retrieve the gene for user defined datatypes.
 * If no user defined datatypes are in use, raise a ValueError.
 * This gets p, pop and returns the gene, a PyObject (an instance of the
 * user defined datatype)
 */
static PyObject *PGA_get_gene (PyObject *self, PyObject *args)
{
    PGAContext *ctx;
    int p, pop;
    PGAIndividual *ind = NULL;
    if (!PyArg_ParseTuple (args, "ii", &p, &pop)) {
        return NULL;
    }
    if (!(ctx = get_context (self))) {
        return NULL;
    }
    CHECK_VALUE_EXCEPTION
        ( ctx->ga.datatype == PGA_DATATYPE_USER
        , "This method is used only for user-defined datatypes"
        , PyExc_ValueError
        , NULL
        );
    ind = PGAGetIndividual (ctx, p, pop);
    CHECK_VALUE_EXCEPTION
        ( ind->chrom != NULL
        , "This gene is not set"
        , PyExc_ValueError
        , NULL
        );
    Py_INCREF (ind->chrom);
    return ind->chrom;
}


/*
 * Get and Set methods.
 */
static PyObject *PGA_get_allele (PyObject *self, PyObject *args)
{
    PGAContext *ctx = NULL;
    int p, pop, i;

    if (!(ctx = get_context (self))) {
        return NULL;
    }
    ERR_CHECK_OCCURRED (ctx, NULL);
    if (!PyArg_ParseTuple (args, "iii", &p, &pop, &i)) {
        return NULL;
    }
    if (!check_allele (ctx, p, pop, i)) {
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
    case PGA_DATATYPE_USER:
    {
        PyErr_SetString (PyExc_ValueError, "No allele for user data type");
        return NULL;
    }
    default:
        assert (0);
    }
    Py_INCREF (Py_None);
    return Py_None;
}

static PyObject *PGA_get_best_index (PyObject *self, PyObject *args)
{
    PGAContext *ctx = NULL;
    int pop;

    if (!PyArg_ParseTuple (args, "i", &pop)) {
        return NULL;
    }
    if (!(ctx = get_context (self))) {
        return NULL;
    }
    return Py_BuildValue ("i", PGAGetBestIndex (ctx, pop));
}

static PyObject *PGA_get_best_report (PyObject *self, PyObject *args)
{
    PGAContext *ctx = NULL;
    int pop, idx;

    if (!PyArg_ParseTuple (args, "ii", &pop, &idx)) {
        return NULL;
    }
    if (!(ctx = get_context (self))) {
        return NULL;
    }
    return Py_BuildValue ("d", PGAGetBestReport (ctx, pop, idx));
}

static PyObject *PGA_get_best_report_index (PyObject *self, PyObject *args)
{
    PGAContext *ctx = NULL;
    int pop, idx;

    if (!PyArg_ParseTuple (args, "ii", &pop, &idx)) {
        return NULL;
    }
    if (!(ctx = get_context (self))) {
        return NULL;
    }
    return Py_BuildValue ("i", PGAGetBestReportIndex (ctx, pop, idx));
}

static PyObject *PGA_get_evaluation (PyObject *self, PyObject *args)
{
    PGAContext *ctx = NULL;
    PyObject *tuple = NULL;
    PyObject *ele = NULL;
    int p, pop, i;
    const double *aux;

    if (!PyArg_ParseTuple (args, "ii", &p, &pop))
    {
        return NULL;
    }
    if (!(ctx = get_context (self))) {
        return NULL;
    }
    CHECK_VALUE_EXCEPTION
        ( PGAGetEvaluationUpToDateFlag (ctx, p, pop)
        , "Evaluation not up to date"
        , PyExc_ValueError
        , NULL
        );
    if (ctx->ga.NumAuxEval == 0) {
        return Py_BuildValue ("d", _PGAGetEvaluation (ctx, p, pop, NULL));
    }
    tuple = PyTuple_New (ctx->ga.NumAuxEval + 1);
    if (tuple == NULL) {
        return NULL;
    }
    ele = Py_BuildValue ("d", PGAGetEvaluation (ctx, p, pop, &aux));
    ERR_DECREF_RET (ele != NULL, tuple, NULL);
    PyTuple_SetItem (tuple, 0, ele);
    for (i=0; i<ctx->ga.NumAuxEval; i++) {
        ele = Py_BuildValue ("d", aux [i]);
        ERR_DECREF_RET (ele != NULL, tuple, NULL);
        PyTuple_SetItem (tuple, i + 1, ele);
    }
    return tuple;
}

static PyObject *PGA_get_evaluation_up_to_date (PyObject *self, PyObject *args)
{
    PGAContext *ctx = NULL;
    int p, pop;

    if (!PyArg_ParseTuple (args, "ii", &p, &pop)) {
        return NULL;
    }
    if (!(ctx = get_context (self))) {
        return NULL;
    }
    return Py_BuildValue ("i", PGAGetEvaluationUpToDateFlag (ctx, p, pop));
}

static PyObject *PGA_get_fitness (PyObject *self, PyObject *args)
{
    PGAContext *ctx = NULL;
    int p, pop;

    if (!PyArg_ParseTuple (args, "ii", &p, &pop)) {
        return NULL;
    }
    if (!(ctx = get_context (self))) {
        return NULL;
    }
    return Py_BuildValue ("d", PGAGetFitness (ctx, p, pop));
}

static PyObject *PGA_get_int_from_binary (PyObject *self, PyObject *args)
{
    PGAContext *ctx = NULL;
    int p, pop, frm, to;

    if (!PyArg_ParseTuple (args, "iiii", &p, &pop, &frm, &to)) {
        return NULL;
    }
    if (!(ctx = get_context (self))) {
        return NULL;
    }
    if (!check_allele (ctx, p, pop, frm)) {
        return NULL;
    }
    if (!check_allele (ctx, p, pop, to)) {
        return NULL;
    }
    CHECK_VALUE_EXCEPTION
        ( PGAGetDataType (ctx) == PGA_DATATYPE_BINARY
        , "Only valid for binary allele"
        , PyExc_ValueError
        , NULL
        );
    return Py_BuildValue ("i", PGAGetIntegerFromBinary (ctx, p, pop, frm, to));
}

static PyObject *PGA_get_int_from_gray_code (PyObject *self, PyObject *args)
{
    PGAContext *ctx = NULL;
    int p, pop, frm, to;

    if (!PyArg_ParseTuple (args, "iiii", &p, &pop, &frm, &to)) {
        return NULL;
    }
    if (!(ctx = get_context (self))) {
        return NULL;
    }
    if (!check_allele (ctx, p, pop, frm)) {
        return NULL;
    }
    if (!check_allele (ctx, p, pop, to)) {
        return NULL;
    }
    CHECK_VALUE_EXCEPTION
        ( PGAGetDataType (ctx) == PGA_DATATYPE_BINARY
        , "Only valid for binary allele"
        , PyExc_ValueError
        , NULL
        );
    return Py_BuildValue
        ("i", PGAGetIntegerFromGrayCode (ctx, p, pop, frm, to));
}

static PyObject *PGA_get_iteration (PyObject *self, PyObject *args)
{
    PGAContext *ctx = NULL;

    if (!(ctx = get_context (self))) {
        return NULL;
    }
    return Py_BuildValue ("i", PGAGetGAIterValue (ctx));
}

static PyObject *PGA_get_real_from_binary (PyObject *self, PyObject *args)
{
    PGAContext *ctx = NULL;
    int p, pop, frm, to;
    double l, u;

    if (!PyArg_ParseTuple (args, "iiiidd", &p, &pop, &frm, &to, &l, &u))
    {
        return NULL;
    }
    if (!(ctx = get_context (self))) {
        return NULL;
    }
    if (!check_allele (ctx, p, pop, frm)) {
        return NULL;
    }
    if (!check_allele (ctx, p, pop, to)) {
        return NULL;
    }
    CHECK_VALUE_EXCEPTION
        ( PGAGetDataType (ctx) == PGA_DATATYPE_BINARY
        , "Only valid for binary allele"
        , PyExc_ValueError
        , NULL
        );
    return Py_BuildValue
        ("d", PGAGetRealFromBinary (ctx, p, pop, frm, to, l, u));
}

static PyObject *PGA_get_real_from_gray_code (PyObject *self, PyObject *args)
{
    PGAContext *ctx = NULL;
    int p, pop, frm, to;
    double l, u;

    if (!PyArg_ParseTuple (args, "iiiidd", &p, &pop, &frm, &to, &l, &u))
    {
        return NULL;
    }
    if (!(ctx = get_context (self))) {
        return NULL;
    }
    if (!check_allele (ctx, p, pop, frm)) {
        return NULL;
    }
    if (!check_allele (ctx, p, pop, to)) {
        return NULL;
    }
    CHECK_VALUE_EXCEPTION
        ( PGAGetDataType (ctx) == PGA_DATATYPE_BINARY
        , "Only valid for binary allele"
        , PyExc_ValueError
        , NULL
        );
    return Py_BuildValue
        ("d", PGAGetRealFromGrayCode (ctx, p, pop, frm, to, l, u));
}

static PyObject *PGA_get_worst_index (PyObject *self, PyObject *args)
{
    PGAContext *ctx = NULL;
    int pop;

    if (!PyArg_ParseTuple (args, "i", &pop)) {
        return NULL;
    }
    if (!(ctx = get_context (self))) {
        return NULL;
    }
    return Py_BuildValue ("i", PGAGetWorstIndex (ctx, pop));
}

/*
 * Python gene print function -- can be overridden in descendent class
 */
static PyObject *PGA_print_string (PyObject *self, PyObject *args)
{
    PyFileObject *file = NULL;
    PGAContext   *ctx = NULL;
    int           p, pop;
    FILE         *fp = NULL;
    int           do_close = 0;
    char         *errmsg = NULL;

    if (!PyArg_ParseTuple (args, "Oii", &file, &p, &pop)) {
        return NULL;
    }
    if (!(ctx = get_context (self))) {
        return NULL;
    }
    if (!(fp = get_fp (self, file))) {
        PyErr_Clear ();
        fp = get_fp_from_file (file);
        if (fp == NULL) {
            return NULL;
        }
        do_close = 1;
    }
    switch (PGAGetDataType (ctx)) {
    case PGA_DATATYPE_BINARY:
        PGABinaryPrintString    (ctx, fp, p, pop);
        break;
    case PGA_DATATYPE_CHARACTER:
        PGACharacterPrintString (ctx, fp, p, pop);
        break;
    case PGA_DATATYPE_INTEGER:
        PGAIntegerPrintString   (ctx, fp, p, pop);
        break;
    case PGA_DATATYPE_REAL:
        PGARealPrintString      (ctx, fp, p, pop);
        break;
    case PGA_DATATYPE_USER:
    {
        PGAIndividual *ind = PGAGetIndividual (ctx, p, pop);
        CHECK_VALUE_EXCEPTION
            ( ind->chrom != NULL
            , "This gene is not set"
            , PyExc_ValueError
            , NULL
            );
        PyObject_Print (ind->chrom, fp, 0);
        break;
    }
    default:
        assert (0);
    }
    if (fflush (fp) != 0) {
        errmsg = "Cannot flush FILE";
    }
    if (do_close) {
        if (fclose (fp) != 0) {
            errmsg = "Cannot close FILE";
        }
    }
    if (errmsg) {
        PyErr_SetString (PyExc_ValueError, errmsg);
        return NULL;
    }
    Py_INCREF (Py_None);
    return Py_None;
}

/*
 * Python context print function -- used for debugging purposes
 * Currently we always print to stderr
 */
static PyObject *PGA_print_context (PyObject *self, PyObject *args)
{
    PGAContext   *ctx = NULL;

    if (!PyArg_ParseTuple (args, "")) {
        return NULL;
    }
    if (!(ctx = get_context (self))) {
        return NULL;
    }
    PGAPrintContextVariable (ctx, stderr);
    Py_INCREF (Py_None);
    return Py_None;
}

static PyObject *PGA_random_01 (PyObject *self, PyObject *args)
{
    PGAContext *ctx = NULL;

    if (!(ctx = get_context (self))) {
        return NULL;
    }
    return Py_BuildValue ("d", PGARandom01 (ctx, 0));
}

static PyObject *PGA_random_flip (PyObject *self, PyObject *args)
{
    PGAContext *ctx = NULL;
    double     probability;

    if (!PyArg_ParseTuple (args, "d", &probability)) {
        return NULL;
    }
    if (!check_probability (probability)) {
        return NULL;
    }
    if (!(ctx = get_context (self))) {
        return NULL;
    }
    return Py_BuildValue ("i", PGARandomFlip (ctx, probability));
}

static PyObject *PGA_random_gaussian (PyObject *self, PyObject *args)
{
    PGAContext *ctx;
    double     l, r;

    if (!PyArg_ParseTuple (args, "dd", &l, &r)) {
        return NULL;
    }
    if (!(ctx = get_context (self))) {
        return NULL;
    }
    return Py_BuildValue ("d", PGARandomGaussian (ctx, l, r));
}

static PyObject *PGA_random_interval (PyObject *self, PyObject *args)
{
    PGAContext *ctx = NULL;
    int        l, r;

    if (!(ctx = get_context (self))) {
        return NULL;
    }
    ERR_CHECK_OCCURRED (ctx, NULL);
    if (!PyArg_ParseTuple (args, "ii", &l, &r)) {
        return NULL;
    }
    check_interval (l, r);
    return Py_BuildValue ("i", PGARandomInterval (ctx, l, r));
}

static PyObject *PGA_random_uniform (PyObject *self, PyObject *args)
{
    PGAContext *ctx;
    double     l, r;

    if (!PyArg_ParseTuple (args, "dd", &l, &r)) {
        return NULL;
    }
    check_interval (l, r);
    if (!(ctx = get_context (self))) {
        return NULL;
    }
    return Py_BuildValue ("d", PGARandomUniform (ctx, l, r));
}

static PyObject *PGA_run (PyObject *self, PyObject *args)
{
    PGAContext *ctx = NULL;

    if (!(ctx = get_context (self))) {
        return NULL;
    }
    PGARun (ctx, evaluate);
    ERR_CHECK_OCCURRED (ctx, NULL);
    Py_INCREF (Py_None);
    return Py_None;
}

static PyObject *PGA_select_next_index (PyObject *self, PyObject *args)
{
    PGAContext *ctx = NULL;
    int pop;

    if (!PyArg_ParseTuple (args, "i", &pop)) {
        return NULL;
    }
    if (!(ctx = get_context (self))) {
        return NULL;
    }
    return Py_BuildValue ("i", PGASelectNextIndex (ctx, pop));
}

static PyObject *PGA_set_allele (PyObject *self, PyObject *args)
{
    PyObject *val = NULL;
    PGAContext *ctx = NULL;
    int p, pop, i;

    if (!(ctx = get_context (self))) {
        return NULL;
    }
    ERR_CHECK_OCCURRED (ctx, NULL);
    if (!PyArg_ParseTuple (args, "iiiO", &p, &pop, &i, &val)) {
        return NULL;
    }
    if (!check_allele (ctx, p, pop, i)) {
        return NULL;
    }

    switch (PGAGetDataType (ctx)) {
    case PGA_DATATYPE_BINARY:
    {
        int allele;
        if (!PyArg_Parse (val, "i", &allele)) {
            return NULL;
        }
        PGASetBinaryAllele (ctx, p, pop, i, allele);
        break;
    }
    case PGA_DATATYPE_CHARACTER:
    {
        char allele;
        if (!PyArg_Parse (val, "c", &allele)) {
            return NULL;
        }
        PGASetCharacterAllele (ctx, p, pop, i, allele);
        break;
    }
    case PGA_DATATYPE_INTEGER:
    {
        int allele;
        if (!PyArg_Parse (val, "i", &allele)) {
            return NULL;
        }
        PGASetIntegerAllele (ctx, p, pop, i, allele);
        break;
    }
    case PGA_DATATYPE_REAL:
    {
        double allele;
        if (!PyArg_Parse (val, "d", &allele)) {
            return NULL;
        }
        PGASetRealAllele (ctx, p, pop, i, allele);
        break;
    }
    case PGA_DATATYPE_USER:
    {
        PyErr_SetString (PyExc_ValueError, "No allele for user data type");
        return NULL;
    }
    default:
        assert (0);
    }
    Py_INCREF   (Py_None);
    return Py_None;
}

static PyObject *PGA_set_evaluation (PyObject *self, PyObject *args)
{
    PGAContext *ctx = NULL;
    int p, pop;
    double val;
    Py_ssize_t nargs, i;
    double *values;
    PyObject *res;

    if (!(ctx = get_context (self))) {
        return NULL;
    }
    nargs = PGAGetNumAuxEval (ctx);
    CHECK_VALUE_EXCEPTION
        ( PySequence_Length (args) == nargs + 3
        , "set_evaluation needs p, pop and the values"
        , PyExc_ValueError
        , NULL
        );
    if ((values = malloc (sizeof (double) * nargs)) == NULL) {
        return PyErr_NoMemory ();
    }
    res = PySequence_GetItem (args, 0);
    if (res == NULL) {
        return NULL;
    }
    if (!PyArg_Parse (res, "i", &p)) {
        return NULL;
    }
    Py_DECREF (res);
    res = PySequence_GetItem (args, 1);
    if (res == NULL) {
        return NULL;
    }
    if (!PyArg_Parse (res, "i", &pop)) {
        return NULL;
    }
    Py_DECREF (res);
    res = PySequence_GetItem (args, 2);
    if (res == NULL) {
        return NULL;
    }
    if (!PyArg_Parse (res, "d", &val)) {
        return NULL;
    }
    Py_DECREF (res);
    for (i=0; i<nargs; i++) {
        res = PySequence_GetItem (args, i + 3);
        if (res == NULL) {
            return NULL;
        }
        if (!PyArg_Parse (res, "d", values + i)) {
            return NULL;
        }
        Py_DECREF (res);
    }
    _PGASetEvaluation (ctx, p, pop, val, values);
    Py_INCREF (Py_None);
    return Py_None;
}

static PyObject *PGA_set_evaluation_up_to_date (PyObject *self, PyObject *args)
{
    PGAContext *ctx = NULL;
    int p, pop, status;

    if (!PyArg_ParseTuple (args, "iii", &p, &pop, &status)) {
        return NULL;
    }
    if (!(ctx = get_context (self))) {
        return NULL;
    }

    PGASetEvaluationUpToDateFlag (ctx, p, pop, status);
    Py_INCREF   (Py_None);
    return Py_None;
}

static PyMethodDef PGA_methods [] =
{ { "check_stopping_conditions", PGA_check_stopping_conditions, METH_VARARGS
  , "Return original stop condition check"
  }
, { "encode_int_as_binary",      PGA_encode_int_as_binary,      METH_VARARGS
  , "Encode int as BCD in binary string"
  }
, { "encode_int_as_gray_code",   PGA_encode_int_as_gray_code,   METH_VARARGS
  , "Encode int as gray code in binary string"
  }
, { "encode_real_as_binary",     PGA_encode_real_as_binary,     METH_VARARGS
  , "Encode real as BCD in binary string"
  }
, { "encode_real_as_gray_code",  PGA_encode_real_as_gray_code,  METH_VARARGS
  , "Encode real as gray code in binary string"
  }
, { "euclidian_distance",        PGA_euclidian_distance,        METH_VARARGS
  , "Get euclidian distance betwee two strings"
  }
, { "evaluate",                  PGA_evaluate,                  METH_VARARGS
  , "Evaluate"
  }
, { "fitness",                   PGA_fitness,                   METH_VARARGS
  , "(Re) compute fitness from evaluations"
  }
, { "get_allele",                PGA_get_allele,                METH_VARARGS
  , "Get allele"
  }
, { "get_best_index",            PGA_get_best_index,            METH_VARARGS
  , "Get best index in population pop"
  }
, { "get_best_report",           PGA_get_best_report,           METH_VARARGS
  , "Get best evaluation for function with given index"
  }
, { "get_best_report_index",     PGA_get_best_report_index,     METH_VARARGS
  , "Get best index for evaluation function with given index"
  }
, { "get_evaluation",            PGA_get_evaluation,            METH_VARARGS
  , "Get evaluation"
  }
, { "get_evaluation_up_to_date", PGA_get_evaluation_up_to_date, METH_VARARGS
  , "Get evaluation up-to-date info"
  }
, { "get_fitness",               PGA_get_fitness,               METH_VARARGS
  , "Get fitness of an individual"
  }
, { "get_gene",                  PGA_get_gene,                  METH_VARARGS
  , "Get gene for user defined datatype"
  }
, { "get_int_from_binary",       PGA_get_int_from_binary,       METH_VARARGS
  , "Get integer value from binary string encoded in BCD"
  }
, { "get_int_from_gray_code",    PGA_get_int_from_gray_code,    METH_VARARGS
  , "Get integer value from binary string encoded in gray code"
  }
, { "get_iteration",             PGA_get_iteration,             METH_VARARGS
  , "Current iteration (GA iter)"
  }
, { "get_real_from_binary",      PGA_get_real_from_binary,      METH_VARARGS
  , "Get real value from binary string encoded in BCD"
  }
, { "get_real_from_gray_code",   PGA_get_real_from_gray_code,   METH_VARARGS
  , "Get real value from binary string encoded in gray code"
  }
, { "get_worst_index",           PGA_get_worst_index,           METH_VARARGS
  , "Get worst index in population pop"
  }
, { "print_context",             PGA_print_context,             METH_VARARGS
  , "Python context print, debug info about PGApack context"
  }
, { "print_string",              PGA_print_string,              METH_VARARGS
  , "Python gene print function -- can be overridden in descendent class."
  }
, { "random01",                  PGA_random_01,                 METH_VARARGS
  , "Random float 0 <= f <= 1"
  }
, { "random_flip",               PGA_random_flip,               METH_VARARGS
  , "Random int 0/1 with probability p"
  }
, { "random_gaussian",           PGA_random_gaussian,           METH_VARARGS
  , "Random value from gaussian distribution with mean, std_deviation"
  }
, { "random_interval",           PGA_random_interval,           METH_VARARGS
  , "Random int [l,r]"
  }
, { "random_uniform",            PGA_random_uniform,            METH_VARARGS
  , "Random float [l,r]"
  }
, { "run",                       PGA_run,                       METH_VARARGS
  , "Run optimization"
  }
, { "select_next_index",         PGA_select_next_index,         METH_VARARGS
  , "Get index of next individual after selection"
  }
, { "set_allele",                PGA_set_allele,                METH_VARARGS
  , "Set allele"
  }
, { "set_evaluation",            PGA_set_evaluation,            METH_VARARGS
  , "Set evaluation"
  }
, { "set_evaluation_up_to_date", PGA_set_evaluation_up_to_date, METH_VARARGS
  , "Set evaluation up to date or not up to date (to True or False)"
  }
, { "set_gene",                  PGA_set_gene,                  METH_VARARGS
  , "Set gene for user defined datatype"
  }
, { NULL } /* EMPTY VALUE AS END-MARKER */
};

/**********************************
 * Getters (and optionally setters)
 **********************************/

#define XSTR(x) #x
#define GETTER_FUNCTION(pganame, attrname, arg) \
    static PyObject *PGA_ ## attrname (PyObject *self, void *closure) \
    {                                                                 \
        PGAContext *ctx;                                              \
        if (!(ctx = get_context (self))) {                            \
            return NULL;                                              \
        }                                                             \
        return Py_BuildValue (XSTR(arg), pganame (ctx));              \
    }

#define ST_DATATYPE_d double
#define ST_DATATYPE_i int
#define ST_CONV_d(x) PyFloat_AsDouble(x)
#define ST_CONV_i(x) PyLong_AsLong(x)
#define ST_DECLARE_TMPVAR(name, arg) ST_DATATYPE_##arg name
#define ST_DATATYPE_CONV(arg, value) ST_CONV_##arg(value)

/* Only to be used with arg=d or arg=i */
#define SETTER_FUNCTION(pganame, attrname, arg, maxval) \
    static int PGA_SET_ ## attrname (PyObject *self, PyObject *v, void *c) \
    {                                                                      \
        PGAContext *ctx;                                                   \
        ST_DECLARE_TMPVAR(tmp, arg);                                       \
        if (!(ctx = get_context (self))) {                                 \
            return -1;                                                     \
        }                                                                  \
        tmp = ST_DATATYPE_CONV(arg, v);                                    \
	if (PyErr_Occurred ()) {                                           \
	    return -1;                                                     \
	}                                                                  \
        CHECK_VALUE_EXCEPTION                                              \
            ( !maxval || tmp <= maxval                                     \
            , "Too large"                                                  \
            , PyExc_ValueError                                             \
            , -1                                                           \
            );                                                             \
        pganame (ctx, tmp);                                                \
        return 0;                                                          \
    }


/* These do *NOT* end in semicolon */
GETTER_FUNCTION (PGAGetCrossoverProb,            crossover_prob,         d)
GETTER_FUNCTION (PGAGetCrossoverBounceBackFlag,  crossover_bounce_back,  i)
GETTER_FUNCTION (PGAGetCrossoverBoundedFlag,     crossover_bounded,      i)
GETTER_FUNCTION (PGAGetCrossoverSBXEta,          crossover_SBX_eta,      d)
GETTER_FUNCTION (PGAGetCrossoverSBXOncePerString,
                 crossover_SBX_once_per_string,                          i)
GETTER_FUNCTION (PGAGetDEVariant,                DE_variant,             i)
GETTER_FUNCTION (PGAGetDENumDiffs,               DE_num_diffs,           i)
GETTER_FUNCTION (PGAGetDEScaleFactor,            DE_scale_factor,        d)
GETTER_FUNCTION (PGAGetDEAuxFactor,              DE_aux_factor,          d)
GETTER_FUNCTION (PGAGetDECrossoverProb,          DE_crossover_prob,      d)
GETTER_FUNCTION (PGAGetDEJitter,                 DE_jitter,              d)
GETTER_FUNCTION (PGAGetDEDither,                 DE_dither,              d)
GETTER_FUNCTION (PGAGetDEDitherPerIndividual,
                 DE_dither_per_individual,                               i)
GETTER_FUNCTION (PGAGetDEProbabilityEO,          DE_probability_EO,      d)
GETTER_FUNCTION (PGAGetEpsilonExponent,          epsilon_exponent,       d)
GETTER_FUNCTION (PGAGetEpsilonGeneration,        epsilon_generation,     i)
GETTER_FUNCTION (PGAGetEpsilonTheta,             epsilon_theta,          i)
GETTER_FUNCTION (PGAGetEvalCount,                eval_count,             i)
GETTER_FUNCTION (PGAGetFitnessCmaxValue,         fitness_cmax,           d)
GETTER_FUNCTION (PGAGetFitnessMinType,           fitness_min_type,       i)
GETTER_FUNCTION (PGAGetFitnessType,              fitness_type,           i)
GETTER_FUNCTION (PGAGetGAIterValue,              GA_iter,                i)
GETTER_FUNCTION (PGAGetMaxFitnessRank,           max_fitness_rank,       d)
GETTER_FUNCTION (PGAGetMaxSimilarityValue,       max_similarity,         i)
GETTER_FUNCTION (PGAGetMaxGAIterValue,           max_GA_iter,            i)
GETTER_FUNCTION (PGAGetMultiObjPrecision,        multi_obj_precision,    i)
GETTER_FUNCTION (PGAGetMutationAndCrossoverFlag, mutation_and_crossover, i)
GETTER_FUNCTION (PGAGetMutationBounceBackFlag,   mutation_bounce_back,   i)
GETTER_FUNCTION (PGAGetMutationBoundedFlag,      mutation_bounded,       i)
GETTER_FUNCTION (PGAGetMutationPolyEta,          mutation_poly_eta,      d)
GETTER_FUNCTION (PGAGetMutationOrCrossoverFlag,  mutation_or_crossover,  i)
GETTER_FUNCTION (PGAGetMutationOnlyFlag,         mutation_only,          i)
GETTER_FUNCTION (PGAGetMutationPolyValue,        mutation_poly_value,    d)
GETTER_FUNCTION (PGAGetMutationProb,             mutation_prob,          d)
GETTER_FUNCTION (PGAGetMutationType,             mutation_type,          i)
GETTER_FUNCTION (PGAGetNumConstraint,            num_constraint,         i)
GETTER_FUNCTION (PGAGetNumReplaceValue,          num_replace,            i)
GETTER_FUNCTION (PGAGetPopSize,                  pop_size,               i)
GETTER_FUNCTION (PGAGetPrintFrequencyValue,      print_frequency,        i)
GETTER_FUNCTION (PGAGetPTournamentProb,          p_tournament_prob,      d)
GETTER_FUNCTION (PGAGetRandomizeSelect,          randomize_select,       i)
GETTER_FUNCTION (PGAGetRandomSeed,               random_seed,            i)
GETTER_FUNCTION (PGAGetRestartFlag,              restart,                i)
GETTER_FUNCTION (PGAGetRestartFrequencyValue,    restart_frequency,      i)
GETTER_FUNCTION (PGAGetRTRWindowSize,            rtr_window_size,        i)
GETTER_FUNCTION (PGAGetStringLength,             string_length,          i)
GETTER_FUNCTION (PGAGetSumConstraintsFlag,       sum_constraints,        i)
GETTER_FUNCTION (PGAGetTournamentSize,           tournament_size,        d)
GETTER_FUNCTION (PGAGetTournamentWithReplacement,tournament_with_replacement,i)
GETTER_FUNCTION (PGAGetTruncationProportion,     truncation_proportion,  d)
GETTER_FUNCTION (PGAGetUniformCrossoverProb,     uniform_crossover_prob, d)

/* These do *NOT* end in semicolon */
/*               pgapack name                pgapy name           type  max */
SETTER_FUNCTION (PGASetCrossoverProb,        crossover_prob,         d,  1)
SETTER_FUNCTION (PGASetEpsilonExponent,      epsilon_exponent,       d, 10)
SETTER_FUNCTION (PGASetMultiObjPrecision,    multi_obj_precision,    i, 14)
SETTER_FUNCTION (PGASetPTournamentProb,      p_tournament_prob,      d,  1)
SETTER_FUNCTION (PGASetUniformCrossoverProb, uniform_crossover_prob, d,  1)

/* This is a special getter because it doesn't have a single data type */
static PyObject *PGA_mutation_value (PyObject *self, void *closure)
{
    PGAContext *ctx;
    if (!(ctx = get_context (self))) {
        return NULL;
    }
    if (PGAGetDataType (ctx) == PGA_DATATYPE_REAL) {
        return Py_BuildValue ("d", PGAGetMutationRealValue (ctx));
    } else if (PGAGetDataType (ctx) == PGA_DATATYPE_INTEGER) {
	return Py_BuildValue ("i", PGAGetMutationIntegerValue (ctx));
    }
    Py_INCREF (Py_None);
    return Py_None;
}

/* Another special getter because we need to increment NumAuxEval */
static PyObject *PGA_num_eval (PyObject *self, void *closure)
{
    PGAContext *ctx;
    if (!(ctx = get_context (self))) {
        return NULL;
    }
    return Py_BuildValue ("i", PGAGetNumAuxEval (ctx) + 1);
}

/* Special getter for maximize */
static PyObject *PGA_maximize (PyObject *self, void *closure)
{
    PGAContext *ctx;
    if (!(ctx = get_context (self))) {
        return NULL;
    }
    return Py_BuildValue ("i", PGAGetOptDirFlag (ctx) == PGA_MAXIMIZE ? 1:0);
}

/* Yet another explicit getter which uses the default mpi communicator */
static PyObject *PGA_mpi_rank (PyObject *self, void *closure)
{
    PGAContext *ctx;
    if (!(ctx = get_context (self))) {
        return NULL;
    }
    return Py_BuildValue ("i", PGAGetRank (ctx, MPI_COMM_WORLD));
}

/* Yet another explicit getter which uses the default mpi communicator */
static PyObject *PGA_mpi_n_proc (PyObject *self, void *closure)
{
    PGAContext *ctx;
    if (!(ctx = get_context (self))) {
        return NULL;
    }
    return Py_BuildValue ("i", PGAGetNumProcs (ctx, MPI_COMM_WORLD));
}


#define GETTER_ENTRY(name) \
    { XSTR(name), PGA_ ## name }
#define GETSET_ENTRY(name) \
    { XSTR(name), PGA_ ## name, PGA_SET_ ## name }

static PyGetSetDef PGA_getset [] =
/*  name      .get                  .set   .doc .closure */
{ GETSET_ENTRY (crossover_prob)
, GETTER_ENTRY (crossover_bounce_back)
, GETTER_ENTRY (crossover_bounded)
, GETTER_ENTRY (crossover_SBX_eta)
, GETTER_ENTRY (crossover_SBX_once_per_string)
, GETTER_ENTRY (DE_variant)
, GETTER_ENTRY (DE_num_diffs)
, GETTER_ENTRY (DE_scale_factor)
, GETTER_ENTRY (DE_aux_factor)
, GETTER_ENTRY (DE_crossover_prob)
, GETTER_ENTRY (DE_jitter)
, GETTER_ENTRY (DE_dither)
, GETTER_ENTRY (DE_dither_per_individual)
, GETTER_ENTRY (DE_probability_EO)
, GETSET_ENTRY (epsilon_exponent)
, GETTER_ENTRY (epsilon_generation)
, GETTER_ENTRY (epsilon_theta)
, GETTER_ENTRY (eval_count)
, GETTER_ENTRY (fitness_cmax)
, GETTER_ENTRY (fitness_min_type)
, GETTER_ENTRY (fitness_type)
, GETTER_ENTRY (GA_iter)
, GETTER_ENTRY (max_fitness_rank)
, GETTER_ENTRY (max_GA_iter)
, GETTER_ENTRY (max_similarity)
, GETTER_ENTRY (maximize)
, GETTER_ENTRY (mpi_n_proc)
, GETTER_ENTRY (mpi_rank)
, GETSET_ENTRY (multi_obj_precision)
, GETTER_ENTRY (mutation_and_crossover)
, GETTER_ENTRY (mutation_or_crossover)
, GETTER_ENTRY (mutation_only)
, GETTER_ENTRY (mutation_poly_eta)
, GETTER_ENTRY (mutation_poly_value)
, GETTER_ENTRY (mutation_bounded)
, GETTER_ENTRY (mutation_bounce_back)
, GETTER_ENTRY (mutation_prob)
, GETTER_ENTRY (mutation_type)
, GETTER_ENTRY (mutation_value)
, GETTER_ENTRY (num_constraint)
, GETTER_ENTRY (num_eval)
, GETTER_ENTRY (num_replace)
, GETTER_ENTRY (pop_size)
, GETTER_ENTRY (print_frequency)
, GETSET_ENTRY (p_tournament_prob)
, GETTER_ENTRY (randomize_select)
, GETTER_ENTRY (random_seed)
, GETTER_ENTRY (restart)
, GETTER_ENTRY (restart_frequency)
, GETTER_ENTRY (rtr_window_size)
, GETTER_ENTRY (string_length)
, GETTER_ENTRY (sum_constraints)
, GETTER_ENTRY (tournament_size)
, GETTER_ENTRY (tournament_with_replacement)
, GETTER_ENTRY (truncation_proportion)
, GETSET_ENTRY (uniform_crossover_prob)
, { NULL }
};

/***************************
 * Creation of the PGA Class
 ***************************/

typedef struct {
    PyObject_HEAD
    PyObject *inst_dict;
    PyObject *weakreflist;
} PGAObject;

static void PGA_dealloc (PyObject *self)
{
    PGAContext *ctx;
    PyObject *Py_MPI_i = NULL;
    PyObject *PGA_ctx = NULL;
    int mpi_initialized = 0;

    ctx = get_context (self);
    #if 0
    fprintf (stderr, "After ctx in PGA_dealloc: %p\n", ctx);
    fflush  (stderr);
    #endif
    if (ctx != NULL) {
        PGA_ctx = PyObject_GetAttrString (self, "context");
        if (PGA_ctx != NULL) {
            /* Ignore return code here, can't do anything if this fails */
            PyObject_DelItem (contexts, PGA_ctx);
            Py_DECREF (PGA_ctx);
        }
        PGADestroy (ctx);
    }

    Py_MPI_i = PyObject_GetAttrString (self, "mpi_initialized");
    if (!Py_MPI_i) {
        return;
    }
    if (!PyArg_Parse (Py_MPI_i, "i", &mpi_initialized)) {
        Py_DECREF   (Py_MPI_i);
        return;
    }
    Py_DECREF (Py_MPI_i);
    /* Only call exitfunc if MPI wasn't initialized externally */
    if (!mpi_initialized) {
        exitfunc ();
    }
    self->ob_type->tp_free (self);
}

static PyObject *PGA_new (PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    PGAObject *self = NULL;
    self = (PGAObject *) type->tp_alloc (type, 0);
    return (PyObject *) self;
}

static PyMappingMethods PGA_mm = {
    .mp_length         = PGA_len,
};

static PyTypeObject PGA_Type = {
    PyVarObject_HEAD_INIT(NULL, 0)
    .tp_name           = "pga.PGA",
    .tp_doc            = "PGApack wrapper class",
    .tp_as_mapping     = &PGA_mm,
    .tp_weaklistoffset = offsetof (PGAObject, weakreflist),
    .tp_dictoffset     = offsetof (PGAObject, inst_dict),
    .tp_basicsize      = sizeof (PGAObject),
    .tp_itemsize       = 0,
    .tp_flags          = ( Py_TPFLAGS_DEFAULT
                         | Py_TPFLAGS_BASETYPE
                         ),
    .tp_new            = PGA_new,
    .tp_init           = (initproc)   PGA_init,
    .tp_dealloc        = (destructor) PGA_dealloc,
    .tp_getset         = PGA_getset,
    .tp_methods        = PGA_methods,
};

/*****************
 * Module creation
 *****************/

#if IS_PY3
static struct PyModuleDef module_definition =
{ PyModuleDef_HEAD_INIT
, "pga"
, NULL
, -1
, Module_Methods
};
// python3 has different convention for module name
#define initpga PyInit_pga
#endif

PyMODINIT_FUNC initpga (void)
{
    PyObject *module = NULL;
    constdef_t  *cd;
    PyObject *module_Dict;
    PyObject *version = PyString_FromString_Compat (VERSION);
    PyObject *weakref = PyImport_ImportModule ("weakref");

    if (!version || !weakref) {
        return FAIL;
    }
    /* We do not want to keep a reference to all PGA objects
     * But we still need to look them up via the context pointer.
     * So use a WeakValueDictionary to store them, this doesn't prevent
     * the destructor from being called.
     */
    contexts = PyObject_CallMethod (weakref, "WeakValueDictionary", "");
    if (!contexts) {
        return FAIL;
    }
    Py_DECREF (weakref);

    if (PyType_Ready (&PGA_Type) < 0) {
        return FAIL;
    }
#if IS_PY3
    module = PyModule_Create (&module_definition);
#else
    module = Py_InitModule   ("pga", Module_Methods);
#endif
    if (module == NULL) {
        return FAIL;
    }
    Py_INCREF (&PGA_Type);
    if (PyModule_AddObject (module, "PGA", (PyObject *) &PGA_Type) < 0) {
        Py_DECREF (&PGA_Type);
        return FAIL;
    }
    if (all_symbols (module) != 0) {
        Py_DECREF (&PGA_Type);
        return FAIL;
    }
    module_Dict = PyModule_GetDict (module);
    PyDict_SetItemString (module_Dict, "contexts", contexts);
    PyDict_SetItemString (module_Dict, "VERSION", version);
    for (cd = constdef; cd->cd_name; cd++) {
        PyObject *constant = Py_BuildValue    ("i", cd->cd_value);
        PyDict_SetItemString (module_Dict, cd->cd_name, constant);
        Py_DECREF (constant);
    }

#if IS_PY3
    return module;
#endif
}
