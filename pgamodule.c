/* Copyright (C) 2005 Dr. Ralf Schlatterbeck Open Source Consulting.
 * Reichergasse 131, A-3411 Weidling.
 * Web: http://www.runtux.com Email: office@runtux.com
 * All rights reserved
 * ****************************************************************************
 *
 * This library is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Library General Public License as
 * published by the Free Software Foundation; either version 2 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Library General Public License for more details.
 *
 * You should have received a copy of the GNU Library General Public
 * License along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 * ****************************************************************************
 */

#include <Python.h>
#include <pgapack.h>
#include <stdio.h>
#undef NDEBUG
#include <assert.h>
#include <Version.h>

static PyObject *context        = NULL;
static int       error_occurred = 0;

typedef struct
{
    char *cd_name;
    int   cd_value;
} constdef_t;

static constdef_t constdef [] =
    { {"PGA_NEWPOP",               PGA_NEWPOP                }
    , {"PGA_OLDPOP",               PGA_OLDPOP                }
    , {"PGA_POPREPL_BEST",         PGA_POPREPL_BEST          }
    , {"PGA_POPREPL_RANDOM_NOREP", PGA_POPREPL_RANDOM_NOREP  }
    , {"PGA_POPREPL_RANDOM_REP",   PGA_POPREPL_RANDOM_REP    }
    , {"PGA_REPORT_AVERAGE",       PGA_REPORT_AVERAGE        }
    , {"PGA_REPORT_HAMMING",       PGA_REPORT_HAMMING        }
    , {"PGA_REPORT_OFFLINE",       PGA_REPORT_OFFLINE        }
    , {"PGA_REPORT_ONLINE",        PGA_REPORT_ONLINE         }
    , {"PGA_REPORT_STRING",        PGA_REPORT_STRING         }
    , {"PGA_REPORT_WORST",         PGA_REPORT_WORST          }
    , {"PGA_STOP_MAXITER",         PGA_STOP_MAXITER          }
    , {"PGA_STOP_NOCHANGE",        PGA_STOP_NOCHANGE         }
    , {"PGA_STOP_TOOSIMILAR",      PGA_STOP_TOOSIMILAR       }
    , {NULL,                       0                         }
    };

int compare_constdef (const void *v1, const void *v2)
{
    const constdef_t *p1 = v1, *p2 = v2;
    return strcmp (p1->cd_name, p2->cd_name);
}

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
        error_occurred = 1;      \
        return r;                \
    }                            \
} while (0)

/*
 * Used if the calling object has an endofgen method.
 * Otherwise use built-in default for the datatype.
 * Called after each generation.
 * Do something useful, display the population on a graphics output,
 * let the user adjust the population, etc.
 */
static void endofgen (PGAContext *ctx)
{
    PyObject *self, *r;
    ERR_CHECK (!error_occurred,);
    self = get_self (ctx);
    ERR_CHECK (self, );
    r    = PyObject_CallMethod (self, "endofgen", "");
    ERR_CHECK (r, );
    Py_DECREF (r);
    Py_DECREF (self);
    return;
}

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

/*
 * Used if the calling object has a duplicate-check method.
 * Otherwise use built-in default for the datatype.
 * Perform duplicate checking, compare strings p1 and p2 to see if they
 * are different. If they are, return non-zero, else return 0.
 */
static int check_duplicate (PGAContext *ctx, int p1, int pop1, int p2, int pop2)
{
    PyObject *self, *r;
    int rr, retval;
    ERR_CHECK (!error_occurred, 0);
    self = get_self (ctx);
    ERR_CHECK (self, 0);
    r    = PyObject_CallMethod
        (self, "check_duplicate", "iiii", p1, pop1, p2, pop2);
    ERR_CHECK (r, 0);
    rr = PyArg_Parse (r, "i", &retval);
    ERR_CHECK (rr, 0);
    Py_DECREF (r);
    Py_DECREF (self);
    return !!retval;
}

/*
 * Check stopping criteria, this is always active.
 * User can set a stop_cond method to add stopping criteria.
 * We perform the magic here that if during one of the callback
 * functions (calling into python, e.g. evaluate) an error occurrs, we
 * check the error flag here and stop. This way we can return an error
 * in one of the callback functions to python and raise the appropriate
 * exception there.
 */
static int check_stop (PGAContext *ctx)
{
    PyObject *self;
    ERR_CHECK (!error_occurred, PGA_TRUE);
    self = get_self (ctx);
    ERR_CHECK (self, PGA_TRUE);
    if (PyObject_HasAttrString (self, "stop_cond"))
    {
        int retval, rr;
        PyObject *r = PyObject_CallMethod (self, "stop_cond", NULL);
        ERR_CHECK (r, PGA_TRUE);
        rr = PyArg_Parse (r, "i", &retval);
        ERR_CHECK (rr, PGA_TRUE);
        Py_DECREF (r);
        Py_DECREF (self);
        return !!retval;
    }
    Py_DECREF (self);
    return PGACheckStoppingConditions (ctx);
}

/*
 * Used if the calling object has a crossover method.
 * Otherwise use built-in default for the datatype.
 * Perform crossover from p1 and p2 into c1 and c2
 */
static void crossover
    (PGAContext *ctx, int p1, int p2, int p_pop, int c1, int c2, int c_pop)
{
    PyObject *self, *r;
    ERR_CHECK (!error_occurred, );
    self = get_self (ctx);
    ERR_CHECK (self, );
    r    = PyObject_CallMethod
        (self, "crossover", "iiiiii", p1, p2, p_pop, c1, c2, c_pop);
    ERR_CHECK (r, );
    Py_DECREF (r);
    Py_DECREF (self);
    return;
}

/*
 * Used if the calling object has an initstring method.
 * Otherwise use built-in default for the datatype.
 */
static void initstring (PGAContext *ctx, int p, int pop)
{
    PyObject *self, *r;
    ERR_CHECK (!error_occurred, );
    self = get_self (ctx);
    ERR_CHECK (self, );
    r    = PyObject_CallMethod (self, "initstring", "ii", p, pop);
    ERR_CHECK (r, );
    Py_DECREF (r);
    Py_DECREF (self);
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
    PyObject *self, *r;
    int retval, rr;
    ERR_CHECK (!error_occurred, 0);
    self = get_self (ctx);
    ERR_CHECK (self, 0);
    r    = PyObject_CallMethod (self, "mutation", "iid", p, pop, mr);
    ERR_CHECK (r, 0);
    rr = PyArg_Parse (r, "i", &retval);
    ERR_CHECK (rr, 0);
    Py_DECREF (r);
    Py_DECREF (self);
    return retval;
}

/*
 * Low-level gene print function
 */
static void print_gene (PGAContext *ctx, FILE *fp, int p, int pop)
{
    PyObject *self, *file, *r;
    ERR_CHECK (!error_occurred,);
    self    = get_self (ctx);
    ERR_CHECK (self,);
    file = PyFile_FromFile (fp, "<PGA_file>", "w", fclose);
    ERR_CHECK (file,);
    r    = PyObject_CallMethod (self, "print_string", "Oii", file, p, pop);
    ERR_CHECK (r,);
    Py_DECREF (r);
    Py_DECREF (self);
}

/*
 * Python gene print function -- can be overridden in descendent class
 */
static PyObject *PGA_print_string (PyObject *self0, PyObject *args)
{
    PyObject     *self;
    PyFileObject *file;
    PGAContext   *ctx;
    int           p, pop;

    if (!PyArg_ParseTuple(args, "OOii", &self, &file, &p, &pop))
        return NULL;
    if (!(ctx = get_context (self)))
        return NULL;
    if (!PyFile_Check (file))
        return NULL;
    assert (file->f_fp);
    switch (PGAGetDataType (ctx)) {
    case PGA_DATATYPE_BINARY :
        PGABinaryPrintString    (ctx, file->f_fp, p, pop);
        break;
    case PGA_DATATYPE_CHARACTER :
        PGACharacterPrintString (ctx, file->f_fp, p, pop);
        break;
    case PGA_DATATYPE_INTEGER :
        PGAIntegerPrintString   (ctx, file->f_fp, p, pop);
        break;
    case PGA_DATATYPE_REAL :
        PGARealPrintString      (ctx, file->f_fp, p, pop);
        break;
    default :
        assert (0);
    }
    Py_INCREF (Py_None);
    return Py_None;
}

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
    if (sequence)
    {
        int i, len = PySequence_Length (sequence);
        if (len < 0)
        {
            return 0;
        }
        for (i = 0; i < len; i++)
        {
            PyObject *x = PySequence_GetItem (sequence, i);
            int val, ci;
            if (!x)
                return 0;
            if (!PyArg_Parse (x, "i", &val))
                return 0;
            for (ci = 0; ci < constlen; ci++)
            {
                if (val == constcheck [ci].cd_value)
                    break;
            }
            if (ci == constlen)
            {
                PyErr_SetString 
                    (PyExc_ValueError, name);
                return 0;
            }
            if (checkfun && !checkfun (ctx, val, name))
                return 0;
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
    PyObject *print_options = NULL;
    PyObject *no_duplicates = NULL;
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
        , "print_options"
        , "no_duplicates"
        , NULL
        };

    if  (!PyArg_ParseTupleAndKeywords 
            ( args
            , kw
            , "OOi|OiOOiiiidOiiOO"
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
            , &print_options
            , &no_duplicates
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
    if (PyObject_HasAttrString (self, "check_duplicate"))
    {
        PGASetUserFunction
            (ctx, PGA_USERFUNCTION_DUPLICATE, (void *)check_duplicate);
    }
    if (PyObject_HasAttrString (self, "crossover"))
    {
        PGASetUserFunction (ctx, PGA_USERFUNCTION_CROSSOVER, (void *)crossover);
    }
    if (PyObject_HasAttrString (self, "endofgen"))
    {
        PGASetUserFunction (ctx, PGA_USERFUNCTION_ENDOFGEN, (void *)endofgen);
    }
    if (PyObject_HasAttrString (self, "initstring"))
    {
        PGASetUserFunction
            (ctx, PGA_USERFUNCTION_INITSTRING, (void *)initstring);
    }
    if (PyObject_HasAttrString (self, "mutation"))
    {
        PGASetUserFunction (ctx, PGA_USERFUNCTION_MUTATION, (void *)mutation);
    }
    PGASetUserFunction (ctx, PGA_USERFUNCTION_PRINTSTRING, (void *)print_gene);
    PGASetUserFunction (ctx, PGA_USERFUNCTION_STOPCOND,    (void *)check_stop);
    if (max_no_change)
    {
        PGASetMaxNoChangeValue (ctx, max_no_change);
    }
    if (max_similarity)
    {
        PGASetMaxSimilarityValue (ctx, max_similarity);
    }
    if (mutation_prob)
    {
        PGASetMutationProb (ctx, mutation_prob);
    }
    if (no_duplicates && PyObject_IsTrue (no_duplicates))
    {
        PGASetNoDuplicatesFlag (ctx, PGA_TRUE);
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
    if (random_seed)
    {
        PGASetRandomSeed   (ctx, random_seed);
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
    if (max_GA_iter) 
    {
        if (max_GA_iter <= 2)
        {
            PyErr_SetString \
                ( PyExc_ValueError
                , "Iteration count must be at least 2"
                );
            return NULL;
        }
        PGASetMaxGAIterValue (ctx, max_GA_iter);
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
        return NULL;
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
        return NULL;
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
    case PGA_DATATYPE_BINARY :
    {
        int allele = PGAGetBinaryAllele (ctx, p, pop, i);
        return Py_BuildValue ("i", allele);
        break;
    }
    case PGA_DATATYPE_CHARACTER :
    {
        char allele = PGAGetCharacterAllele (ctx, p, pop, i);
        return Py_BuildValue ("c", allele);
        break;
    }
    case PGA_DATATYPE_INTEGER :
    {
        int allele = PGAGetIntegerAllele (ctx, p, pop, i);
        return Py_BuildValue ("i", allele);
        break;
    }
    case PGA_DATATYPE_REAL :
    {
        double allele = PGAGetRealAllele (ctx, p, pop, i);
        return Py_BuildValue ("d", allele);
        break;
    }
    default :
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
    case PGA_DATATYPE_BINARY :
    {
        int allele;
        if (!PyArg_Parse (val, "i", &allele))
            return NULL;
        PGASetBinaryAllele (ctx, p, pop, i, allele);
        break;
    }
    case PGA_DATATYPE_CHARACTER :
    {
        char allele;
        if (!PyArg_Parse (val, "c", &allele))
            return NULL;
        PGASetCharacterAllele (ctx, p, pop, i, allele);
        break;
    }
    case PGA_DATATYPE_INTEGER :
    {
        int allele;
        if (!PyArg_Parse (val, "i", &allele))
            return NULL;
        PGASetIntegerAllele (ctx, p, pop, i, allele);
        break;
    }
    case PGA_DATATYPE_REAL :
    {
        double allele;
        if (!PyArg_Parse (val, "d", &allele))
            return NULL;
        PGASetRealAllele (ctx, p, pop, i, allele);
        break;
    }
    default :
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
, { "print_string",              PGA_print_string,              METH_VARARGS
  , "Python gene print function -- can be overridden in descendent class."
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

static PyMethodDef Module_Methods[] = { {NULL, NULL, 0, NULL} };

PyMODINIT_FUNC initpga (void)
{
    PyMethodDef *def;
    constdef_t  *cd;
    PyObject *module      = Py_InitModule       ("pga", Module_Methods);
    PyObject *module_Dict = PyModule_GetDict    (module);
    PyObject *class_Dict  = PyDict_New          ();
    PyObject *class_Name  = PyString_FromString ("PGA");
    PyObject *version     = PyString_FromString (VERSION);
    PyObject *pga         = PyString_FromString ("pga");
    PyObject *pga_Class   = PyClass_New         (NULL, class_Dict, class_Name);
    context               = Py_BuildValue       ("{}");
    PyDict_SetItemString (module_Dict, "PGA",        pga_Class);
    PyDict_SetItemString (module_Dict, "context",    context);
    PyDict_SetItemString (module_Dict, "VERSION",    version);
    PyDict_SetItemString (class_Dict,  "__module__", pga);

    for (cd = constdef; cd->cd_name; cd++) {
        PyObject *constant = Py_BuildValue    ("i", cd->cd_value);
        PyDict_SetItemString (module_Dict, cd->cd_name, constant);
        Py_DECREF (constant);
    }

    Py_DECREF(class_Dict);
    Py_DECREF(class_Name);
    Py_DECREF(version);
    Py_DECREF(pga);
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
