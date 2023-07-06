/* Patched version by Berry S. */
/*
 * Copyright 2010 VIFF Development Team.
 *
 * This file is part of VIFF, the Virtual Ideal Functionality Framework.
 *
 * VIFF is free software: you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License (LGPL) as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * VIFF is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 * or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General
 * Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with VIFF. If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include <Python.h>
#include <structmember.h>

#include <gmp.h>
#include <longintrepr.h>

/* be compatible with Python < 2.6 */

#ifndef Py_TYPE
#define Py_TYPE(ob)             (((PyObject*)(ob))->ob_type)
#endif

// dummy type FieldElement
static PyTypeObject cfield_FieldElementType = {
    PyObject_HEAD_INIT(NULL)
    0,                          /*ob_size*/
    "cfield.FieldElement",      /*tp_name*/
    sizeof(PyObject),           /*tp_basicsize*/
    0,                          /*tp_itemsize*/
    0,                          /*tp_dealloc*/
    0,                          /*tp_print*/
    0,                          /*tp_getattr*/
    0,                          /*tp_setattr*/
    0,                          /*tp_compare*/
    0,                          /*tp_repr*/
    0,                          /*tp_as_number*/
    0,                          /*tp_as_sequence*/
    0,                          /*tp_as_mapping*/
    0,                          /*tp_hash */
    0,                          /*tp_call*/
    0,                          /*tp_str*/
    0,                          /*tp_getattro*/
    0,                          /*tp_setattro*/
    0,                          /*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT|Py_TPFLAGS_BASETYPE, /*tp_flags*/
    0,                          /*tp_doc*/
    0,                          /*tp_traverse*/
    0,                          /*tp_clear*/
    0,                          /*tp_richcompare*/
    0,                          /*tp_weaklistoffset*/
    0,                          /*tp_iter*/
    0,                          /*tp_iternext*/
    0,                          /*tp_methods*/
    0,                          /*tp_members*/
    0,                          /*tp_getset*/
    0,                          /*tp_base*/
    0,                          /*tp_dict*/
    0,                          /*tp_descr_get*/
    0,                          /*tp_descr_set*/
    0,                          /*tp_dictoffset*/
    0,                          /*tp_init*/
    0,                          /*tp_alloc*/
    0,                          /*tp_new*/
    0,                          /*tp_free*/
    0,                          /*tp_is_gc*/
    0,                          /*tp_bases*/
    0,                          /*tp_mro*/
    0,                          /*tp_cache*/
    0,                          /*tp_subclasses*/
    0,                          /*tp_weaklist*/
};

typedef struct {
    PyObject_HEAD
    int value;
} cfield_GF256;

static struct PyMemberDef cfield_GF256_members[] = {
  {"value", T_INT, offsetof(cfield_GF256, value), READONLY, 0},
  {0, 0, 0, 0, 0}
};

static cfield_GF256 cfield_GF256_inst[256];
static PyObject* cfield_GF256_mult_table[256][256];
static PyObject* cfield_GF256_inv_table[256];

static void cfield_GF256_generate_tables(void) {
    int log_table[256];
    int exp_table[256];
    unsigned char c;
    unsigned char a = 1;
    unsigned char d;
    int x, y;
    unsigned char z;

    // Code from http://www.samiam.org/galois.html
    for (c = 0; c < 255; c++) {
	exp_table[c] = a;
	/* Multiply by three */
	d = a & 0x80;
	a <<= 1;

	if (d == 0x80)
	    a ^= 0x1b;

	a ^= exp_table[c];
	/* Set the log table value */
	log_table[exp_table[c]] = c;
    }

    exp_table[255] = exp_table[0];
    log_table[0] = 0;

    for (x = 0; x < 256; x++) {
	for (y = 0; y < 256; y++) {
	    if (x == 0 || y == 0)
		z = 0;
	    else
		z = exp_table[(log_table[x] + log_table[y]) % 255];
	    cfield_GF256_mult_table[x][y] =
		(PyObject*)&cfield_GF256_inst[z];
	}

	cfield_GF256_inv_table[x] = (PyObject*)
	    &cfield_GF256_inst[exp_table[255 - log_table[x]]];
    }
}

static PyObject* cfield_GF256__new(unsigned char value) {
    Py_INCREF(&cfield_GF256_inst[value]);
    return (PyObject*)&cfield_GF256_inst[value];
}

static PyObject* cfield_GF256_new(PyTypeObject* type, PyObject* args,
				  PyObject* kwargs)
{
    unsigned char value;

    if (!PyArg_ParseTuple(args, "B", &value))
	return NULL;

    return cfield_GF256__new(value);
}

#define GET_VALUES(IN1, IN2, OUT1, OUT2, OP)				\
    if (Py_TYPE(IN1)->tp_as_number != NULL &&				\
	Py_TYPE(IN1)->tp_as_number->nb_##OP == cfield_GF256_##OP) {	\
	OUT1 = (unsigned char) ((cfield_GF256*)IN1)->value;		\
	if (Py_TYPE(IN2)->tp_as_number != NULL &&			\
	    Py_TYPE(IN2)->tp_as_number->nb_##OP == cfield_GF256_##OP) { \
	    OUT2 = (unsigned char) ((cfield_GF256*)IN2)->value;		\
	} else if (PyInt_Check(IN2))					\
	    OUT2 = (unsigned char) PyInt_AS_LONG(IN2);			\
	else if (PyLong_Check(IN2))					\
	    OUT2 = (unsigned char) PyLong_AsLong(IN2);			\
	else {								\
	    Py_INCREF(Py_NotImplemented);				\
	    return Py_NotImplemented;					\
	}								\
    } else {								\
	if (PyInt_Check(IN1))						\
	    OUT1 = (unsigned char) PyInt_AS_LONG(IN1);			\
	else if (PyLong_Check(IN1))					\
	    OUT1 = (unsigned char) PyLong_AsLong(IN1);			\
	else {								\
	    Py_INCREF(Py_NotImplemented);				\
	    return Py_NotImplemented;					\
	}								\
	OUT2 = (unsigned char) ((cfield_GF256*)IN2)->value;		\
    }

static PyObject* cfield_GF256_add(PyObject* self, PyObject* other) {
    unsigned char a, b;
    GET_VALUES(self, other, a, b, add);
    return cfield_GF256__new(a ^ b);
}

static PyObject* cfield_GF256_multiply(PyObject* self, PyObject* other) {
    unsigned char a, b;
    PyObject* result;
    GET_VALUES(self, other, a, b, multiply);
    result = cfield_GF256_mult_table[a][b];
    Py_INCREF(result);
    return result;
}

static PyObject* cfield_GF256_divide(PyObject* self, PyObject* other) {
    unsigned char a, b;
    PyObject* result;
    GET_VALUES(self, other, a, b, divide);

    if (b == 0) {
	PyErr_SetString(PyExc_ZeroDivisionError, "Cannot divide by zero.");
	return NULL;
    }

    result = cfield_GF256_mult_table[a]
	[((cfield_GF256*)cfield_GF256_inv_table[b])->value];
    Py_INCREF(result);
    return result;
}

static PyObject* cfield_GF256_inv(cfield_GF256* self) {
    PyObject* result;

    if (self->value == 0) {
	PyErr_SetString(PyExc_ZeroDivisionError, "Cannot invert zero.");
	return NULL;
    }

    result = cfield_GF256_inv_table[self->value];
    Py_INCREF(result);
    return result;
}

static cfield_GF256* cfield_GF256_square_multiply(cfield_GF256* base,
						  unsigned char exponent)
{
    cfield_GF256* tmp;

    if (exponent == 0)
	return &cfield_GF256_inst[1];
    else if (exponent % 2 == 0) {
	base = cfield_GF256_square_multiply(base, exponent / 2);
	return (cfield_GF256*)cfield_GF256_mult_table[base->value][base->value];
    } else {
	tmp = cfield_GF256_square_multiply(base, exponent - 1);
	return (cfield_GF256*)cfield_GF256_mult_table[base->value][tmp->value];
    }
}

static PyObject* cfield_GF256_pow(PyObject* self, PyObject* other,
				  PyObject* modulus)
{
    unsigned char exponent;
    cfield_GF256* result;

    if (!PyInt_Check(other)) {
	PyErr_SetString(PyExc_TypeError, "Exponent must be integer.");
	return NULL;
    }

    if (modulus != Py_None) {
	PyErr_SetString(PyExc_TypeError,
			"Exponentiation with modulus not possible.");
	return NULL;
    }

    exponent = (unsigned char)PyInt_AS_LONG(other);
    result = cfield_GF256_square_multiply((cfield_GF256*)self, exponent);
    Py_INCREF(result);
    return (PyObject*)result;
}

static PyObject* cfield_GF256_richcompare(PyObject* self, PyObject* other,
					  int op)
{
    unsigned char a, b;

    switch (op) {
    case Py_EQ:
    case Py_NE:
	if (Py_TYPE(self)->tp_richcompare == cfield_GF256_richcompare) {
	    a = (unsigned char) ((cfield_GF256*)self)->value;
	} else {
	    if (PyInt_Check(self))
		a = (unsigned char) PyInt_AS_LONG(self);
	    else if (PyLong_Check(self))
		a = (unsigned char) PyLong_AsLong(self);
	    else {
		Py_INCREF(Py_NotImplemented);
		return Py_NotImplemented;
	    }

	    b = (unsigned char) ((cfield_GF256*)other)->value;
	}

	if (Py_TYPE(other)->tp_richcompare == cfield_GF256_richcompare) {
	    b = (unsigned char) ((cfield_GF256*)other)->value;
	} else if (PyInt_Check(other))
	    b = (unsigned char) PyInt_AS_LONG(other);
	else if (PyLong_Check(other))
	    b = (unsigned char) PyLong_AsLong(other);
 	else {
	    Py_INCREF(Py_NotImplemented);
	    return Py_NotImplemented;
	}

	if ((a == b) ^ (op == Py_EQ)) {
	    Py_INCREF(Py_False);
	    return Py_False;
	} else {
	    Py_INCREF(Py_True);
	    return Py_True;
	}

	break;
    default:
	PyErr_SetString(PyExc_TypeError, "GF256 elements are not ordered.");
	return NULL;
    }
}

static long cfield_GF256_hash(PyObject* self) {
    return (long)self;
}

static PyObject* cfield_GF256_repr(PyObject* self) {
    return PyString_FromFormat("[%d]", ((cfield_GF256*)self)->value);
}

static int cfield_GF256_nonzero(cfield_GF256* self) {
    if (self->value == 0)
	return 0;
    else
	return 1;
}

static PyObject* cfield_GF256_int(cfield_GF256* self) {
    return PyInt_FromLong(self->value);
}

static PyNumberMethods cfield_GF256_as_number = {
        cfield_GF256_add,               /*nb_add*/
        cfield_GF256_add,               /*nb_subtract*/
        cfield_GF256_multiply,          /*nb_multiply*/
        cfield_GF256_divide,            /*nb_divide*/
        0,                              /*nb_remainder*/
        0,                              /*nb_divmod*/
        cfield_GF256_pow,               /*nb_power*/
        0,                              /*nb_negative*/
        0,                              /*tp_positive*/
        0,                              /*tp_absolute*/
        (inquiry)cfield_GF256_nonzero,  /*tp_nonzero*/
        (unaryfunc)cfield_GF256_inv,    /*nb_invert*/
        0,                              /*nb_lshift*/
        0,                              /*nb_rshift*/
        0,                              /*nb_and*/
	cfield_GF256_add,               /*nb_xor*/
        0,                              /*nb_or*/
        0,                              /*nb_coerce*/
        (unaryfunc)cfield_GF256_int,    /*nb_int*/
        0,                              /*nb_long*/
        0,                              /*nb_float*/
        0,                              /*nb_oct*/
        0,                              /*nb_hex*/
        0,                              /* nb_inplace_add */
        0,                              /* nb_inplace_subtract */
        0,                              /* nb_inplace_multiply */
        0,                              /* nb_inplace_divide */
        0,                              /* nb_inplace_remainder */
        0,                              /* nb_inplace_power */
        0,                              /* nb_inplace_lshift */
        0,                              /* nb_inplace_rshift */
        0,                              /* nb_inplace_and */
        0,                              /* nb_inplace_xor */
        0,                              /* nb_inplace_or */
	cfield_GF256_divide,            /* nb_floor_divide */
        cfield_GF256_divide,            /* nb_true_divide */
        0,                              /* nb_inplace_floor_divide */
        0                               /* nb_inplace_true_divide */
};

static PyTypeObject cfield_GF256Type = {
    PyObject_HEAD_INIT(NULL)
    0,                          /*ob_size*/
    "cfield.GF256",             /*tp_name*/
    sizeof(cfield_GF256),       /*tp_basicsize*/
    0,                          /*tp_itemsize*/
    0,                          /*tp_dealloc*/
    0,                          /*tp_print*/
    0,                          /*tp_getattr*/
    0,                          /*tp_setattr*/
    0,                          /*tp_compare*/
    cfield_GF256_repr,          /*tp_repr*/
    &cfield_GF256_as_number,    /*tp_as_number*/
    0,                          /*tp_as_sequence*/
    0,                          /*tp_as_mapping*/
    cfield_GF256_hash,          /*tp_hash */
    0,                          /*tp_call*/
    cfield_GF256_repr,          /*tp_str*/
    0,                          /*tp_getattro*/
    0,                          /*tp_setattro*/
    0,                          /*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT|Py_TPFLAGS_BASETYPE|Py_TPFLAGS_CHECKTYPES, /*tp_flags*/
    0,                          /*tp_doc*/
    0,                          /*tp_traverse*/
    0,                          /*tp_clear*/
    cfield_GF256_richcompare,   /*tp_richcompare*/
    0,                          /*tp_weaklistoffset*/
    0,                          /*tp_iter*/
    0,                          /*tp_iternext*/
    0,                          /*tp_methods*/
    cfield_GF256_members,       /*tp_members*/
    0,                          /*tp_getset*/
    &cfield_FieldElementType,   /*tp_base*/
    0,                          /*tp_dict*/
    0,                          /*tp_descr_get*/
    0,                          /*tp_descr_set*/
    0,                          /*tp_dictoffset*/
    0,                          /*tp_init*/
    0,                          /*tp_alloc*/
    cfield_GF256_new,           /*tp_new*/
    0,                          /*tp_free*/
    0,                          /*tp_is_gc*/
    0,                          /*tp_bases*/
    0,                          /*tp_mro*/
    0,                          /*tp_cache*/
    0,                          /*tp_subclasses*/
    0,                          /*tp_weaklist*/
};


// Code stolen from gmpy.
static void cfield_long2mpz(PyObject* obj, mpz_t z) {
    mpz_t digit;
    int len, i, negative;
    PyLongObject *l = (PyLongObject *) obj;

    assert(PyLong_Check(obj));

    mpz_init_set_si(z, 0);
    mpz_init(digit);

    if(l->ob_size < 0) {
        len = - l->ob_size;
        negative = 1;
    } else {
        len = l->ob_size;
        negative = 0;
    }
    for(i=0; i<len; i++) {
        mpz_set_ui(digit, l->ob_digit[i]);
        mpz_mul_2exp(digit, digit, i * SHIFT);
        mpz_ior(z, z, digit);
    }
    if(negative)
        mpz_neg(z, z);
    mpz_clear(digit);
}

static PyObject* cfield_mpz2long(mpz_t z) {
    int negative;
    int size;
    int i;
    PyLongObject *newob;
    mpz_t temp;

    /* Assume gmp uses limbs as least as large as the builtin longs do */
    assert(mp_bits_per_limb >= SHIFT);

    mpz_init(temp);
    if(mpz_sgn(z) < 0) {
        negative = 1;
        mpz_neg(temp, z);
    } else {
        negative = 0;
        mpz_set(temp, z);
    }

    size = (mpz_sizeinbase(temp, 2) + SHIFT - 1) / SHIFT;

    if(!(newob = _PyLong_New(size))) {
        mpz_clear(temp);
        return NULL;
    }
    for(i=0; i<size; i++) {
        newob->ob_digit[i] = (unsigned int)(mpz_get_ui(temp) & MASK);
        mpz_fdiv_q_2exp(temp, temp, SHIFT);
    }
    assert(mpz_sgn(temp) == 0);
    mpz_clear(temp);

    /* long_normalize() is file-static so we must reimplement it */
    /* longobjp = long_normalize(longobjp); */
    i = size;
    while ( (i>0) && (newob->ob_digit[i-1] == 0))
        i--;
    newob->ob_size = i;

    if(negative)
        newob->ob_size = - newob->ob_size;
    return (PyObject *) newob;
}

// This breaks binary compatibility with Python versions that change
// PyTypeObject.
typedef struct {
    PyTypeObject type;
    mpz_t modulus;
} cfield_GFType;

typedef struct {
    PyObject_HEAD
    mpz_t value;
} cfield_GFElement;

static int cfield_py2mpz(PyObject* o, mpz_t z) {
    if (PyInt_Check(o)) {
	mpz_init_set_si(z, PyInt_AS_LONG(o));
    } else if (PyLong_Check(o))
	cfield_long2mpz(o, z);
    else
	return -1;

    return 0;
}

static PyObject* cfield_GFElement_new(PyTypeObject* type, PyObject* args,
				      PyObject* kwargs)
{
    cfield_GFElement* self;

    if (PyTuple_GET_SIZE(args) != 1) {
	PyErr_Format(PyExc_TypeError,
		     "GFELement takes exactly one argument (%zd given)",
		     PyTuple_GET_SIZE(args));
	return NULL;
    }

    self = (cfield_GFElement*)type->tp_alloc(type, 0);

    if(cfield_py2mpz(PyTuple_GET_ITEM(args, 0), self->value) < 0) {
	PyErr_SetString(PyExc_TypeError,
			"GFElement parameter must be a number.");
	return NULL;
    }

    mpz_mod(self->value, self->value, ((cfield_GFType*)Py_TYPE(self))->modulus);
    return (PyObject*)self;
}

static void cfield_GFElement_dealloc(PyObject* o) {
    cfield_GFElement* self = (cfield_GFElement*)o;
    mpz_clear(self->value);
    Py_TYPE(o)->tp_free(o);
}

#define GF_ELEMENT_GET_VALUES(IN1, IN2, OUT1, OUT2, OP, TYPE)		\
    mpz_t tmp1, tmp2;							\
    if (Py_TYPE(IN1)->tp_as_number != NULL &&				\
	Py_TYPE(IN1)->tp_as_number->nb_##OP == cfield_GFElement_##OP) { \
	OUT1 = &((cfield_GFElement*)IN1)->value;			\
	TYPE = Py_TYPE(IN1);						\
	if (Py_TYPE(IN2)->tp_as_number != NULL &&			\
	    Py_TYPE(IN2)->tp_as_number->nb_##OP == cfield_GFElement_##OP) { \
	    OUT2 = &((cfield_GFElement*)IN2)->value;			\
	    if (Py_TYPE(IN1) != Py_TYPE(IN2)) {				\
		Py_INCREF(Py_NotImplemented);				\
		return Py_NotImplemented;				\
	    }								\
	} else {							\
	    OUT2 = &tmp2;						\
	    if (cfield_py2mpz(IN2, *OUT2) < 0) {			\
		Py_INCREF(Py_NotImplemented);				\
		return Py_NotImplemented;				\
	    }								\
	}								\
    } else {								\
	OUT1 = &tmp1;							\
	if (cfield_py2mpz(IN1, *OUT1) < 0) {				\
	    Py_INCREF(Py_NotImplemented);				\
	    return Py_NotImplemented;					\
	}								\
	OUT2 = &((cfield_GFElement*)IN2)->value;			\
	TYPE = Py_TYPE(IN2);						\
    }

#define GF_ELEMENT_VALUE_FREE(OUT1, OUT2)	\
    if (OUT1 == &tmp1)				\
	mpz_clear(tmp1);			\
    else if (OUT2 == &tmp2)			\
	mpz_clear(tmp2);

static PyObject* cfield_GFElement_add(PyObject* self, PyObject* other) {
    mpz_t *a, *b;
    cfield_GFElement* result;
    PyTypeObject* type;
    GF_ELEMENT_GET_VALUES(self, other, a, b, add, type);
    result = PyObject_New(cfield_GFElement, type);
    mpz_init(result->value);
    mpz_add(result->value, *a, *b);
    mpz_mod(result->value, result->value, ((cfield_GFType*)type)->modulus);
    GF_ELEMENT_VALUE_FREE(a, b);
    return (PyObject*)result;
}

static PyObject* cfield_GFElement_subtract(PyObject* self, PyObject* other) {
    mpz_t *a, *b;
    cfield_GFElement* result;
    PyTypeObject* type;
    GF_ELEMENT_GET_VALUES(self, other, a, b, subtract, type);
    result = PyObject_New(cfield_GFElement, type);
    mpz_init(result->value);
    mpz_sub(result->value, *a, *b);
    mpz_mod(result->value, result->value, ((cfield_GFType*)type)->modulus);
    GF_ELEMENT_VALUE_FREE(a, b);
    return (PyObject*)result;
}

static PyObject* cfield_GFElement_multiply(PyObject* self, PyObject* other) {
    mpz_t *a, *b;
    cfield_GFElement* result;
    PyTypeObject* type;
    GF_ELEMENT_GET_VALUES(self, other, a, b, multiply, type);
    result = PyObject_New(cfield_GFElement, type);
    mpz_init(result->value);
    mpz_mul(result->value, *a, *b);
    mpz_mod(result->value, result->value, ((cfield_GFType*)type)->modulus);
    GF_ELEMENT_VALUE_FREE(a, b);
    return (PyObject*)result;
}

static PyObject* cfield_GFElement_divide(PyObject* self, PyObject* other) {
    mpz_t *a, *b;
    cfield_GFElement* result;
    PyTypeObject* type;
    GF_ELEMENT_GET_VALUES(self, other, a, b, divide, type);
    result = PyObject_New(cfield_GFElement, type);
    mpz_init(result->value);

    if(mpz_invert(result->value, *b, ((cfield_GFType*)type)->modulus) == 0) {
	PyErr_SetString(PyExc_ZeroDivisionError,
			"Cannot divide by zero.");
	return NULL;
    }

    mpz_mul(result->value, *a, result->value);
    mpz_mod(result->value, result->value, ((cfield_GFType*)type)->modulus);
    GF_ELEMENT_VALUE_FREE(a, b);
    return (PyObject*)result;
}

static PyObject* cfield_GFElement_inv(PyObject* self) {
    cfield_GFElement* result = PyObject_New(cfield_GFElement, Py_TYPE(self));
    mpz_init(result->value);

    if(mpz_invert(result->value, ((cfield_GFElement*)self)->value,
		  ((cfield_GFType*)Py_TYPE(self))->modulus) == 0) {
	PyErr_SetString(PyExc_ZeroDivisionError,
			"Cannot invert zero.");
	return NULL;
    }

    return (PyObject*)result;
}

static PyObject* cfield_GFElement_pow(PyObject* self, PyObject* other,
				      PyObject* modulus)
{
    cfield_GFElement* result;

    if (modulus != Py_None) {
	PyErr_SetString(PyExc_TypeError,
			"Exponentiation with modulus not possible.");
	return NULL;
    }

    result = PyObject_New(cfield_GFElement, Py_TYPE(self));

    if (PyInt_Check(other)) {
	if (PyInt_AS_LONG(other) < 0) {
	    PyErr_SetString(PyExc_AssertionError,
			    "Cannot take power to negative exponent.");
	    return NULL;
	}

	mpz_init(result->value);
	mpz_powm_ui(result->value, ((cfield_GFElement*)self)->value,
		    PyInt_AS_LONG(other),
		    ((cfield_GFType*)Py_TYPE(self))->modulus);
    } else if (PyLong_Check(other)) {
	cfield_long2mpz(other, result->value);
	mpz_powm(result->value, ((cfield_GFElement*)self)->value,
		 result->value, ((cfield_GFType*)Py_TYPE(self))->modulus);
    } else {
	PyErr_SetString(PyExc_TypeError, "Exponent must be integer.");
	return NULL;	
    }

    return (PyObject*)result;
}

static PyObject* cfield_GFElement_xor(PyObject* self, PyObject* other) {
    mpz_t *a, *b;
    cfield_GFElement* result;
    PyTypeObject* type;
    GF_ELEMENT_GET_VALUES(self, other, a, b, xor, type);
    result = PyObject_New(cfield_GFElement, type);
    mpz_init(result->value);
    mpz_xor(result->value, *a, *b);
    GF_ELEMENT_VALUE_FREE(a, b);
    return (PyObject*)result;
}

static PyObject* cfield_GFElement_neg(PyObject* self) {
    cfield_GFElement* result = PyObject_New(cfield_GFElement, Py_TYPE(self));
    mpz_init(result->value);
    mpz_neg(result->value, ((cfield_GFElement*)self)->value);
    return (PyObject*)result;
}

static PyObject* cfield_GFElement_int(PyObject* self) {
    return cfield_mpz2long(((cfield_GFElement*)self)->value);
}

static int cfield_GFElement_nonzero(PyObject* self) {
    int result = mpz_sgn(((cfield_GFElement*)self)->value);
    return result * result;
}

static PyObject* cfield_GFElement_richcompare(PyObject* self, PyObject* other,
					      int op)
{
    mpz_t *a, *b;
    mpz_t tmp_a, tmp_b;
    int c;
    PyObject* result;

    if (Py_TYPE(self)->tp_richcompare == cfield_GFElement_richcompare) {
	a = &((cfield_GFElement*)self)->value;

	if (Py_TYPE(other)->tp_richcompare == cfield_GFElement_richcompare) {
	    b = &((cfield_GFElement*)other)->value;

	    if (Py_TYPE(self) != Py_TYPE(other)) {
		PyErr_SetString(PyExc_AssertionError,
				"Cannot compare elements of different fields.");
		return NULL;
	    }
	} else {
	    b = &tmp_b;

	    if (cfield_py2mpz(other, *b) < 0) {
		Py_INCREF(Py_NotImplemented);
		return Py_NotImplemented;
	    }
	}
    } else {
	a = &tmp_a;

	if (cfield_py2mpz(self, *a) < 0) {
	    Py_INCREF(Py_NotImplemented);
	    return Py_NotImplemented;
	}

	b = &((cfield_GFElement*)other)->value;
    }

    c = mpz_cmp(*a, *b);

    if (a == &tmp_a)
	mpz_clear(*a);
    else if (b == &tmp_b)
	mpz_clear(*b);

    // code from object.c
    switch (op) {
    case Py_LT: c = c <  0; break;
    case Py_LE: c = c <= 0; break;
    case Py_EQ: c = c == 0; break;
    case Py_NE: c = c != 0; break;
    case Py_GT: c = c >  0; break;
    case Py_GE: c = c >= 0; break;
    default:
	Py_INCREF(Py_NotImplemented);
	return Py_NotImplemented;
    }

    result = c ? Py_True : Py_False;
    Py_INCREF(result);
    return result;
}

static long cfield_GFElement_hash(PyObject* self) {
    long result = mpz_get_si(((cfield_GFElement*)self)->value);
    result += mpz_get_si(((cfield_GFType*)Py_TYPE(self))->modulus) * 1000003L;

    if (result == -1)
	return -2;
    else
	return result;
}

static PyObject* cfield_GFElement_repr(PyObject* self) {
    char* str;
    PyObject* result;
    str = mpz_get_str(NULL, 10, ((cfield_GFElement*)self)->value);
    result = PyString_FromFormat("{%s}", str);
    free(str);
    return result;
}

static PyObject* cfield_GFElement_value(PyObject* self, void* noarg) {
    return cfield_mpz2long(((cfield_GFElement*)self)->value);
}

static PyObject* cfield_GFElement_sqrt(PyObject* self, PyObject* noarg) {
    mpz_t* modulus = &((cfield_GFType*)Py_TYPE(self))->modulus;
    char* str;
    cfield_GFElement* result;

    if (mpz_congruent_ui_p(*modulus, 3, 4) == 0) {
	str = mpz_get_str(NULL, 10, *modulus);
	PyErr_Format(PyExc_AssertionError, "Cannot compute square root "
		     "with modulus %s.", str);
	free(str);
	return NULL;
    }

    result = PyObject_New(cfield_GFElement, Py_TYPE(self));
    mpz_init(result->value);
    mpz_add_ui(result->value, *modulus, 1);
    mpz_divexact_ui(result->value, result->value, 4);
    mpz_powm(result->value, ((cfield_GFElement*)self)->value,
	     result->value, *modulus);
    return (PyObject*)result;
}

static PyObject* cfield_GFElement_bit(PyObject* self, PyObject* index) {
    if (!PyInt_Check(index)) {
	PyErr_SetString(PyExc_TypeError, "Index must be an integer");
	return NULL;
    }

    return PyInt_FromLong(mpz_tstbit(((cfield_GFElement*)self)->value,
				     PyInt_AS_LONG(index)));
}

static PyObject* cfield_GFElement_signed(PyObject* self, PyObject* noarg) {
    mpz_t* value = &((cfield_GFElement*)self)->value;
    mpz_t* modulus = &((cfield_GFType*)Py_TYPE(self))->modulus;
    mpz_t tmp;
    PyObject* result;

    mpz_init(tmp);
    mpz_sub_ui(tmp, *modulus, 1);
    mpz_divexact_ui(tmp, tmp, 2);

    if (mpz_cmp(*value, tmp) > 0) {
	mpz_sub(tmp, *value, *modulus);
	result = cfield_mpz2long(tmp);
    } else
	result = cfield_GFElement_value(self, NULL);

    mpz_clear(tmp);
    return result;
}

static PyMethodDef cfield_GFElement_methods[] = {
    {"sqrt", cfield_GFElement_sqrt, METH_NOARGS, 0},
    {"bit", cfield_GFElement_bit, METH_O, 0},
    {"signed", cfield_GFElement_signed, METH_NOARGS, 0},
    {"unsigned", (PyCFunction)cfield_GFElement_value, METH_NOARGS, 0},
    {0, 0, 0, 0}
};

static PyGetSetDef cfield_GFElement_getset[] = {
    {"value", cfield_GFElement_value, 0, 0, 0},
    {0, 0, 0, 0, 0}
};

static PyNumberMethods cfield_GFElement_as_number = {
        cfield_GFElement_add,           /*nb_add*/
        cfield_GFElement_subtract,      /*nb_subtract*/
        cfield_GFElement_multiply,      /*nb_multiply*/
        cfield_GFElement_divide,        /*nb_divide*/
        0,                              /*nb_remainder*/
        0,                              /*nb_divmod*/
        cfield_GFElement_pow,           /*nb_power*/
        cfield_GFElement_neg,           /*nb_negative*/
        0,                              /*tp_positive*/
        0,                              /*tp_absolute*/
        cfield_GFElement_nonzero,       /*tp_nonzero*/
        cfield_GFElement_inv,           /*nb_invert*/
        0,                              /*nb_lshift*/
        0,                              /*nb_rshift*/
        0,                              /*nb_and*/
	cfield_GFElement_xor,           /*nb_xor*/
        0,                              /*nb_or*/
        0,                              /*nb_coerce*/
        cfield_GFElement_int,           /*nb_int*/
        cfield_GFElement_int,           /*nb_long*/
        0,                              /*nb_float*/
        0,                              /*nb_oct*/
        0,                              /*nb_hex*/
        0,                              /* nb_inplace_add */
        0,                              /* nb_inplace_subtract */
        0,                              /* nb_inplace_multiply */
        0,                              /* nb_inplace_divide */
        0,                              /* nb_inplace_remainder */
        0,                              /* nb_inplace_power */
        0,                              /* nb_inplace_lshift */
        0,                              /* nb_inplace_rshift */
        0,                              /* nb_inplace_and */
        0,                              /* nb_inplace_xor */
        0,                              /* nb_inplace_or */
        cfield_GFElement_divide,        /* nb_floor_divide */
        cfield_GFElement_divide,        /* nb_true_divide */
        0,                              /* nb_inplace_floor_divide */
        0                               /* nb_inplace_true_divide */
};

static PyTypeObject cfield_GFTypeType = {
    PyObject_HEAD_INIT(NULL)
    0,                          /*ob_size*/
    "cfield.GFType",            /*tp_name*/
    sizeof(cfield_GFType),      /*tp_basicsize*/
    0,                          /*tp_itemsize*/
    0,                          /*tp_dealloc*/
    0,                          /*tp_print*/
    0,                          /*tp_getattr*/
    0,                          /*tp_setattr*/
    0,                          /*tp_compare*/
    0,                          /*tp_repr*/
    0,                          /*tp_as_number*/
    0,                          /*tp_as_sequence*/
    0,                          /*tp_as_mapping*/
    0,                          /*tp_hash */
    0,                          /*tp_call*/
    0,                          /*tp_str*/
    0,                          /*tp_getattro*/
    0,                          /*tp_setattro*/
    0,                          /*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT,         /*tp_flags*/
    0,                          /*tp_doc*/
    0,                          /*tp_traverse*/
    0,                          /*tp_clear*/
    0,                          /*tp_richcompare*/
    0,                          /*tp_weaklistoffset*/
    0,                          /*tp_iter*/
    0,                          /*tp_iternext*/
    0,                          /*tp_methods*/
    0,                          /*tp_members*/
    0,                          /*tp_getset*/
    0,  /*   &PyType_Type,  Berry S   */          /*tp_base*/
    0,                          /*tp_dict*/
    0,                          /*tp_descr_get*/
    0,                          /*tp_descr_set*/
    0,                          /*tp_dictoffset*/
    0,                          /*tp_init*/
    0,                          /*tp_alloc*/
    PyType_GenericNew,          /*tp_new*/
    0,                          /*tp_free*/
    0,                          /*tp_is_gc*/
    0,                          /*tp_bases*/
    0,                          /*tp_mro*/
    0,                          /*tp_cache*/
    0,                          /*tp_subclasses*/
    0,                          /*tp_weaklist*/
};

static PyObject* cfield_field_cache;

static PyObject* cfield_GF(PyObject* self, PyObject* modulus) {
    cfield_GFType* field;

    if (PyDict_Contains(cfield_field_cache, modulus)) {
	field = (cfield_GFType*)PyDict_GetItem(cfield_field_cache, modulus);
	Py_INCREF((PyObject*)field);
	return (PyObject*)field;
    }

    field = PyObject_New(cfield_GFType, &cfield_GFTypeType);

    memset(&field->type.ob_size, 0, sizeof(PyTypeObject) - sizeof(PyObject));
    field->type.tp_name = "cfield.GFElement";
    field->type.tp_basicsize = sizeof(cfield_GFElement);
    field->type.tp_dealloc = cfield_GFElement_dealloc;
    field->type.tp_repr = cfield_GFElement_repr;
    field->type.tp_as_number = &cfield_GFElement_as_number;
    field->type.tp_hash = cfield_GFElement_hash;
    field->type.tp_str = cfield_GFElement_repr;
    field->type.tp_flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE |
	Py_TPFLAGS_CHECKTYPES;
    field->type.tp_richcompare = cfield_GFElement_richcompare;
    field->type.tp_base = &cfield_FieldElementType;
    field->type.tp_methods = cfield_GFElement_methods;
    field->type.tp_getset = cfield_GFElement_getset;
    field->type.tp_new = cfield_GFElement_new;

    // Have to change type for PyType_Ready() to make mro_internal() happy.
    field->type.ob_type = &PyType_Type;
    if (PyType_Ready(&field->type) < 0)
	return NULL;
    field->type.ob_type = &cfield_GFTypeType;

    PyDict_SetItemString(field->type.tp_dict, "modulus", modulus);
    PyDict_SetItemString(field->type.tp_dict, "field", (PyObject*)field);

    cfield_py2mpz(modulus, field->modulus);

    if (!field->modulus)
	return NULL;

    if (mpz_probab_prime_p(field->modulus, 25) == 0) {
	PyErr_SetString(PyExc_ValueError, "Modulus must be prime.");
	return NULL;
    }

    if (PyDict_SetItem(cfield_field_cache, modulus, (PyObject*)field) < 0)
	return NULL;

    return (PyObject*)field;
}

static PyMethodDef cfield_methods[] = {
    {"GF", cfield_GF, METH_O, 0},
    {0, 0, 0, 0}
};

PyMODINIT_FUNC initcfield(void) {
    PyObject* m;
    int i;

    if (PyType_Ready(&cfield_GF256Type) < 0)
	return;

    m = Py_InitModule3("cfield", cfield_methods,
		       "C implemementation of viff.field");

    if (!m)
	return;

    Py_INCREF(&cfield_FieldElementType);
    PyModule_AddObject(m, "FieldElement", (PyObject*)&cfield_FieldElementType);
    Py_INCREF(&cfield_GF256Type);
    PyModule_AddObject(m, "GF256", (PyObject*)&cfield_GF256Type);

    PyDict_SetItemString(cfield_GF256Type.tp_dict, "modulus",
			 PyInt_FromLong(256));
    PyDict_SetItemString(cfield_GF256Type.tp_dict, "field",
			 (PyObject*)&cfield_GF256Type);

    for (i = 0; i < 256; i++) {
	PyObject_Init((PyObject*)(&cfield_GF256_inst[i]),
		      (PyTypeObject*)&cfield_GF256Type);
	cfield_GF256_inst[i].value = i;
    }

    cfield_GF256_generate_tables();

    cfield_GFTypeType.tp_base = &PyType_Type;  /* Berry S */
    PyType_Ready(&cfield_GFTypeType);

    cfield_field_cache = PyDict_New();
    PyDict_SetItem(cfield_field_cache, PyInt_FromLong(256),
		   (PyObject*)&cfield_GF256Type);
    PyModule_AddObject(m, "_field_cache", cfield_field_cache);
}
