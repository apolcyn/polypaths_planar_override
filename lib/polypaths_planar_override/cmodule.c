/***************************************************************************
* Copyright (c) 2010 by Casey Duncan
* All rights reserved.
*
* This software is subject to the provisions of the BSD License
* A copy of the license should accompany this distribution.
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
* FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
****************************************************************************/
#include "Python.h"
#include "polypaths_planar_override.h"

double polypaths_planar_override_EPSILON = 1e-5;
double polypaths_planar_override_EPSILON2 = 1e-5 * 1e-5;

static PyObject *
_set_epsilon_func(PyObject *self, PyObject *epsilon)
{
    epsilon = PyObject_ToFloat(epsilon);
    if (epsilon == NULL) {
        return NULL;
    }

    polypaths_planar_override_EPSILON = PyFloat_AS_DOUBLE(epsilon);
    polypaths_planar_override_EPSILON2 = polypaths_planar_override_EPSILON * polypaths_planar_override_EPSILON;
    Py_DECREF(epsilon);
    Py_INCREF(Py_None);
    return Py_None;
}

PyObject *polypaths_planar_overrideTransformNotInvertibleError;

static PyMethodDef module_functions[] = {
    {"_set_epsilon", (PyCFunction) _set_epsilon_func, METH_O,
     "PRIVATE: Set epsilon value used by C extension"},
    {NULL}
};

PyDoc_STRVAR(module_doc, "polypaths_planar_override native code classes");

#define INIT_TYPE(type, name) {                                         \
    if ((type).tp_new == 0) {                                           \
		(type).tp_new = PyType_GenericNew;                              \
    }                                                                   \
    if (PyType_Ready(&(type)) < 0) {                                    \
        goto fail;                                                      \
    }                                                                   \
    if (PyModule_AddObject(module, (name), (PyObject *)&(type)) < 0) {  \
        goto fail;                                                      \
    }                                                                   \
}

#if PY_MAJOR_VERSION >= 3

static struct PyModuleDef moduledef = {
        PyModuleDef_HEAD_INIT,
        "cvector",
        module_doc,
        -1,                 /* m_size */
        module_functions,   /* m_methods */
        NULL,               /* m_reload (unused) */
        NULL,               /* m_traverse */
        NULL,               /* m_clear */
        NULL                /* m_free */
};

#define INITERROR return NULL

PyObject *
PyInit_c(void)

#else
#define INITERROR return

void
initc(void)
#endif
{
#if PY_MAJOR_VERSION >= 3
    PyObject *module = PyModule_Create(&moduledef);
#else
    PyObject *module = Py_InitModule3("c", module_functions, module_doc);
#endif
    Py_INCREF((PyObject *)&polypaths_planar_overrideVec2Type);
    Py_INCREF((PyObject *)&polypaths_planar_overrideSeq2Type);
    Py_INCREF((PyObject *)&polypaths_planar_overrideVec2ArrayType);
    Py_INCREF((PyObject *)&polypaths_planar_overrideAffineType);
    Py_INCREF((PyObject *)&polypaths_planar_overrideBBoxType);
    Py_INCREF((PyObject *)&polypaths_planar_overrideLineType);
    Py_INCREF((PyObject *)&polypaths_planar_overrideRayType);
    Py_INCREF((PyObject *)&polypaths_planar_overrideSegmentType);
    Py_INCREF((PyObject *)&polypaths_planar_overridePolygonType);

    INIT_TYPE(polypaths_planar_overrideVec2Type, "Vec2");
    INIT_TYPE(polypaths_planar_overrideSeq2Type, "Seq2");
    INIT_TYPE(polypaths_planar_overrideVec2ArrayType, "Vec2Array");
	/* Override inheritance of tp_itemsize, ugly */
	polypaths_planar_overrideVec2ArrayType.tp_itemsize = 0;
    INIT_TYPE(polypaths_planar_overrideAffineType, "Affine");
    INIT_TYPE(polypaths_planar_overrideBBoxType, "BoundingBox");
    INIT_TYPE(polypaths_planar_overrideLineType, "Line");
    INIT_TYPE(polypaths_planar_overrideRayType, "Ray");
    INIT_TYPE(polypaths_planar_overrideSegmentType, "LineSegment");
    INIT_TYPE(polypaths_planar_overridePolygonType, "Polygon");

	polypaths_planar_overrideTransformNotInvertibleError = PyErr_NewException(
		"polypaths_planar_override.TransformNotInvertibleError", NULL, NULL);
	if (polypaths_planar_overrideTransformNotInvertibleError == NULL) {
		goto fail;
	}
    if (PyModule_AddObject(
        module, "TransformNotInvertibleError", 
		polypaths_planar_overrideTransformNotInvertibleError) < 0) {
        Py_DECREF(polypaths_planar_overrideTransformNotInvertibleError);
        goto fail;
    }
#if PY_MAJOR_VERSION >= 3
    return module;
#else
    return;
#endif
fail:
    Py_DECREF((PyObject *)&polypaths_planar_overrideVec2Type);
    Py_DECREF((PyObject *)&polypaths_planar_overrideSeq2Type);
    Py_DECREF((PyObject *)&polypaths_planar_overrideVec2ArrayType);
    Py_DECREF((PyObject *)&polypaths_planar_overrideAffineType);
    Py_DECREF((PyObject *)&polypaths_planar_overrideBBoxType);
    Py_DECREF((PyObject *)&polypaths_planar_overrideLineType);
    Py_DECREF((PyObject *)&polypaths_planar_overrideRayType);
    Py_DECREF((PyObject *)&polypaths_planar_overrideSegmentType);
    Py_DECREF((PyObject *)&polypaths_planar_overridePolygonType);
    Py_DECREF(module);
    INITERROR;
}

