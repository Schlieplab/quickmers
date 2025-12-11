/*
Filename: quickmers/_cbindings.c
Author: Kian Jalilian
Copyright: 2025, Alexander Schliep
Version: 0.1.0
Description: Python C extension providing Hamming and Levenshtein distance functions.
License: LGPL-3.0-or-later
*/
#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <numpy/arrayobject.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "hamming.h"
#include "levenshtein.h"

// ---------------- Helper: convert list or numpy array to const uint8_t** ----------------
static int convert_to_kmer_array(PyObject* obj, const uint8_t*** kmers_out, int64_t* n_out, int64_t* kmer_len_out) {
    Py_ssize_t n = 0;
    const uint8_t** kmers = NULL;
    int64_t kmer_len = 0;

    if (PyList_Check(obj)) {
        n = PyList_Size(obj);

        if (n == 0) {
            *kmers_out = NULL;
            *n_out = 0;
            *kmer_len_out = 0;
            return 0;
        }

        kmers = malloc(sizeof(uint8_t*) * n);
        if (!kmers) return -1;

        for (Py_ssize_t i = 0; i < n; i++) {
            PyObject* item = PyList_GetItem(obj, i);
            if (!PyUnicode_Check(item)) { free(kmers); return -1; }
            kmers[i] = (const uint8_t*)PyUnicode_AsUTF8(item);
            if (i == 0) kmer_len = strlen((const char*)kmers[0]);
        }

    } 
    else if (PyArray_Check(obj)) {
        PyArrayObject* arr = (PyArrayObject*)obj;
        if (PyArray_TYPE(arr) != NPY_OBJECT || PyArray_NDIM(arr) != 1) return -1;

        n = PyArray_DIM(arr, 0);

        if (n == 0) {
            *kmers_out = NULL;
            *n_out = 0;
            *kmer_len_out = 0;
            return 0;
        }

        kmers = malloc(sizeof(uint8_t*) * n);
        if (!kmers) return -1;

        for (Py_ssize_t i = 0; i < n; i++) {
            PyObject* item = *(PyObject**)PyArray_GETPTR1(arr, i);
            if (!PyUnicode_Check(item)) { free(kmers); return -1; }
            kmers[i] = (const uint8_t*)PyUnicode_AsUTF8(item);
            if (i == 0) kmer_len = strlen((const char*)kmers[0]);
        }
    } 
    else {
        return -1;
    }

    *kmers_out = kmers;
    *n_out = n;
    *kmer_len_out = kmer_len;
    return 0;
}

// ---------------- 64-bit Hamming Wrappers ------------------

static PyObject* py_hamming_distance_array(PyObject* self, PyObject* args) {
    const char* query;
    PyObject* kmers_obj;
    if (!PyArg_ParseTuple(args, "sO", &query, &kmers_obj)) return NULL;

    const uint8_t** kmers = NULL;
    int64_t n = 0, kmer_len = 0;
    if (convert_to_kmer_array(kmers_obj, &kmers, &n, &kmer_len) != 0) {
        PyErr_SetString(PyExc_TypeError, "Invalid list or numpy array");
        return NULL;
    }

    int* distances = malloc(n * sizeof(int));
    if (!distances) { free(kmers); return PyErr_NoMemory(); }

    hamming_distance_array_64bit(query, (const char**)kmers, n, distances);

    npy_intp dims[1] = { n };
    PyObject* np_array = PyArray_SimpleNewFromData(1, dims, NPY_INT32, distances);
    PyObject* result_copy = PyArray_NewCopy((PyArrayObject*)np_array, NPY_CORDER);

    free(distances);
    free(kmers);
    Py_DECREF(np_array);
    return result_copy;
}

static PyObject* py_hamming_distance(PyObject* self, PyObject* args) {
    const char* s1;
    const char* s2;
    if (!PyArg_ParseTuple(args, "ss", &s1, &s2)) return NULL;

    int64_t dist = hamming_distance_64bit(s1, s2);
    return PyLong_FromLongLong(dist);
}

// ---------------- Levenshtein Wrappers ---------------------

static PyObject* py_levenshtein(PyObject* self, PyObject* args) {
    const char* s1;
    const char* s2;
    if (!PyArg_ParseTuple(args, "ss", &s1, &s2)) return NULL;

    int64_t dist = myers((uint8_t*)s1, strlen(s1), (uint8_t*)s2, strlen(s2));
    return PyLong_FromLongLong(dist);
}

// ---------------- Levenshtein List -------------------------

static PyObject* py_levenshtein_array(PyObject* self, PyObject* args) {
    const char* query;
    PyObject* kmers_obj;
    if (!PyArg_ParseTuple(args, "sO", &query, &kmers_obj)) return NULL;

    const uint8_t** kmers = NULL;
    int64_t n = 0, kmer_len = 0;
    if (convert_to_kmer_array(kmers_obj, &kmers, &n, &kmer_len) != 0) {
        PyErr_SetString(PyExc_TypeError, "Invalid list or numpy array");
        return NULL;
    }

    int64_t* out = malloc(sizeof(int64_t) * n);
    if (!out) { free(kmers); return PyErr_NoMemory(); }

    myers_dispatch((const uint8_t*)query, strlen(query), kmers, n, kmer_len, out);

    npy_intp dims[1] = { n };
    PyObject* np_array = PyArray_SimpleNew(1, dims, NPY_INT64);
    memcpy(PyArray_DATA((PyArrayObject*)np_array), out, n * sizeof(int64_t));

    free(kmers);
    free(out);
    return np_array;
}

// ---------------- Levenshtein List with Min Dist -----------------

static PyObject* py_levenshtein_array_with_min_dist(PyObject* self, PyObject* args) {
    const char* kmer;
    PyObject* kmers_obj;
    long min_distance;

    // Parse arguments: string, object (list or numpy array), integer
    if (!PyArg_ParseTuple(args, "sOl", &kmer, &kmers_obj, &min_distance)) return NULL;

    const uint8_t** kmers = NULL;
    int64_t n = 0, kmer_len = 0;
    if (convert_to_kmer_array(kmers_obj, &kmers, &n, &kmer_len) != 0) {
        PyErr_SetString(PyExc_TypeError, "Invalid list or numpy array");
        return NULL;
    }

    if (n == 0) {
        PyObject* result_list = PyList_New(0);
        PyObject* result_tuple = PyTuple_New(2);
        PyTuple_SetItem(result_tuple, 0, PyBool_FromLong(1));
        PyTuple_SetItem(result_tuple, 1, result_list);
        free(kmers);
        return result_tuple;
    }

    PyObject* result_list = PyList_New(n);
    if (!result_list) { free(kmers); return NULL; }

    int early_exit = 0;
    int64_t len_kmer = (int64_t)strlen(kmer);

    for (Py_ssize_t i = 0; i < n; i++) {
        const char* current_kmer = (const char*)kmers[i];
        int64_t len_current = (int64_t)strlen(current_kmer);

        int64_t dist = myers((uint8_t*)kmer, len_kmer, (uint8_t*)current_kmer, len_current);

        PyObject* py_dist = PyLong_FromLongLong(dist);
        PyList_SetItem(result_list, i, py_dist);

        if (dist < min_distance) {
            early_exit = 1;
            for (Py_ssize_t j = i + 1; j < n; j++) {
                Py_INCREF(Py_None);
                PyList_SetItem(result_list, j, Py_None);
            }
            break;
        }
    }

    PyObject* result_tuple = PyTuple_New(2);
    PyTuple_SetItem(result_tuple, 0, PyBool_FromLong(!early_exit));
    PyTuple_SetItem(result_tuple, 1, result_list);

    free(kmers);
    return result_tuple;
}

// ---------------- Module Definition -----------------------

static PyMethodDef QuickmersMethods[] = {
    {"hamming_distance_array", py_hamming_distance_array, METH_VARARGS, "64-bit Hamming distance for list."},
    {"hamming_distance", py_hamming_distance, METH_VARARGS, "64-bit Hamming distance scalar."},
    {"levenshtein_distance", py_levenshtein, METH_VARARGS, "Levenshtein scalar."},
    {"levenshtein_distance_array", py_levenshtein_array, METH_VARARGS, "Levenshtein list (list or numpy array)."},
    {"levenshtein_distance_array_with_min_dist", py_levenshtein_array_with_min_dist, METH_VARARGS, "Levenshtein list with min distance (list or numpy array)."},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef quickmersmodule = {
    PyModuleDef_HEAD_INIT,
    "_cbindings",
    NULL,
    -1,
    QuickmersMethods
};

PyMODINIT_FUNC PyInit__cbindings(void) {
    PyObject* module = PyModule_Create(&quickmersmodule);
    if (!module) return NULL;
    import_array();
    return module;
}
