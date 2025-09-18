#define PY_SSIZE_T_CLEAN
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include "hamming.h"
#include "levenshtein.h"
#include <Python.h>
#include <numpy/arrayobject.h>

static PyObject* py_hamming_distance_array_32bit(PyObject* self, PyObject* args) {
    const char* query;
    PyObject* seq_list_obj;

    if (!PyArg_ParseTuple(args, "sO!", &query, &PyList_Type, &seq_list_obj))
        return NULL;

    size_t n = PyList_Size(seq_list_obj);

    // Allocate C array of char* for passing to C function
    const char** seq_list = malloc(n * sizeof(char*));
    if (!seq_list) return PyErr_NoMemory();

    for (size_t i = 0; i < n; i++) {
        PyObject* item = PyList_GetItem(seq_list_obj, i);  // borrowed ref
        if (!PyUnicode_Check(item)) {
            free(seq_list);
            PyErr_SetString(PyExc_TypeError, "List elements must be strings");
            return NULL;
        }
        seq_list[i] = PyUnicode_AsUTF8(item);
    }

    // Allocate output array
    int* distances = malloc(n * sizeof(int));
    if (!distances) {
        free(seq_list);
        return PyErr_NoMemory();
    }

    // Call the actual C function
    hamming_distance_array_32bit(query, seq_list, n, distances);

    // Wrap into NumPy array
    npy_intp dims[1] = {n};
    PyObject* result = PyArray_SimpleNewFromData(1, dims, NPY_INT32, distances);
    PyObject* result_copy = PyArray_NewCopy((PyArrayObject*)result, NPY_CORDER);

    // Clean up
    free(distances);
    free(seq_list);
    Py_DECREF(result);

    return result_copy;
}

static PyObject* py_hamming_distance_array_64bit(PyObject* self, PyObject* args) {
    const char* query;
    PyObject* seq_list_obj;

    if (!PyArg_ParseTuple(args, "sO!", &query, &PyList_Type, &seq_list_obj))
        return NULL;

    size_t n = PyList_Size(seq_list_obj);

    // Allocate C array of char* for passing to C function
    const char** seq_list = malloc(n * sizeof(char*));
    if (!seq_list) return PyErr_NoMemory();

    for (size_t i = 0; i < n; i++) {
        PyObject* item = PyList_GetItem(seq_list_obj, i);  // borrowed ref
        if (!PyUnicode_Check(item)) {
            free(seq_list);
            PyErr_SetString(PyExc_TypeError, "List elements must be strings");
            return NULL;
        }
        seq_list[i] = PyUnicode_AsUTF8(item);
    }

    // Allocate output array
    int* distances = malloc(n * sizeof(int));
    if (!distances) {
        free(seq_list);
        return PyErr_NoMemory();
    }

    // Call the actual C function
    hamming_distance_array_64bit(query, seq_list, n, distances);

    // Wrap into NumPy array
    npy_intp dims[1] = {n};
    PyObject* result = PyArray_SimpleNewFromData(1, dims, NPY_INT32, distances);
    PyObject* result_copy = PyArray_NewCopy((PyArrayObject*)result, NPY_CORDER);

    // Clean up
    free(distances);
    free(seq_list);
    Py_DECREF(result);

    return result_copy;
}

static PyObject* py_hamming_distance_encoded_array_32bit(PyObject* self, PyObject* args) {
    uint32_t kmer;
    PyObject* arrayObj;

    if (!PyArg_ParseTuple(args, "IO!", &kmer, &PyArray_Type, &arrayObj))
        return NULL;

    PyArrayObject* arr = (PyArrayObject*)arrayObj;
    uint32_t* kmers_list = (uint32_t*)PyArray_DATA(arr);
    size_t n = PyArray_SIZE(arr);

    int* distances = malloc(n * sizeof(int));
    hamming_distance_encoded_array_32bit(kmer, kmers_list, n, distances);

    npy_intp dims[1] = {n};
    PyObject* result = PyArray_SimpleNewFromData(1, dims, NPY_INT32, distances);
    // Optionally make a copy so we can free distances
    PyObject* result_copy = PyArray_NewCopy((PyArrayObject*)result, NPY_CORDER);
    free(distances);
    Py_DECREF(result);
    return result_copy;
}

static PyObject* py_hamming_distance_encoded_array_64bit(PyObject* self, PyObject* args) {
    uint64_t kmer;
    PyObject* arrayObj;

    if (!PyArg_ParseTuple(args, "KO!", &kmer, &PyArray_Type, &arrayObj))
        return NULL;

    PyArrayObject* arr = (PyArrayObject*)arrayObj;
    uint64_t* kmers_list = (uint64_t*)PyArray_DATA(arr);
    size_t n = PyArray_SIZE(arr);

    int* distances = malloc(n * sizeof(int));
    hamming_distance_encoded_array_64bit(kmer, kmers_list, n, distances);

    npy_intp dims[1] = {n};
    PyObject* result = PyArray_SimpleNewFromData(1, dims, NPY_INT32, distances);
    // Optionally make a copy so we can free distances
    PyObject* result_copy = PyArray_NewCopy((PyArrayObject*)result, NPY_CORDER);
    free(distances);
    Py_DECREF(result);
    return result_copy;
}

static PyObject* py_hamming_distance_encoded_32bit(PyObject* self, PyObject* args) {
    uint32_t kmer1, kmer2;

    // Parse two unsigned ints from Python
    if (!PyArg_ParseTuple(args, "II", &kmer1, &kmer2))
        return NULL;

    // Call the C function
    uint32_t dist = hamming_distance_encoded_32bit(kmer1, kmer2);

    // Return as Python int
    return PyLong_FromUnsignedLong(dist);
}

static PyObject* py_hamming_distance_32bit(PyObject* self, PyObject* args) {
    const char* kmer1;
    const char* kmer2;

    // Parse two strings from Python
    if (!PyArg_ParseTuple(args, "ss", &kmer1, &kmer2))
        return NULL;

    uint32_t dist = hamming_distance_32bit(kmer1, kmer2);

    return PyLong_FromUnsignedLong(dist);
}

static PyObject* py_hamming_distance_encoded_64bit(PyObject* self, PyObject* args) {
    PyObject *kmer1_obj, *kmer2_obj;

    // Accept any Python object
    if (!PyArg_ParseTuple(args, "OO", &kmer1_obj, &kmer2_obj))
        return NULL;

    // Convert to uint64_t
    uint64_t kmer1 = PyLong_AsUnsignedLongLong(kmer1_obj);
    if (PyErr_Occurred()) return NULL;
    uint64_t kmer2 = PyLong_AsUnsignedLongLong(kmer2_obj);
    if (PyErr_Occurred()) return NULL;

    uint64_t dist = hamming_distance_encoded_64bit(kmer1, kmer2);

    return PyLong_FromUnsignedLongLong(dist);
}

static PyObject* py_hamming_distance_64bit(PyObject* self, PyObject* args) {
    const char* kmer1;
    const char* kmer2;

    // Parse two strings from Python
    if (!PyArg_ParseTuple(args, "ss", &kmer1, &kmer2))
        return NULL;

    uint64_t dist = hamming_distance_64bit(kmer1, kmer2);

    return PyLong_FromUnsignedLongLong(dist);
}

static PyObject* py_levenshtein(PyObject* self, PyObject* args) {
    const char* s1;
    const char* s2;

    // Parse two strings from Python
    if (!PyArg_ParseTuple(args, "ss", &s1, &s2))
        return NULL;

    int64_t len1 = (int64_t)strlen(s1);
    int64_t len2 = (int64_t)strlen(s2);

    // Call Myers bit-parallel Levenshtein function
    int64_t dist = myers((uint8_t*)s1, len1, (uint8_t*)s2, len2);

    // Return as Python int
    return PyLong_FromLongLong(dist);
}

static PyObject* py_levenshtein_list(PyObject* self, PyObject* args) {
    const char* kmer;
    PyObject* kmer_list;

    // Parse a string and a Python list
    if (!PyArg_ParseTuple(args, "sO!", &kmer, &PyList_Type, &kmer_list))
        return NULL;

    Py_ssize_t list_size = PyList_Size(kmer_list);
    if (list_size < 0) return NULL; // error checking

    int64_t len_kmer = (int64_t)strlen(kmer);

    // Create Python list to store results
    PyObject* result_list = PyList_New(list_size);
    if (!result_list) return NULL;

    for (Py_ssize_t i = 0; i < list_size; i++) {
        PyObject* item = PyList_GetItem(kmer_list, i);  // borrowed reference
        if (!PyUnicode_Check(item)) {
            Py_DECREF(result_list);
            PyErr_SetString(PyExc_TypeError, "List items must be strings");
            return NULL;
        }

        const char* current_kmer = PyUnicode_AsUTF8(item);
        int64_t len_current = (int64_t)strlen(current_kmer);

        // Call Myers function
        int64_t dist = myers((uint8_t*)kmer, len_kmer, (uint8_t*)current_kmer, len_current);

        // Convert distance to Python integer and set in list
        PyObject* py_dist = PyLong_FromLongLong(dist);
        PyList_SetItem(result_list, i, py_dist);  // steals reference
    }

    return result_list;
}

static PyMethodDef QuickmersMethods[] = {
    {"hamming_distance_array_32bit", py_hamming_distance_array_32bit, METH_VARARGS, "Compute Hamming distances between query string and list of kmer strings. maximum k = 16."},
    {"hamming_distance_array_64bit", py_hamming_distance_array_64bit, METH_VARARGS, "Compute Hamming distances between query string and list of kmer strings. maximum k = 32."},
    {"hamming_distance_encoded_array_32bit", py_hamming_distance_encoded_array_32bit, METH_VARARGS, "Compute Hamming distances between encoded kmer and list of kmers. maximum k = 16."},
    {"hamming_distance_encoded_32bit", py_hamming_distance_encoded_32bit, METH_VARARGS, "Compute Hamming distance between two encoded kmers. maximum k = 16."},
    {"hamming_distance_32bit", py_hamming_distance_32bit, METH_VARARGS, "Compute Hamming distance between two kmer strings. maximum k = 16."},
    {"hamming_distance_encoded_array_64bit", py_hamming_distance_encoded_array_64bit, METH_VARARGS, "Compute Hamming distances between encoded kmer and list of kmers. maximum k = 32."},
    {"hamming_distance_encoded_64bit", py_hamming_distance_encoded_64bit, METH_VARARGS, "Compute Hamming distance between two encoded kmers. maximum k = 32."},
    {"hamming_distance_64bit", py_hamming_distance_64bit, METH_VARARGS, "Compute Hamming distance between two kmer strings. maximum k = 32."},
    {"levenshtein", py_levenshtein, METH_VARARGS, "Compute Levenshtein edit distance using Myers bit-parallel algorithm. maximum length = 63."},
    {"levenshtein_list", py_levenshtein_list, METH_VARARGS, "Compute Levenshtein edit distances between query string and list of kmer strings using Myers bit-parallel algorithm. maximum length = 63."},
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
    PyObject *module = PyModule_Create(&quickmersmodule);
    if (!module) return NULL;
    import_array();  // initialize NumPy C API
    return module;
}

