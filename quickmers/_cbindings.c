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

static PyObject* py_levenshtein_list_avx2(PyObject* self, PyObject* args) {
    const char* query;
    PyObject* kmer_list;

    if (!PyArg_ParseTuple(args, "sO!", &query, &PyList_Type, &kmer_list))
        return NULL;

    Py_ssize_t n = PyList_Size(kmer_list);
    if (n <= 0) return NULL;

    // Convert Python list to array of pointers
    const uint8_t **kmers = (const uint8_t**)malloc(sizeof(uint8_t*)*n);
    for (Py_ssize_t i = 0; i < n; i++) {
        PyObject* item = PyList_GetItem(kmer_list, i);
        if (!PyUnicode_Check(item)) {
            free(kmers);
            PyErr_SetString(PyExc_TypeError, "List items must be strings");
            return NULL;
        }
        kmers[i] = (const uint8_t*)PyUnicode_AsUTF8(item);
    }

    int64_t kmer_len = strlen((const char*)kmers[0]);
    int64_t *out = (int64_t*)malloc(sizeof(int64_t)*n);

    myers_batch_avx2((const uint8_t*)query, strlen(query), kmers, n, kmer_len, out);

    PyObject* result_list = PyList_New(n);
    for (Py_ssize_t i = 0; i < n; i++)
        PyList_SetItem(result_list, i, PyLong_FromLongLong(out[i]));

    free(kmers);
    free(out);
    return result_list;
}

static PyObject* py_levenshtein_list_avx2_numpy(PyObject* self, PyObject* args) {
    const char* query;
    PyObject* kmer_list;

    if (!PyArg_ParseTuple(args, "sO!", &query, &PyList_Type, &kmer_list))
        return NULL;

    Py_ssize_t n = PyList_Size(kmer_list);
    if (n <= 0) {
        // return an empty NumPy array
        npy_intp dims[1] = {0};
        return PyArray_SimpleNew(1, dims, NPY_INT64);
    }

    // Convert Python list to array of pointers
    const uint8_t **kmers = (const uint8_t**)malloc(sizeof(uint8_t*) * n);
    if (!kmers) return PyErr_NoMemory();

    for (Py_ssize_t i = 0; i < n; i++) {
        PyObject* item = PyList_GetItem(kmer_list, i);
        if (!PyUnicode_Check(item)) {
            free(kmers);
            PyErr_SetString(PyExc_TypeError, "List items must be strings");
            return NULL;
        }
        kmers[i] = (const uint8_t*)PyUnicode_AsUTF8(item);
    }

    int64_t kmer_len = strlen((const char*)kmers[0]);
    int64_t* out = (int64_t*)malloc(sizeof(int64_t) * n);
    if (!out) {
        free(kmers);
        return PyErr_NoMemory();
    }

    myers_batch_avx2((const uint8_t*)query, strlen(query), kmers, n, kmer_len, out);

    // Create NumPy array
    npy_intp dims[1] = { n };
    PyObject* np_array = PyArray_SimpleNew(1, dims, NPY_INT64);
    if (!np_array) {
        free(kmers);
        free(out);
        return NULL;
    }

    // Copy results into NumPy array buffer
    int64_t* arr_data = (int64_t*)PyArray_DATA((PyArrayObject*)np_array);
    memcpy(arr_data, out, n * sizeof(int64_t));

    free(kmers);
    free(out);

    return np_array;
}

static PyObject* py_levenshtein_list_with_min_dist(PyObject* self, PyObject* args) {
    const char* kmer;
    PyObject* kmer_list;
    long min_distance;  // new arg

    // Parse: string, list, integer
    if (!PyArg_ParseTuple(args, "sO!l", &kmer, &PyList_Type, &kmer_list, &min_distance))
        return NULL;

    Py_ssize_t list_size = PyList_Size(kmer_list);
    if (list_size <= 0) {
        // return (True, empty list) instead of error
        PyObject* result_list = PyList_New(0);
        if (!result_list) return NULL; // malloc failure
        PyObject* result_tuple = PyTuple_New(2);
        if (!result_tuple) {
            Py_DECREF(result_list);
            return NULL;
        }
        PyTuple_SetItem(result_tuple, 0, PyBool_FromLong(1)); // True
        PyTuple_SetItem(result_tuple, 1, result_list);
        return result_tuple;
    }

    int64_t len_kmer = (int64_t)strlen(kmer);

    // Create Python list to store results
    PyObject* result_list = PyList_New(list_size);
    if (!result_list) return NULL;

    int early_exit = 0;

    for (Py_ssize_t i = 0; i < list_size; i++) {
        PyObject* item = PyList_GetItem(kmer_list, i);
        if (!PyUnicode_Check(item)) {
            Py_DECREF(result_list);
            PyErr_SetString(PyExc_TypeError, "List items must be strings");
            return NULL;
        }

        const char* current_kmer = PyUnicode_AsUTF8(item);
        int64_t len_current = (int64_t)strlen(current_kmer);

        int64_t dist = myers((uint8_t*)kmer, len_kmer, (uint8_t*)current_kmer, len_current);

        // Store distance
        PyObject* py_dist = PyLong_FromLongLong(dist);
        PyList_SetItem(result_list, i, py_dist);

        if (dist < min_distance) {
            early_exit = 1;
            // Truncate list (set remaining items to None)
            for (Py_ssize_t j = i+1; j < list_size; j++) {
                Py_INCREF(Py_None);
                PyList_SetItem(result_list, j, Py_None);
            }
            break;
        }
    }

    // Return tuple: (bool, list)
    PyObject* result_tuple = PyTuple_New(2);
    PyTuple_SetItem(result_tuple, 0, PyBool_FromLong(!early_exit));
    PyTuple_SetItem(result_tuple, 1, result_list);

    return result_tuple;
}

static PyObject* py_levenshtein_list_avx2_with_min_dist(PyObject* self, PyObject* args) {
    const char* kmer;
    PyObject* kmer_list;
    long min_distance;

    if (!PyArg_ParseTuple(args, "sO!l", &kmer, &PyList_Type, &kmer_list, &min_distance))
        return NULL;

    Py_ssize_t list_size = PyList_Size(kmer_list);
    if (list_size <= 0) {
        // return (True, empty list) instead of error
        PyObject* result_list = PyList_New(0);
        if (!result_list) return NULL; // malloc failure
        PyObject* result_tuple = PyTuple_New(2);
        if (!result_tuple) {
            Py_DECREF(result_list);
            return NULL;
        }
        PyTuple_SetItem(result_tuple, 0, PyBool_FromLong(1)); // True
        PyTuple_SetItem(result_tuple, 1, result_list);
        return result_tuple;
    }

    const uint8_t **kmers = (const uint8_t**)malloc(sizeof(uint8_t*)*list_size);
    for (Py_ssize_t i = 0; i < list_size; i++) {
        PyObject* item = PyList_GetItem(kmer_list, i);
        if (!PyUnicode_Check(item)) {
            free(kmers);
            PyErr_SetString(PyExc_TypeError, "List items must be strings");
            return NULL;
        }
        kmers[i] = (const uint8_t*)PyUnicode_AsUTF8(item);
    }

    int64_t kmer_len = strlen((const char*)kmers[0]);
    int64_t *out = (int64_t*)malloc(sizeof(int64_t)*list_size);

    myers_batch_avx2((const uint8_t*)kmer, strlen(kmer), kmers, list_size, kmer_len, out);

    // Build result list
    PyObject* result_list = PyList_New(list_size);
    int early_exit = 0;
    for (Py_ssize_t i = 0; i < list_size; i++) {
        PyObject* py_dist = PyLong_FromLongLong(out[i]);
        PyList_SetItem(result_list, i, py_dist);

        if (out[i] < min_distance) {
            early_exit = 1;
            // Truncate rest
            for (Py_ssize_t j = i+1; j < list_size; j++) {
                Py_INCREF(Py_None);
                PyList_SetItem(result_list, j, Py_None);
            }
            break;
        }
    }

    free(kmers);
    free(out);

    PyObject* result_tuple = PyTuple_New(2);
    PyTuple_SetItem(result_tuple, 0, PyBool_FromLong(!early_exit));
    PyTuple_SetItem(result_tuple, 1, result_list);

    return result_tuple;
}

static PyObject* py_levenshtein_list_avx2_with_min_dist_numpy(PyObject* self, PyObject* args) {
    const char* kmer;
    PyObject* kmer_list;
    long min_distance;

    if (!PyArg_ParseTuple(args, "sO!l", &kmer, &PyList_Type, &kmer_list, &min_distance))
        return NULL;

    Py_ssize_t list_size = PyList_Size(kmer_list);
    if (list_size <= 0) {
        // return (True, empty np.array([]))
        npy_intp dims[1] = {0};
        PyObject* empty_array = PyArray_SimpleNew(1, dims, NPY_INT64);
        if (!empty_array) return NULL;

        PyObject* result_tuple = PyTuple_New(2);
        if (!result_tuple) {
            Py_DECREF(empty_array);
            return NULL;
        }
        PyTuple_SetItem(result_tuple, 0, PyBool_FromLong(1)); // True
        PyTuple_SetItem(result_tuple, 1, empty_array);
        return result_tuple;
    }

    const uint8_t **kmers = (const uint8_t**)malloc(sizeof(uint8_t*) * list_size);
    if (!kmers) return PyErr_NoMemory();

    for (Py_ssize_t i = 0; i < list_size; i++) {
        PyObject* item = PyList_GetItem(kmer_list, i);
        if (!PyUnicode_Check(item)) {
            free(kmers);
            PyErr_SetString(PyExc_TypeError, "List items must be strings");
            return NULL;
        }
        kmers[i] = (const uint8_t*)PyUnicode_AsUTF8(item);
    }

    int64_t kmer_len = strlen((const char*)kmers[0]);
    int64_t* out = (int64_t*)malloc(sizeof(int64_t) * list_size);
    if (!out) {
        free(kmers);
        return PyErr_NoMemory();
    }

    myers_batch_avx2((const uint8_t*)kmer, strlen(kmer), kmers, list_size, kmer_len, out);

    // Create NumPy array
    npy_intp dims[1] = { list_size };
    PyObject* np_array = PyArray_SimpleNew(1, dims, NPY_INT64);
    if (!np_array) {
        free(kmers);
        free(out);
        return NULL;
    }
    int64_t* arr_data = (int64_t*)PyArray_DATA((PyArrayObject*)np_array);

    int early_exit = 0;
    for (Py_ssize_t i = 0; i < list_size; i++) {
        arr_data[i] = out[i];
        if (out[i] < min_distance) {
            early_exit = 1;
            // Fill remaining with -1 sentinel
            for (Py_ssize_t j = i+1; j < list_size; j++) {
                arr_data[j] = -1;
            }
            break;
        }
    }

    free(kmers);
    free(out);

    PyObject* result_tuple = PyTuple_New(2);
    PyTuple_SetItem(result_tuple, 0, PyBool_FromLong(!early_exit));
    PyTuple_SetItem(result_tuple, 1, np_array);

    return result_tuple;
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
    {"levenshtein_list_avx2", py_levenshtein_list_avx2, METH_VARARGS, "Compute Levenshtein edit distances between query string and list of kmer strings using Myers bit-parallel algorithm. maximum length = 63. Uses avx2 operations and is faster than levenshtein_list if cpu supports the avx2 commands."},
    {"levenshtein_list_avx2_numpy", py_levenshtein_list_avx2_numpy, METH_VARARGS, "Compute Levenshtein edit distances between query string and list of kmer strings using Myers bit-parallel algorithm. maximum length = 63. Uses avx2 operations and is faster than levenshtein_list if cpu supports the avx2 commands."},
    {"levenshtein_list_with_min_dist", py_levenshtein_list_with_min_dist, METH_VARARGS, "Compute Levenshtein edit distances between query string and list of kmer strings unless one of the pairs violate min distance constraint. maximum length = 63."},
    {"levenshtein_list_avx2_with_min_dist", py_levenshtein_list_avx2_with_min_dist, METH_VARARGS, "Compute Levenshtein edit distances between query string and list of kmer strings unless one of the pairs violate min distance constraint. maximum length = 63. Uses avx2 operations and is faster than levenshtein_list if cpu supports the avx2 commands."},
    {"levenshtein_list_avx2_with_min_dist_numpy", py_levenshtein_list_avx2_with_min_dist_numpy,  METH_VARARGS, "Compute Levenshtein edit distances between query string and list of kmer strings unless one of the pairs violate min distance constraint. maximum length = 63. Uses avx2 operations and is faster than levenshtein_list if cpu supports the avx2 commands."},
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
    import_array();
    return module;
}

