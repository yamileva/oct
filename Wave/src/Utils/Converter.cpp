#include <Utils/Converter.h>

namespace Convert
{
	Converter::Converter() {

	}
	// =====
	// LISTS
	// =====
	PyObject* Converter::vectorToList_Double(const std::vector<double>& data) {
		PyObject* listObj = PyList_New(data.size());
		if (!listObj) throw std::logic_error("Unable to allocate memory for Python list");
		for (unsigned int i = 0; i < data.size(); i++) {
			PyObject* num = PyFloat_FromDouble(data[i]);
			if (!num) {
				Py_DECREF(listObj);
				throw std::logic_error("Unable to allocate memory for Python list");
			}
			PyList_SET_ITEM(listObj, i, num);
		}
		return listObj;
	}

	// ======
	// TUPLES
	// ======

	PyObject* Converter::vectorToTuple_Double(const std::vector<double>& data) {
		PyObject* tuple = PyTuple_New(data.size());
		if (!tuple) throw std::logic_error("Unable to allocate memory for Python tuple");
		for (unsigned int i = 0; i < data.size(); i++) {
			PyObject* num = PyFloat_FromDouble(data[i]);
			if (!num) {
				Py_DECREF(tuple);
				throw std::logic_error("Unable to allocate memory for Python tuple");
			}
			PyTuple_SET_ITEM(tuple, i, num);
		}

		return tuple;
	}

	PyObject* Converter::vectorVectorToTuple_Double(const std::vector< std::vector< double > >& data) {
		PyObject* tuple = PyTuple_New(data.size());
		if (!tuple) throw std::logic_error("Unable to allocate memory for Python tuple");
		for (unsigned int i = 0; i < data.size(); i++) {
			PyObject* subTuple = NULL;
			try {
				subTuple = vectorToTuple_Double(data[i]);
			}
			catch (std::logic_error& e) {
				throw e;
			}
			if (!subTuple) {
				Py_DECREF(tuple);
				throw std::logic_error("Unable to allocate memory for Python tuple of tuples");
			}
			PyTuple_SET_ITEM(tuple, i, subTuple);
		}

		return tuple;
	}


	PyObject* Converter::vectorVectorVectorToTuple_Double(const std::vector<std::vector< std::vector< double > > >& data) {
		PyObject* tuple = PyTuple_New(data.size());
		if (!tuple) throw std::logic_error("Unable to allocate memory for Python tuple");
		for (unsigned int i = 0; i < data.size(); i++) {
			PyObject* subTuple = NULL;
			try {
				subTuple = vectorVectorToTuple_Double(data[i]);
			}
			catch (std::logic_error& e) {
				throw e;
			}
			if (!subTuple) {
				Py_DECREF(tuple);
				throw std::logic_error("Unable to allocate memory for Python tuple of tuples");
			}
			PyTuple_SET_ITEM(tuple, i, subTuple);
		}

		return tuple;
	}
	// PyObject -> Vector
	std::vector<double> Converter::listTupleToVector_Double(PyObject* incoming) {
		std::vector < double > data;
		if (PyTuple_Check(incoming)) {
			for (Py_ssize_t i = 0; i < PyTuple_Size(incoming); i++) {
				PyObject* value = PyTuple_GetItem(incoming, i);
				data.push_back(PyFloat_AsDouble(value));
			}
		}
		else {
			if (PyList_Check(incoming)) {
				for (Py_ssize_t i = 0; i < PyList_Size(incoming); i++) {
					PyObject* value = PyList_GetItem(incoming, i);
					data.push_back(PyFloat_AsDouble(value));
				}
			}
			else {
				throw std::logic_error("Passed PyObject pointer was not a list or tuple!");
			}
		}
		return data;
	}

	std::vector<std::vector<double> > Converter::listListTupleTupleToVectorVector_Double(PyObject* incoming) {
		std::vector<std::vector<double> > data;
		if (PyTuple_Check(incoming)) {
			for (Py_ssize_t i = 0; i < PyTuple_Size(incoming); i++) {
				PyObject* value = PyTuple_GetItem(incoming, i);
				data.push_back(listTupleToVector_Double(value));
			}
		}
		else {
			if (PyList_Check(incoming)) {
				for (Py_ssize_t i = 0; i < PyList_Size(incoming); i++) {
					PyObject* value = PyList_GetItem(incoming, i);
					data.push_back(listTupleToVector_Double(value));
				}
			}
			else {
				throw std::logic_error("Passed PyObject pointer was not a list or tuple!");
			}
		}
		return data;
	}

	// PyObject -> Vector
	std::vector<int> Converter::listTupleToVector_Int(PyObject* incoming) {
		std::vector<int> data;
		if (PyTuple_Check(incoming)) {
			for (Py_ssize_t i = 0; i < PyTuple_Size(incoming); i++) {
				PyObject* value = PyTuple_GetItem(incoming, i);
				data.push_back(PyFloat_AsDouble(value));
			}
		}
		else {
			if (PyList_Check(incoming)) {
				for (Py_ssize_t i = 0; i < PyList_Size(incoming); i++) {
					PyObject* value = PyList_GetItem(incoming, i);
					data.push_back(PyFloat_AsDouble(value));
				}
			}
			else {
				throw std::logic_error("Passed PyObject pointer was not a list or tuple!");
			}
		}
		return data;
	}
}