#include <Cpython/Call_Cpython.h>
namespace CPython
{
	Call_Cpython::Call_Cpython(const std::string& module, const std::string& function, const std::string& path) : Converter()
	{
		_module = module;
		_function = function;
		_path = path;
	}


	void Call_Cpython::dielectric_params(const double& lamda, std::vector<std::vector<double> >& data)
	{
		int t_end = 0;
#ifndef NDEBUG
		std::cout << "run plot\n";
#endif //  NDEBUG

		PyObject* pName, * pModule, * pFunc;
		PyObject* pArgTuple, * pValue;

		Py_Initialize();

		PyObject* sys = PyImport_ImportModule("sys");
		PyObject* path = PyObject_GetAttrString(sys, "path");
		PyList_Append(path, PyUnicode_FromString(_path.c_str()));

		pName = PyUnicode_FromString(_module.c_str());   //Get the name of the module
		pModule = PyImport_Import(pName);     //Get the module

		Py_DECREF(pName);

		if (pModule != NULL) {
			pFunc = PyObject_GetAttrString(pModule, _function.c_str());   //Get the function by its name
			/* pFunc is a new reference */

			if (pFunc && PyCallable_Check(pFunc)) {

				pArgTuple = PyTuple_New(1);

				//Set the argument tuple to contain the two input tuples
				PyTuple_SetItem(pArgTuple, 0, PyFloat_FromDouble(lamda));

				//Call the python function
				pValue = PyObject_CallObject(pFunc, pArgTuple);
				Py_DECREF(pArgTuple);


				if (pValue != NULL) {
					data = listListTupleTupleToVectorVector_Double(pValue);
					//std::ranges::copy(buf, std::ostream_iterator<double>(std::cout, " "));
					Py_DECREF(pValue);
				}

				//Some error catching
				else {
					Py_DECREF(pFunc);
					Py_DECREF(pModule);
					PyErr_Print();
					std::cout << "Call failed\n";
					return;
				}
			}
			else {
				if (PyErr_Occurred())
					PyErr_Print();
				std::cout << "Cannot find function " << _function << "\n";
				return;
			}
			Py_XDECREF(pFunc);
			Py_DECREF(pModule);
		}
		else {
			PyErr_Print();
			std::cout << "Failed to load " << _module << "\n";
			return;
		}
		Py_Finalize();
	
	}
}
