#pragma once
#include <Python.h>
#include <vector>
#include <stdexcept>
namespace Convert {
	/**
	 * @class Converter
	 * Converts vectors to Python objects
	 */
	class Converter
	{
	public:
		Converter();
	protected:
		PyObject* vectorToList_Double(const std::vector<double>& data);

		PyObject* vectorToTuple_Double(const std::vector<double>& data);

		PyObject* vectorVectorToTuple_Double(const std::vector< std::vector< double > >& data);
		PyObject* vectorVectorVectorToTuple_Double(const std::vector<std::vector< std::vector< double > > >& data);

		std::vector<double> listTupleToVector_Double(PyObject* incoming);
		std::vector < std::vector<double> > listListTupleTupleToVectorVector_Double(PyObject * incoming);

		std::vector<int> listTupleToVector_Int(PyObject* incoming);
	};


}

