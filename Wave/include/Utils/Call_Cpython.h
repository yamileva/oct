#pragma once
#include <Utils/Converter.h>
#include <numpy/arrayobject.h>
#include <Utils/Grid.h>
#include <Utils/Defines.h>
#include <ctime>
#include <sstream>
#include <iostream>
#include <string>
#include <algorithm>

namespace CPython
{
	/**
	 * @class CGraph
	 * Sends vector to the numpy module
	 */
	class Call_Cpython : public Convert::Converter
	{
	public:
		/*!
		Sends vector to the numpy module
		\param[in] module name 
		\param[in] function name 
		\param[in] module path
		*/
		Call_Cpython(const std::string& module, const std::string& function, const std::string& path);

		/*
		int plot(Dimension& dim, double* mtrx, std::vector<double>& grid_settings, std::vector<std::vector<double>>& pml_setting, double* xi);
		int plot(Dimension& dim, std::vector<std::vector<double>>& mtrx, std::vector<double>& grid_settings, std::vector<std::vector<double>>& pml_setting, std::vector<double>& xi);*/

		void dielectric_params(const double& lamda, std::vector<std::vector<double> >& data);
		
	private:

		std::string _module;
		std::string _function;
		std::string _path;
	};
}

