#pragma once
#include <Utils/Defines.h>
#include <fstream>
#include <sstream>
#include <ctime>
#include <cmath>
#include <Utils/Grid.h>
#include <filesystem>
#include <omp.h>
#include "H5Cpp.h"

namespace fs = std::filesystem; // c++17 и более
/**
 * @class Solver
 * .
 */
template<typename value_t, typename ind_t, typename steps_t>
class Solver : public Grid< ind_t, steps_t>
{
protected:
	/// relative reflection
	steps_t _r_reflection;					
	/// pml borders COORD. First -- left coord ; Second --> right coord
	std::vector<std::pair<steps_t, steps_t> > _pml_crd;
	/// vacuum + skin layers borders COORD. 
	std::vector<std::vector<value_t> > _skin_crd;
	/// pml borders INDEX. First -- left coord ; Second --> right coord
	std::vector<std::pair<ind_t, ind_t> > _pml_ind;
	/// [pml,vacuum + skin, pml]
	std::vector<std::vector<ind_t> > _layers_ind;
	/// pml border WIDTH.
	std::vector<steps_t> _pml_width;
private:

	/// write file flag
	bool _fwrite;
	bool _all_layers;

	/// file path
	std::string _PATH = "../../../Data";

	/// file name
	std::string _FILE_NAME;
	std::vector<value_t> _b;
	std::vector<value_t> _a;
	std::vector<value_t> _sigma;
	std::vector<value_t> _epsilon;

	value_t _lamda;
	value_t _omega_min;
	value_t _omega_max;
	value_t _tau;



	ind_t _layerN_t = 1;
	ind_t _layerN_x = 1;
	ind_t _layerN_y = 1;

	int _n_threads;

	/// size of vector for X-space
	const ind_t _x_vector_size = this->_numberOfStepsX + 1;

	/// size of vector for Y-space
	const ind_t _y_vector_size = this->_numberOfStepsY + 1;

	/// size of vector for Z-space
	const ind_t _z_vector_size = this->_numberOfStepsZ + 1;

	/// size of vector for Time
	const ind_t _t_vector_size = this->_numberOfStepsT + 1;

	/*!
	If folder created - return true, else false
	\param[in] folder path
	*/
	bool demo_exists(const fs::path& p, fs::file_status s = fs::file_status{})
	{
		if (fs::status_known(s) ? fs::exists(s) : fs::exists(p))
		{
			return true;
		}
		else
		{
			return false;
		}
	}

public:
	Solver(std::map<std::string, steps_t>& lborders, std::map<std::string, steps_t>& rborders, std::map<std::string, steps_t>& h,
		steps_t& r_reflection,
		std::vector<std::pair<value_t, value_t> >& pml_crd,
		std::vector<std::vector<value_t> >& skin_crd,
		std::vector<value_t>& a,
		std::vector<value_t>& b,
		std::vector<value_t>& epsilon,
		std::vector<value_t>& sigma,
		Dimension dim,bool fwrite, bool parallel, std::string file_name = "wave", int num_threads = 1);

	/*!
	Calculation of the pml layer index and pml layer width
	*/
	void CalcPMLParam();


	/*!
	Wave for 1D
	\param[in,out] mtrx memory that stores the result in time and x-space
	*/

	void Wave1d(value_t* mtrx, value_t tau, value_t omega, value_t(*func_t0)(steps_t), value_t(*func_u_t0)(steps_t, steps_t, steps_t), value_t(*impulse)(steps_t, steps_t, steps_t, steps_t));
	/*!
	Wave for 2D
	\param[in,out] mtrx memory that stores the result in time and x-space
	*/
	void Wave2d(value_t* mtrx, value_t lamda, value_t omega_min, value_t omega_max, value_t tau,
		value_t(*func_t0)(steps_t, steps_t), value_t(*func_u_t0)(steps_t, steps_t), value_t(*impulse)(steps_t, steps_t, steps_t, steps_t, steps_t, steps_t),
		ind_t layerN_t, ind_t layerN_x, ind_t layerN_y, steps_t centr_gauss, steps_t radius_gauss);

	void WriteParallel(H5::H5File group, std::string dataSetName, value_t* data, ind_t nrows, ind_t ncols);
	void Write_to_extendible_H5(const char* FILENAME, const char* dset_name, value_t* data, ind_t rows, ind_t cols);
	void Write_param(std::string filenane, bool hdf5_flag);
	ind_t idx(ind_t p)
	{
		if (_fwrite && !_all_layers)
		{
			return p % 3;
		}
		else
		{
			return p;
		}
	};
	ind_t idx_phi(ind_t p)
	{
		if (_fwrite && !_all_layers)
		{
			return p % 2;
		}
		else
		{
			return p;
		}
	};
};
