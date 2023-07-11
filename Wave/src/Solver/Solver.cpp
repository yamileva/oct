#include <Solver/Solver.h>
template<typename value_t, typename ind_t, typename steps_t>
Solver<value_t, ind_t, steps_t>::Solver(std::map<std::string, steps_t>& lborders, std::map<std::string, steps_t>& rborders,
	std::map<std::string, steps_t>& h,
	steps_t& r_reflection,
	std::vector<std::pair<value_t, value_t> >& pml_crd,
	std::vector<std::vector<value_t> >& skin_crd,
	std::vector<value_t>& a,
	std::vector<value_t>& b,
	std::vector<value_t>& epsilon,
	std::vector<value_t>& sigma,
	Dimension dim,
	bool fwrite,bool parallel, std::string file_name, int num_threads) : Grid< ind_t, steps_t>(lborders, rborders, h, dim)
{
	_r_reflection = r_reflection;
	_pml_crd = pml_crd;
	_skin_crd = skin_crd;
	_a = a;
	_b = b;
	_epsilon = epsilon;
	_sigma = sigma;
	_fwrite = fwrite;
	_all_layers = parallel;
	_FILE_NAME = file_name;
	_n_threads = num_threads;
	CalcPMLParam();
	if (!demo_exists(_PATH)) {
		fs::create_directories(_PATH);
	}
}
//----------------------------------------------------------------------------------------------------------------------------------
template<typename value_t, typename ind_t, typename steps_t>
void Solver<value_t, ind_t, steps_t>::CalcPMLParam() {
	std::vector < std::string > tmp = { "x","y","z" };
	for (ind_t i = 0; i < _pml_crd.size(); i++)
	{
		_pml_ind.push_back(std::make_pair((_pml_crd[i].first - this->_lborders[tmp[i]]) / this->_h[tmp[i]],
			(_pml_crd[i].second - this->_lborders[tmp[i]]) / this->_h[tmp[i]]));
		// PML SECOND
		std::vector<ind_t> buf(static_cast<ind_t>(_skin_crd[i].size() + 3));
		buf[0] = _pml_ind[i].first + 1;
		buf[buf.size() - 2] = _pml_ind[i].second + 1;
		buf[buf.size() - 1] = (this->_rborders[tmp[i]] - this->_lborders[tmp[i]]) / this->_h[tmp[i]];

		for (ind_t j = 1; j < buf.size() - 2; j++)
		{
			buf[j] = (_skin_crd[i][j - 1] - this->_lborders[tmp[i]]) / this->_h[tmp[i]];
		}

		_layers_ind.push_back(buf);
		buf.clear();
		_pml_width.push_back(_pml_crd[i].first - this->_lborders[tmp[i]]);
	}
};
//----------------------------------------------------------------------------------------------------------------------------------
template<typename value_t, typename ind_t, typename steps_t>
void Solver<value_t, ind_t, steps_t>::Wave1d(value_t* mtrx, value_t tau, value_t omega,
	value_t(*func_t0)(steps_t), value_t(*func_u_t0)(steps_t, steps_t, steps_t), value_t(*impulse)(steps_t, steps_t, steps_t, steps_t)) {

	std::string filename = _PATH + "/" + _FILE_NAME + ".h5";
	hid_t file = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
	std::cout << _FILE_NAME << "\n";
	//запись параметров сетки
	std::string output = "Without write file";
	if (_fwrite)
	{
		output = "With write file";
		Write_param(filename, true);
	}

	ind_t xi_left_size = (_pml_ind.front().first) * 2 + 1;
	ind_t xi_right_size = ((_x_vector_size - 1) - _pml_ind.front().second) * 2 + 1;

	ind_t phi_left_size = (_pml_ind.front().first);
	ind_t phi_right_size = ((_x_vector_size - 1) - _pml_ind.front().second);

	value_t* _xi_left = (value_t*)calloc(xi_left_size, sizeof(value_t));
	value_t* _xi_right = (value_t*)calloc(xi_right_size, sizeof(value_t));

	ind_t phi_time;
	if (_all_layers) phi_time = this->_numberOfStepsT + 1;
	else phi_time = 2;
	value_t* _phi_left = (value_t*)calloc(phi_time * phi_left_size, sizeof(value_t));
	value_t* _phi_right = (value_t*)calloc(phi_time * phi_right_size, sizeof(value_t));

	auto base = [&](value_t a, value_t b) -> value_t {
		return ((2. * this->_h_time * this->_h_time) / (2. * b + a * this->_h_time));
	};
	auto alpha = [&](value_t a, value_t b) -> value_t {
		return base(a, b) * (1. / (this->_h_x * this->_h_x));
	};
	auto gamma = [&](value_t a, value_t b) -> value_t {
		return  ((4. * b) / (2. * b + a * this->_h_time)) - 2. * alpha(a, b);
	};
	auto betta = [&](value_t a, value_t b) -> value_t {
		return (a * this->_h_time - 2. * b) / (2. * b + a * this->_h_time);
	};
	value_t _A;
	ind_t i_shift;

	double t_total = omp_get_wtime();
	for (ind_t i = 0; i < _x_vector_size; i++)
	{
		mtrx[i] = func_t0(i * this->_h_x + this->_lborders_x);
		//mtrx[1 * _x_vector_size + i] = mtrx[i] + func_u_t0(i * this->_h_x + this->_lborders_x,0,omega);
	}
	
	ind_t l_prev = static_cast<ind_t>(1.);
	ind_t l_ind = static_cast<ind_t>(0.);
	for(const auto &l : _layers_ind.front())
	{
		for (ind_t i = l_prev; i < l; i++)
		{
			value_t x = i * this->_h_x + this->_lborders_x;
			value_t t = 0.;
			/*
			if (l_ind == 0) {
				_A = _xi_left[2 * i] * _b[l_ind];
				mtrx[1 * _x_vector_size + i] = mtrx[i] + func_u_t0(x,t,omega) * this->_h_time
					+ (this->_h_time * this->_h_time) / 2
					* (-(_A / _b[l_ind]) * func_u_t0(x, t, omega));
			}
			else if (l_ind == _layers_ind.front().size() - 1)
			{
				_A = _xi_right[2 * i_shift] * _b[l_ind];
				mtrx[1 * _x_vector_size + i] = mtrx[i] + func_u_t0(x, t, omega) * this->_h_time
					+ (this->_h_time * this->_h_time) / 2
					* (-(_A / _b[l_ind]) * func_u_t0(x, t, omega));
			}
			else*/
			if(i == _pml_ind.front().first)
			{
				mtrx[1 * _x_vector_size + i] = mtrx[i] + func_u_t0(x, t, omega) * this->_h_time
					+ ((this->_h_time * this->_h_time) / 2)
					* (-(_a[l_ind - 1] / _b[l_ind]) * func_u_t0(x, t, omega));
	 		}
			else
			{
				mtrx[1 * _x_vector_size + i] = mtrx[i];
			}
		}
		l_prev = l;
		l_ind++;
	}
	

	// XI FOR PML
	for (ind_t i = 0; i < xi_left_size; i++)
	{
		_xi_left[i] = (1. / (sqrt(_b.front()) * _pml_width.front()))
			* std::log10(1. / _r_reflection)
			* (
				abs((i * this->_h_x / 2. + this->_lborders_x) - _pml_crd.front().first) / _pml_width.front()
				- sin(2. * M_PI * abs((i * this->_h_x / 2. + this->_lborders_x)
					- _pml_crd.front().first) / _pml_width.front()) / (2. * M_PI)
				);
	}
	for (ind_t i = 0; i < xi_right_size; i++)
	{
		_xi_right[i] = (1 / (sqrt(_b.back()) * _pml_width.front()))
			* std::log10(1. / _r_reflection)
			* (
				abs(((i + 2. * _pml_ind.front().second) * this->_h_x / 2. + this->_lborders_x) - _pml_crd.front().second) / _pml_width.front()
				- sin(2. * M_PI * abs(((i + 2. * _pml_ind.front().second) * this->_h_x / 2. + this->_lborders_x)
					- _pml_crd.front().second) / _pml_width.front()) / (2. * M_PI)
				);
	}
	//
	std::cout << "Start sheme\n";
	for (ind_t p = 1; p < _t_vector_size - 1; p++)
	{
		//запись в файл трех слоев по времени
		if (_fwrite && !_all_layers && p % 10 == 0 ) {
			Write_to_extendible_H5(filename.c_str(), "wave", mtrx, 1, _x_vector_size);
		}

		double s_iter = omp_get_wtime();
		mtrx[idx(p + 1) * _x_vector_size] = static_cast<value_t>(0.);
		

		ind_t l_prev = static_cast<ind_t>(1.);
		ind_t l_ind = static_cast<ind_t>(0.);
		for(const auto &l : _layers_ind.front())
		{
			for (ind_t i = l_prev; i < l; i++)
			{	
				// испульс некоторое время tau // p * this->_h_time <= tau && 
				if (p * this->_h_time <= tau && i == _pml_ind.front().first) {
					value_t x = i * this->_h_x + this->_lborders_x;
					value_t t = p * this->_h_time;
					mtrx[idx(p + 1) * _x_vector_size + i] = impulse(x, t, omega, tau);
				}
				/*else if (p * this->_h_time <= 2*tau && i * this->_h_x + this->_lborders_x == _pml_crd.front().first) {
					mtrx[idx(p + 1) * _x_vector_size + i] = static_cast<value_t>(0.);
				}*/
				else
				{
					//pml слева
					if (l_ind == 0)
					{
						_phi_left[idx_phi(p) * phi_left_size + i] = _phi_left[idx_phi(p - 1) * phi_left_size + i] * ((2. - _xi_left[2 * i + 1] * this->_h_time) / (2. + _xi_left[2 * i + 1] * this->_h_time))
							- ((_xi_left[2 * i + 1] * this->_h_time * (1. / _b[l_ind])) / (2. + _xi_left[2 * i + 1] * this->_h_time))
							* ((mtrx[idx(p) * _x_vector_size + (i + 1)] - mtrx[idx(p) * _x_vector_size + i]) / this->_h_x
								+ (mtrx[idx(p - 1) * _x_vector_size + (i + 1)] - mtrx[idx(p - 1) * _x_vector_size + i]) / this->_h_x
								);

						_A = _xi_left[2 * i] * _b[l_ind];

						mtrx[idx(p + 1) * _x_vector_size + i] = alpha(_A, _b[l_ind]) * mtrx[idx(p) * _x_vector_size + (i - 1)]
							+ gamma(_A, _b[l_ind]) * mtrx[idx(p) * _x_vector_size + i]
							+ alpha(_A, _b[l_ind]) * mtrx[idx(p) * _x_vector_size + (i + 1)]
							+ betta(_A, _b[l_ind]) * mtrx[idx(p - 1) * _x_vector_size + i]
							+ base(_A, _b[l_ind]) * _b[l_ind]
							* ((_phi_left[idx_phi(p) * phi_left_size + i] - _phi_left[idx_phi(p) * phi_left_size + (i - 1)]) / this->_h_x);
						
					}
					// pml справа
					else if (l_ind == _layers_ind.front().size() - 1)
					{
						i_shift = i - _pml_ind.front().second;

						_phi_right[idx_phi(p) * phi_right_size + (i_shift + 1)] = _phi_right[idx_phi(p - 1) * phi_right_size + (i_shift + 1)] * ((2. - _xi_right[2 * i_shift + 1] * this->_h_time) / (2. + _xi_right[2 * i_shift + 1] * this->_h_time))
							- ((_xi_right[2 * i_shift + 1] * this->_h_time * (1. / _b[l_ind])) / (2. + _xi_right[2 * i_shift + 1] * this->_h_time))
							* ((mtrx[idx(p) * _x_vector_size + (i + 1)] - mtrx[idx(p) * _x_vector_size + i]) / this->_h_x
								+ (mtrx[idx(p - 1) * _x_vector_size + (i + 1)] - mtrx[idx(p - 1) * _x_vector_size + i]) / this->_h_x
								);

						_A = _xi_right[2 * i_shift] * _b[l_ind];

						mtrx[idx(p + 1) * _x_vector_size + i] = alpha(_A, _b[l_ind]) * mtrx[idx(p) * _x_vector_size + (i - 1)]
							+ gamma(_A, _b[l_ind]) * mtrx[idx(p) * _x_vector_size + i]
							+ alpha(_A, _b[l_ind]) * mtrx[idx(p) * _x_vector_size + (i + 1)]
							+ betta(_A, _b[l_ind]) * mtrx[idx(p - 1) * _x_vector_size + i]
							+ base(_A, _b[l_ind]) * _b[l_ind]
							* ((_phi_right[idx_phi(p) * phi_right_size + (i_shift + 1)] - _phi_right[idx_phi(p) * phi_right_size + i_shift]) / this->_h_x);
					}
					// wave
					else
					{
						mtrx[idx(p + 1) * _x_vector_size + i] = alpha(_a[l_ind - 1], _b[l_ind]) * mtrx[idx(p) * _x_vector_size + (i - 1)]
							+ gamma(_a[l_ind - 1], _b[l_ind]) * mtrx[idx(p) * _x_vector_size + i]
							+ alpha(_a[l_ind - 1], _b[l_ind]) * mtrx[idx(p) * _x_vector_size + (i + 1)]
							+ betta(_a[l_ind - 1], _b[l_ind]) * mtrx[idx(p - 1) * _x_vector_size + i];
					}
				}
			}
			if (l_ind != 0 && l_ind != _layers_ind.front().size() - 2 && l_ind != _layers_ind.front().size() - 1)
			{
				l_prev = l + 1;
			}
			else
			{
				l_prev = l;
			}
			l_ind++;
		}

		for (ind_t l = 0; l < _layers_ind.front().size(); l++)
		{
			ind_t i = _layers_ind.front()[l];
			if (l !=0 && l != _layers_ind.front().size() - 2 && l != _layers_ind.front().size() - 1)
			{
				mtrx[idx(p + 1) * _x_vector_size + i] = (1. / 6) * (-mtrx[idx(p + 1) * _x_vector_size + (i - 2)]
					+ 4. * mtrx[idx(p + 1) * _x_vector_size + (i - 1)]
					+ 4. * mtrx[idx(p + 1) * _x_vector_size + (i + 1)]
					- mtrx[idx(p + 1) * _x_vector_size + (i + 2)]);
			}
		}

		mtrx[idx(p + 1) * _x_vector_size + (_x_vector_size - 1)] = static_cast<value_t>(0.);

		if (p == 1)
		{
			double t_iter = omp_get_wtime() - s_iter;
			std::cout << "Time for one iteration: " << t_iter <<" seconds.\n";
			std::cout << "Number of iterations:" << this->_numberOfStepsT << "\n";
			std::cout << "Time for all iteration: " << t_iter * this->_numberOfStepsT <<" seconds.\n";
		}
		if (p % (this->_numberOfStepsT / 100) == 0) std::cout << "End calculation for t = " << (p * 100) / (this->_numberOfStepsT - 1) << "%\n";
	}
	//запись в файл остальных слоев
	if (_fwrite && !_all_layers) {
		Write_to_extendible_H5(filename.c_str(), "wave", mtrx, idx(_t_vector_size), _x_vector_size);
	}

	t_total = omp_get_wtime() - t_total;
	std::cout << "Time Sheme : " << t_total <<" seconds.\n";

	// запись всех слоев по времени
	if (_all_layers)
	{
		t_total = omp_get_wtime();
		WriteParallel(file, "wave", mtrx, _t_vector_size, _x_vector_size);
		t_total = omp_get_wtime() - t_total;
		std::cout << "Write file : " << t_total << " seconds.\n";
	}
	H5Fclose(file);
	free(_phi_right);
	free(_xi_right);
	free(_xi_left);
	free(_phi_left);
}
//----------------------------------------------------------------------------------------------------------------------------------
template<typename value_t, typename ind_t, typename steps_t>
void Solver<value_t, ind_t, steps_t>::Wave2d(value_t* mtrx, value_t lamda , value_t omega_min, value_t omega_max, value_t tau,
	value_t(*func_t0)(steps_t, steps_t), value_t(*func_u_t0)(steps_t, steps_t), value_t(*impulse)(steps_t, steps_t, steps_t, steps_t, steps_t, steps_t),
	ind_t layerN_t, ind_t layerN_x, ind_t layerN_y, steps_t centr_gauss, steps_t radius_gauss) {

	_layerN_t = layerN_t;
	_layerN_x = layerN_x;
	_layerN_y = layerN_y;

	_lamda = lamda;
	_omega_min = omega_min;
	_omega_max = omega_max;
	_tau = tau;
	
	std::string filename = _PATH + "/" + _FILE_NAME + ".h5";
	hid_t file = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
	std::cout << _FILE_NAME << "\n";

	//запись параметров сетки
	std::string output = "Without write file";
	if (_fwrite)
	{
		output = "With write file";
		Write_param(filename, true);
	}

	//------------------------------------------------------------------------------
	ind_t xi1_size = (_x_vector_size - 1) * 2 + 1;

	ind_t xi2_size = (_y_vector_size - 1) * 2 + 1;
	//------------------------------------------------------------------------------
	ind_t phi1_size = ((_x_vector_size - 1) * (_y_vector_size - 1));
	ind_t phi2_size = ((_x_vector_size - 1) * (_y_vector_size - 1));
	//------------------------------------------------------------------------------
	value_t* _xi1 = (value_t*)calloc(xi1_size, sizeof(value_t));
	value_t* _xi2 = (value_t*)calloc(_layers_ind.front().size() * xi2_size, sizeof(value_t));
	//------------------------------------------------------------------------------

	ind_t phi_time;
	if (_all_layers) phi_time = this->_numberOfStepsT + 1;
	else phi_time = 2;
	value_t* _phi1 = (value_t*)calloc(phi_time * phi1_size, sizeof(value_t));
	value_t* _phi2 = (value_t*)calloc(phi_time * phi2_size, sizeof(value_t));
	//------------------------------------------------------------------------------

	auto base = [&](value_t a, value_t b) -> value_t {
		return ((2. * this->_h_time * this->_h_time) / (2. * b + a * this->_h_time));
	};
	auto alpha1 = [&](value_t a, value_t b) -> value_t {
		return  base(a, b) * (1. / (this->_h_x * this->_h_x));
	};
	auto alpha2 = [&](value_t a, value_t b) -> value_t {
		return base(a, b) * (1. / (this->_h_y * this->_h_y));
	};
	auto gamma = [&](value_t a, value_t b) -> value_t {
		return ((4. * b) / (2. * b + a * this->_h_time)) - 2. * alpha1(a, b) - 2. * alpha2(a, b);
	};
	auto betta = [&](value_t a, value_t b) -> value_t {
		return (a * this->_h_time - 2. * b) / (2. * b + a * this->_h_time);
	};
	auto omega_func = [&](value_t t) ->value_t {
		return ((t / _tau) * (_omega_max - _omega_min) + _omega_min);
	};

	ind_t width = _x_vector_size * _y_vector_size;
	ind_t tmp_width = (this->_numberOfStepsX / _layerN_x + 1) * (this->_numberOfStepsY / _layerN_y + 1);
	value_t* tmp_mtrx;
	ind_t counter = 0;
	if (!_all_layers) {
		tmp_mtrx = (value_t*)calloc(tmp_width, sizeof(value_t));
	}

	ind_t l_prev_x;
	ind_t l_ind_x;
	ind_t l_prev_y;
	ind_t l_ind_y;
	value_t _A;
	double t_total = omp_get_wtime();
	double t_write;

	// Xi_1 FOR PML
	for (ind_t i = 0; i <= 2 * _pml_ind.front().first; i++)
	{
		_xi1[i] = (1. / (sqrt(_b.front()) * _pml_width.front()))
			* std::log10(1. / _r_reflection)
			* (
				abs((i * this->_h_x / 2. + this->_lborders_x) - _pml_crd.front().first) / _pml_width.front()
				- sin(2. * M_PI * abs((i * this->_h_x / 2. + this->_lborders_x)
					- _pml_crd.front().first) / _pml_width.front()) / (2. * M_PI)
				);
	}
	for (ind_t i = 2 * _pml_ind.front().second; i < xi1_size; i++)
	{
		_xi1[i] = (1. / (sqrt(_b.back()) * _pml_width.front()))
			* std::log10(1. / _r_reflection)
			* (
				abs((i * this->_h_x / 2. + this->_lborders_x) - _pml_crd.front().second) / _pml_width.front()
				- sin(2. * M_PI * abs((i * this->_h_x / 2. + this->_lborders_x)
					- _pml_crd.front().second) / _pml_width.front()) / (2. * M_PI)
				);
	}
	// Xi_2 FOR PML
	for (ind_t l = 0; l < _layers_ind.front().size(); l++)
	{
		for (ind_t j = 0; j <= 2 * _pml_ind[1].first; j++)
		{
			_xi2[l * xi2_size + j] = (1./ (sqrt(_b[l]) * _pml_width[1]))
				* std::log10(1. / _r_reflection)
				* (
					abs((j * this->_h_y / 2. + this->_lborders_y) - _pml_crd[1].first) / _pml_width[1]
					- sin(2. * M_PI * abs((j * this->_h_y / 2. + this->_lborders_y)
						- _pml_crd[1].first) / _pml_width[1]) / (2. * M_PI)
					);
		}
		for (ind_t j = 2 * _pml_ind[1].second; j < xi2_size; j++)
		{
			_xi2[l * xi2_size + j] = (1. / (sqrt(_b[l]) * _pml_width[1]))
				* std::log10(1. / _r_reflection)
				* (
					abs((j * this->_h_y / 2. + this->_lborders_y) - _pml_crd[1].second) / _pml_width[1]
					- sin(2. * M_PI * abs((j * this->_h_y / 2. + this->_lborders_y)
						- _pml_crd[1].second) / _pml_width[1]) / (2. * M_PI)
					);
		}
	}
	//
	for (ind_t i = 0; i < _x_vector_size; i++)
	{
		for (ind_t j = 0; j < _y_vector_size; j++)
		{
			mtrx[i * _y_vector_size + j] = func_t0(i * this->_h_x + this->_lborders_x, j * this->_h_y + this->_lborders_y);
			mtrx[1 * width + i * _y_vector_size + j] = mtrx[i * _y_vector_size + j] + func_u_t0(i * this->_h_x + this->_lborders_x, j * this->_h_y + this->_lborders_y);
			
		}
	}
	
	l_prev_x = static_cast<ind_t>(1.);
	l_ind_x = static_cast<ind_t>(0.);
	for (const auto& lx : _layers_ind.front())
	{
		l_prev_y = static_cast<ind_t>(1.);
		l_ind_y = static_cast<ind_t>(0.);
		for (const auto& ly : _layers_ind[1])
		{
			if (l_ind_x == 0 || l_ind_y == 0
				|| l_ind_x == _layers_ind.front().size() - 1 || l_ind_y == _layers_ind[1].size() - 1)
			{
#pragma omp parallel for num_threads(_n_threads)
				for (ind_t i = l_prev_x; i < lx; i++)
				{
					for (ind_t j = l_prev_y; j < ly; j++)
					{

						_phi1[idx_phi(1) * phi1_size + i * (_y_vector_size - 1) + j] = _phi1[idx_phi(0) * phi1_size + i * (_y_vector_size - 1) + j] * ((2. - _xi1[2 * i + 1] * this->_h_time) / (2. + _xi1[2 * i + 1] * this->_h_time))
							+ (
								(((_xi2[l_ind_x * xi2_size + (2 * j + 1)] - _xi1[2 * i + 1]) * this->_h_time * (1. / _b[l_ind_x])) / (2. + this->_h_time * _xi1[2 * i + 1]))
								*
								(
									(mtrx[idx(1) * width + (i + 1) * _y_vector_size + j]
										+ mtrx[idx(1) * width + (i + 1) * _y_vector_size + (j + 1)]
										- mtrx[idx(1) * width + i * _y_vector_size + j]
										- mtrx[idx(1) * width + i * _y_vector_size + (j + 1)]) / (2. * this->_h_x)

									+ (mtrx[idx(0) * width + (i + 1) * _y_vector_size + j]
										+ mtrx[idx(0) * width + (i + 1) * _y_vector_size + (j + 1)]
										- mtrx[idx(0) * width + i * _y_vector_size + j]
										- mtrx[idx(0) * width + i * _y_vector_size + (j + 1)]) / (2. * this->_h_x)
									)
								);

						_phi2[idx_phi(1) * phi2_size + i * (_y_vector_size - 1) + j] = _phi2[idx_phi(0) * phi2_size + i * (_y_vector_size - 1) + j] * ((2. - _xi2[l_ind_x * xi2_size + (2 * j + 1)] * this->_h_time) / (2. + _xi2[l_ind_x * xi2_size + (2 * j + 1)] * this->_h_time))
							+ (
								(((_xi1[2 * i + 1] - _xi2[l_ind_x * xi2_size + (2 * j + 1)]) * this->_h_time * (1. / _b[l_ind_x])) / (2. + this->_h_time * _xi2[l_ind_x * xi2_size + (2 * j + 1)]))
								*
								(
									(mtrx[idx(1) * width + i * _y_vector_size + (j + 1)]
										+ mtrx[idx(1) * width + (i + 1) * _y_vector_size + (j + 1)]
										- mtrx[idx(1) * width + i * _y_vector_size + j]
										- mtrx[idx(1) * width + (i + 1) * _y_vector_size + j]) / (2. * this->_h_y)

									+ (mtrx[idx(0) * width + i * _y_vector_size + (j + 1)]
										+ mtrx[idx(0) * width + (i + 1) * _y_vector_size + (j + 1)]
										- mtrx[idx(0) * width + i * _y_vector_size + j]
										- mtrx[idx(0) * width + (i + 1) * _y_vector_size + j]) / (2. * this->_h_y)
									)
								);
					}	
				}
			}
			l_prev_y = ly;
			l_ind_y++;
		}
		l_prev_x = lx;
		l_ind_x++;
	}
	
	///////////////////////////////////////////////////////////////////////////
	std::cout << "Start sheme\n";
	for (ind_t p = 1; p < _t_vector_size - 1; p++)
	{
		//запись в файл трех слоев по времени
		if (_fwrite && !_all_layers && p % _layerN_t == 0) {
			//t_write = omp_get_wtime();
			if (_layerN_x > 1 || _layerN_y > 1) {
				counter = 0;
				for (ind_t i = 0; i < _x_vector_size; i++)
				{
					for (ind_t j = 0; j < _y_vector_size; j++)
					{
						if (i % _layerN_x == 0 && j % _layerN_y == 0) {
							tmp_mtrx[counter] = mtrx[idx(p) * width + i * _y_vector_size + j];
							counter++;
						}
					}
				}
				Write_to_extendible_H5(filename.c_str(), "wave", tmp_mtrx, 1, tmp_width);
			}
			else
			{
				Write_to_extendible_H5(filename.c_str(), "wave", mtrx, 1, width);
			}
			//t_write = omp_get_wtime() - t_write;
			//std::cout << "Write file : " << t_write << " seconds.\n";
		}
		double s_iter = omp_get_wtime();
		// граница
#pragma omp parallel for num_threads(_n_threads)
		for (ind_t j = 1; j < _y_vector_size - 1; j++)
		{
			mtrx[idx(p + 1) * width + 0 * _y_vector_size + j] = static_cast<value_t>(0.);
			mtrx[idx(p + 1) * width + (_x_vector_size - 1) * _y_vector_size + j] = static_cast<value_t>(0.);
		}

		//wave
		l_prev_x = static_cast<ind_t>(1.);
		l_ind_x = static_cast<ind_t>(0.);
		for (const auto& lx : _layers_ind.front())
		{
			l_prev_y = static_cast<ind_t>(1.);
			l_ind_y = static_cast<ind_t>(0.);
			for (const auto& ly : _layers_ind[1])
			{
#pragma omp parallel for num_threads(_n_threads)
				for (ind_t i = l_prev_x; i < lx; i++)
				{
					for (ind_t j = l_prev_y; j < ly; j++)
					{
						/*// импульс с временем "испускания" tau 
						if (i == _pml_ind.front().first 
							&& j * this->_h_y + this->_lborders_y > centr_gauss - radius_gauss
							&& j * this->_h_y + this->_lborders_y < centr_gauss + radius_gauss)
						{
							value_t x = i * this->_h_x + this->_lborders_x;
							value_t y = j * this->_h_y + this->_lborders_y;
							value_t t = p * this->_h_time;

							value_t omega_t = omega_func(t);
							mtrx[idx(p + 1) * width + i * _y_vector_size + j] = impulse(x, y, t, omega_t, centr_gauss, radius_gauss);
						}
						else {*/
							//pml
							if (l_ind_x == 0 || l_ind_y == 0
								|| l_ind_x == _layers_ind.front().size() - 1 || l_ind_y == _layers_ind[1].size() - 1)
							{

								_A = (_xi1[2 * i] + _xi2[l_ind_x * xi2_size + (2 * j)]) * _b[l_ind_x];
								//pml
								mtrx[idx(p + 1) * width + i * _y_vector_size + j] = alpha1(_A, _b[l_ind_x]) * mtrx[idx(p) * width + (i - 1) * _y_vector_size + j]
									+ alpha1(_A, _b[l_ind_x]) * mtrx[idx(p) * width + (i + 1) * _y_vector_size + j]
									+ alpha2(_A, _b[l_ind_x]) * mtrx[idx(p) * width + i * _y_vector_size + (j - 1)]
									+ alpha2(_A, _b[l_ind_x]) * mtrx[idx(p) * width + i * _y_vector_size + (j + 1)]
									+ gamma(_A, _b[l_ind_x]) * mtrx[idx(p) * width + i * _y_vector_size + j]
									- mtrx[idx(p) * width + i * _y_vector_size + j] * base(_A, _b[l_ind_x]) * _xi1[2 * i] * _xi2[l_ind_x * xi2_size + (2 * j)] * _b[l_ind_x]
									+ betta(_A, _b[l_ind_x]) * mtrx[idx(p - 1) * width + i * _y_vector_size + j]
									+ base(_A, _b[l_ind_x]) * _b[l_ind_x]
									* (
										((_phi1[idx_phi(p) * phi1_size + i * (_y_vector_size - 1) + (j - 1)]
											+ _phi1[idx_phi(p) * phi1_size + i * (_y_vector_size - 1) + j]
											- _phi1[idx_phi(p) * phi1_size + (i - 1) * (_y_vector_size - 1) + (j - 1)]
											- _phi1[idx_phi(p) * phi1_size + (i - 1) * (_y_vector_size - 1) + j])
											/ (2. * this->_h_x))
										+
										((_phi2[idx_phi(p) * phi2_size + (i - 1) * (_y_vector_size - 1) + j]
											+ _phi2[idx_phi(p) * phi2_size + i * (_y_vector_size - 1) + j]
											- _phi2[idx_phi(p) * phi2_size + (i - 1) * (_y_vector_size - 1) + (j - 1)]
											- _phi2[idx_phi(p) * phi2_size + i * (_y_vector_size - 1) + (j - 1)])
											/ (2. * this->_h_y))
										);
							}
							//wave
							else
							{
								mtrx[idx(p + 1) * width + i * _y_vector_size + j] = alpha1(_a[l_ind_x - 1], _b[l_ind_x]) * mtrx[idx(p) * width + (i - 1) * _y_vector_size + j]
									+ alpha1(_a[l_ind_x - 1], _b[l_ind_x]) * mtrx[idx(p) * width + (i + 1) * _y_vector_size + j]
									+ alpha2(_a[l_ind_x - 1], _b[l_ind_x]) * mtrx[idx(p) * width + i * _y_vector_size + (j - 1)]
									+ alpha2(_a[l_ind_x - 1], _b[l_ind_x]) * mtrx[idx(p) * width + i * _y_vector_size + (j + 1)]
									+ gamma(_a[l_ind_x - 1], _b[l_ind_x]) * mtrx[idx(p) * width + i * _y_vector_size + j]
									+ betta(_a[l_ind_x - 1], _b[l_ind_x]) * mtrx[idx(p - 1) * width + i * _y_vector_size + j];
							}
						//}
					}
				}
				l_prev_y = ly;
				l_ind_y++;
			}
			// сдвиг левой границы на одну. Для пропуска i=границе раздела двух слоев.
			if (l_ind_x != 0 && l_ind_x != _layers_ind.front().size() - 2 && l_ind_x != _layers_ind.front().size() - 1)
			{
				l_prev_x = lx;
			}
			else
			{
				l_prev_x = lx;
			}
			l_ind_x++;
		}

		// условие на границе раздела
		for (ind_t lx = 0; lx < _layers_ind.front().size(); lx++)
		{
			ind_t i = _layers_ind.front()[lx];
			if (lx != 0 && lx != _layers_ind.front().size() - 2 && lx != _layers_ind.front().size() - 1)
			{
#pragma omp parallel for num_threads(_n_threads)
				for (ind_t j = 1; j < _y_vector_size - 1; j++)
				{
					mtrx[idx(p + 1) * width + i * _y_vector_size + j] = (1. / 6) * (-mtrx[idx(p + 1) * width + (i - 2) * _y_vector_size + j]
						+ 4. * mtrx[idx(p + 1) * width + (i - 1) * _y_vector_size + j]
						+ 4. * mtrx[idx(p + 1) * width + (i + 1) * _y_vector_size + j]
						- mtrx[idx(p + 1) * width + (i + 2) * _y_vector_size + j]);
				}
			}
		}
#pragma omp parallel for num_threads(_n_threads)
		// граница
		for (ind_t i = 1; i < _x_vector_size - 1; i++)
		{
			mtrx[idx(p + 1) * width + i * _y_vector_size + 0] = static_cast<value_t>(0.);
			mtrx[idx(p + 1) * width + i * _y_vector_size + (_y_vector_size - 1)] = static_cast<value_t>(0.);
		}
		//phi
		l_prev_x = static_cast<ind_t>(1.);
		l_ind_x = static_cast<ind_t>(0.);
		for (const auto& lx : _layers_ind.front())
		{
			l_prev_y = static_cast<ind_t>(1.);
			l_ind_y = static_cast<ind_t>(0.);
			for (const auto& ly : _layers_ind[1])
			{
				if (l_ind_x == 0 || l_ind_y == 0
					|| l_ind_x == _layers_ind.front().size() - 1 || l_ind_y == _layers_ind[1].size() - 1)
				{
#pragma omp parallel for num_threads(_n_threads)
					for (ind_t i = l_prev_x; i < lx; i++)
					{
						for (ind_t j = l_prev_y; j < ly; j++)
						{
							_phi1[idx_phi(p + 1) * phi1_size + i * (_y_vector_size - 1) + j] = _phi1[idx_phi(p) * phi1_size + i * (_y_vector_size - 1) + j] * ((2. - _xi1[2 * i + 1] * this->_h_time) / (2. + _xi1[2 * i + 1] * this->_h_time))
								+ (
									(((_xi2[l_ind_x * xi2_size + (2 * j + 1)] - _xi1[2 * i + 1]) * this->_h_time * (1. / _b[l_ind_x])) / (2. + this->_h_time * _xi1[2 * i + 1]))
									*
									(
										(mtrx[idx(p + 1) * width + (i + 1) * _y_vector_size + j]
											+ mtrx[idx(p + 1) * width + (i + 1) * _y_vector_size + (j + 1)]
											- mtrx[idx(p + 1) * width + i * _y_vector_size + j]
											- mtrx[idx(p + 1) * width + i * _y_vector_size + (j + 1)]) / (2. * this->_h_x)

										+ (mtrx[idx(p) * width + (i + 1) * _y_vector_size + j]
											+ mtrx[idx(p) * width + (i + 1) * _y_vector_size + (j + 1)]
											- mtrx[idx(p) * width + i * _y_vector_size + j]
											- mtrx[idx(p) * width + i * _y_vector_size + (j + 1)]) / (2. * this->_h_x)
										)
									);

							_phi2[idx_phi(p + 1) * phi2_size + i * (_y_vector_size - 1) + j] = _phi2[idx_phi(p) * phi2_size + i * (_y_vector_size - 1) + j] * ((2. - _xi2[l_ind_x * xi2_size + (2 * j + 1)] * this->_h_time) / (2. + _xi2[l_ind_x * xi2_size + (2 * j + 1)] * this->_h_time))
								+ (
									(((_xi1[2 * i + 1] - _xi2[l_ind_x * xi2_size + (2 * j + 1)]) * this->_h_time * (1. / _b[l_ind_x])) / (2. + this->_h_time * _xi2[l_ind_x * xi2_size + (2 * j + 1)]))
									*
									(
										(mtrx[idx(p + 1) * width + i * _y_vector_size + (j + 1)]
											+ mtrx[idx(p + 1) * width + (i + 1) * _y_vector_size + (j + 1)]
											- mtrx[idx(p + 1) * width + i * _y_vector_size + j]
											- mtrx[idx(p + 1) * width + (i + 1) * _y_vector_size + j]) / (2. * this->_h_y)

										+ (mtrx[idx(p) * width + i * _y_vector_size + (j + 1)]
											+ mtrx[idx(p) * width + (i + 1) * _y_vector_size + (j + 1)]
											- mtrx[idx(p) * width + i * _y_vector_size + j]
											- mtrx[idx(p) * width + (i + 1) * _y_vector_size + j]) / (2. * this->_h_y)
										)
									);

						}
					}
				}
				l_prev_y = ly;
				l_ind_y++;
			}
			l_prev_x = lx;
			l_ind_x++;
		}
		//

		if (p == 1)
		{
			double t_iter = omp_get_wtime() - s_iter;
			std::cout << "Time for one iteration: " << t_iter << " seconds.\n";
			std::cout << "Number of iterations:" << this->_numberOfStepsT << "\n";
			std::cout << "Time for all iteration: " << t_iter * this->_numberOfStepsT << " seconds.\n";
		}
		if (p % (this->_numberOfStepsT/100) == 0) std::cout << "End calculation for t = " << (p * 100) / (this->_numberOfStepsT - 1) << "%\n";
	}

	////запись в файл остальных слоев по времени
	//if (_fwrite && !_all_layers) {
	//	Write_to_extendible_H5(filename.c_str(), "wave", mtrx, idx(_t_vector_size), width);
	//}
	////
	t_total = omp_get_wtime() - t_total;
	std::cout << "Time Sheme : " << t_total << " seconds.\n";

	// запись всех слоев по времени
	if (_all_layers)
	{
		t_total = omp_get_wtime();
		WriteParallel(file, "wave", mtrx, _t_vector_size, width);
		t_total = omp_get_wtime() - t_total;
		std::cout << "Write file : " << t_total << " seconds.\n";
	}
	H5Fclose(file);
	free(tmp_mtrx);
	free(_xi1);
	free(_xi2);
	free(_phi1);
	free(_phi2);
}
//----------------------------------------------------------------------------------------------------------------------------------
template<typename value_t, typename ind_t, typename steps_t>
void Solver<value_t, ind_t, steps_t>::Write_to_extendible_H5(const char* FILENAME, const char* dset_name, value_t* data, ind_t rows, ind_t cols)
{
	hsize_t ndims = 2;
	hsize_t nrows = rows;
	hsize_t ncols = cols;

	hsize_t* mbuff_dims = (hsize_t*)malloc(ndims * sizeof(hsize_t));
	mbuff_dims[0] = nrows;
	mbuff_dims[1] = ncols;
	hid_t mem_space = H5Screate_simple(ndims, mbuff_dims, NULL);

	hid_t file = H5Fopen(FILENAME, H5F_ACC_RDWR, H5P_DEFAULT);

	if (file < 0)
	{
		exit(EXIT_FAILURE);
	}

	else
	{
		if (!H5Lexists(file, dset_name, H5P_DEFAULT))
		{

			hsize_t* dims = (hsize_t*)malloc(ndims * sizeof(hsize_t));
			dims[0] = nrows;
			dims[1] = ncols;
			
			hsize_t* max_dims = (hsize_t*)malloc(ndims * sizeof(hsize_t));
			max_dims[0] = H5S_UNLIMITED;
			max_dims[1] = ncols;

			hid_t file_space = H5Screate_simple(ndims, dims, max_dims);

			hid_t plist = H5Pcreate(H5P_DATASET_CREATE);
			H5Pset_layout(plist, H5D_CHUNKED);

			hsize_t* chunk_dims = (hsize_t*)malloc(ndims * sizeof(hsize_t));
			chunk_dims[0] = nrows;
			chunk_dims[1] = ncols;

			H5Pset_chunk(plist, ndims, chunk_dims);

			hid_t dset = H5Dcreate(file, dset_name, H5T_NATIVE_DOUBLE,
				file_space, H5P_DEFAULT, plist, H5P_DEFAULT);

			H5Pclose(plist);

			file_space = H5Dget_space(dset);
			hsize_t start[2] = { 0, 0 };
			hsize_t count[2] = { nrows, ncols };
			H5Sselect_hyperslab(file_space, H5S_SELECT_SET, start,
				NULL, count, NULL);
			for (hsize_t i = 0; i < nrows; ++i)
			{
				H5Dwrite(dset, H5T_NATIVE_DOUBLE, mem_space, file_space,
					H5P_DEFAULT, data);
			}

			H5Sclose(mem_space);
			H5Dclose(dset);
			H5Sclose(file_space);
			H5Fclose(file);

			free(dims);
			free(max_dims);
		}
		else
		{
			hid_t dset = H5Dopen(file, dset_name, H5P_DEFAULT);
			hid_t file_space = H5Dget_space(dset);

			hsize_t* dims = (hsize_t*)calloc(ndims, sizeof(hsize_t));

			H5Sget_simple_extent_dims(file_space, dims, NULL);

			dims[0] += nrows;
			dims[1] = ncols;
			H5Dset_extent(dset, dims);
			file_space = H5Dget_space(dset);

			hsize_t start[2] = { dims[0] - nrows, 0 };
			hsize_t count[2] = { nrows, ncols };
			H5Sselect_hyperslab(file_space, H5S_SELECT_SET, start,
				NULL, count, NULL);

			for (hsize_t i = 0; i < nrows; ++i)
			{
				H5Dwrite(dset, H5T_NATIVE_DOUBLE, mem_space, file_space,
					H5P_DEFAULT, data);
			}

			H5Sclose(mem_space);
			H5Dclose(dset);
			H5Sclose(file_space);
			H5Fclose(file);
			free(dims);
		}
	}
	free(mbuff_dims);
}
//----------------------------------------------------------------------------------------------------------------------------------
template<typename value_t, typename ind_t, typename steps_t>
void Solver<value_t, ind_t, steps_t>::Write_param(std::string filename, bool hdf5_flag) {
	std::string param_filename = _PATH + "/" + _FILE_NAME + "_params" + ".h5";
	hid_t param_file = H5Fcreate(param_filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
	steps_t* grid_settings = (steps_t*)calloc(this->_h.size() + this->_lborders.size() + this->_rborders.size(), sizeof(steps_t));
	//steps
	auto it = this->_h.begin();
	for (ind_t i = 0; i < this->_h.size(); i++)
	{
		grid_settings[i] = (it->second);
		it++;
	}
	//left borders
	it = this->_lborders.begin();
	for (ind_t i = this->_h.size(); i < this->_h.size() + this->_lborders.size(); i++)
	{
		grid_settings[i] = (it->second);
		it++;
	}
	it = this->_rborders.begin();
	for (ind_t i = this->_h.size() + this->_lborders.size(); i < this->_h.size() + this->_lborders.size() + this->_rborders.size(); i++)
	{
		grid_settings[i] = (it->second);
		it++;
	}
	// pml crd
	steps_t* pml_left = (steps_t*)calloc(_pml_crd.size(), sizeof(steps_t));
	steps_t* pml_right = (steps_t*)calloc(_pml_crd.size(), sizeof(steps_t));
	for (ind_t i = 0; i < _pml_crd.size(); i++)
	{
		pml_left[i] = _pml_crd[i].first;
		pml_right[i] = _pml_crd[i].second;
	}
	// skin crd
	steps_t* skin = (steps_t*)calloc(_skin_crd.size() * _skin_crd[0].size(), sizeof(steps_t));
	for (ind_t i = 0; i < _skin_crd.size(); i++)
	{
		for (ind_t j = 0; j < _skin_crd[i].size(); j++)
		{	
			skin[i * _skin_crd[0].size() + j] = _skin_crd[i][j];
		}
	}
	steps_t* N_l = (steps_t*)calloc(3, sizeof(steps_t));
	N_l[0] = _layerN_t;
	N_l[1] = _layerN_x;
	N_l[2] = _layerN_y;
	steps_t* light_param = (steps_t*)calloc(4, sizeof(steps_t));
	light_param[0] = _lamda;
	light_param[1] = _omega_min;
	light_param[2] = _omega_max;
	light_param[3] = _tau;


	Write_to_extendible_H5(param_filename.c_str(), "grid_settings", grid_settings, 1, this->_h.size() + this->_lborders.size() + this->_rborders.size());
	Write_to_extendible_H5(param_filename.c_str(), "pml_left", pml_left, 1, _pml_crd.size());
	Write_to_extendible_H5(param_filename.c_str(), "pml_right", pml_right, 1, _pml_crd.size());
	Write_to_extendible_H5(param_filename.c_str(), "skin_crd", skin, 1, _skin_crd.size() * _skin_crd[0].size());
	Write_to_extendible_H5(param_filename.c_str(), "write_layer_N", N_l, 1, 3);
	Write_to_extendible_H5(param_filename.c_str(), "light_param", light_param, 1, 4);

	free(grid_settings);
	free(pml_left);
	free(pml_right);
	free(skin);
	H5Fclose(param_file);
}
//----------------------------------------------------------------------------------------------------------------------------------
template<typename value_t, typename ind_t, typename steps_t>
void Solver<value_t, ind_t, steps_t>::WriteParallel(H5::H5File group, std::string dataSetName, value_t* data, ind_t nrows, ind_t ncols)
{
	hsize_t dims[2] = { nrows, ncols };
	H5::DataSpace dspace(2, dims);
	H5::DataSet dset(group.createDataSet(dataSetName, H5::PredType::NATIVE_DOUBLE, dspace));
	dset.write(data, H5::PredType::NATIVE_DOUBLE);
	dset.close();
}