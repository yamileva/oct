#include "CMakeWave.h"

double mu_0 = 1.25663706212e-06;
double eps_0 = 8.8541878128e-12;
double lightspeed = 2.998e8; // lightspeed. measured in [meter/seconds]

double Amplitude(double t, double tau, double delta) {
	if (t <= delta)
		return (t / delta);
	else if (t > delta && t <= tau - delta)
		return 1;
	else if (t > tau - delta and t <= tau)
		return (1 - (t - tau + delta) * 1 / delta);
	else return 0;
}
double func_t0(double x) {
	/*double delta_x = 30.;
	double sigma = 0.8;
	return (10. / (sqrt(M_PI) * sigma)) * exp(-((delta_x - x) * (delta_x - x)) / (sigma * sigma));*/
	return 0;
}
double impulse(double x, double t, double omega, double tau) {
	return Amplitude(t, tau, t / 4) * sin(omega * t);
}
double func_u_t0(double x, double t, double omega) {
	return omega;
}

double func2_t0(double x, double y) {
	double delta_x = (35. / 1300.) * 1e3, delta_y = (35. / 1300.) * 1e3;
	double sigma_x = 0.5, sigma_y = 0.2;
	//return exp(-((delta_x - x) * (delta_x - x) + (delta_y - y) * (delta_y - y)) / (sigma_x * sigma_y));
	return 0;
}
double impulse2(double x, double y, double t, double omega, double delta_y, double radius) {
	return exp(-((y - delta_y) * (y - delta_y)) / (radius * radius)) * sin(omega * t);
}
double func2_u_t0(double x, double y) {
	return 0.;
}


void set_get_params(const double& lamda, const double& period,
	std::string& name, std::map<std::string, double>& left_borders, std::map<std::string, double>& right_borders,
	std::map<std::string, double>& h,
	std::vector<std::pair<double, double> >& pml_crd,
	std::vector<std::vector<double> >& skin_crd,
	Dimension& dim,
	std::vector<double>& grid_settings,
	std::vector<double>& a,
	std::vector<double>& b,
	std::vector<double>& epsilon,
	std::vector<double>& sigma,
	std::vector<std::vector<double>>& pml_setting)
{

	/*CPython::Call_Cpython py("appr", "cpp_dielectric_params", "D:/USATU/OCT/Wave/src/Cpython");
	std::vector<std::vector<double> > diel_params;
	std::cout << "Get Params\n";
	py.dielectric_params(lamda, diel_params);*/
	//std::ranges::copy(diel_params.front(), std::ostream_iterator<double>(std::cout, " "));
	//std::cout << "\n";
	//std::ranges::copy(diel_params.back(), std::ostream_iterator<double>(std::cout, " "));
	//std::cout << "\n";


	if (dim == Dimension::One)
	{
		//std::vector<std::vector<double> > diel_params({ {2.3408999999983005,1.795599999954565,1.9599999999910123,1.9320999999531239,1.959999999987876}
		//										,{0.0665302480707485,0.3012867941819113,0.13999913997329247,0.3174449114081058,0.16260148493089227} });

		//sigma.push_back(0.);//VACUUM
		//sigma.insert(sigma.end(), diel_params.back().begin(), diel_params.back().end());

		//epsilon.push_back(1.); //PML
		//epsilon.push_back(1.); //VACUUM
		//epsilon.insert(epsilon.end(), diel_params.front().begin(), diel_params.front().end());
		//epsilon.push_back(diel_params.front().back()); //PML
		//
		//a = sigma;
		//for (std::size_t i = 0; i < sigma.size(); i++)
		//{
		//	a[i] = (mu_0 * sigma[i] * (lamda * 1e-9) * (lamda * 1e-9)) / period;
		//}
		//b = epsilon;
		//for (std::size_t i = 0; i < epsilon.size(); i++)
		//{
		//	b[i] = (mu_0 * eps_0 * epsilon[i] * (lamda * 1e-9) * (lamda * 1e-9)) / (period * period);
		//}

		//name = "wave1d";
		//// TODO: Поправить задание параметров, некрасиво и запутанно!!!
		//left_borders["x"] = 0.;
		//left_borders["time"] = 0.;

		//double lx = 460; //[mkm] 460
		//right_borders["x"] = (lx / lamda) * 1e3; // lx * 10^-6 [m] / lamda * 10^-9 [m]
		//right_borders["time"] = (2 * lx * 1e-6) / (lightspeed * period);

		//h["x"] = 1. / 40.;
		//h["time"] = 1. / 50;
		//
		//pml_crd.push_back(std::make_pair((20. / lamda) * 1e3, (440. / lamda) * 1e3)); //440
		//std::vector<double> layers = { 40., 60., 100., 280., 370. };
		//
		//for (std::size_t i = 0; i < layers.size(); i++)
		//{
		//	layers[i] = (layers[i] / lamda) * 1e3;
		//}
		//skin_crd.push_back(layers);
		//layers.clear();

		std::vector<std::vector<double> > diel_params({ {1}
												,{0} });


		sigma.insert(sigma.end(), diel_params.back().begin(), diel_params.back().end());

		epsilon.push_back(diel_params.front().back()); //PML
		epsilon.insert(epsilon.end(), diel_params.front().begin(), diel_params.front().end());
		epsilon.push_back(diel_params.front().back()); //PML

		a = sigma;
		for (std::size_t i = 0; i < sigma.size(); i++)
		{
			a[i] = (mu_0 * sigma[i] * (lamda * 1e-9) * (lamda * 1e-9)) / period;
		}
		b = epsilon;
		for (std::size_t i = 0; i < epsilon.size(); i++)
		{
			b[i] = (mu_0 * eps_0 * epsilon[i] * (lamda * 1e-9) * (lamda * 1e-9)) / (period * period);
		}

		name = "1d_test";
		// TODO: Поправить задание параметров, некрасиво и запутанно!!!
		left_borders["x"] = 0.;
		left_borders["time"] = 0.;

		double lx = 80; //[mkm] 460
		right_borders["x"] = (lx / lamda) * 1e3; // lx * 10^-6 [m] / lamda * 10^-9 [m]
		right_borders["time"] = (lx * 1e-6) / (lightspeed * period);

		h["x"] = 1. / 40.;
		h["time"] = 1. / 50;

		pml_crd.push_back(std::make_pair((20. / lamda) * 1e3, (60. / lamda) * 1e3)); //440
		std::vector<double> layers = { 0 };

		for (std::size_t i = 0; i < layers.size(); i++)
		{
			layers[i] = (layers[i] / lamda) * 1e3;
		}
		skin_crd.push_back(layers);
		layers.clear();

	}
	else if (dim == Dimension::Two)
	{
		std::vector<std::vector<double> > diel_params({ {1.7955990820}
												,{32.9423183078} });

		sigma.push_back(0.);//VACUUM
		sigma.insert(sigma.end(), diel_params.back().begin(), diel_params.back().end());

		epsilon.push_back(1.); //PML
		epsilon.push_back(1.); //VACUUM
		epsilon.insert(epsilon.end(), diel_params.front().begin(), diel_params.front().end());
		epsilon.push_back(diel_params.front().back()); //PML

		a = sigma;
		for (std::size_t i = 0; i < sigma.size(); i++)
		{
			a[i] = (mu_0 * sigma[i] * (lamda * 1e-9) * (lamda * 1e-9)) / period;
		}
		b = epsilon;
		for (std::size_t i = 0; i < epsilon.size(); i++)
		{
			b[i] = (mu_0 * eps_0 * epsilon[i] * (lamda * 1e-9) * (lamda * 1e-9)) / (period * period);
		}
		time_t now = time(0);
		tm* ltm = localtime(&now);

		name = "Save/hm" + std::to_string(ltm->tm_hour) + "_" + std::to_string(ltm->tm_min);
		left_borders["x"] = 0.;
		left_borders["y"] = 0.;
		left_borders["time"] = 0.;

		double lx = 110; //[mkm] 460
		double ly = 70; //[mkm] 60
		right_borders["x"] = (lx / lamda) * 1e3; // lx * 10^-6 [m] / lamda * 10^-9 [m]
		right_borders["y"] = (ly / lamda) * 1e3; // lx * 10^-6 [m] / lamda * 10^-9 [m]
		right_borders["time"] = (3. * lx * 1e-6) / (lightspeed * period);

		h["x"] = 1. / 40.;
		h["y"] = 1. / 40.;
		h["time"] = 1. / 60.;

		pml_crd.push_back(std::make_pair((10. / lamda) * 1e3, (100. / lamda) * 1e3)); //440
		pml_crd.push_back(std::make_pair((10. / lamda) * 1e3, (60. / lamda) * 1e3));
		//std::vector<double> layers = { 40., 60., 100., 280., 370. };
		std::vector<double> layers = { 70. };  // distance to 2nd layer

		for (std::size_t i = 0; i < layers.size(); i++)
		{
			layers[i] = (layers[i] / lamda) * 1e3;
		}
		skin_crd.push_back(layers);
		skin_crd.push_back({});
		layers.clear();
	}

	for (const auto& el : h) {
		grid_settings.push_back(el.second);
	}
	for (const auto& el : left_borders) {
		grid_settings.push_back(el.second);
	}
	for (const auto& el : right_borders) {
		grid_settings.push_back(el.second);
	}
	for (const auto& el : pml_crd) {
		pml_setting.push_back({ el.first,el.second });
	}
}

int main(int argc, char* argv[]) {

	/*std::cout << "Avialable threads " << omp_get_max_threads() << "\n" << omp_in_parallel() << std::endl;
#pragma omp parallel num_threads(12)
	{
		std::cout << omp_get_thread_num() << " " <<  omp_in_parallel() << std::endl;
	}*/
	int threads_N = 4;  //4
	bool solve = true;
	bool write_file = true;
	bool write_all_layers = false;
	Dimension dim;
	dim = Dimension::Two;

	// parametrs
	double lamda = 1000;//lambda - wavelength. measured in [nanometer]
	double delta_lamda = 0; //lambda change on sourse (+-)
	double period = (lamda * 1e-9) / lightspeed;
	double lamda_min = lamda - delta_lamda, lamda_max = lamda + delta_lamda;
	double omega_min = period * 2. * M_PI * lightspeed / (lamda_max * 1e-9);
	double omega_max = period * 2. * M_PI * lightspeed / (lamda_min * 1e-9);
	int Nt = 5, Nx = 5, Ny = 5;
	double tau = 0.2 * 1e-12 / period;

	double omega_0 = 0.5;
	double delta_y = (30. / lamda) * 1e3;//4

	std::string name;
	std::map<std::string, double> left_borders;
	std::map<std::string, double> right_borders;
	std::map<std::string, double> h;
	double r_reflection = 1e-10;
	std::vector<std::pair<double, double> > pml_crd;
	std::vector<std::vector<double> > skin_crd;
	std::vector<double> grid_settings;
	std::vector<std::vector<double>> pml_setting;
	std::vector<double> a;
	std::vector<double> b;
	std::vector<double> epsilon;
	std::vector<double> sigma;
	//

	set_get_params(lamda, period, name, left_borders, right_borders, h, pml_crd, skin_crd, dim, grid_settings, a, b, epsilon, sigma, pml_setting);
	Solver<double, int, double> wave(left_borders, right_borders, h, r_reflection, pml_crd,
		skin_crd, a, b, epsilon, sigma, dim, write_file, write_all_layers, name, threads_N);


	if (dim == Dimension::One)
	{
		if (solve) {
			std::size_t T, X;
			std::cout << "Start resize. Number of steps T: " << wave.getNumberOfStepsT()
				<< ". Number of steps X: " << wave.getNumberOfStepsX() << "\n";
			if (write_file && !write_all_layers)
			{
				T = 3;
			}
			else
			{
				T = (wave.getNumberOfStepsT() + 1);
			}
			X = wave.getNumberOfStepsX() + 1;
			std::cout << "Size mtrx = " << T * X * sizeof(double) * 1e-9 << " GB" << "\n";
			double* mtrx = (double*)calloc(T * X, sizeof(double));
			double tau = (5. * lamda * 1e-9) / (lightspeed * period);
			double omega = period * 2. * M_PI * lightspeed / (lamda * 1e-9);
			wave.Wave1d(mtrx, tau, omega, func_t0, func_u_t0, impulse);
			free(mtrx);
		}
	}
	else if (dim == Dimension::Two)
	{
		if (solve) {
			std::size_t T, X, Y;
			std::cout << "Start resize. Number of steps T: " << wave.getNumberOfStepsT()
				<< ". Number of steps X: " << wave.getNumberOfStepsX()
				<< ". Number of steps Y: " << wave.getNumberOfStepsY() << "\n";
			if (write_file && !write_all_layers)
			{
				T = 3;
			}
			else
			{
				T = (wave.getNumberOfStepsT() + 1);
			}
			X = wave.getNumberOfStepsX() + 1;
			Y = wave.getNumberOfStepsY() + 1;


			std::cout << "Size file = " << (wave.getNumberOfStepsT() + 1) / Nt * X / Nx * Y / Ny * sizeof(double) * 1e-9 << " GB" << "\n";
			std::cout << "Size mtrx = " << T * X * Y * sizeof(double) * 1e-9 << " GB" << "\n";

			double* mtrx = (double*)calloc(T * X * Y, sizeof(double));

			wave.Wave2d(mtrx, lamda, omega_min, omega_max, tau, func2_t0, func2_u_t0, impulse2, Nt, Nx, Ny, delta_y, omega_0);
			free(mtrx);
		}
	}

	/*Graphics::CGraph graph("pyqt", "cpp_app", "D:/USATU/OCT/PML/src/Graphics");
	if (write_file)
	{

		graph.plot(dim, name + ".h5");
	}
	else
	{

		graph.plot(dim, mtrx, grid_settings, pml_setting, xi);
	}*/


	return EXIT_SUCCESS;
}