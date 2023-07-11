#pragma once
#include <iostream>
#include <map>
#include <vector>
#include <string>
enum class Dimension { One, Two, Three };
/**
 * @class Grid
 * The Grid class contains information about space discretization.
 */
template<typename ind_t, typename steps_t>
class Grid
{
public:
	Grid(std::map<std::string, steps_t>& lborders, std::map<std::string, steps_t>& rborders, std::map<std::string, steps_t>& h,
		Dimension dim)
	{
		_lborders = lborders;
		_rborders = rborders;
		_h = h;
		_dim = dim;

		if (_dim == Dimension::One)
		{
			_h_x = _h["x"];
			_h_y = NULL;
			_h_z = NULL;
			_h_time = _h["time"];

			_lborders_x = _lborders["x"];
			_lborders_y = NULL;
			_lborders_z = NULL;
			_lborders_time = _lborders["time"];

			_rborders_x = _rborders["x"];
			_rborders_y = NULL;
			_rborders_z = NULL;
			_rborders_time = _rborders["time"];
		}
		else if (_dim == Dimension::Two)
		{
			_h_x = _h["x"];
			_h_y = _h["y"];
			_h_z = NULL;
			_h_time = _h["time"];

			_lborders_x = _lborders["x"];
			_lborders_y = _lborders["y"];
			_lborders_z = NULL;
			_lborders_time = _lborders["time"];

			_rborders_x = _rborders["x"];
			_rborders_y = _rborders["y"];
			_rborders_z = NULL;
			_rborders_time = _rborders["time"];
		}
		else
		{
			_h_x = _h["x"];
			_h_y = _h["y"];
			_h_z = _h["z"];
			_h_time = _h["time"];

			_lborders_x = _lborders["x"];
			_lborders_y = _lborders["y"];
			_lborders_z = _lborders["z"];
			_lborders_time = _lborders["time"];

			_rborders_x = _rborders["x"];
			_rborders_y = _rborders["y"];
			_rborders_z = _rborders["z"];
			_rborders_time = _rborders["time"];
		}
		CalcStepsNumber();
	};
	/*!
	Calculation of the number of points/iterations depending on the dimension
	\param[in] dim dimension of space
	*/
	void CalcStepsNumber() 
	{
		if (_dim == Dimension::One)
		{
			_numberOfStepsX = (_rborders["x"] - _lborders["x"]) / _h["x"];
			_numberOfStepsT = (_rborders["time"] - _lborders["time"]) / _h["time"];
			_numberOfStepsY = -1;
			_numberOfStepsZ = -1;
		}
		else if (_dim == Dimension::Two)
		{
			_numberOfStepsX = (_rborders["x"] - _lborders["x"]) / _h["x"];
			_numberOfStepsY = (_rborders["y"] - _lborders["y"]) / _h["y"];
			_numberOfStepsT = (_rborders["time"] - _lborders["time"]) / _h["time"];
			_numberOfStepsZ = -1;
		}
		else
		{
			_numberOfStepsX = (_rborders["x"] - _lborders["x"]) / _h["x"];
			_numberOfStepsY = (_rborders["y"] - _lborders["y"]) / _h["y"];
			_numberOfStepsT = (_rborders["time"] - _lborders["time"]) / _h["time"];
			_numberOfStepsZ = (_rborders["z"] - _lborders["z"]) / _h["z"];
		}

	};
	/*!
	Returns the number of points/iteration by X-space
	*/
	ind_t getNumberOfStepsX() { return _numberOfStepsX; };
	/*!
	Returns the number of points/iteration by Y-space
	*/
	ind_t getNumberOfStepsY() { return _numberOfStepsY; };
	/*!
	Returns the number of points/iteration by Z-space
	*/
	ind_t getNumberOfStepsZ() { return _numberOfStepsZ; };
	/*!
	Returns the number of points/iteration by T-space
	*/
	ind_t getNumberOfStepsT() { return _numberOfStepsT; };
protected:
	/// steps _h["x"] == hx ..
	std::map<std::string, steps_t> _h;
	/// _lborders["x"] == left borders x
	std::map<std::string, steps_t> _lborders;
	/// _rborders["x"] == right borders x
	std::map<std::string, steps_t> _rborders;
	//
	steps_t _h_x;
	steps_t _h_y;
	steps_t _h_z;
	steps_t _h_time;
	///
	steps_t _lborders_x;
	steps_t _lborders_y;
	steps_t _lborders_z;
	steps_t _lborders_time;
	//
	steps_t _rborders_x;
	steps_t _rborders_y;
	steps_t _rborders_z;
	steps_t _rborders_time;
	///
	/// iterations count of X-space
	ind_t _numberOfStepsX;
	/// iterations count of Y-space
	ind_t _numberOfStepsY;
	/// iterations count of Z-space
	ind_t _numberOfStepsZ;
	/// iterations count of Time
	ind_t _numberOfStepsT;
	Dimension _dim;
};

