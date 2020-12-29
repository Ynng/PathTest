// PathTest.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <vector>
#include <cstdlib>
#include <any>
#include <map>
#include <stdexcept>
#include <string>
#include <cmath>
#include <algorithm>
#include <iomanip>
#include <sstream>
#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

class Pos2d
{
private:
	const std::any& getID(const std::string& iid) const;
protected:
	std::map<std::string, std::any> data{};
public:
	Pos2d();

	Pos2d(double x, double y);

	double magnitude() const;

	double distance(const Pos2d& other) const;

	double distance(double x, double y) const;

	void addIndex(int index, double change);

	double norm();

	//TODO: Pos2d RotateBy(Rotation? & other);

	Pos2d operator+(const Pos2d& other) const;

	Pos2d& operator+=(const Pos2d& other);

	Pos2d operator-(const Pos2d& other) const;

	Pos2d& operator-=(const Pos2d& other);

	Pos2d operator-();

	Pos2d operator*(const double scalar) const;

	Pos2d& operator*=(double scalar);

	Pos2d operator/(const double scalar) const;

	bool operator==(const Pos2d& other) const;

	bool operator!=(const Pos2d& other) const;

	Pos2d& operator/=(double scalar);

	double operator[](int index);

	/**
	 * Set the point data.
	 *
	 * @param iid   The data name
	 * @param idata The data
	 */
	void setData(const std::string& iid, const std::any& idata);

	/**
	 * Get the point data.
	 *
	 * @param  iid The data name
	 * @tparam T   The data type
	 * @return The data
	 */
	template <typename T> T getData(const std::string& iid) const {
		const std::any& idata = getID(iid);
		try {
			return std::any_cast<T>(idata);
		}
		catch (const std::bad_any_cast& e) {
			throw std::runtime_error("Pos2d::getData:: \"" + iid + "\" contains wrong type \"" +
				idata.type().name() + "\"");
		}
	}

	double x = 0;
	double y = 0;
};
/*-----Pos2d-------*/
Pos2d::Pos2d()
{
	x = 0;
	y = 0;
}

Pos2d::Pos2d(double x, double y)
{
	this->x = x;
	this->y = y;
}

double Pos2d::magnitude() const {
	return std::sqrt(x * x + y * y);
}

double Pos2d::distance(const Pos2d& other) const
{
	return std::sqrt(pow(x - other.x, 2) + pow(y - other.y, 2));
}

double Pos2d::distance(double x, double y) const
{
	return std::sqrt(pow(this->x - x, 2) + pow(this->y - y, 2));
}

double Pos2d::norm()
{
	return std::sqrt(pow(x, 2) + pow(y, 2));
}

Pos2d Pos2d::operator+(const Pos2d& other) const
{
	return { x + other.x, y + other.y };
}

Pos2d& Pos2d::operator+=(const Pos2d& other)
{
	x += other.x;
	y += other.y;
	return *this;
}

Pos2d Pos2d::operator-(const Pos2d& other)const
{
	return { x - other.x, y - other.y };
}

Pos2d& Pos2d::operator-=(const Pos2d& other)
{
	x -= other.x;
	y -= other.y;
	return *this;
}

Pos2d Pos2d::operator-() { return { -x, -y }; }

Pos2d Pos2d::operator*(const double scalar) const
{
	return { scalar * x, scalar * y };
}

Pos2d& Pos2d::operator*=(double scalar)
{
	x *= scalar;
	y *= scalar;
	return *this;
}

Pos2d Pos2d::operator/(const double scalar) const
{
	return *this * (1.0 / scalar);
}

bool Pos2d::operator==(const Pos2d& other) const
{
	return abs(x - other.x) < 1e-9 &&
		abs(y - other.y) < 1e-9;
}

bool Pos2d::operator!=(const Pos2d& other) const
{
	return !operator==(other);
}

Pos2d& Pos2d::operator/=(double scalar)
{
	*this *= (1.0 / scalar);
	return *this;
}

double Pos2d::operator[](int index) {
	if (index == 0)return x;
	if (index == 1)return y;

	throw std::runtime_error("Vector::at():: \"" + std::to_string(index) + "\" is invalid index");
}

void Pos2d::addIndex(int index, double change) {
	if (index == 0)x += change;
	if (index == 1)y += change;
}

void Pos2d::setData(const std::string& iid, const std::any& idata) {
	data[iid] = idata;
}

const std::any& Pos2d::getID(const std::string& iid) const {
	try {
		return data.at(iid);
	}
	catch (const std::out_of_range& e) {
		throw std::runtime_error("DataPoint::getID:: \"" + iid + "\" does not exist in point");
	}
}

class QuinticPolynomial {
private:
	double coeffs[6];
public:
	QuinticPolynomial(double xstart, double vstart, double xend, double vend) {
		double u = xend - xstart - vstart;
		double v = vend - vstart;

		double a3 = 10 * u - 4 * v;
		double a4 = -15 * u + 7 * v;
		double a5 = 6 * u - 3 * v;

		coeffs[0] = xstart;
		coeffs[1] = vstart;
		coeffs[2] = 0;
		coeffs[3] = a3;
		coeffs[4] = a4;
		coeffs[5] = a5;
	}

	double calcPoint(double t) {
		double xt = 0;
		for (int power = 0; power < 6; power++) {
			xt += coeffs[power] * pow(t, power);
		}
		return xt;
	}
};

double minVel=5;
double maxAccel=60;
double maxVel=45;
double maxDecel=40;
double finalVel=1;
double k=10;

std::vector<Pos2d> quinticSegment(Pos2d s, Pos2d g, int steps, bool end) {
	std::vector<double> rx;
	std::vector<double> ry;
	rx.reserve(end ? steps + 1 : steps);
	ry.reserve(end ? steps + 1 : steps);

	double vxs = s.getData<double>("vel") * sin(s.getData<double>("theta"));
	double vys = s.getData<double>("vel") * cos(s.getData<double>("theta"));
	double vxg = g.getData<double>("vel") * sin(g.getData<double>("theta"));
	double vyg = g.getData<double>("vel") * cos(g.getData<double>("theta"));

	QuinticPolynomial xqp = QuinticPolynomial(s.x, vxs, g.x, vxg);
	QuinticPolynomial yqp = QuinticPolynomial(s.y, vys, g.y, vyg);

	for (double i = 0; i <= (end ? steps : steps - 1); i++) {
		rx.push_back(xqp.calcPoint(i / steps));
		ry.push_back(yqp.calcPoint(i / steps));
	}

	std::vector<Pos2d> path;
	path.reserve(rx.size());

	for (size_t i = 0; i < rx.size(); i++) {
		path.push_back({rx[i], ry[i]});
	}

	return path;
}

double wrapAngle(double angle) {
	return angle - (M_PI * 2) * std::floor((angle + M_PI) / (M_PI * 2));
};

double angleBetweenPointsSpline(Pos2d current, Pos2d target) {
	//return wrapAngle(atan2(target.y - current.y, target.x - current.x));
	return wrapAngle(atan2(target.x - current.x, target.y - current.y));
}

void computeCurvature(std::vector<Pos2d>& ipath) {
	//adding additional information
	ipath[0].setData("curvature", 0.0);
	for (size_t i = 1; i < ipath.size() - 1; i++)
	{
		double distOne = ipath[i].distance(ipath[i - 1]);
		double distTwo = ipath[i].distance(ipath[i + 1]);
		double distThree = ipath[i + 1].distance(ipath[i - 1]);

		double productOfSides = distOne * distTwo * distThree;
		double semiPerimeter = (distOne + distTwo + distThree) / 2.0;

		double triangleArea = sqrt(semiPerimeter * (semiPerimeter - distOne) * (semiPerimeter - distTwo) * (semiPerimeter - distThree));

		double r = productOfSides / (4.0 * triangleArea);
		double curvature = isnormal(1.0 / r) ? 1.0 / r : 0;
		curvature = curvature * curvature;

		ipath[i].setData("curvature", curvature);
	}
	ipath.back().setData("curvature", 0.0);
}

void computeVelocity(std::vector<Pos2d>& ipath) {
	ipath.back().setData("velocity", finalVel);
	for (int i = ipath.size() - 1; i > 0; i--)
	{
		// k / curvature, limited to max
		double wantedVel = std::min(maxVel, k / ipath[i].getData<double>("curvature"));

		//distance from last point
		double distance = ipath[i].distance(ipath[i - 1]);

		// maximum velocity given distance respecting acceleration
	// vf = sqrt(vi2 + 2ad)
		double maxIncrement = sqrt(pow(ipath[i].getData<double>("velocity"), 2) + (2.0 * maxDecel * distance));

		double newVel = std::min(wantedVel, maxIncrement);
		ipath[i - 1].setData("velocity", newVel);
	}
}

std::vector<Pos2d> withPathQuintic(std::vector<Pos2d> ipath, int steps = 50, double slopeScalar = 0.8) {
	//calculating angles
	ipath[0].setData("theta", angleBetweenPointsSpline(ipath[0], ipath[1]));
	for (size_t i = 1; i < ipath.size()-1; i++)
	{
		ipath[i].setData("theta", angleBetweenPointsSpline(ipath[i-1], ipath[i+1]));
	}
	ipath[ipath.size()-1].setData("theta", angleBetweenPointsSpline(ipath[ipath.size() - 2], ipath[ipath.size() - 1]));
	

	//generating velocities
	for (size_t i = 0; i < ipath.size()-1; i++)
	{
		double vel = slopeScalar * ipath[i].distance(ipath[i + 1]);
		ipath[i].setData("vel", vel);

		if(i==ipath.size()-2)
			ipath[i+1].setData("vel", vel);
	}

	//debuggin
	for (size_t i = 0; i < ipath.size(); i++) {
		printf("theta: %.2f\n", ipath[i].getData<double>("theta"));
	}


	//generating path
	std::vector<Pos2d> path(0);
	if (ipath.size() == 2) {
		path = quinticSegment(ipath[0], ipath[1], steps, true);
	}
	else {
		for (int i = 0; i < ipath.size() - 1; i++) {
			std::vector<Pos2d> segment = quinticSegment(ipath[i], ipath[i+1], steps, i >= ipath.size() - 2);
			for (size_t j = 0; j < segment.size(); j++)
			{
				segment[j].setData("segIndex", i);
			}
			path.insert(path.end(), segment.begin(), segment.end());
		}
	}

	computeCurvature(path);
	computeVelocity(path);

	return path;
}

std::vector<Pos2d> withPathBezier(std::vector<Pos2d> ipath) {
	double resolution = 0.75;
	double iweight = 0.001;
	double iTolerance = 1e-9;
	//generating path
	std::vector<Pos2d> rawpath;
	for (size_t i = 0; i < ipath.size() - 1; i++) {
		Pos2d start = ipath[i];
		Pos2d end = ipath[i + 1];

		Pos2d diff = end - start;
		int steps = ceil(diff.magnitude() / resolution);
		Pos2d step = diff / steps;
		rawpath.reserve(rawpath.size() + steps);
		for (int j = 0; j < steps; j++)
		{
			rawpath.push_back(start + step * j);
		}
	}
	//end on the last point of the original path
	if (ipath.size() > 0)  rawpath.push_back(ipath.back());

	//smoothing path
	std::vector<Pos2d> path;
	//copying path
	path.reserve(rawpath.size());
	for (Pos2d point : rawpath)
		path.push_back(point);

	double smoothWeight = 1.0 - iweight;
	double change = iTolerance;
	while (change >= iTolerance) {
		change = 0;
		for (size_t i = 1; i < rawpath.size() - 1; i++)
		{
			for (int j = 0; j < 2; j++)
			{
				double aux = path[i][j];
				double dataFac = iweight * (rawpath[i][j] - aux);
				double smoothFac = smoothWeight * (path[i - 1][j] + path[i + 1][j] - 2.0 * aux);
				path[i].addIndex(j, dataFac + smoothFac);
				change = abs(aux - path[i][j]);
			}
		}
	}

	computeCurvature(path);
	computeVelocity(path);

	return path;
}
/* For matplot lib only */
std::string int_to_hex(int i)
{
	std::stringstream stream;
	stream << "0x"
		<< std::setfill('0') << std::setw(sizeof(int) * 2)
		<< std::hex << i;
	return stream.str();
}

std::string perc2multcolor(double perc, double min, double max) {
	double base = max - min;
	if (base == 0) {
		perc = 0;
	}
	else {
		perc = (perc - min) / base * 100;
	}
	int r, g, b = 0;
	if (perc >= 0 && perc <= 20) {
		r = 255;
		g = round(12.75 * perc);
		b = 0;
	}
	else if (perc > 20 && perc <= 40) {
		r = round(-12.75 * perc + 510);
		g = 255;
		b = 0;
	}
	else if (perc > 40 && perc <= 60) {
		r = 0;
		g = 255;
		b = round(12.75 * perc) - 510;
	}
	else if (perc > 60 && perc <= 80) {
		r = 0;
		g = round(-12.75 * perc + 1020);
		b = 255;
	}
	else {
		r = round(12.75 * perc - 1020);
		g = 0;
		b = 255;
	}
	int h = r * 0x10000 + g * 0x100 + b * 0x1;
	std::string temp = "000000" + int_to_hex(h);
	temp = temp.substr(temp.size() - 6);
	return "#" + temp;
}

int main()
{
	plt::figure_size(700, 700);
	
	//std::vector<Pos2d> rawPath = { {24,72} , {24,108},  {36,144}, {72, 120}, {72, 96}, {48, 72} };
	std::vector<Pos2d> rawPath = { {0,0} , {24,24},  {26,-24}, {48,0} };
	//std::vector<Pos2d> rawPath = { {24,72} , {24,108},  {36,144}};
	int n = rawPath.size();
	std::vector<double> xr(n), yr(n);
	for (int i = 0; i < rawPath.size(); i++)
	{
		xr[i] = rawPath[i].x;
		yr[i] = rawPath[i].y;
	}
	plt::plot(xr, yr, "go--");

	std::vector<Pos2d> path = withPathBezier(rawPath);
	std::vector<Pos2d> quinticPath = withPathQuintic(rawPath);

	n = path.size();
	std::vector<double> x(1), y(1), c(1), v(1);

	for (int i = 0; i < path.size(); i++)
	{
		x[0] = path[i].x;
		y[0] = path[i].y;
		c[0] = path[i].getData<double>("curvature");
		v[0] = path[i].getData<double>("velocity");
		plt::scatter(x, y, 20, { {"color", perc2multcolor(v[0], minVel, maxVel)} });
		printf("%d:\t%.2f %.2f\t%.2f\t|\t%s\n", i, x[0], y[0], v[0], perc2multcolor(v[0], minVel, maxVel));
	}

	for (int i = 0; i < quinticPath.size(); i++)
	{
		x[0] = quinticPath[i].x;
		y[0] = quinticPath[i].y;
		c[0] = quinticPath[i].getData<double>("curvature");
		v[0] = quinticPath[i].getData<double>("velocity");
		plt::scatter(x, y, 20, { {"color", perc2multcolor(v[0], minVel, maxVel)} });
		printf("%d:\t%.2f %.2f\t%.2f\t|\t%s\n", i, x[0], y[0], v[0], perc2multcolor(v[0], minVel, maxVel));
	}

	plt::axis("equal");
	plt::show();
}


// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
