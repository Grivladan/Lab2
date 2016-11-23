#include "stdafx.h"

inline float Vector2f::length() const
{
	return sqrt(x * x + y * y);
}

inline float Vector2f::operator*(Vector2f v) const
{
	return x * v.x + y * v.y;
}

inline Vector2f Vector2f::operator+(Vector2f v) const
{
	return Vector2f(x + v.x, y + v.y);
}

Vector2f Vector2f::operator-(Vector2f v) const
{
	return Vector2f(x - v.x, y - v.y);
}

inline Vector2f Vector2f::operator*(float f) const
{
	return Vector2f(x * f, y * f);
}

Vector2f Vector2f::normalize() const
{
	Vector2f v = *(this);
	
	float l = v.length();

	v.x /= l;
	v.y /= l;

	return v;
}

inline float Vector2f::area(Vector2f a, Vector2f b, Vector2f c)
{
	return (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
}

inline bool Vector2f::intersect_1(float a, float b, float c, float d)
{
	if (a > b) std::swap(a, b);
	if (c > d) std::swap(c, d);
	return std::max(a, c) <= std::min(b, d);
}

inline bool Vector2f::intersect(Vector2f a, Vector2f b, Vector2f c, Vector2f d)
{
	return intersect_1(a.x, b.x, c.x, d.x)
		&& intersect_1(a.y, b.y, c.y, d.y)
		&& area(a, b, c) * area(a, b, d) <= 0
		&& area(c, d, a) * area(c, d, b) <= 0;
}

void FlowSolver::splitCurve()
{
	if (M < 4)
		throw "Not enough singularities";

	singularities.clear();
	colocations.clear();
	normals.clear();

	float l1 = 8.0f / 9.0f, l2 = 8.0f / 9.0f;
	int m1, m2;
	float L = l1 + l2;
	m1 = std::max(1, int(l1 * (M - 1) / L));
	m2 = (M - 1) - m1;

	big_delta = std::min(l1 / m1, l2 / m2);
	delta = big_delta * 0.45f;

	singularities.push_back(Vector2f(-2.0f / 9.0f, 2.0f / 9.0f));

	for (int i = 0; i < m1; i++)
	{
		float x, y;

		x = -2.0f / 9.0 + (i + 1) * (2.0f / 9.0f) / m1;
		y = 2.0f / 9.0f - (i + 1) * (4.0f / 9.0f) / m1;;
		singularities.push_back(Vector2f(x, y));

		x = -2.0f / 9.0f + (i + 0.5f) * (2.0f / 9.0f) / m1;
		y = 2.0f / 9.0f - (i + 0.5f) * (4.0f / 9.0f) / m1;
		colocations.push_back(Vector2f(x, y));
		normals.push_back(Vector2f(2.0f/sqrt(5), 1.0f/sqrt(5)).normalize());
	}

	for (int i = 0; i < m2; i++)
	{
		float x, y;

		x = 0.0f + (i + 1) * (2.0f / 9.0f) / m2;
		y = -2.0f / 9.0f + (i + 1) * (4.0f / 9.0f) / m2;
		singularities.push_back(Vector2f(x, y));

		x = 0.0f + (i + 0.5f) * (2.0f / 9.0f) / m2;
		y = -2.0f / 9.0f + (i + 0.5f) * (4.0f / 9.0f) / m2;
		colocations.push_back(Vector2f(x, y));
		normals.push_back(Vector2f(2.0f / sqrt(5), -1.0f / sqrt(5)).normalize());
		//normals.push_back(Vector2f(-1.0f, 0.0f).normalize()); 
	}


	tailSources.clear();
	tailSourceGammaNum.clear();
	tailSources.push_back(Vector2f(-2.0f / 9.0f, 2.0f / 9.0f));
	tailSourceGammaNum.push_back(0);
	tailSources.push_back(Vector2f(0.0f, -2.0f / 9.0f));
	tailSourceGammaNum.push_back(M/2);
	tailSources.push_back(Vector2f(2.0f / 9.0f, 2.0f / 9.0f));
	tailSourceGammaNum.push_back(M - 1);

	tails.clear();
	tails.resize(tailSources.size());

	tailGamma.clear();
	tailGamma.resize(tailSources.size());

	stage = 1;
	return;
}

float FlowSolver::u_j(int j, float x, float y) const
{
	if (stage < 1)
		throw "Curve is not splitted";
	if (j >= M)
		throw "Point index (j) out of range";

	float x_j = singularities[j].x;
	float y_j = singularities[j].y;
	
	float dist_sq = (x - x_j) * (x - x_j) + (y - y_j) * (y - y_j);
	dist_sq = std::max(dist_sq, delta * delta);

	float res = 1.0f / (2.0f * Pi) * (y_j - y) / dist_sq;

	return res;
}

float FlowSolver::u_p(float x, float y, float x_p, float y_p) const
{
	if (stage < 1)
		throw "Curve is not splitted";

	float x_j = x_p;
	float y_j = y_p;

	float dist_sq = (x - x_j) * (x - x_j) + (y - y_j) * (y - y_j);
	dist_sq = std::max(dist_sq, delta * delta);

	float res = 1.0f / (2.0f * Pi) * (y_j - y) / dist_sq;

	return res;
}

float FlowSolver::v_j(int j, float x, float y) const
{
	if (stage < 1)
		throw "Curve is not splitted";
	if (j >= M)
		throw "Point index (j) out of range";

	float x_j = singularities[j].x;
	float y_j = singularities[j].y;
	
	float dist_sq = (x - x_j) * (x - x_j) + (y - y_j) * (y - y_j);
	dist_sq = std::max(dist_sq, delta * delta);

	float res = 1.0f / (2.0f * Pi) * (x - x_j) / dist_sq;

	return res;
}

float FlowSolver::v_p(float x, float y, float x_p, float y_p) const
{
	if (stage < 1)
		throw "Curve is not splitted";

	float x_j = x_p;
	float y_j = y_p;

	float dist_sq = (x - x_j) * (x - x_j) + (y - y_j) * (y - y_j);
	dist_sq = std::max(dist_sq, delta * delta);

	float res = 1.0f / (2.0f * Pi) * (x - x_j) / dist_sq;

	return res;
}

Vector2f FlowSolver::V_j(int j, float x, float y) const
{
	if (stage < 1)
		throw "Curve is not splitted";
	if (j >= M)
		throw "Point index (j) out of range";

	return Vector2f(u_j(j, x, y), v_j(j, x, y));
}

Vector2f FlowSolver::V_p(float x, float y, float x_p, float y_p) const
{
	if (stage < 1)
		throw "Curve is not splitted";

	return Vector2f(u_p(x, y, x_p, y_p), v_p(x, y, x_p, y_p));
}

float FlowSolver::phi_j(int j, float x, float y) const
{
	if (stage < 1)
		throw "Curve is not splitted";
	if (j >= M)
		throw "Point index (j) out of range";

	float x_j = singularities[j].x;
	float y_j = singularities[j].y;
	float dist = sqrt((x - x_j) * (x - x_j) + (y - y_j) * (y - y_j));

	if (dist < delta)
		return 0.0f;

	return std::atan2(y - y_j, x - x_j);
}

float FlowSolver::phi_p(float x, float y, float x_p, float y_p) const
{
	if (stage < 1)
		throw "Curve is not splitted";

	float x_j = x_p;
	float y_j = y_p;
	float dist = sqrt((x - x_j) * (x - x_j) + (y - y_j) * (y - y_j));

	if (dist < delta)
		return 0.0f;

	return std::atan2(y - y_j, x - x_j);
}

inline float FlowSolver::u_Inf() const
{
	return cos(alpha);
}

inline float FlowSolver::v_Inf() const
{
	return sin(alpha);
}

inline Vector2f FlowSolver::V_Inf() const
{
	return Vector2f(cos(alpha), sin(alpha));
}

Vector2f FlowSolver::V(float x, float y) const
{
	if (stage < 3)
		throw "System of equations is not solved";
	
	Vector2f res = V_Inf();
	for (int j = 0; j < M; j++)
		res = res + V_j(j, x, y) * gamma[j];

	for (int tail = 0; tail < int(tails.size()); tail++)
		for (int i = 0; i < int(tails[tail].size()); i++)
			res = res + V_p(x, y, tails[tail][i].x, tails[tail][i].y) * tailGamma[tail][i];

	return res;
}

float FlowSolver::phiVihr(float x, float y) const
{
	if (stage < 3)
		throw "System of equations is not solved";

	float res = x * u_Inf() + y * v_Inf();
	for (int j = 0; j < M; j++)
		res += gamma[j] * phi_j(j, x, y) / 2.0f / Pi;

	return res;
}

float FlowSolver::phi(float x, float y) const
{
	/* TODO: modify this function to avoid stripped output*/
	if (stage < 3)
		throw "System of equations is not solved";

		
	float res = x * u_Inf() + y * v_Inf();
	for (int j = 0; j < M; j++)
		res += gamma[j] * phi_j(j, x, y) / 2.0f / Pi;

	for (int tail = 0; tail < int(tails.size()); tail++)
		for (int i = 0; i < int(tails[tail].size()); i++)
			res = res + tailGamma[tail][i] *
				phi_p(x, y, tails[tail][i].x, tails[tail][i].y) / 2.0f / Pi;
	
	return res;

}

float FlowSolver::psi(float x, float y) const
{
	if (stage < 3)
		throw "System of equations is not solved";
	
	float res = exp(y * u_Inf() - x * v_Inf());;
	for (int j = 0; j < M; j++)
	{
		float x_j = singularities[j].x;
		float y_j = singularities[j].y;
	
		float dist = sqrt((x - x_j) * (x - x_j) + (y - y_j) * (y - y_j));
		dist = std::max(dist, delta);

		res *= pow(dist, - gamma[j] / (2.0f * Pi));
	}

	for (int tail = 0; tail < int(tails.size()); tail++)
		for (int i = 0; i < int(tails[tail].size()); i++)
		{
			float x_j = tails[tail][i].x;
			float y_j = tails[tail][i].y;

			float dist = sqrt((x - x_j) * (x - x_j) + (y - y_j) * (y - y_j));
			dist = std::max(dist, delta);

			res *= pow(dist, -tailGamma[tail][i] / (2.0f * Pi));
		}

	res = log(res);

	return res;
}

float FlowSolver::phiContinuous(float x, float y) const
{
	float res = 0.0f;

	float gamma_sum1 = 0.0f, gamma_sum2 = 0.0f;
	compf tmp_res(0.0f, 0.0f);
	compf z(x, y);
	compf z01, z02, z_star, A;

	// starting
	gamma_sum1 = 0.0f;

	if (int(tails[0].size()) > 0)
	{
		// bottom tail
		for (int p = 0; p < int(tails[0].size()) - 1; p++)
		{
			gamma_sum1 += tailGamma[0][p];
			z01 = compf(tails[0][p].x, tails[0][p].y);
			z02 = compf(tails[0][p + 1].x, tails[0][p + 1].y);
			z_star = (z01 + z02) / 2.0f;
			A = gamma_sum1 * (z02 - z01);
			if (abs(z - z_star) > delta)
				tmp_res += A / (2.0f * Pi * compf(0.0, 1.0) * (z - z_star));
		}

		// bottom tail link
		gamma_sum1 += tailGamma[0][int(tails[0].size()) - 1];
		z01 = compf(tails[0][int(tails[0].size()) - 1].x, tails[0][int(tails[0].size()) - 1].y);
		z02 = compf(singularities[0].x, singularities[0].y);
		z_star = (z01 + z02) / 2.0f;
		A = gamma_sum1 * (z02 - z01);
		if (abs(z - z_star) > delta)
			tmp_res += A / (2.0f * Pi * compf(0.0, 1.0) * (z - z_star));
	}

	// going from the other end now
	gamma_sum2 = 0.0f;

	if ((tails[1].size()) > 0)
	{
		// the top tail
		for (int p = 0; p < int(tails[1].size()) - 1; p++)
		{
			gamma_sum2 += tailGamma[1][p];
			z01 = compf(tails[1][p].x, tails[1][p].y);
			z02 = compf(tails[1][p + 1].x, tails[1][p + 1].y);
			z_star = (z01 + z02) / 2.0f;
			A = gamma_sum2 * (z02 - z01);
			if (abs(z - z_star) > delta)
				tmp_res += A / (2.0f * Pi * compf(0.0, 1.0) * (z - z_star));
		}

		// top tail link
		gamma_sum2 += tailGamma[1][int(tails[1].size()) - 1];
		z01 = compf(tails[1][int(tails[1].size()) - 1].x, tails[1][int(tails[1].size()) - 1].y);
		z02 = compf(singularities[M - 1].x, singularities[M - 1].y);
		z_star = (z01 + z02) / 2.0f;
		A = gamma_sum2 * (z02 - z01);
		if (abs(z - z_star) > delta)
			tmp_res += A / (2.0f * Pi * compf(0.0, 1.0) * (z - z_star));
	}

	// the obstacle, reversed traversal
	for (int j = M - 1; j > 0; j--)
	{
		gamma_sum2 += gamma[j];
		z01 = compf(singularities[j].x, singularities[j].y);
		z02 = compf(singularities[j - 1].x, singularities[j - 1].y);
		z_star = (z01 + z02) / 2.0f;
		A = gamma_sum2 * (z02 - z01);
		if (abs(z - z_star) > delta)
			tmp_res += A / (2.0f * Pi * compf(0.0, 1.0) * (z - z_star));
	}

	res = tmp_res.real();
	res += x * u_Inf() + y * v_Inf();

	// // we don't have to add anything else here
	// // since there is zero circulation
	// res += gamma[0] * phi_j(0, x, y) / 2.0f / Pi;

	return res;
}

float FlowSolver::Cp(float x, float y) const
{
	if (stage < 3)
		throw "System of equations is not solved";

	Vector2f v1 = V(x, y);
	Vector2f v2 = V_Inf();

	float v1_sq = v1.x * v1.x + v1.y * v1.y;
	float v2_sq = v2.x * v2.x + v2.y * v2.y;
	float res = 1.0f - v1_sq / v2_sq;

	return res;
}

void FlowSolver::setGamma0(float value)
{
	gamma0 = value;
	stage = std::min(stage, 1);
}

void FlowSolver::setAlpha(float value)
{
	if (!(- Pi - 0.0001 < value && value < Pi + 0.0001))
		throw "New value of alpha must be from (-Pi/2; Pi/2)";

	alpha = value;
	stage = std::min(stage, 1);
}

void FlowSolver::setM(int value)
{
	if (!(4 <= value && value <= 100))
		throw "New value of M must be from [4; 100]";

	M = value - value % 4 + 1;
	stage = std::min(stage, 0);
}

void FlowSolver::makeGammaSystem()
{
	if (stage < 1)
		throw "Curve is not splitted";

	A.clear();
	A.resize(M, std::vector<float>(M + 1));
	
	for (int k = 0; k < M - 1; k++)
	{
		for (int j = 0; j < M; j++)
		{
			A[k][j] = V_j(j, colocations[k].x, colocations[k].y) * normals[k];
		}
		A[k][M] = - (V_Inf() * normals[k]);

		for (int tail = 0; tail < int(tails.size()); tail++)
			for (int i = 0; i < int(tails[tail].size()); i++)
				A[k][M] -= tailGamma[tail][i] *
				(V_p(colocations[k].x, colocations[k].y,
				tails[tail][i].x, tails[tail][i].y) * normals[k]);

	}

	for (int j = 0; j < M; j++)
	{
		A[M - 1][j] = 1.0f;
	}
	//A[M - 1][M] = gamma0;
	A[M - 1][M] = 0.0f;
	for (int tail = 0; tail < int(tails.size()); tail++)
		for (int i = 0; i < int(tails[tail].size()); i++)
			A[M - 1][M] -= tailGamma[tail][i];

	stage = 2;
	return;
}

void FlowSolver::solveGamma()
{
	if (stage < 2)
		throw "System of equations is not created";

	gamma.clear();

	for (int k = 0; k < M; k++)
	{
		double mx = -1;
		int v = -1;
		for (int l = k; l < M; l++)
		{
			if (abs(A[l][k]) > mx) {
				mx = abs(A[l][k]);
				v = l;
			}
		}
		if (v != k) swap(A[v], A[k]);

		if (abs(A[k][k]) <= Eps)
			throw "System of equations is singular";

		for (int j = k + 1; j <= M; j++)
			A[k][j] /= A[k][k];
		A[k][k] = 1.0f;

		for (int l = 0; l < M; l++)
			if (l != k)
			{
				for (int j = k + 1; j <= M; j++)
					A[l][j] -= A[k][j] * A[l][k];
				A[l][k] = 0.0f;
			}
	}

	for (int k = 0; k < M; k++)
		gamma.push_back(A[k][M]);

	stage = 3;
	return;
}

void FlowSolver::init()
{
	splitCurve();
	makeGammaSystem();
	solveGamma();
}

void FlowSolver::nextStep()
{
	if (stage < 3)
		throw "Stage is not initialized. Use init() to initialize";

	float W_max = 0.0f;

	std::vector<Vector2f> controlPoints;
	controlPoints.push_back(Vector2f(0.0f, -0.5f));
	controlPoints.push_back(Vector2f(0.0f, 0.5f));
	controlPoints.push_back(Vector2f(0.5f, 0.5f));

	for (int i = 0; i < int(controlPoints.size()); ++i)
		W_max = std::max(W_max, V(controlPoints[i].x, controlPoints[i].y).length());

	float timeDelta = 0.5f * big_delta / W_max;

	for (int tail = 0; tail < int(tails.size()); tail++)
	{
		tails[tail].push_back(tailSources[tail]);
		tailGamma[tail].push_back(gamma[tailSourceGammaNum[tail]]);

		for (int i = 0; i < int(tails[tail].size()); i++)
		{
			// the points are not allowed to cross the obstacle
			Vector2f s = tails[tail][i];
			Vector2f step = V(s.x, s.y) * timeDelta;
			Vector2f f = s + step;
			float shift = 0.f;

			if (i == tails[tail].size() - 1)
			{
				if (tail == 3)
					f.y += 3.0f * step.y;
				tails[tail][i] = f;
				continue;
			}
			
			// checking top horizontal segment
			//if the dot crossed the segment - reflecting its trajectory
			if (Vector2f::intersect(s, f, Vector2f(0.0f, 0.5f), Vector2f(0.5f, 0.5f)))
			{
				f.y += 2.0f * (0.5f - f.y);
			}
			// if the dot is withing the strip - shifting it out
			if (0.0f < f.x && f.x < 0.5f && abs(f.y - 0.5f) < big_delta)
			{
				if (f.y < 0.5f)
					f.y = 0.5f - big_delta;
				else
					f.y = 0.5f + big_delta;
			}

			// checking left vertical segment
			// if the dot crossed the segment - reflecting its trajectory
			if (Vector2f::intersect(s, f, Vector2f(0.0f, -0.5f), Vector2f(0.0f, 0.5f)))
			{
				f.x += 2.0f * (0.0f - f.x);
			}
			// if the dot is withing the strip - shifting it out
			if (-0.5f <= f.y && f.y <= 0.5f && abs(f.x - 0.0f) < big_delta)
			{
				if (f.x < 0.0f)
					f.x = 0.0f - big_delta;
				else
					f.x = 0.0f + big_delta;
			}

			tails[tail][i] = f;
		}
	}

	time += timeDelta;
	makeGammaSystem();
	solveGamma();
}

std::vector<Vector2f> FlowSolver::getColocations() const
{
	return colocations;
}

std::vector<Vector2f> FlowSolver::getSingularities() const
{
	return singularities;
}

std::vector< std::vector<Vector2f> > FlowSolver::getTails() const
{
	return tails;
}

std::vector<Vector2f> FlowSolver::getTailSources() const
{
	return tailSources;
}

std::vector< std::vector<float> > FlowSolver::getTailGamma() const
{
	return tailGamma;
}


float FlowSolver::getDelta() const
{
	return delta;
}

float FlowSolver::V_scal(float x, float y) const
{
	Vector2f v = V(x, y);
	return sqrt(v.x * v.x + v.y * v.y);
}

float FlowSolver::getTime() const
{
	return time;
}

void FlowSolver::restart()
{
	time = 0.0f;
	init();
}

float FlowSolver::getMaxGamma() const
{
	return maxGamma;
}