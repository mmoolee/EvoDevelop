#ifndef WEIGHTEDL1DR_HPP
#define WEIGHTEDL1DR_HPP

#include "DRSplitting.hpp"
#include "OMPHelper.h"

inline double lpNorm(VectorX x, double lp = 2)
{
	double rt = 0;
	for (int i = 0; i < x.size(); ++i)
	{
		rt += pow(fabs(x(i)), lp);
	}

	return pow(rt, 1.0 / lp);
}

/*
solve the following problem:
min ||wy||_1, s.t.y=Ax+b, ||x||_p<sigma
*/
template <bool useAA, bool wNorm> class WeightedL1DR : public DRSplitting<MatrixXX, useAA>
{
private:
	MatrixXX proxSubjectFun(const MatrixXX& t)
	{
		MatrixXX u = t;
		auto ux = u.block(0, 0, n, 1);
		auto uy = u.block(0, 1, n, 1);
		auto tx = t.block(0, 0, n, 1);
		auto ty = t.block(0, 1, n, 1);

		if (wNorm) ux = ATWAI.solve(tx + AT * (w.cwiseProduct(ty - b)));
		else ux = ATWAI.solve(tx + AT * (ty - b));
		uy = A * ux + b;

		return u;
	}

	MatrixXX proxObjectFun(const MatrixXX& s)
	{
		MatrixXX v = s;

		OMP_PARALLEL
		{
			OMP_SECTIONS
			{
				OMP_SECTION
				{
					auto vx = v.block(0, 0, n, 1);
					vx = sigma / fmax(sigma, lpNorm(vx, lp)) * vx;
				}

				OMP_SECTION
				{
					auto vy = v.block(0, 1, n, 1);
					if (wNorm)
						vy = (VectorX::Ones(n) - (vy / this->gamma).cwiseAbs().cwiseMax(1).cwiseInverse()).cwiseProduct(vy);
					else
						vy = (VectorX::Ones(n) -
							(vy.cwiseQuotient(w) / this->gamma).cwiseAbs().cwiseMax(1).cwiseInverse()).cwiseProduct(vy);
				}
			}
		}

		return v;
	}

	double norm(const MatrixXX& s)
	{
		auto x = s.block(0, 0, n, 1);
		auto y = s.block(0, 1, n, 1);
		if (wNorm) return (x.squaredNorm() + w.cwiseProduct(y.cwiseAbs2()).sum()) / n;
		else return (x.squaredNorm() + y.squaredNorm()) / n;
	}

	void updateGamma(int count, double res)
	{
		if (count % gammaN == 0) lres = res;
		if ((count + 1) % gammaN == 0 && res >= lres - 1e-3)
		{
			this->gamma /= 2;
			this->gamma = fmax(this->gamma, gammaL);
			this->resetAnderson(AA_trait<useAA>());
		}
	}

public:
	WeightedL1DR(const ColMajorSparseMatrix& A_, const VectorX& b_, double sigma_, int lp_)
		: A(A_), b(b_), n(b_.size()), sigma(sigma_), lp(lp_) {
		AT = A.transpose();
	}

	void init(double gamma_, double thres_, int maxIters_, int andersM_) {
		DRSplitting<MatrixXX, useAA>::init(gamma_, thres_, maxIters_, andersM_);
		gammaN = 2 * andersM_;
	}

	VectorX run(const VectorX& w_, const VectorX& x_)
	{
		w = w_;

		MatrixXX s(n, 2);
		s.block(0, 0, n, 1) = x_;
		s.block(0, 1, n, 1) = A * x_ + b;

		ColMajorSparseMatrix I(n, n);
		I.setIdentity();
		if (wNorm) ATWAI.compute(AT * w.asDiagonal() * A + I);
		else ATWAI.compute(AT * A + I);
		if (ATWAI.info() != Eigen::Success)
		{
			printf("WL1 LDLT factory failed\n");
			exit(EXIT_FAILURE);
		}

		s = DRSplitting<MatrixXX, useAA>::run(s);
		return s.block(0, 0, n, 1);
	}

private:
	int n, lp;
	ColMajorSparseMatrix A, AT;
	VectorX b;

	int gammaN = 20;
	double gammaL = 1e-3;

	double lres;

	VectorX w;
	Eigen::SimplicialLDLT<ColMajorSparseMatrix> ATWAI;

public:
	double sigma;
};

#endif // !WEIGHTEDL1DR_HPP
