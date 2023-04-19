#pragma once
#include "../MeshDefinition.h"
#include <Eigen/Sparse>

typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::Triplet<double> T;

extern const double W_SMALL_PATCH;

struct DistortionMetric 
{
	double max_u_;
	double min_u_;
	double distortion_;
};

class Distortion
{
public:
	Distortion(const Mesh& mesh);

	void Compute(const std::vector<bool>& seam_v, 
		DistortionMetric& distortion);

private:
	void Init();

	void VertArea();
	void FaceAngle();
	void VertGauss();

	void MatrixC();

private:
	const Mesh& mesh_;

	Eigen::VectorXd f_angle;

	Eigen::VectorXd m_M_inv; // Areas of vertices
	Eigen::VectorXd m_K; // Curvatures of vertices

	SpMat m_L;
	SpMat m_C; // C = L*A^-1*L^T
	SpMat m_sub_C;
	SpMat P;
	SpMat PT;

	Eigen::SimplicialLDLT<SpMat> solver;
	Eigen::VectorXd m_b;

	Eigen::VectorXd chol_solution;
	Eigen::VectorXd solve_u;
};

