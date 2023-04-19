#include "ConesFlattening.h"
#include "Opt/WeightedL1DR.hpp"

//namespace ConesFlattening
//{
//	ColMajorSparseMatrix L, A, P;
//	VectorX K;
//	
//	int lp = 2;
//	double sigma;
//
//	void areaMat(const Mesh& mesh)
//	{
//		A.resize(mesh.n_vertices(), mesh.n_vertices());
//
//		double sumArea = 0;
//		std::vector<double> fArea(mesh.n_faces());
//
//		for (auto f : mesh.faces())
//		{
//			OpenMesh::Vec3d p[3];
//			auto cfv_it = mesh.cfv_begin(f);
//			p[0] = mesh.point(*cfv_it);		++cfv_it;
//			p[1] = mesh.point(*cfv_it);		++cfv_it;
//			p[2] = mesh.point(*cfv_it);
//
//			fArea[f.idx()] = ((p[0] - p[1]) % (p[1] - p[2])).norm();
//			sumArea += fArea[f.idx()];
//		}
//
//		std::vector<Eigen::Triplet<double>> trips;
//		trips.reserve(mesh.n_vertices());
//
//		for (auto v : mesh.vertices())
//		{
//			double vArea = 0;
//			for (auto vf : mesh.vf_range(v))
//			{
//				vArea += fArea[vf.idx()];
//			}
//
//			vArea /= (3 * sumArea);
//			trips.emplace_back(v.idx(), v.idx(), pow(vArea, 1.0 / lp));
//		}
//
//		A.setFromTriplets(trips.begin(), trips.end());
//	}
//
//	void interiorMat(const Mesh& mesh)
//	{
//		std::vector<Eigen::Triplet<double>> trips;
//		trips.reserve(mesh.n_vertices());
//
//		int rowId = 0;
//		for (auto v : mesh.vertices())
//		{
//			if (mesh.is_boundary(v)) continue;
//			trips.emplace_back(rowId, v.idx(), 1);
//			rowId++;
//		}
//
//		P.resize(rowId, mesh.n_vertices());
//		P.setFromTriplets(trips.begin(), trips.end());
//	}
//
//	void initYamabeCoef(const Mesh& mesh)
//	{
//		std::vector<OpenMesh::Vec3i> fv_idx(mesh.n_faces());
//		std::vector<OpenMesh::Vec3d> fv_theta(mesh.n_faces()), hf_fvCot(mesh.n_faces());
//
//		for (auto f : mesh.faces())
//		{
//			Mesh::VertexHandle v[3];
//			auto cfv_it = mesh.cfv_begin(f);
//			v[0] = *cfv_it;		++cfv_it;
//			v[1] = *cfv_it;		++cfv_it;
//			v[2] = *cfv_it;
//
//			OpenMesh::Vec3d e[3];
//			e[0] = (mesh.point(v[2]) - mesh.point(v[1])).normalized();
//			e[1] = (mesh.point(v[0]) - mesh.point(v[2])).normalized();
//			e[2] = (mesh.point(v[1]) - mesh.point(v[0])).normalized();
//
//			fv_idx[f.idx()][0] = v[0].idx();
//			fv_idx[f.idx()][1] = v[1].idx();
//			fv_idx[f.idx()][2] = v[2].idx();
//
//			fv_theta[f.idx()][0] = acos(-e[1] | e[2]);
//			fv_theta[f.idx()][1] = acos(-e[2] | e[0]);
//			fv_theta[f.idx()][2] = acos(-e[0] | e[1]);
//
//			hf_fvCot[f.idx()][0] = -0.5 * (e[1] | e[2]) / (e[1] % e[2]).norm();
//			hf_fvCot[f.idx()][1] = -0.5 * (e[2] | e[0]) / (e[2] % e[0]).norm();
//			hf_fvCot[f.idx()][2] = -0.5 * (e[0] | e[1]) / (e[0] % e[1]).norm();
//		}
//
//		K.resize(mesh.n_vertices());
//		for (auto v : mesh.vertices())
//		{
//			if (mesh.is_boundary(v))
//				K[v.idx()] = M_PI;
//			else
//				K[v.idx()] = 2 * M_PI;
//		}
//
//		for (auto f : mesh.faces())
//		{
//			K[fv_idx[f.idx()][0]] -= fv_theta[f.idx()][0];
//			K[fv_idx[f.idx()][1]] -= fv_theta[f.idx()][1];
//			K[fv_idx[f.idx()][2]] -= fv_theta[f.idx()][2];
//		}
//
//		std::vector<Eigen::Triplet<double>> trips;
//		trips.reserve(9 * mesh.n_faces());
//
//		for (auto f : mesh.faces())
//		{
//			trips.emplace_back(fv_idx[f.idx()][0], fv_idx[f.idx()][0], hf_fvCot[f.idx()][1] + hf_fvCot[f.idx()][2]);
//			trips.emplace_back(fv_idx[f.idx()][1], fv_idx[f.idx()][1], hf_fvCot[f.idx()][0] + hf_fvCot[f.idx()][2]);
//			trips.emplace_back(fv_idx[f.idx()][2], fv_idx[f.idx()][2], hf_fvCot[f.idx()][0] + hf_fvCot[f.idx()][1]);
//			trips.emplace_back(fv_idx[f.idx()][0], fv_idx[f.idx()][1], -hf_fvCot[f.idx()][2]);
//			trips.emplace_back(fv_idx[f.idx()][1], fv_idx[f.idx()][0], -hf_fvCot[f.idx()][2]);
//			trips.emplace_back(fv_idx[f.idx()][1], fv_idx[f.idx()][2], -hf_fvCot[f.idx()][0]);
//			trips.emplace_back(fv_idx[f.idx()][2], fv_idx[f.idx()][1], -hf_fvCot[f.idx()][0]);
//			trips.emplace_back(fv_idx[f.idx()][2], fv_idx[f.idx()][0], -hf_fvCot[f.idx()][1]);
//			trips.emplace_back(fv_idx[f.idx()][0], fv_idx[f.idx()][2], -hf_fvCot[f.idx()][1]);
//		}
//
//		L.resize(mesh.n_vertices(), mesh.n_vertices());
//		L.setFromTriplets(trips.begin(), trips.end());
//	}
//
//	void initCoef(const Mesh& mesh, double sigma_)
//	{
//		//printf("Initialization\n");
//		sigma = sigma_;
//
//		areaMat(mesh);
//		interiorMat(mesh);
//		initYamabeCoef(mesh);
//	}
//
//	void geneCone(VectorX& conesK)
//	{
//		ColMajorSparseMatrix LN = P * L * A.cwiseInverse() * P.transpose();
//		VectorX KN = P * K;
//		VectorX Aphi = VectorX::Zero(LN.cols());
//
//		double lambda = 1;
//		double thres = 1e-3;
//		double gamma = 1;
//
//		VectorX w;
//		WeightedL1DR<true, true> WL1(LN * sigma, KN, 1.0, lp);
//		VectorX conesK_last = KN;
//		conesK = KN;
//
//		//printf("Compute non-integer-constrained cones\n");
//
//		while (thres >= 1e-10)
//		{
//			WL1.init(gamma, thres, 2e4, 10);
//
//			w = (conesK.cwiseAbs() + lambda * VectorX::Ones(KN.size())).cwiseInverse();
//			//printf("------------Weighted L1------------");
//			Aphi = WL1.run(w, Aphi);
//			conesK = LN * Aphi * sigma + KN;
//				
//			if (thres > 1e-6 && (conesK - conesK_last).norm() < conesK.size()*1e-8) break;
//			conesK_last = conesK;
//
//			lambda /= 2;
//			thres /= 2;
//		}
//
//		conesK = P.transpose() * (LN * Aphi * sigma + KN);
//	}
//}

ConesFlattening::ConesFlattening(const Mesh& m)
	:mesh(m)
{
}

void ConesFlattening::initCoef(double sigma_)
{
	//printf("Initialization\n");
	sigma = sigma_;

	areaMat();
	interiorMat();
	initYamabeCoef();
}

void ConesFlattening::geneCone(VectorX& conesK)
{
	ColMajorSparseMatrix LN = P * L * A.cwiseInverse() * P.transpose();
	VectorX KN = P * K;
	VectorX Aphi = VectorX::Zero(LN.cols());

	double lambda = 1;
	double thres = 1e-3;
	double gamma = 1;

	VectorX w;
	WeightedL1DR<true, true> WL1(LN * sigma, KN, 1.0, lp);
	VectorX conesK_last = KN;
	conesK = KN;

	//printf("Compute non-integer-constrained cones\n");

	while (thres >= 1e-10)
	{
		WL1.init(gamma, thres, 2000, 10);

		w = (conesK.cwiseAbs() + lambda * VectorX::Ones(KN.size())).cwiseInverse();
		//printf("------------Weighted L1------------");
		Aphi = WL1.run(w, Aphi);
		conesK = LN * Aphi * sigma + KN;

		if (thres > 1e-6 && (conesK - conesK_last).norm() < conesK.size() * 1e-8) break;
		conesK_last = conesK;

		lambda /= 2;
		thres /= 2;
	}

	conesK = P.transpose() * (LN * Aphi * sigma + KN);
}

void ConesFlattening::areaMat()
{
	A.resize(mesh.n_vertices(), mesh.n_vertices());

	double sumArea = 0;
	std::vector<double> fArea(mesh.n_faces());

	for (auto f : mesh.faces())
	{
		OpenMesh::Vec3d p[3];
		auto cfv_it = mesh.cfv_begin(f);
		p[0] = mesh.point(*cfv_it);		++cfv_it;
		p[1] = mesh.point(*cfv_it);		++cfv_it;
		p[2] = mesh.point(*cfv_it);

		fArea[f.idx()] = ((p[0] - p[1]) % (p[1] - p[2])).norm();
		sumArea += fArea[f.idx()];
	}

	std::vector<Eigen::Triplet<double>> trips;
	trips.reserve(mesh.n_vertices());

	for (auto v : mesh.vertices())
	{
		double vArea = 0;
		for (auto vf : mesh.vf_range(v))
		{
			vArea += fArea[vf.idx()];
		}

		vArea /= (3 * sumArea);
		trips.emplace_back(v.idx(), v.idx(), pow(vArea, 1.0 / lp));
	}

	A.setFromTriplets(trips.begin(), trips.end());
}

void ConesFlattening::initYamabeCoef()
{
	std::vector<OpenMesh::Vec3i> fv_idx(mesh.n_faces());
	std::vector<OpenMesh::Vec3d> fv_theta(mesh.n_faces()), hf_fvCot(mesh.n_faces());

	for (auto f : mesh.faces())
	{
		Mesh::VertexHandle v[3];
		auto cfv_it = mesh.cfv_begin(f);
		v[0] = *cfv_it;		++cfv_it;
		v[1] = *cfv_it;		++cfv_it;
		v[2] = *cfv_it;

		OpenMesh::Vec3d e[3];
		e[0] = (mesh.point(v[2]) - mesh.point(v[1])).normalized();
		e[1] = (mesh.point(v[0]) - mesh.point(v[2])).normalized();
		e[2] = (mesh.point(v[1]) - mesh.point(v[0])).normalized();

		fv_idx[f.idx()][0] = v[0].idx();
		fv_idx[f.idx()][1] = v[1].idx();
		fv_idx[f.idx()][2] = v[2].idx();

		fv_theta[f.idx()][0] = acos(-e[1] | e[2]);
		fv_theta[f.idx()][1] = acos(-e[2] | e[0]);
		fv_theta[f.idx()][2] = acos(-e[0] | e[1]);

		hf_fvCot[f.idx()][0] = -0.5 * (e[1] | e[2]) / (e[1] % e[2]).norm();
		hf_fvCot[f.idx()][1] = -0.5 * (e[2] | e[0]) / (e[2] % e[0]).norm();
		hf_fvCot[f.idx()][2] = -0.5 * (e[0] | e[1]) / (e[0] % e[1]).norm();
	}

	K.resize(mesh.n_vertices());
	for (auto v : mesh.vertices())
	{
		if (mesh.is_boundary(v))
			K[v.idx()] = M_PI;
		else
			K[v.idx()] = 2 * M_PI;
	}

	for (auto f : mesh.faces())
	{
		K[fv_idx[f.idx()][0]] -= fv_theta[f.idx()][0];
		K[fv_idx[f.idx()][1]] -= fv_theta[f.idx()][1];
		K[fv_idx[f.idx()][2]] -= fv_theta[f.idx()][2];
	}

	std::vector<Eigen::Triplet<double>> trips;
	trips.reserve(9 * mesh.n_faces());

	for (auto f : mesh.faces())
	{
		trips.emplace_back(fv_idx[f.idx()][0], fv_idx[f.idx()][0], hf_fvCot[f.idx()][1] + hf_fvCot[f.idx()][2]);
		trips.emplace_back(fv_idx[f.idx()][1], fv_idx[f.idx()][1], hf_fvCot[f.idx()][0] + hf_fvCot[f.idx()][2]);
		trips.emplace_back(fv_idx[f.idx()][2], fv_idx[f.idx()][2], hf_fvCot[f.idx()][0] + hf_fvCot[f.idx()][1]);
		trips.emplace_back(fv_idx[f.idx()][0], fv_idx[f.idx()][1], -hf_fvCot[f.idx()][2]);
		trips.emplace_back(fv_idx[f.idx()][1], fv_idx[f.idx()][0], -hf_fvCot[f.idx()][2]);
		trips.emplace_back(fv_idx[f.idx()][1], fv_idx[f.idx()][2], -hf_fvCot[f.idx()][0]);
		trips.emplace_back(fv_idx[f.idx()][2], fv_idx[f.idx()][1], -hf_fvCot[f.idx()][0]);
		trips.emplace_back(fv_idx[f.idx()][2], fv_idx[f.idx()][0], -hf_fvCot[f.idx()][1]);
		trips.emplace_back(fv_idx[f.idx()][0], fv_idx[f.idx()][2], -hf_fvCot[f.idx()][1]);
	}

	L.resize(mesh.n_vertices(), mesh.n_vertices());
	L.setFromTriplets(trips.begin(), trips.end());
}

void ConesFlattening::interiorMat()
{
	std::vector<Eigen::Triplet<double>> trips;
	trips.reserve(mesh.n_vertices());

	int rowId = 0;
	for (auto v : mesh.vertices())
	{
		if (mesh.is_boundary(v)) continue;
		trips.emplace_back(rowId, v.idx(), 1);
		rowId++;
	}

	P.resize(rowId, mesh.n_vertices());
	P.setFromTriplets(trips.begin(), trips.end());
}
