#include "Distortion.h"

Distortion::Distortion(const Mesh& mesh)
	:mesh_(mesh)
{
	Init();
}

void Distortion::Compute(const std::vector<bool>& seam_v,
	DistortionMetric& distortion)
{
	int n_v = mesh_.n_vertices();

	//std::vector<T> row_trip;
	//std::vector<T> col_trip;
	//row_trip.reserve(n_v);
	//col_trip.reserve(n_v);
	//
	//int row_id = 0;
	//for (int i = 0; i < n_v; i++)
	//{
	//	if (!seam_v[i])
	//	{
	//		row_trip.emplace_back(row_id, i, 1);
	//		col_trip.emplace_back(i, row_id, 1);
	//		row_id++;
	//	}
	//}

	////std::cout << "inner size : " << row_id << std::endl;
	//P.resize(row_id, n_v);
	//PT.resize(n_v, row_id);

	//P.setFromTriplets(row_trip.begin(), row_trip.end());
	//PT.setFromTriplets(col_trip.begin(), col_trip.end());

	//m_sub_C = P * m_C * PT;

	std::vector<T> row_trip;
	std::vector<T> col_trip;
	row_trip.reserve(n_v);
	col_trip.reserve(n_v);
	
	int row_id = 0;
	std::vector<int> glo_loc(n_v, -1);
	for (int i = 0; i < n_v; i++)
	{
		if (!seam_v[i])
		{
			row_trip.emplace_back(row_id, i, 1);
			col_trip.emplace_back(i, row_id, 1);
			glo_loc[i] = row_id;
			row_id++;
		}
	}

	//std::cout << "inner size : " << row_id << std::endl;
	P.resize(row_id, n_v);
	PT.resize(n_v, row_id);

	P.setFromTriplets(row_trip.begin(), row_trip.end());
	PT.setFromTriplets(col_trip.begin(), col_trip.end());

	std::vector<T> C_trip;
	for (int k = 0; k < m_C.outerSize(); ++k)
	{
		for (SpMat::InnerIterator it(m_C, k); it; ++it)
		{
			if (glo_loc[it.row()] == -1
				|| glo_loc[it.col()] == -1) continue;
			
			C_trip.emplace_back(glo_loc[it.row()],
				glo_loc[it.col()],
				it.value());
		}
	}
	m_sub_C.resize(row_id, row_id);
	m_sub_C.setFromTriplets(C_trip.begin(), C_trip.end());




	solver.compute(m_sub_C);

	m_b = P * m_K;

	chol_solution = solver.solve(m_b);

	solve_u = m_M_inv.cwiseProduct(m_L.transpose() * PT * chol_solution);

	distortion.max_u_ = solve_u.maxCoeff();
	distortion.min_u_ = solve_u.minCoeff();

	distortion.distortion_ = m_b.transpose() * chol_solution;
}

void Distortion::Init()
{
	VertArea();

	FaceAngle();
	VertGauss();

	MatrixC();
}

void Distortion::VertArea()
{
	Eigen::VectorXd m_M;
	m_M.setConstant(mesh_.n_vertices(), 0);
	double fA;
	for (const FH& f_h : mesh_.faces())
	{
		fA = mesh_.calc_face_area(f_h) / 3;
		for (const VH& fv_h : mesh_.fv_range(f_h))
		{
			m_M[fv_h.idx()] += fA;
		}
	}

	m_M = m_M / m_M.sum();

	m_M_inv = m_M.cwiseInverse();
}

void Distortion::FaceAngle()
{
	f_angle.setConstant(mesh_.n_faces() * 3, 0);
	for (const FH& f_h : mesh_.faces())
	{
		std::vector<OpenMesh::Vec3d> he_v(3);
		int cout = 0;
		for (const HEH& he_h : mesh_.fh_range(f_h))
		{
			he_v[cout] = mesh_.calc_edge_vector(he_h).normalized();
			cout++;
		}

		double cos_0 = -dot(he_v[1], he_v[2]);
		double cos_1 = -dot(he_v[2], he_v[0]);
		double cos_2 = -dot(he_v[0], he_v[1]);

		f_angle[3 * f_h.idx() + 0] = acos(cos_0);
		f_angle[3 * f_h.idx() + 1] = acos(cos_1);
		f_angle[3 * f_h.idx() + 2] = acos(cos_2);
	}
}

void Distortion::VertGauss()
{
	m_K.setConstant(mesh_.n_vertices(), 2 * M_PI);

	for (const VH& v_h : mesh_.vertices())
	{
		if (mesh_.is_boundary(v_h)) {
			m_K[v_h.idx()] = M_PI;
		}
	}

	for (const FH& f_h : mesh_.faces())
	{
		std::vector<int> v_ids(3);
		int cout = 0;
		for (const HEH& he_h : mesh_.fh_range(f_h))
		{
			v_ids[cout] = mesh_.to_vertex_handle(mesh_.next_halfedge_handle(he_h)).idx();
			cout++;
		}

		m_K[v_ids[0]] -= f_angle[3 * f_h.idx() + 0];
		m_K[v_ids[1]] -= f_angle[3 * f_h.idx() + 1];
		m_K[v_ids[2]] -= f_angle[3 * f_h.idx() + 2];
	}
}


void Distortion::MatrixC()
{
	double cot_weight;
	VH to_v, from_v;
	std::vector<T> trip;
	trip.reserve(12 * mesh_.n_faces());
	for (const FH& f_h : mesh_.faces())
	{
		int cout = 0;
		for (const HEH& fh_h : mesh_.fh_range(f_h))
		{
			to_v = mesh_.to_vertex_handle(fh_h);
			from_v = mesh_.from_vertex_handle(fh_h);

			cot_weight = 0.5 / tan(f_angle[3 * f_h.idx() + cout]);
			cout++;

			trip.emplace_back(to_v.idx(), from_v.idx(), -cot_weight);
			trip.emplace_back(from_v.idx(), to_v.idx(), -cot_weight);
			trip.emplace_back(from_v.idx(), from_v.idx(), cot_weight);
			trip.emplace_back(to_v.idx(), to_v.idx(), cot_weight);
		}
	}

	m_L.resize(mesh_.n_vertices(), mesh_.n_vertices());
	m_L.setFromTriplets(trip.begin(), trip.end());

	m_C = m_L * m_M_inv.asDiagonal() * m_L.transpose();
}
