#include "EnergyUpdate.h"
#include <Eigen/Dense>

const double P_SMALL_PATCH = 0.005;
const double MAX_DISTORTION = 1e4;

const double W_DISTORTION_EXP = 1000.0;

const double W_BOUND_LENGTH = 1.0;
const double W_DISTORTION = 1.0;
const double W_SEAM_SMOOTH = 1e3;
const double W_SMALL_PATCH = 100;
const double W_NARROW_PATCH = 10;

//// 47094
//const double W_SMALL_PATCH = 0;
//const double W_NARROW_PATCH = 0;

//// Flower, CAD
//const double W_SMALL_PATCH = 1;
//const double W_NARROW_PATCH = 1;

EnergyUpdate::EnergyUpdate(const Mesh& mesh,
	const Eigen::VectorXd& e_l,
	const Eigen::VectorXd& e_w,
	const Eigen::VectorXd& f_area,
	const std::vector<std::vector<int>>& f_f_adj,
	const std::vector<std::vector<int>>& v_f_adj)
	:mesh_(mesh), e_l_(e_l), e_w_(e_w), 
	f_f_adj_(f_f_adj), v_f_adj_(v_f_adj),
	f_area_(f_area), distortion_(mesh_)
{
	sum_e_l_ = e_l_.sum();

	mesh_area = MeshTools::Area(mesh);
}

void EnergyUpdate::SeamUpdate(const std::vector<bool>& seam,
	SeamEnergy& energy, double& fit_energy)
{
	UpdateSegId(seam);

	ComputeSmallPatch(seam, energy.seg_num_, energy.small_cont_);
	RateLength(seam, energy.rate_length_);
	SeamSmooth(seam, energy.seam_smooth_);
	SeamDistortion(seam, energy.dis_metric_);

	
	NarrowRegion(seam, energy.narrow_cont_);
	int narrow_energy = energy.narrow_cont_;
	//int narrow_energy = energy.narrow_cont_ > 0 ? 1 : 0;

	double dis_energy = exp(W_DISTORTION_EXP *
		(energy.dis_metric_.distortion_ / pow(DISTORTION_BOUND, 2) - 1.0));
	dis_energy = std::min(dis_energy, MAX_DISTORTION);

	fit_energy
		= W_BOUND_LENGTH * energy.rate_length_
		+ W_SEAM_SMOOTH * energy.seam_smooth_
		+ W_DISTORTION * dis_energy
		+ W_SMALL_PATCH * energy.small_cont_
		+ W_NARROW_PATCH * narrow_energy
		;
}

void EnergyUpdate::SetBound(double bound)
{
	DISTORTION_BOUND = bound;
}

void EnergyUpdate::ComputeSmallPatch(const std::vector<bool>& seam,
	int& seg_num, int& small_cont) const
{
	seg_num = seg_num_;

	std::vector<double> f_num(seg_num_, 0);
	for (int i = 0; i < mesh_.n_faces(); i++)
	{
		f_num[seg_id_[i]] += f_area_[i];
	}

	small_cont = 0;
	for (size_t i = 0; i < f_num.size(); i++)
	{
		if (f_num[i] < P_SMALL_PATCH * mesh_area)
		{
			small_cont++;
		}
	}

	//small_cont = *std::min_element(f_num.begin(), f_num.end());
	//small_cont = small_cont > P_SMALL_PATCH * mesh_area ? 0 : 1;
}

void EnergyUpdate::RateLength(const std::vector<bool>& seam,
	double& rate_length) const
{
	// Length
	double total_seam_length = 0;
	for (const EH& e_h : mesh_.edges())
	{
		if (mesh_.is_boundary(e_h)) continue;

		if (seam[e_h.idx()]) {
			//total_seam_length += e_l_[e_h.idx()];
			total_seam_length +=
				e_w_[e_h.idx()] * e_l_[e_h.idx()];
		}
	}

	rate_length = total_seam_length / sum_e_l_;
}

void EnergyUpdate::SeamSmooth(const std::vector<bool>& seam, 
	double& seam_energy) const
{
	std::vector<int> v_cont(mesh_.n_vertices(), 0);
	for (const EH& e_h : mesh_.edges())
	{
		if (!seam[e_h.idx()]) continue;

		const HEH& he_0 = mesh_.halfedge_handle(e_h, 0);
		const HEH& he_1 = mesh_.halfedge_handle(e_h, 1);

		v_cont[mesh_.to_vertex_handle(he_0).idx()]++;
		v_cont[mesh_.to_vertex_handle(he_1).idx()]++;
	}

	seam_energy = 0;
	for (const VH& v_h : mesh_.vertices())
	{
		if (v_cont[v_h.idx()] != 2) continue;

		OpenMesh::Vec3d temp_vec(0, 0, 0);
		for (const HEH& adj_he : mesh_.voh_range(v_h))
		{
			if (seam[mesh_.edge_handle(adj_he).idx()])
			{
				temp_vec += mesh_.calc_edge_vector(adj_he);
			}
		}

		seam_energy += temp_vec.sqrnorm();
	}

	seam_energy /= sum_e_l_;
}

void EnergyUpdate::SeamDistortion(const std::vector<bool>& seam,
	DistortionMetric& dis_energy)
{
	// Vertices on seams
	VH to_v, from_v;
	std::vector<bool> is_on_seam(mesh_.n_vertices(), false);
	for (const EH& e_h : mesh_.edges())
	{
		if (mesh_.is_boundary(e_h) || seam[e_h.idx()])
		{
			to_v = mesh_.to_vertex_handle(mesh_.halfedge_handle(e_h, 0));
			from_v = mesh_.from_vertex_handle(mesh_.halfedge_handle(e_h, 0));

			is_on_seam[to_v.idx()] = true;
			is_on_seam[from_v.idx()] = true;
		}
	}

	distortion_.Compute(is_on_seam, dis_energy);
}

void EnergyUpdate::NarrowRegion(const std::vector<bool>& seam, 
	int& cont)
{
	cont = 0;

	std::vector<int> v_cont(mesh_.n_vertices(), 0);
	for (const EH& e_h : mesh_.edges())
	{
		if (!seam[e_h.idx()]) continue;

		const HEH& he_0 = mesh_.halfedge_handle(e_h, 0);
		const HEH& he_1 = mesh_.halfedge_handle(e_h, 1);

		v_cont[mesh_.to_vertex_handle(he_0).idx()]++;
		v_cont[mesh_.to_vertex_handle(he_1).idx()]++;
	}

	std::vector<bool> f_status(mesh_.n_faces(), true);
	for (const VH& v_h : mesh_.vertices())
	{
		if (v_cont[v_h.idx()] > 2)
		{
			for (const int& adj_f_id : v_f_adj_[v_h.idx()])
			{
				f_status[adj_f_id] = false;
			}
		}
	}

	for (const FH& f_h : mesh_.faces())
	{
		if (!f_status[f_h.idx()]) continue;

		std::set<int> seg_set;
		for (const int& adj_f : f_f_adj_[f_h.idx()])
		{
			seg_set.insert(seg_id_[adj_f]);
		}

		if (seg_set.size() < 2)
		{
			f_status[f_h.idx()] = false;
		}
	}

	std::vector<std::vector<HEH>> he_path;
	std::vector<std::vector<HEH>>::iterator path_it;
	std::vector<HEH> temp_path;
	for (const FH& f_h : mesh_.faces())
	{
		if (!f_status[f_h.idx()]) continue;

		std::set<EH> temp_e_set;
		std::vector<bool> visited_f_status(mesh_.n_faces(), false);

		for (const int& adj_f : f_f_adj_[f_h.idx()])
		{
			visited_f_status[adj_f] = true;
		}

		for (const int& adj_f : f_f_adj_[f_h.idx()])
		{
			for (const EH& e_h : mesh_.fe_range(FH(adj_f)))
			{
				temp_e_set.insert(e_h);
			}
		}

		he_path.clear();
		for (const EH& e_h : temp_e_set)
		{
			if (!seam[e_h.idx()]) continue;

			HEH he_h = mesh_.halfedge_handle(e_h, 0);
			const FH& f_0 = mesh_.face_handle(he_h);
			const FH& f_1 = mesh_.opposite_face_handle(he_h);

			if (!f_0.is_valid() ||
				seg_id_[f_0.idx()] < seg_id_[f_1.idx()])
			{
				he_h = mesh_.opposite_halfedge_handle(he_h);
			}

			// Separate the moved boundary
			temp_path.clear();
			temp_path.emplace_back(he_h);

			path_it = he_path.begin();
			while (path_it != he_path.end())
			{
				// Is on the seam?
				// Can be merge with other pathes?
				if (v_cont[mesh_.to_vertex_handle(temp_path.back()).idx()] == 2
					&& mesh_.to_vertex_handle(temp_path.back())
					== mesh_.from_vertex_handle((*path_it).front()))
				{
					temp_path.insert(temp_path.end(), (*path_it).begin(), (*path_it).end());
					path_it = he_path.erase(path_it);
				}
				else if (v_cont[mesh_.from_vertex_handle(temp_path.front()).idx()] == 2
					&& mesh_.from_vertex_handle(temp_path.front())
					== mesh_.to_vertex_handle((*path_it).back()))
				{
					temp_path.insert(temp_path.begin(), (*path_it).begin(), (*path_it).end());
					path_it = he_path.erase(path_it);
				}
				else {
					path_it++;
				}
			}

			he_path.push_back(temp_path);
		}

		int inner_path_cont = 0;
		for (size_t i = 0; i < he_path.size(); i++)
		{
			bool all_bound = true;
			for (const HEH& he_h : he_path[i])
			{
				const FH& f_0 = mesh_.face_handle(he_h);
				const FH& f_1 = mesh_.opposite_face_handle(he_h);
				if (f_0.is_valid() && visited_f_status[f_0.idx()]
					&& f_1.is_valid() && visited_f_status[f_1.idx()])
				{
					all_bound = false;
					break;
				}
			}

			if (!all_bound)
			{
				inner_path_cont++;
			}
		}
		if (inner_path_cont > 1)
		{
			cont++;
		}
	}
}

void EnergyUpdate::UpdateSegId(const std::vector<bool>& seam)
{
	// Seg Num
	seg_id_.clear();
	seg_id_.resize(mesh_.n_faces(), -1);
	seg_num_ = 0;

	for (int i = 0; i < mesh_.n_faces(); i++)
	{
		if (seg_id_[i] != -1) continue;

		std::vector<int> f_stack;
		f_stack.push_back(i);
		do
		{
			int top_idx = f_stack.back(); f_stack.pop_back();
			seg_id_[top_idx] = seg_num_;

			for (const HEH& fh_h : mesh_.fh_range(FH(top_idx)))
			{
				FH oppo_f_h = mesh_.opposite_face_handle(fh_h);

				if (!oppo_f_h.is_valid()
					|| seg_id_[oppo_f_h.idx()] != -1) continue;

				int adj_e_idx = mesh_.edge_handle(fh_h).idx();

				if (!seam[adj_e_idx])
				{
					f_stack.push_back(oppo_f_h.idx());
				}
			}
		} while (f_stack.size() != 0);

		++seg_num_;
	}
}
