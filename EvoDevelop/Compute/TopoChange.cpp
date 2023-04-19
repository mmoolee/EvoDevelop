#include "TopoChange.h"

const int N_RANDOM = 10000;
const double LENGTH_RATE = 0.05;
const double SMALL_PERCENT = 0.005;
const double ANCHOR_COS = cos(M_PI / 180 * 150);

TopoChange::TopoChange(const Mesh& mesh)
	:SegMesh(mesh)
{
}

TopoChange::TopoChange(const Mesh& mesh,
	const std::vector<bool>& seam,
	const std::vector<bool>& cone)
	:SegMesh(mesh, seam), cone_status_(cone)
{
}



bool TopoChange::RandomSeamRemove(const Eigen::VectorXd& e_l,
	const Eigen::VectorXd& e_w,
	const Eigen::VectorXd& v_K,
	bool len_fixed,
	bool cone_fixed)
{
	// Cut seam
	std::vector<std::vector<HEH>> he_path;
	InnerSeamSegment(he_path);

	//for (size_t i = 0; i < he_path.size(); i++)
	//{
	//	std::ofstream out_f("out_f_" + std::to_string(i) + ".txt");
	//	for (const HEH& he_h:he_path[i])
	//	{
	//		out_f << mesh_.face_handle(he_h).idx() << std::endl;
	//	}
	//	out_f.close();
	//}

	if (he_path.size() == 0) return false;

	std::vector<int> edge_seam_id(mesh_.n_edges(), -1); 	// Boundary index for each edge
	std::map<std::pair<int, int>, std::vector<int>> patch_id_to_seam_id;
	EdgeAndPatchToSeam(he_path, edge_seam_id, patch_id_to_seam_id);

	std::vector<bool> picked_boundaries(he_path.size(), true); // Archive of boundaries

	// ============= Remove the cases with more than two adjacent boundaries ===============
	std::map<std::pair<int, int>, std::vector<int>>::iterator
		iter = patch_id_to_seam_id.begin();
	while (iter != patch_id_to_seam_id.end())
	{
		if (iter->second.size() != 1) {
			for (int path_id : iter->second) {
				picked_boundaries[path_id] = false;
			}
		}

		iter++;
	}

	std::vector<double> boundaries_length(he_path.size(), 0);
	std::vector<double> w_boundaries_length(he_path.size(), 0);
	std::vector<double> patches_length(seg_num_, 0);

	// Find each boundary
	for (const EH& e_h : mesh_.edges())
	{
		if (!seam_status_[e_h.idx()]) continue;

		const HEH& he_h = mesh_.halfedge_handle(e_h, 0);
		const FH& f_0 = mesh_.face_handle(he_h);
		const FH& f_1 = mesh_.opposite_face_handle(he_h);

		if (f_0.is_valid())
		{
			patches_length[seg_id_[f_0.idx()]] += e_l[e_h.idx()];
		}

		if (f_1.is_valid())
		{
			patches_length[seg_id_[f_1.idx()]] += e_l[e_h.idx()];
		}

		if (edge_seam_id[e_h.idx()] == -1) continue;

		boundaries_length[edge_seam_id[e_h.idx()]] += e_l[e_h.idx()];
		w_boundaries_length[edge_seam_id[e_h.idx()]] +=
			e_l[e_h.idx()] * e_w[e_h.idx()];
	}

	for (int path_id = 0; path_id < boundaries_length.size(); path_id++)
	{
		if (!picked_boundaries[path_id]) continue;

		double path_length = boundaries_length[path_id];

		// Indices of patches and their boundaries
		const HEH he_h = he_path[path_id].front();
		const FH f_0 = mesh_.face_handle(he_h);
		const FH f_1 = mesh_.opposite_face_handle(he_h);

		int patch_id_0 = seg_id_[f_0.idx()];
		int patch_id_1 = seg_id_[f_1.idx()];

		if (len_fixed)
		{
			if (1.0 * path_length < LENGTH_RATE * patches_length[patch_id_0]
				|| 1.0 * path_length < LENGTH_RATE * patches_length[patch_id_1])
			{
				picked_boundaries[path_id] = false;
			}
		}
		else
		{
			if (1.0 * path_length < LENGTH_RATE * patches_length[patch_id_0]
				&& 1.0 * path_length < LENGTH_RATE * patches_length[patch_id_1])
			{
				picked_boundaries[path_id] = false;
			}
		}
	}

	int remove_id = -1;
	double sum_length = 0;
	std::vector<int> picked_ids;
	std::vector<double> percent_len;
	if (cone_fixed)
	{
		for (int path_id = 0; path_id < he_path.size(); path_id++)
		{
			if (!picked_boundaries[path_id]) continue;

			for (int i = 1; i < he_path[path_id].size(); i++)
			{
				const VH from_v = mesh_.from_vertex_handle(he_path[path_id][i]);

				if (cone_status_[from_v.idx()])
				{
					picked_boundaries[path_id] = false;
					break;
				}
			}
		}

		for (size_t i = 0; i < picked_boundaries.size(); i++)
		{
			if (picked_boundaries[i])
			{
				const HEH he_h = he_path[i].front();
				double temp_percent = w_boundaries_length[i];

				sum_length += temp_percent;
				percent_len.emplace_back(temp_percent);

				picked_ids.push_back(i);
			}
		}
	}
	else
	{
		//// Numbers of cones on boundaries
		//std::vector<int> cone_num(he_path.size(), 0);
		//for (size_t i = 0; i < he_path.size(); i++)
		//{
		//	for (size_t j = 1; j < he_path[i].size(); j++)
		//	{
		//		const VH& from_v = mesh_.from_vertex_handle(he_path[i][j]);

		//		if (cone_status_[from_v.idx()])
		//		{
		//			cone_num[i]++;
		//		}
		//	}
		//}

		//int min_num = 1;
		////int min_num = DBL_MAX;
		////for (size_t i = 0; i < cone_num.size(); i++)
		////{
		////	if (cone_num[i] == 0) continue;
		////	
		////	min_num = std::min(min_num, cone_num[i]);
		////}


		for (size_t i = 0; i < he_path.size(); i++)
		{
			//if (picked_boundaries[i] && cone_num[i] == min_num)
			if (picked_boundaries[i])
			{
				//picked_ids.emplace_back(i);

				//// Prabability
				//double temp_len = 0;
				//for (size_t j = 1; j < he_path[i].size(); j++)
				//{
				//	const VH& from_v = mesh_.from_vertex_handle(he_path[i][j]);

				//	temp_len += v_K[from_v.idx()];
				//}

				//temp_len /= he_path[i].size();

				//percent_len.emplace_back(temp_len);
				//sum_length += temp_len;


				//double temp_percent = w_boundaries_length[i];
				double temp_percent = 1.0;

				sum_length += temp_percent;
				percent_len.emplace_back(temp_percent);

				picked_ids.push_back(i);
			}
		}
	}

	if (picked_ids.size() == 0) return false;

	std::vector<double> total_len(picked_ids.size() + 1, 0);
	for (int i = 0; i < picked_ids.size(); i++)
	{
		total_len[i + 1] = total_len[i]
			+ percent_len[i] / sum_length;
	}

	double random_prab = 1.0 * (rand() % N_RANDOM) / N_RANDOM;

	for (int i = 0; i < picked_ids.size(); i++)
	{
		if (random_prab >= total_len[i] && random_prab < total_len[i + 1])
		{
			remove_id = picked_ids[i];
			break;
		}
	}

	if (remove_id == -1) return false;

	std::vector<bool> ori_seam = seam_status_;

	for (const HEH& he_h : he_path[remove_id])
	{
		seam_status_[mesh_.edge_handle(he_h).idx()] = false;
	}

	SeamUpdate();
	IdxUpdate();

	if (HighGenus())
	{
		seam_status_ = ori_seam;
		SeamUpdate();
		IdxUpdate();
		return false;
	}

	if (!cone_fixed)
	{
		// Cone
		for (size_t j = 1; j < he_path[remove_id].size(); j++)
		{
			const VH& from_v = mesh_.from_vertex_handle(he_path[remove_id][j]);
			cone_status_[from_v.idx()] = false;
		}
	}

	return true;
}

bool TopoChange::RandomSmallPatch(const Eigen::VectorXd& e_l,
	const Eigen::VectorXd& e_w,
	bool cone_fixed)
{
	// Find the small patchs' Id
	std::vector<double> patch_area(seg_num_, 0);
	for (int i = 0; i < mesh_.n_faces(); i++)
	{
		patch_area[seg_id_[i]] += mesh_.calc_face_area(FH(i));
	}

	double mesh_area = MeshTools::Area(mesh_);
	std::vector<int> small_id;
	for (size_t i = 0; i < patch_area.size(); i++)
	{
		if (patch_area[i] < SMALL_PERCENT * mesh_area)
		{
			small_id.push_back(i);
		}
	}

	if (small_id.size() == 0) return false;

	std::vector<bool> patch_status(seg_num_, false);
	for (int p_id : small_id)
	{
		patch_status[p_id] = true;
	}

	// Cut seam
	std::vector<std::vector<HEH>> he_path;
	InnerSeamSegment(he_path);

	if (he_path.size() == 0) return false;

	//for (size_t i = 0; i < he_path.size(); i++)
	//{
	//	std::ofstream out_f("out_f_" + std::to_string(i) + ".txt");
	//	for (const HEH& he_h : he_path[i])
	//	{
	//		out_f << mesh_.face_handle(he_h).idx() << std::endl;
	//	}
	//	out_f.close();
	//}

	std::vector<int> edge_seam_id(mesh_.n_edges(), -1); 	// Boundary index for each edge
	std::map<std::pair<int, int>, std::vector<int>> patch_id_to_seam_id;
	EdgeAndPatchToSeam(he_path, edge_seam_id, patch_id_to_seam_id);

	std::vector<bool> picked_boundaries(he_path.size(), true); // Archive of boundaries

	// ============= Remove the cases with more than two adjacent boundaries ===============
	std::map<std::pair<int, int>, std::vector<int>>::iterator
		iter = patch_id_to_seam_id.begin();
	while (iter != patch_id_to_seam_id.end())
	{
		if (iter->second.size() != 1) {
			for (int path_id : iter->second) {
				picked_boundaries[path_id] = false;
			}
		}

		iter++;
	}

	if (cone_fixed)
	{
		for (int path_id = 0; path_id < he_path.size(); path_id++)
		{
			if (!picked_boundaries[path_id]) continue;

			for (int i = 1; i < he_path[path_id].size(); i++)
			{
				const VH& from_v = mesh_.from_vertex_handle(he_path[path_id][i]);

				if (cone_status_[from_v.idx()])
				{
					picked_boundaries[path_id] = false;
					break;
				}
			}
		}
	}

	for (int path_id = 0; path_id < he_path.size(); path_id++)
	{
		if (!picked_boundaries[path_id]) continue;

		for (int i = 1; i < he_path[path_id].size(); i++)
		{
			const VH& from_v = mesh_.from_vertex_handle(he_path[path_id][i]);

			if (cone_status_[from_v.idx()])
			{
				picked_boundaries[path_id] = false;
				break;
			}
		}
	}

	std::vector<double> boundaries_length(he_path.size(), 0);
	std::vector<double> w_boundaries_length(he_path.size(), 0);
	std::vector<double> patches_length(seg_num_, 0);

	// Find each boundary
	for (const EH& e_h : mesh_.edges())
	{
		if (!seam_status_[e_h.idx()]) continue;

		const HEH& he_h = mesh_.halfedge_handle(e_h, 0);
		const FH& f_0 = mesh_.face_handle(he_h);
		const FH& f_1 = mesh_.opposite_face_handle(he_h);

		if (f_0.is_valid())
		{
			patches_length[seg_id_[f_0.idx()]] += e_l[e_h.idx()];
		}

		if (f_1.is_valid())
		{
			patches_length[seg_id_[f_1.idx()]] += e_l[e_h.idx()];
		}
		
		if (edge_seam_id[e_h.idx()] == -1) continue;

		boundaries_length[edge_seam_id[e_h.idx()]] += e_l[e_h.idx()];
		w_boundaries_length[edge_seam_id[e_h.idx()]] +=
			e_l[e_h.idx()] * e_w[e_h.idx()];
	}

	for (int path_id = 0; path_id < boundaries_length.size(); path_id++)
	{
		if (!picked_boundaries[path_id]) continue;

		double path_length = boundaries_length[path_id];

		// Indices of patches and their boundaries
		const HEH& he_h = he_path[path_id].front();
		const FH& f_0 = mesh_.face_handle(he_h);
		const FH& f_1 = mesh_.opposite_face_handle(he_h);

		int patch_id_0 = seg_id_[f_0.idx()];
		int patch_id_1 = seg_id_[f_1.idx()];

		if (1.0 * path_length < LENGTH_RATE * patches_length[patch_id_0]
			&& 1.0 * path_length < LENGTH_RATE * patches_length[patch_id_1])
		{
			picked_boundaries[path_id] = false;
		}
	}

	// ============= Small patches ===============
	for (int path_id = 0; path_id < boundaries_length.size(); path_id++)
	{
		if (!picked_boundaries[path_id]) continue;

		// Indices of patches and their boundaries
		const HEH& he_h = he_path[path_id].front();
		const FH& f_0 = mesh_.face_handle(he_h);
		const FH& f_1 = mesh_.opposite_face_handle(he_h);

		int patch_id_0 = seg_id_[f_0.idx()];
		int patch_id_1 = seg_id_[f_1.idx()];

		if (!patch_status[patch_id_0] && !patch_status[patch_id_1])
		{
			picked_boundaries[path_id] = false;
		}
	}

	std::vector<int> picked_ids;
	for (size_t i = 0; i < picked_boundaries.size(); i++)
	{
		if (picked_boundaries[i])
		{
			picked_ids.push_back(i);
		}
	}

	if (picked_ids.size() == 0) return false;

	double sum_length = 0;
	std::vector<double> percent_len;
	for (int path_id : picked_ids)
	{
		const HEH& he_h = he_path[path_id].front();
		const FH& f_0 = mesh_.face_handle(he_h);
		const FH& f_1 = mesh_.opposite_face_handle(he_h);

		//double max_len = std::max(patches_length[seg_id_[f_0.idx()]],
		//	patches_length[seg_id_[f_1.idx()]]);

		//double temp_percent = w_boundaries_length[path_id] / max_len;
		double temp_percent = w_boundaries_length[path_id];

		sum_length += temp_percent;
		percent_len.emplace_back(temp_percent);
	}

	std::vector<double> total_len(picked_ids.size() + 1, 0);
	for (int i = 0; i < picked_ids.size(); i++)
	{
		total_len[i + 1] = total_len[i]
			+ percent_len[i] / sum_length;
	}

	double random_prab = 1.0 * (rand() % N_RANDOM) / N_RANDOM;

	int remove_id = -1;
	for (int i = 0; i < picked_ids.size(); i++)
	{
		if (random_prab >= total_len[i] && random_prab < total_len[i + 1])
		{
			remove_id = picked_ids[i];
			break;
		}
	}

	std::vector<bool> ori_seam = seam_status_;

	for (const HEH& he_h : he_path[remove_id])
	{
		seam_status_[mesh_.edge_handle(he_h).idx()] = false;
	}

	SeamUpdate();
	IdxUpdate();

	if (HighGenus())
	{
		seam_status_ = ori_seam;
		SeamUpdate();
		IdxUpdate();
		return false;
	}

	if (!cone_fixed)
	{
		// Cone
		for (size_t j = 1; j < he_path[remove_id].size(); j++)
		{
			const VH& from_v = mesh_.from_vertex_handle(he_path[remove_id][j]);
			cone_status_[from_v.idx()] = false;
		}
	}

	return true;
}

bool TopoChange::RandomNarrowRemove(const Eigen::VectorXd& e_l, 
	const Eigen::VectorXd& e_w,
	const std::vector<std::vector<int>>& f_f_adj,
	const std::vector<std::vector<int>>& v_f_adj,
	bool cone_fixed)
{
	std::vector<int> narrow_set;

	std::vector<int> v_cont(mesh_.n_vertices(), 0);
	VertexAdjSeam(v_cont);

	std::vector<bool> f_status(mesh_.n_faces(), true);
	for (const VH& v_h : mesh_.vertices())
	{
		if (v_cont[v_h.idx()] > 2)
		{
			for (const int& adj_f_id : v_f_adj[v_h.idx()])
			{
				f_status[adj_f_id] = false;
			}
		}
		//else if (cone_fixed && cone_status_[v_h.idx()])
		//{
		//	for (const int& adj_f_id : v_f_adj[v_h.idx()])
		//	{
		//		f_status[adj_f_id] = false;
		//	}
		//}
	}

	for (const FH& f_h:mesh_.faces())
	{
		if (!f_status[f_h.idx()]) continue;
		
		std::set<int> seg_set;
		for (const int& adj_f : f_f_adj[f_h.idx()])
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

		for (const int& adj_f : f_f_adj[f_h.idx()])
		{
			visited_f_status[adj_f] = true;
		}

		for (const int& adj_f : f_f_adj[f_h.idx()])
		{
			for (const EH& e_h:mesh_.fe_range(FH(adj_f)))
			{
				temp_e_set.insert(e_h);
			}
		}

		he_path.clear();
		for (const EH& e_h: temp_e_set)
		{
			if (!seam_status_[e_h.idx()]) continue;

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
			for (const HEH& he_h:he_path[i])
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
			//for (size_t i = 0; i < he_path.size(); i++)
			//{
			//	std::ofstream file_out("f_" + std::to_string(f_h.idx()) 
			//		+ "_" + std::to_string(i) + ".txt");
			//	for (const HEH& temp_he : he_path[i])
			//	{
			//		file_out << mesh_.face_handle(temp_he).idx() << std::endl;
			//	}
			//	file_out.close();
			//}
			//
			narrow_set.emplace_back(f_h.idx());
		}


		//const std::vector<int>& temp_adj_f = f_f_adj[f_h.idx()];

		//std::vector<int> temp_seg(mesh_.n_faces(), -2);

		//for (const int& temp_f_id : temp_adj_f)
		//{
		//	temp_seg[temp_f_id] = -1;
		//}

		//int temp_seg_num = 0;

		//for (const int& temp_f_id : temp_adj_f)
		//{
		//	if (temp_seg[temp_f_id] != -1) continue;

		//	std::vector<int> f_stack;
		//	f_stack.push_back(temp_f_id);
		//	do
		//	{
		//		int top_idx = f_stack.back(); f_stack.pop_back();
		//		temp_seg[top_idx] = temp_seg_num;

		//		for (const HEH& fh_h : mesh_.fh_range(FH(top_idx)))
		//		{
		//			FH oppo_f_h = mesh_.opposite_face_handle(fh_h);

		//			if (oppo_f_h.is_valid()
		//				&& temp_seg[oppo_f_h.idx()] == -1)
		//			{
		//				int adj_e_idx = mesh_.edge_handle(fh_h).idx();

		//				if (!seam_status_[adj_e_idx])
		//				{
		//					f_stack.push_back(oppo_f_h.idx());
		//				}
		//			}
		//		}
		//	} while (f_stack.size() != 0);

		//	++temp_seg_num;
		//}

		//if (temp_seg_num > 2)
		//{
		//	narrow_set.emplace_back(f_h.idx());
		//}
	}

	if (narrow_set.size() == 0) return false;
	
	//std::ofstream file_out("narrow_f.txt");
	//for (int f_id:narrow_set)
	//{
	//	file_out << f_id << std::endl;
	//}
	//file_out.close();

	int pick_f = narrow_set[rand() % narrow_set.size()];

	//std::cout <<"pick: " << pick_f << std::endl;

	int cur_seg_id = seg_id_[pick_f];

	f_status.clear();
	f_status.resize(mesh_.n_faces(), false);
	for (const int& adj_f : f_f_adj[pick_f])
	{
		f_status[adj_f] = true;
	}

	he_path.clear();
	for (const int& adj_f:f_f_adj[pick_f])
	{
		if (seg_id_[adj_f] != cur_seg_id) continue;
		
		for (const HEH& tri_he:mesh_.fh_range(FH(adj_f)))
		{
			if (seam_status_[mesh_.edge_handle(tri_he).idx()]) continue;

			const FH& f_0 = mesh_.face_handle(tri_he);
			const FH& f_1 = mesh_.opposite_face_handle(tri_he);

			if (!f_0.is_valid() || !f_status[f_0.idx()]) continue;

			if (!f_1.is_valid() || !f_status[f_1.idx()])
			{
				//std::cout << f_0.idx() << " " << f_1.idx() << std::endl;

				// Separate the moved boundary
				temp_path.clear();
				temp_path.emplace_back(tri_he);

				path_it = he_path.begin();
				while (path_it != he_path.end())
				{
					// Is on the seam?
					// Can be merge with other pathes?
					if (v_cont[mesh_.to_vertex_handle(temp_path.back()).idx()] == 0
						&& mesh_.to_vertex_handle(temp_path.back())
						== mesh_.from_vertex_handle((*path_it).front()))
					{
						temp_path.insert(temp_path.end(), (*path_it).begin(), (*path_it).end());
						path_it = he_path.erase(path_it);
					}
					else if (v_cont[mesh_.from_vertex_handle(temp_path.front()).idx()] == 0
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
		}
	}

	if (he_path.size() == 0) return false;
	
	//for (int i = 0; i < he_path.size(); ++i)
	//{
	//	std::ofstream file_path("narrow_f_" + std::to_string(i) + ".txt");
	//	for (const HEH& he_h : he_path[i])
	//	{
	//		file_path << mesh_.face_handle(he_h) << std::endl;
	//	}
	//	file_path.close();
	//}

	//// ======================== Extended region =============================== //
	//std::vector<std::set<int>> narrow_region_set;
	//for (const auto& single_path:he_path)
	//{
	//	FH start_f = mesh_.face_handle(single_path.front());

	//	std::vector<bool> temp_seam = seam_status_;
	//	for (const HEH& he_h:single_path)
	//	{
	//		temp_seam[mesh_.edge_handle(he_h).idx()] = true;
	//	}

	//	std::vector<bool> visited_f(mesh_.n_faces(), false);
	//	std::vector<int> f_stack;
	//	std::set<int> temp_region;

	//	f_stack.push_back(start_f.idx());
	//	temp_region.insert(start_f.idx());
	//	do
	//	{
	//		int top_idx = f_stack.back(); f_stack.pop_back();
	//		visited_f[top_idx] = true;

	//		for (const HEH& fh_h : mesh_.fh_range(FH(top_idx)))
	//		{
	//			FH oppo_f_h = mesh_.opposite_face_handle(fh_h);

	//			if (oppo_f_h.is_valid()
	//				&& !visited_f[oppo_f_h.idx()])
	//			{
	//				int adj_e_idx = mesh_.edge_handle(fh_h).idx();

	//				if (!temp_seam[adj_e_idx])
	//				{
	//					f_stack.push_back(oppo_f_h.idx());
	//					temp_region.insert(oppo_f_h.idx());
	//				}
	//			}
	//		}
	//	} while (f_stack.size() != 0);

	//	narrow_region_set.emplace_back(temp_region);		
	//}

	////for (int i = 0; i < narrow_region_set.size(); ++i)
	////{
	////	std::ofstream region_out("narrow_region" + std::to_string(i) +".txt");
	////	for (const int& f_id:narrow_region_set[i])
	////	{
	////		region_out << f_id << std::endl;
	////	}
	////	region_out.close();
	////}
	//
	//std::set<int> min_region;
	//std::vector<HEH> new_he_set;
	//int min_size = INT_MAX;
	//for (int i = 0; i < narrow_region_set.size(); ++i)
	//{
	//	if (min_size > narrow_region_set[i].size())
	//	{
	//		min_region = narrow_region_set[i];
	//		min_size = min_region.size();
	//		new_he_set = he_path[i];
	//	}
	//}

	////std::ofstream min_region_out("min_region.txt");
	////for (const int& f_id : min_region)
	////{
	////	min_region_out << f_id << std::endl;
	////}
	////min_region_out.close();

	std::set<int> min_region;
	for (const int& adj_f:f_f_adj[pick_f])
	{
		if (seg_id_[adj_f] == cur_seg_id)
		{
			min_region.insert(adj_f);
		}
	}

	// ==========================================================================

	if (min_region.size() == 0) return false;

	// Cut seam
	he_path.clear();
	for (const int& f_id:min_region)
	{
		for (const EH& e_h:mesh_.fe_range(FH(f_id)))
		{
			if (!seam_status_[e_h.idx()]) continue;
			if (mesh_.is_boundary(e_h)) continue;

			HEH he_h = mesh_.halfedge_handle(e_h, 0);

			const FH& f_0 = mesh_.face_handle(he_h);
			const FH& f_1 = mesh_.opposite_face_handle(he_h);

			// Pick one halfedge of each edge
			if (seg_id_[f_0.idx()] < seg_id_[f_1.idx()])
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
	}

	if (he_path.size() == 0) return false;

	//for (int i = 0; i < he_path.size(); ++i)
	//{
	//	std::ofstream file_path("path_" + std::to_string(i) + ".txt");
	//	for (const HEH& he_h : he_path[i])
	//	{
	//		file_path << mesh_.face_handle(he_h) << std::endl;
	//	}
	//	file_path.close();
	//}

	std::vector<int> edge_seam_id(mesh_.n_edges(), -1); // Boundary index for each edge
	std::map<std::pair<int, int>, std::vector<int>> patch_id_to_seam_id;
	EdgeAndPatchToSeam(he_path, edge_seam_id, patch_id_to_seam_id);

	std::vector<bool> picked_boundaries(he_path.size(), true); // Archive of boundaries

	// ============= Remove the cases with more than two adjacent boundaries ===============
	std::map<std::pair<int, int>, std::vector<int>>::iterator
		iter = patch_id_to_seam_id.begin();
	while (iter != patch_id_to_seam_id.end())
	{
		if (iter->second.size() != 1) {
			for (int path_id : iter->second) {
				picked_boundaries[path_id] = false;
			}
		}

		iter++;
	}

	std::vector<double> boundaries_length(he_path.size(), 0);
	std::vector<double> w_boundaries_length(he_path.size(), 0);
	std::vector<double> patches_length(seg_num_, 0);

	// Find each boundary
	for (const EH& e_h : mesh_.edges())
	{
		if (!seam_status_[e_h.idx()]) continue;

		const HEH& he_h = mesh_.halfedge_handle(e_h, 0);
		const FH& f_0 = mesh_.face_handle(he_h);
		const FH& f_1 = mesh_.opposite_face_handle(he_h);

		if (f_0.is_valid())
		{
			patches_length[seg_id_[f_0.idx()]] += e_l[e_h.idx()];
		}

		if (f_1.is_valid())
		{
			patches_length[seg_id_[f_1.idx()]] += e_l[e_h.idx()];
		}

		if (edge_seam_id[e_h.idx()] == -1) continue;

		boundaries_length[edge_seam_id[e_h.idx()]] += e_l[e_h.idx()];
		w_boundaries_length[edge_seam_id[e_h.idx()]] +=
			e_l[e_h.idx()] * e_w[e_h.idx()];
	}

	//for (int path_id = 0; path_id < boundaries_length.size(); path_id++)
	//{
	//	if (!picked_boundaries[path_id]) continue;

	//	double path_length = boundaries_length[path_id];

	//	// Indices of patches and their boundaries
	//	const HEH he_h = he_path[path_id].front();
	//	const FH f_0 = mesh_.face_handle(he_h);
	//	const FH f_1 = mesh_.opposite_face_handle(he_h);

	//	int patch_id_0 = seg_id_[f_0.idx()];
	//	int patch_id_1 = seg_id_[f_1.idx()];

	//	if (1.0 * path_length < LENGTH_RATE * patches_length[patch_id_0]
	//			&& 1.0 * path_length < LENGTH_RATE * patches_length[patch_id_1])
	//	{
	//		picked_boundaries[path_id] = false;
	//	}
	//}

	int remove_id = -1;
	double sum_length = 0;
	std::vector<int> picked_ids;
	std::vector<double> percent_len;
	if (cone_fixed)
	{
		for (int path_id = 0; path_id < he_path.size(); path_id++)
		{
			if (!picked_boundaries[path_id]) continue;

			for (int i = 1; i < he_path[path_id].size(); i++)
			{
				const VH from_v = mesh_.from_vertex_handle(he_path[path_id][i]);

				if (cone_status_[from_v.idx()])
				{
					picked_boundaries[path_id] = false;
					break;
				}
			}
		}
	}
	
	for (size_t i = 0; i < picked_boundaries.size(); i++)
	{
		if (picked_boundaries[i])
		{
			double temp_percent = w_boundaries_length[i];

			sum_length += temp_percent;
			percent_len.emplace_back(temp_percent);

			picked_ids.push_back(i);
		}
	}

	if (picked_ids.size() == 0) return false;

	std::vector<double> total_len(picked_ids.size() + 1, 0);
	for (int i = 0; i < picked_ids.size(); i++)
	{
		total_len[i + 1] = total_len[i]
			+ percent_len[i] / sum_length;
	}

	double random_prab = 1.0 * (rand() % N_RANDOM) / N_RANDOM;

	for (int i = 0; i < picked_ids.size(); i++)
	{
		if (random_prab >= total_len[i] && random_prab < total_len[i + 1])
		{
			remove_id = picked_ids[i];
			break;
		}
	}

	if (remove_id == -1) return false;

	std::vector<int> ori_seg_id = seg_id_;
	
	int new_seg_id = seg_id_[mesh_.face_handle(he_path[remove_id].front()).idx()];
	if (new_seg_id == cur_seg_id)
	{
		new_seg_id = seg_id_[mesh_.opposite_face_handle(he_path[remove_id].front()).idx()];
	}

	for (const int& f_id: min_region)
	{
		seg_id_[f_id] = new_seg_id;
	}

	IdxUpdate();
	IdxUpdate();

	if (HighGenus())
	{
		seg_id_ = ori_seg_id;
		IdxUpdate();
		IdxUpdate();
		return false;
	}

	//std::vector<bool> ori_seam = seam_status_;

	//for (const HEH& he_h : he_path[remove_id])
	//{
	//	seam_status_[mesh_.edge_handle(he_h).idx()] = false;
	//}

	//for (const HEH& he_h : new_he_set)
	//{
	//	seam_status_[mesh_.edge_handle(he_h).idx()] = true;
	//}

	//SeamUpdate();
	//IdxUpdate();

	//if (HighGenus())
	//{
	//	seam_status_ = ori_seam;
	//	SeamUpdate();
	//	IdxUpdate();
	//	return false;
	//}

	if (!cone_fixed)
	{
		// Cone
		for (size_t j = 1; j < he_path[remove_id].size(); j++)
		{
			const VH& from_v = mesh_.from_vertex_handle(he_path[remove_id][j]);
			cone_status_[from_v.idx()] = false;
		}
	}

	return true;
}

void TopoChange::VertexAdjSeam(std::vector<int>& v_cont)
{
	v_cont.assign(mesh_.n_vertices(), 0);
	for (const EH& e_h : mesh_.edges())
	{
		if (!seam_status_[e_h.idx()]) continue;

		const HEH& he_0 = mesh_.halfedge_handle(e_h, 0);
		const HEH& he_1 = mesh_.halfedge_handle(e_h, 1);

		v_cont[mesh_.to_vertex_handle(he_0).idx()]++;
		v_cont[mesh_.to_vertex_handle(he_1).idx()]++;
	}
}

void TopoChange::InnerSeamSegment(std::vector<std::vector<HEH>>& he_path)
{
	he_path.clear();

	std::vector<int> v_cont;
	VertexAdjSeam(v_cont);

	// Cut the seams
	std::vector<std::vector<HEH>>::iterator path_it;
	std::vector<HEH> temp_path;
	for (const EH& e_h : mesh_.edges())
	{
		if (!seam_status_[e_h.idx()]) continue;
		if (mesh_.is_boundary(e_h)) continue;

		HEH he_h = mesh_.halfedge_handle(e_h, 0);

		const FH& f_0 = mesh_.face_handle(he_h);
		const FH& f_1 = mesh_.opposite_face_handle(he_h);

		// Pick one halfedge of each edge
		if (seg_id_[f_0.idx()] < seg_id_[f_1.idx()])
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
}

bool TopoChange::RandomOneTriangleChange(bool cone_fixed)
{
	//std::vector<int> v_cont;
	//VertexAdjSeam(v_cont);

	//std::vector<HEH> he_stack;
	//for (const HEH& he_h : mesh_.halfedges())
	//{
	//	if (mesh_.is_boundary(he_h)) continue;

	//	const HEH& next_he = mesh_.next_halfedge_handle(he_h);
	//	const HEH& prev_he = mesh_.prev_halfedge_handle(he_h);

	//	const EH& next_e = mesh_.edge_handle(next_he);
	//	const EH& prev_e = mesh_.edge_handle(prev_he);

	//	const FH& oppo_next_f = mesh_.opposite_face_handle(next_he);
	//	const FH& oppo_prev_f = mesh_.opposite_face_handle(prev_he);

	//	if (seam_status_[prev_e.idx()]
	//		&& seam_status_[next_e.idx()])
	//	{
	//		if (!oppo_prev_f.is_valid()
	//			&& !oppo_next_f.is_valid()) continue;

	//		const VH& to_v = mesh_.to_vertex_handle(next_he);

	//		if (cone_fixed && 
	//			v_cont[to_v.idx()] == 2 && cone_status_[to_v.idx()]) continue;
	//		

	//		int adj_seg_0 = oppo_prev_f.is_valid() ?
	//			seg_id_[oppo_prev_f.idx()] : -1;

	//		int adj_seg_1 = oppo_next_f.is_valid() ?
	//			seg_id_[oppo_next_f.idx()] : -1;

	//		if (v_cont[to_v.idx()] != 2 && adj_seg_0 == adj_seg_1) continue;

	//		he_stack.push_back(he_h);

	//		if (v_cont[to_v.idx()] == 2)
	//		{
	//			he_stack.push_back(he_h);
	//		}
	//	}
	//}

	//if (he_stack.size() == 0) return false;

	//HEH pick_he = he_stack[rand() % he_stack.size()];

	////std::cout << "Tri change: " << pick_f.idx() << std::endl;

	//FH oppo_cur_f = mesh_.opposite_face_handle(pick_he);
	//FH oppo_next_f = mesh_.opposite_face_handle(mesh_.next_halfedge_handle(pick_he));
	//FH oppo_prev_f = mesh_.opposite_face_handle(mesh_.prev_halfedge_handle(pick_he));
	//std::vector<int> adj_seg_id;
	//if (oppo_next_f.is_valid() && 
	//	(v_cont[mesh_.to_vertex_handle(pick_he).idx()] == 2 
	//		|| seg_id_[oppo_cur_f.idx()] != seg_id_[oppo_next_f.idx()]))
	//{
	//	adj_seg_id.emplace_back(seg_id_[oppo_next_f.idx()]);
	//}

	//if (oppo_prev_f.is_valid() && 
	//	(v_cont[mesh_.from_vertex_handle(pick_he).idx()] == 2
	//		|| seg_id_[oppo_cur_f.idx()] != seg_id_[oppo_prev_f.idx()]))
	//{
	//	adj_seg_id.emplace_back(seg_id_[oppo_prev_f.idx()]);
	//}

	//if (adj_seg_id.size() == 0) return false;

	//int pick_seg_id = adj_seg_id[rand() % adj_seg_id.size()];

	//seg_id_[mesh_.face_handle(pick_he).idx()] = pick_seg_id;

	//if (!cone_fixed && 
	//	cone_status_[mesh_.opposite_vh(pick_he).idx()]) {
	//	
	//	cone_status_[mesh_.opposite_vh(pick_he).idx()] = false;

	//	if (seg_id_[mesh_.face_handle(pick_he).idx()] 
	//		== seg_id_[oppo_next_f.idx()])
	//	{
	//		cone_status_[mesh_.from_vertex_handle(pick_he).idx()] = true;
	//	}
	//	else
	//	{
	//		cone_status_[mesh_.to_vertex_handle(pick_he).idx()] = true;
	//	}
	//}

	//IdxUpdate();

	//return true;

	std::vector<int> ori_seg = seg_id_;

	std::vector<int> v_cont;
	VertexAdjSeam(v_cont);

	std::vector<HEH> he_stack;
	for (const HEH& he_h : mesh_.halfedges())
	{
		if (mesh_.is_boundary(he_h)) continue;

		const HEH& next_he = mesh_.next_halfedge_handle(he_h);
		const HEH& prev_he = mesh_.prev_halfedge_handle(he_h);

		const EH& next_e = mesh_.edge_handle(next_he);
		const EH& prev_e = mesh_.edge_handle(prev_he);

		const FH& oppo_next_f = mesh_.opposite_face_handle(next_he);
		const FH& oppo_prev_f = mesh_.opposite_face_handle(prev_he);

		if (seam_status_[prev_e.idx()]
			&& seam_status_[next_e.idx()])
		{
			if (!oppo_prev_f.is_valid()
				&& !oppo_next_f.is_valid()) continue;

			const VH& to_v = mesh_.to_vertex_handle(next_he);

			if (cone_fixed &&
				v_cont[to_v.idx()] == 2 && cone_status_[to_v.idx()]) continue;


			int adj_seg_0 = oppo_prev_f.is_valid() ?
				seg_id_[oppo_prev_f.idx()] : -1;

			int adj_seg_1 = oppo_next_f.is_valid() ?
				seg_id_[oppo_next_f.idx()] : -1;

			if (v_cont[to_v.idx()] != 2 && adj_seg_0 == adj_seg_1) continue;

			he_stack.push_back(he_h);

			if (v_cont[to_v.idx()] == 2)
			{
				he_stack.push_back(he_h);
			}
		}
	}

	if (he_stack.size() == 0) return false;

	HEH pick_he = he_stack[rand() % he_stack.size()];

	//std::cout << "Tri change: " << mesh_.face_handle(pick_he).idx() << std::endl;

	FH oppo_cur_f = mesh_.opposite_face_handle(pick_he);
	FH oppo_next_f = mesh_.opposite_face_handle(mesh_.next_halfedge_handle(pick_he));
	FH oppo_prev_f = mesh_.opposite_face_handle(mesh_.prev_halfedge_handle(pick_he));
	std::vector<int> adj_seg_id;
	if (oppo_next_f.is_valid() &&
		(v_cont[mesh_.to_vertex_handle(pick_he).idx()] == 2
			|| seg_id_[oppo_cur_f.idx()] != seg_id_[oppo_next_f.idx()]))
	{
		adj_seg_id.emplace_back(seg_id_[oppo_next_f.idx()]);
	}

	if (oppo_prev_f.is_valid() &&
		(v_cont[mesh_.from_vertex_handle(pick_he).idx()] == 2
			|| seg_id_[oppo_cur_f.idx()] != seg_id_[oppo_prev_f.idx()]))
	{
		adj_seg_id.emplace_back(seg_id_[oppo_prev_f.idx()]);
	}

	if (adj_seg_id.size() == 0) return false;

	int pick_seg_id = adj_seg_id[rand() % adj_seg_id.size()];

	seg_id_[mesh_.face_handle(pick_he).idx()] = pick_seg_id;

	IdxUpdate();

	if (HighGenus())
	{
		seg_id_ = ori_seg;
		IdxUpdate();
		return false;
	}


	if (!cone_fixed &&
		cone_status_[mesh_.opposite_vh(pick_he).idx()]) {

		cone_status_[mesh_.opposite_vh(pick_he).idx()] = false;

		if (seg_id_[mesh_.face_handle(pick_he).idx()]
			== seg_id_[oppo_next_f.idx()])
		{
			cone_status_[mesh_.from_vertex_handle(pick_he).idx()] = true;
		}
		else
		{
			cone_status_[mesh_.to_vertex_handle(pick_he).idx()] = true;
		}
	}

	return true;
}

bool TopoChange::RandomAdjRemove(
	const Eigen::VectorXd& he_ang)
{
	//std::vector<int> ori_seg = seg_id_;

	//std::vector<OpenMesh::Vec3d> he_vec_n(mesh_.n_halfedges());
	//for (const HEH& he_h : mesh_.halfedges())
	//{
	//	he_vec_n[he_h.idx()] = mesh_.calc_edge_vector(he_h).normalized();
	//}

	//std::vector<int> v_cont;
	//VertexAdjSeam(v_cont);

	//std::vector<bool> anchor_status(mesh_.n_vertices(), false);
	//for (const VH& v_h : mesh_.vertices())
	//{
	//	if (v_cont[v_h.idx()] < 2) continue;

	//	if (v_cont[v_h.idx()] > 2)
	//	{
	//		anchor_status[v_h.idx()] = true;
	//		continue;
	//	}

	//	std::vector<int> adj_hes;
	//	for (const HEH& he_h : mesh_.voh_range(v_h))
	//	{
	//		if (seam_status_[mesh_.edge_handle(he_h).idx()])
	//		{
	//			adj_hes.emplace_back(he_h.idx());
	//		}
	//	}

	//	if (dot(he_vec_n[adj_hes[0]], he_vec_n[adj_hes[1]]) > ANCHOR_COS)
	//	{
	//		anchor_status[v_h.idx()] = true;
	//	}
	//}

	//// Cut the seams
	//std::vector<std::vector<HEH>> he_path;
	//for (const EH& e_h : mesh_.edges())
	//{
	//	if (!seam_status_[e_h.idx()]) continue;
	//	if (mesh_.is_boundary(e_h)) continue;

	//	HEH he_h = mesh_.halfedge_handle(e_h, 0);

	//	const FH& f_0 = mesh_.face_handle(he_h);
	//	const FH& f_1 = mesh_.opposite_face_handle(he_h);

	//	// Pick one halfedge of each edge
	//	if (seg_id_[f_0.idx()] < seg_id_[f_1.idx()])
	//	{
	//		he_h = mesh_.opposite_halfedge_handle(he_h);
	//	}

	//	// Separate the moved boundary
	//	std::vector<HEH> temp_path;
	//	temp_path.emplace_back(he_h);

	//	std::vector<std::vector<HEH>>::iterator path_it = he_path.begin();
	//	while (path_it != he_path.end())
	//	{
	//		// Is on the seam?
	//		// Can be merge with other pathes?
	//		if (!anchor_status[mesh_.to_vertex_handle(temp_path.back()).idx()]
	//			&& mesh_.to_vertex_handle(temp_path.back())
	//			== mesh_.from_vertex_handle((*path_it).front()))
	//		{
	//			temp_path.insert(temp_path.end(), (*path_it).begin(), (*path_it).end());
	//			path_it = he_path.erase(path_it);
	//		}
	//		else if (!anchor_status[mesh_.from_vertex_handle(temp_path.front()).idx()]
	//			&& mesh_.from_vertex_handle(temp_path.front())
	//			== mesh_.to_vertex_handle((*path_it).back()))
	//		{
	//			temp_path.insert(temp_path.begin(), (*path_it).begin(), (*path_it).end());
	//			path_it = he_path.erase(path_it);
	//		}
	//		else {
	//			path_it++;
	//		}
	//	}

	//	he_path.push_back(temp_path);
	//}

	//if (he_path.size() == 0) return false;

	//std::vector<bool> picked_boundaries(he_path.size(), true); // Archive of boundaries

	//int cont_0 = 0, cont_1 = 0;
	//std::vector<int> adj_id;
	//for (int path_id = 0; path_id < he_path.size(); path_id++)
	//{
	//	const HEH& start_he = he_path[path_id].front();
	//	const HEH& end_he = he_path[path_id].back();

	//	int seg_0 = seg_id_[mesh_.face_handle(start_he).idx()];
	//	int seg_1 = seg_id_[mesh_.opposite_face_handle(start_he).idx()];

	//	const VH& start_v = mesh_.from_vertex_handle(start_he);
	//	const VH& end_v = mesh_.to_vertex_handle(end_he);

	//	// Start V
	//	adj_id.clear();
	//	for (const FH& adj_f : mesh_.vf_range(start_v))
	//	{
	//		if (adj_id.size() == 0 || seg_id_[adj_f.idx()] != adj_id.back())
	//		{
	//			adj_id.emplace_back(seg_id_[adj_f.idx()]);
	//		}
	//	}
	//	if (adj_id.front() == adj_id.back()) adj_id.pop_back();

	//	cont_0 = 0;	cont_1 = 0;
	//	for (size_t i = 0; i < adj_id.size(); i++)
	//	{
	//		if (adj_id[i] == seg_0)
	//		{
	//			cont_0++;
	//		}

	//		if (adj_id[i] == seg_1)
	//		{
	//			cont_1++;
	//		}
	//	}

	//	if (cont_0 > 1 || cont_1 > 1)
	//	{
	//		picked_boundaries[path_id] = false;
	//		continue;
	//	}

	//	adj_id.clear();
	//	for (const FH& adj_f : mesh_.vf_range(end_v))
	//	{
	//		if (adj_id.size() == 0 || seg_id_[adj_f.idx()] != adj_id.back())
	//		{
	//			adj_id.emplace_back(seg_id_[adj_f.idx()]);
	//		}
	//	}
	//	if (adj_id.front() == adj_id.back()) adj_id.pop_back();

	//	cont_0 = 0;	cont_1 = 0;
	//	for (size_t i = 0; i < adj_id.size(); i++)
	//	{
	//		if (adj_id[i] == seg_0)
	//		{
	//			cont_0++;
	//		}

	//		if (adj_id[i] == seg_1)
	//		{
	//			cont_1++;
	//		}
	//	}

	//	if (cont_0 > 1 || cont_1 > 1)
	//	{
	//		picked_boundaries[path_id] = false;
	//		continue;
	//	}
	//}

	//for (int path_id = 0; path_id < he_path.size(); path_id++)
	//{
	//	if (!picked_boundaries[path_id]) continue;

	//	const HEH& front_he = he_path[path_id].front();
	//	const HEH& back_he = he_path[path_id].back();

	//	const VH& start_v = mesh_.from_vertex_handle(front_he);
	//	const VH& end_v = mesh_.to_vertex_handle(back_he);

	//	//if (v_cont[start_v.idx()] > 2 || v_cont[end_v.idx()] > 2)
	//	//{
	//	//	picked_boundaries[path_id] = false;
	//	//	continue;
	//	//}

	//	int cont = 0;

	//	double max_start = -DBL_MAX, min_start = DBL_MAX;
	//	double sum_start = 0;
	//	HEH start_he = front_he;
	//	HEH cur_he = start_he;
	//	do
	//	{
	//		sum_start += he_ang[mesh_.next_halfedge_handle(cur_he).idx()];

	//		cur_he = mesh_.opposite_halfedge_handle(
	//			mesh_.prev_halfedge_handle(cur_he));

	//		if (cur_he != start_he && seam_status_[mesh_.edge_handle(cur_he).idx()])
	//		{
	//			min_start = std::min(sum_start, min_start);
	//			max_start = std::max(sum_start, max_start);
	//		}

	//	} while (cont < 2 * mesh_.valence(start_v) && cur_he != start_he);

	//	double max_end = -DBL_MAX, min_end = DBL_MAX;
	//	double sum_end = 0;
	//	start_he = back_he;
	//	cur_he = start_he;
	//	do
	//	{
	//		sum_end += he_ang[mesh_.prev_halfedge_handle(cur_he).idx()];

	//		cur_he = mesh_.opposite_halfedge_handle(
	//			mesh_.next_halfedge_handle(cur_he));

	//		if (cur_he != start_he && seam_status_[mesh_.edge_handle(cur_he).idx()])
	//		{
	//			min_end = std::min(sum_end, min_end);
	//			max_end = std::max(sum_end, max_end);
	//		}

	//	} while (cont < 2 * mesh_.valence(end_v) && cur_he != start_he);

	//	if ((min_start < sum_start / 2
	//		&& min_end < sum_end / 2) ||
	//		(max_start > sum_start / 2
	//			&& max_end > sum_end / 2)) continue;

	//	picked_boundaries[path_id] = false;
	//}


	//std::vector<int> picked_ids;

	//// Numbers of cones on boundaries
	//std::vector<int> cone_num(he_path.size(), 0);
	//for (size_t i = 0; i < he_path.size(); i++)
	//{
	//	for (size_t j = 1; j < he_path[i].size(); j++)
	//	{
	//		const VH& from_v = mesh_.from_vertex_handle(he_path[i][j]);

	//		if (cone_status_[from_v.idx()])
	//		{
	//			cone_num[i]++;
	//		}
	//	}
	//}

	//int min_num = 1;
	////int min_num = DBL_MAX;
	////for (size_t i = 0; i < cone_num.size(); i++)
	////{
	////	if (cone_num[i] == 0) continue;
	////	
	////	min_num = std::min(min_num, cone_num[i]);
	////}

	//for (size_t i = 0; i < he_path.size(); i++)
	//{
	//	if (cone_num[i] != min_num)
	//	{
	//		picked_boundaries[i] = false;
	//	}
	//}

	//for (size_t i = 0; i < picked_boundaries.size(); i++)
	//{
	//	if (picked_boundaries[i])
	//	{
	//		picked_ids.push_back(i);
	//	}
	//}

	//if (picked_ids.size() == 0) return false;

	//const std::vector<HEH>& pick_path =
	//	he_path[picked_ids[rand() % picked_ids.size()]];

	//int f_0 = mesh_.face_handle(pick_path.front()).idx();
	//int f_1 = mesh_.opposite_face_handle(pick_path.front()).idx();

	////OutputSeg("seg.txt");
	////std::cout << f_0 << "  " << f_1 << std::endl;

	//int seg_0 = std::min(seg_id_[f_0], seg_id_[f_1]);
	//int seg_1 = std::max(seg_id_[f_0], seg_id_[f_1]);

	//std::vector<int> set_0, set_1;
	//for (int i = 1; i < pick_path.size(); i++)
	//{
	//	for (const FH& adj_f : mesh_.vf_range(mesh_.from_vertex_handle(pick_path[i])))
	//	{
	//		if (seg_id_[adj_f.idx()] == seg_0)
	//		{
	//			set_0.emplace_back(adj_f.idx());
	//		}
	//		else if (seg_id_[adj_f.idx()] == seg_1)
	//		{
	//			set_1.emplace_back(adj_f.idx());
	//		}
	//	}

	//}
	//for (const HEH& he_h : pick_path)
	//{
	//	for (const FH& adj_f : mesh_.vf_range(mesh_.from_vertex_handle(he_h)))
	//	{
	//		if (seg_id_[adj_f.idx()] == seg_0)
	//		{
	//			set_0.emplace_back(adj_f.idx());
	//		}
	//		else if (seg_id_[adj_f.idx()] == seg_1)
	//		{
	//			set_1.emplace_back(adj_f.idx());
	//		}
	//	}
	//}

	//for (const FH& adj_f : mesh_.vf_range(mesh_.to_vertex_handle(pick_path.back())))
	//{
	//	if (seg_id_[adj_f.idx()] == seg_0)
	//	{
	//		set_0.emplace_back(adj_f.idx());
	//	}
	//	else if (seg_id_[adj_f.idx()] == seg_1)
	//	{
	//		set_1.emplace_back(adj_f.idx());
	//	}
	//}

	//std::sort(set_0.begin(), set_0.end());
	//set_0.erase(std::unique(set_0.begin(), set_0.end()), set_0.end());

	//std::sort(set_1.begin(), set_1.end());
	//set_1.erase(std::unique(set_1.begin(), set_1.end()), set_1.end());

	//double area_0 = 0, area_1 = 0;
	//for (const int& f_id : set_0)
	//{
	//	area_0 += mesh_.calc_face_area(FH(f_id));
	//}

	//for (const int& f_id : set_1)
	//{
	//	area_1 += mesh_.calc_face_area(FH(f_id));
	//}

	//if (area_0 > area_1)
	//{
	//	for (const int& f_id : set_1)
	//	{
	//		seg_id_[f_id] = seg_0;
	//	}
	//}
	//else if (area_0 < area_1)
	//{
	//	for (const int& f_id : set_0)
	//	{
	//		seg_id_[f_id] = seg_1;
	//	}
	//}
	//else
	//{
	//	return false;
	//}

	//IdxUpdate();

	//std::vector<bool> new_seam_v(mesh_.n_vertices(), false);
	//for (const EH& e_h : mesh_.edges())
	//{
	//	const HEH& he_0 = mesh_.halfedge_handle(e_h, 0);

	//	const VH& to_v = mesh_.to_vertex_handle(he_0);
	//	const VH& from_v = mesh_.from_vertex_handle(he_0);

	//	new_seam_v[to_v.idx()] = true;
	//	new_seam_v[from_v.idx()] = true;
	//}

	//for (const VH& v_h : mesh_.vertices())
	//{
	//	if (!cone_status_[v_h.idx()]) continue;

	//	if (cone_status_[v_h.idx()] && !new_seam_v[v_h.idx()])
	//	{
	//		bool find_status = false;
	//		for (const VH& adj_v : mesh_.vv_range(v_h))
	//		{
	//			if (new_seam_v[adj_v.idx()])
	//			{
	//				find_status = true;
	//				cone_status_[v_h.idx()] = false;
	//				cone_status_[adj_v.idx()] = true;
	//				break;
	//			}
	//		}

	//		if (find_status)
	//		{
	//			continue;
	//		}
	//		else
	//		{
	//			return false;
	//		}
	//	}
	//}
	//
	//std::vector<int> new_v_cont;
	//VertexAdjSeam(new_v_cont);
	//for (int i = 0; i < cone_status_.size(); i++)
	//{
	//	if (cone_status_[i] && new_v_cont[i] == 0)
	//	{
	//		seg_id_ = ori_seg;
	//		IdxUpdate();
	//		return false;
	//	}
	//}

	//return true;

	std::vector<int> ori_seg = seg_id_;

	std::vector<OpenMesh::Vec3d> he_vec_n(mesh_.n_halfedges());
	for (const HEH& he_h : mesh_.halfedges())
	{
		he_vec_n[he_h.idx()] = mesh_.calc_edge_vector(he_h).normalized();
	}

	std::vector<int> v_cont;
	VertexAdjSeam(v_cont);

	std::vector<bool> anchor_status(mesh_.n_vertices(), false);
	for (const VH& v_h : mesh_.vertices())
	{
		if (v_cont[v_h.idx()] < 2) continue;

		if (v_cont[v_h.idx()] > 2)
		{
			anchor_status[v_h.idx()] = true;
			continue;
		}

		std::vector<int> adj_hes;
		for (const HEH& he_h : mesh_.voh_range(v_h))
		{
			if (seam_status_[mesh_.edge_handle(he_h).idx()])
			{
				adj_hes.emplace_back(he_h.idx());
			}
		}

		if (dot(he_vec_n[adj_hes[0]], he_vec_n[adj_hes[1]]) > ANCHOR_COS)
		{
			anchor_status[v_h.idx()] = true;
		}
	}

	// Cut the seams
	std::vector<std::vector<HEH>> he_path;
	for (const EH& e_h : mesh_.edges())
	{
		if (!seam_status_[e_h.idx()]) continue;
		if (mesh_.is_boundary(e_h)) continue;

		HEH he_h = mesh_.halfedge_handle(e_h, 0);

		const FH& f_0 = mesh_.face_handle(he_h);
		const FH& f_1 = mesh_.opposite_face_handle(he_h);

		// Pick one halfedge of each edge
		if (seg_id_[f_0.idx()] < seg_id_[f_1.idx()])
		{
			he_h = mesh_.opposite_halfedge_handle(he_h);
		}

		// Separate the moved boundary
		std::vector<HEH> temp_path;
		temp_path.emplace_back(he_h);

		std::vector<std::vector<HEH>>::iterator path_it = he_path.begin();
		while (path_it != he_path.end())
		{
			// Is on the seam?
			// Can be merge with other pathes?
			if (!anchor_status[mesh_.to_vertex_handle(temp_path.back()).idx()]
				&& mesh_.to_vertex_handle(temp_path.back())
				== mesh_.from_vertex_handle((*path_it).front()))
			{
				temp_path.insert(temp_path.end(), (*path_it).begin(), (*path_it).end());
				path_it = he_path.erase(path_it);
			}
			else if (!anchor_status[mesh_.from_vertex_handle(temp_path.front()).idx()]
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

	if (he_path.size() == 0) return false;

	std::vector<bool> picked_boundaries(he_path.size(), true); // Archive of boundaries

	int cont_0 = 0, cont_1 = 0;
	std::vector<int> adj_id;
	for (int path_id = 0; path_id < he_path.size(); path_id++)
	{
		const HEH& start_he = he_path[path_id].front();
		const HEH& end_he = he_path[path_id].back();

		int seg_0 = seg_id_[mesh_.face_handle(start_he).idx()];
		int seg_1 = seg_id_[mesh_.opposite_face_handle(start_he).idx()];

		const VH& start_v = mesh_.from_vertex_handle(start_he);
		const VH& end_v = mesh_.to_vertex_handle(end_he);

		// Start V
		adj_id.clear();
		for (const FH& adj_f : mesh_.vf_range(start_v))
		{
			if (adj_id.size() == 0 || seg_id_[adj_f.idx()] != adj_id.back())
			{
				adj_id.emplace_back(seg_id_[adj_f.idx()]);
			}
		}
		if (adj_id.front() == adj_id.back()) adj_id.pop_back();

		cont_0 = 0;	cont_1 = 0;
		for (size_t i = 0; i < adj_id.size(); i++)
		{
			if (adj_id[i] == seg_0)
			{
				cont_0++;
			}

			if (adj_id[i] == seg_1)
			{
				cont_1++;
			}
		}

		if (cont_0 > 1 || cont_1 > 1)
		{
			picked_boundaries[path_id] = false;
			continue;
		}

		adj_id.clear();
		for (const FH& adj_f : mesh_.vf_range(end_v))
		{
			if (adj_id.size() == 0 || seg_id_[adj_f.idx()] != adj_id.back())
			{
				adj_id.emplace_back(seg_id_[adj_f.idx()]);
			}
		}
		if (adj_id.front() == adj_id.back()) adj_id.pop_back();

		cont_0 = 0;	cont_1 = 0;
		for (size_t i = 0; i < adj_id.size(); i++)
		{
			if (adj_id[i] == seg_0)
			{
				cont_0++;
			}

			if (adj_id[i] == seg_1)
			{
				cont_1++;
			}
		}

		if (cont_0 > 1 || cont_1 > 1)
		{
			picked_boundaries[path_id] = false;
			continue;
		}
	}

	//std::ofstream temp_status_0("status_0.txt");
	//for (int i = 0; i < picked_boundaries.size(); i++)
	//{
	//	temp_status_0 << picked_boundaries[i] << std::endl;
	//}
	//temp_status_0.close();


	for (int path_id = 0; path_id < he_path.size(); path_id++)
	{
		if (!picked_boundaries[path_id]) continue;

		const HEH& front_he = he_path[path_id].front();
		const HEH& back_he = he_path[path_id].back();

		const VH& start_v = mesh_.from_vertex_handle(front_he);
		const VH& end_v = mesh_.to_vertex_handle(back_he);

		//if (v_cont[start_v.idx()] > 2 || v_cont[end_v.idx()] > 2)
		//{
		//	picked_boundaries[path_id] = false;
		//	continue;
		//}

		int cont = 0;

		double max_start = -DBL_MAX, min_start = DBL_MAX;
		double sum_start = 0;
		HEH start_he = front_he;
		HEH cur_he = start_he;
		do
		{
			sum_start += he_ang[mesh_.next_halfedge_handle(cur_he).idx()];

			cur_he = mesh_.opposite_halfedge_handle(
				mesh_.prev_halfedge_handle(cur_he));

			if (cur_he != start_he && seam_status_[mesh_.edge_handle(cur_he).idx()])
			{
				min_start = std::min(sum_start, min_start);
				max_start = std::max(sum_start, max_start);
			}

		} while (cont < 2 * mesh_.valence(start_v) && cur_he != start_he);

		double max_end = -DBL_MAX, min_end = DBL_MAX;
		double sum_end = 0;
		start_he = back_he;
		cur_he = start_he;
		do
		{
			sum_end += he_ang[mesh_.prev_halfedge_handle(cur_he).idx()];

			cur_he = mesh_.opposite_halfedge_handle(
				mesh_.next_halfedge_handle(cur_he));

			if (cur_he != start_he && seam_status_[mesh_.edge_handle(cur_he).idx()])
			{
				min_end = std::min(sum_end, min_end);
				max_end = std::max(sum_end, max_end);
			}

		} while (cont < 2 * mesh_.valence(end_v) && cur_he != start_he);

		if ((min_start < sum_start / 2
			&& min_end < sum_end / 2) ||
			(max_start > sum_start / 2
				&& max_end > sum_end / 2)) continue;

		picked_boundaries[path_id] = false;
	}

	//std::ofstream temp_status_1("status_1.txt");
	//for (int i = 0; i < picked_boundaries.size(); i++)
	//{
	//	temp_status_1 << picked_boundaries[i] << std::endl;
	//}
	//temp_status_1.close();


	//for (int i = 0; i < he_path.size(); i++)
	//{
	//	if (picked_boundaries[i])
	//	{
	//		std::ofstream out_file("out_f_" + std::to_string(i) + ".txt");
	//		for (const HEH& he_h : he_path[i])
	//		{
	//			out_file << mesh_.face_handle(he_h).idx() << std::endl;
	//		}
	//		out_file.close();
	//	}

	//}

	//// Numbers of cones on boundaries
	//std::vector<int> cone_num(he_path.size(), 0);
	//for (size_t i = 0; i < he_path.size(); i++)
	//{
	//	for (size_t j = 1; j < he_path[i].size(); j++)
	//	{
	//		const VH& from_v = mesh_.from_vertex_handle(he_path[i][j]);

	//		if (cone_status_[from_v.idx()])
	//		{
	//			cone_num[i]++;
	//		}
	//	}
	//}

	//int min_num = 1;
	////int min_num = DBL_MAX;
	////for (size_t i = 0; i < cone_num.size(); i++)
	////{
	////	if (cone_num[i] == 0) continue;
	////	
	////	min_num = std::min(min_num, cone_num[i]);
	////}

	//for (size_t i = 0; i < he_path.size(); i++)
	//{
	//	if (cone_num[i] != min_num)
	//	{
	//		picked_boundaries[i] = false;
	//	}
	//}

	std::vector<int> picked_ids;
	for (size_t i = 0; i < picked_boundaries.size(); i++)
	{
		if (picked_boundaries[i])
		{
			picked_ids.push_back(i);
		}
	}

	if (picked_ids.size() == 0) return false;

	const std::vector<HEH>& pick_path =
		he_path[picked_ids[rand() % picked_ids.size()]];

	//std::ofstream out_pick("pick_f.txt");
	//for (const HEH& he_h : pick_path)
	//{
	//	out_pick << mesh_.face_handle(he_h).idx() << std::endl;
	//}
	//out_pick.close();

	int f_0 = mesh_.face_handle(pick_path.front()).idx();
	int f_1 = mesh_.opposite_face_handle(pick_path.front()).idx();

	//OutputSeg("seg.txt");
	//std::cout << f_0 << "  " << f_1 << std::endl;

	int seg_0 = std::min(seg_id_[f_0], seg_id_[f_1]);
	int seg_1 = std::max(seg_id_[f_0], seg_id_[f_1]);

	std::vector<int> set_0, set_1;
	for (const HEH& he_h : pick_path)
	{
		for (const FH& adj_f : mesh_.vf_range(mesh_.from_vertex_handle(he_h)))
		{
			if (seg_id_[adj_f.idx()] == seg_0)
			{
				set_0.emplace_back(adj_f.idx());
			}
			else if (seg_id_[adj_f.idx()] == seg_1)
			{
				set_1.emplace_back(adj_f.idx());
			}
		}
	}

	for (const FH& adj_f : mesh_.vf_range(mesh_.to_vertex_handle(pick_path.back())))
	{
		if (seg_id_[adj_f.idx()] == seg_0)
		{
			set_0.emplace_back(adj_f.idx());
		}
		else if (seg_id_[adj_f.idx()] == seg_1)
		{
			set_1.emplace_back(adj_f.idx());
		}
	}

	std::sort(set_0.begin(), set_0.end());
	set_0.erase(std::unique(set_0.begin(), set_0.end()), set_0.end());

	std::sort(set_1.begin(), set_1.end());
	set_1.erase(std::unique(set_1.begin(), set_1.end()), set_1.end());

	double area_0 = 0, area_1 = 0;
	for (const int& f_id : set_0)
	{
		area_0 += mesh_.calc_face_area(FH(f_id));
	}

	for (const int& f_id : set_1)
	{
		area_1 += mesh_.calc_face_area(FH(f_id));
	}

	if (area_0 > area_1)
	{
		for (const int& f_id : set_1)
		{
			seg_id_[f_id] = seg_0;
		}
	}
	else if (area_0 < area_1)
	{
		for (const int& f_id : set_0)
		{
			seg_id_[f_id] = seg_1;
		}
	}
	else
	{
		return false;
	}

	IdxUpdate();

	if (HighGenus())
	{
		seg_id_ = ori_seg;
		IdxUpdate();
		return false;
	}

	std::vector<bool> new_seam_v(mesh_.n_vertices(), false);
	for (const EH& e_h : mesh_.edges())
	{
		if (!seam_status_[e_h.idx()]) continue;

		const HEH& he_0 = mesh_.halfedge_handle(e_h, 0);

		const VH& to_v = mesh_.to_vertex_handle(he_0);
		const VH& from_v = mesh_.from_vertex_handle(he_0);

		new_seam_v[to_v.idx()] = true;
		new_seam_v[from_v.idx()] = true;
	}

	for (const VH& v_h : mesh_.vertices())
	{
		if (!cone_status_[v_h.idx()]) continue;

		if (cone_status_[v_h.idx()] && !new_seam_v[v_h.idx()])
		{
			bool find_status = false;
			for (const VH& adj_v : mesh_.vv_range(v_h))
			{
				if (new_seam_v[adj_v.idx()])
				{
					find_status = true;
					cone_status_[v_h.idx()] = false;
					cone_status_[adj_v.idx()] = true;
					break;
				}
			}

			if (find_status)
			{
				continue;
			}
			else
			{
				return false;
			}
		}
	}

	std::vector<int> new_v_cont;
	VertexAdjSeam(new_v_cont);
	for (int i = 0; i < cone_status_.size(); i++)
	{
		if (cone_status_[i] && new_v_cont[i] == 0)
		{
			//std::cout << i << std::endl;

			//std::ofstream cone_file("cone_file.txt");
			//std::ofstream cont_file("cont_file.txt");
			//for (int cone_id = 0; cone_id < cone_status_.size(); cone_id++)
			//{
			//	if (cone_status_[cone_id])
			//	{
			//		cone_file << cone_id << std::endl;
			//	}
			//	cont_file << new_v_cont[cone_id] << std::endl;
			//}
			//cone_file.close(); cont_file.close();

			seg_id_ = ori_seg;
			IdxUpdate();
			return false;
		}
	}

	return true;
}

bool TopoChange::RandomDijkstra(const Eigen::VectorXd& e_l, const Eigen::VectorXd& e_w)
{
	Eigen::VectorXd e_weight = e_l.cwiseProduct(e_w);

	std::vector<int> v_cont;
	VertexAdjSeam(v_cont);

	std::vector<bool> anchor_status(mesh_.n_vertices(), false);
	for (const VH& v_h : mesh_.vertices())
	{
		if (v_cont[v_h.idx()] > 2 || cone_status_[v_h.idx()])
		{
			anchor_status[v_h.idx()] = true;
		}
	}

	// Cut the seams
	std::vector<std::vector<HEH>> archor_cut_path;
	std::vector<std::vector<HEH>> cross_cut_path;
	std::vector<HEH> temp_path;
	std::vector<std::vector<HEH>>::iterator path_it;
	for (const EH& e_h : mesh_.edges())
	{
		if (!seam_status_[e_h.idx()]) continue;

		HEH he_h = mesh_.halfedge_handle(e_h, 0);

		const FH& f_0 = mesh_.face_handle(he_h);
		const FH& f_1 = mesh_.opposite_face_handle(he_h);

		// Pick one halfedge of each edge
		if (!f_0.is_valid() 
			|| (f_1.is_valid() && seg_id_[f_0.idx()] < seg_id_[f_1.idx()]))
		{
			he_h = mesh_.opposite_halfedge_handle(he_h);
		}

		// Separate the moved boundary
		temp_path.clear();
		temp_path.emplace_back(he_h);

		path_it = cross_cut_path.begin();
		while (path_it != cross_cut_path.end())
		{
			// Is on the seam?
			// Can be merge with other pathes?
			//if (v_cont[mesh_.to_vertex_handle(temp_path.back()).idx()] == 2
			if (!anchor_status[mesh_.to_vertex_handle(temp_path.back()).idx()]
				&& mesh_.to_vertex_handle(temp_path.back())
				== mesh_.from_vertex_handle((*path_it).front()))
			{
				temp_path.insert(temp_path.end(), (*path_it).begin(), (*path_it).end());
				path_it = cross_cut_path.erase(path_it);
			}
			//else if (v_cont[mesh_.from_vertex_handle(temp_path.front()).idx()] == 2
			else if (!anchor_status[mesh_.from_vertex_handle(temp_path.front()).idx()]
				&& mesh_.from_vertex_handle(temp_path.front())
				== mesh_.to_vertex_handle((*path_it).back()))
			{
				temp_path.insert(temp_path.begin(), (*path_it).begin(), (*path_it).end());
				path_it = cross_cut_path.erase(path_it);
			}
			else {
				path_it++;
			}
		}

		cross_cut_path.push_back(temp_path);

		if (mesh_.is_boundary(e_h)) continue;

		// Separate the moved boundary
		temp_path.clear();
		temp_path.emplace_back(he_h);

		path_it = archor_cut_path.begin();
		while (path_it != archor_cut_path.end())
		{
			// Is on the seam?
			// Can be merge with other pathes?
			if (!anchor_status[mesh_.to_vertex_handle(temp_path.back()).idx()]
				&& mesh_.to_vertex_handle(temp_path.back())
				== mesh_.from_vertex_handle((*path_it).front()))
			{
				temp_path.insert(temp_path.end(), (*path_it).begin(), (*path_it).end());
				path_it = archor_cut_path.erase(path_it);
			}
			else if (!anchor_status[mesh_.from_vertex_handle(temp_path.front()).idx()]
				&& mesh_.from_vertex_handle(temp_path.front())
				== mesh_.to_vertex_handle((*path_it).back()))
			{
				temp_path.insert(temp_path.begin(), (*path_it).begin(), (*path_it).end());
				path_it = archor_cut_path.erase(path_it);
			}
			else {
				path_it++;
			}
		}

		archor_cut_path.push_back(temp_path);
	}

	if (archor_cut_path.size() == 0) return false;

	Eigen::VectorXd cut_len;;
	cut_len.setConstant(archor_cut_path.size(), 0);
	for (size_t i = 0; i < archor_cut_path.size(); i++)
	{
		for (const HEH& he_h : archor_cut_path[i])
		{
			//cut_len[i] += e_weight[mesh_.edge_handle(he_h).idx()];
			cut_len[i] += e_l[mesh_.edge_handle(he_h).idx()];
		}
	}

	double max_len = cut_len.maxCoeff();
	std::vector<int> picked_ids;
	for (int path_id = 0; path_id < archor_cut_path.size(); path_id++)
	{
		if (archor_cut_path[path_id].size() < 2 
			|| cut_len[path_id] < max_len * LENGTH_RATE) continue;

		picked_ids.emplace_back(path_id);
	}

	if (picked_ids.size() == 0) return false;

	const std::vector<HEH>& picked_path = archor_cut_path[picked_ids[rand() % picked_ids.size()]];
	const VH& start_v = mesh_.from_vertex_handle(picked_path.front());
	const VH& end_v = mesh_.to_vertex_handle(picked_path.back());

	//std::ofstream out_f("out_f.txt");
	//for (const HEH& he_h : picked_path)
	//{
	//	out_f << mesh_.face_handle(he_h).idx() << std::endl;
	//}
	//out_f.close();

	//for (size_t i = 0; i < cross_cut_path.size(); i++)
	//{
	//	std::ofstream out_fff("out_fff" + std::to_string(i) + ".txt");
	//	for (const HEH& he_h : cross_cut_path[i])
	//	{
	//		out_fff << mesh_.face_handle(he_h).idx() << std::endl;
	//	}
	//	out_fff.close();
	//}

	std::vector<int> start_v_set;
	if (v_cont[start_v.idx()] == 2)
	{
		start_v_set.emplace_back(start_v.idx());
	}
	else
	{
		HEH left_he = picked_path.front();
		do
		{
			left_he = mesh_.opposite_halfedge_handle(
				mesh_.prev_halfedge_handle(left_he));
		} while (!seam_status_[mesh_.edge_handle(left_he).idx()]);

		HEH right_he = picked_path.front();
		do
		{
			right_he = mesh_.next_halfedge_handle(
				mesh_.opposite_halfedge_handle(right_he));
		} while (!seam_status_[mesh_.edge_handle(right_he).idx()]);

		if (left_he.is_valid() && right_he.is_valid() 
			&& mesh_.face_handle(left_he).is_valid()
			&& mesh_.face_handle(right_he).is_valid()
			&& seg_id_[mesh_.face_handle(left_he).idx()]
			== seg_id_[mesh_.face_handle(right_he).idx()])
		{
			left_he.reset();
		}
		if (left_he.is_valid() && right_he.is_valid()
			&& mesh_.opposite_face_handle(left_he).is_valid()
			&& mesh_.opposite_face_handle(right_he).is_valid()
			&& seg_id_[mesh_.opposite_face_handle(left_he).idx()]
			== seg_id_[mesh_.opposite_face_handle(right_he).idx()])
		{
			right_he.reset();
		}

		for (const auto& single_path : cross_cut_path)
		{
			//std::cout << mesh_.face_handle(left_he).idx() 
			//	<< " " << mesh_.face_handle(right_he).idx()
			//	<< " " << mesh_.opposite_face_handle(left_he).idx()
			//	<< " " << mesh_.opposite_face_handle(right_he).idx()
			//	<< std::endl;

			if (mesh_.opposite_halfedge_handle(left_he) == single_path.back()
				|| left_he == single_path.front()
				|| right_he == single_path.front()
				|| mesh_.opposite_halfedge_handle(right_he) == single_path.back())
			{
				for (const HEH& he_h : single_path)
				{
					const VH& to_v = mesh_.to_vertex_handle(he_h);
					const VH& from_v = mesh_.from_vertex_handle(he_h);

					start_v_set.emplace_back(to_v.idx());
					start_v_set.emplace_back(from_v.idx());
				}

				//if (left_he == mesh_.opposite_halfedge_handle(single_path.back())
				//	|| left_he == single_path.front())
				//{
				//	std::ofstream out_f_left_start("out_f_left_start.txt");
				//	for (const HEH& he_h : single_path)
				//	{
				//		out_f_left_start << mesh_.face_handle(he_h).idx() << std::endl;
				//	}
				//	out_f_left_start.close();
				//}
				//else
				//{
				//	std::ofstream out_f_right_start("out_f_right_start.txt");
				//	for (const HEH& he_h : single_path)
				//	{
				//		out_f_right_start << mesh_.face_handle(he_h).idx() << std::endl;
				//	}
				//	out_f_right_start.close();
				//}
			}
		}

		std::sort(start_v_set.begin(), start_v_set.end());
		start_v_set.erase(unique(start_v_set.begin(), start_v_set.end()), start_v_set.end());
	}

	//std::ofstream start_v_out("out_v_start.txt");
	//for (const int& v_id : start_v_set)
	//{
	//	start_v_out << v_id << std::endl;
	//}
	//start_v_out.close();

	std::vector<int> end_v_set;
	if (v_cont[end_v.idx()] == 2)
	{
		end_v_set.emplace_back(end_v.idx());
	}
	else
	{
		HEH left_he = picked_path.back();
		do
		{
			left_he = mesh_.prev_halfedge_handle(
				mesh_.opposite_halfedge_handle(left_he));
		} while (!seam_status_[mesh_.edge_handle(left_he).idx()]);

		HEH right_he = picked_path.back();
		do
		{
			right_he = mesh_.opposite_halfedge_handle(
				mesh_.next_halfedge_handle(right_he));
		} while (!seam_status_[mesh_.edge_handle(right_he).idx()]);

		if (left_he.is_valid() && right_he.is_valid()
			&& seg_id_[mesh_.face_handle(left_he).idx()]
			== seg_id_[mesh_.face_handle(right_he).idx()])
		{
			right_he.reset();
		}
		
		if (left_he.is_valid() && right_he.is_valid()
			&& seg_id_[mesh_.opposite_face_handle(left_he).idx()]
			== seg_id_[mesh_.opposite_face_handle(right_he).idx()])
		{
			left_he.reset();
		}

		for (const auto& single_path : cross_cut_path)
		{
			if (left_he == single_path.back()
				|| left_he == mesh_.opposite_halfedge_handle(single_path.front())
				|| right_he == mesh_.opposite_halfedge_handle(single_path.front())
				|| right_he == single_path.back())
			{
				for (const HEH& he_h : single_path)
				{
					const VH& to_v = mesh_.to_vertex_handle(he_h);
					const VH& from_v = mesh_.from_vertex_handle(he_h);

					end_v_set.emplace_back(to_v.idx());
					end_v_set.emplace_back(from_v.idx());
				}

				//if (left_he == mesh_.opposite_halfedge_handle(single_path.back())
				//	|| left_he == single_path.front())
				//{
				//	std::ofstream out_f_left_end("out_f_left_end.txt");
				//	for (const HEH& he_h : single_path)
				//	{
				//		out_f_left_end << mesh_.face_handle(he_h).idx() << std::endl;
				//	}
				//	out_f_left_end.close();
				//}
				//else
				//{
				//	std::ofstream out_f_right_end("out_f_right_end.txt");
				//	for (const HEH& he_h : single_path)
				//	{
				//		out_f_right_end << mesh_.face_handle(he_h).idx() << std::endl;
				//	}
				//	out_f_right_end.close();
				//}
			}
		}

		std::sort(end_v_set.begin(), end_v_set.end());
		end_v_set.erase(unique(end_v_set.begin(), end_v_set.end()), end_v_set.end());
	}

	//std::ofstream end_v_out("out_v_end.txt");
	//for (const int& v_id : end_v_set)
	//{
	//	end_v_out << v_id << std::endl;
	//}
	//end_v_out.close();

	int seg_0 = seg_id_[mesh_.face_handle(picked_path.front()).idx()];
	int seg_1 = seg_id_[mesh_.opposite_face_handle(picked_path.front()).idx()];

	std::vector<bool> inner_v_status(mesh_.n_edges(), false);
	for (const EH& e_h : mesh_.edges())
	{
		if (seam_status_[e_h.idx()]) continue;
		if (mesh_.is_boundary(e_h)) continue;

		const HEH he_0 = mesh_.halfedge_handle(e_h, 0);
		const FH f_0 = mesh_.face_handle(he_0);

		if (seg_id_[f_0.idx()] == seg_0
			|| seg_id_[f_0.idx()] == seg_1)
		{
			const VH to_v = mesh_.to_vertex_handle(he_0);
			const VH from_v = mesh_.from_vertex_handle(he_0);
			inner_v_status[to_v.idx()] = true;
			inner_v_status[from_v.idx()] = true;
		}
	}

	for (const EH& e_h : mesh_.edges())
	{
		const HEH he_0 = mesh_.halfedge_handle(e_h, 0);
		if (seam_status_[e_h.idx()] || mesh_.is_boundary(e_h))
		{
			const VH to_v = mesh_.to_vertex_handle(he_0);
			const VH from_v = mesh_.from_vertex_handle(he_0);

			inner_v_status[to_v.idx()] = false;
			inner_v_status[from_v.idx()] = false;
		}
	}

	double dis, min_dis = DBL_MAX;
	std::vector<HEH> out_temp_path, min_path;
	if (start_v_set.size() < end_v_set.size())
	{
		for (int i = 0; i < picked_path.size(); i++)
		{
			inner_v_status[mesh_.to_vertex_handle(picked_path[i]).idx()] = true;
		}

		for (int start_id : start_v_set)
		{
			//std::cout << " ========== Dij ===========" << std::endl;
			Dijkstra(start_id, end_v_set, inner_v_status, e_weight, dis, out_temp_path);

			//std::ofstream temp_path_f("out_f_" + std::to_string(start_id) + ".txt");
			//for (const HEH& he_h : out_temp_path)
			//{
			//	temp_path_f << mesh_.face_handle(he_h).idx() << std::endl;
			//}
			//temp_path_f.close();

			if (dis < min_dis)
			{
				min_dis = dis;
				min_path = out_temp_path;
			}
		}
	}
	else
	{
		for (int i = 0; i < picked_path.size(); i++)
		{
			inner_v_status[mesh_.from_vertex_handle(picked_path[i]).idx()] = true;
		}

		for (int end_id : end_v_set)
		{
			//std::cout << " ========== Dij ===========" << std::endl;
			Dijkstra(end_id, start_v_set, inner_v_status, e_weight, dis, out_temp_path);

			//std::ofstream temp_path_f("out_f_" + std::to_string(end_id) + ".txt");
			//for (const HEH& he_h : out_temp_path)
			//{
			//	temp_path_f << mesh_.face_handle(he_h).idx() << std::endl;
			//}
			//temp_path_f.close();

			if (dis < min_dis)
			{
				min_dis = dis;
				min_path = out_temp_path;
			}
		}
	}

	std::vector<bool> ori_seam = seam_status_;

	for (const HEH& he_h : picked_path)
	{
		seam_status_[mesh_.edge_handle(he_h).idx()] = false;
	}

	//std::ofstream picked_path_f("out_f_picked_path.txt");
	//for (const HEH& he_h : picked_path)
	//{
	//	picked_path_f << mesh_.face_handle(he_h).idx() << std::endl;
	//}
	//picked_path_f.close();

	for (const HEH& he_h : min_path)
	{
		seam_status_[mesh_.edge_handle(he_h).idx()] = true;
	}

	//std::ofstream min_path_f("out_f_min_path.txt");
	//for (const HEH& he_h : min_path)
	//{
	//	min_path_f << mesh_.face_handle(he_h).idx() << std::endl;
	//}
	//min_path_f.close();

	SeamUpdate();
	IdxUpdate();

	if (HighGenus())
	{
		seam_status_ = ori_seam;
		SeamUpdate();
		IdxUpdate();
		return false;
	}

	return true;
}

void TopoChange::EdgeAndPatchToSeam(const std::vector<std::vector<HEH>>& path_he,
	std::vector<int>& edge_seam_id,
	std::map<std::pair<int, int>, std::vector<int>>& patch_id_to_seam_id)
{
	edge_seam_id.assign(mesh_.n_edges(), -1);
	patch_id_to_seam_id.clear();

	// Compute index for each boundaries
	for (int path_id = 0; path_id < path_he.size(); path_id++)
	{
		// Index for boundaries
		for (const HEH& he_h : path_he[path_id])
		{
			edge_seam_id[mesh_.edge_handle(he_h).idx()] = path_id;
		}

		// Indices of patches and their boundaries
		const HEH& he_h = path_he[path_id].front();
		const FH& f_0 = mesh_.face_handle(he_h);
		const FH& f_1 = mesh_.opposite_face_handle(he_h);

		int min_id = std::min(seg_id_[f_0.idx()], seg_id_[f_1.idx()]);
		int max_id = std::max(seg_id_[f_0.idx()], seg_id_[f_1.idx()]);

		std::pair<int, int> temp_pair = { min_id, max_id };

		if (patch_id_to_seam_id.find(temp_pair) != patch_id_to_seam_id.end())
		{
			patch_id_to_seam_id[temp_pair].push_back(path_id);
		}
		else
		{
			patch_id_to_seam_id[temp_pair] = { path_id };
		}
	}
}

bool TopoChange::Dijkstra(
	const int& start_v,
	const std::vector<int>& end_v,
	const std::vector<bool>& v_status,
	const Eigen::VectorXd& e_weight,
	double& min_dis,
	std::vector<HEH>& dij_path)
{
	min_dis = INFINITY;

	std::vector<bool> end_v_status(mesh_.n_vertices(), false);
	for (const int& v_id:end_v)
	{
		end_v_status[v_id] = true;
	}

	std::vector<double> dis(mesh_.n_vertices(), INFINITY);
	std::vector<int> pre(mesh_.n_vertices(), -1);
	dis[start_v] = 0;
	std::vector<int> Q;
	std::vector<int> U;
	std::vector<bool> Qf(mesh_.n_vertices(), false);
	std::vector<bool> Uf(mesh_.n_vertices(), false);
	Q.push_back(start_v);
	Qf[start_v] = true;
	int Index;
	int end = -1;

	while (!Q.empty())
	{
		double min = INFINITY;
		for (int i = 0; i < Q.size(); i++)
		{
			if (dis[Q[i]] < min)
			{
				min = dis[Q[i]];
				Index = i;
			}
		}

		int pointIndex = Q[Index];
		Q[Index] = Q[Q.size() - 1];
		Q.pop_back();
		U.push_back(pointIndex);
		Uf[pointIndex] = true;

		for (const HEH& adj_he : mesh_.voh_range(mesh_.vertex_handle(pointIndex)))
		{
			const VH& to_v = mesh_.to_vertex_handle(adj_he);

			if (!Uf[to_v.idx()])
			{
				if (end_v_status[pointIndex] || 
					v_status[to_v.idx()])
				{
					//std::cout << to_v.idx() << std::endl;
					double distance = min + e_weight[mesh_.edge_handle(adj_he).idx()];
					if (dis[to_v.idx()] > distance)
					{
						dis[to_v.idx()] = distance;
						pre[to_v.idx()] = pointIndex;
						if (!Qf[to_v.idx()])
						{
							Q.push_back(to_v.idx());
							Qf[to_v.idx()] = true;
						}
					}
				}

			}
		}

		//std::cout <<pointIndex << "  " << end_v_status[pointIndex] << std::endl;

		if (end_v_status[pointIndex])
		{
			end = pointIndex;
			break;
		}
	}

	if (end != -1)
	{
		int temp = end;
		dij_path.clear();
		min_dis = 0;
		while (temp != start_v)
		{
			auto heh = mesh_.find_halfedge(mesh_.vertex_handle(temp), mesh_.vertex_handle(pre[temp]));
			dij_path.emplace_back(heh);
			min_dis += e_weight[mesh_.edge_handle(heh).idx()];
			temp = pre[temp];
		}
	}
	else
	{
		return false;
	}

	return true;
}

bool TopoChange::RandomPatchInherit(TopoChange& better_topo,
	const Eigen::VectorXd& e_l,
	const Eigen::VectorXd& e_w)
{
	const std::vector<bool>& better_seam = better_topo.GetSeam();
	const std::vector<int>& better_seg = better_topo.GetSegId();

	std::vector<int> seg_id(mesh_.n_faces(), -1);
	int seg_num = 0;
	int cont = 0;
	std::vector<int> f_stack;
	for (int f_id = 0; f_id < mesh_.n_faces(); f_id++)
	{
		if (seg_id[f_id] != -1) continue;

		f_stack.clear();
		f_stack.push_back(f_id);
		do
		{
			int top_idx = f_stack.back(); f_stack.pop_back();
			seg_id[top_idx] = seg_num;

			//std::cout << cont << std::endl;
			//cont++;

			for (const HEH& fh_h : mesh_.fh_range(FH(top_idx)))
			{
				if (better_seam[mesh_.edge_handle(fh_h).idx()]
					|| seam_status_[mesh_.edge_handle(fh_h).idx()]) continue;

				FH oppo_f_h = mesh_.opposite_face_handle(fh_h);

				if (!oppo_f_h.is_valid() || seg_id[oppo_f_h.idx()] != -1) continue;

				f_stack.push_back(oppo_f_h.idx());
			}
		} while (f_stack.size() != 0);

		seg_num++;
	}

	std::vector<int> ori_v_cont(mesh_.n_vertices(), 0);
	for (const EH& e_h : mesh_.edges())
	{
		if (!seam_status_[e_h.idx()]) continue;

		HEH he_h = mesh_.halfedge_handle(e_h, 0);

		VH to_v = mesh_.to_vertex_handle(he_h);
		VH from_v = mesh_.from_vertex_handle(he_h);

		ori_v_cont[to_v.idx()]++;
		ori_v_cont[from_v.idx()]++;
	}

	std::map<std::pair<int, int>, int> adj_num;
	std::vector<int> new_num(seg_num, 0);
	for (const HEH& he_h : mesh_.halfedges())
	{
		if (mesh_.is_boundary(he_h)) continue;
		
		EH cur_e = mesh_.edge_handle(he_h);

		if (!better_seam[cur_e.idx()] &&
			!seam_status_[cur_e.idx()]) continue;

		HEH next_he = mesh_.opposite_halfedge_handle(he_h);
		int cont = 0;
		do
		{
			next_he = mesh_.next_halfedge_handle(
				mesh_.opposite_halfedge_handle(next_he));
			cont++;

			if (cont > 100)
			{
				std::cout << mesh_.face_handle(he_h).idx() << std::endl;
				std::cout << mesh_.face_handle(next_he).idx() << std::endl;
				return false;
			}
		} while (!better_seam[mesh_.edge_handle(next_he).idx()]
			&& !seam_status_[mesh_.edge_handle(next_he).idx()]);

		VH to_v = mesh_.to_vertex_handle(he_h);
		EH next_e = mesh_.edge_handle(next_he);
		FH cur_f = mesh_.face_handle(he_h);
		FH cur_oppo_f = mesh_.opposite_face_handle(he_h);
		FH next_oppo_f = mesh_.opposite_face_handle(next_he);

		int cur_sub_id = seg_id[cur_f.idx()];
		std::pair<int, int> temp_pair =
		{ seg_id[cur_f.idx()], seg_id_[cur_oppo_f.idx()] };
		if (!seam_status_[cur_e.idx()] && seam_status_[next_e.idx()])
		{
			new_num[cur_sub_id]++;
		}
		else if (seg_id_[cur_oppo_f.idx()] != seg_id_[next_oppo_f.idx()])
		{
			if (adj_num.find(temp_pair) != adj_num.end())
			{
				if (adj_num[temp_pair] != -1)
				{
					adj_num[temp_pair] += 1;
				}
			}
			else
			{
				adj_num.insert({ temp_pair, 1 });
			}
		}
		else if (seg_id_[cur_oppo_f.idx()] == seg_id_[next_oppo_f.idx()] 
			&& (ori_v_cont[to_v.idx()] != 2 || cone_status_[to_v.idx()]))
		{
			if (adj_num.find(temp_pair) != adj_num.end())
			{
				adj_num[temp_pair] = -1;
			}
			else
			{
				adj_num.insert({ temp_pair, -1 });
			}
		}
	}

	Eigen::VectorXd e_weight = e_l.cwiseProduct(e_w);
	std::map<std::pair<int, int>, double> sub_ori_length; // (sub_id, ori_id) ---> boundary length
	std::vector<double> sub_new_length(seg_num, 0);
	for (const HEH& he_h : mesh_.halfedges())
	{
		EH cur_e = mesh_.edge_handle(he_h);
		if (mesh_.is_boundary(cur_e)) continue;

		if (!better_seam[cur_e.idx()]
			&& !seam_status_[cur_e.idx()]) continue;

		int cur_seg_id = seg_id[mesh_.face_handle(he_h).idx()];
		int adj_ori_id = seg_id_[mesh_.opposite_face_handle(he_h).idx()];

		//if (!seg_status[cur_seg_id]) continue;

		if (seam_status_[cur_e.idx()])
		{
			if (sub_ori_length.find({ cur_seg_id , adj_ori_id }) != sub_ori_length.end())
			{
				//// Have cones on the boundary
				//if (sub_ori_length[{ cur_seg_id, adj_ori_id }] < 0) continue;

				//VH to_v = mesh_.to_vertex_handle(he_h);
				//if (union_cut_cont[to_v.idx()] == 2 && cone_status_[to_v.idx()])
				//{
				//	sub_ori_length[{ cur_seg_id, adj_ori_id }] = -1;
				//}
				//else if (sub_ori_length[{ cur_seg_id, adj_ori_id }] > 0)
				//{
					sub_ori_length[{ cur_seg_id, adj_ori_id }] += e_weight[cur_e.idx()];
				//}
			}
			else
			{
				sub_ori_length.insert({ { cur_seg_id, adj_ori_id }, e_weight[cur_e.idx()] });
			}
			
		}
		else if (better_seam[mesh_.edge_handle(he_h).idx()])
		{
			sub_new_length[cur_seg_id] += e_weight[cur_e.idx()];
		}
	}

	//for (const HEH& he_h : mesh_.halfedges())
	//{
	//	if (!seam_status_[mesh_.edge_handle(he_h).idx()]) continue;

	//	HEH next_he = mesh_.opposite_halfedge_handle(he_h);
	//	int cont = 0;
	//	do
	//	{
	//		next_he = mesh_.next_halfedge_handle(
	//			mesh_.opposite_halfedge_handle(next_he));
	//		cont++;

	//		if (cont > 100)
	//		{
	//			std::cout << mesh_.face_handle(he_h).idx() << std::endl;
	//			std::cout << mesh_.face_handle(next_he).idx() << std::endl;
	//		}
	//	} while (!better_seam[mesh_.edge_handle(next_he).idx()]
	//		&& !seam_status_[mesh_.edge_handle(next_he).idx()]);

	//	if (!seam_status_[mesh_.edge_handle(next_he).idx()]) continue;

	//	//FH f_0 = mesh_.opposite_face_handle(he_h);
	//	//FH f_1 = mesh_.opposite_face_handle(next_he);
	//	//if (seg_id_[f_0.idx()] == seg_id_[f_1.idx()] 
	//	//	&& union_cut_cont[mesh_.to_vertex_handle(he_h).idx()] != 2)
	//	//{
	//	//	sub_ori_length[{ seg_id[mesh_.face_handle(he_h).idx()], seg_id_[f_0.idx()] }] = -1;
	//	//}

	//	FH f_0 = mesh_.opposite_face_handle(he_h);
	//	FH f_1 = mesh_.opposite_face_handle(next_he);

	//	std::pair<int, int> temp_pair =
	//	{ seg_id[mesh_.face_handle(he_h).idx()], seg_id_[f_0.idx()] };
	//	if (adj_num[temp_pair] != 1)
	//	{
	//		sub_ori_length[temp_pair] = -1;
	//	}
	//}

	//std::ofstream sub_seg_out("sub_seg_out.txt");
	//for (int i = 0; i < mesh_.n_faces(); i++)
	//{
	//	sub_seg_out << seg_id[i] << std::endl;
	//}
	//sub_seg_out.close();

	//std::ofstream sub_ori_out("sub_ori_out.txt");
	//for (auto sub_id : sub_ori_length)
	//{
	//	sub_ori_out << sub_id.first.first
	//		<< " " << sub_id.first.second
	//		<< " " << sub_id.second
	//		<< " " << adj_num[sub_id.first]
	//		<< " " << new_num[sub_id.first.first]
	//		<< " " << sub_id.second
	//		<< " " << sub_new_length[sub_id.first.second]
	//		<< std::endl;
	//}sub_ori_out.close();

	std::vector<std::pair<int, int>> picked_sub;
	std::vector<double> percent_len;
	double sum_length = 0;
	for (auto sub_id:sub_ori_length) 
	{
		if (adj_num[sub_id.first] == 1 && new_num[sub_id.first.first] == 1)
		{
			double temp_len = sub_id.second - sub_new_length[sub_id.first.second];

			if (temp_len > 0)
			{
				percent_len.emplace_back(temp_len);
				sum_length += temp_len;
				picked_sub.emplace_back(sub_id.first);
			}
		}
	}

	if (picked_sub.size() == 0) return false;

	std::vector<double> total_len(percent_len.size() + 1, 0);
	for (int i = 0; i < percent_len.size(); i++)
	{
		total_len[i + 1] = total_len[i]
			+ percent_len[i] / sum_length;
	}

	double random_prab = 1.0 * (rand() % N_RANDOM) / N_RANDOM;

	std::pair<int, int> remove_sub_ori = {-1, -1};
	for (int i = 0; i < percent_len.size(); i++)
	{
		if (random_prab >= total_len[i] && random_prab < total_len[i + 1])
		{
			remove_sub_ori = picked_sub[i];
			break;
		}
	}

	if (remove_sub_ori.first == -1) return false;

	std::vector<int> temp_seg_id = seg_id_;

	for (const FH& f_h:mesh_.faces())
	{
		if (seg_id[f_h.idx()] == remove_sub_ori.first)
		{
			seg_id_[f_h.idx()] = remove_sub_ori.second;
		}
	}

	IdxUpdate();
	IdxUpdate();

	if (HighGenus())
	{
		seg_id_ = temp_seg_id;
		IdxUpdate();
		IdxUpdate();
		return false;
	}

	return true;
}

void TopoChange::SetCone(const std::vector<bool>& cone)
{
	cone_status_ = cone;
}

bool TopoChange::HighGenus()
{
	std::vector<int> seg_len(seg_num_, 0);
	for (const HEH& he_h : mesh_.halfedges())
	{
		if (mesh_.is_boundary(he_h)) continue;

		if (!seam_status_[
			mesh_.edge_handle(he_h).idx()]) continue;

		seg_len[seg_id_[mesh_.face_handle(he_h).idx()]]++;
	}

	std::vector<int> seg_status(seg_num_, true);
	for (const HEH& he_h : mesh_.halfedges())
	{
		if (mesh_.is_boundary(he_h)) continue;

		if (!seam_status_[mesh_.edge_handle(he_h).idx()]) continue;

		const FH& f_h = mesh_.face_handle(he_h);
		int cur_seg = seg_id_[f_h.idx()];

		int seam_cont = 0;
		if (seg_status[cur_seg])
		{
			seg_status[cur_seg] = false;

			std::vector<HEH> he_set;
			HEH cur_he = he_h; seam_cont++;

			he_set.emplace_back(cur_he);

			HEH next_he = mesh_.next_halfedge_handle(cur_he);
			for (size_t i = 0; i < seg_len[cur_seg]; i++)
			{
				int cont = 0;
				while (!seam_status_[mesh_.edge_handle(next_he).idx()] && cont < 1000)
				{
					next_he = mesh_.next_halfedge_handle(
						mesh_.opposite_halfedge_handle(next_he));

					cont++;
				}

				cur_he = next_he;
				he_set.emplace_back(cur_he);
				if (cur_he == he_h)
				{
					break;
				}

				seam_cont++;
				next_he = mesh_.next_halfedge_handle(cur_he);

				//if (cont > 999)
				//{
				//	std::cout << "Error" << std::endl;
				//	system("pause");
				//}
			}

			if (seam_cont != seg_len[cur_seg])
			{
				//std::cout << seam_cont << "  " << seg_len[cur_seg]  << std::endl;

				//std::ofstream f_file("torus.txt");
				//for (const HEH& he_h:he_set)
				//{
				//	f_file << mesh_.face_handle(he_h).idx() << std::endl;
				//}
				//f_file.close();

				return true;
			}
		}
	}

	return false;
}

std::vector<bool> TopoChange::ReturnCone()
{
	return cone_status_;
}