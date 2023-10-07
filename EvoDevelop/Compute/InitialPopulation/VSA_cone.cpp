#include "VSA_cone.h"
#ifdef WIN32
#include <Windows.h>
#endif
#include <ctime>
#include <format>

const double CURVATURE_THESHOLD = 1e-4;
const double SMALL_ARGCONE = 0.025;
const double LARGE_ARGCONE = 0.1;
const int LARGE_CONE_NUM = 10000;
const int VSA_INIT_NUM = 20;
const int VSA_ITER_NUM = 500;
//const int VSA_ITER_NUM = 250;

VSA_cone::VSA_cone(const Mesh& mesh)
	:mesh_(mesh)
{
}

std::vector<int> VSA_cone::GenerateCone(double cone_coeff)
{
	VectorX conesK;

	// Cone curvature
	ConesFlattening cone_flat(mesh_);
	cone_flat.initCoef(cone_coeff);
	cone_flat.geneCone(conesK);

	//std::cout << std::format("cones K size: {}", conesK.size()) << std::endl;

	cone_list.clear();
	for (int i = 0; i < conesK.size(); i++)
	{
		if (abs(conesK[i]) > CURVATURE_THESHOLD)
		{
			cone_list.push_back(i);
		}
	}

	//std::cout << std::format("initial cone list size: {}", cone_list.size()) << std::endl;

	// Merge the adjacent cones
	for (size_t i = 0; i < cone_list.size(); i++)
	{
		//std::cout << std::format("merge adjacent cones: {}/{}...", i, cone_list.size()) << std::endl;
		for (size_t j = i + 1; j < cone_list.size(); j++)
		{
			for (const VH& vv : mesh_.vv_range(mesh_.vertex_handle(cone_list[i])))
			{
				if (vv.idx() == cone_list[j])
				{
					cone_list[j] = cone_list[i];
				}
			}
		}
	}

	//std::cout << std::format("merged cone list size: {}", cone_list.size()) << std::endl;

	// Remove same idx
	std::sort(cone_list.begin(), cone_list.end());
	cone_list.erase(std::unique(cone_list.begin(), cone_list.end()),
		cone_list.end());

	//std::ofstream cconefile("cone_" + std::to_string(cone_coeff) + ".txt");
	//for (int i = 0; i < cone_list.size(); i++)
	//{
	//	cconefile << cone_list[i] << "  " << conesK[cone_list[i]] << std::endl;
	//}
	//cconefile.close();

	return cone_list;
}

void VSA_cone::UpdateCone(const std::vector<int>& cone_set)
{
	cone_list = cone_set;
}

std::vector<bool> VSA_cone::Generate(int vsa_num)
{
	is_vertex_cone.setConstant(mesh_.n_vertices(), -1);
	for (size_t i = 0; i < cone_list.size(); i++)
	{
		is_vertex_cone[cone_list[i]] = 1;
	}

	vsa_num_ = vsa_num;
	
	is_edge_cone_seam.setConstant(mesh_.n_edges(), -1);
	if (cone_list.size()!=0)
	{
		SpanningTree();
	}

	VSA();
	UpdateTopoByPartition();
	CutInnerCone();
	UpdateTopoByPartition();

	//std::string segname = "seg_num_000.txt";
	//std::ofstream pfile(segname);
	//for (int i = 0; i < mesh_.n_faces(); i++)
	//{
	//	pfile << partition[i] << std::endl;
	//
	//}
	//pfile.close();

	CutHighGenus();

	std::vector<bool> seam_status(mesh_.n_edges(), false);
	for (const EH& e_h : mesh_.edges())
	{
		if (mesh_.is_boundary(e_h))
		{
			seam_status[e_h.idx()] = true;
		}

		HEH he_h = mesh_.halfedge_handle(e_h, 0);
		int f_0_idx = mesh_.face_handle(he_h).idx();
		int f_1_idx = mesh_.opposite_face_handle(he_h).idx();

		if (partition[f_0_idx] != partition[f_1_idx])
		{
			seam_status[e_h.idx()] = true;
		}
	}

	return seam_status;
}

void VSA_cone::OutputCone(const std::string& file_name)
{
	// Find parents
	std::ofstream cone_cout(file_name);
	for (size_t i = 0; i < cone_list.size(); i++)
	{
		cone_cout << cone_list[i] << std::endl;
	}
	cone_cout.close();
}

void VSA_cone::SpanningTree()
{
	std::vector<int> SamplePoints;
	for (int i = 0; i < cone_list.size(); i++)
	{
		SamplePoints.push_back(cone_list[i]);
	}

	MeshCache MC(mesh_);
	Algorithm::Dijkstra_group(MC, SamplePoints);
	std::vector<int> cutVertex, cutEdge;
	Algorithm::Kruskal(MC, SamplePoints, cutVertex, cutEdge);

	// Edge on the seam
	for (size_t i = 0; i < cutEdge.size(); i++)
	{
		is_edge_cone_seam[cutEdge[i]] = 1;
	}

	//time_t now_time;
	//time(&now_time);

	//std::ofstream seamfile("cone_seam_s" + std::to_string(now_time) + ".txt");
	//for (int i = 0; i < cutEdge.size(); i++)
	//{
	//	HEH he_h = mesh_.halfedge_handle(EH(cutEdge[i]), 0);
	//	VH to_v = mesh_.to_vertex_handle(he_h);
	//	VH from_v = mesh_.from_vertex_handle(he_h);
	//	seamfile << to_v << " " << from_v << std::endl;
	//}
	//seamfile.close();

	//Sleep(1000);
}

void VSA_cone::UpdateTopoByPartition()
{
	is_edge_seam.setConstant(mesh_.n_edges(), -1);
	for (auto eh : mesh_.edges())
	{
		auto heh = mesh_.halfedge_handle(eh, 0);
		if (!mesh_.is_boundary(eh))
		{
			int f1 = mesh_.face_handle(heh).idx();
			int f2 = mesh_.opposite_face_handle(heh).idx();
			if (partition[f1] != partition[f2])
			{
				is_edge_seam[eh.idx()] = 1;
			}
		}
		else
		{
			is_edge_seam[eh.idx()] = 1;
		}
	}
}

void VSA_cone::CutInnerCone()
{
	struct en
	{
		int fid;
		double cost;
		int pid;
	};

	struct cmp
	{
		bool operator()(en a, en b)
		{
			return a.cost > b.cost;
		}
	};

	// Check cone
	std::vector<int> proxy;
	int ccount = 0;
	for (size_t i = 0; i < mesh_.n_vertices(); i++)
	{
		if (is_vertex_cone[i] > 0)
		{
			ccount++;
			auto vh = mesh_.vertex_handle(i);
			int fp = -1;
			bool isconein = true;
			for (auto vf : mesh_.vf_range(vh))
			{
				if (fp < 0)
				{
					fp = partition[vf.idx()];
				}
				if (fp > -1 && fp != partition[vf.idx()])
				{
					isconein = false;
					break;
				}
			}
			if (isconein)
			{
				// Find cone seam edge
				//std::cout << i << std::endl;
				int eid = -1;
				//std::cout << " v " << vh.idx() << std::endl;
				for (auto ve : mesh_.ve_range(vh))
				{
					//std::cout << " dk " << ve.idx()
					//	<< "  " << is_edge_cone_seam[ve.idx()] << std::endl;
					if (is_edge_cone_seam[ve.idx()] > 0)
					{
						eid = ve.idx();
						break;
					}
				}
				std::priority_queue<en, std::vector<en>, cmp> q;
				while (!q.empty())
				{
					q.pop();
				}

				// Update proxy
				proxy.clear();
				auto heh = mesh_.halfedge_handle(mesh_.edge_handle(eid), 0);

				//std::cout << " dkkk " << eid << std::endl;

				if (mesh_.face_handle(heh).is_valid())
				{
					proxy.push_back(mesh_.face_handle(heh).idx());
				}
				
				if (mesh_.opposite_face_handle(heh).is_valid())
				{
					proxy.push_back(mesh_.opposite_face_handle(heh).idx());
				}
				
				for (int i = 0; i < 2; i++)
				{
					en temp;
					temp.fid = proxy[i];
					temp.pid = proxy_num + i;
					temp.cost = 0;
					q.push(temp);

					//std::cout << " proxy " << proxy[i] << std::endl;
				}

				int iter = 0;
				// Grow regions
				while (!q.empty())
				{
					//std::cout << iter << "-----" << q.size() << std::endl;
					en temp = q.top();
					q.pop();

					//std::cout << temp.fid << "  " << partition.rows() << "  " << partition.cols() << std::endl;
 
					if (partition[temp.fid] == fp)
					{
						partition[temp.fid] = temp.pid;
						for (auto ffh : mesh_.ff_range(mesh_.face_handle(temp.fid)))
						{
							if (partition[ffh.idx()] == fp)
							{
								en temp1;
								temp1.fid = ffh.idx();
								temp1.pid = temp.pid;
								temp1.cost = temp.cost + 1;

								//std::cout << " pr " << temp1.fid << "  " << temp.pid << "  " << temp1.cost << std::endl;

								q.push(temp1);
							}
						}
					}
					iter++;

				}
				//std::cout << "  split over  "
				//	<< std::endl;
				for (size_t i = 0; i < mesh_.n_faces(); i++)
				{
					if (partition[i] == proxy_num + 1)
					{
						partition[i] = fp;
					}
				}
				proxy_num++;
			}
		}


		//std::ofstream pfile(tpath + "partition_vsa.txt");
		//for (int i = 0; i < mesh.n_faces(); i++)
		//{
		//	pfile << partition[i] << std::endl;
		//}
		//pfile.close();

	}
}

void VSA_cone::CutHighGenus()
{
	struct en
	{
		int fid;
		double cost;
		int pid;
	};

	struct cmp
	{
		bool operator()(en a, en b)
		{
			return a.cost > b.cost;
		}
	};

	int k_num = 0;
	for (size_t i = 0; i < mesh_.n_faces(); i++)
	{
		k_num = fmax(k_num, partition[i]);
	}
	k_num++;

	int seg_num;
	seg_num = k_num;

	bool cut_status;
	do
	{
		k_num = seg_num;

		cut_status = false;

		VectorX seamlabel;
		VectorX seamseam;
		VectorX vlabel1;
		VectorX vlabel2;

		std::vector<int> seg_f_num(k_num, 0);
		for (const FH& f_h:mesh_.faces())
		{
			seg_f_num[partition[f_h.idx()]]++;
		}

		for (size_t p_id = 0; p_id < k_num; p_id++)
		{
			if (seg_f_num[p_id] == 0) continue;

			vlabel1.setConstant(mesh_.n_vertices(), -1);
			vlabel2.setConstant(mesh_.n_vertices(), -1);

			seamlabel.setConstant(mesh_.n_edges(), -1);
			seamseam.setConstant(mesh_.n_edges(), -1);
			int first_e = -1;
			int seamcount = 0;
			int bdcount = 1;
			for (auto eh : mesh_.edges())
			{
				auto heh1 = mesh_.halfedge_handle(eh, 0);
				auto heh2 = mesh_.opposite_halfedge_handle(heh1);
				auto fh1 = mesh_.face_handle(heh1);
				auto fh2 = mesh_.face_handle(heh2);
				if (((fh1.is_valid() && partition[fh1.idx()] == p_id)
					||(fh2.is_valid() && partition[fh2.idx()] == p_id))
					&& ((partition[fh1.idx()] != partition[fh2.idx()])))
				{
					seamcount++;
					seamlabel[eh.idx()] = 1;
					vlabel1[mesh_.from_vertex_handle(heh1).idx()] = 1;
					vlabel1[mesh_.to_vertex_handle(heh1).idx()] = 1;
					if (first_e < 0)
					{
						first_e = eh.idx();
					}
				}
			}

			if (seamcount == 0) continue;

			HEH start_he;
			for (const HEH & he_h:mesh_.halfedges())
			{
				if (mesh_.is_boundary(he_h)) continue;
				
				if (seamlabel[mesh_.edge_handle(he_h).idx()] != 1) continue;

				if (partition[mesh_.face_handle(he_h).idx()] == p_id)
				{
					start_he = he_h; break;
				}
			}

			std::vector<int> start_v_set;
			int seam_cont = 0;
			HEH cur_he = start_he;
			do
			{
				cur_he = mesh_.opposite_halfedge_handle(cur_he);

				do
				{
					cur_he = mesh_.next_halfedge_handle(
						mesh_.opposite_halfedge_handle(cur_he));
				} while (seamlabel[mesh_.edge_handle(cur_he).idx()] != 1);
				
				seam_cont++;
				start_v_set.emplace_back(mesh_.to_vertex_handle(cur_he).idx());
				vlabel2[mesh_.to_vertex_handle(cur_he).idx()] = 1;

			} while (cur_he != start_he);

			//std::cout << "cont: " << seam_cont
			//	<< " " << seamcount << std::endl;

			//std::ofstream v_file("min_start_v.txt");
			//for (int v_id:start_v_set)
			//{
			//	v_file << v_id << std::endl;
			//}
			//v_file.close();

			if (seam_cont != seamcount)
			{
				cut_status = true;

				// Cut 
				double min_dis = DBL_MAX;
				double temp_dis;
				std::vector<HEH> temp_seam;
				std::vector<HEH> min_seam;
				for (int v_id:start_v_set)
				{
					Dijkstra(mesh_, v_id, vlabel1, vlabel2, temp_seam, temp_dis);

					if (temp_dis < min_dis)
					{
						min_seam = temp_seam;
						min_dis = temp_dis;
					}
				}
			
				for (const HEH& he_h:min_seam)
				{
					seamlabel[mesh_.edge_handle(he_h).idx()] = 1;
				}

				HEH start_he = min_seam.front();

				FH f_0 = mesh_.face_handle(start_he);
				FH f_1 = mesh_.opposite_face_handle(start_he);

				std::priority_queue<en, std::vector<en>, cmp> q;
				while (!q.empty())
				{
					q.pop();
				}

				en temp;
				temp.cost = 0;

				temp.fid = f_0.idx();
				temp.pid = seg_num;
				q.push(temp);

				temp.fid = f_1.idx();
				temp.pid = seg_num + 1;
				q.push(temp);

				// Grow regions
				while (!q.empty())
				{
					en temp = q.top();
					q.pop();

					if (partition[temp.fid] == p_id)
					{
						partition[temp.fid] = temp.pid;
						for (auto fh : mesh_.fh_range(mesh_.face_handle(temp.fid)))
						{
							if (seamlabel[mesh_.edge_handle(fh).idx()] < 0
								&& partition[mesh_.opposite_face_handle(fh).idx()] == p_id)
							{
								en temp1;
								temp1.fid = mesh_.opposite_face_handle(fh).idx();
								temp1.pid = temp.pid;
								temp1.cost = temp.cost + 1;

								q.push(temp1);
							}
						}
					}
				}

				seg_num += 2;
			}
		}

	} while (cut_status);
}

void VSA_cone::Dijkstra(const Mesh& mesh, 
	const int& start,
	const VectorX& vlabel1,
	const VectorX& vlabel2,
	std::vector<HEH>& seam,
	double& min_dis)
{
	min_dis = INFINITY;
	seam.clear();

	int vertices_num = mesh.n_vertices();
	std::vector<double> dis(vertices_num, INFINITY);
	std::vector<int> pre(vertices_num, -1); 
	dis[start] = 0;
	std::vector<int> Q;
	std::vector<int> U;
	std::vector<bool> Qf(vertices_num, false);
	std::vector<bool> Uf(vertices_num, false);
	Q.push_back(start);
	Qf[start] = true;
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

		for (const auto vh : mesh.vv_range(mesh.vertex_handle(pointIndex)))
		{
			if (!Uf[vh.idx()])
			{
				if (vlabel1[vh.idx()] < 0 || (vlabel1[vh.idx()] > 0 && vlabel2[vh.idx()] < 0))
				{
					double distance = min + (mesh.point(mesh.vertex_handle(vh.idx()))
						- mesh.point(mesh.vertex_handle(pointIndex))).norm();
					if (dis[vh.idx()] > distance)
					{
						dis[vh.idx()] = distance;
						pre[vh.idx()] = pointIndex;
						if (!Qf[vh.idx()])
						{
							Q.push_back(vh.idx());
							Qf[vh.idx()] = true;
						}
					}
				}

			}
		}
		if (vlabel1[pointIndex] > 0 && vlabel2[pointIndex] < 0)
		{
			end = pointIndex;
			break;
		}
	}

	int temp = end;

	if (temp != -1)
	{
		min_dis = 0;
		while (temp != start)
		{
			auto heh = mesh.find_halfedge(mesh.vertex_handle(temp), mesh.vertex_handle(pre[temp]));
			seam.emplace_back(heh);
			min_dis += mesh.calc_edge_length(heh);
			temp = pre[temp];
		}
	}
}

//void VSA_cone::OutputCone(const std::vector<int>& cone_list,
//	const std::string& file_name)
//{
//	std::ofstream conefile(file_name);
//	for (const int& i:cone_list)
//	{
//		conefile << i << std::endl;
//	}
//	conefile.close();
//}

void VSA_cone::VSA()
{
	struct Tri
	{
		int id;
		int tag;
		double weight;

		Tri(int id, int tag, double weight) :tag(tag), id(id), weight(weight) {};

		bool operator<(const Tri& a) const
		{
			return weight > a.weight;
		}
	};

	int k_num = VSA_INIT_NUM;
	int f_num = mesh_.n_faces();

	std::vector<OpenMesh::Vec3d> plane;
	std::vector<OpenMesh::Vec3d> face_normal(f_num);

	// Area
	std::vector<double> face_area(f_num);
	double total_area = 0;
	for (auto f : mesh_.faces())
	{
		face_normal[f.idx()] = mesh_.calc_face_normal(f);
		auto he = mesh_.halfedge_handle(f);
		auto p0 = mesh_.point(mesh_.from_vertex_handle(he));
		auto p1 = mesh_.point(mesh_.to_vertex_handle(he));
		auto p2 = mesh_.point(mesh_.to_vertex_handle(mesh_.next_halfedge_handle(he)));
		face_area[f.idx()] = fabs(((p0 - p1) % (p2 - p1)).norm());
		total_area += face_area[f.idx()];
	}
	total_area = sqrt(total_area);
	//std::cout << "total area" << total_area << std::endl;

	plane.resize(k_num);
	partition.setConstant(f_num, -1);

	std::vector<int> seed_set(k_num);
	for (size_t i = 0; i < k_num; i++)
	{
		int r_id = i * f_num / k_num;
		partition[r_id] = i;
		plane[i] = face_normal[r_id];
		seed_set[i] = 0;
	}

	for (int itnum = 0; ; itnum++)
	{
		std::vector<double> small_energy(k_num, DBL_MAX);
		std::vector<double> max_energy(k_num, DBL_EPSILON);
		std::vector<int> max_seed(k_num);
		double total_energy = 0;
		std::priority_queue<Tri> que_add;
		int max_id = -1;
		double max_energy_d = DBL_EPSILON;

		// Find proxy for each patch
		for (int i = 0; i < f_num; i++)
		{
			if (partition[i] != -1)
			{
				double energy = (plane[partition[i]] - face_normal[i]).sqrnorm();
				//double energy = (plane[partition[i]] - face_normal[i]).norm();
				total_energy += energy * face_area[i];
				if (energy < small_energy[partition[i]])
				{
					small_energy[partition[i]] = energy;
					seed_set[partition[i]] = i;
				}
				if (energy > max_energy[partition[i]])
				{
					max_energy[partition[i]] = energy;
					max_seed[partition[i]] = i;
					que_add.push(Tri(i, 0, -energy));
				}
				if (energy > max_energy_d)
				{
					max_energy_d = energy;
					max_id = i;
				}
			}
		}

		// Find the unassigned triangle 
		bool is_valid = true;
		for (int i = 0; i < f_num; i++)
		{
			if (partition[i] < 0)
			{
				seed_set.push_back(i);
				partition[i] = k_num;
				plane.push_back(face_normal[i]);
				k_num++;
				is_valid = false;
				break;
			}
		}

		if (is_valid)
		{
			// Add new plane
			int k_num_temp = k_num;
			if (k_num < vsa_num_ && true && !que_add.empty())
			{
				int id = que_add.top().id;
				//std::cout << "max id " << id << std::endl;
				if (id > -1)
				{
					seed_set.push_back(id);
					partition[id] = k_num_temp;
					plane.push_back(face_normal[id]);
					k_num_temp++;
				}
			}
			k_num = k_num_temp;
		}

		//if (itnum > 0) {
		//	std::cout << "Iteration " << itnum
		//		<< " plane number " << k_num
		//		<< " max angle " << -que_add.top().weight
		//		<< " Total energy : " << total_energy << std::endl;
		//}

		if (/*itnum > 1 && abs(total_energy - last_energy) / total_energy < 1e-7 ||*/ itnum > VSA_ITER_NUM) break;
		
		double last_energy = total_energy;
		
		std::vector<bool> is_conq(f_num, false);

		std::priority_queue<Tri> que;
		for (int i = 0; i < k_num; i++)
		{
			partition[seed_set[i]] = i;
			is_conq[seed_set[i]] = true;
			for (auto fheh : mesh_.fh_range(mesh_.face_handle(seed_set[i])))
			{
				auto eh = mesh_.edge_handle(fheh);
				if (!mesh_.is_boundary(eh))
				{
					auto fh = mesh_.opposite_face_handle(fheh);
					if (is_edge_cone_seam[eh.idx()] < 0)
					{
						double energy = (face_normal[fh.idx()] - plane[i]).norm();
						que.push(Tri(fh.idx(), i, energy));
					}
				}
			}
		}

		while (que.size() != 0)
		{
			int fid = que.top().id;
			int pid = que.top().tag;
			que.pop();
			if (!is_conq[fid])
			{
				partition[fid] = pid;
				is_conq[fid] = true;
				for (auto fheh : mesh_.fh_range(mesh_.face_handle(fid)))
				{
					auto eh = mesh_.edge_handle(fheh);
					if (!mesh_.is_boundary(eh))
					{
						auto fh = mesh_.opposite_face_handle(fheh);
						if (is_edge_cone_seam[eh.idx()] < 0 && !is_conq[fh.idx()])
						{
							double energy = (face_normal[fh.idx()] - plane[pid]).norm();
							que.push(Tri(fh.idx(), pid, energy));
						}
					}
				}
			}
		}

		std::vector<OpenMesh::Vec3d> new_normal(k_num, OpenMesh::Vec3d(0, 0, 0));
		VectorX testlm;
		testlm.setConstant(k_num, -1);
		int wc = 0;
		for (int j = 0; j < f_num; j++)
		{
			if (partition[j] > -1)
			{
				new_normal[partition[j]] += face_normal[j] * face_area[j];
			}
		}
		//std::cout << k_num << " -- " << new_normal.size() << " -- " << plane.size() << std::endl;
		for (int i = 0; i < k_num; i++)
		{
			new_normal[i].normalize();
			plane[i] = new_normal[i];
		}
	}

	//std::ofstream pfile("partition_vsa.txt");
	//for (int i = 0; i < mesh_.n_faces(); i++)
	//{
	//	pfile << partition[i] << std::endl;
	//}
	//pfile.close();
}

