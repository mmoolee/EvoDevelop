#pragma once
#include "ConesFlattening.h"
#include <iostream>
#include <fstream>
#include <vector>
void Dijkstra(Mesh mesh, int start, VectorX &vlabel1, VectorX &vlabel2, VectorX &seam)
{
	int vertices_num = mesh.n_vertices();
	std::vector<double> dis(vertices_num, INFINITY); //记录最短路径长度，初始都为无穷大
	std::vector<int> pre(vertices_num, -1);  //记录每个点在最短路径上的前一个点
	dis[start] = 0;
	std::vector<int> Q;       //距离dis不为INFINITY，且不知道最短路径的点
	std::vector<int> U;       //已知最短路径的点
	std::vector<bool> Qf(vertices_num, false);    //判断是否在Q和U中
	std::vector<bool> Uf(vertices_num, false);
	Q.push_back(start);
	Qf[start] = true;
	int Index;
	int end;
	while (!Q.empty())
	{
		double min = INFINITY;
		for (int i = 0; i < Q.size(); i++)   //找到Q中dis最小的点的序号pointIndex，并加入到U中
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
		for (const auto vh : mesh.vv_range(mesh.vertex_handle(pointIndex)))     //更新pointIndex的1-邻域的点的dis值
		{
			if (!Uf[vh.idx()])
			{
				if (vlabel1[vh.idx()] < 0 || (vlabel1[vh.idx()] > 0 && vlabel2[vh.idx()] > 0))
				{
					double distance = min + (mesh.point(mesh.vertex_handle(vh.idx())) - mesh.point(mesh.vertex_handle(pointIndex))).norm();
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
		if (vlabel1[pointIndex] > 0 && vlabel2[pointIndex] > 0)   //当end已经在U中时，停止循环
		{
			end = pointIndex;
			break;
		}
	}
	int temp = end;
	vlabel2[end] = -1;
	while (temp != start)
	{
		auto heh = mesh.find_halfedge(mesh.vertex_handle(temp), mesh.vertex_handle(pre[temp]));
		seam[mesh.edge_handle(heh).idx()] = 1;
		vlabel1[temp] = 1;
		temp = pre[temp];
	}
}


void highgenus()
{
	// read mesh
	std::string objpath = "obj.obj";
	std::string parpath = "seg.txt";
	Mesh mesh;
	std::cout << "load mesh from " << objpath << std::endl;
	if (!MeshTools::ReadMesh(mesh, objpath))
	{
		std::cout << "load failed!\n";
		exit(EXIT_FAILURE);
	}
	// read segment
	VectorX partition;
	int k_num = 0;
	partition.resize(mesh.n_faces());
	std::ifstream seg_file(parpath);
	for (size_t i = 0; i < mesh.n_faces(); i++)
	{
		int p_num;
		seg_file >> p_num;
		k_num = fmax(k_num, p_num);
		partition[i] = p_num;
	}
	seg_file.close();
	k_num++;
	int seg_num;
	seg_num = k_num;

	VectorX seamlabel;
	VectorX seamseam;
	VectorX vlabel1;
	VectorX vlabel2;
	vlabel1.setConstant(mesh.n_vertices(), -1);
	vlabel2.setConstant(mesh.n_vertices(), -1);
	for (size_t p_id = 0; p_id < k_num; p_id++)
	{
		seamlabel.setConstant(mesh.n_edges(), -1);
		seamseam.setConstant(mesh.n_edges(), -1);
		int first_e = -1;
		int seamcount = 0;
		int bdcount = 1;
		for (auto eh : mesh.edges())
		{
			auto heh1 = mesh.halfedge_handle(eh, 0);
			auto heh2 = mesh.opposite_halfedge_handle(heh1);
			auto fh1 = mesh.face_handle(heh1);
			auto fh2 = mesh.face_handle(heh2);
			if ((partition[fh1.idx()] != partition[fh2.idx()]) && (partition[fh1.idx()] == p_id || partition[fh2.idx()] == p_id))
			{
				seamcount++;
				seamlabel[eh.idx()] = 1;
				vlabel1[mesh.from_vertex_handle(heh1).idx()] = 1;
				vlabel1[mesh.to_vertex_handle(heh1).idx()] = 1;
				if (first_e < 0)
				{
					first_e = eh.idx();
				}
			}
		}
		seamseam[first_e] = 1;
		int iter_e = first_e;
		int iter_he = mesh.halfedge_handle(mesh.edge_handle(iter_e), 0).idx();
		for (size_t iter = 0; iter < seamcount; iter++)
		{
			auto vh = mesh.to_vertex_handle(mesh.halfedge_handle(iter_he));
			for (auto vhe : mesh.voh_range(vh))
			{
				auto ve = mesh.edge_handle(vhe);
				if (seamseam[ve.idx()] < 0 && seamlabel[ve.idx()] > 0)
				{
					bdcount++;
					seamseam[ve.idx()] = 1;
					vlabel2[mesh.from_vertex_handle(vhe).idx()] = 1;
					vlabel2[mesh.to_vertex_handle(vhe).idx()] = 1;
					iter_he = vhe.idx();
					break;
				}
			}
		}
		if (bdcount != seamcount)
		{
			std::cout << " high genus, patch : " << p_id << std::endl;
			std::string segname = "seg_num" + std::to_string(p_id) + ".txt";
			std::ofstream pfile(segname);
			for (int i = 0; i < mesh.n_faces(); i++)
			{
				if (partition[i] == p_id)
				{
					pfile << 1 << std::endl;
				}
				else
				{
					pfile << 0 << std::endl;
				}
			}
			pfile.close();

			// cut 
			// find two init vertex
			std::vector<int> init_list;
			init_list.clear();
			for (size_t i = 0; i < mesh.n_vertices(); i++)
			{
				if (vlabel1[i] > 0 && vlabel2[i] < 0)
				{
					init_list.push_back(i);
				}
			}
			double max_dis = 0;
			int init_point2 = -1;
			for (size_t i = 1; i < init_list.size(); i++)
			{
				double dis = (mesh.point(mesh.vertex_handle(init_list[0])) - mesh.point(mesh.vertex_handle(init_list[i]))).norm();
				if (dis > max_dis)
				{
					max_dis = dis;
					init_point2 = init_list[i];
				}
			}

			std::string seamname = "seam" + std::to_string(p_id) + "-1-s.txt";
			std::ofstream sfile(seamname);
			for (int i = 0; i < mesh.n_edges(); i++)
			{
				sfile << seamlabel[i] << std::endl;
			}
			sfile.close();
			Dijkstra(mesh, init_list[0], vlabel1, vlabel2, seamlabel);
			std::string seamname1 = "seam" + std::to_string(p_id) + "-2-s.txt";
			std::ofstream sfile1(seamname1);
			for (int i = 0; i < mesh.n_edges(); i++)
			{
				sfile1 << seamlabel[i] << std::endl;
			}
			sfile1.close();
			Dijkstra(mesh, init_point2, vlabel1, vlabel2, seamlabel);
			std::string seamname3 = "seam" + std::to_string(p_id) + "-3-s.txt";
			std::ofstream sfile3(seamname3);
			for (int i = 0; i < mesh.n_edges(); i++)
			{
				sfile3 << seamlabel[i] << std::endl;
			}
			sfile3.close();

			// update partition
			int seed = -1;
			for (size_t i = 0; i < mesh.n_faces(); i++)
			{
				if (partition[i] == p_id)
				{
					seed = i;
					break;
				}
			}
			std::queue<int> pq;
			pq.push(seed);
			while (!pq.empty())
			{
				int temp = pq.front();
				pq.pop();
				partition[temp] = seg_num;
				for (auto fhe : mesh.fh_range(mesh.face_handle(temp)))
				{
					if (seamlabel[mesh.edge_handle(fhe).idx()] < 0 && partition[mesh.opposite_face_handle(fhe).idx()]==p_id)
					{
						pq.push(mesh.opposite_face_handle(fhe).idx());
					}
				}
			}
			seg_num++;

			std::string segnamen = "seg_num_new" + std::to_string(p_id) + ".txt";
			std::ofstream pfilen(segnamen);
			for (int i = 0; i < mesh.n_faces(); i++)
			{
				pfilen << partition[i] << std::endl;
			}
			pfilen.close();
		}
		/* code */
	}
}


