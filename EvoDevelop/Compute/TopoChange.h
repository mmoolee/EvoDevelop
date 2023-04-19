#pragma once
#include "SegMesh.h"
#include <Eigen/Dense>

extern const int N_RANDOM;

class TopoChange :
	public SegMesh
{
public:
	TopoChange(const Mesh&);
	TopoChange(const Mesh& mesh, 
		const std::vector<bool>& seam,
		const std::vector<bool>& cone);

	bool RandomSeamRemove(const Eigen::VectorXd& e_l,
		const Eigen::VectorXd& e_w,
		const Eigen::VectorXd& v_K,
		bool len_fixed = true,
		bool cone_fixed = true);

	bool RandomSmallPatch(const Eigen::VectorXd& e_l,
		const Eigen::VectorXd& e_w,
		bool cone_fixed = true);

	bool RandomNarrowRemove(
		const Eigen::VectorXd& e_l,
		const Eigen::VectorXd& e_w,
		const std::vector<std::vector<int>>& f_f_adj,
		const std::vector<std::vector<int>>& v_f_adj, 
		bool cone_fixed = true);

	bool RandomOneTriangleChange(bool cone_fixed = true);

	bool RandomAdjRemove(const Eigen::VectorXd& h_ang);

	bool RandomDijkstra(const Eigen::VectorXd& e_l,
		const Eigen::VectorXd& e_w);

	bool RandomPatchInherit(TopoChange& better_topo,
		const Eigen::VectorXd& e_l,
		const Eigen::VectorXd& e_w);

	void SetCone(const std::vector<bool>& cone);
	std::vector<bool> ReturnCone();

	bool HighGenus();

private:
	void VertexAdjSeam(std::vector<int>& v_cont);

	void InnerSeamSegment(std::vector<std::vector<HEH>>& he_path);

	void EdgeAndPatchToSeam(const std::vector<std::vector<HEH>>& path_he,
		std::vector<int>& edge_seam_id,
		std::map<std::pair<int, int>, std::vector<int>>& patch_id_to_seam_id);

	bool Dijkstra(const int& start_v,
		const std::vector<int>& end_v, 
		const std::vector<bool>& v_status, 
		const Eigen::VectorXd& e_weight,
		double &min_dis,
		std::vector<HEH>& dij_path);

private:
	std::vector<bool> cone_status_;


};

