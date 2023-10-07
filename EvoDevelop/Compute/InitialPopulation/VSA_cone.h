#pragma once
#include <fstream>
#include "../MeshDefinition.h"
#include "ConesFlattening.h"
#include "Kruskal/Auxiliary.h"

class VSA_cone
{
public:
	VSA_cone(const Mesh& mesh);

	std::vector<int> GenerateCone(double cone_coeff);

	void UpdateCone(const std::vector<int>& cone_set);

	std::vector<bool> Generate(int vsa_num);

	void OutputCone(const std::string& file_name);

private:
	void SpanningTree();

	void VSA();

	void UpdateTopoByPartition();

	void CutInnerCone();

	void CutHighGenus();

	void Dijkstra(const Mesh& mesh, 
		const int& start, 
		const VectorX& vlabel1, 
		const VectorX& vlabel2, 
		std::vector<HEH>& seam,
		double& dis);

private:
	int vsa_num_ = 250;

	const Mesh& mesh_;

	std::vector<int> cone_list;

	VectorX is_vertex_cone;
	VectorX partition;
	VectorX is_edge_cone_seam;
	VectorX is_edge_seam;

	int proxy_num = 200;
};

