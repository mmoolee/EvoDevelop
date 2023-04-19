#pragma once
#include <fstream>
#include "../MeshDefinition.h"

//typedef std::map<std::pair<int, int>, std::vector<EH>> Idx_Boundary;

class SegMesh
{
public:
	SegMesh(const Mesh& mesh);
	SegMesh(const Mesh& mesh, std::vector<bool> seam);
	SegMesh(const Mesh& mesh, std::vector<int> seg_id);

	void Init();

	void ReadSeg(const std::string&);

	const Mesh& GetMesh();

	void SetSeam(const std::vector<bool>&);
	std::vector<bool>& GetSeam();
	void SeamUpdate(); // Update when seams change
	std::vector<bool> ReturnSeam() const;

	void SetIdx(std::vector<int>);
	std::vector<int>& GetSegId();
	void IdxUpdate(); // Update when ids change

	int& GetSegNum();

	void OutputSeg(const std::string&) const;

	void OutputMesh(const std::string&);

	// Face number of each patch
	void PatchFaceNum(std::vector<int>&);

protected:
	const Mesh& mesh_;

	int seg_num_;
	std::vector<bool> seam_status_;
	std::vector<int> seg_id_;
};

