#include "SegMesh.h"

SegMesh::SegMesh(const Mesh& mesh)
	:mesh_(mesh)
{
	Init();
}

SegMesh::SegMesh(const Mesh& mesh, std::vector<bool> seam)
	:mesh_(mesh)
{
	SetSeam(seam);
	SeamUpdate();
}

SegMesh::SegMesh(const Mesh& mesh, std::vector<int> seg_id)
	:mesh_(mesh)
{
	SetIdx(seg_id);
	IdxUpdate();
}

void SegMesh::ReadSeg(const std::string& file_name)
{
	std::ifstream seg_file(file_name);
	for (size_t i = 0; i < mesh_.n_faces(); i++)
	{
		FH f_h = mesh_.face_handle(i);
		int temp_id;
		seg_file >> temp_id;

		seg_id_[i] = temp_id;
	}
	seg_file.close();

	IdxUpdate();
}

void SegMesh::Init()
{
	seg_num_ = 1;

	seam_status_.assign(mesh_.n_edges(), false);
	
	// Boundary
	for (const EH& e_h:mesh_.edges())
	{
		if (!mesh_.is_boundary(e_h)) continue;
		
		seam_status_[e_h.idx()] = true;
	}

	seg_id_.assign(mesh_.n_faces(), 0);
}

const Mesh& SegMesh::GetMesh()
{
	return mesh_;
}

void SegMesh::SetSeam(const std::vector<bool>& seam)
{
	seam_status_ = seam;
	SeamUpdate();
}

std::vector<bool>& SegMesh::GetSeam()
{
	return seam_status_;
}

void SegMesh::SetIdx(std::vector<int> seg_id)
{
	seg_id_ = seg_id;
	IdxUpdate();
}

std::vector<int>& SegMesh::GetSegId()
{
	return seg_id_;
}

int& SegMesh::GetSegNum()
{
	return seg_num_;
}

void SegMesh::SeamUpdate()
{
	int n_f = mesh_.n_faces();

	seg_id_.assign(n_f, -1);

	seg_num_ = 0;

	std::vector<FH> f_stack;
	for (const FH& f_h:mesh_.faces())
	{
		if (seg_id_[f_h.idx()] != -1) continue;

		f_stack.clear();
		f_stack.push_back(f_h);
		do
		{
			FH top_f = f_stack.back(); f_stack.pop_back();
			seg_id_[top_f.idx()] = seg_num_;

			for (const HEH& fh_h : mesh_.fh_range(top_f))
			{
				if (seam_status_[mesh_.edge_handle(fh_h).idx()]) continue;

				FH oppo_f_h = mesh_.opposite_face_handle(fh_h);

				if (!oppo_f_h.is_valid()) continue;
				if (seg_id_[oppo_f_h.idx()] != -1) continue;

				f_stack.push_back(oppo_f_h);
			}
		} while (f_stack.size() != 0);

		++seg_num_;
	}
}

std::vector<bool> SegMesh::ReturnSeam() const
{
	return seam_status_;
}

void SegMesh::IdxUpdate()
{
	seam_status_.assign(mesh_.n_edges(), false);

	int idx_0, idx_1;
	for (const EH& e_h : mesh_.edges())
	{
		if (mesh_.is_boundary(e_h))
		{
			seam_status_[e_h.idx()] = true;
			continue;
		}

		idx_0 = mesh_.face_handle(mesh_.halfedge_handle(e_h, 0)).idx();
		idx_1 = mesh_.face_handle(mesh_.halfedge_handle(e_h, 1)).idx();

		if (seg_id_[idx_0] == seg_id_[idx_1]) continue;

		seam_status_[e_h.idx()] = true;
	}

	SeamUpdate();
}

void SegMesh::OutputSeg(const std::string& file_name) const
{
	std::ofstream seg_cout(file_name);
	for (size_t i = 0; i < mesh_.n_faces(); i++)
	{
		if (i != mesh_.n_faces() - 1)
		{
			seg_cout << seg_id_[i] << '\n';
		}
		else
		{
			seg_cout << seg_id_[i];
		}

	}
	seg_cout.close();
}

void SegMesh::OutputMesh(const std::string& file_name)
{
	MeshTools::WriteMesh(mesh_, file_name);
}

void SegMesh::PatchFaceNum(std::vector<int>& f_num)
{
	f_num.assign(seg_num_, 0);
	for (const int& temp_id:seg_id_)
	{
		++f_num[temp_id];
	}
}