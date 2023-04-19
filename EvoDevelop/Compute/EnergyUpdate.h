#pragma once
#include "Distortion.h"
#include "Individual.h"

extern const double MAX_DISTORTION;
extern const double W_SMALL_PATCH;

class EnergyUpdate
{
public:
	EnergyUpdate(const Mesh& mesh, 
		const Eigen::VectorXd& e_l,
		const Eigen::VectorXd& e_w,
		const Eigen::VectorXd& f_area,
		const std::vector<std::vector<int>>& f_f_adj,
		const std::vector<std::vector<int>>& v_f_adj);

	void SeamUpdate(const std::vector<bool>& seam,
		SeamEnergy& energy, double& fit_energy);

	void SetBound(double bound);

private:
	void ComputeSmallPatch(const std::vector<bool>&, 
		int& seg_num, int& small_cont) const;
	void RateLength(const std::vector<bool>&, double& rate_length) const;
	void SeamSmooth(const std::vector<bool>&, double& seam_energy) const;
	void SeamDistortion(const std::vector<bool>&, DistortionMetric& dis_energy);

	void NarrowRegion(const std::vector<bool>& seam,
		int& cont);

	void UpdateSegId(const std::vector<bool>& seam);

private:
	const Mesh& mesh_;
	const Eigen::VectorXd& e_l_;
	const Eigen::VectorXd& e_w_;
	const Eigen::VectorXd& f_area_;

	const std::vector<std::vector<int>>& f_f_adj_;
	const std::vector<std::vector<int>>& v_f_adj_;
	
	std::vector<int> seg_id_;
	int seg_num_ = 0;

	Distortion distortion_;

	double sum_e_l_;
	double mesh_area;

	double DISTORTION_BOUND = 0.025;
};

