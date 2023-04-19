#pragma once
#include "Distortion.h"

struct SeamEnergy
{
	double rate_length_ = DBL_MAX;
	double seam_smooth_ = DBL_MAX;
	int seg_num_ = 0;

	// Status
	int small_cont_ = INT_MAX;
	int narrow_cont_ = INT_MAX;

	DistortionMetric dis_metric_;

	SeamEnergy& operator=(const SeamEnergy& other)
	{
		this->rate_length_ = other.rate_length_;
		this->dis_metric_ = other.dis_metric_;
		this->small_cont_ = other.small_cont_;
		this->seam_smooth_ = other.seam_smooth_;
		this->seg_num_ = other.seg_num_;
		this->narrow_cont_ = other.narrow_cont_;

		return *this;
	}
};

class Individual
{
public:
	Individual();

	Individual& operator=(const Individual& other_ind);

	bool operator==(const Individual& rhs);

public:
	// Cone
	std::vector<bool> cone_status_;

	// Seam
	std::vector<bool> seam_status_;
	SeamEnergy seam_energy_;
	double fitness_energy_ = 0;
};