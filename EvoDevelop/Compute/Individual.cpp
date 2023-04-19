#include "Individual.h"

Individual::Individual()
{
}

Individual& Individual::operator=(const Individual& other_ind)
{
	this->seam_status_ = other_ind.seam_status_;
	this->cone_status_ = other_ind.cone_status_;

	this->seam_energy_ = other_ind.seam_energy_;
	this->fitness_energy_ = other_ind.fitness_energy_;

	return *this;
}

bool Individual::operator==(const Individual& rhs)
{
	if (this->seam_status_.size() != rhs.seam_status_.size())
		return false;

	for (size_t i = 0; i < this->seam_status_.size(); i++)
	{
		if (this->seam_status_[i] != rhs.seam_status_[i])
		{
			return false;
		}
	}

	return true;
}
