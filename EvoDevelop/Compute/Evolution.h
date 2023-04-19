#pragma once
#include "Individual.h"
#include "TopoChange.h"
#include "InitialPopulation\VSA_cone.h"
#include "EnergyUpdate.h"

// Segment
class Evolution
{
public:
	Evolution(Mesh mesh);

	void DebugPath(const std::string& path = "debug_output");

	~Evolution();

	void Run(double bound = 0.025, int cone_size = 300);

private:
	void VertGauss();
	void FaceAngle();
	void FaceArea();
	void EdgeLength();
	void EdgeWeight();
	void InitAdjRegion();
	void FaceAdjRegion(const FH& f_h,
		const Eigen::MatrixX3d& f_c_,
		std::vector<int>& adj_f_ids);

	void InitialPopulation();
	void ComputeFixedCone();
	void InitialArchive();
	void UpdateArchive();

	void StochasticSelection();
	void RandomMutation(const int& cont);
	void Crossover(const int& cont);
	bool ElitistReinsert(const int& cont);

	void ClearVector(std::vector<Individual*>& vec);
	void PopVector(std::vector<Individual*>& vec);
	
	void OutputInitial();
	void OutputArchive();
	void OutputSelection();
	void OutputMutation();
	void OutputCrossover();
	void OutputPopulation();

private:
	Mesh mesh_; // Normalized mesh;

	double DISTORTION_BOUND = 0.025;
	//double DISTORTION_BOUND = 0.005;

	int N_CONE_THRESHOLD = 300;
	//int N_CONE_THRESHOLD = 150;

	std::vector<EnergyUpdate*> energy_array_;

	std::vector<Individual*> initial_;
	std::vector<Individual*> archive_;

	std::vector<TopoChange> pop_seg_;
	std::vector<int> parent_;

	//std::vector<std::pair<Individual*, Individual*>> island_;

	enum change_type
	{
		C_NON, C_SEAM, C_SMALL, C_NARROW, C_DIJ, C_TRI, C_ADJ, C_CROSS
	};
	std::vector<change_type> change_type_;

	std::vector<Individual*> population_;

	Eigen::VectorXd v_K_;
	Eigen::VectorXd f_ang_;
	Eigen::VectorXd f_area_;
	Eigen::VectorXd e_l_;
	Eigen::VectorXd e_w_;
	Eigen::MatrixX3d f_c_;
	double region_dis_;
	std::vector<std::vector<int>> f_f_adj_;
	std::vector<std::vector<int>> v_f_adj_;

	DistortionMetric dis_bound_;

	std::string ROOT_PATH = "debug_output";

	std::string PATH_INITIAL = "\\initial";
	std::string PATH_ARCHIVE = "\\archive";
	std::string PATH_SELECTION = "\\selection";
	std::string PATH_CROSSOVER = "\\crossover";
	std::string PATH_MUTATION = "\\mutation";
	std::string PATH_POPULATION = "\\population";

	int N_TRI_REMOVE = 5;
	int INIT_MIN_NUM = 300;
	int N_GENERATION = 1000;

	const int N_VSA_COEFF = 10;
	int N_INIT_VSA = 100;
	//int N_VSA_STEP = 40; // 47094
	int N_VSA_STEP = 10;

	const int N_CORE = 6;
	const int N_CONE_COEFF = 2 * N_CORE;

	double DECREASE_STEP = 1.1;
	double N_CONE_STEP = DISTORTION_BOUND / N_CORE;
	int N_ITER_STEP = 200;

	std::vector<int> constraint_list_;
};

