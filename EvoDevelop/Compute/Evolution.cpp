#include "Evolution.h"
#include <io.h>

const int N_ARCHIVE = 50;
const int N_ARCHIVE_REPLACE = 0;

const int N_CONE_NUM = 30;

const int N_ARCHIVE_FIX = 20;
const int N_RANDOM_TRY = 5;
const int N_TRI_REMOVE_STEP = 5000;

const int N_ELITIST_NUM = 10;
const int N_RANDOM_PICK = 10;
const int N_ELITIST_PICK = 10;
const int N_CROSSOVER = 10;
const int N_MAX_ITER = 100;
//const double P_ADJ_PERCENT = 0.001;
const double P_SMALL_REGION = 0.0001;
//const double MIN_DIS_THRESHOLD = 0.0025;

const double CHANGE_THRESHOLD_0 = 0.1;
const double CHANGE_THRESHOLD_1 = 0.05;
const double CHANGE_THRESHOLD_2 = 0.01;

bool IndCmp(Individual* ind_0, Individual* ind_1)
{
	return ind_0->fitness_energy_
		< ind_1->fitness_energy_;
}

bool IndCmp2(std::pair<Individual*, Individual*> ind_0,
	std::pair<Individual*, Individual*> ind_1)
{
	return ind_0.first->fitness_energy_
		< ind_1.first->fitness_energy_;
}

Evolution::Evolution(Mesh mesh)
	:mesh_(mesh)
{
	FaceArea();
	FaceAngle();
	VertGauss();
	EdgeLength();
	EdgeWeight();
	InitAdjRegion();

	int core_num = omp_get_num_procs();
	int fit_compute_num = std::min(core_num, N_CONE_COEFF);

	std::cout << "Fit Compute Num: " << fit_compute_num << std::endl;

	for (size_t i = 0; i < fit_compute_num; i++)
	{
		EnergyUpdate* fit = new EnergyUpdate(mesh_, e_l_, e_w_, f_area_, 
			f_f_adj_, v_f_adj_);
		energy_array_.emplace_back(fit);
	}

	if (mesh_.n_faces() > N_TRI_REMOVE_STEP)
	{
		N_TRI_REMOVE *= ((int) mesh_.n_faces() / N_TRI_REMOVE_STEP);
	}
}

void Evolution::DebugPath(const std::string& path)
{
	std::cout << "Create Output Dir" << std::endl;

	ROOT_PATH = path;
	PATH_INITIAL  = path + "\\initial";
	PATH_ARCHIVE = path + "\\archive";
	PATH_SELECTION = path + "\\selection";
	PATH_CROSSOVER = path + "\\crossover";
	PATH_MUTATION = path + "\\mutation";
	PATH_POPULATION = path + "\\population";


	if (_access(PATH_INITIAL.c_str(), 0) == -1)
	{
		std::string cmd = "mkdir " + PATH_INITIAL;
		system(cmd.c_str());
	}

	if (_access(PATH_ARCHIVE.c_str(), 0) == -1)
	{
		std::string cmd = "mkdir " + PATH_ARCHIVE;
		system(cmd.c_str());
	}

	//if (_access(PATH_SELECTION.c_str(), 0) == -1)
	//{
	//	std::string cmd = "mkdir " + PATH_SELECTION;
	//	system(cmd.c_str());
	//}


	//if (_access(PATH_CROSSOVER.c_str(), 0) == -1)
	//{
	//	std::string cmd = "mkdir " + PATH_CROSSOVER;
	//	system(cmd.c_str());
	//}


	//if (_access(PATH_MUTATION.c_str(), 0) == -1)
	//{
	//	std::string cmd = "mkdir " + PATH_MUTATION;
	//	system(cmd.c_str());
	//}

	//if (_access(PATH_POPULATION.c_str(), 0) == -1)
	//{
	//	std::string cmd = "mkdir " + PATH_POPULATION;
	//	system(cmd.c_str());
	//}
}

Evolution::~Evolution()
{
	ClearVector(archive_);
	ClearVector(population_);
	ClearVector(initial_);

	while (energy_array_.size() != 0) {
		EnergyUpdate* temp_energy = energy_array_.back();
		energy_array_.pop_back();

		if (temp_energy) {
			delete temp_energy;
			temp_energy = NULL;
		}
	}
}

void Evolution::Run(double bound, int cone_size)
{
	N_CONE_THRESHOLD = cone_size;

	DISTORTION_BOUND = bound;
	N_CONE_STEP = DISTORTION_BOUND / N_CORE;

	InitialArchive();

	//system("pause");

	std::ofstream output_fitness(ROOT_PATH + "\\fitness_output.txt");

	// Evolution
	int cont = 0;
	int archive_cont = 0; // How many time the archive is fixed 
	while (cont < N_GENERATION && archive_cont < N_ARCHIVE_FIX)
	{
		std::cout << "Generation: " << cont << " th." << std::endl;

		output_fitness << cont
			<< " " << archive_[0]->fitness_energy_
			<< " " << archive_[0]->seam_energy_.rate_length_
			<< " " << archive_[0]->seam_energy_.seg_num_
			<< " " << archive_[0]->seam_energy_.seam_smooth_
			<< " " << archive_[0]->seam_energy_.dis_metric_.distortion_
			<< " " << archive_[0]->seam_energy_.small_cont_
			<< " " << archive_[0]->seam_energy_.narrow_cont_
			<< std::endl;
		
		//SegMesh seg_mesh(mesh_);
		//seg_mesh.SetSeam(archive_[0]->seam_status_);
		//seg_mesh.OutputSeg(ROOT_PATH + "\\iter_" + std::to_string(cont) + ".txt");

		StochasticSelection();

		RandomMutation(cont);

		Crossover(cont);

		bool update_status = ElitistReinsert(cont);

		//bool update_status = true;

		if (update_status)
		{
			archive_cont = 0;
			std::cout << "Archive is updated!" << std::endl;
		}
		else
		{
			archive_cont++;
			std::cout << "Archive is fixed " << archive_cont << " time(s)!" << std::endl;
		}

		cont++;
	}

	output_fitness << cont
		<< " " << archive_[0]->fitness_energy_
		<< " " << archive_[0]->seam_energy_.rate_length_
		<< " " << archive_[0]->seam_energy_.seg_num_
		<< " " << archive_[0]->seam_energy_.seam_smooth_
		<< " " << archive_[0]->seam_energy_.dis_metric_.distortion_
		<< " " << archive_[0]->seam_energy_.small_cont_
		<< " " << archive_[0]->seam_energy_.narrow_cont_
		<< std::endl;
	output_fitness.close();

	SegMesh seg_mesh(mesh_);
	seg_mesh.SetSeam(archive_[0]->seam_status_);
	seg_mesh.OutputSeg(ROOT_PATH + "\\output_seg.txt");
	seg_mesh.OutputMesh(ROOT_PATH + "\\output_mesh.obj");
}

void Evolution::FaceAngle()
{
	f_ang_.setConstant(mesh_.n_halfedges(), 0);
	for (const FH& f_h : mesh_.faces())
	{
		std::vector<int> he_set(3);
		std::vector<OpenMesh::Vec3d> he_v(3);
		int cont = 0;
		for (const HEH& he_h : mesh_.fh_range(f_h))
		{
			he_set[cont] = he_h.idx();
			he_v[cont] = mesh_.calc_edge_vector(he_h).normalized();
			cont++;
		}

		double cos_0 = -dot(he_v[1], he_v[2]);
		double cos_1 = -dot(he_v[2], he_v[0]);
		double cos_2 = -dot(he_v[0], he_v[1]);

		f_ang_[he_set[0]] = acos(cos_0);
		f_ang_[he_set[1]] = acos(cos_1);
		f_ang_[he_set[2]] = acos(cos_2);
	}
}

void Evolution::FaceArea()
{
	f_area_.setConstant(mesh_.n_faces(), 0);
	for (int i = 0; i < mesh_.n_faces(); i++)
	{
		f_area_[i] = mesh_.calc_face_area(FH(i));
	}
}

void Evolution::VertGauss()
{
	v_K_.setConstant(mesh_.n_vertices(), 2 * M_PI);

	for (const VH& v_h : mesh_.vertices())
	{
		if (mesh_.is_boundary(v_h)) {
			v_K_[v_h.idx()] = M_PI;
		}
	}

	for (const FH& f_h : mesh_.faces())
	{
		std::vector<int> he_set(3);
		std::vector<int> v_ids(3);
		int cont = 0;
		for (const HEH& he_h : mesh_.fh_range(f_h))
		{
			he_set[cont] = he_h.idx();
			v_ids[cont] = mesh_.opposite_vh(he_h).idx();
			cont++;
		}

		v_K_[v_ids[0]] -= f_ang_[he_set[0]];
		v_K_[v_ids[1]] -= f_ang_[he_set[1]];
		v_K_[v_ids[2]] -= f_ang_[he_set[2]];
	}
}

void Evolution::UpdateArchive()
{
	std::sort(archive_.begin(), archive_.end(), IndCmp);
	while (archive_.size() > N_ARCHIVE)
	{
		PopVector(archive_);
	}
}

void Evolution::RandomMutation(const int& cont)
{
	std::cout << "Random Mutation" << std::endl;

	std::vector<TopoChange> copy_set;
	for (int i = 0; i < pop_seg_.size(); i++)
	{
		copy_set.push_back(pop_seg_[i]);
	}

	change_type_.clear();
	change_type_.resize(pop_seg_.size(), C_NON);
	for (int pop_id = 0; pop_id < pop_seg_.size(); ++pop_id)
	{
		TopoChange& seg_mesh = pop_seg_[pop_id];
		const TopoChange& copy_seg = copy_set[pop_id];

		change_type& temp_type = change_type_[pop_id];

		if (cont < N_ITER_STEP)
		{
			int rand_case = rand() % 2;

			if (rand_case < 1)
			{
				// Remove boundary
				for (int i = 0; i < N_RANDOM_TRY; i++)
				{
					if (seg_mesh.RandomSeamRemove(e_l_, e_w_, v_K_))
					{
						temp_type = C_SEAM;
						break;
					}
				}
			}
			else if (rand_case < 2)
			{
				if (archive_[parent_[pop_id]]->seam_energy_.small_cont_ > 0)
				{
					for (int i = 0; i < N_RANDOM_TRY; i++)
					{
						if (seg_mesh.RandomSmallPatch(e_l_, e_w_))
						{
							temp_type = C_SMALL;
							break;
						}
					}
				}
				else
				{
					for (int i = 0; i < N_RANDOM_TRY; i++)
					{
						if (seg_mesh.RandomDijkstra(e_l_, e_w_))
						{
							temp_type = C_DIJ;
							break;
						}
					}
				}
			}

			if (seg_mesh.HighGenus())
			{
				copy_seg.OutputSeg(ROOT_PATH + "\\mutate_" + std::to_string(cont)
					+ "_ori_seg.txt");
				seg_mesh.OutputSeg(ROOT_PATH + "\\mutate_" + std::to_string(cont)
					+ "_high_genus_"
					+ std::to_string(temp_type) + ".txt");
			}
		}
		else if (cont < 2 * N_ITER_STEP)
		{
			int rand_case = rand() % 3;
			
			if (rand_case < 1)
			{
				// Remove boundary
				for (int i = 0; i < N_RANDOM_TRY; i++)
				{
					if (seg_mesh.RandomSeamRemove(e_l_, e_w_, v_K_, true, false))
					{
						temp_type = C_SEAM;
						break;
					}
				}
			}
			else if (rand_case < 2)
			{
				for (int i = 0; i < N_RANDOM_TRY; i++)
				{
					if (seg_mesh.RandomDijkstra(e_l_, e_w_))
					{
						temp_type = C_DIJ;
						break;
					}
				}
			}
			else
			{
				if (archive_[parent_[pop_id]]->seam_energy_.small_cont_ > 0)
				{
					for (int i = 0; i < N_RANDOM_TRY; i++)
					{
						if (seg_mesh.RandomSmallPatch(e_l_, e_w_, false))
						{
							temp_type = C_SMALL;
							break;
						}
					}
				}
				else
				{
					for (int i = 0; i < N_RANDOM_TRY; i++)
					{
						bool remove_status = false;
						int remove_iter = rand() % N_TRI_REMOVE + 1;
						//int remove_iter = 1;
						for (int j = 0; j < remove_iter; j++)
						{
							bool cur_status = seg_mesh.RandomOneTriangleChange(false);
							remove_status = remove_status
								|| cur_status;

							//std::cout << cur_status << std::endl;
						}

						if (remove_status)
						{
							temp_type = C_TRI;
							break;
						}
					}
				}
			}

			if (seg_mesh.HighGenus())
			{
				copy_seg.OutputSeg(ROOT_PATH + "\\mutate_" + std::to_string(cont)
					+ "_ori_seg.txt");
				seg_mesh.OutputSeg(ROOT_PATH + "\\mutate_" + std::to_string(cont)
					+ "_high_genus_"
					+ std::to_string(temp_type) + ".txt");
			}
		}
		else if (cont < 3 * N_ITER_STEP)
		{
			int rand_case = rand() % 3;

			if (rand_case < 1)
			{
				// Remove boundary
				for (int i = 0; i < N_RANDOM_TRY; i++)
				{
					//if (seg_mesh.RandomSeamRemove(e_l_, e_w_, v_K_, false, false))
					if (seg_mesh.RandomSeamRemove(e_l_, e_w_, v_K_, true, true))
					{
						temp_type = C_SEAM;
						break;
					}
				}
			}
			else if (rand_case < 2)
			{
				for (int i = 0; i < N_RANDOM_TRY; i++)
				{
					if (seg_mesh.RandomDijkstra(e_l_, e_w_))
					{
						temp_type = C_DIJ;
						break;
					}
				}
			}
			else
			{
				if (archive_[parent_[pop_id]]->seam_energy_.small_cont_ > 0)
				{
					for (int i = 0; i < N_RANDOM_TRY; i++)
					{
						if (seg_mesh.RandomSmallPatch(e_l_, e_w_, false))
						{
							temp_type = C_SMALL;
							break;
						}
					}
				}
				else if (archive_[parent_[pop_id]]->seam_energy_.narrow_cont_ > 0)
				{
					//copy_seg.OutputSeg(ROOT_PATH + "\\mutate_" + std::to_string(cont)
					//	+ "_narrow.txt");

					for (int i = 0; i < N_RANDOM_TRY; i++)
					{
						if (seg_mesh.RandomNarrowRemove(e_l_, e_w_,
							f_f_adj_, v_f_adj_, false))
						{
							temp_type = C_NARROW;
							break;
						}
					}

					//seg_mesh.OutputSeg(ROOT_PATH + "\\mutate_" + std::to_string(cont)
					//	+ "_after_narrow.txt");
				}
				else
				{
					for (int i = 0; i < N_RANDOM_TRY; i++)
					{
						bool remove_status = false;
						int remove_iter = rand() % N_TRI_REMOVE + 1;
						//int remove_iter = 1;
						for (int j = 0; j < remove_iter; j++)
						{
							bool cur_status = seg_mesh.RandomOneTriangleChange(false);
							remove_status = remove_status
								|| cur_status;

							//std::cout << cur_status << std::endl;
						}

						if (remove_status)
						{
							temp_type = C_TRI;
							break;
						}
					}
				}
			}

			if (seg_mesh.HighGenus())
			{
				copy_seg.OutputSeg(ROOT_PATH + "\\mutate_" + std::to_string(cont)
					+ "_ori_seg.txt");
				seg_mesh.OutputSeg(ROOT_PATH + "\\mutate_" + std::to_string(cont)
					+ "_high_genus_"
					+ std::to_string(temp_type) + ".txt");
			}
		}
		else
		{
			int rand_case = rand() % 3;

			if (rand_case < 1)
			{
				for (int i = 0; i < N_RANDOM_TRY; i++)
				{
					if (seg_mesh.RandomDijkstra(e_l_, e_w_))
					{
						temp_type = C_DIJ;
						break;
					}
				}
			}
			else if (rand_case < 2)
			{
				//if (archive_[parent_[pop_id]]->seam_energy_.small_cont_ > 0)
				//{
				//	for (int i = 0; i < N_RANDOM_TRY; i++)
				//	{
				//		if (seg_mesh.RandomSmallPatch(e_l_, e_w_))
				//		{
				//			temp_type = C_SMALL;
				//			break;
				//		}
				//	}
				//}
				//else if (archive_[parent_[pop_id]]->seam_energy_.narrow_cont_ > 0)
				//{
				//	//copy_seg.OutputSeg(ROOT_PATH + "\\mutate_" + std::to_string(cont)
				//	//	+ "_narrow.txt");

				//	for (int i = 0; i < N_RANDOM_TRY; i++)
				//	{
				//		if (seg_mesh.RandomNarrowRemove(e_l_, e_w_,
				//			f_f_adj_, v_f_adj_, false))
				//		{
				//			temp_type = C_NARROW;
				//			break;
				//		}
				//	}

				//	//seg_mesh.OutputSeg(ROOT_PATH + "\\mutate_" + std::to_string(cont)
				//	//	+ "_after_narrow.txt");
				//}
				//else
				//{
					for (int i = 0; i < N_RANDOM_TRY; i++)
					{
						bool remove_status = false;
						//int remove_iter = rand() % N_TRI_REMOVE + 1;
						int remove_iter = 1;
						for (int j = 0; j < remove_iter; j++)
						{
							bool cur_status = seg_mesh.RandomAdjRemove(f_ang_);
							remove_status = remove_status
								|| cur_status;

							//std::cout << cur_status << std::endl;
						}

						if (remove_status)
						{
							temp_type = C_ADJ;
							break;
						}
					}
				//}
			}
			else
			{
				for (int i = 0; i < N_RANDOM_TRY; i++)
				{
					bool remove_status = false;
					int remove_iter = rand() % N_TRI_REMOVE + 1;
					//int remove_iter = 1;
					for (int j = 0; j < remove_iter; j++)
					{
						bool cur_status = seg_mesh.RandomOneTriangleChange(false);
						remove_status = remove_status
							|| cur_status;

						//std::cout << cur_status << std::endl;
					}

					if (remove_status)
					{
						temp_type = C_TRI;
						break;
					}
				}
			}

			if (seg_mesh.HighGenus())
			{
				copy_seg.OutputSeg(ROOT_PATH + "\\mutate_" + std::to_string(cont)
					+ "_ori_seg.txt");
				seg_mesh.OutputSeg(ROOT_PATH + "\\mutate_" + std::to_string(cont)
					+ "_high_genus_"
					+ std::to_string(temp_type) + ".txt");
			}
		}
	}

	//OutputMutation();
}

void Evolution::StochasticSelection()
{
	std::cout << "Stochastic Selection" << std::endl;

	parent_.clear();
	pop_seg_.clear();

	int n_archive = archive_.size();

	// Pick
	std::vector<int> temp_vec;
	temp_vec.reserve(n_archive + 1);
	std::vector<int> picked_ids;
	for (int i = 0; i < n_archive; i++)
	{
		temp_vec.assign(n_archive - i, i);
		picked_ids.insert(picked_ids.end(), temp_vec.begin(), temp_vec.end());
	}

	for (size_t i = 0; i < N_RANDOM_PICK; i++)
	{
		int ind_id = picked_ids[rand() % picked_ids.size()];
		parent_.push_back(ind_id);
	}

	TopoChange seg_mesh(mesh_);
	for (const int& id : parent_)
	{
		seg_mesh.SetSeam(archive_[id]->seam_status_);
		seg_mesh.SetCone(archive_[id]->cone_status_);
		pop_seg_.push_back(seg_mesh);
	}

	for (int id = 0; id < N_ELITIST_PICK; id++)
	{
		seg_mesh.SetSeam(archive_[id]->seam_status_);
		seg_mesh.SetCone(archive_[id]->cone_status_);
		pop_seg_.push_back(seg_mesh);
		parent_.push_back(id);
	}

	//OutputSelection();
}

void Evolution::Crossover(const int& cont)
{
	std::cout << "Crossover" << std::endl;

	int n_population = pop_seg_.size();
	int n_archive = archive_.size();

	std::vector<std::pair<int, int>> cross_pair;

	// Crossover with the fitness order
	int id_0, id_1, better_id, worse_id;
	int cross_cont = 0;
	TopoChange new_seg(mesh_);
	TopoChange better_seg(mesh_);
	while (parent_.size() < N_CROSSOVER + n_population && cross_cont < N_MAX_ITER)
	{
		//std::cout << cont << std::endl;
		// 
		// From the mutated individuals
		
		id_0 = rand() % N_ELITIST_NUM;

		do {
			id_1 = rand() % N_ELITIST_NUM;
		} while (id_1 == id_0);

		//better_id = std::min(id_0, id_1);
		//worse_id = std::max(id_0, id_1);

		better_id = id_0;
		worse_id = id_1;

		cross_pair.emplace_back(worse_id, better_id);

		// Patch Inherit
		// Move good segmentation to worse individuals
		new_seg.SetSeam(archive_[worse_id]->seam_status_);
		new_seg.SetCone(archive_[worse_id]->cone_status_);
		better_seg.SetSeam(archive_[better_id]->seam_status_);
		better_seg.SetCone(archive_[better_id]->cone_status_);

		for (int i = 0; i < N_RANDOM_TRY; i++)
		{
			//std::cout << i << std::endl;

			if (new_seg.RandomPatchInherit(better_seg, e_l_, e_w_))
			{
				if (new_seg.HighGenus())
				{
					SegMesh worse_seg(mesh_, archive_[worse_id]->seam_status_);

					worse_seg.OutputSeg(ROOT_PATH + "\\cross_" + std::to_string(cont)
						+ "_worse_seg.txt");
					better_seg.OutputSeg(ROOT_PATH + "\\cross_" + std::to_string(cont)
						+ "_better_seg.txt");
					new_seg.OutputSeg(ROOT_PATH + "\\cross_" + std::to_string(cont)
						+ "_high_genus.txt");
				}

				pop_seg_.push_back(new_seg);
				parent_.push_back(worse_id);
				change_type_.emplace_back(C_CROSS);

				break;
			}
		}

		cross_cont++;
	}


	//// Selected ids
	//std::ofstream select_cout(PATH_CROSSOVER + "\\0_cross.txt");
	//for (size_t i = 0; i < cross_pair.size(); i++)
	//{
	//	select_cout << cross_pair[i].first << " " << cross_pair[i].second << std::endl;
	//}
	//select_cout.close();


	//// Segmentation
	//for (size_t i = n_population; i < pop_seg_.size(); i++)
	//{
	//	pop_seg_[i].OutputSeg(PATH_CROSSOVER + "\\seg_" + std::to_string(i - n_population) + ".txt");
	//}

	//OutputCrossover();
}

bool Evolution::ElitistReinsert(const int& cont)
{
	std::cout << "Elitist Reinsert" << std::endl;
	bool update_status = false;

	// Update fitness
	ClearVector(population_);
	population_.resize(pop_seg_.size());
#pragma omp parallel for
	for (int i = 0; i < pop_seg_.size(); i++)
	{
		EnergyUpdate& energy_update = *energy_array_[omp_get_thread_num()];
		Individual* ind = new Individual();

		ind->cone_status_ = pop_seg_[i].ReturnCone();
		ind->seam_status_ = pop_seg_[i].ReturnSeam();
		energy_update.SeamUpdate(ind->seam_status_,
			ind->seam_energy_, ind->fitness_energy_);

		population_[i] = ind;
	}

	//OutputPopulation();

	for (int i = 0; i < population_.size(); i++)
	{
		int parent_id = parent_[i];
		Individual* parent_ind = archive_[parent_id];
		Individual* child_ind = population_[i];

		double child_dis = child_ind->seam_energy_.dis_metric_.distortion_;
		double parent_dis = archive_[parent_[i]]->seam_energy_.dis_metric_.distortion_;
		//std::cout << i << " " << parent_[i] 
		// << " " << child_dis
		// << " " << parent_dis
		// << std::endl;

		TopoChange before_seg(mesh_, parent_ind->seam_status_, parent_ind->cone_status_);
		TopoChange after_seg(mesh_, child_ind->seam_status_, child_ind->cone_status_);

		if (parent_ind->fitness_energy_ > child_ind->fitness_energy_)
		{
			if (cont < N_ITER_STEP)
			{
				*parent_ind = *child_ind;
				update_status = true;
			}
			else if (cont < 2 * N_ITER_STEP
				&& (child_dis - parent_dis) < CHANGE_THRESHOLD_0 * parent_dis)
			{
				*parent_ind = *child_ind;
				update_status = true;
			}
			else if (cont < 3 * N_ITER_STEP 
				&& (child_dis - parent_dis) < CHANGE_THRESHOLD_1 * parent_dis)
			{
				*parent_ind = *child_ind;
				update_status = true;
			}
			else if ( (child_dis - parent_dis) < CHANGE_THRESHOLD_2 * parent_dis)
			{
				*parent_ind = *child_ind;
				update_status = true;
			}

			if (before_seg.HighGenus() || after_seg.HighGenus())
			{
				before_seg.OutputSeg(ROOT_PATH + "\\archive_before_" + std::to_string(change_type_[i]) + ".txt");
				after_seg.OutputSeg(ROOT_PATH + "\\archive_after" + std::to_string(cont)
					+ "_high_genus.txt");
			}
		}
	}

	UpdateArchive();

	if (cont % 100 == 99)
	{
		OutputArchive();
	}

	return update_status;
}

void Evolution::ClearVector(std::vector<Individual*>& vec)
{
	while (vec.size() != 0) {
		PopVector(vec);
	}
}

void Evolution::PopVector(std::vector<Individual*>& vec)
{
	Individual* ind = vec.back();
	vec.pop_back();

	if (ind){
		delete ind;
		ind = NULL;
	}
}

int FilesRead(std::string root, std::vector<std::string>& fileVec)
{
	int Nums = 0;
	long long handle = 0;
	struct _finddata_t fileinfo;
	std::string temp_str;
	if ((handle = _findfirst(temp_str.assign(root).append("\\*").c_str(), &fileinfo)) != -1)
	{
		do
		{
			if ((fileinfo.attrib & _A_SUBDIR))
			{
				if (strcmp(fileinfo.name, ".") != 0 && strcmp(fileinfo.name, "..") != 0)
					FilesRead(temp_str.assign(root).append(fileinfo.name).c_str(), fileVec);
			}
			else
			{
				try
				{
					if (fileinfo.size == 0)
						throw - 1;
					else
						fileVec.push_back(temp_str.assign(root).append("\\").append(fileinfo.name));
				}
				catch (int e)
				{
					if (e == -1)
						std::cout << "file is empty!" << std::endl;
				}
			}
		} while (_findnext(handle, &fileinfo) == 0);
		_findclose(handle);
	}

	std::cout << "Nums: " << fileVec.size() << std::endl;
	if (Nums > 0)
		return Nums;
	else
		return 0;
}

void Evolution::EdgeLength()
{
	e_l_.setConstant(mesh_.n_edges(), 0);
	for (const EH& e_h : mesh_.edges())
	{
		//e_l_[e_h.idx()] = 1;
		e_l_[e_h.idx()] = mesh_.calc_edge_length(e_h);
	}
}

void Evolution::EdgeWeight()
{
	e_w_.setConstant(mesh_.n_edges(), 0);
	OpenMesh::Vec3d n_0, n_1;
	for (const EH& e_h : mesh_.edges())
	{
		if (mesh_.is_boundary(e_h)) continue;

		n_0 = mesh_.normal(mesh_.face_handle(mesh_.halfedge_handle(e_h, 0)));
		n_1 = mesh_.normal(mesh_.face_handle(mesh_.halfedge_handle(e_h, 1)));

		e_w_[e_h.idx()] = pow(1 + (n_0 | n_1), 3);
	}
}

void Evolution::InitAdjRegion()
{
	OpenMesh::Vec3d max_box, min_box;
	max_box[0] = -DBL_MAX; max_box[1] = -DBL_MAX; max_box[2] = -DBL_MAX;
	min_box[0] = DBL_MAX; min_box[1] = DBL_MAX; min_box[2] = DBL_MAX;
	for (const VH& v_h : mesh_.vertices())
	{
		const OpenMesh::Vec3d& temp_pos = mesh_.point(v_h);
		max_box.maximize(temp_pos);
		min_box.minimize(temp_pos);
	}
	region_dis_ = P_SMALL_REGION * (max_box - min_box).norm();

	f_c_.resize(mesh_.n_faces(), 3);
	for (const FH& f_h : mesh_.faces())
	{
		OpenMesh::Vec3d temp_c = mesh_.calc_face_centroid(f_h);
		f_c_(f_h.idx(), 0) = temp_c[0];
		f_c_(f_h.idx(), 1) = temp_c[1];
		f_c_(f_h.idx(), 2) = temp_c[2];
	}

	f_f_adj_.clear();
	f_f_adj_.resize(mesh_.n_faces());
	v_f_adj_.clear();
	v_f_adj_.resize(mesh_.n_vertices());
	for (const FH& f_h : mesh_.faces())
	{
		std::vector<int>& temp_adj_f = f_f_adj_[f_h.idx()];

		FaceAdjRegion(f_h, f_c_, temp_adj_f);

		for (const int& adj_f_id : temp_adj_f)
		{
			for (const VH& adj_v : mesh_.fv_range(
				mesh_.face_handle(adj_f_id)))
			{
				v_f_adj_[adj_v.idx()].emplace_back(f_h.idx());
			}
		}
	}
}

void Evolution::FaceAdjRegion(const FH& f_h, const Eigen::MatrixX3d& f_c_, std::vector<int>& adj_f_ids)
{
	adj_f_ids.clear();
	for (const VH& adj_v : mesh_.fv_range(f_h))
	{
		for (const FH& adj_f : mesh_.vf_range(adj_v))
		{
			adj_f_ids.emplace_back(adj_f.idx());
		}
	}

	std::vector<double> dis(mesh_.n_faces(), INFINITY);
	dis[f_h.idx()] = 0;
	std::vector<int> Q;
	std::vector<int> U;
	std::vector<bool> Qf(mesh_.n_faces(), false);
	std::vector<bool> Uf(mesh_.n_faces(), false);
	Q.push_back(f_h.idx());
	Qf[f_h.idx()] = true;
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

		int faceIndex = Q[Index];
		Q[Index] = Q[Q.size() - 1];
		Q.pop_back();
		U.push_back(faceIndex);
		Uf[faceIndex] = true;

		for (const FH& adj_f : mesh_.ff_range(mesh_.face_handle(faceIndex)))
		{
			if (!Uf[adj_f.idx()])
			{
				//std::cout << to_v.idx() << std::endl;
				double distance = min +
					(f_c_.row(faceIndex) - f_c_.row(adj_f.idx())).norm();
				if (dis[adj_f.idx()] > distance && distance < region_dis_)
				{
					dis[adj_f.idx()] = distance;
					if (!Qf[adj_f.idx()])
					{
						Q.push_back(adj_f.idx());
						adj_f_ids.push_back(adj_f.idx());
						Qf[adj_f.idx()] = true;
					}
				}
			}
		}
	}

	std::sort(adj_f_ids.begin(), adj_f_ids.end());
	adj_f_ids.erase(std::unique(adj_f_ids.begin(), adj_f_ids.end()), adj_f_ids.end());
}

void Evolution::InitialPopulation()
{
	ComputeFixedCone();

	std::cout << "Initial Population" << std::endl;

	std::vector<VSA_cone> vsa_set(N_CONE_COEFF, mesh_);
	std::vector<std::vector<int>> cone_list_list(N_CONE_COEFF);
#pragma omp parallel for
	for (int i = 1; i <= N_CONE_COEFF; i++)
	{
		VSA_cone& vsa_cone = vsa_set[i - 1];
		cone_list_list[i - 1] = vsa_cone.GenerateCone(N_CONE_STEP * i);
	}

	// Energy of the fixed cones
	std::vector<bool> seam_v_status(mesh_.n_vertices(), false);
	for (int v_id : constraint_list_)
	{
		seam_v_status[v_id] = true;
	}

	for (std::vector<int>& cone_list : cone_list_list)
	{
		if (cone_list.size() == 0) continue;

		// Remove cones adjacent to the fixed cones
		for (size_t i = 0; i < cone_list.size(); i++)
		{
			for (const VH& vv : mesh_.vv_range(VH(cone_list[i])))
			{
				if (seam_v_status[vv.idx()])
				{
					cone_list[i] = vv.idx();
				}
			}
		}

		std::vector<bool> temp_status(mesh_.n_vertices(), false);
		for (int v_id : cone_list)
		{
			temp_status[v_id] = true;
		}

		// Merge the adjacent cones
		for (size_t i = 0; i < cone_list.size(); i++)
		{
			for (const VH& vv : mesh_.vv_range(VH(cone_list[i])))
			{
				if (temp_status[vv.idx()])
				{
					cone_list[i] = vv.idx();
				}
			}
		}

		cone_list.insert(cone_list.end(), constraint_list_.begin(), constraint_list_.end());

		std::sort(cone_list.begin(), cone_list.end());
		cone_list.erase(std::unique(cone_list.begin(), cone_list.end()), cone_list.end());
	}

	// Initial
	ClearVector(initial_);
	initial_.resize(N_CONE_COEFF * N_VSA_COEFF);
#pragma omp parallel for
	for (int i = 1; i <= N_CONE_COEFF; i++)
	{
		EnergyUpdate& energy_update = *energy_array_[omp_get_thread_num()];
		VSA_cone& vsa_cone = vsa_set[i - 1];

		vsa_cone.UpdateCone(cone_list_list[i - 1]);

		//vsa_cone.OutputCone(ROOT_PATH + "\\cone_" + std::to_string(N_CONE_STEP * i) + ".txt");

		for (int j = 0; j < N_VSA_COEFF; j++)
		{
			//std::cout << "VSA: " << (i - 1) * N_VSA_COEFF + j << std::endl;

			std::vector<bool> seam_status = vsa_cone.Generate(N_INIT_VSA + N_VSA_STEP * j);

			Individual* ini_ind = new Individual();
			ini_ind->cone_status_ = seam_v_status;
			ini_ind->seam_status_ = seam_status;
			energy_update.SeamUpdate(ini_ind->seam_status_,
				ini_ind->seam_energy_, ini_ind->fitness_energy_);

			initial_[(i - 1) * N_VSA_COEFF + j] = (ini_ind);
		}
	}

	for (size_t i = 0; i < initial_.size(); i++)
	{
		TopoChange topo_mesh(mesh_);
		topo_mesh.SetSeam(initial_[i]->seam_status_);
		if (topo_mesh.HighGenus())
		{
			topo_mesh.OutputSeg(ROOT_PATH + "\\init_" + std::to_string(i)
				+ "_high_genus.txt");
		}
	}

	std::sort(initial_.begin(), initial_.end(), IndCmp);
	OutputInitial();
}

void Evolution::ComputeFixedCone()
{
	std::cout << "Compute Fixed Cones: " << DISTORTION_BOUND
		<< " " << N_CONE_THRESHOLD
		<< std::endl;

	// Determine the step of GenerateCone()
	VSA_cone vsa_cone(mesh_);
	std::vector<int> cone_list;

	constraint_list_ = vsa_cone.GenerateCone(DISTORTION_BOUND);
	if (constraint_list_.size() < N_CONE_THRESHOLD)
	{
		int cont = 0;
		do
		{
			cone_list = vsa_cone.GenerateCone(N_CONE_STEP * N_CONE_COEFF);
			constraint_list_ = vsa_cone.GenerateCone(DISTORTION_BOUND);

			std::cout << N_CONE_STEP * N_CONE_COEFF << std::endl;
			std::cout << cone_list.size() << std::endl;
			std::cout << constraint_list_.size() << std::endl;

			if (cone_list.size() > 0)
			{
				break;
			}
			else
			{
				DISTORTION_BOUND /= DECREASE_STEP;
				N_CONE_STEP = DISTORTION_BOUND / N_CORE;
			}
			cont++;
		} while (cont < N_MAX_ITER);
	}
	else
	{
		//do
		//{
		//	cone_list = vsa_cone.GenerateCone(N_CONE_STEP * N_CONE_COEFF);
		//	constraint_list_ = vsa_cone.GenerateCone(DISTORTION_BOUND);

		//	std::cout << N_CONE_STEP * N_CONE_COEFF << std::endl;
		//	std::cout << cone_list.size() << std::endl;
		//	std::cout << constraint_list_.size() << std::endl;

		//	if (cone_list.size() == 0 || constraint_list_.size() < N_CONE_THRESHOLD)
		//	{
		//		break;
		//	}
		//	else
		//	{
		//		DISTORTION_BOUND *= DECREASE_STEP;
		//		N_CONE_STEP = DISTORTION_BOUND / N_CORE;
		//	}
		//} while (true);

		std::ofstream output_large(ROOT_PATH + "\\large_cone_size.txt");
		output_large.close();
		
		exit(EXIT_FAILURE);
	}

	// Merge the adjacent cones
	for (size_t i = 0; i < constraint_list_.size(); i++)
	{
		for (size_t j = i + 1; j < constraint_list_.size(); j++)
		{
			for (const VH& vv : mesh_.vv_range(
				mesh_.vertex_handle(constraint_list_[i])))
			{
				if (vv.idx() == constraint_list_[j])
				{
					constraint_list_[j] = constraint_list_[i];
				}
			}
		}
	}

	std::sort(constraint_list_.begin(), constraint_list_.end());
	constraint_list_.erase(std::unique(constraint_list_.begin(),
		constraint_list_.end()), constraint_list_.end());

	for (size_t i = 0; i < energy_array_.size(); i++) {
		energy_array_[i]->SetBound(DISTORTION_BOUND);
	}

	std::ofstream max_cone_file(ROOT_PATH + "\\max_cone_" 
		+ std::to_string(N_CONE_STEP * N_CONE_COEFF) + ".txt");
	for (size_t i = 0; i < cone_list.size(); i++)
	{
		max_cone_file << cone_list[i] << std::endl;
	}
	max_cone_file.close();

	std::ofstream cone_file(ROOT_PATH + "\\cone_" 
		+ std::to_string(DISTORTION_BOUND) + ".txt");
	for (size_t i = 0; i < constraint_list_.size(); i++)
	{
		cone_file << constraint_list_[i] << std::endl;
	}
	cone_file.close();
}

void Evolution::InitialArchive()
{
	InitialPopulation();

	std::cout << "Initial Archive" << std::endl;

	ClearVector(archive_);
	for (int i = 0; i < initial_.size(); i++)
	{
		if (initial_[i]->seam_energy_.dis_metric_.distortion_ > DISTORTION_BOUND) break;

		if (archive_.size() != 0 && 
			*archive_.back() == *initial_[i]) continue;

		Individual* new_ind = new Individual();

		*new_ind = *initial_[i];
		archive_.emplace_back(new_ind);
	}

	int cont = 0;
	int max_iter = 10000;
	while (archive_.size() < N_ARCHIVE)
	{
		for (int i = 0; i < initial_.size(); i++)
		{
			if (initial_[i]->seam_energy_.dis_metric_.distortion_ > DISTORTION_BOUND) break;

			Individual* new_ind = new Individual();

			*new_ind = *initial_[i];
			archive_.emplace_back(new_ind);
		}

		cont++;

		if (cont > max_iter)
		{
			exit(EXIT_FAILURE);
		}
	}

	UpdateArchive();

	//N_ITER_STEP = std::max(N_ITER_STEP, 2 * archive_.front()->seam_energy_.seg_num_);
	N_GENERATION = 5 * N_ITER_STEP;
	//N_GENERATION = 4 * N_ITER_STEP;

	OutputArchive();
}

void Evolution::OutputInitial()
{
	// Segmentation
	SegMesh seg_mesh(mesh_);
	for (size_t i = 0; i < initial_.size(); i++)
	{
		seg_mesh.SetSeam(initial_[i]->seam_status_);
		seg_mesh.OutputSeg(PATH_INITIAL + "\\seg_" + std::to_string(i) + ".txt");
	}

	// Fitness
	std::ofstream fit_cout(PATH_INITIAL + "\\0_initial_fitness.txt");
	for (size_t i = 0; i < initial_.size(); i++)
	{
		fit_cout << i
			<< " " << initial_[i]->fitness_energy_
			<< " " << initial_[i]->seam_energy_.rate_length_
			<< " " << initial_[i]->seam_energy_.seg_num_
			<< " " << initial_[i]->seam_energy_.seam_smooth_
			<< " " << initial_[i]->seam_energy_.dis_metric_.distortion_
			<< std::endl;
	}
	fit_cout.close();
}

void Evolution::OutputArchive()
{
	// Segmentation
	SegMesh seg_mesh(mesh_);
	for (size_t i = 0; i < archive_.size(); i++)
	{
		seg_mesh.SetSeam(archive_[i]->seam_status_);
		seg_mesh.OutputSeg(PATH_ARCHIVE + "\\seg_" + std::to_string(i) + ".txt");

		std::ofstream cone_cout(PATH_ARCHIVE + "\\z_cone_" + std::to_string(i) + ".txt");
		const std::vector<bool>& cone_status = archive_[i]->cone_status_;
		for (int i = 0; i < cone_status.size(); i++)
		{
			if (cone_status[i])
			{
				cone_cout << i << std::endl;
			}
		}
		cone_cout.close();
	}

	// Fitness
	std::ofstream fit_cout(PATH_ARCHIVE + "\\0_archive_fitness.txt");
	for (size_t i = 0; i < archive_.size(); i++)
	{
		fit_cout << i
			<< " " << archive_[i]->fitness_energy_
			<< " " << archive_[i]->seam_energy_.rate_length_
			<< " " << archive_[i]->seam_energy_.seg_num_
			<< " " << archive_[i]->seam_energy_.seam_smooth_
			<< " " << archive_[i]->seam_energy_.dis_metric_.distortion_
			<< " " << archive_[i]->seam_energy_.small_cont_
			<< " " << archive_[i]->seam_energy_.narrow_cont_
			<< std::endl;
	}
	fit_cout.close();
}

void Evolution::OutputSelection()
{
	// Selected ids
	std::ofstream select_cout(PATH_SELECTION + "\\0_selection.txt");
	for (size_t i = 0; i < parent_.size(); i++)
	{
		select_cout << parent_[i] <<std::endl;
	}
	select_cout.close();

	// Segmentation
	for (size_t i = 0; i < pop_seg_.size(); i++)
	{
		pop_seg_[i].OutputSeg(
			PATH_SELECTION + "\\seg_" + std::to_string(i) + ".txt");

		std::ofstream cone_cout(PATH_SELECTION + "\\z_cone_" + std::to_string(i) + ".txt");
		const std::vector<bool>& cone_status = pop_seg_[i].ReturnCone();
		for (int i = 0; i < cone_status.size(); i++)
		{
			if (cone_status[i])
			{
				cone_cout << i << std::endl;
			}
		}
		cone_cout.close();
	}
}

void Evolution::OutputMutation()
{
	// Selected ids
	std::ofstream select_cout(PATH_MUTATION + "\\0_mutation.txt");
	for (size_t i = 0; i < parent_.size(); i++)
	{
		select_cout << parent_[i] << std::endl;
	}
	select_cout.close();

	// Segmentation
	for (size_t i = 0; i < pop_seg_.size(); i++)
	{
		pop_seg_[i].OutputSeg(PATH_MUTATION + "\\seg_" + std::to_string(i) + ".txt");

		std::ofstream cone_cout(PATH_MUTATION + "\\z_cone_" + std::to_string(i) + ".txt");
		const std::vector<bool>& cone_status = pop_seg_[i].ReturnCone();
		for (int i = 0; i < cone_status.size(); i++)
		{
			if (cone_status[i])
			{
				cone_cout << i << std::endl;
			}
		}
		cone_cout.close();
	}

	std::ofstream change_cout(PATH_MUTATION + "\\0_change_type.txt");
	for (size_t i = 0; i < change_type_.size(); i++)
	{
		change_cout << i << " " << change_type_[i] << std::endl;
	}
	change_cout.close();
}

void Evolution::OutputCrossover()
{
	// Segmentation
	for (size_t i = 0; i < pop_seg_.size(); i++)
	{
		pop_seg_[i].OutputSeg(PATH_CROSSOVER + "\\seg_" + std::to_string(i) + ".txt");
	}
}

void Evolution::OutputPopulation()
{
	// Segmentation
	SegMesh seg_mesh(mesh_);
	for (size_t i = 0; i < population_.size(); i++)
	{
		if (!population_[i]) continue;

		seg_mesh.SetSeam(population_[i]->seam_status_);
		seg_mesh.OutputSeg(PATH_POPULATION + "\\seg_" + std::to_string(i) + ".txt");
	}

	// Fitness
	std::ofstream fit_cout(PATH_POPULATION 
		+ "\\0_population_fitness.txt");
	for (size_t i = 0; i < population_.size(); i++)
	{
		if (!population_[i]) continue;

		fit_cout << i
			<< " " << population_[i]->fitness_energy_
			<< " " << population_[i]->seam_energy_.rate_length_
			<< " " << population_[i]->seam_energy_.seg_num_
			<< " " << population_[i]->seam_energy_.seam_smooth_
			<< " " << population_[i]->seam_energy_.dis_metric_.distortion_
			<< std::endl;
	}
	fit_cout.close();

	// Selected ids
	std::ofstream change_cout(PATH_POPULATION + "\\0_change_type.txt");
	for (size_t i = 0; i < change_type_.size(); i++)
	{
		change_cout << i << " " << change_type_[i] << std::endl;
	}
	change_cout.close();
}