#include "Compute\Evolution.h"
#include <chrono>
#include <io.h>

const int INPUT_NUM = 4;
const std::string INPUT_OBJ = "[INPUT_OBJ]";
const std::string OUTPUT_PATH = "[OUTPUT_PATH]";
const std::string DIS_THRESHOLD = "[DISTORTION_THRESHOLD]";
const std::string CONE_BOUND = "[CONE_BOUND]";

void printHelp(const std::string& programName)
{
	std::cout << "Help: "
		<< programName
		<< " " << INPUT_OBJ << " "
		<< " " << OUTPUT_PATH << " "
		<< " " << DIS_THRESHOLD << " "
		<< " " << CONE_BOUND << " "
		<< std::endl;
}

bool doesArgExist(const std::string& arg, const std::string& searchStr)
{
	return arg.find(searchStr) != std::string::npos;
}

bool parseArg(const std::string& arg, const std::string& searchStr, std::string& value)
{
	if (doesArgExist(arg, searchStr)) {
		value = arg.substr(arg.find_first_of(searchStr[searchStr.size() - 1]) + 1);
		return true;
	}
	return false;
}

void parseArgs(int argc, const char* argv[], std::string& obj, std::string& output_path)
{
	if (argc < INPUT_NUM) {
		// Input and/or output path not specified
		printHelp(argv[0]);
		exit(EXIT_FAILURE);
	}
	else {
		// parse arguments
		obj = argv[1];
		output_path = argv[2];

		//std::string value;
		//for (int i = 3; i < argc; i++) {
		//	if (parseArg(argv[i], "--dist=", value)) thres = std::stod(value);
		//	if (parseArg(argv[i], "--maxIters=", value)) iters = std::stod(value);
		//}
	}
}

int FilePathsRead(std::string root, std::vector<std::string>& fileVec)
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
				{
					fileVec.push_back(temp_str.assign("\\").append(fileinfo.name));
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

void AreaScaling(Mesh& mesh)
{
	double Area = MeshTools::Area(mesh);
	double sqrtArea = sqrt(Area)  * 10;

	for (auto vh : mesh.vertices())
	{
		mesh.set_point(vh,
			mesh.point(vh) / sqrtArea);
	}
}

int main(int argc, const char* argv[])
{
	if (argc < INPUT_NUM)
	{
		printHelp(argv[0]);
		exit(EXIT_FAILURE);
	}

	auto start_time = std::chrono::system_clock::now();

	std::string obj = "";
	std::string out_path = "";

	double dis_threshold = 0.025;
	int cone_size = 300;

	obj = argv[1];
	out_path = argv[2];

	dis_threshold = atof(argv[3]);
	cone_size = atoi(argv[4]);

	std::cout << "=======================================" << std::endl;
	std::cout << "Model: " << obj << std::endl;

	Mesh input_mesh;
	MeshTools::ReadMesh(input_mesh, obj);
	input_mesh.update_normals();

	// Normalized
	AreaScaling(input_mesh);

	// Set seed
	srand(0);

	// Output path
	if (_access(out_path.c_str(), 0) == -1)
	{
		std::string cmd = "mkdir " + out_path;
		system(cmd.c_str());
	}

	Evolution evo_dev(input_mesh);
	evo_dev.DebugPath(out_path);
	evo_dev.Run(dis_threshold, cone_size);

	auto end_time = std::chrono::system_clock::now();
	auto duration = std::chrono::duration_cast
		<std::chrono::microseconds>(end_time - start_time);

	std::ofstream output_time(out_path + "\\time.txt");
	output_time << double(duration.count())
		* std::chrono::microseconds::period::num
		/ std::chrono::microseconds::period::den << std::endl;
	output_time.close();

	return 0;
}