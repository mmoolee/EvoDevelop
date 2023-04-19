#include "Opt/types.h"
#include "../../MeshDefinition.h"

class ConesFlattening
{
public:
	ConesFlattening(const Mesh& mesh);

	void initCoef(double sigma);

	void geneCone(VectorX& conesK);

private:
	void areaMat();

	void initYamabeCoef();

	void interiorMat();

private:
	const Mesh& mesh;
	ColMajorSparseMatrix L, A, P;
	VectorX K;

	int lp = 2;
	double sigma;
};

//namespace ConesFlattening
//{
//	void initCoef(const Mesh& mesh, double sigma_);
//	void geneCone(VectorX& conesK);
//}



