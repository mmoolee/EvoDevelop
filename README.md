EvoDevelop
===

This is a C++ implementation of the piecewise developable approximation in the following paper:

[Evolutionary Piecewise Developable Approximations](https://doi.org/10.1145/3592140).
Zheng-Yu Zhao\*, Mo Li\*, Zheng Zhang, [Qing Fang](https://qingfang1208.github.io/)\*, [Ligang Liu](http://staff.ustc.edu.cn/~lgliu/), [Xiao-Ming Fu](https://ustc-gcl-f.github.io/).
*ACM Transactions on Graphics (SIGGRAPH)*, 42(4), 2023. (*Joint first authors)

The code is written by using Microsoft Visual Studio 2019, with OpenMP support for parallel programming.

## Dependencies
* [OpenMesh](https://www.graphics.rwth-aachen.de/software/openmesh/)
* [Eigen](http://eigen.tuxfamily.org/)
* OpenMP support

## Usage
* A generated .exe file is included in ```x64/Release```, and run with:
```
EvoDevelop.exe [INPUT_OBJ] [OUTPUT_PATH] [DISTORTION_THRESHOLD] [CONE_BOUND]
```
Command arguments:
* [INPUT_OBJ]: a input mesh (.obj).
* [OUTPUT_PATH]: the output path.
* [DISTORTION_THRESHOLD]: the distortion threshold (a range [0.02, 0.03] is recommended, default: 0.025).
* [CONE_BOUND]: the cone bound to control the number of cone singularities (default: 1000).

## Output
* fitness_output.txt: text file including all terms in the fitness energy in each iteration.
* output_seg.txt: the partition represented by different face labels.