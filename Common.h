#pragma once
#ifndef Common_h
#define Common_h

//论文中提到的无穷大
#define INF 1e8 
#define PI 3.14159265358979323846

#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>

#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <Eigen/Sparse> //稀疏矩阵的代数运算
#include <Eigen/Dense> //稠密矩阵的代数运算

#include "lsd.h"

using namespace std;

typedef cv::Mat CVMat;
typedef cv::Vec3b colorPixel;//RGB

typedef Eigen::SparseMatrix<double> SparseMatrixD;
typedef Eigen::SparseMatrix<double, Eigen::RowMajor> SpareseMatrixD_Row;//行主序矩阵
typedef Eigen::VectorXd VectorXd;
typedef Eigen::MatrixXd MatrixXd;
typedef Eigen::Vector2d Vector2d;
typedef Eigen::Vector2i Vector2i;
typedef Eigen::MatrixXi MatrixXi;
typedef Eigen::Matrix2d Matrix2d;
typedef Eigen::Triplet<double> T;
typedef Eigen::SimplicialCholesky<SparseMatrixD> CSolve;


class Config
{
public:
	//运行前配置
	int rows;
	int cols;
	int meshNumRow;
	int meshNumCol;
	int meshQuadRow;
	int meshQuadCol;
	double rowPermesh;
	double colPermesh;
	Config() {};
	Config(int rows, int cols, int meshNumRow, int meshNumCol)
	{
		this->rows = rows;
		this->cols = cols;
		this->meshNumRow = meshNumRow;
		this->meshNumCol = meshNumCol;
		this->meshQuadCol = meshNumCol - 1;
		this->meshQuadRow = meshNumRow - 1;
		this->rowPermesh = (double(rows) - 1) / (double(meshNumRow) - 1);
		this->colPermesh = double(double(cols) - 1) / (double(meshNumCol) - 1);
	}
};

//坐标类
class Coordinate
{
public:
	int row;
	int col;

	//重载运算符 ==
	bool operator==(const Coordinate& rhs) const
	{
		return (row == rhs.row && col == rhs.col);
	}

	//重载运算符 <
	bool operator<(const Coordinate& rhs) const
	{
		// this operator is used to determine equality, so it must use both x and y
		if (row < rhs.row)
			return true;
		if (row > rhs.row)
			return false;
		return col < rhs.col;
	}

	//构造函数
	Coordinate() { row = 0; col = 0; };
	Coordinate(int setRow, int setCol) { row = setRow; col = setCol; };
};

//double坐标类
//Coordinate，row和col为double
class CoordinateDouble
{
public:
	double row;
	double col;

	//重载运算符 == 
	bool operator==(const CoordinateDouble& rhs) const
	{
		return (row == rhs.row && col == rhs.col);
	}

	//重载运算符 < 
	bool operator<(const CoordinateDouble& rhs) const
	{
		// this operator is used to determine equality, so it must use both x and y
		if (row < rhs.row)
			return true;
		if (row > rhs.row)
			return false;
		return col < rhs.col;
	}

	//重载输出流
	friend ostream& operator<<(ostream& stream, const CoordinateDouble& p)
	{
		stream << "(" << p.col << "," << p.row << ")";
		return stream;
	}

	CoordinateDouble() { row = 0; col = 0; };
	CoordinateDouble(double setRow, double setCol) { row = setRow; col = setCol; };
};

//直线
class LineD
{
public:
	double row1, col1;
	double row2, col2;

	//构造函数
	LineD(double row1, double col1, double row2, double col2)
	{
		this->row1 = row1;
		this->row2 = row2;
		this->col1 = col1;
		this->col2 = col2;
	}
	//构造函数
	LineD()
	{
		row1 = 0; col1 = 0; row2 = 0; col2 = 0;
	}
	LineD(CoordinateDouble p1, CoordinateDouble p2)
	{
		row1 = p1.row; row2 = p2.row; col1 = p1.col; col2 = p2.col;
	}
};


CVMat Mask_contour(CVMat src);

VectorXd mesh_to_vector(vector<vector<CoordinateDouble>> mesh, Config config);

vector<vector<CoordinateDouble>> vector_to_mesh(VectorXd x, Config config);

void print_sparse_mat(SparseMatrixD Q);

SpareseMatrixD_Row row_stack(SparseMatrixD origin, SpareseMatrixD_Row diag);

SpareseMatrixD_Row row_stack(SpareseMatrixD_Row origin, SpareseMatrixD_Row diag);

MatrixXd row_stack(MatrixXd mat1, MatrixXd mat2);

MatrixXd col_stack(MatrixXd mat1, MatrixXd mat2);

void DrawLine(CVMat& img, CoordinateDouble coordstart, CoordinateDouble coordend);

void DrawLine(CVMat& img, LineD line);

CVMat drawmesh(CVMat src, vector<vector<CoordinateDouble>> mesh, Config config);

void enlarge_mesh(vector<vector<CoordinateDouble>>& mesh, double enlarge_factor, Config config);


#endif