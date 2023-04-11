#pragma once
#pragma once

#ifndef GlobalWarpping_h
#define GlobalWarpping_h

#include"Common.h"

#define clamp(x,a,b)    (  ((a)<(b))  \
? ((x)<(a))?(a):(((x)>(b))?(b):(x))	\
: ((x)<(b))?(b):(((x)>(a))?(a):(x))	\
)


struct Line_rotate 
{
	Vector2d pstart = Vector2d::Zero();//开始
	Vector2d pend = Vector2d::Zero();//结束
	double angle = 0;//旋转角度
	Line_rotate(Vector2d pstart, Vector2d pend, double angle)
	{
		this->pstart = pstart;
		this->pend = pend;
		this->angle = angle;
	}
};

struct BilinearWeights 
{
	double s;
	double t;
};

class GlobalWarpping
{
private:
	Config config;//Config类对象,保存一些网格参数s
	CVMat src;//原图像

protected:
	//双线性权重转化为Matrix
	MatrixXd BilinearWeightsToMatrix(BilinearWeights w);
	//检测line1与line2是否相交,返回交点坐标值
	Vector2d detectIntersect(Matrix2d line1, Matrix2d line2, bool& isintersection);
	//修改mask边缘
	void revise_mask_for_lines(CVMat& mask);
	//判断点point是否在topLeft,topRight,bottomLeft,bottomRight组成地四边形内部
	bool is_in_quad(CoordinateDouble point, CoordinateDouble topLeft, CoordinateDouble topRight, CoordinateDouble bottomLeft, CoordinateDouble bottomRight);
	//判断直线Line是否在mask中
	bool line_in_mask(CVMat mask, LineD line);
	//调用lsd进行直线检测,返回检测到的直线
	vector<LineD> lsd_detect(CVMat mask);
	//根据网格的k和b来判断网格边缘是否和直线相交
	bool does_segment_intersect_line(LineD lineSegment, double slope, double intersect, bool vertical, CoordinateDouble& intersectPoint);
	//用四边形网格对线段进行分割，返回交点集
	vector<CoordinateDouble> intersections_with_quad(LineD lineSegment, CoordinateDouble topLeft, CoordinateDouble topRight, CoordinateDouble bottomLeft, CoordinateDouble bottomRight);
	//将三维的线段分割结果展开成一维的vector<LineD>
	vector<LineD> flatten(vector<vector<vector<LineD>>> lineSeg);

	BilinearWeights get_bilinear_weights(CoordinateDouble point, Coordinate upperLeftIndices, vector<vector<CoordinateDouble>> mesh);

	vector<vector<vector<pair<int, double>>>> cal_theta(vector<vector<vector<Line_rotate>>> lineSeg, Config config);

	SpareseMatrixD_Row block_diag(SpareseMatrixD_Row origin, MatrixXd addin, int QuadID, Config config);

	

public:
	//构造函数
	GlobalWarpping(Config& conf, CVMat& src);
	//获取论文中的Shape Preservations能量矩阵
	SpareseMatrixD_Row get_shape_mat(vector<vector<CoordinateDouble>> mesh);
	//获取mesh中位于(row, col)的四邻域网格顶点坐标坐标
	VectorXd get_vertices(int row, int col, vector<vector<CoordinateDouble>>& mesh);
	//获取论文中的Boundary Constraints能量矩阵
	pair<SpareseMatrixD_Row, VectorXd> get_boundary_mat(vector<vector<CoordinateDouble>> mesh);
	//在整个网格中对所有直线进行分割
	vector<vector<vector<LineD>>> segment_line_in_quad(vector<LineD> lines, vector<vector<CoordinateDouble>> mesh);


	SpareseMatrixD_Row get_vertex_to_shape_mat(vector<vector<CoordinateDouble>> mesh);

	
	//获取论文中的Line Preservations能量矩阵
	SpareseMatrixD_Row get_line_mat(CVMat mask, vector<vector<CoordinateDouble>> mesh, vector<double>rotate_theta, vector<vector<vector<LineD>>> lineSeg, vector<pair<MatrixXd, MatrixXd>>& BilinearVec, int& linenum, vector<bool>& bad);
	
	//
	vector<vector<vector<LineD>>> init_line_seg(CVMat mask, vector < LineD >& lineSeg_flatten, vector<vector<CoordinateDouble>> mesh, vector<pair<int, double>>& id_theta, vector<double>& rotate_theta);
	

};

#endif