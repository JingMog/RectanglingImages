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
	Vector2d pstart = Vector2d::Zero();//��ʼ
	Vector2d pend = Vector2d::Zero();//����
	double angle = 0;//��ת�Ƕ�
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
	Config config;//Config�����,����һЩ�������s
	CVMat src;//ԭͼ��

protected:
	//˫����Ȩ��ת��ΪMatrix
	MatrixXd BilinearWeightsToMatrix(BilinearWeights w);
	//���line1��line2�Ƿ��ཻ,���ؽ�������ֵ
	Vector2d detectIntersect(Matrix2d line1, Matrix2d line2, bool& isintersection);
	//�޸�mask��Ե
	void revise_mask_for_lines(CVMat& mask);
	//�жϵ�point�Ƿ���topLeft,topRight,bottomLeft,bottomRight��ɵ��ı����ڲ�
	bool is_in_quad(CoordinateDouble point, CoordinateDouble topLeft, CoordinateDouble topRight, CoordinateDouble bottomLeft, CoordinateDouble bottomRight);
	//�ж�ֱ��Line�Ƿ���mask��
	bool line_in_mask(CVMat mask, LineD line);
	//����lsd����ֱ�߼��,���ؼ�⵽��ֱ��
	vector<LineD> lsd_detect(CVMat mask);
	//���������k��b���ж������Ե�Ƿ��ֱ���ཻ
	bool does_segment_intersect_line(LineD lineSegment, double slope, double intersect, bool vertical, CoordinateDouble& intersectPoint);
	//���ı���������߶ν��зָ���ؽ��㼯
	vector<CoordinateDouble> intersections_with_quad(LineD lineSegment, CoordinateDouble topLeft, CoordinateDouble topRight, CoordinateDouble bottomLeft, CoordinateDouble bottomRight);
	//����ά���߶ηָ���չ����һά��vector<LineD>
	vector<LineD> flatten(vector<vector<vector<LineD>>> lineSeg);

	BilinearWeights get_bilinear_weights(CoordinateDouble point, Coordinate upperLeftIndices, vector<vector<CoordinateDouble>> mesh);

	vector<vector<vector<pair<int, double>>>> cal_theta(vector<vector<vector<Line_rotate>>> lineSeg, Config config);

	SpareseMatrixD_Row block_diag(SpareseMatrixD_Row origin, MatrixXd addin, int QuadID, Config config);

	

public:
	//���캯��
	GlobalWarpping(Config& conf, CVMat& src);
	//��ȡ�����е�Shape Preservations��������
	SpareseMatrixD_Row get_shape_mat(vector<vector<CoordinateDouble>> mesh);
	//��ȡmesh��λ��(row, col)�����������񶥵���������
	VectorXd get_vertices(int row, int col, vector<vector<CoordinateDouble>>& mesh);
	//��ȡ�����е�Boundary Constraints��������
	pair<SpareseMatrixD_Row, VectorXd> get_boundary_mat(vector<vector<CoordinateDouble>> mesh);
	//�����������ж�����ֱ�߽��зָ�
	vector<vector<vector<LineD>>> segment_line_in_quad(vector<LineD> lines, vector<vector<CoordinateDouble>> mesh);


	SpareseMatrixD_Row get_vertex_to_shape_mat(vector<vector<CoordinateDouble>> mesh);

	
	//��ȡ�����е�Line Preservations��������
	SpareseMatrixD_Row get_line_mat(CVMat mask, vector<vector<CoordinateDouble>> mesh, vector<double>rotate_theta, vector<vector<vector<LineD>>> lineSeg, vector<pair<MatrixXd, MatrixXd>>& BilinearVec, int& linenum, vector<bool>& bad);
	
	//
	vector<vector<vector<LineD>>> init_line_seg(CVMat mask, vector < LineD >& lineSeg_flatten, vector<vector<CoordinateDouble>> mesh, vector<pair<int, double>>& id_theta, vector<double>& rotate_theta);
	

};

#endif