#pragma once
#ifndef LocalWrapping_h
#define LocalWrapping_h

#include"Common.h"
#include <algorithm>
using namespace std;

enum class Border
{
	BORDER_TOP = 0,//��
	BORDER_BOTTOM = 1,//��
	BORDER_LEFT = 2,//��
	BORDER_RIGHT = 3//��
};

enum class SeamDirection
{
	SEAM_VERTICAL = 0,//��ֱ
	SEAM_HORIZONTAL = 1//ˮƽ
};

class LocalWarpping
{
private:
	bool Is_transparent(CVMat& mask, int row, int col);//�ж�mask��(row,col)���Ƿ�Ϊ��
	void init_displacement(vector<vector<Coordinate>>& displacement, int rows, int cols);//��ʼ��displacement
	
public:
	//Sobel��Ե��⣬���ڻ�ȡ��������
	CVMat Sobel_img(CVMat src);
	//��ͼ���Ե�ҵ���Ŀհ�,������ʼ����������,����borderѡ����ͼ������������ѡ��seam
	pair<int, int> Choose_longest_border(CVMat src, CVMat mask, Border& direction);
	//��Borderȷ������ͼ��ʹ��DP��ȡ������С��seam
	int* Get_local_seam(CVMat src, CVMat mask, SeamDirection seamdirection, pair<int, int> begin_end);
	//���ݻ�ȡ��seam���������ص�ƽ��
	CVMat Insert_local_seam(CVMat src, CVMat& mask, int* seam, SeamDirection seamdirection, pair<int, int> begin_end, bool shiftToend);
	//��ʾ��ѡ����border
	CVMat Show_longest_border(CVMat src, pair<int, int>begin_end, Border direction);
	//������������ӳ�䵽ԭͼ����
	void warp_mesh_back(vector<vector<CoordinateDouble>>& mesh, vector<vector<Coordinate>> displacementMap, Config config);
	//�ھ��λ�֮���ͼ���Ϸ��ñ�׼��������(20*20)
	vector<vector<CoordinateDouble>> get_rectangle_mesh(CVMat src, Config config);
	//��ȡloaclwarp��ת������
	vector<vector<Coordinate>> Get_Local_warp_displacement(CVMat src, CVMat mask);
	//Local Warp������
	vector<vector<Coordinate>> Local_warp(CVMat src, CVMat& wrap_img, CVMat mask);
};


#endif