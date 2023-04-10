#pragma once
#ifndef LocalWrapping_h
#define LocalWrapping_h

#include"Common.h"
#include <algorithm>
using namespace std;

enum class Border
{
	BORDER_TOP = 0,//上
	BORDER_BOTTOM = 1,//下
	BORDER_LEFT = 2,//左
	BORDER_RIGHT = 3//右
};

enum class SeamDirection
{
	SEAM_VERTICAL = 0,//垂直
	SEAM_HORIZONTAL = 1//水平
};

class LocalWarpping
{
private:
	bool Is_transparent(CVMat& mask, int row, int col);//判断mask的(row,col)处是否为空
	void init_displacement(vector<vector<Coordinate>>& displacement, int rows, int cols);//初始化displacement
	
public:
	//Sobel边缘检测，便于获取能量函数
	CVMat Sobel_img(CVMat src);
	//在图像边缘找到最长的空白,返回起始与最终坐标,根据border选择子图计算能量函数选择seam
	pair<int, int> Choose_longest_border(CVMat src, CVMat mask, Border& direction);
	//从Border确定的子图中使用DP求取能量最小的seam
	int* Get_local_seam(CVMat src, CVMat mask, SeamDirection seamdirection, pair<int, int> begin_end);
	//根据获取的seam来进行像素的平移
	CVMat Insert_local_seam(CVMat src, CVMat& mask, int* seam, SeamDirection seamdirection, pair<int, int> begin_end, bool shiftToend);
	//显示所选择的最长border
	CVMat Show_longest_border(CVMat src, pair<int, int>begin_end, Border direction);
	//将矩形网格反向映射到原图像中
	void warp_mesh_back(vector<vector<CoordinateDouble>>& mesh, vector<vector<Coordinate>> displacementMap, Config config);
	//在矩形化之后的图像上放置标准矩形网格(20*20)
	vector<vector<CoordinateDouble>> get_rectangle_mesh(CVMat src, Config config);
	//获取loaclwarp的转化矩阵
	vector<vector<Coordinate>> Get_Local_warp_displacement(CVMat src, CVMat mask);
	//Local Warp主函数
	vector<vector<Coordinate>> Local_warp(CVMat src, CVMat& wrap_img, CVMat mask);
};


#endif