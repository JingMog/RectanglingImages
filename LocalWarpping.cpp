#include "LocalWarpping.h"

bool cmp(const pair<int, float> a, const pair<int, float> b)
{
	return a.second < b.second;
}

/*
	判断mask的(row,col)处是否为空(即黑色)
*/
bool LocalWarpping::Is_transparent(CVMat& mask, int row, int col)
{
	if (mask.at<uchar>(row, col) == 0) 
		return false;
	else 
		return true;
}

/*
	初始化displacement
*/
void LocalWarpping::init_displacement(vector<vector<Coordinate>>& displacement, int rows, int cols)
{
	for (int row = 0; row < rows; row++) 
	{
		vector<Coordinate> displacement_row;
		for (int col = 0; col < cols; col++)
		{
			Coordinate c;
			displacement_row.push_back(c);
		}
		displacement.push_back(displacement_row);
	}
}

/*
	使用Sobel对图像进行边缘检测
*/
CVMat LocalWarpping::Sobel_img(CVMat src)
{
	CVMat gray;
	//将RGB图转化为灰度图进行边缘检测
	cv::cvtColor(src, gray, CV_BGR2GRAY);
	CVMat grad_x, grad_y, dst;
	//Sobel(输入的灰度图像,输出的梯度,输出梯度的数据类型,
	//dx=1,dy=0,计算x方向上的梯度
	cv::Sobel(gray, grad_x, CV_8U, 1, 0);
	//dx=0,dy=1,计算y方向上的梯度
	cv::Sobel(gray, grad_y, CV_8U, 0, 1);
	addWeighted(grad_x, 0.5, grad_y, 0.5, 0, dst);//将gray_x与gray_y等权组合
	return dst;
}

/*
	在矩形化之后的图像上放置一个标准网格,
	然后需要进行wrap back返回到原始不规则图像中,得到网格
*/
void LocalWarpping::warp_mesh_back(vector<vector<CoordinateDouble>>& mesh, vector<vector<Coordinate>> displacementMap, Config config)
{
	int meshnum_row = config.meshNumRow;//网格的行数
	int meshnum_col = config.meshNumCol;//网格的列数

	for (int row_mesh = 0; row_mesh < meshnum_row; row_mesh++)
	{
		for (int col_mesh = 0; col_mesh < meshnum_col; col_mesh++)
		{
			//如果位于边缘
			if (row_mesh == meshnum_row - 1 && col_mesh == meshnum_col - 1)
			{
				CoordinateDouble& meshVertexPoint = mesh[row_mesh][col_mesh];//网格坐标点
				Coordinate vertexDisplacement = displacementMap[floor(meshVertexPoint.row) - 1][floor(meshVertexPoint.col) - 1];
				meshVertexPoint.row += vertexDisplacement.row;
				meshVertexPoint.col += vertexDisplacement.col;
			}
			//非边缘,直接将原坐标加上displacement
			CoordinateDouble& meshVertexPoint = mesh[row_mesh][col_mesh];
			Coordinate vertexDisplacement = displacementMap[(int)floor(meshVertexPoint.row)][(int)floor(meshVertexPoint.col)];
			meshVertexPoint.row += vertexDisplacement.row;
			meshVertexPoint.col += vertexDisplacement.col;
		}
	}
}

/*
	在图像边缘找到最长的空白,返回起始与最终坐标
	之后根据该坐标选择子图,计算子图的能量函数,即Sobel边缘检测
	找到能量最短的一条线,然后执行Seam Carving
*/
pair<int, int> LocalWarpping::Choose_longest_border(CVMat src, CVMat mask, Border& direction)
{
	//返回最长border的[begin，end]
	int maxLength = 0;//记录最大长度
	int rows = src.rows;//rows
	int cols = src.cols;//cols

	int final_startIndex = 0;//结果起始坐标
	int final_endIndex = 0;//结果最终坐标

	//上下左右四个方向找最长的border,然后对其进行seam carving
	//left
	int tmp_maxLength, tmp_startIndex, tmp_endIndex;
	tmp_maxLength = tmp_startIndex = tmp_endIndex = 0;
	bool isCounting = false;
	for (int row = 0; row < rows; row++)
	{
		//如果mask中(row, 0)点不为空(空白中断),或者为最后一行,则计数中止
		if (!Is_transparent(mask, row, 0) || row == rows - 1)
		{
			if (isCounting)
			{	//如果还在计数，到此为止
				if (Is_transparent(mask, row, 0))
				{
					tmp_endIndex++;
					tmp_maxLength++;
					isCounting = true;
				}
				if (tmp_maxLength > maxLength)
				{
					//更新最大值
					maxLength = tmp_maxLength;
					final_startIndex = tmp_startIndex;
					final_endIndex = tmp_endIndex;
					direction = Border::BORDER_LEFT;
				}
			}
			isCounting = false;//中止计数
			tmp_startIndex = tmp_endIndex = row;
			tmp_maxLength = 0;
		}
		else
		{
			//该点透明,继续计数
			tmp_endIndex++;
			tmp_maxLength++;
			isCounting = true;
		}
	}

	//右侧
	tmp_maxLength = tmp_startIndex = tmp_endIndex = 0;
	isCounting = false;
	for (int row = 0; row < rows; row++) 
	{
		//如果mask中(row, cols - 1)点不为空(空白中断),或者为最后一行,则计数中止
		if (!Is_transparent(mask, row, cols - 1) || row == rows - 1) 
		{
			if (isCounting) 
			{
				//如果还在计数，到此为止
				if (Is_transparent(mask, row, cols - 1)) 
				{
					tmp_endIndex++;
					tmp_maxLength++;
					isCounting = true;
				}
				if (tmp_maxLength > maxLength) 
				{
					//更新最大值
					maxLength = tmp_maxLength;
					final_startIndex = tmp_startIndex;
					final_endIndex = tmp_endIndex;
					direction = Border::BORDER_RIGHT;
				}
			}
			isCounting = false;//中止计数
			tmp_startIndex = tmp_endIndex = row;
			tmp_maxLength = 0;
		}
		else 
		{
			//该点透明，开始计数
			tmp_endIndex++;
			tmp_maxLength++;
			isCounting = true;
		}
	}

	//上
	tmp_maxLength = tmp_startIndex = tmp_endIndex = 0;
	isCounting = false;
	for (int col = 0; col < cols; col++) 
	{
		//如果mask中(0, col)点不为空(空白中断),或者为最后一行,则计数中止
		if (!Is_transparent(mask, 0, col) || col == cols - 1) 
		{
			if (isCounting) 
			{	//如果还在计数，到此为止
				if (Is_transparent(mask, 0, col))
				{
					tmp_endIndex++;
					tmp_maxLength++;
					isCounting = true;
				}
				if (tmp_maxLength > maxLength)
				{
					//更新最大值
					maxLength = tmp_maxLength;
					final_startIndex = tmp_startIndex;
					final_endIndex = tmp_endIndex;
					direction = Border::BORDER_TOP;
				}
			}
			isCounting = false;//中止计数
			tmp_startIndex = tmp_endIndex = col;
			tmp_maxLength = 0;
		}
		else 
		{	//该点透明，开始计数
			tmp_endIndex++;
			tmp_maxLength++;
			isCounting = true;
		}
	}

	//下方
	tmp_maxLength = tmp_startIndex = tmp_endIndex = 0;
	isCounting = false;
	for (int col = 0; col < cols; col++) 
	{
		//如果mask中(row - 1, col)点不为空(空白中断),或者为最后一行,则计数中止
		if (!Is_transparent(mask, rows - 1, col) || col == cols - 1)
		{
			if (isCounting)
				{//如果还在计数，到此为止
				if (Is_transparent(mask, rows - 1, col)) 
				{
					tmp_endIndex++;
					tmp_maxLength++;
					isCounting = true;
				}
				if (tmp_maxLength > maxLength) 
				{
					//更新最大值
					maxLength = tmp_maxLength;
					final_startIndex = tmp_startIndex;
					final_endIndex = tmp_endIndex;
					direction = Border::BORDER_BOTTOM;
				}
			}
			isCounting = false;//中止计数
			tmp_startIndex = tmp_endIndex = col;
			tmp_maxLength = 0;
		}
		else 
		{	//该点透明，开始计数
			tmp_endIndex++;
			tmp_maxLength++;
			isCounting = true;
		}
	}
	//cout << "Border maxLength: " << maxLength << endl;

	//cout << Is_transparent(src.at<Vec3b>(0,final_endIndex));
	//system("pause");
	return make_pair(final_startIndex, final_endIndex - 1);
}

/*
	显示所选择的最长border
	根据choose_longest_border选择最长的border，然后将该Border设置为蓝色
*/
CVMat LocalWarpping::Show_longest_border(CVMat src, pair<int, int>begin_end, Border direction)
{
	CVMat tmpsrc;
	src.copyTo(tmpsrc);//tmpsrc = src
	int rows = src.rows;
	int cols = src.cols;
	switch (direction) 
	{
	case Border::BORDER_LEFT:
		for (int row = begin_end.first; row < begin_end.second; row++) 
			tmpsrc.at<cv::Vec3b>(row, 0) = cv::Vec3b(0, 0, 255);
		break;
	case Border::BORDER_RIGHT:
		for (int row = begin_end.first; row < begin_end.second; row++)
			tmpsrc.at<cv::Vec3b>(row, cols - 1) = cv::Vec3b(0, 0, 255);
		break;
	case Border::BORDER_TOP:
		for (int col = begin_end.first; col < begin_end.second; col++)
			tmpsrc.at<cv::Vec3b>(0, col) = cv::Vec3b(0, 0, 255);
		break;
	case Border::BORDER_BOTTOM:
		for (int col = begin_end.first; col < begin_end.second; col++)
			tmpsrc.at<cv::Vec3b>(rows - 1, col) = cv::Vec3b(0, 0, 255);
		break;
	default:
		break;
	}

	cv::namedWindow("Border", CV_WINDOW_AUTOSIZE);
	cv::imshow("Border", tmpsrc);
	cv::waitKey(0);
	return tmpsrc;
}


/*
	根据获取的seam来进行像素的平移
	插入一条seam(Seam Carving算法)
*/
CVMat LocalWarpping::Insert_local_seam(CVMat src, CVMat& mask, int* seam, SeamDirection seamdirection, pair<int, int> begin_end, bool shiftToend)
{
	//如果seam垂直,则将src和mask转置,转化为水平
	if (seamdirection == SeamDirection::SEAM_HORIZONTAL) 
	{
		transpose(src, src);
		transpose(mask, mask);
	}

	CVMat resimg;
	src.copyTo(resimg);

	int begin = begin_end.first;//最长border所在的local row范围
	int end = begin_end.second;
	int rows = src.rows;
	int cols = src.cols;

	for (int row = begin; row <= end; row++)
	{
		//对于每一行，seam右侧（左侧）的所有像素均往边缘平移
		int local_row = row - begin;//在border子图中的位置
		if (!shiftToend)
		{
			//seam左侧所有像素点，往左移动一格
			for (int col = 0; col < seam[local_row]; col++)
			{
				colorPixel pixel = src.at<colorPixel>(row, col + 1);
				resimg.at<colorPixel>(row, col) = pixel;
				mask.at<uchar>(row, col) = mask.at<uchar>(row, col + 1);
			}
		}
		else
		{
			//seam右侧所有像素点，往右移动一格
			for (int col = cols - 1; col > seam[local_row]; col--)
			{
				colorPixel pixel = src.at<colorPixel>(row, col - 1);
				resimg.at<colorPixel>(row, col) = pixel;
				mask.at<uchar>(row, col) = mask.at<uchar>(row, col - 1);
			}
		}

		//对缝隙两边进行平过渡（左右像素平均值）
		mask.at<uchar>(row, seam[local_row]) = 0;
		if (seam[local_row] == 0) 
		{
			//如果缝隙位于第一列,则将右侧像素填充至第一列
			resimg.at<colorPixel>(row, seam[local_row]) = src.at<colorPixel>(row, seam[local_row] + 1);
		}
		else if (seam[local_row] == cols - 1)
		{
			//如果缝隙位于最后一列，则将左侧像素填充至最后一列
			resimg.at<colorPixel>(row, seam[local_row]) = src.at<colorPixel>(row, seam[local_row] - 1);
		}
		else 
		{
			//否则将seam像素左右进行平滑过度
			colorPixel pixel1 = src.at<colorPixel>(row, seam[local_row] + 1);
			colorPixel pixel2 = src.at<colorPixel>(row, seam[local_row] - 1);
			resimg.at<colorPixel>(row, seam[local_row]) = pixel1 / 2 + pixel2 / 2;
		}
	}

	//重新转置
	if (seamdirection == SeamDirection::SEAM_HORIZONTAL)
	{
		cv::transpose(resimg, resimg);
		cv::transpose(mask, mask);
	}

	return resimg;
}

/*
	从Border确定的子图中使用DP求取能量最小的seam
	使用动态规划算法，首席按按获取图像的能量函数，然后进行路径回溯，最终路径保存在int* seam中
*/
int* LocalWarpping::Get_local_seam(CVMat src, CVMat mask, SeamDirection seamdirection, pair<int, int> begin_end)
{
	//如果为水平,则进行转置将其转化为垂直
	if (seamdirection == SeamDirection::SEAM_HORIZONTAL)
	{
		cv::transpose(src, src);
		cv::transpose(mask, mask);
	}

	//统一寻找竖直的seam
	int rows = src.rows;
	int cols = src.cols;

	//行，列的范围
	int row_start = begin_end.first;
	int row_end = begin_end.second;
	int col_start = 0;
	int col_end = cols - 1;

	int range = row_end - row_start + 1;//空白Border的长度为end - start + 1

	int outputWidth = cols;
	int outputHeight = range;

	CVMat displayimg;
	src.copyTo(displayimg);

	//获取子图,对子图进行DP求解seam
	CVMat local_img = displayimg(cv::Range::Range(row_start, row_end + 1), cv::Range::Range(col_start, col_end + 1));//获取子图
	CVMat local_mask = mask(cv::Range::Range(row_start, row_end + 1), cv::Range::Range(col_start, col_end + 1));//获取mask子图，Range左闭右开，所有右边界需要+1

	CVMat local_energy = Sobel_img(local_img);//Sobel边缘检测获取图像的能量
	CVMat local_energy_32f;
	local_energy.convertTo(local_energy_32f, CV_32F);

	//对图像边缘(非矩形部分)处设置能量为INF,防止seam穿过这一部分
	//论文中的3.1Mesh-free Local Warping
	for (int row = 0; row < range; row++)
		for (int col = col_start; col <= col_end; col++)
			if ((int)local_mask.at<uchar>(row, col) == 255)
				local_energy_32f.at<float>(row, col) = INF;


	//动态规划求解最小能量的seam，（最短路径的优化）
	CVMat tmpenergy;
	local_energy_32f.copyTo(tmpenergy);
	//1.遍历所有子图中的像素点，求解能量图
	for (int row = 1; row < range; row++)
	{
		for (int col = col_start; col <= col_end; col++)
		{
			//如果是第一列,考虑上方,和右上方
			if (col == col_start)
				tmpenergy.at<float>(row, col) += min(tmpenergy.at<float>(row - 1, col), tmpenergy.at<float>(row - 1, col + 1));
			//如果是最后一列,考虑上方,左上方
			else if (col == col_end)
				tmpenergy.at<float>(row, col) += min(tmpenergy.at<float>(row - 1, col - 1), tmpenergy.at<float>(row - 1, col));
			//正常point,考虑左上方,上方和右上方
			else
				tmpenergy.at<float>(row, col) += min(tmpenergy.at<float>(row - 1, col), min(tmpenergy.at<float>(row - 1, col - 1), tmpenergy.at<float>(row - 1, col + 1)));
		}
	}

	vector<pair<int, float>> last_row;//存储最后一行的列数和对应的能量值
	for (int col = col_start; col <= col_end; col++)
	{
		last_row.push_back(make_pair(col, tmpenergy.at<float>(range - 1, col)));
	}

	sort(last_row.begin(), last_row.end(), cmp);//对能量进行排序
	int* seam = new int[range];
	seam[range - 1] = last_row[0].first;

	//2.路径回溯，根据得到的能量图求取最优seam
	for (int row = range - 2; row >= 0; row--)
	{
		if (seam[row + 1] == col_start)
		{
			//如果seam下一行为子图第一列
			//回溯上，右上
			if (tmpenergy.at<float>(row, seam[row + 1] + 1) < tmpenergy.at<float>(row, seam[row + 1]))
				seam[row] = seam[row + 1] + 1;
			else
				seam[row] = seam[row + 1];
		}
		else if (seam[row + 1] == col_end)
		{
			//seam[row + 1]为子图的最后一列
			//回溯左上，上
			if (tmpenergy.at<float>(row, seam[row + 1] - 1) < tmpenergy.at<float>(row, seam[row + 1]))
			{
				seam[row] = seam[row + 1] - 1;
			}
			else
			{
				seam[row] = seam[row + 1];
			}
		}
		else
		{
			//(row, seam[row+1]-1), (row, seam[row+1]), (row, seam[row+1]+1)三者中的最小值
			//回溯左上，上，右上三个
			float min_energy = min(tmpenergy.at<float>(row, seam[row + 1] - 1), min(tmpenergy.at<float>(row, seam[row + 1]), tmpenergy.at<float>(row, seam[row + 1] + 1)));
			if (min_energy == tmpenergy.at<float>(row, seam[row + 1] - 1))
			{
				seam[row] = seam[row + 1] - 1;
			}
			else if (min_energy == tmpenergy.at<float>(row, seam[row + 1] + 1))
			{
				seam[row] = seam[row + 1] + 1;
			}
			else
			{
				seam[row] = seam[row + 1];
			}
			
		}
	}

	/*
	//显示找到的这条Seam
	for (int i = 0; i < range; i++)
		cout << seam[i] << " ";
	for (int row = 0; row < range; row++)
		local_img.at<colorPixel>(row, seam[row]) = colorPixel(255, 255, 0);
	cv::namedWindow("local_seam", CV_WINDOW_AUTOSIZE);
	cv::imshow("local_seam", local_img);
	cv::waitKey(0);
	system("pause");
	*/
	return seam;
}

/*
	在矩形化之后的图像上放置标准矩形网格(20*20)
*/
vector<vector<CoordinateDouble>> LocalWarpping::get_rectangle_mesh(CVMat src, Config config) 
{
	int rows = config.rows;
	int cols = config.cols;
	int meshnum_row = config.meshNumRow;//网格的行数,论文中为20
	int meshnum_col = config.meshNumCol;//网格的列数,论文中为20
	double row_per_mesh = config.rowPermesh;//每一个网格所占的行数
	double col_per_mesh = config.colPermesh;//每一个网格所占的列数

	vector<vector<CoordinateDouble>> mesh;
	for (int row_mesh = 0; row_mesh < meshnum_row; row_mesh++) 
	{
		vector<CoordinateDouble> meshrow;
		for (int col_mesh = 0; col_mesh < meshnum_col; col_mesh++) 
		{
			CoordinateDouble coord;
			coord.row = row_mesh * row_per_mesh;//网格交点的横坐标
			coord.col = col_mesh * col_per_mesh;//网格交点的纵坐标
			meshrow.push_back(coord);
		}
		mesh.push_back(meshrow);
	}
	return mesh;
}


/*
	获取图像的转化矩阵,即每个像素的移动向量
*/
vector<vector<Coordinate>> LocalWarpping::Get_Local_warp_displacement(CVMat src, CVMat mask)
{
	int rows = src.rows;
	int cols = src.cols;

	vector<vector<Coordinate>> displacementMap;//转换矩阵
	vector<vector<Coordinate>> finaldisplacementMap;
	init_displacement(finaldisplacementMap, rows, cols);//初始化displacement
	init_displacement(displacementMap, rows, cols);//初始化displacement

	while (true)
	{
		Border direction;
		pair<int, int> begin_end = Choose_longest_border(src, mask, direction);//选择mask中最长的空白边,方向保存在direction中
		//Show_longest_border(src, begin_end, direction);
		
		if (begin_end.first == begin_end.second)
		{
			//所有像素全部移动完成，直接返回displacement
			/*cv::namedWindow("Border", CV_WINDOW_AUTOSIZE);
			cv::imshow("Border", src);
			cv::waitKey(0);*/
			return displacementMap;
		}
		else
		{
			//仍有空白边界未平移完成
			bool shift_to_end = false;
			SeamDirection seamdirection;//缝隙的方向（水平 or 竖直）
			switch (direction)
			{
			case Border::BORDER_LEFT:
				seamdirection = SeamDirection::SEAM_VERTICAL;//border位于左,即seam为垂直
				shift_to_end = false;
				break;
			case Border::BORDER_RIGHT:
				seamdirection = SeamDirection::SEAM_VERTICAL;//border位于右,即seam为垂直
				shift_to_end = true;
				break;
			case Border::BORDER_TOP:
				seamdirection = SeamDirection::SEAM_HORIZONTAL;//border位于上,即seam为水平
				shift_to_end = false;
				break;
			case Border::BORDER_BOTTOM:
				seamdirection = SeamDirection::SEAM_HORIZONTAL;//border位于下,即seam为水平
				shift_to_end = true;
				break;
			default:
				break;
			}

			int* seam = Get_local_seam(src, mask, seamdirection, begin_end);//使用DP获取能量最小seam

			src = Insert_local_seam(src, mask, seam, seamdirection, begin_end, shift_to_end);//插入一条Seam，即平移像素，消除空白

			//更新置换矩阵
			for (int row = 0; row < rows; row++) 
			{
				for (int col = 0; col < cols; col++) 
				{
					Coordinate tmpdisplacement;
					if (seamdirection == SeamDirection::SEAM_VERTICAL && row >= begin_end.first && row <= begin_end.second)
					{
						//如果为垂直缝,并且row在缝隙内
						int local_row = row - begin_end.first;//子图中的row
						if (col > seam[local_row] && shift_to_end) 
						{
							tmpdisplacement.col = -1;//向右移动
						}
						else if (col < seam[local_row] && !shift_to_end) 
						{
							tmpdisplacement.col = 1;//向左移动
						}
					}
					else if (seamdirection == SeamDirection::SEAM_HORIZONTAL && col >= begin_end.first && col <= begin_end.second) 
					{
						//如果为水平缝,并且row在缝隙内
						int local_col = col - begin_end.first;//子图中的col
						if (row > seam[local_col] && shift_to_end) 
						{
							tmpdisplacement.row = -1;//向下移动
						}
						else if (row < seam[local_col] && !shift_to_end) 
						{
							tmpdisplacement.row = 1;//向上移动
						}
					}

					//更新置换矩阵中的[row][col]
					Coordinate& finaldisplacement = finaldisplacementMap[row][col];
					int tmpdisplace_row = row + tmpdisplacement.row;
					int tmpdisplace_col = col + tmpdisplacement.col;
					Coordinate displacementOftarget = displacementMap[tmpdisplace_row][tmpdisplace_col];
					int rowInOrigin = tmpdisplace_row + displacementOftarget.row;//
					int colInOrigin = tmpdisplace_col + displacementOftarget.col;
					finaldisplacement.row = rowInOrigin - row;//tmpdisplacement.row + displacementOftarget.row
					finaldisplacement.col = colInOrigin - col;//tmpdisplacement.col + displacementOftarget.col
				}
			}
			//displacement = finalDisplacement
			//便于下次继续更新矩阵
			for (int row = 0; row < rows; row++) 
			{
				for (int col = 0; col < cols; col++) 
				{
					Coordinate& displacement = displacementMap[row][col];
					Coordinate finalDisplacement = finaldisplacementMap[row][col];
					displacement.row = finalDisplacement.row;
					displacement.col = finalDisplacement.col;
				}
			}

		}
	}
	return displacementMap;
}


/*
	LocalWarp主函数,返回转化矩阵,warp_img为拉伸之后的图像
*/
vector<vector<Coordinate>> LocalWarpping::Local_warp(CVMat src, CVMat& wrap_img, CVMat mask) 
{
	vector<vector<Coordinate>> displacementMap = Get_Local_warp_displacement(src, mask);//获取置换矩阵
	for (int row = 0; row < src.rows; row++)
	{
		for (int col = 0; col < src.cols; col++)
		{
			Coordinate displacement = displacementMap[row][col];
			colorPixel pixel = src.at<colorPixel>(row + displacement.row, col + displacement.col);
			wrap_img.at<colorPixel>(row, col) = pixel;//根据转化矩阵生成拉伸之后的图像
		}
	}
	return displacementMap;
}