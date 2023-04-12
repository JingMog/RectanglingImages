#include"GlobalWarpping.h"
#pragma warning(disable: 4244)

/*
	默认构造函数,初始化config
*/
GlobalWarpping::GlobalWarpping(Config& conf, CVMat& src)
{
	this->config = conf;//初始化config
	this->src = src;//初始化src
}

/*
	将BilinearWeights转化为MatrixXd
*/
MatrixXd GlobalWarpping::BilinearWeightsToMatrix(BilinearWeights w)
{
	MatrixXd mat(2, 8);
	double v1w = 1 - w.s - w.t + w.s * w.t;
	double v2w = w.s - w.s * w.t;
	double v3w = w.t - w.s * w.t;
	double v4w = w.s * w.t;
	mat << v1w, 0, v2w, 0, v3w, 0, v4w, 0,
		0, v1w, 0, v2w, 0, v3w, 0, v4w;
	return mat;
}

/*
	获取双线性权值
*/
BilinearWeights GlobalWarpping::get_bilinear_weights(CoordinateDouble point, Coordinate upperLeftIndices, vector<vector<CoordinateDouble>> mesh)
{
	//获取四边形网格的四个坐标
	CoordinateDouble p1 = mesh[upperLeftIndices.row][upperLeftIndices.col];//左上角
	CoordinateDouble p2 = mesh[upperLeftIndices.row][upperLeftIndices.col + 1];//右上角
	CoordinateDouble p3 = mesh[upperLeftIndices.row + 1][upperLeftIndices.col];//左下角
	CoordinateDouble p4 = mesh[upperLeftIndices.row + 1][upperLeftIndices.col + 1];//右下角

	double slopeTop    = (p2.row - p1.row) / (p2.col - p1.col);//上边界斜率
	double slopeBottom = (p4.row - p3.row) / (p4.col - p3.col);//下边界斜率
	double slopeLeft   = (p1.row - p3.row) / (p1.col - p3.col);//左边界斜率
	double slopeRight  = (p2.row - p4.row) / (p2.col - p4.col);//右边界斜率

	double quadraticEpsilon = 0.1;//

	//上下边界平行并且左右边界平行
	if (slopeTop == slopeBottom && slopeLeft == slopeRight)
	{
		// method 3
		Matrix2d mat1;
		//上边界的长,左边界的长,上边界的宽,左边界的宽
		mat1 << p2.col - p1.col, p3.col - p1.col,
			p2.row - p1.row, p3.row - p1.row;

		MatrixXd mat2(2, 1);
		mat2 << point.col - p1.col, point.row - p1.row;//point到p1的长,point到p1的宽

		MatrixXd matsolution = mat1.inverse() * mat2;//mat1^-1 * mat2

		BilinearWeights weights;
		weights.s = matsolution(0, 0);
		weights.t = matsolution(1, 0);
		return weights;
	}
	//只有左右边界平行
	else if (slopeLeft == slopeRight)
	{
		//method 2
		//a = (x2-x1)*(y4-y3) - (y2-y1)*(x4-x3)
		double a = (p2.col - p1.col) * (p4.row - p3.row) - (p2.row - p1.row) * (p4.col - p3.col);
		//b = y*((x4-x3)-(x2-x1)) - x*((y4-y3)-(y2-y1)) + x1*(y4-y3) - y1*(x4-x3) + y3*(x2-x1)-x3*(y2-y1)
		double b = point.row * ((p4.col - p3.col) - (p2.col - p1.col)) - point.col * ((p4.row - p3.row) - (p2.row - p1.row)) + p1.col * (p4.row - p3.row) - p1.row * (p4.col - p3.col) + (p2.col - p1.col) * (p3.row) - (p2.row - p1.row) * (p3.col);
		//c = y*(x3-x1) - x*(y3-y1) + x1*y3 - x3*y1;
		double c = point.row * (p3.col - p1.col) - point.col * (p3.row - p1.row) + p1.col * p3.row - p3.col * p1.row;

		//求方程ax^2 + bx + c = 0的根
		double s1 = (-1 * b + sqrt(b * b - 4 * a * c)) / (2 * a);//第一个根
		double s2 = (-1 * b - sqrt(b * b - 4 * a * c)) / (2 * a);//第二个根
		//cout << "s1: " << s1 << " s2: " << s2 << endl;
		double s;
		if (0 <= s1 && s1 <= 1) 
		{
			s = s1;
		}
		else if (0 <= s2 && s2 <= 1) 
		{
			s = s2;
		}
		else 
		{
			//
			if ((s1 > 1 && s1 - quadraticEpsilon < 1) || (s2 > 1 && s2 - quadraticEpsilon < 1)) 
			{
				s = 1;
			}
			else if ((s1 < 0 && s1 + quadraticEpsilon > 0) || (s2 < 0 && s2 + quadraticEpsilon > 0)) 
			{
				s = 0;
			}
			else 
			{
				// this case should not happen
				//cout << "Could not interpolate s weight for coordinate (" << point.col << "," << point.row << ")." << endl;
				s = 0;
			}
		}

		double val = (p3.row + (p4.row - p3.row) * s - p1.row - (p2.row - p1.row) * s);//s weight
		double t = (point.row - p1.row - (p2.row - p1.row) * s) / val;
		double valEpsilon = 0.01; // 0.1 and 0.01 appear identical
		if (fabs(val) < valEpsilon)
		{
			// Py ~= Cy because Dy - Cy ~= 0. So, instead of interpolating with y, we use x.
			t = (point.col - p1.col - (p2.col - p1.col) * s) / (p3.col + (p4.col - p3.col) * s - p1.col - (p2.col - p1.col) * s);
		}

		BilinearWeights weights;
		weights.s = s;
		weights.t = t;
		return weights;
	}
	//只有上下边界平行或者都不平行
	else 
	{
		// method 1
		double a = (p3.col - p1.col) * (p4.row - p2.row) - (p3.row - p1.row) * (p4.col - p2.col);
		double b = point.row * ((p4.col - p2.col) - (p3.col - p1.col)) - point.col * ((p4.row - p2.row) - (p3.row - p1.row)) + (p3.col - p1.col) * (p2.row) - (p3.row - p1.row) * (p2.col) + (p1.col) * (p4.row - p2.row) - (p1.row) * (p4.col - p2.col);
		double c = point.row * (p2.col - p1.col) - (point.col) * (p2.row - p1.row) + p1.col * p2.row - p2.col * p1.row;

		double t1 = (-1 * b + sqrt(b * b - 4 * a * c)) / (2 * a);
		double t2 = (-1 * b - sqrt(b * b - 4 * a * c)) / (2 * a);
		double t;
		//cout << "s1: " << t1 << " s2: " << t2 << endl;
		if (t1 >= 0 && t1 <= 1) 
		{
			t = t1;
		}
		else if (t2 >= 0 && t2 <= 1) 
		{
			t = t2;
		}
		else 
		{
			if ((t1 > 1 && t1 - quadraticEpsilon < 1) || (t2 > 1 && t2 - quadraticEpsilon < 1)) 
			{
				t = 1;
			}
			else if ((t1 < 0 && t1 + quadraticEpsilon > 0) || (t2 < 0 && t2 + quadraticEpsilon > 0)) 
			{
				t = 0;
			}
			else 
			{
				// this case should not happen
				//cout << "Could not interpolate t weight for coordinate (" << point.col << "," << point.row << ")." << endl;
				t = 0;
			}
		}

		double val = (p2.row + (p4.row - p2.row) * t - p1.row - (p3.row - p1.row) * t);
		double s = (point.row - p1.row - (p3.row - p1.row) * t) / val;
		double valEpsilon = 0.1; // 0.1 and 0.01 appear identical
		if (fabs(val) < valEpsilon) 
		{
			// Py ~= Ay because By - Ay ~= 0. So, instead of interpolating with y, we use x.
			s = (point.col - p1.col - (p3.col - p1.col) * t) / (p2.col + (p4.col - p2.col) * t - p1.col - (p3.col - p1.col) * t);
		}
		BilinearWeights weights;
		weights.s = clamp(s, 0, 1);
		weights.t = clamp(t, 0, 1);
		return weights;
	}
}

/*
	获取边界能量矩阵
	用来保证边界的形状一定是矩形，计算上下左右四个边界上点的能量
*/
pair<SpareseMatrixD_Row, VectorXd> GlobalWarpping::get_boundary_mat(vector<vector<CoordinateDouble>> mesh)
{
	//Vq=[x0 y0,x1,y1...]
	int rows = config.rows;
	int cols = config.cols;
	int numMeshRow = config.meshNumRow;//网格的行数
	int numMeshCol = config.meshNumCol;//网格的列数
	int vertexnum = numMeshRow * numMeshCol;//顶点的数量为行数乘以列数，(400)

	VectorXd dvec = VectorXd::Zero(double(vertexnum) * 2);//创建vertexnum*2维的向量devc
	VectorXd B = VectorXd::Zero(double(vertexnum) * 2);//B向量
	for (int i = 0; i < vertexnum * 2; i += numMeshCol * 2)
	{
		//左
		dvec(i) = 1;
		B(i) = 0;
	}
	for (int i = numMeshCol * 2 - 2; i < vertexnum * 2; i += numMeshCol * 2) 
	{
		//右
		dvec(i) = 1;
		B(i) = double(cols) - 1;
	}
	for (int i = 1; i < 2 * numMeshCol; i += 2) 
	{
		//上
		dvec(i) = 1;
		B(i) = 0;
	}
	for (int i = 2 * vertexnum - 2 * numMeshCol + 1; i < vertexnum * 2; i += 2) 
	{
		//下
		dvec(i) = 1;
		B(i) = double(rows) - 1;
	}

	//diag sparse;
	SpareseMatrixD_Row diag(dvec.size(), dvec.size());//用稀疏矩阵来存储对角元素
	for (int i = 0; i < dvec.size(); i++)
		diag.insert(i, i) = dvec(i);
	diag.makeCompressed();//对稀疏矩阵进行压缩
	return make_pair(diag, B);
};

/*
	获取mesh中位于(row, col)的四邻域网格顶点坐标坐标
*/
VectorXd GlobalWarpping::get_vertices(int _row, int _col, vector<vector<CoordinateDouble>>& mesh) 
{	
	//y0,x0,y1,x1...
	double row = _row;
	double col = _col;
	VectorXd Vq = VectorXd::Zero(8);//长度为8的向量
	//四边形网格的四个顶点坐标(左上角为row, col)
	CoordinateDouble p0 = mesh[row][col];//左上
	CoordinateDouble p1 = mesh[row][col + 1];//右上
	CoordinateDouble p2 = mesh[row + 1][col];//左下
	CoordinateDouble p3 = mesh[row + 1][col + 1];//右下
	Vq << p0.col, p0.row, p1.col, p1.row, p2.col, p2.row, p3.col, p3.row;//返回VectorXd Vq
	return Vq;
}

/*
	获取Shape能量矩阵
	保证每个四边形进行相似的变换，引自chengmingming论文
*/
SpareseMatrixD_Row GlobalWarpping::get_shape_mat(vector<vector<CoordinateDouble>> mesh) 
{
	double numMeshRow = config.meshNumRow;
	double numMeshCol = config.meshNumCol;
	double numQuadRow = config.meshQuadRow;//meshNumRow - 1
	double numQuadCol = config.meshQuadCol;//meshNumCol - 1
	//Shape能量矩阵,
	//用一个大的稀疏矩阵来存储shape energy,相当于每个单元为8*8的coeff，行数列数为numQuadRow * numQuadCol
	SpareseMatrixD_Row Shape_energy(8.0 * numQuadRow * numQuadCol, 8.0 * numQuadRow * numQuadCol);
	for (int row = 0; row < numQuadRow; row++) 
	{
		for (int col = 0; col < numQuadCol; col++) 
		{
			double _row = row;
			double _col = col;
			//网格四边形的四个顶点
			CoordinateDouble p0 = mesh[_row][_col];//左上
			CoordinateDouble p1 = mesh[_row][_col + 1];//右上
			CoordinateDouble p2 = mesh[_row + 1][_col];//左下
			CoordinateDouble p3 = mesh[_row + 1][_col + 1];//右下
			MatrixXd Aq(8, 4);//论文中的Aq矩阵
			Aq << p0.col, -p0.row, 1, 0,
				  p0.row, p0.col, 0, 1,
				  p1.col, -p1.row, 1, 0,
				  p1.row, p1.col, 0, 1,
				  p2.col, -p2.row, 1, 0,
				  p2.row, p2.col, 0, 1,
				  p3.col, -p3.row, 1, 0,
				  p3.row, p3.col, 0, 1;

			//计算系数
			MatrixXd Aq_trans = Aq.transpose();//Aq^T
			MatrixXd Aq_trans_mul_Aq_reverse = (Aq_trans * Aq).inverse();//(Aq^T * Aq)^-1
			MatrixXd I = MatrixXd::Identity(8, 8);//单位矩阵
			MatrixXd coeff = (Aq * Aq_trans_mul_Aq_reverse * Aq_trans - I);//矩阵规模为8*8

			int left_top_x = (row * numQuadCol + col) * 8;//大稀疏矩阵的左上角坐标
			//将coeff矩阵保存到大稀疏矩阵中
			for (int i = 0; i < 8; i++) 
			{
				for (int j = 0; j < 8; j++) 
				{
					Shape_energy.insert(double(left_top_x) + i, double(left_top_x) + j) = coeff(i, j);
				}
			}
		}
	}
	Shape_energy.makeCompressed();//对稀疏矩阵进行压缩
	return Shape_energy;
}

/*
	获取Q矩阵?
*/
SpareseMatrixD_Row GlobalWarpping::get_vertex_to_shape_mat(vector<vector<CoordinateDouble>> mesh)
{
	int numMeshRow = config.meshNumRow;
	int numMeshCol = config.meshNumCol;
	int numQuadRow = config.meshQuadRow;
	int numQuadCol = config.meshQuadCol;
	//大稀疏矩阵,行数和列数都为numQuadRow * numQuadCol,每个元素为8*2的矩阵
	SpareseMatrixD_Row Q(8.0 * numQuadRow * numQuadCol, 2.0 * numMeshRow * numMeshCol);
	for (int row = 0; row < numQuadRow; row++)
	{
		for (int col = 0; col < numQuadCol; col++)
		{
			int quadid = 8 * (row * numQuadCol + col);//大矩阵的左上角行坐标
			int topleftvertexId = 2 * (row * numMeshCol + col);//大矩阵的左上角列坐标
			Q.insert(quadid, topleftvertexId) = 1;
			Q.insert(quadid + 1.0, topleftvertexId + 1.0) = 1;
			Q.insert(quadid + 2.0, topleftvertexId + 2.0) = 1;
			Q.insert(quadid + 3.0, topleftvertexId + 3.0) = 1;
			Q.insert(quadid + 4.0, topleftvertexId + 2.0 * numMeshCol) = 1;
			Q.insert(quadid + 5.0, topleftvertexId + 2.0 * numMeshCol + 1.0) = 1;
			Q.insert(quadid + 6.0, topleftvertexId + 2.0 * numMeshCol + 2.0) = 1;
			Q.insert(quadid + 7.0, topleftvertexId + 2.0 * numMeshCol + 3.0) = 1;
		}
	}
	Q.makeCompressed();//Q矩阵进行压缩
	
	/*for (int k = 0; k < Q.outerSize(); ++k)
	{
		for (SparseMatrix<double>::InnerIterator it(Q, k); it; ++it)
		{
			cout << it.row() << " " << it.col() << " : " << it.value() << endl;
		}
		cout << std::endl;
	}*/
	
	return Q;
}

/*
	判断是否满足X0<a<X1或者X1<a<X0，同时距离小于1e-8
*/
bool between(double a, double X0, double X1) 
{
	double temp1 = a - X0;
	double temp2 = a - X1;
	if ((temp1<1e-8 && temp2>-1e-8) || (temp2<1e-6 && temp1>-1e-8))
		return true;
	else
		return false;
}

/*
	检测两条直线是否相交,如果是,返回两条直线的交点
	如果相交,将isintersection设为true,否则为false
*/
Vector2d GlobalWarpping::detectIntersect(Matrix2d line1, Matrix2d line2, bool& isintersection) 
{
	double line_x = 0, line_y = 0; //交点坐标
	//第一条直线的两个端点坐标
	double p1_x = line1(0, 1), p1_y = line1(0, 0), p2_x = line1(1, 1), p2_y = line1(1, 0);
	//第二条直线的两个端点坐标
	double p3_x = line2(0, 1), p3_y = line2(0, 0), p4_x = line2(1, 1), p4_y = line2(1, 0);

	//如果两条直线均平行于y轴,则不相交
	if ((fabs(p1_x - p2_x) < 1e-6) && (fabs(p3_x - p4_x) < 1e-6)) 
	{
		isintersection = false;
	}
	//只有第一条直线平行于y轴
	else if ((fabs(p1_x - p2_x) < 1e-6)) 
	{
		//如果p1_x位于p3_x与p4_x之间,则有可能相交,否则必然无交点
		if (between(p1_x, p3_x, p4_x)) 
		{
			double k = (p4_y - p3_y) / (p4_x - p3_x);//第二条直线斜率
			line_x = p1_x;//交点x坐标
			line_y = k * (line_x - p3_x) + p3_y;//交点y坐标
			//如果交点位于p1_y与p2_y之间,则说明相交
			if (between(line_y, p1_y, p2_y)) 
				isintersection = true;
			else 
				isintersection = false;
		}
		else 
			isintersection = false;
	}
	//只有第二条直线平行于y轴
	else if ((fabs(p3_x - p4_x) < 1e-6)) 
	{ 
		//如果p3_x位于p1_x与p2_x之间,则有可能相交,否则必然无交点
		if (between(p3_x, p1_x, p2_x)) 
		{
			double k = (p2_y - p1_y) / (p2_x - p1_x);//第一条直线斜率
			line_x = p3_x;//交点x坐标
			line_y = k * (line_x - p2_x) + p2_y;//交点y坐标
			if (between(line_y, p3_y, p4_y)) 
				isintersection = true;
			else 
				isintersection = false;
		}
		else 
			isintersection = false;
	}
	//两条直线都不平行与y轴
	else 
	{
		double k1 = (p2_y - p1_y) / (p2_x - p1_x);//第一条直线斜率
		double k2 = (p4_y - p3_y) / (p4_x - p3_x);//第二条直线斜率
		//如果两条直线斜率相同,则必无交点
		if (fabs(k1 - k2) < 1e-6)
			isintersection = false;
		else 
		{
			line_x = ((p3_y - p1_y) - (k2 * p3_x - k1 * p1_x)) / (k1 - k2);//"直线"对应交点x坐标
			line_y = k1 * (line_x - p1_x) + p1_y;//"直线"对应交点y坐标
			if (between(line_x, p1_x, p2_x) && between(line_x, p3_x, p4_x))
				isintersection = true;
			else
				isintersection = false;
		}
	}
	Vector2d p;
	p << line_y, line_x;
	return p;
}

/*
	修改mask边缘
*/
void GlobalWarpping::revise_mask_for_lines(CVMat& mask)
{
	//边缘的检测不予关注
	//对mask腐蚀
	int rows = mask.rows;
	int cols = mask.cols;
	//将左右两边缘设置为空
	for (int row = 0; row < rows; row++)
	{
		mask.at<uchar>(row, 0) = 255;
		mask.at<uchar>(row, cols - 1) = 255;
	}
	//将上下两边缘设置为空
	for (int col = 0; col < cols; col++)
	{
		mask.at<uchar>(0, col) = 255;
		mask.at<uchar>(rows - 1, col) = 255;
	}
	CVMat element = cv::getStructuringElement(cv::MORPH_ELLIPSE, cv::Size(8, 8));//椭圆形
	cv::dilate(mask, mask, element);//膨胀
	cv::dilate(mask, mask, element);//膨胀
}

/*
	判断点point是否在topLeft,topRight,bottomLeft,bottomRight组成地四边形内部
*/
bool GlobalWarpping::is_in_quad(CoordinateDouble point, CoordinateDouble topLeft, CoordinateDouble topRight, CoordinateDouble bottomLeft, CoordinateDouble bottomRight) 
{
	//点必须在topleft右下，topRight左下，bottomLeft右上，bottomRight左上

	//点必须在left右侧
	if (topLeft.col == bottomLeft.col) 
	{
		if (point.col < topLeft.col)
		{
			return false;
		}
	}
	else
	{
		//判断是否和四边形左边界有交点
		double leftSlope = (topLeft.row - bottomLeft.row) / (topLeft.col - bottomLeft.col);//左边界斜率
		double leftIntersect = topLeft.row - leftSlope * topLeft.col;//b
		double yOnLineX = (point.row - leftIntersect) / leftSlope;//(y-b)/k
		if (point.col < yOnLineX)
		{
			return false;
		}
	}
	//点必须在right左侧
	if (topRight.col == bottomRight.col) 
	{
		if (point.col > topRight.col)
		{
			return false;
		}
	}
	else 
	{
		//判断是否和四边形右边界有交点
		double rightSlope = (topRight.row - bottomRight.row) / (topRight.col - bottomRight.col);
		double rightIntersect = topRight.row - rightSlope * topRight.col;
		double yOnLineX = (point.row - rightIntersect) / rightSlope;
		if (point.col > yOnLineX)
		{
			return false;
		}
			
	}
	// must be below top line
	if (topLeft.row == topRight.row) 
	{
		if (point.row < topLeft.row)
		{
			return false;
		}
	}
	else 
	{
		//判断是否和四边形上边界有交点
		double topSlope = (topRight.row - topLeft.row) / (topRight.col - topLeft.col);
		double topIntersect = topLeft.row - topSlope * topLeft.col;
		double xOnLineY = topSlope * point.col + topIntersect;
		if (point.row < xOnLineY)
		{
			return false;
		}
	}
	// must be above bottom line
	if (bottomLeft.row == bottomRight.row) 
	{
		if (point.row > bottomLeft.row)
		{
			return false;
		}
	}
	else 
	{
		//判断是否和四边形下边界有交点
		double bottomSlope = (bottomRight.row - bottomLeft.row) / (bottomRight.col - bottomLeft.col);
		double bottomIntersect = bottomLeft.row - bottomSlope * bottomLeft.col;
		double xOnLineY = bottomSlope * point.col + bottomIntersect;
		if (point.row > xOnLineY)
		{
			return false;
		}
	}
	//如果上述条件均不满足,则说明点在四边形内部
	return true;
}

/*
	判断线段Line是否在mask中
*/
bool GlobalWarpping::line_in_mask(CVMat mask, LineD line)
{
	//如果线段起点和终点都在mask中,则说明线段在mask中
	if (mask.at<uchar>(line.row1, line.col1) == 0 && mask.at<uchar>(line.row2, line.col2) == 0) 
	{
		return true;
	}
	return false;
}

/*
	调用lsd方法进行直线检测
*/
vector<LineD> GlobalWarpping::lsd_detect(CVMat mask)
{
	int rows = src.rows;
	int cols = src.cols;
	CVMat gray_img;
	cv::cvtColor(src, gray_img, CV_BGR2GRAY);//灰度图
	double* image = new double[double(gray_img.rows) * double(gray_img.cols)];//用double数组存储gray_img(lsd要求)
	for (int row = 0; row < gray_img.rows; row++) 
		for (int col = 0; col < gray_img.cols; col++) 
			image[row * gray_img.cols + col] = gray_img.at<uchar>(row, col);

	vector<LineD> lines;
	double* out;//输出
	int num_lines;//检测到线地个数个数
	out = lsd(&num_lines, image, gray_img.cols, gray_img.rows);//调用lsd进行直线检测
	for (int i = 0; i < num_lines; i++) 
	{
		LineD line(out[i * 7 + 1], out[i * 7 + 0], out[i * 7 + 3], out[i * 7 + 2]);
		//CoordinateDouble start(out[i * 7 + 1], out[i * 7 + 0]);
		//CoordinateDouble end(out[i * 7 + 3], out[i * 7 + 2]);
		if (line_in_mask(mask, line)) 
		{
			lines.push_back(line);
		}
		//DrawLine(src, start, end);
	}

	/*cv::namedWindow("Border");
	cv::imshow("Border", src);
	cv::waitKey(0);*/

	return lines;
}

/*
	根据网格的k和b来判断网格边缘是否和直线相交
*/
bool GlobalWarpping::does_segment_intersect_line(LineD lineSegment, double slope, double intersect, bool vertical, CoordinateDouble& intersectPoint)
{
	double lineSegmentSlope = INF;
	if (lineSegment.col1 != lineSegment.col2) 
	{
		//如果这条线段不平行于y轴
		lineSegmentSlope = (lineSegment.row2 - lineSegment.row1) / (lineSegment.col2 - lineSegment.col1);//计算线段的斜率
	}
	double lineSegmentIntersect = lineSegment.row1 - lineSegmentSlope * lineSegment.col1;//b = y - kx
	//如果线段与网格线平行
	if (lineSegmentSlope == slope)
	{
		if (lineSegmentIntersect == intersect)
		{
			//如果线段的b等于网格线的b，说明这两条线段一样
			intersectPoint.col = lineSegment.col1;//交点横坐标即为线段横坐标
			intersectPoint.row = lineSegment.row1;//交点纵坐标即为线段纵坐标
			return true;
		}
		else 
		{
			return false;
		}
	}

	//如果线段与网格线不平行
	//交点x坐标 x = (b2 - b1) / (k1 - k2)
	double intersectX = (intersect - lineSegmentIntersect) / (lineSegmentSlope - slope);
	//交点y坐标 y = k1 * x + b
	double intersectY = lineSegmentSlope * intersectX + lineSegmentIntersect;
	//检查交点是否在线段上
	if (vertical)
	{
		//垂直，检查纵坐标，判断是否相交
		if ((intersectY <= lineSegment.row1 && intersectY >= lineSegment.row2) || (intersectY <= lineSegment.row2 && intersectY >= lineSegment.row1)) 
		{
			intersectPoint.col = intersectX;
			intersectPoint.row = intersectY;
			return true;
		}
		else 
		{
			return false;
		}
	}
	else 
	{
		//水平，检查横坐标，判断是否相交
		if ((intersectX <= lineSegment.col1 && intersectX >= lineSegment.col2) || (intersectX <= lineSegment.col2 && intersectX >= lineSegment.col1)) 
		{
			intersectPoint.col = intersectX;
			intersectPoint.row = intersectY;
			return true;
		}
		else 
		{
			return false;
		}
	}
}

/*
	用四边形网格对线段进行分割，返回交点集
	通过调用does_segment_intersect_line实现
*/
vector<CoordinateDouble> GlobalWarpping::intersections_with_quad(LineD lineSegment, CoordinateDouble topLeft, CoordinateDouble topRight, CoordinateDouble bottomLeft, CoordinateDouble bottomRight)
{
	vector<CoordinateDouble> intersections;//交点集

	//左
	double leftSlope = INF;
	if (topLeft.col != bottomLeft.col) 
	{
		leftSlope = (topLeft.row - bottomLeft.row) / (topLeft.col - bottomLeft.col);//左边界k
	}
	double leftIntersect = topLeft.row - leftSlope * topLeft.col;//左边界b
	//判断线段是否和该边界相交
	CoordinateDouble leftIntersectPoint;//存放相交的点
	if (does_segment_intersect_line(lineSegment, leftSlope, leftIntersect, true, leftIntersectPoint)) 
	{
		if (leftIntersectPoint.row >= topLeft.row && leftIntersectPoint.row <= bottomLeft.row) 
		{
			intersections.push_back(leftIntersectPoint);
		}
	}

	//右
	double rightSlope = INF;
	if (topRight.col != bottomRight.col) 
	{
		rightSlope = (topRight.row - bottomRight.row) / (topRight.col - bottomRight.col);//右边界k
	}
	double rightIntersect = topRight.row - rightSlope * topRight.col;//右边界b
	//判断线段是否和该边界相交
	CoordinateDouble rightIntersectPoint;
	if (does_segment_intersect_line(lineSegment, rightSlope, rightIntersect, true, rightIntersectPoint)) 
	{
		if (rightIntersectPoint.row >= topRight.row && rightIntersectPoint.row <= bottomRight.row) 
		{
			intersections.push_back(rightIntersectPoint);
		}
	}

	//上
	double topSlope = INF;
	if (topLeft.col != topRight.col) 
	{
		topSlope = (topRight.row - topLeft.row) / (topRight.col - topLeft.col);//上边界k
	}
	double topIntersect = topLeft.row - topSlope * topLeft.col;//上边界b
	//判断线段是否和该边界相交
	CoordinateDouble topIntersectPoint;
	if (does_segment_intersect_line(lineSegment, topSlope, topIntersect, false, topIntersectPoint)) 
	{
		if (topIntersectPoint.col >= topLeft.col && topIntersectPoint.col <= topRight.col) 
		{
			intersections.push_back(topIntersectPoint);
		}
	}

	//下
	double bottomSlope = INF;
	if (bottomLeft.col != bottomRight.col) 
	{
		bottomSlope = (bottomRight.row - bottomLeft.row) / (bottomRight.col - bottomLeft.col);//下边界k
	}
	double bottomIntersect = bottomLeft.row - bottomSlope * bottomLeft.col;//下边界b
	//判断线段是否和该边界相交
	CoordinateDouble bottomIntersectPoint;
	if (does_segment_intersect_line(lineSegment, bottomSlope, bottomIntersect, false, bottomIntersectPoint)) 
	{
		if (bottomIntersectPoint.col >= bottomLeft.col && bottomIntersectPoint.col <= bottomRight.col)
		{
			intersections.push_back(bottomIntersectPoint);
		}
	}
	return intersections;//返回交点集
}

/*
	在整个网格中对所有直线进行分割
*/
vector<vector<vector<LineD>>> GlobalWarpping::segment_line_in_quad(vector<LineD> lines, vector<vector<CoordinateDouble>> mesh)
{
	int QuadnumRow = config.meshQuadRow;
	int QuadnumCol = config.meshQuadCol;
	vector<vector<vector<LineD>>> quad_line_seg;//用一个三维vector保存分割结果,分别为行,列,每个网格内线段集
	CVMat src2;
	src.copyTo(src2);

	for (int row = 0; row < QuadnumRow; row++) 
	{
		vector<vector<LineD>> vec_row;
		for (int col = 0; col < QuadnumCol; col++) 
		{
			//每一个网格[row][col]
			//获取网格的四个顶点
			CoordinateDouble lefttop = mesh[double(row)][double(col)];
			CoordinateDouble righttop = mesh[double(row)][double(col) + 1];
			CoordinateDouble leftbottom = mesh[double(row) + 1][double(col)];
			CoordinateDouble rightbottom = mesh[double(row) + 1][double(col) + 1];

			vector<LineD> lineInQuad;//存储网格对所有线段分割的结果
			//遍历每一条直线
			for (int i = 0; i < lines.size(); i++) 
			{
				LineD line = lines[i];
				CoordinateDouble point1(line.row1, line.col1);//线段端点1
				CoordinateDouble point2(line.row2, line.col2);//线段端点2
				bool p1InQuad = is_in_quad(point1, lefttop, righttop, leftbottom, rightbottom);
				bool p2InQuad = is_in_quad(point2, lefttop, righttop, leftbottom, rightbottom);
				if (p1InQuad && p2InQuad) 
				{
					//如果p1,p2都在四边形中,则线段必在四边形中
					lineInQuad.push_back(line);
				}
				else if (p1InQuad) 
				{
					//只有p1在四边形中
					vector<CoordinateDouble> intersections = intersections_with_quad(line, lefttop, righttop, leftbottom, rightbottom);
					//如果有交点,则对应线段两个端点分别为point1, intersections[0]
					if (intersections.size() != 0) 
					{
						LineD cutLine(point1, intersections[0]);
						lineInQuad.push_back(cutLine);
					}
				}
				else if (p2InQuad) 
				{
					//只有p2在四边形中
					vector<CoordinateDouble> intersections = intersections_with_quad(line, lefttop, righttop, leftbottom, rightbottom);
					if (intersections.size() != 0) 
					{
						LineD cutLine(point2, intersections[0]);
						lineInQuad.push_back(cutLine);
					}
				}
				else 
				{
					//两个端点都不在四边形中
					vector<CoordinateDouble> intersections = intersections_with_quad(line, lefttop, righttop, leftbottom, rightbottom);
					//如果有两个交点,则切割线段两个端点分别为intersections[0], intersections[1]
					if (intersections.size() == 2)
					{
						LineD cutLine(intersections[0], intersections[1]);
						lineInQuad.push_back(cutLine);
					}
				}
			}
			vec_row.push_back(lineInQuad);
			//TEST
			//DRAW BORDER
			//DrawLine(src, lefttop, righttop);
			//DrawLine(src, lefttop, leftbottom);
			//DrawLine(src, righttop, rightbottom);
			//DrawLine(src, leftbottom, rightbottom);
			/*for (int i = 0; i < lineInQuad.size(); i++) {
				LineD line = lineInQuad[i];
				DrawLine(src2, line);
			}
			cv::namedWindow("quad", CV_WINDOW_AUTOSIZE);
			cv::imshow("quad", src);
			cv::namedWindow("line", CV_WINDOW_AUTOSIZE);
			cv::imshow("line", src2);

			cv::waitKey(0);
			*/
		}
		quad_line_seg.push_back(vec_row);
	}
	return quad_line_seg;
}

/*
	将三维的线段分割结果展开成一维的vector<LineD>
*/
vector<LineD> GlobalWarpping::flatten(vector<vector<vector<LineD>>>& lineSeg)
{
	vector<LineD> line_vec;//结果vector
	int numQuadRow = config.meshQuadRow;
	int numQuadCol = config.meshQuadCol;
	for (int row = 0; row < numQuadRow; row++) 
	{
		for (int col = 0; col < numQuadCol; col++) 
		{
			for (int k = 0; k < lineSeg[row][col].size(); k++) 
			{
				line_vec.push_back(lineSeg[row][col][k]);
			}
		}
	}
	return line_vec;
}

/*
	在origin矩阵上进行扩展,即矩阵的合并,返回扩展之后的矩阵
*/
SpareseMatrixD_Row GlobalWarpping::block_diag_extend(SpareseMatrixD_Row origin, MatrixXd addin, int QuadID)
{
	int cols_total = 8 * config.meshQuadRow * config.meshQuadCol;//总的列数,每一个网格单元为8*8(Shape Preservation Mat)
	SpareseMatrixD_Row res(origin.rows() + addin.rows(), cols_total);//按行扩展,矩阵上半部分为origin,下半部分为addin
	res.topRows(origin.rows()) = origin;//矩阵上半部分

	int lefttop_row = origin.rows();//addin开始的行数
	int lefttop_col = 8 * QuadID;//addin开始的列数
	for (int row = 0; row < addin.rows(); row++) 
	{
		for (int col = 0; col < addin.cols(); col++) 
		{
			res.insert(lefttop_row + row, lefttop_col + col) = addin(row, col);//将addin矩阵插入到res矩阵
		}
	}
	res.makeCompressed();//对稀疏矩阵进行压缩
	return res;
}

/*
	初始化线段分割,将具有相近倾斜角度的线段分配到一个集合中
*/
vector<vector<vector<LineD>>> GlobalWarpping::init_line_seg(CVMat mask, vector<LineD>& lineSeg_flatten, vector<vector<CoordinateDouble>> mesh, vector<pair<int,double>>& id_theta, vector<double>& rotate_theta)
{
	double thetaPerbin = PI / 49;
	revise_mask_for_lines(mask);
	//第一步:使用LSD进行直线检测
	vector<LineD> lines = lsd_detect(mask);
	//第二步:用网格分割所有线段到各自的网格内部
	vector<vector<vector<LineD>>> lineSeg = segment_line_in_quad(lines, mesh);
	//第三步:根据每个线段的倾斜角度,将线分配到预先设置的格子内部
	lineSeg_flatten = flatten(lineSeg);//将得到的三维vector展开
	for (int i = 0; i < lineSeg_flatten.size(); i++) 
	{
		LineD line = lineSeg_flatten[i];//每条直线
		double theta = atan((line.row1 - line.row2) / (line.col1 - line.col2));//计算线段的倾斜角度arctan(y1-y1)/(x1-x2)
		int lineSegmentBucket = (int)round((theta + PI / 2) / thetaPerbin);//计算线段对应的bin id
		assert(lineSegmentBucket < 50);//保证lineSegmentBucket < 50
		id_theta.push_back(make_pair(lineSegmentBucket, theta));//将该线段加入到id_theta中
		rotate_theta.push_back(0);//
	}
	return lineSeg;
}

/*
	获取直线能量矩阵
*/
SpareseMatrixD_Row GlobalWarpping::get_line_mat(CVMat mask, vector<vector<CoordinateDouble>> mesh, vector<double>rotate_theta, vector<vector<vector<LineD>>> lineSeg, vector<pair<MatrixXd, MatrixXd>>& BilinearVec, int& linenum, vector<bool>& bad) 
{
	int linetmpnum = -1;
	int rows = config.rows;
	int cols = config.cols;
	int QuadnumRow = config.meshQuadRow;
	int QuadnumCol = config.meshQuadCol;
	double gridcols = config.colPermesh;//每个网格所占的列数
	double gridrows = config.rowPermesh;//每个网格所占的行数

	/*
	for (int i = 0; i < lines.size(); i++) {
		LineD line = lines[i];

		CoordinateDouble start(line.row1, line.col1);
		CoordinateDouble end(line.row2,line.col2 );
		if (line_in_mask(mask, line)) {
			DrawLine(src, start, end);
			cv::namedWindow("Border", CV_WINDOW_AUTOSIZE);
			cv::imshow("Border", src);
			cv::waitKey(0);
		}
	}*/
	SpareseMatrixD_Row energy_line;//line能量矩阵
	for (int row = 0; row < QuadnumRow; row++) 
	{
		for (int col = 0; col < QuadnumCol; col++) 
		{
			//遍历每一个网格
			vector<LineD> linesegInquad = lineSeg[row][col];//网格mesh[row][col]中分割得到的线段集合
			int QuadID = row * QuadnumCol + col;
			if (linesegInquad.size() == 0) 
			{
				//如果当前网格没有线段,则直接下一个网格
				continue;
			}
			
			Coordinate topleft(row, col);//左上角顶点
			MatrixXd C_row_stack(0, 8);//用来保存每一个网格中所有线段对应的C矩阵

			for (int k = 0; k < linesegInquad.size(); k++) 
			{
				//遍历每一条线段
				linetmpnum++;//记录当前线的数量,用来进行索引
				LineD line = linesegInquad[k];//第k条线
				CoordinateDouble linestart(line.row1, line.col1);//线段起点
				CoordinateDouble lineend(line.row2, line.col2);//线段终点

				BilinearWeights startWeight = get_bilinear_weights(linestart, topleft, mesh);//s t2n t t1n
				MatrixXd start_W_mat = BilinearWeightsToMatrix(startWeight);
				BilinearWeights endWeight = get_bilinear_weights(lineend, topleft, mesh);
				MatrixXd end_W_mat = BilinearWeightsToMatrix(endWeight);
				//cout << startWeight.s << " " << startWeight.t << endl;//test
				//test
				VectorXd S = get_vertices(row, col, mesh);
				Vector2d ans = start_W_mat * S - Vector2d(linestart.col, linestart.row);
				Vector2d ans2 = end_W_mat * S - Vector2d(lineend.col, lineend.row);

				//如果ans2的2范数>=0.0001或者ans的2范数>=0.0001
				if (ans2.norm() >= 0.0001 || ans.norm() >= 0.0001) 
				{
					//错误情况
					bad.push_back(true);
					BilinearVec.push_back(make_pair(MatrixXd::Zero(2, 8), MatrixXd::Zero(2, 8)));
					continue;
				}
				bad.push_back(false);
				//assert(ans.norm() < 0.0001);
				//assert(ans2.norm() < 0.0001);
				
				//end test
				//system("pause");

				double theta = rotate_theta[linetmpnum];//论文中的theta_m
				BilinearVec.push_back(make_pair(start_W_mat, end_W_mat));
				Matrix2d R;//论文中的旋转矩阵R
				R << cos(theta), -sin(theta),
					sin(theta), cos(theta);
				MatrixXd ehat(2, 1);//ehat
				ehat << line.col1 - line.col2, line.row1 - line.row2;//计算ehat
				MatrixXd tmp = (ehat.transpose() * ehat).inverse();//(ehat^T * ehat)^-1
				Matrix2d I = Matrix2d::Identity();//单位矩阵
				MatrixXd C = R * ehat * tmp * (ehat.transpose()) * (R.transpose()) - I;//论文(6)中的C
				//线段方向向量e = start_W_mat - end_W_mat,
				MatrixXd CT = C * (start_W_mat - end_W_mat); //C*e
				C_row_stack = row_stack(C_row_stack, CT);//将CT加到C_row_stack下方,按行合并
			}
			energy_line = block_diag_extend(energy_line, C_row_stack, QuadID);//将C_row_stack加到energy_line矩阵下方,按行扩展
		}
	}
	linenum = linetmpnum;
	return energy_line;
}

