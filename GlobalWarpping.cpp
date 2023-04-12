#include"GlobalWarpping.h"
#pragma warning(disable: 4244)

/*
	Ĭ�Ϲ��캯��,��ʼ��config
*/
GlobalWarpping::GlobalWarpping(Config& conf, CVMat& src)
{
	this->config = conf;//��ʼ��config
	this->src = src;//��ʼ��src
}

/*
	��BilinearWeightsת��ΪMatrixXd
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
	��ȡ˫����Ȩֵ
*/
BilinearWeights GlobalWarpping::get_bilinear_weights(CoordinateDouble point, Coordinate upperLeftIndices, vector<vector<CoordinateDouble>> mesh)
{
	//��ȡ�ı���������ĸ�����
	CoordinateDouble p1 = mesh[upperLeftIndices.row][upperLeftIndices.col];//���Ͻ�
	CoordinateDouble p2 = mesh[upperLeftIndices.row][upperLeftIndices.col + 1];//���Ͻ�
	CoordinateDouble p3 = mesh[upperLeftIndices.row + 1][upperLeftIndices.col];//���½�
	CoordinateDouble p4 = mesh[upperLeftIndices.row + 1][upperLeftIndices.col + 1];//���½�

	double slopeTop    = (p2.row - p1.row) / (p2.col - p1.col);//�ϱ߽�б��
	double slopeBottom = (p4.row - p3.row) / (p4.col - p3.col);//�±߽�б��
	double slopeLeft   = (p1.row - p3.row) / (p1.col - p3.col);//��߽�б��
	double slopeRight  = (p2.row - p4.row) / (p2.col - p4.col);//�ұ߽�б��

	double quadraticEpsilon = 0.1;//

	//���±߽�ƽ�в������ұ߽�ƽ��
	if (slopeTop == slopeBottom && slopeLeft == slopeRight)
	{
		// method 3
		Matrix2d mat1;
		//�ϱ߽�ĳ�,��߽�ĳ�,�ϱ߽�Ŀ�,��߽�Ŀ�
		mat1 << p2.col - p1.col, p3.col - p1.col,
			p2.row - p1.row, p3.row - p1.row;

		MatrixXd mat2(2, 1);
		mat2 << point.col - p1.col, point.row - p1.row;//point��p1�ĳ�,point��p1�Ŀ�

		MatrixXd matsolution = mat1.inverse() * mat2;//mat1^-1 * mat2

		BilinearWeights weights;
		weights.s = matsolution(0, 0);
		weights.t = matsolution(1, 0);
		return weights;
	}
	//ֻ�����ұ߽�ƽ��
	else if (slopeLeft == slopeRight)
	{
		//method 2
		//a = (x2-x1)*(y4-y3) - (y2-y1)*(x4-x3)
		double a = (p2.col - p1.col) * (p4.row - p3.row) - (p2.row - p1.row) * (p4.col - p3.col);
		//b = y*((x4-x3)-(x2-x1)) - x*((y4-y3)-(y2-y1)) + x1*(y4-y3) - y1*(x4-x3) + y3*(x2-x1)-x3*(y2-y1)
		double b = point.row * ((p4.col - p3.col) - (p2.col - p1.col)) - point.col * ((p4.row - p3.row) - (p2.row - p1.row)) + p1.col * (p4.row - p3.row) - p1.row * (p4.col - p3.col) + (p2.col - p1.col) * (p3.row) - (p2.row - p1.row) * (p3.col);
		//c = y*(x3-x1) - x*(y3-y1) + x1*y3 - x3*y1;
		double c = point.row * (p3.col - p1.col) - point.col * (p3.row - p1.row) + p1.col * p3.row - p3.col * p1.row;

		//�󷽳�ax^2 + bx + c = 0�ĸ�
		double s1 = (-1 * b + sqrt(b * b - 4 * a * c)) / (2 * a);//��һ����
		double s2 = (-1 * b - sqrt(b * b - 4 * a * c)) / (2 * a);//�ڶ�����
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
	//ֻ�����±߽�ƽ�л��߶���ƽ��
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
	��ȡ�߽���������
	������֤�߽����״һ���Ǿ��Σ��������������ĸ��߽��ϵ������
*/
pair<SpareseMatrixD_Row, VectorXd> GlobalWarpping::get_boundary_mat(vector<vector<CoordinateDouble>> mesh)
{
	//Vq=[x0 y0,x1,y1...]
	int rows = config.rows;
	int cols = config.cols;
	int numMeshRow = config.meshNumRow;//���������
	int numMeshCol = config.meshNumCol;//���������
	int vertexnum = numMeshRow * numMeshCol;//���������Ϊ��������������(400)

	VectorXd dvec = VectorXd::Zero(double(vertexnum) * 2);//����vertexnum*2ά������devc
	VectorXd B = VectorXd::Zero(double(vertexnum) * 2);//B����
	for (int i = 0; i < vertexnum * 2; i += numMeshCol * 2)
	{
		//��
		dvec(i) = 1;
		B(i) = 0;
	}
	for (int i = numMeshCol * 2 - 2; i < vertexnum * 2; i += numMeshCol * 2) 
	{
		//��
		dvec(i) = 1;
		B(i) = double(cols) - 1;
	}
	for (int i = 1; i < 2 * numMeshCol; i += 2) 
	{
		//��
		dvec(i) = 1;
		B(i) = 0;
	}
	for (int i = 2 * vertexnum - 2 * numMeshCol + 1; i < vertexnum * 2; i += 2) 
	{
		//��
		dvec(i) = 1;
		B(i) = double(rows) - 1;
	}

	//diag sparse;
	SpareseMatrixD_Row diag(dvec.size(), dvec.size());//��ϡ��������洢�Խ�Ԫ��
	for (int i = 0; i < dvec.size(); i++)
		diag.insert(i, i) = dvec(i);
	diag.makeCompressed();//��ϡ��������ѹ��
	return make_pair(diag, B);
};

/*
	��ȡmesh��λ��(row, col)�����������񶥵���������
*/
VectorXd GlobalWarpping::get_vertices(int _row, int _col, vector<vector<CoordinateDouble>>& mesh) 
{	
	//y0,x0,y1,x1...
	double row = _row;
	double col = _col;
	VectorXd Vq = VectorXd::Zero(8);//����Ϊ8������
	//�ı���������ĸ���������(���Ͻ�Ϊrow, col)
	CoordinateDouble p0 = mesh[row][col];//����
	CoordinateDouble p1 = mesh[row][col + 1];//����
	CoordinateDouble p2 = mesh[row + 1][col];//����
	CoordinateDouble p3 = mesh[row + 1][col + 1];//����
	Vq << p0.col, p0.row, p1.col, p1.row, p2.col, p2.row, p3.col, p3.row;//����VectorXd Vq
	return Vq;
}

/*
	��ȡShape��������
	��֤ÿ���ı��ν������Ƶı任������chengmingming����
*/
SpareseMatrixD_Row GlobalWarpping::get_shape_mat(vector<vector<CoordinateDouble>> mesh) 
{
	double numMeshRow = config.meshNumRow;
	double numMeshCol = config.meshNumCol;
	double numQuadRow = config.meshQuadRow;//meshNumRow - 1
	double numQuadCol = config.meshQuadCol;//meshNumCol - 1
	//Shape��������,
	//��һ�����ϡ��������洢shape energy,�൱��ÿ����ԪΪ8*8��coeff����������ΪnumQuadRow * numQuadCol
	SpareseMatrixD_Row Shape_energy(8.0 * numQuadRow * numQuadCol, 8.0 * numQuadRow * numQuadCol);
	for (int row = 0; row < numQuadRow; row++) 
	{
		for (int col = 0; col < numQuadCol; col++) 
		{
			double _row = row;
			double _col = col;
			//�����ı��ε��ĸ�����
			CoordinateDouble p0 = mesh[_row][_col];//����
			CoordinateDouble p1 = mesh[_row][_col + 1];//����
			CoordinateDouble p2 = mesh[_row + 1][_col];//����
			CoordinateDouble p3 = mesh[_row + 1][_col + 1];//����
			MatrixXd Aq(8, 4);//�����е�Aq����
			Aq << p0.col, -p0.row, 1, 0,
				  p0.row, p0.col, 0, 1,
				  p1.col, -p1.row, 1, 0,
				  p1.row, p1.col, 0, 1,
				  p2.col, -p2.row, 1, 0,
				  p2.row, p2.col, 0, 1,
				  p3.col, -p3.row, 1, 0,
				  p3.row, p3.col, 0, 1;

			//����ϵ��
			MatrixXd Aq_trans = Aq.transpose();//Aq^T
			MatrixXd Aq_trans_mul_Aq_reverse = (Aq_trans * Aq).inverse();//(Aq^T * Aq)^-1
			MatrixXd I = MatrixXd::Identity(8, 8);//��λ����
			MatrixXd coeff = (Aq * Aq_trans_mul_Aq_reverse * Aq_trans - I);//�����ģΪ8*8

			int left_top_x = (row * numQuadCol + col) * 8;//��ϡ���������Ͻ�����
			//��coeff���󱣴浽��ϡ�������
			for (int i = 0; i < 8; i++) 
			{
				for (int j = 0; j < 8; j++) 
				{
					Shape_energy.insert(double(left_top_x) + i, double(left_top_x) + j) = coeff(i, j);
				}
			}
		}
	}
	Shape_energy.makeCompressed();//��ϡ��������ѹ��
	return Shape_energy;
}

/*
	��ȡQ����?
*/
SpareseMatrixD_Row GlobalWarpping::get_vertex_to_shape_mat(vector<vector<CoordinateDouble>> mesh)
{
	int numMeshRow = config.meshNumRow;
	int numMeshCol = config.meshNumCol;
	int numQuadRow = config.meshQuadRow;
	int numQuadCol = config.meshQuadCol;
	//��ϡ�����,������������ΪnumQuadRow * numQuadCol,ÿ��Ԫ��Ϊ8*2�ľ���
	SpareseMatrixD_Row Q(8.0 * numQuadRow * numQuadCol, 2.0 * numMeshRow * numMeshCol);
	for (int row = 0; row < numQuadRow; row++)
	{
		for (int col = 0; col < numQuadCol; col++)
		{
			int quadid = 8 * (row * numQuadCol + col);//���������Ͻ�������
			int topleftvertexId = 2 * (row * numMeshCol + col);//���������Ͻ�������
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
	Q.makeCompressed();//Q�������ѹ��
	
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
	�ж��Ƿ�����X0<a<X1����X1<a<X0��ͬʱ����С��1e-8
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
	�������ֱ���Ƿ��ཻ,�����,��������ֱ�ߵĽ���
	����ཻ,��isintersection��Ϊtrue,����Ϊfalse
*/
Vector2d GlobalWarpping::detectIntersect(Matrix2d line1, Matrix2d line2, bool& isintersection) 
{
	double line_x = 0, line_y = 0; //��������
	//��һ��ֱ�ߵ������˵�����
	double p1_x = line1(0, 1), p1_y = line1(0, 0), p2_x = line1(1, 1), p2_y = line1(1, 0);
	//�ڶ���ֱ�ߵ������˵�����
	double p3_x = line2(0, 1), p3_y = line2(0, 0), p4_x = line2(1, 1), p4_y = line2(1, 0);

	//�������ֱ�߾�ƽ����y��,���ཻ
	if ((fabs(p1_x - p2_x) < 1e-6) && (fabs(p3_x - p4_x) < 1e-6)) 
	{
		isintersection = false;
	}
	//ֻ�е�һ��ֱ��ƽ����y��
	else if ((fabs(p1_x - p2_x) < 1e-6)) 
	{
		//���p1_xλ��p3_x��p4_x֮��,���п����ཻ,�����Ȼ�޽���
		if (between(p1_x, p3_x, p4_x)) 
		{
			double k = (p4_y - p3_y) / (p4_x - p3_x);//�ڶ���ֱ��б��
			line_x = p1_x;//����x����
			line_y = k * (line_x - p3_x) + p3_y;//����y����
			//�������λ��p1_y��p2_y֮��,��˵���ཻ
			if (between(line_y, p1_y, p2_y)) 
				isintersection = true;
			else 
				isintersection = false;
		}
		else 
			isintersection = false;
	}
	//ֻ�еڶ���ֱ��ƽ����y��
	else if ((fabs(p3_x - p4_x) < 1e-6)) 
	{ 
		//���p3_xλ��p1_x��p2_x֮��,���п����ཻ,�����Ȼ�޽���
		if (between(p3_x, p1_x, p2_x)) 
		{
			double k = (p2_y - p1_y) / (p2_x - p1_x);//��һ��ֱ��б��
			line_x = p3_x;//����x����
			line_y = k * (line_x - p2_x) + p2_y;//����y����
			if (between(line_y, p3_y, p4_y)) 
				isintersection = true;
			else 
				isintersection = false;
		}
		else 
			isintersection = false;
	}
	//����ֱ�߶���ƽ����y��
	else 
	{
		double k1 = (p2_y - p1_y) / (p2_x - p1_x);//��һ��ֱ��б��
		double k2 = (p4_y - p3_y) / (p4_x - p3_x);//�ڶ���ֱ��б��
		//�������ֱ��б����ͬ,����޽���
		if (fabs(k1 - k2) < 1e-6)
			isintersection = false;
		else 
		{
			line_x = ((p3_y - p1_y) - (k2 * p3_x - k1 * p1_x)) / (k1 - k2);//"ֱ��"��Ӧ����x����
			line_y = k1 * (line_x - p1_x) + p1_y;//"ֱ��"��Ӧ����y����
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
	�޸�mask��Ե
*/
void GlobalWarpping::revise_mask_for_lines(CVMat& mask)
{
	//��Ե�ļ�ⲻ���ע
	//��mask��ʴ
	int rows = mask.rows;
	int cols = mask.cols;
	//����������Ե����Ϊ��
	for (int row = 0; row < rows; row++)
	{
		mask.at<uchar>(row, 0) = 255;
		mask.at<uchar>(row, cols - 1) = 255;
	}
	//����������Ե����Ϊ��
	for (int col = 0; col < cols; col++)
	{
		mask.at<uchar>(0, col) = 255;
		mask.at<uchar>(rows - 1, col) = 255;
	}
	CVMat element = cv::getStructuringElement(cv::MORPH_ELLIPSE, cv::Size(8, 8));//��Բ��
	cv::dilate(mask, mask, element);//����
	cv::dilate(mask, mask, element);//����
}

/*
	�жϵ�point�Ƿ���topLeft,topRight,bottomLeft,bottomRight��ɵ��ı����ڲ�
*/
bool GlobalWarpping::is_in_quad(CoordinateDouble point, CoordinateDouble topLeft, CoordinateDouble topRight, CoordinateDouble bottomLeft, CoordinateDouble bottomRight) 
{
	//�������topleft���£�topRight���£�bottomLeft���ϣ�bottomRight����

	//�������left�Ҳ�
	if (topLeft.col == bottomLeft.col) 
	{
		if (point.col < topLeft.col)
		{
			return false;
		}
	}
	else
	{
		//�ж��Ƿ���ı�����߽��н���
		double leftSlope = (topLeft.row - bottomLeft.row) / (topLeft.col - bottomLeft.col);//��߽�б��
		double leftIntersect = topLeft.row - leftSlope * topLeft.col;//b
		double yOnLineX = (point.row - leftIntersect) / leftSlope;//(y-b)/k
		if (point.col < yOnLineX)
		{
			return false;
		}
	}
	//�������right���
	if (topRight.col == bottomRight.col) 
	{
		if (point.col > topRight.col)
		{
			return false;
		}
	}
	else 
	{
		//�ж��Ƿ���ı����ұ߽��н���
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
		//�ж��Ƿ���ı����ϱ߽��н���
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
		//�ж��Ƿ���ı����±߽��н���
		double bottomSlope = (bottomRight.row - bottomLeft.row) / (bottomRight.col - bottomLeft.col);
		double bottomIntersect = bottomLeft.row - bottomSlope * bottomLeft.col;
		double xOnLineY = bottomSlope * point.col + bottomIntersect;
		if (point.row > xOnLineY)
		{
			return false;
		}
	}
	//�������������������,��˵�������ı����ڲ�
	return true;
}

/*
	�ж��߶�Line�Ƿ���mask��
*/
bool GlobalWarpping::line_in_mask(CVMat mask, LineD line)
{
	//����߶������յ㶼��mask��,��˵���߶���mask��
	if (mask.at<uchar>(line.row1, line.col1) == 0 && mask.at<uchar>(line.row2, line.col2) == 0) 
	{
		return true;
	}
	return false;
}

/*
	����lsd��������ֱ�߼��
*/
vector<LineD> GlobalWarpping::lsd_detect(CVMat mask)
{
	int rows = src.rows;
	int cols = src.cols;
	CVMat gray_img;
	cv::cvtColor(src, gray_img, CV_BGR2GRAY);//�Ҷ�ͼ
	double* image = new double[double(gray_img.rows) * double(gray_img.cols)];//��double����洢gray_img(lsdҪ��)
	for (int row = 0; row < gray_img.rows; row++) 
		for (int col = 0; col < gray_img.cols; col++) 
			image[row * gray_img.cols + col] = gray_img.at<uchar>(row, col);

	vector<LineD> lines;
	double* out;//���
	int num_lines;//��⵽�ߵظ�������
	out = lsd(&num_lines, image, gray_img.cols, gray_img.rows);//����lsd����ֱ�߼��
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
	���������k��b���ж������Ե�Ƿ��ֱ���ཻ
*/
bool GlobalWarpping::does_segment_intersect_line(LineD lineSegment, double slope, double intersect, bool vertical, CoordinateDouble& intersectPoint)
{
	double lineSegmentSlope = INF;
	if (lineSegment.col1 != lineSegment.col2) 
	{
		//��������߶β�ƽ����y��
		lineSegmentSlope = (lineSegment.row2 - lineSegment.row1) / (lineSegment.col2 - lineSegment.col1);//�����߶ε�б��
	}
	double lineSegmentIntersect = lineSegment.row1 - lineSegmentSlope * lineSegment.col1;//b = y - kx
	//����߶���������ƽ��
	if (lineSegmentSlope == slope)
	{
		if (lineSegmentIntersect == intersect)
		{
			//����߶ε�b���������ߵ�b��˵���������߶�һ��
			intersectPoint.col = lineSegment.col1;//��������꼴Ϊ�߶κ�����
			intersectPoint.row = lineSegment.row1;//���������꼴Ϊ�߶�������
			return true;
		}
		else 
		{
			return false;
		}
	}

	//����߶��������߲�ƽ��
	//����x���� x = (b2 - b1) / (k1 - k2)
	double intersectX = (intersect - lineSegmentIntersect) / (lineSegmentSlope - slope);
	//����y���� y = k1 * x + b
	double intersectY = lineSegmentSlope * intersectX + lineSegmentIntersect;
	//��齻���Ƿ����߶���
	if (vertical)
	{
		//��ֱ����������꣬�ж��Ƿ��ཻ
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
		//ˮƽ���������꣬�ж��Ƿ��ཻ
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
	���ı���������߶ν��зָ���ؽ��㼯
	ͨ������does_segment_intersect_lineʵ��
*/
vector<CoordinateDouble> GlobalWarpping::intersections_with_quad(LineD lineSegment, CoordinateDouble topLeft, CoordinateDouble topRight, CoordinateDouble bottomLeft, CoordinateDouble bottomRight)
{
	vector<CoordinateDouble> intersections;//���㼯

	//��
	double leftSlope = INF;
	if (topLeft.col != bottomLeft.col) 
	{
		leftSlope = (topLeft.row - bottomLeft.row) / (topLeft.col - bottomLeft.col);//��߽�k
	}
	double leftIntersect = topLeft.row - leftSlope * topLeft.col;//��߽�b
	//�ж��߶��Ƿ�͸ñ߽��ཻ
	CoordinateDouble leftIntersectPoint;//����ཻ�ĵ�
	if (does_segment_intersect_line(lineSegment, leftSlope, leftIntersect, true, leftIntersectPoint)) 
	{
		if (leftIntersectPoint.row >= topLeft.row && leftIntersectPoint.row <= bottomLeft.row) 
		{
			intersections.push_back(leftIntersectPoint);
		}
	}

	//��
	double rightSlope = INF;
	if (topRight.col != bottomRight.col) 
	{
		rightSlope = (topRight.row - bottomRight.row) / (topRight.col - bottomRight.col);//�ұ߽�k
	}
	double rightIntersect = topRight.row - rightSlope * topRight.col;//�ұ߽�b
	//�ж��߶��Ƿ�͸ñ߽��ཻ
	CoordinateDouble rightIntersectPoint;
	if (does_segment_intersect_line(lineSegment, rightSlope, rightIntersect, true, rightIntersectPoint)) 
	{
		if (rightIntersectPoint.row >= topRight.row && rightIntersectPoint.row <= bottomRight.row) 
		{
			intersections.push_back(rightIntersectPoint);
		}
	}

	//��
	double topSlope = INF;
	if (topLeft.col != topRight.col) 
	{
		topSlope = (topRight.row - topLeft.row) / (topRight.col - topLeft.col);//�ϱ߽�k
	}
	double topIntersect = topLeft.row - topSlope * topLeft.col;//�ϱ߽�b
	//�ж��߶��Ƿ�͸ñ߽��ཻ
	CoordinateDouble topIntersectPoint;
	if (does_segment_intersect_line(lineSegment, topSlope, topIntersect, false, topIntersectPoint)) 
	{
		if (topIntersectPoint.col >= topLeft.col && topIntersectPoint.col <= topRight.col) 
		{
			intersections.push_back(topIntersectPoint);
		}
	}

	//��
	double bottomSlope = INF;
	if (bottomLeft.col != bottomRight.col) 
	{
		bottomSlope = (bottomRight.row - bottomLeft.row) / (bottomRight.col - bottomLeft.col);//�±߽�k
	}
	double bottomIntersect = bottomLeft.row - bottomSlope * bottomLeft.col;//�±߽�b
	//�ж��߶��Ƿ�͸ñ߽��ཻ
	CoordinateDouble bottomIntersectPoint;
	if (does_segment_intersect_line(lineSegment, bottomSlope, bottomIntersect, false, bottomIntersectPoint)) 
	{
		if (bottomIntersectPoint.col >= bottomLeft.col && bottomIntersectPoint.col <= bottomRight.col)
		{
			intersections.push_back(bottomIntersectPoint);
		}
	}
	return intersections;//���ؽ��㼯
}

/*
	�����������ж�����ֱ�߽��зָ�
*/
vector<vector<vector<LineD>>> GlobalWarpping::segment_line_in_quad(vector<LineD> lines, vector<vector<CoordinateDouble>> mesh)
{
	int QuadnumRow = config.meshQuadRow;
	int QuadnumCol = config.meshQuadCol;
	vector<vector<vector<LineD>>> quad_line_seg;//��һ����άvector����ָ���,�ֱ�Ϊ��,��,ÿ���������߶μ�
	CVMat src2;
	src.copyTo(src2);

	for (int row = 0; row < QuadnumRow; row++) 
	{
		vector<vector<LineD>> vec_row;
		for (int col = 0; col < QuadnumCol; col++) 
		{
			//ÿһ������[row][col]
			//��ȡ������ĸ�����
			CoordinateDouble lefttop = mesh[double(row)][double(col)];
			CoordinateDouble righttop = mesh[double(row)][double(col) + 1];
			CoordinateDouble leftbottom = mesh[double(row) + 1][double(col)];
			CoordinateDouble rightbottom = mesh[double(row) + 1][double(col) + 1];

			vector<LineD> lineInQuad;//�洢����������߶ηָ�Ľ��
			//����ÿһ��ֱ��
			for (int i = 0; i < lines.size(); i++) 
			{
				LineD line = lines[i];
				CoordinateDouble point1(line.row1, line.col1);//�߶ζ˵�1
				CoordinateDouble point2(line.row2, line.col2);//�߶ζ˵�2
				bool p1InQuad = is_in_quad(point1, lefttop, righttop, leftbottom, rightbottom);
				bool p2InQuad = is_in_quad(point2, lefttop, righttop, leftbottom, rightbottom);
				if (p1InQuad && p2InQuad) 
				{
					//���p1,p2�����ı�����,���߶α����ı�����
					lineInQuad.push_back(line);
				}
				else if (p1InQuad) 
				{
					//ֻ��p1���ı�����
					vector<CoordinateDouble> intersections = intersections_with_quad(line, lefttop, righttop, leftbottom, rightbottom);
					//����н���,���Ӧ�߶������˵�ֱ�Ϊpoint1, intersections[0]
					if (intersections.size() != 0) 
					{
						LineD cutLine(point1, intersections[0]);
						lineInQuad.push_back(cutLine);
					}
				}
				else if (p2InQuad) 
				{
					//ֻ��p2���ı�����
					vector<CoordinateDouble> intersections = intersections_with_quad(line, lefttop, righttop, leftbottom, rightbottom);
					if (intersections.size() != 0) 
					{
						LineD cutLine(point2, intersections[0]);
						lineInQuad.push_back(cutLine);
					}
				}
				else 
				{
					//�����˵㶼�����ı�����
					vector<CoordinateDouble> intersections = intersections_with_quad(line, lefttop, righttop, leftbottom, rightbottom);
					//�������������,���и��߶������˵�ֱ�Ϊintersections[0], intersections[1]
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
	����ά���߶ηָ���չ����һά��vector<LineD>
*/
vector<LineD> GlobalWarpping::flatten(vector<vector<vector<LineD>>>& lineSeg)
{
	vector<LineD> line_vec;//���vector
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
	��origin�����Ͻ�����չ,������ĺϲ�,������չ֮��ľ���
*/
SpareseMatrixD_Row GlobalWarpping::block_diag_extend(SpareseMatrixD_Row origin, MatrixXd addin, int QuadID)
{
	int cols_total = 8 * config.meshQuadRow * config.meshQuadCol;//�ܵ�����,ÿһ������ԪΪ8*8(Shape Preservation Mat)
	SpareseMatrixD_Row res(origin.rows() + addin.rows(), cols_total);//������չ,�����ϰ벿��Ϊorigin,�°벿��Ϊaddin
	res.topRows(origin.rows()) = origin;//�����ϰ벿��

	int lefttop_row = origin.rows();//addin��ʼ������
	int lefttop_col = 8 * QuadID;//addin��ʼ������
	for (int row = 0; row < addin.rows(); row++) 
	{
		for (int col = 0; col < addin.cols(); col++) 
		{
			res.insert(lefttop_row + row, lefttop_col + col) = addin(row, col);//��addin������뵽res����
		}
	}
	res.makeCompressed();//��ϡ��������ѹ��
	return res;
}

/*
	��ʼ���߶ηָ�,�����������б�Ƕȵ��߶η��䵽һ��������
*/
vector<vector<vector<LineD>>> GlobalWarpping::init_line_seg(CVMat mask, vector<LineD>& lineSeg_flatten, vector<vector<CoordinateDouble>> mesh, vector<pair<int,double>>& id_theta, vector<double>& rotate_theta)
{
	double thetaPerbin = PI / 49;
	revise_mask_for_lines(mask);
	//��һ��:ʹ��LSD����ֱ�߼��
	vector<LineD> lines = lsd_detect(mask);
	//�ڶ���:������ָ������߶ε����Ե������ڲ�
	vector<vector<vector<LineD>>> lineSeg = segment_line_in_quad(lines, mesh);
	//������:����ÿ���߶ε���б�Ƕ�,���߷��䵽Ԥ�����õĸ����ڲ�
	lineSeg_flatten = flatten(lineSeg);//���õ�����άvectorչ��
	for (int i = 0; i < lineSeg_flatten.size(); i++) 
	{
		LineD line = lineSeg_flatten[i];//ÿ��ֱ��
		double theta = atan((line.row1 - line.row2) / (line.col1 - line.col2));//�����߶ε���б�Ƕ�arctan(y1-y1)/(x1-x2)
		int lineSegmentBucket = (int)round((theta + PI / 2) / thetaPerbin);//�����߶ζ�Ӧ��bin id
		assert(lineSegmentBucket < 50);//��֤lineSegmentBucket < 50
		id_theta.push_back(make_pair(lineSegmentBucket, theta));//�����߶μ��뵽id_theta��
		rotate_theta.push_back(0);//
	}
	return lineSeg;
}

/*
	��ȡֱ����������
*/
SpareseMatrixD_Row GlobalWarpping::get_line_mat(CVMat mask, vector<vector<CoordinateDouble>> mesh, vector<double>rotate_theta, vector<vector<vector<LineD>>> lineSeg, vector<pair<MatrixXd, MatrixXd>>& BilinearVec, int& linenum, vector<bool>& bad) 
{
	int linetmpnum = -1;
	int rows = config.rows;
	int cols = config.cols;
	int QuadnumRow = config.meshQuadRow;
	int QuadnumCol = config.meshQuadCol;
	double gridcols = config.colPermesh;//ÿ��������ռ������
	double gridrows = config.rowPermesh;//ÿ��������ռ������

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
	SpareseMatrixD_Row energy_line;//line��������
	for (int row = 0; row < QuadnumRow; row++) 
	{
		for (int col = 0; col < QuadnumCol; col++) 
		{
			//����ÿһ������
			vector<LineD> linesegInquad = lineSeg[row][col];//����mesh[row][col]�зָ�õ����߶μ���
			int QuadID = row * QuadnumCol + col;
			if (linesegInquad.size() == 0) 
			{
				//�����ǰ����û���߶�,��ֱ����һ������
				continue;
			}
			
			Coordinate topleft(row, col);//���ϽǶ���
			MatrixXd C_row_stack(0, 8);//��������ÿһ�������������߶ζ�Ӧ��C����

			for (int k = 0; k < linesegInquad.size(); k++) 
			{
				//����ÿһ���߶�
				linetmpnum++;//��¼��ǰ�ߵ�����,������������
				LineD line = linesegInquad[k];//��k����
				CoordinateDouble linestart(line.row1, line.col1);//�߶����
				CoordinateDouble lineend(line.row2, line.col2);//�߶��յ�

				BilinearWeights startWeight = get_bilinear_weights(linestart, topleft, mesh);//s t2n t t1n
				MatrixXd start_W_mat = BilinearWeightsToMatrix(startWeight);
				BilinearWeights endWeight = get_bilinear_weights(lineend, topleft, mesh);
				MatrixXd end_W_mat = BilinearWeightsToMatrix(endWeight);
				//cout << startWeight.s << " " << startWeight.t << endl;//test
				//test
				VectorXd S = get_vertices(row, col, mesh);
				Vector2d ans = start_W_mat * S - Vector2d(linestart.col, linestart.row);
				Vector2d ans2 = end_W_mat * S - Vector2d(lineend.col, lineend.row);

				//���ans2��2����>=0.0001����ans��2����>=0.0001
				if (ans2.norm() >= 0.0001 || ans.norm() >= 0.0001) 
				{
					//�������
					bad.push_back(true);
					BilinearVec.push_back(make_pair(MatrixXd::Zero(2, 8), MatrixXd::Zero(2, 8)));
					continue;
				}
				bad.push_back(false);
				//assert(ans.norm() < 0.0001);
				//assert(ans2.norm() < 0.0001);
				
				//end test
				//system("pause");

				double theta = rotate_theta[linetmpnum];//�����е�theta_m
				BilinearVec.push_back(make_pair(start_W_mat, end_W_mat));
				Matrix2d R;//�����е���ת����R
				R << cos(theta), -sin(theta),
					sin(theta), cos(theta);
				MatrixXd ehat(2, 1);//ehat
				ehat << line.col1 - line.col2, line.row1 - line.row2;//����ehat
				MatrixXd tmp = (ehat.transpose() * ehat).inverse();//(ehat^T * ehat)^-1
				Matrix2d I = Matrix2d::Identity();//��λ����
				MatrixXd C = R * ehat * tmp * (ehat.transpose()) * (R.transpose()) - I;//����(6)�е�C
				//�߶η�������e = start_W_mat - end_W_mat,
				MatrixXd CT = C * (start_W_mat - end_W_mat); //C*e
				C_row_stack = row_stack(C_row_stack, CT);//��CT�ӵ�C_row_stack�·�,���кϲ�
			}
			energy_line = block_diag_extend(energy_line, C_row_stack, QuadID);//��C_row_stack�ӵ�energy_line�����·�,������չ
		}
	}
	linenum = linetmpnum;
	return energy_line;
}

