#include"Common.h"

/*
	填补Hole,返回值为CVMat dstBw
*/
CVMat fillHole(const CVMat srcBw)
{
	cv::Size m_Size = srcBw.size();
	CVMat Temp = CVMat::zeros(m_Size.height + 2, m_Size.width + 2, srcBw.type());//延展图像,宽度加2，高度加2
	srcBw.copyTo(Temp(cv::Range::Range(1, m_Size.height + 1), cv::Range::Range(1, m_Size.width + 1)));

	cv::floodFill(Temp, cv::Point(0, 0), cv::Scalar(255)); //使用泛洪填充
	CVMat cutImg;//裁剪延展的图像
	Temp(cv::Range::Range(1, m_Size.height + 1), cv::Range::Range(1, m_Size.width + 1)).copyTo(cutImg);

	return srcBw | (~cutImg);
}

/*
	标记连通域
*/
CVMat Mask_contour(CVMat src)
{
	CVMat bw;
	cv::cvtColor(src, src, CV_BGR2GRAY);//将彩色图转换为灰度图

	uchar thr = 252;
	CVMat mask = CVMat::zeros(src.size(), CV_8UC1);
	//生成mask，将原图中灰度值小于252的均变为255
	for (int row = 0; row < src.rows; row++)
	{
		for (int col = 0; col < src.cols; col++)
		{
			if (src.at<uchar>(row, col) < thr)
				mask.at<uchar>(row, col) = 255;
		}
	}

	bw = fillHole(mask);//根据mask填充
	bw = ~bw;//反转用于后续的腐蚀和膨胀操作

	CVMat element = cv::getStructuringElement(cv::MORPH_ELLIPSE, cv::Size(5, 5));
	CVMat dilate_out;//膨胀结果
	cv::dilate(bw, dilate_out, element);//膨胀
	cv::dilate(dilate_out, dilate_out, element);//二次膨胀
	cv::dilate(dilate_out, dilate_out, element);//三次膨胀

	CVMat erode_out;//腐蚀结果
	erode(dilate_out, erode_out, element);//腐蚀
	return erode_out;
}

/*
	将二维网格转化为一维向量
	向量的存储方式为行序为主,col,row,...
*/
VectorXd mesh_to_vector(vector<vector<CoordinateDouble>> mesh, Config config)
{
	int numMeshRow = config.meshNumRow;//网格的行数,论文中为20
	int numMeshCol = config.meshNumCol;//网格的列数,论文中为20
	VectorXd vec = VectorXd::Zero(numMeshRow * numMeshCol * 2);
	for (int row = 0; row < numMeshRow; row++)
	{
		for (int col = 0; col < numMeshCol; col++)
		{
			CoordinateDouble coord = mesh[row][col];//网格坐标
			vec((row * numMeshCol + col) * 2) = coord.col;//网格点的列数
			vec((row * numMeshCol + col) * 2 + 1) = coord.row;//网格点的行数
		}
	}
	return vec;
}

/*
	将一维向量坐标转换为二维网格坐标,为mesh_to_vector函数的逆过程
*/
vector<vector<CoordinateDouble>> vector_to_mesh(VectorXd x, Config config)
{
	int numMeshRow = config.meshNumRow;//网格的行数
	int numMeshCol = config.meshNumCol;//网格的列数
	vector<vector<CoordinateDouble>> mesh;//二维vector
	for (int row = 0; row < numMeshRow; row++)
	{
		vector<CoordinateDouble> meshRow;
		for (int col = 0; col < numMeshCol; col++)
		{
			int xid = (row * numMeshCol + col) * 2;
			CoordinateDouble coord;
			coord.row = x(xid + 1);//行坐标
			coord.col = x(xid);//列坐标
			meshRow.push_back(coord);//加入到meshRow中
		}
		mesh.push_back(meshRow);//将meshRow加入到mesh中
	}
	return mesh;
}

/*
	输出SparseMatrixD矩阵信息
*/
void print_sparse_mat(SparseMatrixD Q)
{
	for (int k = 0; k < Q.outerSize(); k++)
	{
		for (SparseMatrixD::InnerIterator it(Q, k); it; ++it)
		{
			cout << it.row() << " " << it.col() << " : " << it.value() << endl;
		}
		cout << endl;
	}
}

/*
	矩阵按行合并,上半部分为origin,下半部分为diag
*/
SpareseMatrixD_Row row_stack(SparseMatrixD origin, SpareseMatrixD_Row diag)
{
	SpareseMatrixD_Row res(origin.rows() + diag.rows(), origin.cols());
	res.topRows(origin.rows()) = origin;//将矩阵前origin.rows()行赋值为origin
	res.bottomRows(diag.rows()) = diag;//将矩阵后didg.rows()行复制为diag
	return res;
}

/*
	矩阵按行合并,上半部分为origin,下半部分为diag
*/
SpareseMatrixD_Row row_stack(SpareseMatrixD_Row origin, SpareseMatrixD_Row diag)
{
	SpareseMatrixD_Row res(origin.rows() + diag.rows(), origin.cols());
	res.topRows(origin.rows()) = origin;
	res.bottomRows(diag.rows()) = diag;
	return res;
}

/*
	矩阵按行合并,上半部分为origin,下半部分为diag
*/
MatrixXd row_stack(MatrixXd mat1, MatrixXd mat2)
{
	MatrixXd res(mat1.rows() + mat2.rows(), mat1.cols());
	res.topRows(mat1.rows()) = mat1;
	res.bottomRows(mat2.rows()) = mat2;
	return res;
}

/*
	矩阵按列合并,左半部分为mat1,右半部分为mat2
*/
MatrixXd col_stack(MatrixXd mat1, MatrixXd mat2)
{
	MatrixXd res(mat1.rows(), mat1.cols() + mat2.cols());
	res.leftCols(mat1.cols()) = mat1;
	res.rightCols(mat2.cols()) = mat2;
	return res;
}

/*
	根据startPoint,endPoint绘制直线
*/
void DrawLine(CVMat& img, CoordinateDouble startPoint, CoordinateDouble endPoint)
{
	cv::Point start((int)startPoint.col, (int)startPoint.row);
	cv::Point end((int)endPoint.col, (int)endPoint.row);
	int thickness = 1;//线条的粗细
	int lineType = 1;//线条的类型
	cv::line(img, start, end, cv::Scalar(0, 255, 0), thickness, lineType);
}

/*
	根据LineD绘制直线
*/
void DrawLine(CVMat& img, LineD line)
{
	cv::Point start((int)line.col1, (int)line.row1);//起点
	cv::Point end((int)line.col2, (int)line.row2);//终点
	int thickness = 1;
	int lineType = 1;
	cv::line(img, start, end, cv::Scalar(0, 255, 0), thickness, lineType);
}

/*
	在源图像src中绘制网格mesh
*/
CVMat drawmesh(CVMat src, vector<vector<CoordinateDouble>> mesh, Config config)
{
	int meshNumRow = config.meshNumRow;
	int meshNumCol = config.meshNumCol;

	// 遍历网格,二维
	for (int row = 0; row < meshNumRow; row++)
	{
		for (int col = 0; col < meshNumCol; col++)
		{
			CoordinateDouble now = mesh[row][col];//当前网格点
			if (row == meshNumRow - 1 && col < meshNumCol - 1)
			{
				// 最后一行,画当前点到右侧点的直线
				CoordinateDouble right = mesh[row][col + 1];
				DrawLine(src, now, right);
			}
			else if (row < meshNumRow - 1 && col == meshNumCol - 1)
			{
				// 最后一列,画当前点到下方点的直线
				CoordinateDouble down = mesh[row + 1][col];
				DrawLine(src, now, down);
			}
			else if (row < meshNumRow - 1 && col < meshNumCol - 1)
			{
				// 非最后一行也非最后一列
				// 需要同时画到右侧的线和到下方的线
				CoordinateDouble right = mesh[row][col + 1];
				DrawLine(src, now, right);
				CoordinateDouble down = mesh[row + 1][col];
				DrawLine(src, now, down);
			}
		}
	}
	cv::namedWindow("Mesh", CV_WINDOW_AUTOSIZE);
	cv::imshow("Mesh", src);
	cv::waitKey(0);
	return src;
}

/*
	扩充网格,
	根据扩充因子enlargeFacrtor来相应将坐标扩充enlargeFacrtor倍
*/
void enlarge_mesh(vector<vector<CoordinateDouble>>& mesh, double enlargeFacrtor, Config config)
{
	int numMeshRow = config.meshNumRow;
	int numMeshCol = config.meshNumCol;
	for (int row = 0; row < numMeshRow; row++)
	{
		for (int col = 0; col < numMeshCol; col++)
		{
			CoordinateDouble& coord = mesh[row][col];
			coord.row = coord.row * enlargeFacrtor;//row * enlargeFacrtor
			coord.col = coord.col * enlargeFacrtor;//col * enlargeFacrtor
		}
	}
};
