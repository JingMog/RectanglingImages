#include"Common.h"

/*
	�Hole,����ֵΪCVMat dstBw
*/
CVMat fillHole(const CVMat srcBw)
{
	cv::Size m_Size = srcBw.size();
	CVMat Temp = CVMat::zeros(m_Size.height + 2, m_Size.width + 2, srcBw.type());//��չͼ��,��ȼ�2���߶ȼ�2
	srcBw.copyTo(Temp(cv::Range::Range(1, m_Size.height + 1), cv::Range::Range(1, m_Size.width + 1)));

	cv::floodFill(Temp, cv::Point(0, 0), cv::Scalar(255)); //ʹ�÷������
	CVMat cutImg;//�ü���չ��ͼ��
	Temp(cv::Range::Range(1, m_Size.height + 1), cv::Range::Range(1, m_Size.width + 1)).copyTo(cutImg);

	return srcBw | (~cutImg);
}

/*
	�����ͨ��
*/
CVMat Mask_contour(CVMat src)
{
	CVMat bw;
	cv::cvtColor(src, src, CV_BGR2GRAY);//����ɫͼת��Ϊ�Ҷ�ͼ

	uchar thr = 252;
	CVMat mask = CVMat::zeros(src.size(), CV_8UC1);
	//����mask����ԭͼ�лҶ�ֵС��252�ľ���Ϊ255
	for (int row = 0; row < src.rows; row++)
	{
		for (int col = 0; col < src.cols; col++)
		{
			if (src.at<uchar>(row, col) < thr)
				mask.at<uchar>(row, col) = 255;
		}
	}

	bw = fillHole(mask);//����mask���
	bw = ~bw;//��ת���ں����ĸ�ʴ�����Ͳ���

	CVMat element = cv::getStructuringElement(cv::MORPH_ELLIPSE, cv::Size(5, 5));
	CVMat dilate_out;//���ͽ��
	cv::dilate(bw, dilate_out, element);//����
	cv::dilate(dilate_out, dilate_out, element);//��������
	cv::dilate(dilate_out, dilate_out, element);//��������

	CVMat erode_out;//��ʴ���
	erode(dilate_out, erode_out, element);//��ʴ
	return erode_out;
}

/*
	����ά����ת��Ϊһά����
	�����Ĵ洢��ʽΪ����Ϊ��,col,row,...
*/
VectorXd mesh_to_vector(vector<vector<CoordinateDouble>> mesh, Config config)
{
	int numMeshRow = config.meshNumRow;//���������,������Ϊ20
	int numMeshCol = config.meshNumCol;//���������,������Ϊ20
	VectorXd vec = VectorXd::Zero(numMeshRow * numMeshCol * 2);
	for (int row = 0; row < numMeshRow; row++)
	{
		for (int col = 0; col < numMeshCol; col++)
		{
			CoordinateDouble coord = mesh[row][col];//��������
			vec((row * numMeshCol + col) * 2) = coord.col;//����������
			vec((row * numMeshCol + col) * 2 + 1) = coord.row;//����������
		}
	}
	return vec;
}

/*
	��һά��������ת��Ϊ��ά��������,Ϊmesh_to_vector�����������
*/
vector<vector<CoordinateDouble>> vector_to_mesh(VectorXd x, Config config)
{
	int numMeshRow = config.meshNumRow;//���������
	int numMeshCol = config.meshNumCol;//���������
	vector<vector<CoordinateDouble>> mesh;//��άvector
	for (int row = 0; row < numMeshRow; row++)
	{
		vector<CoordinateDouble> meshRow;
		for (int col = 0; col < numMeshCol; col++)
		{
			int xid = (row * numMeshCol + col) * 2;
			CoordinateDouble coord;
			coord.row = x(xid + 1);//������
			coord.col = x(xid);//������
			meshRow.push_back(coord);//���뵽meshRow��
		}
		mesh.push_back(meshRow);//��meshRow���뵽mesh��
	}
	return mesh;
}

/*
	���SparseMatrixD������Ϣ
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
	�����кϲ�,�ϰ벿��Ϊorigin,�°벿��Ϊdiag
*/
SpareseMatrixD_Row row_stack(SparseMatrixD origin, SpareseMatrixD_Row diag)
{
	SpareseMatrixD_Row res(origin.rows() + diag.rows(), origin.cols());
	res.topRows(origin.rows()) = origin;//������ǰorigin.rows()�и�ֵΪorigin
	res.bottomRows(diag.rows()) = diag;//�������didg.rows()�и���Ϊdiag
	return res;
}

/*
	�����кϲ�,�ϰ벿��Ϊorigin,�°벿��Ϊdiag
*/
SpareseMatrixD_Row row_stack(SpareseMatrixD_Row origin, SpareseMatrixD_Row diag)
{
	SpareseMatrixD_Row res(origin.rows() + diag.rows(), origin.cols());
	res.topRows(origin.rows()) = origin;
	res.bottomRows(diag.rows()) = diag;
	return res;
}

/*
	�����кϲ�,�ϰ벿��Ϊorigin,�°벿��Ϊdiag
*/
MatrixXd row_stack(MatrixXd mat1, MatrixXd mat2)
{
	MatrixXd res(mat1.rows() + mat2.rows(), mat1.cols());
	res.topRows(mat1.rows()) = mat1;
	res.bottomRows(mat2.rows()) = mat2;
	return res;
}

/*
	�����кϲ�,��벿��Ϊmat1,�Ұ벿��Ϊmat2
*/
MatrixXd col_stack(MatrixXd mat1, MatrixXd mat2)
{
	MatrixXd res(mat1.rows(), mat1.cols() + mat2.cols());
	res.leftCols(mat1.cols()) = mat1;
	res.rightCols(mat2.cols()) = mat2;
	return res;
}

/*
	����startPoint,endPoint����ֱ��
*/
void DrawLine(CVMat& img, CoordinateDouble startPoint, CoordinateDouble endPoint)
{
	cv::Point start((int)startPoint.col, (int)startPoint.row);
	cv::Point end((int)endPoint.col, (int)endPoint.row);
	int thickness = 1;//�����Ĵ�ϸ
	int lineType = 1;//����������
	cv::line(img, start, end, cv::Scalar(0, 255, 0), thickness, lineType);
}

/*
	����LineD����ֱ��
*/
void DrawLine(CVMat& img, LineD line)
{
	cv::Point start((int)line.col1, (int)line.row1);//���
	cv::Point end((int)line.col2, (int)line.row2);//�յ�
	int thickness = 1;
	int lineType = 1;
	cv::line(img, start, end, cv::Scalar(0, 255, 0), thickness, lineType);
}

/*
	��Դͼ��src�л�������mesh
*/
CVMat drawmesh(CVMat src, vector<vector<CoordinateDouble>> mesh, Config config)
{
	int meshNumRow = config.meshNumRow;
	int meshNumCol = config.meshNumCol;

	// ��������,��ά
	for (int row = 0; row < meshNumRow; row++)
	{
		for (int col = 0; col < meshNumCol; col++)
		{
			CoordinateDouble now = mesh[row][col];//��ǰ�����
			if (row == meshNumRow - 1 && col < meshNumCol - 1)
			{
				// ���һ��,����ǰ�㵽�Ҳ���ֱ��
				CoordinateDouble right = mesh[row][col + 1];
				DrawLine(src, now, right);
			}
			else if (row < meshNumRow - 1 && col == meshNumCol - 1)
			{
				// ���һ��,����ǰ�㵽�·����ֱ��
				CoordinateDouble down = mesh[row + 1][col];
				DrawLine(src, now, down);
			}
			else if (row < meshNumRow - 1 && col < meshNumCol - 1)
			{
				// �����һ��Ҳ�����һ��
				// ��Ҫͬʱ�����Ҳ���ߺ͵��·�����
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
	��������,
	������������enlargeFacrtor����Ӧ����������enlargeFacrtor��
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
