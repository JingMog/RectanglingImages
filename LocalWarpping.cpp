#include "LocalWarpping.h"

bool cmp(const pair<int, float> a, const pair<int, float> b)
{
	return a.second < b.second;
}

/*
	�ж�mask��(row,col)���Ƿ�Ϊ��(����ɫ)
*/
bool LocalWarpping::Is_transparent(CVMat& mask, int row, int col)
{
	if (mask.at<uchar>(row, col) == 0) 
		return false;
	else 
		return true;
}

/*
	��ʼ��displacement
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
	ʹ��Sobel��ͼ����б�Ե���
*/
CVMat LocalWarpping::Sobel_img(CVMat src)
{
	CVMat gray;
	//��RGBͼת��Ϊ�Ҷ�ͼ���б�Ե���
	cv::cvtColor(src, gray, CV_BGR2GRAY);
	CVMat grad_x, grad_y, dst;
	//Sobel(����ĻҶ�ͼ��,������ݶ�,����ݶȵ���������,
	//dx=1,dy=0,����x�����ϵ��ݶ�
	cv::Sobel(gray, grad_x, CV_8U, 1, 0);
	//dx=0,dy=1,����y�����ϵ��ݶ�
	cv::Sobel(gray, grad_y, CV_8U, 0, 1);
	addWeighted(grad_x, 0.5, grad_y, 0.5, 0, dst);//��gray_x��gray_y��Ȩ���
	return dst;
}

/*
	�ھ��λ�֮���ͼ���Ϸ���һ����׼����,
	Ȼ����Ҫ����wrap back���ص�ԭʼ������ͼ����,�õ�����
*/
void LocalWarpping::warp_mesh_back(vector<vector<CoordinateDouble>>& mesh, vector<vector<Coordinate>> displacementMap, Config config)
{
	int meshnum_row = config.meshNumRow;//���������
	int meshnum_col = config.meshNumCol;//���������

	for (int row_mesh = 0; row_mesh < meshnum_row; row_mesh++)
	{
		for (int col_mesh = 0; col_mesh < meshnum_col; col_mesh++)
		{
			//���λ�ڱ�Ե
			if (row_mesh == meshnum_row - 1 && col_mesh == meshnum_col - 1)
			{
				CoordinateDouble& meshVertexPoint = mesh[row_mesh][col_mesh];//���������
				Coordinate vertexDisplacement = displacementMap[floor(meshVertexPoint.row) - 1][floor(meshVertexPoint.col) - 1];
				meshVertexPoint.row += vertexDisplacement.row;
				meshVertexPoint.col += vertexDisplacement.col;
			}
			//�Ǳ�Ե,ֱ�ӽ�ԭ�������displacement
			CoordinateDouble& meshVertexPoint = mesh[row_mesh][col_mesh];
			Coordinate vertexDisplacement = displacementMap[(int)floor(meshVertexPoint.row)][(int)floor(meshVertexPoint.col)];
			meshVertexPoint.row += vertexDisplacement.row;
			meshVertexPoint.col += vertexDisplacement.col;
		}
	}
}

/*
	��ͼ���Ե�ҵ���Ŀհ�,������ʼ����������
	֮����ݸ�����ѡ����ͼ,������ͼ����������,��Sobel��Ե���
	�ҵ�������̵�һ����,Ȼ��ִ��Seam Carving
*/
pair<int, int> LocalWarpping::Choose_longest_border(CVMat src, CVMat mask, Border& direction)
{
	//�����border��[begin��end]
	int maxLength = 0;//��¼��󳤶�
	int rows = src.rows;//rows
	int cols = src.cols;//cols

	int final_startIndex = 0;//�����ʼ����
	int final_endIndex = 0;//�����������

	//���������ĸ����������border,Ȼ��������seam carving
	//left
	int tmp_maxLength, tmp_startIndex, tmp_endIndex;
	tmp_maxLength = tmp_startIndex = tmp_endIndex = 0;
	bool isCounting = false;
	for (int row = 0; row < rows; row++)
	{
		//���mask��(row, 0)�㲻Ϊ��(�հ��ж�),����Ϊ���һ��,�������ֹ
		if (!Is_transparent(mask, row, 0) || row == rows - 1)
		{
			if (isCounting)
			{	//������ڼ���������Ϊֹ
				if (Is_transparent(mask, row, 0))
				{
					tmp_endIndex++;
					tmp_maxLength++;
					isCounting = true;
				}
				if (tmp_maxLength > maxLength)
				{
					//�������ֵ
					maxLength = tmp_maxLength;
					final_startIndex = tmp_startIndex;
					final_endIndex = tmp_endIndex;
					direction = Border::BORDER_LEFT;
				}
			}
			isCounting = false;//��ֹ����
			tmp_startIndex = tmp_endIndex = row;
			tmp_maxLength = 0;
		}
		else
		{
			//�õ�͸��,��������
			tmp_endIndex++;
			tmp_maxLength++;
			isCounting = true;
		}
	}

	//�Ҳ�
	tmp_maxLength = tmp_startIndex = tmp_endIndex = 0;
	isCounting = false;
	for (int row = 0; row < rows; row++) 
	{
		//���mask��(row, cols - 1)�㲻Ϊ��(�հ��ж�),����Ϊ���һ��,�������ֹ
		if (!Is_transparent(mask, row, cols - 1) || row == rows - 1) 
		{
			if (isCounting) 
			{
				//������ڼ���������Ϊֹ
				if (Is_transparent(mask, row, cols - 1)) 
				{
					tmp_endIndex++;
					tmp_maxLength++;
					isCounting = true;
				}
				if (tmp_maxLength > maxLength) 
				{
					//�������ֵ
					maxLength = tmp_maxLength;
					final_startIndex = tmp_startIndex;
					final_endIndex = tmp_endIndex;
					direction = Border::BORDER_RIGHT;
				}
			}
			isCounting = false;//��ֹ����
			tmp_startIndex = tmp_endIndex = row;
			tmp_maxLength = 0;
		}
		else 
		{
			//�õ�͸������ʼ����
			tmp_endIndex++;
			tmp_maxLength++;
			isCounting = true;
		}
	}

	//��
	tmp_maxLength = tmp_startIndex = tmp_endIndex = 0;
	isCounting = false;
	for (int col = 0; col < cols; col++) 
	{
		//���mask��(0, col)�㲻Ϊ��(�հ��ж�),����Ϊ���һ��,�������ֹ
		if (!Is_transparent(mask, 0, col) || col == cols - 1) 
		{
			if (isCounting) 
			{	//������ڼ���������Ϊֹ
				if (Is_transparent(mask, 0, col))
				{
					tmp_endIndex++;
					tmp_maxLength++;
					isCounting = true;
				}
				if (tmp_maxLength > maxLength)
				{
					//�������ֵ
					maxLength = tmp_maxLength;
					final_startIndex = tmp_startIndex;
					final_endIndex = tmp_endIndex;
					direction = Border::BORDER_TOP;
				}
			}
			isCounting = false;//��ֹ����
			tmp_startIndex = tmp_endIndex = col;
			tmp_maxLength = 0;
		}
		else 
		{	//�õ�͸������ʼ����
			tmp_endIndex++;
			tmp_maxLength++;
			isCounting = true;
		}
	}

	//�·�
	tmp_maxLength = tmp_startIndex = tmp_endIndex = 0;
	isCounting = false;
	for (int col = 0; col < cols; col++) 
	{
		//���mask��(row - 1, col)�㲻Ϊ��(�հ��ж�),����Ϊ���һ��,�������ֹ
		if (!Is_transparent(mask, rows - 1, col) || col == cols - 1)
		{
			if (isCounting)
				{//������ڼ���������Ϊֹ
				if (Is_transparent(mask, rows - 1, col)) 
				{
					tmp_endIndex++;
					tmp_maxLength++;
					isCounting = true;
				}
				if (tmp_maxLength > maxLength) 
				{
					//�������ֵ
					maxLength = tmp_maxLength;
					final_startIndex = tmp_startIndex;
					final_endIndex = tmp_endIndex;
					direction = Border::BORDER_BOTTOM;
				}
			}
			isCounting = false;//��ֹ����
			tmp_startIndex = tmp_endIndex = col;
			tmp_maxLength = 0;
		}
		else 
		{	//�õ�͸������ʼ����
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
	��ʾ��ѡ����border
	����choose_longest_borderѡ�����border��Ȼ�󽫸�Border����Ϊ��ɫ
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
	���ݻ�ȡ��seam���������ص�ƽ��
	����һ��seam(Seam Carving�㷨)
*/
CVMat LocalWarpping::Insert_local_seam(CVMat src, CVMat& mask, int* seam, SeamDirection seamdirection, pair<int, int> begin_end, bool shiftToend)
{
	//���seam��ֱ,��src��maskת��,ת��Ϊˮƽ
	if (seamdirection == SeamDirection::SEAM_HORIZONTAL) 
	{
		transpose(src, src);
		transpose(mask, mask);
	}

	CVMat resimg;
	src.copyTo(resimg);

	int begin = begin_end.first;//�border���ڵ�local row��Χ
	int end = begin_end.second;
	int rows = src.rows;
	int cols = src.cols;

	for (int row = begin; row <= end; row++)
	{
		//����ÿһ�У�seam�Ҳࣨ��ࣩ���������ؾ�����Եƽ��
		int local_row = row - begin;//��border��ͼ�е�λ��
		if (!shiftToend)
		{
			//seam����������ص㣬�����ƶ�һ��
			for (int col = 0; col < seam[local_row]; col++)
			{
				colorPixel pixel = src.at<colorPixel>(row, col + 1);
				resimg.at<colorPixel>(row, col) = pixel;
				mask.at<uchar>(row, col) = mask.at<uchar>(row, col + 1);
			}
		}
		else
		{
			//seam�Ҳ��������ص㣬�����ƶ�һ��
			for (int col = cols - 1; col > seam[local_row]; col--)
			{
				colorPixel pixel = src.at<colorPixel>(row, col - 1);
				resimg.at<colorPixel>(row, col) = pixel;
				mask.at<uchar>(row, col) = mask.at<uchar>(row, col - 1);
			}
		}

		//�Է�϶���߽���ƽ���ɣ���������ƽ��ֵ��
		mask.at<uchar>(row, seam[local_row]) = 0;
		if (seam[local_row] == 0) 
		{
			//�����϶λ�ڵ�һ��,���Ҳ������������һ��
			resimg.at<colorPixel>(row, seam[local_row]) = src.at<colorPixel>(row, seam[local_row] + 1);
		}
		else if (seam[local_row] == cols - 1)
		{
			//�����϶λ�����һ�У������������������һ��
			resimg.at<colorPixel>(row, seam[local_row]) = src.at<colorPixel>(row, seam[local_row] - 1);
		}
		else 
		{
			//����seam�������ҽ���ƽ������
			colorPixel pixel1 = src.at<colorPixel>(row, seam[local_row] + 1);
			colorPixel pixel2 = src.at<colorPixel>(row, seam[local_row] - 1);
			resimg.at<colorPixel>(row, seam[local_row]) = pixel1 / 2 + pixel2 / 2;
		}
	}

	//����ת��
	if (seamdirection == SeamDirection::SEAM_HORIZONTAL)
	{
		cv::transpose(resimg, resimg);
		cv::transpose(mask, mask);
	}

	return resimg;
}

/*
	��Borderȷ������ͼ��ʹ��DP��ȡ������С��seam
	ʹ�ö�̬�滮�㷨����ϯ������ȡͼ�������������Ȼ�����·�����ݣ�����·��������int* seam��
*/
int* LocalWarpping::Get_local_seam(CVMat src, CVMat mask, SeamDirection seamdirection, pair<int, int> begin_end)
{
	//���Ϊˮƽ,�����ת�ý���ת��Ϊ��ֱ
	if (seamdirection == SeamDirection::SEAM_HORIZONTAL)
	{
		cv::transpose(src, src);
		cv::transpose(mask, mask);
	}

	//ͳһѰ����ֱ��seam
	int rows = src.rows;
	int cols = src.cols;

	//�У��еķ�Χ
	int row_start = begin_end.first;
	int row_end = begin_end.second;
	int col_start = 0;
	int col_end = cols - 1;

	int range = row_end - row_start + 1;//�հ�Border�ĳ���Ϊend - start + 1

	int outputWidth = cols;
	int outputHeight = range;

	CVMat displayimg;
	src.copyTo(displayimg);

	//��ȡ��ͼ,����ͼ����DP���seam
	CVMat local_img = displayimg(cv::Range::Range(row_start, row_end + 1), cv::Range::Range(col_start, col_end + 1));//��ȡ��ͼ
	CVMat local_mask = mask(cv::Range::Range(row_start, row_end + 1), cv::Range::Range(col_start, col_end + 1));//��ȡmask��ͼ��Range����ҿ��������ұ߽���Ҫ+1

	CVMat local_energy = Sobel_img(local_img);//Sobel��Ե����ȡͼ�������
	CVMat local_energy_32f;
	local_energy.convertTo(local_energy_32f, CV_32F);

	//��ͼ���Ե(�Ǿ��β���)����������ΪINF,��ֹseam������һ����
	//�����е�3.1Mesh-free Local Warping
	for (int row = 0; row < range; row++)
		for (int col = col_start; col <= col_end; col++)
			if ((int)local_mask.at<uchar>(row, col) == 255)
				local_energy_32f.at<float>(row, col) = INF;


	//��̬�滮�����С������seam�������·�����Ż���
	CVMat tmpenergy;
	local_energy_32f.copyTo(tmpenergy);
	//1.����������ͼ�е����ص㣬�������ͼ
	for (int row = 1; row < range; row++)
	{
		for (int col = col_start; col <= col_end; col++)
		{
			//����ǵ�һ��,�����Ϸ�,�����Ϸ�
			if (col == col_start)
				tmpenergy.at<float>(row, col) += min(tmpenergy.at<float>(row - 1, col), tmpenergy.at<float>(row - 1, col + 1));
			//��������һ��,�����Ϸ�,���Ϸ�
			else if (col == col_end)
				tmpenergy.at<float>(row, col) += min(tmpenergy.at<float>(row - 1, col - 1), tmpenergy.at<float>(row - 1, col));
			//����point,�������Ϸ�,�Ϸ������Ϸ�
			else
				tmpenergy.at<float>(row, col) += min(tmpenergy.at<float>(row - 1, col), min(tmpenergy.at<float>(row - 1, col - 1), tmpenergy.at<float>(row - 1, col + 1)));
		}
	}

	vector<pair<int, float>> last_row;//�洢���һ�е������Ͷ�Ӧ������ֵ
	for (int col = col_start; col <= col_end; col++)
	{
		last_row.push_back(make_pair(col, tmpenergy.at<float>(range - 1, col)));
	}

	sort(last_row.begin(), last_row.end(), cmp);//��������������
	int* seam = new int[range];
	seam[range - 1] = last_row[0].first;

	//2.·�����ݣ����ݵõ�������ͼ��ȡ����seam
	for (int row = range - 2; row >= 0; row--)
	{
		if (seam[row + 1] == col_start)
		{
			//���seam��һ��Ϊ��ͼ��һ��
			//�����ϣ�����
			if (tmpenergy.at<float>(row, seam[row + 1] + 1) < tmpenergy.at<float>(row, seam[row + 1]))
				seam[row] = seam[row + 1] + 1;
			else
				seam[row] = seam[row + 1];
		}
		else if (seam[row + 1] == col_end)
		{
			//seam[row + 1]Ϊ��ͼ�����һ��
			//�������ϣ���
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
			//(row, seam[row+1]-1), (row, seam[row+1]), (row, seam[row+1]+1)�����е���Сֵ
			//�������ϣ��ϣ���������
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
	//��ʾ�ҵ�������Seam
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
	�ھ��λ�֮���ͼ���Ϸ��ñ�׼��������(20*20)
*/
vector<vector<CoordinateDouble>> LocalWarpping::get_rectangle_mesh(CVMat src, Config config) 
{
	int rows = config.rows;
	int cols = config.cols;
	int meshnum_row = config.meshNumRow;//���������,������Ϊ20
	int meshnum_col = config.meshNumCol;//���������,������Ϊ20
	double row_per_mesh = config.rowPermesh;//ÿһ��������ռ������
	double col_per_mesh = config.colPermesh;//ÿһ��������ռ������

	vector<vector<CoordinateDouble>> mesh;
	for (int row_mesh = 0; row_mesh < meshnum_row; row_mesh++) 
	{
		vector<CoordinateDouble> meshrow;
		for (int col_mesh = 0; col_mesh < meshnum_col; col_mesh++) 
		{
			CoordinateDouble coord;
			coord.row = row_mesh * row_per_mesh;//���񽻵�ĺ�����
			coord.col = col_mesh * col_per_mesh;//���񽻵��������
			meshrow.push_back(coord);
		}
		mesh.push_back(meshrow);
	}
	return mesh;
}


/*
	��ȡͼ���ת������,��ÿ�����ص��ƶ�����
*/
vector<vector<Coordinate>> LocalWarpping::Get_Local_warp_displacement(CVMat src, CVMat mask)
{
	int rows = src.rows;
	int cols = src.cols;

	vector<vector<Coordinate>> displacementMap;//ת������
	vector<vector<Coordinate>> finaldisplacementMap;
	init_displacement(finaldisplacementMap, rows, cols);//��ʼ��displacement
	init_displacement(displacementMap, rows, cols);//��ʼ��displacement

	while (true)
	{
		Border direction;
		pair<int, int> begin_end = Choose_longest_border(src, mask, direction);//ѡ��mask����Ŀհױ�,���򱣴���direction��
		//Show_longest_border(src, begin_end, direction);
		
		if (begin_end.first == begin_end.second)
		{
			//��������ȫ���ƶ���ɣ�ֱ�ӷ���displacement
			/*cv::namedWindow("Border", CV_WINDOW_AUTOSIZE);
			cv::imshow("Border", src);
			cv::waitKey(0);*/
			return displacementMap;
		}
		else
		{
			//���пհױ߽�δƽ�����
			bool shift_to_end = false;
			SeamDirection seamdirection;//��϶�ķ���ˮƽ or ��ֱ��
			switch (direction)
			{
			case Border::BORDER_LEFT:
				seamdirection = SeamDirection::SEAM_VERTICAL;//borderλ����,��seamΪ��ֱ
				shift_to_end = false;
				break;
			case Border::BORDER_RIGHT:
				seamdirection = SeamDirection::SEAM_VERTICAL;//borderλ����,��seamΪ��ֱ
				shift_to_end = true;
				break;
			case Border::BORDER_TOP:
				seamdirection = SeamDirection::SEAM_HORIZONTAL;//borderλ����,��seamΪˮƽ
				shift_to_end = false;
				break;
			case Border::BORDER_BOTTOM:
				seamdirection = SeamDirection::SEAM_HORIZONTAL;//borderλ����,��seamΪˮƽ
				shift_to_end = true;
				break;
			default:
				break;
			}

			int* seam = Get_local_seam(src, mask, seamdirection, begin_end);//ʹ��DP��ȡ������Сseam

			src = Insert_local_seam(src, mask, seam, seamdirection, begin_end, shift_to_end);//����һ��Seam����ƽ�����أ������հ�

			//�����û�����
			for (int row = 0; row < rows; row++) 
			{
				for (int col = 0; col < cols; col++) 
				{
					Coordinate tmpdisplacement;
					if (seamdirection == SeamDirection::SEAM_VERTICAL && row >= begin_end.first && row <= begin_end.second)
					{
						//���Ϊ��ֱ��,����row�ڷ�϶��
						int local_row = row - begin_end.first;//��ͼ�е�row
						if (col > seam[local_row] && shift_to_end) 
						{
							tmpdisplacement.col = -1;//�����ƶ�
						}
						else if (col < seam[local_row] && !shift_to_end) 
						{
							tmpdisplacement.col = 1;//�����ƶ�
						}
					}
					else if (seamdirection == SeamDirection::SEAM_HORIZONTAL && col >= begin_end.first && col <= begin_end.second) 
					{
						//���Ϊˮƽ��,����row�ڷ�϶��
						int local_col = col - begin_end.first;//��ͼ�е�col
						if (row > seam[local_col] && shift_to_end) 
						{
							tmpdisplacement.row = -1;//�����ƶ�
						}
						else if (row < seam[local_col] && !shift_to_end) 
						{
							tmpdisplacement.row = 1;//�����ƶ�
						}
					}

					//�����û������е�[row][col]
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
			//�����´μ������¾���
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
	LocalWarp������,����ת������,warp_imgΪ����֮���ͼ��
*/
vector<vector<Coordinate>> LocalWarpping::Local_warp(CVMat src, CVMat& wrap_img, CVMat mask) 
{
	vector<vector<Coordinate>> displacementMap = Get_Local_warp_displacement(src, mask);//��ȡ�û�����
	for (int row = 0; row < src.rows; row++)
	{
		for (int col = 0; col < src.cols; col++)
		{
			Coordinate displacement = displacementMap[row][col];
			colorPixel pixel = src.at<colorPixel>(row + displacement.row, col + displacement.col);
			wrap_img.at<colorPixel>(row, col) = pixel;//����ת��������������֮���ͼ��
		}
	}
	return displacementMap;
}