#include "ImgRecting.h"

GLuint ImgRecting::matToTexture(cv::Mat mat, GLenum minFilter = GL_LINEAR, GLenum magFilter = GL_LINEAR, GLenum wrapFilter = GL_REPEAT)
{
	//cv::flip(mat, mat, 0);
	// Generate a number for our textureID's unique handle
	GLuint textureID;
	glGenTextures(1, &textureID);

	// Bind to our texture handle
	glBindTexture(GL_TEXTURE_2D, textureID);

	// Catch silly-mistake texture interpolation method for magnification
	if (magFilter == GL_LINEAR_MIPMAP_LINEAR ||
		magFilter == GL_LINEAR_MIPMAP_NEAREST ||
		magFilter == GL_NEAREST_MIPMAP_LINEAR ||
		magFilter == GL_NEAREST_MIPMAP_NEAREST)
	{
		//cout << "You can't use MIPMAPs for magnification - setting filter to GL_LINEAR" << endl;
		magFilter = GL_LINEAR;
	}

	// Set texture interpolation methods for minification and magnification
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, minFilter);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, magFilter);

	// Set texture clamping method
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, wrapFilter);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, wrapFilter);

	// Set incoming texture format to:
	// GL_BGR       for CV_CAP_OPENNI_BGR_IMAGE,
	// GL_LUMINANCE for CV_CAP_OPENNI_DISPARITY_MAP,
	// Work out other mappings as required ( there's a list in comments in main() )
	GLenum inputColourFormat = GL_BGR_EXT;
	if (mat.channels() == 1)
	{
		inputColourFormat = GL_LUMINANCE;
	}

	// Create the texture
	glTexImage2D(GL_TEXTURE_2D,     // Type of texture
		0,                 // Pyramid level (for mip-mapping) - 0 is the top level
		GL_RGB,            // Internal colour format to convert to
		mat.cols,          // Image width  i.e. 640 for Kinect in standard mode
		mat.rows,          // Image height i.e. 480 for Kinect in standard mode
		0,                 // Border width in pixels (can either be 1 or 0)
		inputColourFormat, // Input image format (i.e. GL_RGB, GL_RGBA, GL_BGR etc.)
		GL_UNSIGNED_BYTE,  // Image data type
		mat.ptr());        // The actual image data itself

						   // If we're using mipmaps then generate them. Note: This requires OpenGL 3.0 or higher
	return textureID;
}



void ImgRecting::display()
{
	glLoadIdentity();
	// 清除屏幕
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glBindTexture(GL_TEXTURE_2D, texGround);
	for (int row = 0; row < 20; row++) {
		for (int col = 0; col < 20; col++) {
			CoordinateDouble& coord = outputmesh[row][col];
			CoordinateDouble& localcoord = mesh[row][col];
			coord.row /= img.rows;
			coord.col /= img.cols;
			coord.row -= 0.5;
			coord.col -= 0.5;
			coord.row *= 2;
			coord.col *= 2;
			coord.row = clamp(coord.row, -1, 1);
			coord.col = clamp(coord.col, -1, 1);
			//cout << coord << " ";

			localcoord.row /= img.rows;
			localcoord.col /= img.cols;
			localcoord.row = clamp(localcoord.row, 0, 1);
			localcoord.col = clamp(localcoord.col, 0, 1);
		}
		// cout << endl;
	}
	//system("pause");

	for (int row = 0; row < 19; row++) {
		for (int col = 0; col < 19; col++) {
			CoordinateDouble local_left_top = mesh[row][col];
			CoordinateDouble local_right_top = mesh[row][col + 1];
			CoordinateDouble local_left_bottom = mesh[row + 1][col];
			CoordinateDouble local_right_bottom = mesh[row + 1][col + 1];


			CoordinateDouble global_left_top = outputmesh[row][col];
			CoordinateDouble global_right_top = outputmesh[row][col + 1];
			CoordinateDouble global_left_bottom = outputmesh[row + 1][col];
			CoordinateDouble global_right_bottom = outputmesh[row + 1][col + 1];

			//
			glBegin(GL_QUADS);
			glTexCoord2d(local_right_top.col, local_right_top.row); glVertex3d(global_right_top.col, -1 * global_right_top.row, 0.0f);
			glTexCoord2d(local_right_bottom.col, local_right_bottom.row); glVertex3d(global_right_bottom.col, -1 * global_right_bottom.row, 0.0f);
			glTexCoord2d(local_left_bottom.col, local_left_bottom.row);	glVertex3d(global_left_bottom.col, -1 * global_left_bottom.row, 0.0f);
			glTexCoord2d(local_left_top.col, local_left_top.row); glVertex3d(global_left_top.col, -1 * global_left_top.row, 0.0f);
			glEnd();

		}
	}
	/*
	int row = 18;
	int col = 18;
	CoordinateDouble local_left_top = mesh[row][col];
	CoordinateDouble local_right_top = mesh[row][col + 1];
	CoordinateDouble local_left_bottom = mesh[row + 1][col];
	CoordinateDouble local_right_bottom = mesh[row + 1][col + 1];

	CoordinateDouble global_left_top = outputmesh[row][col];
	CoordinateDouble global_right_top = outputmesh[row][col + 1];
	CoordinateDouble global_left_bottom = outputmesh[row + 1][col];
	CoordinateDouble global_right_bottom = outputmesh[row + 1][col + 1];
	glBegin(GL_QUADS);
	glTexCoord2d(local_right_top.col, local_right_top.row);
	glVertex3d(1, 0.7, 0.0f);
	glTexCoord2d(local_right_bottom.col, local_right_bottom.row);
	glVertex3d(1, -0.7, 0.0f);
	glTexCoord2d(local_left_bottom.col, local_left_bottom.row);
	glVertex3d(-1, -1, 0.0f);
	glTexCoord2d(local_left_top.col, local_left_top.row);
	glVertex3d(-1, 1, 0.0f);
	glEnd();
	*/
	glutSwapBuffers();
}


void ImgRecting::runImgRecting(string imgPath)
{
	double Time = (double)cvGetTickCount();//记录程序开始执行的时间

	img = cv::imread(imgPath);//读取图片

	cv::namedWindow("origin Img");
	cv::imshow("origin Img", img);

	CVMat scaled_img;
	cv::resize(img, scaled_img, cv::Size(0, 0), 0.5, 0.5);//对图像进行缩放
	Config config(scaled_img.rows, scaled_img.cols, 20, 20);//设置参数
	CVMat mask = Mask_contour(scaled_img);//获取mask
	CVMat tmpmask;
	mask.copyTo(tmpmask);
	CVMat wrapped_img = CVMat::zeros(scaled_img.size(), CV_8UC3);//保存warp之后得图像
	LocalWarpping localWrap;//实例化LocalWarp对象

	//1.局部变形
	vector<vector<Coordinate>> displacementMap = localWrap.Local_warp(scaled_img, wrapped_img, tmpmask);//localwarp，获取图像变形

	cv::resize(wrapped_img, wrapped_img, cv::Size(0, 0), 2, 2, cv::INTER_AREA);
	cv::namedWindow("warped Img", cv::WINDOW_AUTOSIZE);
	cv::imshow("warped Img", wrapped_img);
	

	cout << "origin img size(): " << scaled_img.size() << endl;
	cout << "warpped_img.size():  " << wrapped_img.size() << endl;

	//2.在变形之后得图像上放置标准矩形网格
	mesh = localWrap.get_rectangle_mesh(scaled_img, config);
	//drawmesh(wrapped_img, mesh, config);
	

	//3.将标准网格根据变形矩阵转化到原图像上
	localWrap.warp_mesh_back(mesh, displacementMap, config);
	//drawmesh(scaled_img, mesh, config);
	cout << "局部变形wrap back完成" << endl;

	//4.全局调整
	GlobalWarpping global(config, scaled_img);
	SpareseMatrixD_Row shape_energy = global.get_shape_mat(mesh);//获取shape Energy Mat
	cout << "get shape energy" << endl;
	SpareseMatrixD_Row Q = global.get_vertex_to_shape_mat(mesh);//获取Q矩阵
	pair<SpareseMatrixD_Row, VectorXd> pair_dvec_B = global.get_boundary_mat(mesh);
	cout << "get border constraint" << endl;

	vector<pair<int, double>>id_theta;//保存每一个bin的theta
	vector<LineD> line_flatten;
	vector<double> rotate_theta;
	vector<vector<vector<LineD>>> LineSeg = global.init_line_seg(mask, line_flatten, mesh, id_theta, rotate_theta);//用网格对所有LSD检测到的直线进行分割

	// 开始优化迭代
	for (int iter = 1; iter <= 10; iter++) 
	{
		cout << "iter: " << iter << endl;
		int linenum = 0;
		vector<pair<MatrixXd, MatrixXd>> BilinearVec;//need to update
		vector<bool> bad;
		SpareseMatrixD_Row line_energy = global.get_line_mat(mask, mesh, rotate_theta, LineSeg, BilinearVec, linenum, bad);//每次迭代需要重新获取line_energy
		cout << "get line energy, line num: " << linenum << endl;

		double Nq = config.meshQuadRow * config.meshQuadCol;//网格点的个数
		double lambdaB = INF;//论文中公式(9)中的lambdaB
		double lambdaL = 100;//论文中公式(9)中的lambdaL
		SpareseMatrixD_Row shape = (1 / Nq) * (shape_energy * Q);//shape
		SpareseMatrixD_Row boundary = lambdaB * pair_dvec_B.first;//boundary
		SpareseMatrixD_Row line = (lambdaL / linenum) * (line_energy * Q);//line

		SpareseMatrixD_Row K = row_stack(shape, line);
		SpareseMatrixD_Row K2 = row_stack(K, boundary);//将shape,line,boundary合并成一个矩阵

		//print_sparse_mat_row(shape_energy,1);
		//system("pause");
		VectorXd B = pair_dvec_B.second;//B
		VectorXd BA = VectorXd::Zero(K2.rows());
		BA.tail(B.size()) = lambdaB * B;//将LambdaB*B放到BA下方

		SparseMatrixD K2_trans = K2.transpose();//K2^T
		SparseMatrixD A = K2_trans * K2;//A为所有能量的平方即||E(V,threta)||^2

		VectorXd b = K2_trans * BA;//K2^T * BA
		VectorXd x;

		//A=K2^T * K2, b=K2^T * BA,求解bx = A
		CSolve* p_A = new CSolve(A);//Eigen库中的数值求解器
		x = p_A->solve(b);

		//update theta
		outputmesh = vector_to_mesh(x, config);//将求解得到的x转化为更新之后的网格

		int tmplinenum = -1;
		VectorXd thetagroup = VectorXd::Zero(50);
		VectorXd thetagroupcnt = VectorXd::Zero(50);
		for (int row = 0; row < config.meshQuadRow; row++)
		{
			for (int col = 0; col < config.meshQuadCol; col++)
			{
				//遍历每一个网格点
				//更新网格点中的直线
				vector<LineD> linesegInquad = LineSeg[row][col];//该网格中的线段
				int QuadID = row * config.meshQuadCol + col;//网格一维ID
				if (linesegInquad.size() == 0) 
				{
					continue;//如果该网格中没有线段,直接下一个网格
				}
				else 
				{
					VectorXd S = global.get_vertices(row, col, outputmesh);
					for (int k = 0; k < linesegInquad.size(); k++) 
					{
						tmplinenum++;
						//cout << tmplinenum<<endl;
						if (bad[tmplinenum] == true) 
						{
							//cout << "bad is true " << endl;
							continue;
						}
						//cout << tmplinenum;
						pair<MatrixXd, MatrixXd> Bstartend = BilinearVec[tmplinenum];//获取逆双线性插值的权重
						MatrixXd start_W_mat = Bstartend.first;//线段起始点Weight mat
						MatrixXd end_W_mat = Bstartend.second;//线段终点Weight mat
						Vector2d newstart = start_W_mat * S;//对线段起点进行更新
						Vector2d newend = end_W_mat * S;//对线段终点进行更新

						double theta = atan((newstart(1) - newend(1)) / (newstart(0) - newend(0)));//新的角度
						double deltatheta = theta - id_theta[tmplinenum].second;//角度的变化值
						if (isnan(id_theta[tmplinenum].second) || isnan(deltatheta)) 
						{
							continue;
						}
						//保证deltatheta 在[-pi,pi]之间
						if (deltatheta > (PI / 2)) 
						{
							deltatheta -= PI;
						}
						if (deltatheta < (-PI / 2)) 
						{
							deltatheta += PI;
						}
						thetagroup(id_theta[tmplinenum].first) += deltatheta;//更新theta
						thetagroupcnt(id_theta[tmplinenum].first) += 1;
					}
				}
			}
		}

		//计算theta的均值
		for (int ii = 0; ii < thetagroup.size(); ii++) 
		{
			thetagroup(ii) /= thetagroupcnt(ii);
		}
		//更新rotate_theta数组
		for (int ii = 0; ii < rotate_theta.size(); ii++) {
			rotate_theta[ii] = thetagroup[id_theta[ii].first];
		}
		/*if (iter % 5 == 0)
		{
			vector<vector<CoordinateDouble>> outputmesh = vector_to_mesh(x, config);
			drawmesh(scaled_img, outputmesh, config);
		}*/
		
		//system("pause");
	}

	//cout << x;
	//system("pause");
	
	enlarge_mesh(mesh, 2, config);
	enlarge_mesh(outputmesh, 2, config);
	CVMat outputimg = CVMat::zeros(img.size(), CV_32FC3);
	CVMat ouputcnt = CVMat::zeros(img.size(), CV_32FC3);


	for (int row = 0; row < config.meshQuadRow; row++) 
	{
		cout << "row: " << row << endl;
		for (int col = 0; col < config.meshQuadCol; col++)
		{
			VectorXd Vq = global.get_vertices(row, col, outputmesh);//x0,y0,x1,y1(Vq)
			VectorXd Vo = global.get_vertices(row, col, mesh);//x0,y0
			double col_len = max(Vq(0), max(Vq(2), max(Vq(4), Vq(6)))) - min(Vq(0), min(Vq(2), min(Vq(4), Vq(6))));//
			double row_len = max(Vq(1), max(Vq(3), max(Vq(5), Vq(7)))) - min(Vq(1), min(Vq(3), min(Vq(5), Vq(7))));//
			double col_step = 1 / (4 * col_len);
			double row_step = 1 / (4 * row_len);
			//cout << "col_step: " << col_step << endl;
			//cout << "row_step: " << row_step << endl;

			//system("pause");
			for (double i = 0; i < 1; i += row_step)
			{
				for (double j = 0; j < 1; j += col_step)
				{
					double v1w = 1 - i - j + i * j;
					double v2w = j - i * j;
					double v3w = i - i * j;
					double v4w = i * j;
					MatrixXd matt(2, 8);
					matt << v1w, 0, v2w, 0, v3w, 0, v4w, 0,
						0, v1w, 0, v2w, 0, v3w, 0, v4w;
					VectorXd pout = matt * Vq;
					VectorXd pref = matt * Vo;
					if (int(pout(1)) >= 0 && int(pout(0)) >= 0 && int(pout(1)) < img.rows && int(pout(0)) < img.cols)
					{
						colorPixel pixel = img.at<colorPixel>(int(pref(1)), int(pref(0)));
						cv::Vec3f pixelf = cv::Vec3f(float(pixel[0]), float(pixel[1]), float(pixel[2]));
						outputimg.at<cv::Vec3f>(int(pout(1)), int(pout(0))) = outputimg.at<cv::Vec3f>(int(pout(1)), int(pout(0))) + pixelf;
						ouputcnt.at<cv::Vec3f>(int(pout(1)), int(pout(0))) += cv::Vec3f(1, 1, 1);
					}
				}
			}
		}
	}

	CVMat finaloutput = outputimg / (255 * ouputcnt);
	CVMat meshimg;
	finaloutput.copyTo(meshimg);

	drawmesh(meshimg, outputmesh, config);
	//fill_image(finaloutput);
	cv::namedWindow("Border", CV_WINDOW_AUTOSIZE);
	cv::imshow("Border", finaloutput);

	Time = (double)cvGetTickCount() - Time	;

	::printf("run time = %gms\n", Time / (cvGetTickFrequency() * 1000));//毫秒
	cv::waitKey(0);

	/*
	//glut
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA);
	glutInitWindowPosition(100, 100);
	glutInitWindowSize(img.cols, img.rows);
	glutCreateWindow("OpenGL纹理");
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_TEXTURE_2D);    // 启用纹理
	texGround = matToTexture(img);
	glutDisplayFunc(&display);   //注册函数
	Time = (double)cvGetTickCount() - Time;

	printf("run time = %gms\n", Time / (cvGetTickFrequency() * 1000));//毫秒
	glutMainLoop(); //循环调用
	*/
	return;
}




