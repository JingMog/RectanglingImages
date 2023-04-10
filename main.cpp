#include <iostream>
#include "LocalWarpping.h"

#include "ImgRecting.h"

using namespace std;

int main()
{
	double Time = (double)cvGetTickCount();

	ImgRecting rect;
	string imgPath = "C:/Users/jingmo/Pictures/CImages/Rectangling/1.jpg";
	rect.runImgRecting(imgPath);//运行ImgRecting

	/*
	LocalWarpping local;
	CVMat img = cv::imread(imgPath);
	cv::namedWindow("img");
	cv::imshow("img", img);

	Border direction;
	CVMat warpped_img = CVMat::zeros(img.size(), CV_8UC3);
	CVMat mask = Mask_contour(img);

	vector<vector<Coordinate>> displacementMap = local.Local_warp(img, warpped_img, mask);

	cv::namedWindow("warpImg");
	cv::imshow("warpImg", warpped_img);


	Time = (double)cvGetTickCount() - Time;
	cout << "Warp Img Time:  " << Time / ((cvGetTickFrequency() * 1000) * 1000) << "s" << endl;
	//printf("Warp Img Time = %gms\n", Time / (cvGetTickFrequency() * 1000));//输出warp Img时间

	cv::waitKey(0);
	*/
	return 0;
}