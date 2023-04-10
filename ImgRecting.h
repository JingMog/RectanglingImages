#pragma once

#include <cmath>
#include "GlobalWarpping.h"
#include "lsd.h"
#include "GL/glut.h"
#include "LocalWarpping.h"

class ImgRecting
{
private:
	GLuint texGround;
	vector<vector<CoordinateDouble>> outputmesh;
	vector<vector<CoordinateDouble>> mesh;
	CVMat img;

public:
	GLuint matToTexture(cv::Mat mat, GLenum minFilter, GLenum magFilter, GLenum wrapFilter);
	void display();
	void runImgRecting(string imgPath);

};

