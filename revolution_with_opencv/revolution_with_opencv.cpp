// revolution_with_opencv.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "vango.h"
#include <iostream>
using namespace std;
void main()
{
	//Mat temp = imread("fire.jpg");
	//RNG g_rng(clock());
	//Mat bg(250,250,CV_8UC4);
	//vector<Point>p{
	//	Point(g_rng.uniform(0, 250), g_rng.uniform(0, 250)),
	//	Point(g_rng.uniform(0, 250), g_rng.uniform(0, 250)),
	//	Point(g_rng.uniform(0, 250), g_rng.uniform(0, 250))
	//};
	//
	//vector<vector<Point>> poi;
	//poi.push_back(p);
	//auto t = clock();
	//for (int i = 0; i < 1000000; i++)
	//{
	//	fillPoly(bg, poi, Scalar(g_rng.uniform(0, 255), g_rng.uniform(0, 255), g_rng.uniform(0, 255), 0));
	//
	//}
	//cout << clock() - t << endl;
	//
	//imshow("", temp);
	//waitKey(0);

	vango temp;
	temp.revolution();
}

