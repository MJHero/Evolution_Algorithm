#pragma once
#include <opencv2/opencv.hpp>
#include <iostream>
#include <vector>

using namespace cv;
using std::vector;


struct trangle
{
	trangle() {};
	trangle(const Point &a_, const Point& b_, const Point& c_, const Scalar &color_)
		:a(a_), b(b_), c(c_), color(color_) {};
	// 点
	Point a, b, c;
	// 颜色
	Scalar color;
};

struct FShell
{
	vector<trangle> DNAs;
	RNG rng;
	// 女娲创世纪初始化
	FShell(int seed):rng(seed)
	{
		for (int i = 0; i < 100; i++)
		{
			DNAs.emplace_back(
				Point(rng.uniform(0, 250), rng.uniform(0, 250)),
				Point(rng.uniform(0, 250), rng.uniform(0, 250)),
				Point(rng.uniform(0, 250), rng.uniform(0, 250)),
				Scalar(rng.uniform(0, 250), rng.uniform(0, 250), rng.uniform(0, 250))
			);
		}
	}

	// 基因组合初始化
	FShell(const vector<trangle>& lh, const vector<trangle> &rh)
	{
		DNAs.clear();
		//DNAs
		for (auto & t : lh)
		{
			DNAs.emplace_back(t);
		}
		for (auto & t : rh)
		{
			DNAs.emplace_back(t);
		}
	}

	// 基因减数分裂
	void division(vector<trangle>& lh, vector<trangle> &rh, float ratio)
	{
		// 分裂时可能变异
		lh.clear();
		rh.clear();
		random_shuffle(DNAs.begin(),DNAs.end());
		for (int i = 0; i < 50; i++)
		{
			auto temp = DNAs[i];
			int rnd = rng.uniform(0, 10000);
			if (rnd > 10000*(1-ratio))
			{
				switch (rnd%2)
				{
				case 0:
					temp.a = Point(rng.uniform(0, 250), rng.uniform(0, 250));
					break;
				case 1:
					temp.b = Point(rng.uniform(0, 250), rng.uniform(0, 250));
					break;
				case 2:
					temp.c = Point(rng.uniform(0, 250), rng.uniform(0, 250));
					break;
				case 3:
					temp.color = Scalar(rng.uniform(0, 250), rng.uniform(0, 250), rng.uniform(0, 250));
					break;
				}
			}
			lh.emplace_back(temp);
		}
		for (int i = 50; i < 100; i++)
		{
			auto temp = DNAs[i];
			int rnd = rng.uniform(0, 10000);
			if (rnd > 10000 * (1 - ratio))
			{
				switch (rnd % 2)
				{
				case 0:
					temp.a = Point(rng.uniform(0, 250), rng.uniform(0, 250));
					break;
				case 1:
					temp.b = Point(rng.uniform(0, 250), rng.uniform(0, 250));
					break;
				case 2:
					temp.c = Point(rng.uniform(0, 250), rng.uniform(0, 250));
					break;
				case 3:
					temp.color = Scalar(rng.uniform(0, 250), rng.uniform(0, 250), rng.uniform(0, 250));
					break;
				}
			}
			rh.emplace_back(DNAs[i]);
		}
	}
};

struct mark_shell
{
	mark_shell(FShell* s_, double mark_) :shell(s_), mark(mark_) {}
	bool operator<(const mark_shell&lh) { return mark < lh.mark; }

	FShell* shell;
	double mark;// 分数
};

class vango
{
public:

	// 自然选择进化
	void revolution();

private:

	Mat target;
private:

	/*	大爆炸伊始
	*	@species				种群
	*	@size					种群数量
	*/
	void bang(vector<FShell*> &species, int size);

	/*	宇宙毁灭
	*/
	void bangbang(vector<FShell*> &species);

	/*	生物物种繁衍
	*	@init_species			初始种群
	*	@round					繁衍代数
	*	@target					目标样子
	*/
	void Running(vector<FShell*>& init_species, const int round, const Mat& target);

	/*	雌雄同体,无性繁殖
	*/
	void birth_live(vector<FShell*>& parents, int ratio, const Mat& target);

	/*	计算一个shell 和target之间的相似度分数
	* 
	*/
	double meansure(FShell* shell, const Mat &target);

	/*	计算图像相似度
	*/ 
	double calc_sim(const Mat& mat1, const Mat& mat2);

	/*	把一个DNA画成一张图
	*/
	void transform_dna(const vector<trangle>& DNA, Mat &result);

	void show_result(const vector<FShell*>& result, int NO,bool wait = false);

};
