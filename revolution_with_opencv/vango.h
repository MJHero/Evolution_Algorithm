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
	// ��
	Point a, b, c;
	// ��ɫ
	Scalar color;
};

struct FShell
{
	vector<trangle> DNAs;
	RNG rng;
	// Ů洴����ͳ�ʼ��
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

	// ������ϳ�ʼ��
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

	// �����������
	void division(vector<trangle>& lh, vector<trangle> &rh, float ratio)
	{
		// ����ʱ���ܱ���
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
	double mark;// ����
};

class vango
{
public:

	// ��Ȼѡ�����
	void revolution();

private:

	Mat target;
private:

	/*	��ը��ʼ
	*	@species				��Ⱥ
	*	@size					��Ⱥ����
	*/
	void bang(vector<FShell*> &species, int size);

	/*	�������
	*/
	void bangbang(vector<FShell*> &species);

	/*	�������ַ���
	*	@init_species			��ʼ��Ⱥ
	*	@round					���ܴ���
	*	@target					Ŀ������
	*/
	void Running(vector<FShell*>& init_species, const int round, const Mat& target);

	/*	����ͬ��,���Է�ֳ
	*/
	void birth_live(vector<FShell*>& parents, int ratio, const Mat& target);

	/*	����һ��shell ��target֮������ƶȷ���
	* 
	*/
	double meansure(FShell* shell, const Mat &target);

	/*	����ͼ�����ƶ�
	*/ 
	double calc_sim(const Mat& mat1, const Mat& mat2);

	/*	��һ��DNA����һ��ͼ
	*/
	void transform_dna(const vector<trangle>& DNA, Mat &result);

	void show_result(const vector<FShell*>& result, int NO,bool wait = false);

};
