#include "stdafx.h"
#include "vango.h"
#include <map>
#include "time.h"
using std::map;
using std::cout;
using std::endl;


void vango::bangbang(vector<FShell*> &species)
{
	for (auto &shell : species)
	{
		delete shell;
		shell = nullptr;
	}
}

void vango::bang(vector<FShell*> &species, int size)
{
	while (size--)
	{
		species.push_back(new FShell(clock() + size));
	}
}

// ����ͼ�����ƶ�
double vango::calc_sim(const Mat& mat1, const Mat& mat2)
{
	//return getMSSIM(mat1, mat2);

	double marks = 0;

	int nr = mat1.rows; // number of rows
	int nc = mat1.cols * mat1.channels(); // total number of elements per line
	int step = mat1.step; // effective width
						  // get the pointer to the image buffer
	uchar *data1 = mat1.data;
	uchar *data2 = mat2.data;
	for (int j = 0; j<nr; j++) {
		for (int i = 0; i<nc; i++) {
			marks += log(1 + abs(*(data1 + i) - *(data2 + i)));
			//*(data1 + i) = *data&mask + div / 2;
		}
		data1 += step;  // next line
		data2 += step;
	}

	return marks;
}

void vango::transform_dna(const vector<trangle>& DNA, Mat &result)
{
	Mat canvas(result.size(),result.type(),Scalar::all(255));
	vector<vector<Point>> poi;
	for (auto &s : DNA)
	{
		poi.clear();
		vector<Point>p{ s.a, s.b, s.c };
		poi.push_back(p);
		fillPoly(canvas, poi, s.color);
		addWeighted(result, 0.985, canvas, 0.015, 0, result);
	}

}

double vango::meansure(FShell* shell, const Mat &target)
{
	Mat bg(target.size(), target.type(), Scalar::all(255));
	transform_dna(shell->DNAs, bg);
	return calc_sim(bg, target);
}

void vango::birth_live(vector<FShell*>& parents, int ratio, const Mat& target)
{
	if (ratio <= 0)return;// û�з���

	// ����һ�¸�������,ȡǰһ��ͺ�һ��ɶ��ӽ�
	random_shuffle(parents.begin(), parents.end());
	int species_size = parents.size() / 2;
	vector<trangle> lh, rh, lh1, rh1;
	vector<FShell> next_genera;
	for (int i = 0; i < species_size; i++)
	{
		// һ��ֻ��һ��, ��ratio��
		for (int j = 0; j < ratio; j++)
		{
			parents[i]->division(lh, rh);
			parents[i+ species_size]->division(lh1, rh1);
			// �ӽ�
			next_genera.emplace_back(FShell(lh, rh1));
			next_genera.emplace_back(FShell(lh1, rh));
		}
	}
	vector<mark_shell> genera_rank;
	for (auto &g : next_genera)
	{
		genera_rank.emplace_back(&g, meansure(&g, target));
	}
	sort(genera_rank.begin(), genera_rank.end());
	for (int i = 0;i<parents.size();i++)
	{
		parents[i]->DNAs = genera_rank[i].shell->DNAs;

	}
}

void vango::Running(vector<FShell*>& init_species, const int round, const Mat& target)
{
	auto t = clock();
	vector<int>output{ 10,100,300,1000,5000,10000 };
	int total = round;
	while (total--)
	{
		cout << "���ܽ��̻���: " << total << " ��" << endl;

		birth_live(init_species, 3, target);
		for (auto & o : output)
		{
			if (o == round - total && total !=0)
			{
				// ���
				show_result(init_species, round - total);
				break;
			}
		}
	}
	// ���������ɺ������
	show_result(init_species, round, true);
	cout << "���!!!"<< "��ʱ��:"<<(clock()-t)/1000.f << endl;
	waitKey(0);
}

void vango::show_result(const vector<FShell*>& result, int NO, bool wait)
{
	Mat bg(target.size(), target.type(), Scalar::all(255));
	for (int i = 0; i < result.size()*0.1; i++)
	{
		bg.setTo(Scalar::all(255));
		transform_dna(result[i]->DNAs, bg);
		if(wait)imshow(format("%d", i), bg);
		imwrite(format("GenareResult/G_%d_NO_%d_.jpg", NO, i), bg);
	}
}

void vango::revolution()
{
	// ����Ŀ����Ⱥ������̬ͼ��
	target = imread("fire.jpg");

	vector<FShell*> species;

	// Ů洴���������Ⱥ
	bang(species, 100);

	// ����
	Running(species, 10, target);

	// ��������
	bangbang(species);
}
