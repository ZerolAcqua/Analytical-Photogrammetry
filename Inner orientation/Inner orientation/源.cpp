#define _CRT_SECURE_NO_WARNINGS 1

#include "stdio.h"
#include "Matrix.h"
#include "math.h"


// 内定向编程
int main()
{
	// 框标坐标x',y'
	Matrix matMatkCoor = { {-106.0010,-106.0040},
							{106.0020,-106.0030},
							{105.9990,106.0020},
							{-106.0000,106.0020} };

	// 量测坐标x,y
	Matrix matMeasureCoor = { {446.625,594.750},
							{10546.625,586.313 },
							{10556.000,10687.438},
							{456.000,10696.500} };

	// 量测了四个框标，可以采用仿射变换进行内定向
	// x'=a_0+a_1x+a_2y
	// y'=b_0+b_1x+b_2y

	Matrix matA = Matrix::zeros(8, 6);
	// 初值全部取为0,因为是线性的，不用迭代
	Matrix matX = Matrix::zeros(6, 1);	// a_0 a_1 a_2 b_0 b_1 b_2
	Matrix matL = Matrix::zeros(8, 1);

	for (int i = 0; i < 4; i++)
	{
		matA[i * 2][0] = 1;
		matA[i * 2][1] = matMeasureCoor[i][0];
		matA[i * 2][2] = matMeasureCoor[i][1];
		matA[i * 2 + 1][3] = 1;
		matA[i * 2 + 1][4] = matMeasureCoor[i][0];
		matA[i * 2 + 1][5] = matMeasureCoor[i][1];

		//						x		-		x^0
		matL[i * 2][0] = matMatkCoor[i][0] - 0;
		matL[i * 2 + 1][0] = matMatkCoor[i][1] - 0;
	}

	matX = (matA.transpose() * matA).inverse() * (matA.transpose() * matL);

	// 计算残差
	Matrix matEpsilon = matA * matX - matL;	//残差
	Matrix matEpsilonX = Matrix::zeros(4, 1);
	Matrix matEpsilonY = Matrix::zeros(4, 1);
	for (int i = 0; i < 4; i++)
	{
		matEpsilonX[i][0] = matEpsilon[i * 2][0];
		matEpsilonY[i][0] = matEpsilon[i * 2 + 1][0];
	}

	FILE* fout = fopen("innerOrientation_result.txt", "w");
	std::cout << "++++++++++++++++++++++++++++++++++++++++++++++\n";
	std::cout <<"内定向正算参数："<< matX;
	fprintf(fout, "++++++++++++++++++++++++++++++++++++++++++++++\n");
	fprintf(fout, "内定向正算参数：\na0\t a1\t a2\t b0\t b1\t b2\n%7e\t%7e\t%7e\t%7e\t%7e\t\n",
		matX[0][0], matX[1][0], matX[2][0], matX[3][0], matX[4][0], matX[5][0]);
	fprintf(fout, "残差：\nx1'\t y1'\t x2'\t x2'\t x3'\t y3'\t x4'\t y4'\t\n%7e\t%7e\t%7e\t%7e\t%7e\t\n",
		matEpsilon[0][0], matEpsilon[1][0], matEpsilon[2][0], matEpsilon[3][0],
		matEpsilon[4][0], matEpsilon[5][0], matEpsilon[6][0], matEpsilon[7][0]);
	printf("rms x(mm)：\n");
	printf("%f\n", sqrt((matEpsilonX.transpose() * matEpsilonX)[0][0] / 4));
	printf("rms y(mm)：\n");
	printf("%f\n", sqrt((matEpsilonY.transpose() * matEpsilonY)[0][0] / 4));

	fprintf(fout, "rms x(mm)：\n");
	fprintf(fout, "%f\n", sqrt((matEpsilonX.transpose() * matEpsilonX)[0][0] / 4));
	fprintf(fout, "rms y(mm)：\n");
	fprintf(fout, "%f\n", sqrt((matEpsilonY.transpose() * matEpsilonY)[0][0] / 4));


	// ------------------------------------------------------------------------------------
	// 反算参数，代码相同，只需要把“框标坐标”和“量测坐标互换”即可
	
	// 框标坐标x',y'
	matMeasureCoor = { {-106.0010,-106.0040},
						{106.0020,-106.0030},
						{105.9990,106.0020},
						{-106.0000,106.0020} };

	// 量测坐标x,y
	 matMatkCoor = { {446.625,594.750},
						{10546.625,586.313 },
						{10556.000,10687.438},
						{456.000,10696.500} };

	// 量测了四个框标，可以采用仿射变换进行内定向
	// x'=a_0+a_1x+a_2y
	// y'=b_0+b_1x+b_2y

	matA = Matrix::zeros(8, 6);
	// 初值全部取为0,因为是线性的，不用迭代
	matX = Matrix::zeros(6, 1);	// a_0 a_1 a_2 b_0 b_1 b_2
	matL = Matrix::zeros(8, 1);

	for (int i = 0; i < 4; i++)
	{
		matA[i * 2][0] = 1;
		matA[i * 2][1] = matMeasureCoor[i][0];
		matA[i * 2][2] = matMeasureCoor[i][1];
		matA[i * 2 + 1][3] = 1;
		matA[i * 2 + 1][4] = matMeasureCoor[i][0];
		matA[i * 2 + 1][5] = matMeasureCoor[i][1];

		//						x		-		x^0
		matL[i * 2][0] = matMatkCoor[i][0] - 0;
		matL[i * 2 + 1][0] = matMatkCoor[i][1] - 0;
	}

	matX = (matA.transpose() * matA).inverse() * (matA.transpose() * matL);

	// 计算残差
	matEpsilon = matA * matX - matL;	//残差
	for (int i = 0; i < 4; i++)
	{
		matEpsilonX[i][0] = matEpsilon[i * 2][0];
		matEpsilonY[i][0] = matEpsilon[i * 2 + 1][0];
	}

	std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++\n";
	std::cout << "内定向反算参数：" << matX;
	fprintf(fout, "++++++++++++++++++++++++++++++++++++++++++++++\n");
	fprintf(fout, "内定向反算参数：\na0'\t a1'\t a2'\t b0'\t b1'\t b2'\n%7e\t%7e\t%7e\t%7e\t%7e\t\n",
		matX[0][0], matX[1][0], matX[2][0], matX[3][0], matX[4][0], matX[5][0]);
	fprintf(fout, "残差：\nx1\t y1\t x2\t x2\t x3\t y3\t x4\t y4\t\n%7e\t%7e\t%7e\t%7e\t%7e\t\n",
		matEpsilon[0][0], matEpsilon[1][0], matEpsilon[2][0], matEpsilon[3][0],
		matEpsilon[4][0], matEpsilon[5][0], matEpsilon[6][0], matEpsilon[7][0]);
	printf("rms x：\n");
	printf("%f\n", sqrt((matEpsilonX.transpose() * matEpsilonX)[0][0] / 4));
	printf("rms y：\n");
	printf("%f\n", sqrt((matEpsilonY.transpose() * matEpsilonY)[0][0] / 4));

	fprintf(fout, "rms x：\n");
	fprintf(fout, "%f\n", sqrt((matEpsilonX.transpose() * matEpsilonX)[0][0] / 4));
	fprintf(fout, "rms y：\n");
	fprintf(fout, "%f\n", sqrt((matEpsilonY.transpose() * matEpsilonY)[0][0] / 4));
}