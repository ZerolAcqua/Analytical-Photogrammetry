#define  _CRT_SECURE_NO_WARNINGS 1
#include <stdio.h>
#include "Matrix.h"
#include "math.h"

#define ITER_TIMES 10


// 后方交会编程
int main()
{
	Matrix::setPrecise(10);
	////// 数据读取部分 //////
	// 读取Data.txt
	FILE* fp = fopen("resection_data.txt", "r");
	if (fp == NULL) {
		perror("打开文件时发生错误");
		return -1;
	}
	char* data = new char[512];
	int n = 0;	// 控制点个数

	fgets(data, 512, fp);
	sscanf(data, "%d", &n);

	int id = 0;
	double* param = new double[n * 5];
	for (int i = 0; fgets(data, 512, fp) && i < n; i++)
	{
		// ID x y X Y Z
		sscanf(data, "%d %lf %lf %lf %lf %lf "
			, &id, param + 5 * i, param + 1 + 5 * i, param + 2 + 5 * i, param + 3 + 5 * i, param + 4 + 5 * i
		);
	}
	fclose(fp);
	delete[]data;

	// 外方位元素 X_s,Y_s,Z_s,phi,omega,kappa
	Matrix matOuterOrientElem = Matrix::zeros(6, 1);
	// 解法方程所得的改正数dX_s,dY_s,dZ_s,dphi,domega,dkappa
	Matrix matX = Matrix::zeros(6, 1);
	// 误差方程系数矩阵
	Matrix matA = Matrix::zeros(n * 2, 6);
	// 误差程常数矩阵
	Matrix matL = Matrix::zeros(n * 2, 1);
	// 旋转矩阵
	//	a_1	a_2	a_3	
	//	b_1	b_2	b_3
	//	c_1	c_2	c_3
	Matrix matR = Matrix::zeros(3, 3);

	double phi = 0;
	double omega = 0;
	double kappa = 0;
	double X_s = 0;
	double Y_s = 0;
	double Z_s = 0;
	double a_1 = 0, a_2 = 0, a_3 = 0;
	double b_1 = 0, b_2 = 0, b_3 = 0;
	double c_1 = 0, c_2 = 0, c_3 = 0;

	const double x_0 = 0;	//mm
	const double y_0 = 0;	//mm
	const double f = 153.24; //焦距mm

	// 像点坐标和地面控制点坐标
	double x = 0, y = 0, X = 0, Y = 0, Z = 0;
	double H = 0;

	// 赋初值 角元素取0，线元素X_s,Y_s取地面控制点,Z_s取7000
	for (int i = 0; i < n; i++)
	{
		matOuterOrientElem[0][0] += param[2 + 5 * i];
		matOuterOrientElem[0][1] += param[3 + 5 * i];
	}
	X_s = matOuterOrientElem[0][0] /= n;	//X_s
	Y_s = matOuterOrientElem[0][1] /= n;	//Y_s
	Z_s = matOuterOrientElem[0][2] = 7000;	//Z_s
	phi = matOuterOrientElem[0][3] = 0;		//phi
	omega = matOuterOrientElem[0][4] = 0;	//omega
	kappa = matOuterOrientElem[0][5] = 0;	//kappa

	a_1 = matR[0][0] = cos(phi) * cos(kappa) - sin(phi) * sin(omega) * sin(kappa);
	a_2 = matR[0][1] = -cos(phi) * sin(kappa) - sin(phi) * sin(omega) * cos(kappa);
	a_3 = matR[0][2] = -sin(phi) * cos(omega);
	b_1 = matR[1][0] = cos(omega) * sin(kappa);
	b_2 = matR[1][1] = cos(omega) * cos(kappa);
	b_3 = matR[1][2] = -sin(omega);
	c_1 = matR[2][0] = sin(phi)* cos(kappa) + cos(phi) * sin(omega) * sin(kappa);
	c_2 = matR[2][1] = -sin(phi) * sin(kappa) + cos(phi) * sin(omega) * cos(kappa);
	c_3 = matR[2][2] = cos(phi) * cos(omega);

	for (int i = 0; i < n; i++)
	{ 
		// 常数项矩阵初值，这里的单位都是mm
		x = param[i * 5];
		y = param[i * 5 + 1];
		X = param[i * 5 + 2];
		Y = param[i * 5 + 3];
		Z = param[i * 5 + 4];
		matL[i * 2][0] = x
			- (-f * (a_1 * (X - X_s) + b_1 * (Y - Y_s) + c_1 * (Z - Z_s))
				/ (a_3 * (X - X_s) + b_3 * (Y - Y_s) + c_3 * (Z - Z_s)) + x_0);
		matL[i * 2 + 1][0] = y
			- (-f * (a_2 * (X - X_s) + b_2 * (Y - Y_s) + c_2 * (Z - Z_s))
				/ (a_3 * (X - X_s) + b_3 * (Y - Y_s) + c_3 * (Z - Z_s)) + y_0);
		// 系数项矩阵初值(垂直摄影使用简化的偏导数)
		H = Z_s - Z;
		matA[i * 2][0] = -f / H * cos(kappa);
		matA[i * 2][1] = -f / H * sin(kappa);
		matA[i * 2][2] = -(x - x_0) / H;
		matA[i * 2][3] = -(f + pow(x - x_0, 2) / f) * cos(kappa) + (x - x_0) * (y - y_0) / f * sin(kappa);
		matA[i * 2][4] = -(x - x_0) * (y - y_0) / f * cos(kappa) - (f + pow(x - x_0, 2) / f) * sin(kappa);
		matA[i * 2][5] = +(y - y_0);
		matA[1 + i * 2][0] = +f / H * sin(kappa);
		matA[1 + i * 2][1] = -f / H * cos(kappa);
		matA[1 + i * 2][2] = -(y - y_0) / H;
		matA[1 + i * 2][3] = -(x - x_0) * (y - y_0) / f * cos(kappa) + (f + pow(y - y_0, 2) / f) * sin(kappa);
		matA[1 + i * 2][4] = -(f + pow(y - y_0, 2) / f) * cos(kappa) - (x - x_0) * (y - y_0) / f * sin(kappa);
		matA[1 + i * 2][5] = -(x - x_0);
	}
	
	// 迭代
	FILE* fout = fopen("resection_result.txt", "w");
	bool flag = false;
	for (int i = 0; i < ITER_TIMES; i++)
	{
		matX = (matA.transpose() * matA).inverse() * (matA.transpose() * matL);
		matOuterOrientElem += matX;

		if (fabs(matX[3][0]) < 3e-5 && fabs(matX[4][0]) < 3e-5 && fabs(matX[5][0]) < 3e-5)
		{
			flag = true;

			std::cout << "完成！迭代次数：" << i+1 << std::endl << "外方位元素：" << matOuterOrientElem << "旋转矩阵：" << matR;
			fprintf(fout, "外方位元素 X_s,Y_s,Z_s,phi,omega,kappa\n %f\t%f\t%f\t%f\t%f\t%f ",
				matOuterOrientElem[0][0], matOuterOrientElem[1][0], matOuterOrientElem[2][0],
				matOuterOrientElem[3][0], matOuterOrientElem[4][0], matOuterOrientElem[5][0]);

			// 计算残差
			Matrix matEpsilon = matA * matX - matL;
			Matrix matEpsilonX = Matrix::zeros(n, 1);
			Matrix matEpsilonY = Matrix::zeros(n, 1);
			fprintf(fout, "残差(mm)：\nid    V\n");
			for (int i = 0; i < n; i++)
			{
				matEpsilonX[i][0] = matEpsilon[i * 2][0];
				matEpsilonY[i][0] = matEpsilon[i * 2 + 1][0];
				fprintf(fout, "x%d\t    %f\n", i + 1, matEpsilon[i * 2][0]);
				fprintf(fout, "y%d\t    %f\n", i + 1, matEpsilon[i * 2 + 1][0]);
			}

			// 中误差
			printf("rms x(mm)：\n");
			printf("%f\n", sqrt((matEpsilonX.transpose() * matEpsilonX)[0][0] / n));
			printf("rms y(mm)：\n");
			printf("%f\n", sqrt((matEpsilonY.transpose() * matEpsilonY)[0][0] / n));

			fprintf(fout, "rms x(mm)：\n");
			fprintf(fout, "%f\n", sqrt((matEpsilonX.transpose() * matEpsilonX)[0][0] / n));
			fprintf(fout, "rms y(mm)：\n");
			fprintf(fout, "%f\n", sqrt((matEpsilonY.transpose() * matEpsilonY)[0][0] / n));

			break;
		}
		// 新值
		X_s = matOuterOrientElem[0][0];		//X_s
		Y_s = matOuterOrientElem[0][1];		//Y_s
		Z_s = matOuterOrientElem[0][2];		//Z_s
		phi = matOuterOrientElem[0][3];		//phi
		omega = matOuterOrientElem[0][4];	//omega
		kappa = matOuterOrientElem[0][5];	//kappa

		a_1 = matR[0][0] = cos(phi) * cos(kappa) - sin(phi) * sin(omega) * sin(kappa);
		a_2 = matR[0][1] = -cos(phi) * sin(kappa) - sin(phi) * sin(omega) * cos(kappa);
		a_3 = matR[0][2] = -sin(phi) * cos(omega);
		b_1 = matR[1][0] = cos(omega) * sin(kappa);
		b_2 = matR[1][1] = cos(omega) * cos(kappa);
		b_3 = matR[1][2] = -sin(omega);
		c_1 = matR[2][0] = sin(phi) * cos(kappa) + cos(phi) * sin(omega) * sin(kappa);
		c_2 = matR[2][1] = -sin(phi) * sin(kappa) + cos(phi) * sin(omega) * cos(kappa);
		c_3 = matR[2][2] = cos(phi) * cos(omega);

		for (int i = 0; i < n; i++)
		{
			// 常数项矩阵初值，这里的单位都是mm
			x = param[i * 5];
			y = param[i * 5 + 1];
			X = param[i * 5 + 2];
			Y = param[i * 5 + 3];
			Z = param[i * 5 + 4];
			matL[i * 2][0] = x
				- (-f * (a_1 * (X - X_s) + b_1 * (Y - Y_s) + c_1 * (Z - Z_s))
					/ (a_3 * (X - X_s) + b_3 * (Y - Y_s) + c_3 * (Z - Z_s)) + x_0);
			matL[i * 2 + 1][0] = y
				- (-f * (a_2 * (X - X_s) + b_2 * (Y - Y_s) + c_2 * (Z - Z_s))
					/ (a_3 * (X - X_s) + b_3 * (Y - Y_s) + c_3 * (Z - Z_s)) + y_0);
			// 系数项矩阵初值(垂直摄影使用简化的偏导数)
			H = Z_s - Z;
			matA[i * 2][0] = -f / H * cos(kappa);
			matA[i * 2][1] = -f / H * sin(kappa);
			matA[i * 2][2] = -(x - x_0) / H;
			matA[i * 2][3] = -(f + pow(x - x_0, 2) / f) * cos(kappa) + (x - x_0) * (y - y_0) / f * sin(kappa);
			matA[i * 2][4] = -(x - x_0) * (y - y_0) / f * cos(kappa) - (f + pow(x - x_0, 2) / f) * sin(kappa);
			matA[i * 2][5] = +(y - y_0);
			matA[1 + i * 2][0] = +f / H * sin(kappa);
			matA[1 + i * 2][1] = -f / H * cos(kappa);
			matA[1 + i * 2][2] = -(y - y_0) / H;
			matA[1 + i * 2][3] = -(x - x_0) * (y - y_0) / f * cos(kappa) + (f + pow(y - y_0, 2) / f) * sin(kappa);
			matA[1 + i * 2][4] = -(f + pow(y - y_0, 2) / f) * cos(kappa) - (x - x_0) * (y - y_0) / f * sin(kappa);
			matA[1 + i * 2][5] = -(x - x_0);
		}

	}
	delete[] param;

	if (!flag)
	{
		std::cout << "结束\n";
		std::cout << "迭代次数达到最大值";
		fprintf(fout, "迭代次数达到最大值!");
	}
	fclose(fout);


}

