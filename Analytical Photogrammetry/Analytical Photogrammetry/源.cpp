#define _CRT_SECURE_NO_WARNINGS 1

#include "stdio.h"
#include "Matrix.h"
#include "math.h"

#ifndef PI
#define PI 3.14159265358979323846
#endif // !PI

#define ITER_TIMES 10

// 后方交会编程
int resection();
// 前方交会编程
int forwardIntersection();
// 内定向编程
int innerOrientation();
// 相对定向编程
int relativeOrientation();
// 绝对定向编程
int absoluteOrientation();

// 辅助函数
double dms2rad(double dms);

int main()
{
	printf("#############################################\n");
	printf("编程一：后方交会\n");
	resection();
	printf("#############################################\n");
	printf("编程二：前方交会\n");
	forwardIntersection();
	printf("#############################################\n");
	printf("编程三：内定向\n");
	innerOrientation();
	printf("#############################################\n");
	printf("编程四：相对定向\n");
	relativeOrientation();
	printf("#############################################\n");
	printf("编程五：绝对定向\n");
	absoluteOrientation();
	printf("#############################################\n");
}

int resection()
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

			std::cout << "完成！迭代次数：" << i + 1 << std::endl << "外方位元素：" << matOuterOrientElem << "旋转矩阵：" << matR;
			fprintf(fout, "外方位元素 X_s,Y_s,Z_s,phi,omega,kappa\n %f\t%f\t%f\t%f\t%f\t%f\n ",
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
	return 0;
}
int forwardIntersection()
{
	// 读取数据(像点坐标只保留了点8031901的数据,故直接使用)
	// 0319上的像点
	FILE* fp = fopen("forwardIntersection_0319.dat", "r");
	if (fp == NULL) {
		perror("打开文件时发生错误");
		return -1;
	}
	char* data = new char[512];

	int n = 0;	// 地面点个数,实际上是1
	char id[50] = { 0 };	// 像点ID，应该是8031901

	fgets(data, 512, fp);
	sscanf(data, "%d", &n);

	double imgCoor0319[2];

	for (int i = 0; fgets(data, 512, fp) && i < n; i++)
	{
		sscanf(data, "%s %*lf %*lf %lf %lf %*d "
			, id, imgCoor0319, imgCoor0319 + 1);
	}
	fclose(fp);

	// 0320上的像点
	fp = fopen("forwardIntersection_0320.dat", "r");
	if (fp == NULL) {
		perror("打开文件时发生错误");
		return -1;
	}

	double imgCoor0320[2];

	fgets(data, 512, fp);
	sscanf(data, "%d", &n);

	fgets(data, 512, fp);
	sscanf(data, "%*s %*lf %*lf %lf %lf %*d ", imgCoor0320, imgCoor0320 + 1);
	fclose(fp);

	// 相机外方位元素
	fp = fopen("forwardIntersection_cmrOuterElem.dat", "r");
	if (fp == NULL) {
		perror("打开文件时发生错误");
		return -1;
	}

	// Phi Omega Kappa Xs Ys Zs
	double cmrElem0319[6] = { 0 };
	double cmrElem0320[6] = { 0 };


	fgets(data, 512, fp);	// 跳过信息头

	fgets(data, 512, fp);
	sscanf(data, "%*d %*s %lf %lf %lf %lf %lf %lf "
		, cmrElem0320, cmrElem0320 + 1, cmrElem0320 + 2, cmrElem0320 + 3, cmrElem0320 + 4, cmrElem0320 + 5);
	fgets(data, 512, fp);
	sscanf(data, "%*d %*s %lf %lf %lf %lf %lf %lf "
		, cmrElem0319, cmrElem0319 + 1, cmrElem0319 + 2, cmrElem0319 + 3, cmrElem0319 + 4, cmrElem0319 + 5);
	fclose(fp);

	// 相机内方位元素
	fp = fopen("forwardIntersection_cmrInnerElem.dat", "r");
	if (fp == NULL) {
		perror("打开文件时发生错误");
		return -1;
	}

	// 内方位元素各相片相同 fk  x0  x0
	double cmrInnerElem[3] = { 0 };

	fgets(data, 512, fp);	// 跳过信息头

	fgets(data, 512, fp);
	sscanf(data, "%lf %lf %lf "
		, cmrInnerElem, cmrInnerElem + 1, cmrInnerElem + 2);
	fclose(fp);
	delete[]data;

	// -------- 方法一：点投影系数（为共线方程法提供初值）--------
	// 定义变量
	// 地面坐标
	double X_A = 0, Y_A = 0, Z_A = 0;
	// 内方位元素
	double x_0 = 0, y_0 = 0, f = 0;
	// 左像片（0320）
	double X_s1 = 0, Y_s1 = 0, Z_s1 = 0;	// 外方位线元素
	double phi_1 = 0, omega_1 = 0, kappa_1 = 0;	// 外方位角元素
	double x_1 = 0, y_1 = 0;	// 像点坐标
	double X_1 = 0, Y_1 = 0, Z_1 = 0;	// 像空间辅助坐标
	// 右像片（0319）
	double X_s2 = 0, Y_s2 = 0, Z_s2 = 0;	// 外方位线元素
	double phi_2 = 0, omega_2 = 0, kappa_2 = 0;	// 外方位角元素
	double x_2 = 0, y_2 = 0;	// 像点坐标
	double X_2 = 0, Y_2 = 0, Z_2 = 0;	// 像空间辅助坐标

	Matrix matR_1 = Matrix::zeros(3, 3);
	Matrix matR_2 = Matrix::zeros(3, 3);

	// 赋值
	f = cmrInnerElem[0], x_0 = cmrInnerElem[1], y_0 = cmrInnerElem[2];
	// 在这里,交换了phi和omega，以及X_S和Y_S
	phi_1 = dms2rad(cmrElem0320[1]);
	omega_1 = dms2rad(cmrElem0320[0]);
	kappa_1 = dms2rad(cmrElem0320[2]);
	X_s1 = cmrElem0320[4];
	Y_s1 = cmrElem0320[3];
	Z_s1 = cmrElem0320[5];
	x_1 = imgCoor0320[0];
	y_1 = imgCoor0320[1];
	phi_2 = dms2rad(cmrElem0319[1]);
	omega_2 = dms2rad(cmrElem0319[0]);
	kappa_2 = dms2rad(cmrElem0319[2]);
	X_s2 = cmrElem0319[4];
	Y_s2 = cmrElem0319[3];
	Z_s2 = cmrElem0319[5];
	x_2 = imgCoor0319[0];
	y_2 = imgCoor0319[1];

	matR_1[0][0] = cos(phi_1) * cos(kappa_1) - sin(phi_1) * sin(omega_1) * sin(kappa_1);
	matR_1[0][1] = -cos(phi_1) * sin(kappa_1) - sin(phi_1) * sin(omega_1) * cos(kappa_1);
	matR_1[0][2] = -sin(phi_1) * cos(omega_1);
	matR_1[1][0] = cos(omega_1) * sin(kappa_1);
	matR_1[1][1] = cos(omega_1) * cos(kappa_1);
	matR_1[1][2] = -sin(omega_1);
	matR_1[2][0] = sin(phi_1) * cos(kappa_1) + cos(phi_1) * sin(omega_1) * sin(kappa_1);
	matR_1[2][1] = -sin(phi_1) * sin(kappa_1) + cos(phi_1) * sin(omega_1) * cos(kappa_1);
	matR_1[2][2] = cos(phi_1) * cos(omega_1);

	matR_2[0][0] = cos(phi_2) * cos(kappa_2) - sin(phi_2) * sin(omega_2) * sin(kappa_2);
	matR_2[0][1] = -cos(phi_2) * sin(kappa_2) - sin(phi_2) * sin(omega_2) * cos(kappa_2);
	matR_2[0][2] = -sin(phi_2) * cos(omega_2);
	matR_2[1][0] = cos(omega_2) * sin(kappa_2);
	matR_2[1][1] = cos(omega_2) * cos(kappa_2);
	matR_2[1][2] = -sin(omega_2);
	matR_2[2][0] = sin(phi_2) * cos(kappa_2) + cos(phi_2) * sin(omega_2) * sin(kappa_2);
	matR_2[2][1] = -sin(phi_2) * sin(kappa_2) + cos(phi_2) * sin(omega_2) * cos(kappa_2);
	matR_2[2][2] = cos(phi_2) * cos(omega_2);

	// 计算其他参数
	double B_X = X_s2 - X_s1;
	double B_Y = Y_s2 - Y_s1;
	double B_Z = Z_s2 - Z_s1;

	Matrix matX_1Y_1Z_1 = matR_1 * (Matrix({ x_1 - x_0,y_1 - y_0,-f }).transpose());
	Matrix matX_2Y_2Z_2 = matR_2 * (Matrix({ x_2 - x_0,y_2 - y_0,-f }).transpose());
	X_1 = matX_1Y_1Z_1[0][0];
	Y_1 = matX_1Y_1Z_1[1][0];
	Z_1 = matX_1Y_1Z_1[2][0];
	X_2 = matX_2Y_2Z_2[0][0];
	Y_2 = matX_2Y_2Z_2[1][0];
	Z_2 = matX_2Y_2Z_2[2][0];

	// 计算点投影系数
	double N_1 = (B_X * Z_2 - B_Z * X_2) / (X_1 * Z_2 - X_2 * Z_1);
	double N_2 = (B_X * Z_1 - B_Z * X_1) / (X_1 * Z_2 - X_2 * Z_1);

	// 计算地面坐标
	X_A = X_s1 + N_1 * X_1;
	Y_A = 0.5 * (Y_s1 + N_1 * Y_1 + Y_s2 + N_2 * Y_2);
	Z_A = Z_s1 + N_1 * Z_1;


	FILE* fout = fopen("forwardIntersection_result.txt", "w");
	// 此处，对X_A和Y_A进行了交换
	std::cout << "+++++++++++++++++++++++++++++++++++++++++++++\n";
	std::cout << "点投影系数:" << N_1 << "\t" << N_2 << std::endl;
	std::cout << "点投影系数法计算的" << id << "地面坐标为：" << std::endl;
	std::cout << "(" << Y_A << "," << X_A << "," << Z_A << ")" << std::endl;

	fprintf(fout, "##############################################\n");
	fprintf(fout, "方法一：点投影系数法计算%s地面坐标\n", id);
	fprintf(fout, "(%f , %f , %f)\n\n", Y_A, X_A, Z_A);



	// -------- 方法二：共线方程严密解法 --------
	std::cout << "+++++++++++++++++++++++++++++++++++++++++++++\n";
	std::cout << "共线方程法\n";
	fprintf(fout, "##############################################\n");
	fprintf(fout, "方法二：共线方程严密解法计算%s地面坐标\n", id);
	int num = 2;	// 像点个数，为2


	// 地面坐标矩阵
	Matrix matGC = Matrix::zeros(3, 1);
	// 解法方程所得的地面坐标
	Matrix matX = Matrix::zeros(3, 1);
	// 误差方程系数矩阵
	Matrix matA = Matrix::zeros(2 * 2, 3);
	// 误差程常数矩阵
	Matrix matL = Matrix::zeros(2 * 2, 1);
	// 旋转矩阵
	//	a_1	a_2	a_3	
	//	b_1	b_2	b_3
	//	c_1	c_2	c_3

	double phi = 0;
	double omega = 0;
	double kappa = 0;
	double X_s = 0;
	double Y_s = 0;
	double Z_s = 0;
	double a_1 = 0, a_2 = 0, a_3 = 0;
	double b_1 = 0, b_2 = 0, b_3 = 0;
	double c_1 = 0, c_2 = 0, c_3 = 0;


	// 像点坐标
	double x = 0, y = 0;
	double H = 0;

	double Xbar;
	double Ybar;
	double Zbar;

	// 赋初值 地面点坐标取点投影系数法结果
	matGC[0][0] = X_A;
	matGC[1][0] = Y_A;
	matGC[2][0] = Z_A;

	X_A = matGC[0][0];
	Y_A = matGC[1][0];
	Z_A = matGC[2][0];
	// 左影像相关参数
	x = x_1;
	y = y_1;
	phi = phi_1;		//在方法一中已经交换过了
	omega = omega_1;	//在方法一中已经交换过了
	kappa = kappa_1;
	X_s = X_s1;			//在方法一中已经交换过了
	Y_s = Y_s1;			//在方法一中已经交换过了
	Z_s = Z_s1;
	a_1 = matR_1[0][0];
	a_2 = matR_1[0][1];
	a_3 = matR_1[0][2];
	b_1 = matR_1[1][0];
	b_2 = matR_1[1][1];
	b_3 = matR_1[1][2];
	c_1 = matR_1[2][0];
	c_2 = matR_1[2][1];
	c_3 = matR_1[2][2];
	H = Z_s - Z_A;
	// x-x^0
	Xbar = a_1 * (X_A - X_s) + b_1 * (Y_A - Y_s) + c_1 * (Z_A - Z_s);
	Ybar = a_2 * (X_A - X_s) + b_2 * (Y_A - Y_s) + c_2 * (Z_A - Z_s);
	Zbar = a_3 * (X_A - X_s) + b_3 * (Y_A - Y_s) + c_3 * (Z_A - Z_s);

	matL[0][0] = x - (-f * Xbar / Zbar + x_0);
	matL[1][0] = y - (-f * Ybar / Zbar + y_0);

	matA[0][0] = +f / H * cos(kappa);
	matA[0][1] = +f / H * sin(kappa);
	matA[0][2] = +x / H;
	matA[1][0] = -f / H * sin(kappa);
	matA[1][1] = +f / H * cos(kappa);
	matA[1][2] = +y / H;


	// 右影像相关参数
	x = x_2;
	y = y_2;
	phi = phi_2;		//在方法一中已经交换过了
	omega = omega_2;	//在方法一中已经交换过了
	kappa = kappa_2;
	X_s = X_s2;			//在方法一中已经交换过了
	Y_s = Y_s2;			//在方法一中已经交换过了
	Z_s = Z_s2;
	a_1 = matR_2[0][0];
	a_2 = matR_2[0][1];
	a_3 = matR_2[0][2];
	b_1 = matR_2[1][0];
	b_2 = matR_2[1][1];
	b_3 = matR_2[1][2];
	c_1 = matR_2[2][0];
	c_2 = matR_2[2][1];
	c_3 = matR_2[2][2];
	H = Z_s - Z_A;

	// x-x^0
	Xbar = a_1 * (X_A - X_s) + b_1 * (Y_A - Y_s) + c_1 * (Z_A - Z_s);
	Ybar = a_2 * (X_A - X_s) + b_2 * (Y_A - Y_s) + c_2 * (Z_A - Z_s);
	Zbar = a_3 * (X_A - X_s) + b_3 * (Y_A - Y_s) + c_3 * (Z_A - Z_s);

	matL[2][0] = x - (-f * Xbar / Zbar + x_0);
	matL[3][0] = y - (-f * Ybar / Zbar + y_0);

	matA[2][0] = +f / H * cos(kappa);
	matA[2][1] = +f / H * sin(kappa);
	matA[2][2] = +x / H;
	matA[3][0] = -f / H * sin(kappa);
	matA[3][1] = +f / H * cos(kappa);
	matA[3][2] = +y / H;

	// 迭代
	bool flag = false;
	Matrix row12Exchange = { {0,1,0} ,{1,0,0}, {0,0,1} };
	for (int i = 0; i < ITER_TIMES; i++)
	{
		matX = (matA.transpose() * matA).inverse() * (matA.transpose() * matL);
		matGC += matX;

		if (fabs(matX[0][0]) < 1e-6 && fabs(matX[1][0]) < 1e-6 && fabs(matX[2][0]) < 1e-6)
		{
			flag = true;
			std::cout << "完成！迭代次数：" << i + 1 << std::endl << "地面点" << id << "坐标：" << row12Exchange * matGC;
			fprintf(fout, "(%f , %f , %f)\n", matGC[1][0], matGC[0][0], matGC[2][0]);

			// 计算残差
			Matrix matEpsilon = matA * matX - matL;
			Matrix matEpsilonX = Matrix::zeros(2, 1);
			Matrix matEpsilonY = Matrix::zeros(2, 1);
			fprintf(fout, "残差(mm)：\nid    V\n");
			for (int i = 0; i < 2; i++)
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
		X_A = matGC[0][0];
		Y_A = matGC[1][0];
		Z_A = matGC[2][0];
		// 左影像相关参数
		x = x_1;
		y = y_1;
		phi = phi_1;		//在方法一中已经交换过了
		omega = omega_1;	//在方法一中已经交换过了
		kappa = kappa_1;
		X_s = X_s1;			//在方法一中已经交换过了
		Y_s = Y_s1;			//在方法一中已经交换过了
		Z_s = Z_s1;
		a_1 = matR_1[0][0];
		a_2 = matR_1[0][1];
		a_3 = matR_1[0][2];
		b_1 = matR_1[1][0];
		b_2 = matR_1[1][1];
		b_3 = matR_1[1][2];
		c_1 = matR_1[2][0];
		c_2 = matR_1[2][1];
		c_3 = matR_1[2][2];
		H = Z_s - Z_A;
		// x-x^0
		Xbar = a_1 * (X_A - X_s) + b_1 * (Y_A - Y_s) + c_1 * (Z_A - Z_s);
		Ybar = a_2 * (X_A - X_s) + b_2 * (Y_A - Y_s) + c_2 * (Z_A - Z_s);
		Zbar = a_3 * (X_A - X_s) + b_3 * (Y_A - Y_s) + c_3 * (Z_A - Z_s);

		matL[0][0] = x - (-f * Xbar / Zbar + x_0);
		matL[1][0] = y - (-f * Ybar / Zbar + y_0);

		matA[0][0] = +f / H * cos(kappa);
		matA[0][1] = +f / H * sin(kappa);
		matA[0][2] = +x / H;
		matA[1][0] = -f / H * sin(kappa);
		matA[1][1] = +f / H * cos(kappa);
		matA[1][2] = +y / H;


		// 右影像相关参数
		x = x_2;
		y = y_2;
		phi = phi_2;		//在方法一中已经交换过了
		omega = omega_2;	//在方法一中已经交换过了
		kappa = kappa_2;
		X_s = X_s2;			//在方法一中已经交换过了
		Y_s = Y_s2;			//在方法一中已经交换过了
		Z_s = Z_s2;
		a_1 = matR_2[0][0];
		a_2 = matR_2[0][1];
		a_3 = matR_2[0][2];
		b_1 = matR_2[1][0];
		b_2 = matR_2[1][1];
		b_3 = matR_2[1][2];
		c_1 = matR_2[2][0];
		c_2 = matR_2[2][1];
		c_3 = matR_2[2][2];
		H = Z_s - Z_A;

		// x-x^0
		Xbar = a_1 * (X_A - X_s) + b_1 * (Y_A - Y_s) + c_1 * (Z_A - Z_s);
		Ybar = a_2 * (X_A - X_s) + b_2 * (Y_A - Y_s) + c_2 * (Z_A - Z_s);
		Zbar = a_3 * (X_A - X_s) + b_3 * (Y_A - Y_s) + c_3 * (Z_A - Z_s);

		matL[2][0] = x - (-f * Xbar / Zbar + x_0);
		matL[3][0] = y - (-f * Ybar / Zbar + y_0);

		matA[2][0] = +f / H * cos(kappa);
		matA[2][1] = +f / H * sin(kappa);
		matA[2][2] = +x / H;
		matA[3][0] = -f / H * sin(kappa);
		matA[3][1] = +f / H * cos(kappa);
		matA[3][2] = +y / H;

	}

	if (!flag)
	{
		std::cout << "结束\n";
		std::cout << "迭代次数达到最大值";
		fprintf(fout, "迭代次数达到最大值!");
	}
	fclose(fout);
	return 0;
}
int innerOrientation()
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
	std::cout << "+++++++++++++++++++++++++++++++++++++++++++++\n";
	std::cout << "内定向正算参数：" << matX;
	fprintf(fout, "+++++++++++++++++++++++++++++++++++++++++++++\n");
	fprintf(fout, "内定向正算参数：\na0\t a1\t a2\t b0\t b1\t b2\n%7e\t%7e\t%7e\t%7e\t%7e\t%7e\n",
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

	std::cout << "+++++++++++++++++++++++++++++++++++++++++++++\n";
	std::cout << "内定向反算参数：" << matX;
	fprintf(fout, "+++++++++++++++++++++++++++++++++++++++++++++\n");
	fprintf(fout, "内定向反算参数：\na0'\t a1'\t a2'\t b0'\t b1'\t b2'\n%7e\t%7e\t%7e\t%7e\t%7e\t%7e\n",
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
	return 0;
}
int relativeOrientation()
{
	// 内方位元素mm
	const double f = 153.8400;
	const double x_0 = 0.0110;
	const double y_0 = 0.0020;

	// 读取像点坐标(文件中的数据已经按顺序整理好)
	// 320 是左像片（1） 319是右像片（2）
	FILE* fp = fopen("relativeOrientation_0320.dat", "r");
	if (fp == NULL) {
		perror("打开文件时发生错误");
		return -1;
	}
	char* data = new char[512];

	int n = 0;	// 同名像点个数

	fgets(data, 512, fp);
	sscanf(data, "%d", &n);

	Matrix imgCoor = Matrix::zeros(n, 5); //id x_1 y_1 x_2 y_2

	for (int i = 0; fgets(data, 512, fp) && i < n; i++)
	{
		sscanf(data, "%lf %*lf %*lf %lf %lf %*d "
			, imgCoor[i] + 0, imgCoor[i] + 1, imgCoor[i] + 2);
	}
	fclose(fp);

	fp = fopen("relativeOrientation_0319.dat", "r");
	if (fp == NULL) {
		perror("打开文件时发生错误");
		return -1;
	}
	fgets(data, 512, fp);

	for (int i = 0; fgets(data, 512, fp) && i < n; i++)
	{
		sscanf(data, "%*lf %*lf %*lf %lf %lf %*d "
			, imgCoor[i] + 3, imgCoor[i] + 4, imgCoor[i]);
	}
	fclose(fp);
	delete[]data;

	// 方法一：连续法相对定向
	// 计算初值
	// 假定基线长
	double B_x = imgCoor[0][1] - imgCoor[0][3];

	// 设定相对定向参数
	// 对应 B_y    B_z
	double mu = 0, nu = 0, phi = 0, omega = 0, kappa = 0;
	Matrix matR = Matrix::zeros(3, 3);


	// 计算像空间辅助坐标
	Matrix matImgAssistCoor1 = Matrix::zeros(3, 1);
	Matrix matImgAssistCoor2 = Matrix::zeros(3, 1);
	double X_1 = 0, Y_1 = 0, Z_1 = 0;
	double X_2 = 0, Y_2 = 0, Z_2 = 0;
	double N_1 = 1, N_2 = 1;
	Matrix matL = Matrix::zeros(n, 1);	// 这里就是Q
	Matrix matX = Matrix::zeros(5, 1);	// 相对定向元素改正数mu, nu, phi, omega, kappa;
	Matrix matRelateOrientElem = Matrix::zeros(5, 1);
	Matrix matA = Matrix::zeros(n, 5);

	double B_y = mu * B_x;
	double B_z = nu * B_x;
	matR[0][0] = cos(phi) * cos(kappa) - sin(phi) * sin(omega) * sin(kappa);
	matR[0][1] = -cos(phi) * sin(kappa) - sin(phi) * sin(omega) * cos(kappa);
	matR[0][2] = -sin(phi) * cos(omega);
	matR[1][0] = cos(omega) * sin(kappa);
	matR[1][1] = cos(omega) * cos(kappa);
	matR[1][2] = -sin(omega);
	matR[2][0] = sin(phi) * cos(kappa) + cos(phi) * sin(omega) * sin(kappa);
	matR[2][1] = -sin(phi) * sin(kappa) + cos(phi) * sin(omega) * cos(kappa);
	matR[2][2] = cos(phi) * cos(omega);
	for (int i = 0; i < n; i++)
	{
		matImgAssistCoor1 = Matrix({ imgCoor[i][1] - x_0,imgCoor[i][2] - y_0,-f }).transpose();
		matImgAssistCoor2 = matR * Matrix({ imgCoor[i][3] - x_0,imgCoor[i][4] - y_0,-f }).transpose();
		X_1 = matImgAssistCoor1[0][0];
		Y_1 = matImgAssistCoor1[1][0];
		Z_1 = matImgAssistCoor1[2][0];
		X_2 = matImgAssistCoor2[0][0];
		Y_2 = matImgAssistCoor2[1][0];
		Z_2 = matImgAssistCoor2[2][0];

		N_1 = (B_x * Z_2 - B_z * X_2) / (X_1 * Z_2 - X_2 * Z_1);
		N_2 = (B_x * Z_1 - B_z * X_1) / (X_1 * Z_2 - X_2 * Z_1);

		// 系数阵
		matA[i][0] = B_x;
		matA[i][1] = -Y_2 / Z_2 * B_x;
		matA[i][2] = -X_2 * Y_2 / Z_2 * N_2;
		matA[i][3] = -(Z_2 + pow(Y_2, 2) / Z_2) * N_2;
		matA[i][4] = X_2 * N_2;
		// 常数项
		matL[i][0] = N_1 * Y_1 - (N_2 * Y_2 + B_y);
	}





	FILE* fout = fopen("relativeOrientation_result.txt", "w");
	printf("+++++++++++++++++++++++++++++++++++++++++++++\n");
	printf("方法一:连续相对定向\n");
	fprintf(fout, "+++++++++++++++++++++++++++++++++++++++++++++\n");
	fprintf(fout, "方法一:连续相对定向\n");

	// 迭代
	bool flag = false;
	for (int i = 0; i < ITER_TIMES; i++)
	{
		matX = (matA.transpose() * matA).inverse() * (matA.transpose() * matL);

		matRelateOrientElem += matX;

		if (fabs(matX[0][0]) < 3e-5 && fabs(matX[1][0]) < 3e-5 && fabs(matX[2][0]) < 3e-5
			&& fabs(matX[3][0]) < 3e-5 && fabs(matX[4][0]) < 3e-5)
		{
			flag = true;
			std::cout << "完成！迭代次数：" << i + 1 << std::endl << "相对定向元素" << matRelateOrientElem;
			fprintf(fout, "相对定向元素:\nmu    nu    phi    omega    kappa\n%lf    %lf    %lf    %lf    %lf\n",
				matRelateOrientElem[0][0], matRelateOrientElem[1][0], matRelateOrientElem[2][0], matRelateOrientElem[3][0], matRelateOrientElem[4][0]);

			// 计算残差
			Matrix matEpsilon = matA * matX - matL;
			fprintf(fout, "残差(mm)：\nid    V_Q\n");
			for (int i = 0; i < n; i++)
			{
				fprintf(fout, "%d\t    %f\n", (int)imgCoor[i][0], matEpsilon[i][0]);
			}

			// 中误差
			fprintf(fout, "中误差(mm)：\n");
			fprintf(fout, "%f\n", sqrt((matEpsilon.transpose() * matEpsilon)[0][0] / n));

			break;
		}

		// 新值
		mu = matRelateOrientElem[0][0];
		nu = matRelateOrientElem[1][0];
		phi = matRelateOrientElem[2][0];
		omega = matRelateOrientElem[3][0];
		kappa = matRelateOrientElem[4][0];

		// 计算像空间辅助坐标
		B_y = mu * B_x;
		B_z = nu * B_x;
		matR[0][0] = cos(phi) * cos(kappa) - sin(phi) * sin(omega) * sin(kappa);
		matR[0][1] = -cos(phi) * sin(kappa) - sin(phi) * sin(omega) * cos(kappa);
		matR[0][2] = -sin(phi) * cos(omega);
		matR[1][0] = cos(omega) * sin(kappa);
		matR[1][1] = cos(omega) * cos(kappa);
		matR[1][2] = -sin(omega);
		matR[2][0] = sin(phi) * cos(kappa) + cos(phi) * sin(omega) * sin(kappa);
		matR[2][1] = -sin(phi) * sin(kappa) + cos(phi) * sin(omega) * cos(kappa);
		matR[2][2] = cos(phi) * cos(omega);
		for (int i = 0; i < n; i++)
		{
			matImgAssistCoor1 = Matrix({ imgCoor[i][1] - x_0,imgCoor[i][2] - y_0,-f }).transpose();
			matImgAssistCoor2 = matR * Matrix({ imgCoor[i][3] - x_0,imgCoor[i][4] - y_0,-f }).transpose();
			X_1 = matImgAssistCoor1[0][0];
			Y_1 = matImgAssistCoor1[1][0];
			Z_1 = matImgAssistCoor1[2][0];
			X_2 = matImgAssistCoor2[0][0];
			Y_2 = matImgAssistCoor2[1][0];
			Z_2 = matImgAssistCoor2[2][0];

			N_1 = (B_x * Z_2 - B_z * X_2) / (X_1 * Z_2 - X_2 * Z_1);
			N_2 = (B_x * Z_1 - B_z * X_1) / (X_1 * Z_2 - X_2 * Z_1);

			// 系数阵
			matA[i][0] = B_x;
			matA[i][1] = -Y_2 / Z_2 * B_x;
			matA[i][2] = -X_2 * Y_2 / Z_2 * N_2;
			matA[i][3] = -(Z_2 + pow(Y_2, 2) / Z_2) * N_2;
			matA[i][4] = X_2 * N_2;
			// 常数项
			matL[i][0] = N_1 * Y_1 - (N_2 * Y_2 + B_y);
		}

	}

	if (!flag)
	{
		std::cout << "结束\n";
		std::cout << "迭代次数达到最大值";
		fprintf(fout, "迭代次数达到最大值!");
	}


	// 方法二：单独法相对定向
	printf("+++++++++++++++++++++++++++++++++++++++++++++\n");
	printf("方法二:单独法相对定向\n");
	fprintf(fout, "+++++++++++++++++++++++++++++++++++++++++++++\n");
	fprintf(fout, "方法二:单独法相对定向\n");
	// 内方位元素
	f;
	x_0;;
	y_0;
	// 计算初值
	// 假定基线长
	double B = imgCoor[0][1] - imgCoor[0][3];

	// 设定相对定向参数
	double phi_1 = 0, kappa_1 = 0, phi_2 = 0, omega_2 = 0, kappa_2 = 0;
	Matrix matR_1 = Matrix::zeros(3, 3);
	Matrix matR_2 = Matrix::zeros(3, 3);


	// 计算像空间辅助坐标
	matImgAssistCoor1 = Matrix::zeros(3, 1);
	matImgAssistCoor2 = Matrix::zeros(3, 1);
	X_1 = 0, Y_1 = 0, Z_1 = 0;
	X_2 = 0, Y_2 = 0, Z_2 = 0;
	matL = Matrix::zeros(n, 1);	// 这里就是Q
	matX = Matrix::zeros(5, 1);	// 相对定向元素改正数mu, nu, phi, omega, kappa;
	matRelateOrientElem = Matrix::zeros(5, 1);
	matA = Matrix::zeros(n, 5);

	matR_1[0][0] = cos(phi_1) * cos(kappa_1);
	matR_1[0][1] = -cos(phi_1) * sin(kappa_1);
	matR_1[0][2] = -sin(phi_1);
	matR_1[1][0] = sin(kappa_1);
	matR_1[1][1] = cos(kappa_1);
	matR_1[1][2] = 0;
	matR_1[2][0] = sin(phi_1) * cos(kappa_1);
	matR_1[2][1] = -sin(phi_1) * sin(kappa_1);
	matR_1[2][2] = cos(phi_1);

	matR_2[0][0] = cos(phi_2) * cos(kappa_2) - sin(phi_2) * sin(omega_2) * sin(kappa_2);
	matR_2[0][1] = -cos(phi_2) * sin(kappa_2) - sin(phi_2) * sin(omega_2) * cos(kappa_2);
	matR_2[0][2] = -sin(phi_2) * cos(omega_2);
	matR_2[1][0] = cos(omega_2) * sin(kappa_2);
	matR_2[1][1] = cos(omega_2) * cos(kappa_2);
	matR_2[1][2] = -sin(omega_2);
	matR_2[2][0] = sin(phi_2) * cos(kappa_2) + cos(phi_2) * sin(omega_2) * sin(kappa_2);
	matR_2[2][1] = -sin(phi_2) * sin(kappa_2) + cos(phi_2) * sin(omega_2) * cos(kappa_2);
	matR_2[2][2] = cos(phi_2) * cos(omega_2);

	for (int i = 0; i < n; i++)
	{
		matImgAssistCoor1 = matR_1 * Matrix({ imgCoor[i][1] - x_0,imgCoor[i][2] - y_0,-f }).transpose();
		matImgAssistCoor2 = matR_2 * Matrix({ imgCoor[i][3] - x_0,imgCoor[i][4] - y_0,-f }).transpose();
		X_1 = matImgAssistCoor1[0][0];
		Y_1 = matImgAssistCoor1[1][0];
		Z_1 = matImgAssistCoor1[2][0];
		X_2 = matImgAssistCoor2[0][0];
		Y_2 = matImgAssistCoor2[1][0];
		Z_2 = matImgAssistCoor2[2][0];

		// 系数阵
		matA[i][0] = X_1 * Y_2 / Z_1;
		matA[i][1] = -X_1;
		matA[i][2] = -X_2 * Y_1 / Z_1;
		matA[i][3] = -(Z_1 + Y_1 * Y_2 / Z_1);
		matA[i][4] = X_2;
		// 常数项q
		matL[i][0] = -f * Y_1 / Z_1 + f * Y_2 / Z_2;
	}



	// 迭代
	flag = false;
	for (int i = 0; i < ITER_TIMES; i++)
	{
		matX = (matA.transpose() * matA).inverse() * (matA.transpose() * matL);

		matRelateOrientElem += matX;

		if (fabs(matX[0][0]) < 3e-5 && fabs(matX[1][0]) < 3e-5 && fabs(matX[2][0]) < 3e-5
			&& fabs(matX[3][0]) < 3e-5 && fabs(matX[4][0]) < 3e-5)
		{
			flag = true;
			std::cout << "完成！迭代次数：" << i + 1 << std::endl << "相对定向元素" << matRelateOrientElem;
			fprintf(fout, "相对定向元素:\nphi_1    kappa_1    phi_2    omega_2    kappa_2\n%lf    %lf    %lf    %lf    %lf\n",
				matRelateOrientElem[0][0], matRelateOrientElem[1][0], matRelateOrientElem[2][0], matRelateOrientElem[3][0], matRelateOrientElem[4][0]);

			// 计算残差
			Matrix matEpsilon = matA * matX - matL;
			fprintf(fout, "残差(mm)：\nid    V_Q\n");
			for (int i = 0; i < n; i++)
			{
				fprintf(fout, "%d\t    %f\n", (int)imgCoor[i][0], matEpsilon[i][0]);
			}

			// 中误差
			fprintf(fout, "中误差(mm)：\n");
			fprintf(fout, "%f\n", sqrt((matEpsilon.transpose() * matEpsilon)[0][0] / n));

			break;
		}

		// 新值
		phi_1 = matRelateOrientElem[0][0];
		kappa_1 = matRelateOrientElem[1][0];
		phi_2 = matRelateOrientElem[2][0];
		omega_2 = matRelateOrientElem[3][0];
		kappa_2 = matRelateOrientElem[4][0];


		matR_1[0][0] = cos(phi_1) * cos(kappa_1);
		matR_1[0][1] = -cos(phi_1) * sin(kappa_1);
		matR_1[0][2] = -sin(phi_1);
		matR_1[1][0] = sin(kappa_1);
		matR_1[1][1] = cos(kappa_1);
		matR_1[1][2] = 0;
		matR_1[2][0] = sin(phi_1) * cos(kappa_1);
		matR_1[2][1] = -sin(phi_1) * sin(kappa_1);
		matR_1[2][2] = cos(phi_1);

		matR_2[0][0] = cos(phi_2) * cos(kappa_2) - sin(phi_2) * sin(omega_2) * sin(kappa_2);
		matR_2[0][1] = -cos(phi_2) * sin(kappa_2) - sin(phi_2) * sin(omega_2) * cos(kappa_2);
		matR_2[0][2] = -sin(phi_2) * cos(omega_2);
		matR_2[1][0] = cos(omega_2) * sin(kappa_2);
		matR_2[1][1] = cos(omega_2) * cos(kappa_2);
		matR_2[1][2] = -sin(omega_2);
		matR_2[2][0] = sin(phi_2) * cos(kappa_2) + cos(phi_2) * sin(omega_2) * sin(kappa_2);
		matR_2[2][1] = -sin(phi_2) * sin(kappa_2) + cos(phi_2) * sin(omega_2) * cos(kappa_2);
		matR_2[2][2] = cos(phi_2) * cos(omega_2);

		for (int i = 0; i < n; i++)
		{
			matImgAssistCoor1 = matR_1 * Matrix({ imgCoor[i][1] - x_0,imgCoor[i][2] - y_0,-f }).transpose();
			matImgAssistCoor2 = matR_2 * Matrix({ imgCoor[i][3] - x_0,imgCoor[i][4] - y_0,-f }).transpose();
			X_1 = matImgAssistCoor1[0][0];
			Y_1 = matImgAssistCoor1[1][0];
			Z_1 = matImgAssistCoor1[2][0];
			X_2 = matImgAssistCoor2[0][0];
			Y_2 = matImgAssistCoor2[1][0];
			Z_2 = matImgAssistCoor2[2][0];

			// 系数阵
			matA[i][0] = X_1 * Y_2 / Z_1;
			matA[i][1] = -X_1;
			matA[i][2] = -X_2 * Y_1 / Z_1;
			matA[i][3] = -(Z_1 + Y_1 * Y_2 / Z_1);
			matA[i][4] = X_2;
			// 常数项q
			matL[i][0] = -f * Y_1 / Z_1 + f * Y_2 / Z_2;
		}

	}

	if (!flag)
	{
		std::cout << "结束\n";
		std::cout << "迭代次数达到最大值";
		fprintf(fout, "迭代次数达到最大值!");
	}

	Matrix temp = matR.inverse() * matR_1.inverse() - matR_2.inverse();

	fclose(fout);
	return 0;
}
int absoluteOrientation()
{
	// 读取文件
	FILE* fp = fopen("absoluteOrientation_data.txt", "r");
	if (fp == NULL) {
		perror("打开文件时发生错误");
		return -1;
	}
	char* data = new char[512];

	int n = 0;	// 模型点个数

	fgets(data, 512, fp);
	sscanf(data, "%d", &n);

	fgets(data, 512, fp); // 跳过信息头

	Matrix modelCoor = Matrix::zeros(n, 3); //x y z
	Matrix groundCoor = Matrix::zeros(n, 3); //X Y Z

	for (int i = 0; fgets(data, 512, fp) && i < n; i++)
	{
		sscanf(data, "%*s %lf %lf %lf %lf %lf %lf ",
			modelCoor[i] + 0, modelCoor[i] + 1, modelCoor[i] + 2,
			groundCoor[i] + 0, groundCoor[i] + 1, groundCoor[i] + 2);
	}
	fclose(fp);
	delete[] data;

	// 方法一：不进行重心化（由于点数很少，计算可行）
	Matrix matA = Matrix::zeros(3 * n, 7);
	Matrix matL = Matrix::zeros(3 * n, 1);
	Matrix matX = Matrix::zeros(7, 1);	// X0 Y0 Z0 lambda Phi Omega Kappa
	Matrix matSevenParam = Matrix::zeros(7, 1);	// X0 Y0 Z0 lambda Phi Omega Kappa
	Matrix matR = Matrix::zeros(3, 3);

	double Phi = 0, Omega = 0, Kappa = 0;
	double X_0 = 0, Y_0 = 0, Z_0 = 0;
	double lambda = 1;


	// 确定初值
	Phi = 0;
	Omega = 0;
	Kappa = 0;

	Matrix modelCoorG = Matrix::ones(1, n) * (modelCoor) / n;
	Matrix groundCoorG = Matrix::ones(1, n) * (groundCoor) / n;
	lambda = (groundCoor[0][0] - groundCoorG[0][0]) / (modelCoor[0][0] - modelCoorG[0][0]);
	X_0 = (groundCoorG - lambda * modelCoorG)[0][0];
	Y_0 = (groundCoorG - lambda * modelCoorG)[0][1];
	Z_0 = (groundCoorG - lambda * modelCoorG)[0][2];


	matSevenParam[0][0] = X_0;
	matSevenParam[0][1] = Y_0;
	matSevenParam[0][2] = Z_0;
	matSevenParam[0][3] = lambda;
	matSevenParam[0][4] = Phi;
	matSevenParam[0][5] = Omega;
	matSevenParam[0][6] = Kappa;

	matR[0][0] = cos(Phi) * cos(Kappa) - sin(Phi) * sin(Omega) * sin(Kappa);
	matR[0][1] = -cos(Phi) * sin(Kappa) - sin(Phi) * sin(Omega) * cos(Kappa);
	matR[0][2] = -sin(Phi) * cos(Omega);
	matR[1][0] = cos(Omega) * sin(Kappa);
	matR[1][1] = cos(Omega) * cos(Kappa);
	matR[1][2] = -sin(Omega);
	matR[2][0] = sin(Phi) * cos(Kappa) + cos(Phi) * sin(Omega) * sin(Kappa);
	matR[2][1] = -sin(Phi) * sin(Kappa) + cos(Phi) * sin(Omega) * cos(Kappa);
	matR[2][2] = cos(Phi) * cos(Omega);

	// [x' y' z']T = R[xp yp zp]T
	Matrix modelCoor_Prime = modelCoor * matR.transpose();

	for (int i = 0; i < n; i++)
	{
		// 系数阵
		matA[i * 3][0] = 1;

		matA[i * 3 + 1][1] = 1;

		matA[i * 3 + 2][2] = 1;

		matA[i * 3][3] = modelCoor_Prime[i][0];
		matA[i * 3 + 1][3] = modelCoor_Prime[i][1];
		matA[i * 3 + 2][3] = modelCoor_Prime[i][2];

		matA[i * 3][4] = -lambda * modelCoor_Prime[i][2];
		matA[i * 3 + 2][4] = lambda * modelCoor_Prime[i][0];

		matA[i * 3 + 1][5] = -lambda * modelCoor_Prime[i][2];
		matA[i * 3 + 2][5] = lambda * modelCoor_Prime[i][1];

		matA[i * 3][6] = -lambda * modelCoor_Prime[i][1];
		matA[i * 3 + 1][6] = lambda * modelCoor_Prime[i][0];

		// 常数项
		matL[i * 3][0] = groundCoor[i][0] - lambda * modelCoor_Prime[i][0] - X_0;
		matL[i * 3 + 1][0] = groundCoor[i][1] - lambda * modelCoor_Prime[i][1] - Y_0;
		matL[i * 3 + 2][0] = groundCoor[i][2] - lambda * modelCoor_Prime[i][2] - Z_0;

	}

	FILE* fout = fopen("absoluteOrientation_result.txt", "w");
	printf("+++++++++++++++++++++++++++++++++++++++++++++\n");
	printf("方法一:不进行重心化\n");
	fprintf(fout, "+++++++++++++++++++++++++++++++++++++++++++++\n");
	fprintf(fout, "方法一:不进行重心化\n");

	// 迭代
	bool flag = false;
	for (int i = 0; i < ITER_TIMES; i++)
	{
		matX = (matA.transpose() * matA).inverse() * (matA.transpose() * matL);
		matSevenParam += matX;

		if (fabs(matX[4][0]) < 3e-5 && fabs(matX[5][0]) < 3e-5 && fabs(matX[6][0]) < 3e-5)
		{
			flag = true;
			std::cout << "完成！迭代次数：" << i + 1 << std::endl << "七参数" << matSevenParam;
			fprintf(fout, "绝对定向元素:\nX0    Y0    Z0    lambda    Phi    Omega    Kappa\n%lf    %lf    %lf    %lf    %lf    %lf    %lf\n",
				matSevenParam[0][0], matSevenParam[1][0], matSevenParam[2][0], matSevenParam[3][0],
				matSevenParam[4][0], matSevenParam[5][0], matSevenParam[6][0]);

			// 计算残差
			Matrix matEpsilon = matA * matX - matL;
			Matrix matEpsilonX = Matrix::zeros(n, 1);
			Matrix matEpsilonY = Matrix::zeros(n, 1);
			Matrix matEpsilonZ = Matrix::zeros(n, 1);
			fprintf(fout, "残差(m)：\nid    V_Q\n");
			for (int i = 0; i < n; i++)
			{
				matEpsilonX = matEpsilon[3 * i][0];
				matEpsilonY = matEpsilon[3 * i + 1][0];
				matEpsilonZ = matEpsilon[3 * i + 2][0];
				fprintf(fout, "Xtp%d\t    %f\n", i + 1, matEpsilon[3 * i][0]);
				fprintf(fout, "Ytp%d\t    %f\n", i + 1, matEpsilon[3 * i + 1][0]);
				fprintf(fout, "Ztp%d\t    %f\n", i + 1, matEpsilon[3 * i + 2][0]);
			}

			// 中误差
			printf("rms X(m)：\n");
			printf("%f\n", sqrt((matEpsilonX.transpose() * matEpsilonX)[0][0] / n));
			printf("rms Y(m)：\n");
			printf("%f\n", sqrt((matEpsilonY.transpose() * matEpsilonY)[0][0] / n));
			printf("rms Z(m)：\n");
			printf("%f\n", sqrt((matEpsilonY.transpose() * matEpsilonZ)[0][0] / n));

			fprintf(fout, "rms X(m)：\n");
			fprintf(fout, "%f\n", sqrt((matEpsilonX.transpose() * matEpsilonX)[0][0] / n));
			fprintf(fout, "rms Y(m)：\n");
			fprintf(fout, "%f\n", sqrt((matEpsilonY.transpose() * matEpsilonY)[0][0] / n));
			fprintf(fout, "rms Z(m)：\n");
			fprintf(fout, "%f\n", sqrt((matEpsilonY.transpose() * matEpsilonZ)[0][0] / n));
			break;
		}

		// 新值
		X_0 = matSevenParam[0][0];
		Y_0 = matSevenParam[0][1];
		Z_0 = matSevenParam[0][2];
		lambda = matSevenParam[0][3];
		Phi = matSevenParam[0][4];
		Omega = matSevenParam[0][5];
		Kappa = matSevenParam[0][6];

		matR[0][0] = cos(Phi) * cos(Kappa) - sin(Phi) * sin(Omega) * sin(Kappa);
		matR[0][1] = -cos(Phi) * sin(Kappa) - sin(Phi) * sin(Omega) * cos(Kappa);
		matR[0][2] = -sin(Phi) * cos(Omega);
		matR[1][0] = cos(Omega) * sin(Kappa);
		matR[1][1] = cos(Omega) * cos(Kappa);
		matR[1][2] = -sin(Omega);
		matR[2][0] = sin(Phi) * cos(Kappa) + cos(Phi) * sin(Omega) * sin(Kappa);
		matR[2][1] = -sin(Phi) * sin(Kappa) + cos(Phi) * sin(Omega) * cos(Kappa);
		matR[2][2] = cos(Phi) * cos(Omega);

		// [x' y' z']T = R[xp yp zp]T
		modelCoor_Prime = modelCoor * matR.transpose();

		for (int i = 0; i < n; i++)
		{
			// 系数阵
			matA[i * 3][0] = 1;

			matA[i * 3 + 1][1] = 1;

			matA[i * 3 + 2][2] = 1;

			matA[i * 3][3] = modelCoor_Prime[i][0];
			matA[i * 3 + 1][3] = modelCoor_Prime[i][1];
			matA[i * 3 + 2][3] = modelCoor_Prime[i][2];

			matA[i * 3][4] = -lambda * modelCoor_Prime[i][2];
			matA[i * 3 + 2][4] = lambda * modelCoor_Prime[i][0];

			matA[i * 3 + 1][5] = -lambda * modelCoor_Prime[i][2];
			matA[i * 3 + 2][5] = lambda * modelCoor_Prime[i][1];

			matA[i * 3][6] = -lambda * modelCoor_Prime[i][1];
			matA[i * 3 + 1][6] = lambda * modelCoor_Prime[i][0];

			// 常数项
			matL[i * 3][0] = groundCoor[i][0] - lambda * modelCoor_Prime[i][0] - X_0;
			matL[i * 3 + 1][0] = groundCoor[i][1] - lambda * modelCoor_Prime[i][1] - Y_0;
			matL[i * 3 + 2][0] = groundCoor[i][2] - lambda * modelCoor_Prime[i][2] - Z_0;

		}
	}
	if (!flag)
	{
		std::cout << "结束\n";
		std::cout << "迭代次数达到最大值";
		fprintf(fout, "迭代次数达到最大值!");
	}



	// 方法二：重心化 
	printf("+++++++++++++++++++++++++++++++++++++++++++++\n");
	printf("方法二:进行重心化\n");
	fprintf(fout, "+++++++++++++++++++++++++++++++++++++++++++++\n");
	fprintf(fout, "方法二:进行重心化\n");
	// 确定初值
	matA = Matrix::zeros(3 * n, 7);
	matL = Matrix::zeros(3 * n, 1);
	matX = Matrix::zeros(7, 1);	// X0 Y0 Z0 lambda Phi Omega Kappa
	matSevenParam = Matrix::zeros(7, 1);	// X0 Y0 Z0 lambda Phi Omega Kappa
	matR = Matrix::zeros(3, 3);

	Phi = 0;
	Omega = 0;
	Kappa = 0;


	modelCoorG = Matrix::ones(1, n) * (modelCoor) / n;
	groundCoorG = Matrix::ones(1, n) * (groundCoor) / n;
	lambda = (groundCoor[0][0] - groundCoorG[0][0]) / (modelCoor[0][0] - modelCoorG[0][0]);

	Matrix modelCoorBar = modelCoor - Matrix::ones(n, 1) * modelCoorG;
	Matrix groundCoorBar = groundCoor - Matrix::ones(n, 1) * groundCoorG;


	// 不需要迭代计算X0 Y0 Z0
	matSevenParam[0][0] = (groundCoorG - modelCoorG)[0][0];
	matSevenParam[0][1] = (groundCoorG - modelCoorG)[0][1];
	matSevenParam[0][2] = (groundCoorG - modelCoorG)[0][2];

	matSevenParam[0][3] = lambda;
	matSevenParam[0][4] = Phi;
	matSevenParam[0][5] = Omega;
	matSevenParam[0][6] = Kappa;

	matR[0][0] = cos(Phi) * cos(Kappa) - sin(Phi) * sin(Omega) * sin(Kappa);
	matR[0][1] = -cos(Phi) * sin(Kappa) - sin(Phi) * sin(Omega) * cos(Kappa);
	matR[0][2] = -sin(Phi) * cos(Omega);
	matR[1][0] = cos(Omega) * sin(Kappa);
	matR[1][1] = cos(Omega) * cos(Kappa);
	matR[1][2] = -sin(Omega);
	matR[2][0] = sin(Phi) * cos(Kappa) + cos(Phi) * sin(Omega) * sin(Kappa);
	matR[2][1] = -sin(Phi) * sin(Kappa) + cos(Phi) * sin(Omega) * cos(Kappa);
	matR[2][2] = cos(Phi) * cos(Omega);


	Matrix modelCoorBar_Prime = modelCoorBar * matR.transpose();
	Matrix A_Degree = Matrix::zeros(3, 3);
	Matrix L_Degree = Matrix::zeros(3, 1);
	double sigma_x2_y2_z2 = 0;
	double sigma_xl_yl_zl = 0;

	for (int i = 0; i < n; i++)
	{
		// 系数阵
		matA[i * 3][0] = 1;

		matA[i * 3 + 1][1] = 1;

		matA[i * 3 + 2][2] = 1;

		matA[i * 3][3] = modelCoorBar[i][0];
		matA[i * 3 + 1][3] = modelCoorBar[i][1];
		matA[i * 3 + 2][3] = modelCoorBar[i][2];

		matA[i * 3][4] = -lambda * modelCoorBar[i][2];
		matA[i * 3 + 2][4] = lambda * modelCoorBar[i][0];

		matA[i * 3 + 1][5] = -lambda * modelCoorBar[i][2];
		matA[i * 3 + 2][5] = lambda * modelCoorBar[i][1];

		matA[i * 3][6] = -lambda * modelCoorBar[i][1];
		matA[i * 3 + 1][6] = lambda * modelCoorBar[i][0];

		// 常数项		
		matL[i * 3][0] = groundCoorBar[i][0] - lambda * modelCoorBar_Prime[i][0];
		matL[i * 3 + 1][0] = groundCoorBar[i][1] - lambda * modelCoorBar_Prime[i][1];
		matL[i * 3 + 2][0] = groundCoorBar[i][2] - lambda * modelCoorBar_Prime[i][2];


		// 法方程系数阵特定项和部分值
		sigma_x2_y2_z2 += pow(modelCoorBar[i][0], 2) + pow(modelCoorBar[i][1], 2) + pow(modelCoorBar[i][2], 2);
		sigma_xl_yl_zl += modelCoorBar[i][0] * matL[i * 3][0] + modelCoorBar[i][1] * matL[i * 3 + 1][0] + modelCoorBar[i][2] * matL[i * 3 + 2][0];

		A_Degree[0][0] += pow(modelCoorBar[i][0], 2) + pow(modelCoorBar[i][2], 2);
		A_Degree[0][1] += modelCoorBar[i][0] * modelCoorBar[i][1];
		A_Degree[0][2] += modelCoorBar[i][1] * modelCoorBar[i][2];

		A_Degree[1][1] += pow(modelCoorBar[i][1], 2) + pow(modelCoorBar[i][2], 2);
		A_Degree[1][2] += modelCoorBar[i][0] * modelCoorBar[i][2];;

		A_Degree[2][2] += pow(modelCoorBar[i][0], 2) + pow(modelCoorBar[i][1], 2);

		L_Degree[0][0] += modelCoorBar[i][0] * matL[i * 3 + 2][0] - modelCoorBar[i][2] * matL[i * 3][0];
		L_Degree[0][1] += modelCoorBar[i][1] * matL[i * 3 + 2][0] - modelCoorBar[i][2] * matL[i * 3 + 1][0];
		L_Degree[0][2] += modelCoorBar[i][0] * matL[i * 3 + 1][0] - modelCoorBar[i][1] * matL[i * 3][0];
	}
	A_Degree[1][0] = A_Degree[0][1];
	A_Degree[2][0] = A_Degree[0][2];
	A_Degree[2][1] = A_Degree[1][2];
	A_Degree = A_Degree * pow(lambda, 2);
	L_Degree = L_Degree * lambda;

	// 迭代
	flag = false;
	Matrix matTemp;
	Matrix a;
	for (int i = 0; i < ITER_TIMES; i++)
	{
		matX[3][0] = sigma_xl_yl_zl / sigma_x2_y2_z2;
		matTemp = A_Degree.inverse() * L_Degree;
		matX[4][0] = matTemp[0][0];
		matX[5][0] = matTemp[1][0];
		matX[6][0] = matTemp[2][0];

		a = matA.transpose() * matA;

		//matX = (matA.transpose() * matA).inverse() * (matA.transpose() * matL);


		matSevenParam += matX;

		if (fabs(matX[4][0]) < 3e-5 && fabs(matX[5][0]) < 3e-5 && fabs(matX[6][0]) < 3e-5)
		{
			flag = true;
			std::cout << "完成！迭代次数：" << i + 1 << std::endl << "七参数" << matSevenParam;
			printf("重心坐标：\n(Xpg,Ypg,Zpg)\n(%lf,%lf,%lf)\n(Xtpg,Ytpg,Ztpg)\n(%lf,%lf,%lf)\n",
				modelCoorG[0][0], modelCoorG[0][1], modelCoorG[0][2], groundCoorG[0][0], groundCoorG[0][1], groundCoorG[0][2]);
			fprintf(fout, "重心坐标：\n(Xpg,Ypg,Zpg)    (Xtpg,Ytpg,Ztpg)\n(%lf,%lf,%lf)    (%lf,%lf,%lf)\n",
				modelCoorG[0][0], modelCoorG[0][1], modelCoorG[0][2], groundCoorG[0][0], groundCoorG[0][1], groundCoorG[0][2]);
			fprintf(fout, "绝对定向元素:\nX0    Y0    Z0    lambda    Phi    Omega    Kappa\n%lf    %lf    %lf    %lf    %lf    %lf    %lf\n",
				matSevenParam[0][0], matSevenParam[1][0], matSevenParam[2][0], matSevenParam[3][0],
				matSevenParam[4][0], matSevenParam[5][0], matSevenParam[6][0]);

			// 计算残差
			Matrix matEpsilon = matA * matX - matL;
			Matrix matEpsilonX = Matrix::zeros(n, 1);
			Matrix matEpsilonY = Matrix::zeros(n, 1);
			Matrix matEpsilonZ = Matrix::zeros(n, 1);
			fprintf(fout, "残差(m)：\nid    V_Q\n");
			for (int i = 0; i < n; i++)
			{
				matEpsilonX = matEpsilon[3 * i][0];
				matEpsilonY = matEpsilon[3 * i + 1][0];
				matEpsilonZ = matEpsilon[3 * i + 2][0];
				fprintf(fout, "Xtp%d\t    %f\n", i + 1, matEpsilon[3 * i][0]);
				fprintf(fout, "Ytp%d\t    %f\n", i + 1, matEpsilon[3 * i + 1][0]);
				fprintf(fout, "Ztp%d\t    %f\n", i + 1, matEpsilon[3 * i + 2][0]);
			}

			// 中误差
			printf("rms X(m)：\n");
			printf("%f\n", sqrt((matEpsilonX.transpose() * matEpsilonX)[0][0] / n));
			printf("rms Y(m)：\n");
			printf("%f\n", sqrt((matEpsilonY.transpose() * matEpsilonY)[0][0] / n));
			printf("rms Z(m)：\n");
			printf("%f\n", sqrt((matEpsilonY.transpose() * matEpsilonZ)[0][0] / n));

			fprintf(fout, "rms X(m)：\n");
			fprintf(fout, "%f\n", sqrt((matEpsilonX.transpose() * matEpsilonX)[0][0] / n));
			fprintf(fout, "rms Y(m)：\n");
			fprintf(fout, "%f\n", sqrt((matEpsilonY.transpose() * matEpsilonY)[0][0] / n));
			fprintf(fout, "rms Z(m)：\n");
			fprintf(fout, "%f\n", sqrt((matEpsilonY.transpose() * matEpsilonZ)[0][0] / n));
			break;
		}

		// 新值
		X_0;
		Y_0;
		Z_0;
		lambda = matSevenParam[0][3];
		Phi = matSevenParam[0][4];
		Omega = matSevenParam[0][5];
		Kappa = matSevenParam[0][6];


		matR[0][0] = cos(Phi) * cos(Kappa) - sin(Phi) * sin(Omega) * sin(Kappa);
		matR[0][1] = -cos(Phi) * sin(Kappa) - sin(Phi) * sin(Omega) * cos(Kappa);
		matR[0][2] = -sin(Phi) * cos(Omega);
		matR[1][0] = cos(Omega) * sin(Kappa);
		matR[1][1] = cos(Omega) * cos(Kappa);
		matR[1][2] = -sin(Omega);
		matR[2][0] = sin(Phi) * cos(Kappa) + cos(Phi) * sin(Omega) * sin(Kappa);
		matR[2][1] = -sin(Phi) * sin(Kappa) + cos(Phi) * sin(Omega) * cos(Kappa);
		matR[2][2] = cos(Phi) * cos(Omega);


		modelCoorBar_Prime = modelCoorBar * matR.transpose();
		A_Degree = Matrix::zeros(3, 3);
		L_Degree = Matrix::zeros(3, 1);
		sigma_x2_y2_z2 = 0;
		sigma_xl_yl_zl = 0;

		for (int i = 0; i < n; i++)
		{
			// 系数阵
			matA[i * 3][0] = 1;

			matA[i * 3 + 1][1] = 1;

			matA[i * 3 + 2][2] = 1;

			matA[i * 3][3] = modelCoorBar[i][0];
			matA[i * 3 + 1][3] = modelCoorBar[i][1];
			matA[i * 3 + 2][3] = modelCoorBar[i][2];

			matA[i * 3][4] = -lambda * modelCoorBar[i][2];
			matA[i * 3 + 2][4] = lambda * modelCoorBar[i][0];

			matA[i * 3 + 1][5] = -lambda * modelCoorBar[i][2];
			matA[i * 3 + 2][5] = lambda * modelCoorBar[i][1];

			matA[i * 3][6] = -lambda * modelCoorBar[i][1];
			matA[i * 3 + 1][6] = lambda * modelCoorBar[i][0];

			// 常数项		
			matL[i * 3][0] = groundCoorBar[i][0] - lambda * modelCoorBar_Prime[i][0];
			matL[i * 3 + 1][0] = groundCoorBar[i][1] - lambda * modelCoorBar_Prime[i][1];
			matL[i * 3 + 2][0] = groundCoorBar[i][2] - lambda * modelCoorBar_Prime[i][2];


			// 法方程系数阵特定项和部分值
			sigma_x2_y2_z2 += pow(modelCoorBar[i][0], 2) + pow(modelCoorBar[i][1], 2) + pow(modelCoorBar[i][2], 2);
			sigma_xl_yl_zl += modelCoorBar[i][0] * matL[i * 3][0] + modelCoorBar[i][1] * matL[i * 3 + 1][0] + modelCoorBar[i][2] * matL[i * 3 + 2][0];

			A_Degree[0][0] += pow(modelCoorBar[i][0], 2) + pow(modelCoorBar[i][2], 2);
			A_Degree[0][1] += modelCoorBar[i][0] * modelCoorBar[i][1];
			A_Degree[0][2] += modelCoorBar[i][1] * modelCoorBar[i][2];

			A_Degree[1][1] += pow(modelCoorBar[i][1], 2) + pow(modelCoorBar[i][2], 2);
			A_Degree[1][2] += modelCoorBar[i][0] * modelCoorBar[i][2];;

			A_Degree[2][2] += pow(modelCoorBar[i][0], 2) + pow(modelCoorBar[i][1], 2);

			L_Degree[0][0] += modelCoorBar[i][0] * matL[i * 3 + 2][0] - modelCoorBar[i][2] * matL[i * 3][0];
			L_Degree[0][1] += modelCoorBar[i][1] * matL[i * 3 + 2][0] - modelCoorBar[i][2] * matL[i * 3 + 1][0];
			L_Degree[0][2] += modelCoorBar[i][0] * matL[i * 3 + 1][0] - modelCoorBar[i][1] * matL[i * 3][0];
		}
		A_Degree[1][0] = A_Degree[0][1];
		A_Degree[2][0] = A_Degree[0][2];
		A_Degree[2][1] = A_Degree[1][2];
		A_Degree = A_Degree * pow(lambda, 2);
		L_Degree = L_Degree * lambda;
	}
	if (!flag)
	{
		std::cout << "结束\n";
		std::cout << "迭代次数达到最大值";
		fprintf(fout, "迭代次数达到最大值!");
	}

	fclose(fout);
	return 0;
}



















// 将外方位元素的特殊格式（度.分秒）化为弧度
double dms2rad(double dms)
{
	double rad = 0;
	double temp = fabs(dms);
	double d = floor(temp);
	temp -= d;
	temp *= 100;
	double m = floor(temp);
	temp -= m;
	temp *= 100;
	double s = temp;
	rad = (d + m / 60.0 + s / 3600.0) / 180 * PI;
	return dms < 0 ? -rad : rad;
}