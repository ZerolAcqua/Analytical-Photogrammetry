#define _CRT_SECURE_NO_WARNINGS 1

#include "stdio.h"
#include "Matrix.h"
#include "math.h"

#ifndef PI
#define PI 3.14159265358979323846
#endif // !PI

#define ITER_TIMES 10


double dms2rad(double dms);

// 前方交会编程
// 由于WUCAPS中存在左右手坐标系的问题的
// 所以需要把外方位元素phi和omega，X_s,Y_s互换
// 方可得到地面坐标
// 同时还需要把地面坐标的X_A和Y_A交换

int main()
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
	std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++\n";
	std::cout << "点投影系数:"<< N_1 << "\t" << N_2 << std::endl;
	std::cout << "点投影系数法计算的" << id << "地面坐标为：" << std::endl;
	std::cout << "(" << Y_A << "," << X_A << "," << Z_A << ")" << std::endl;

	fprintf(fout, "##############################################\n");
	fprintf(fout, "方法一：点投影系数法计算%s地面坐标\n",id);
	fprintf(fout, "(%f , %f , %f)\n\n", Y_A, X_A, Z_A);



	// -------- 方法二：共线方程严密解法 --------
	std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++\n";
	std::cout << "共线方程法\n";
	fprintf(fout, "##############################################\n");
	fprintf(fout, "方法二：共线方程严密解法计算%s地面坐标\n",id);
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