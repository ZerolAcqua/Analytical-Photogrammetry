#define _CRT_SECURE_NO_WARNINGS 1

#include "stdio.h"
#include "Matrix.h"
#include "math.h"

#ifndef PI
#define PI 3.14159265358979323846
#endif // !PI

#define ITER_TIMES 10

// 相对定向编程

int main()
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
	printf("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
	printf("方法一:连续相对定向\n");
	fprintf(fout, "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
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
				fprintf(fout, "%d\t    %f\n",(int)imgCoor[i][0], matEpsilon[i][0]);
			}
			
			// 中误差
			fprintf(fout, "中误差(mm)：\n");
			fprintf(fout, "%f\n",sqrt((matEpsilon.transpose() * matEpsilon)[0][0] / n));

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
	printf("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
	printf("方法二:单独法相对定向\n");
	fprintf(fout, "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
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
		matA[i][3] = -(Z_1 + Y_1*Y_2 / Z_1);
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
}