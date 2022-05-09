#define _CRT_SECURE_NO_WARNINGS 1

#include "stdio.h"
#include "Matrix.h"
#include "math.h"

#ifndef PI
#define PI 3.14159265358979323846
#endif // !PI

#define ITER_TIMES 10

// 绝对定向编程

int main()
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
	printf("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
	printf("方法一:不进行重心化\n");
	fprintf(fout, "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
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
	printf("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
	printf("方法二:进行重心化\n");
	fprintf(fout, "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
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

		a= matA.transpose() * matA;

		//matX = (matA.transpose() * matA).inverse() * (matA.transpose() * matL);


		matSevenParam += matX;

		if (fabs(matX[4][0]) < 3e-5 && fabs(matX[5][0]) < 3e-5 && fabs(matX[6][0]) < 3e-5)
		{
			flag = true;
			std::cout << "完成！迭代次数：" << i + 1 << std::endl << "七参数" << matSevenParam;
			printf("重心坐标：\n(Xpg,Ypg,Zpg)    (Xtpg,Ytpg,Ztpg)\n(%lf,%lf,%lf)    (%lf,%lf,%lf)\n",
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
			printf("%f\n", sqrt((matEpsilonY.transpose()* matEpsilonZ)[0][0] / n));

			fprintf(fout, "rms X(m)：\n");
			fprintf(fout, "%f\n", sqrt((matEpsilonX.transpose() * matEpsilonX)[0][0] / n));
			fprintf(fout, "rms Y(m)：\n");
			fprintf(fout, "%f\n", sqrt((matEpsilonY.transpose() * matEpsilonY)[0][0] / n));
			fprintf(fout, "rms Z(m)：\n");
			fprintf(fout, "%f\n", sqrt((matEpsilonY.transpose()* matEpsilonZ)[0][0] / n));
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