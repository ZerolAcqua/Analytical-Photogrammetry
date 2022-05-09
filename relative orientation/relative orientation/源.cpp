#define _CRT_SECURE_NO_WARNINGS 1

#include "stdio.h"
#include "Matrix.h"
#include "math.h"

#ifndef PI
#define PI 3.14159265358979323846
#endif // !PI

#define ITER_TIMES 10

// ��Զ�����

int main()
{
	// �ڷ�λԪ��mm
	const double f = 153.8400;
	const double x_0 = 0.0110;
	const double y_0 = 0.0020;

	// ��ȡ�������(�ļ��е������Ѿ���˳�������)
	// 320 ������Ƭ��1�� 319������Ƭ��2��
	FILE* fp = fopen("relativeOrientation_0320.dat", "r");
	if (fp == NULL) {
		perror("���ļ�ʱ��������");
		return -1;
	}
	char* data = new char[512];

	int n = 0;	// ͬ��������

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
		perror("���ļ�ʱ��������");
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

	// ����һ����������Զ���
	// �����ֵ
	// �ٶ����߳�
	double B_x = imgCoor[0][1] - imgCoor[0][3];
	
	// �趨��Զ������
	// ��Ӧ B_y    B_z
	double mu = 0, nu = 0, phi = 0, omega = 0, kappa = 0;
	Matrix matR = Matrix::zeros(3, 3);


	// ������ռ丨������
	Matrix matImgAssistCoor1 = Matrix::zeros(3, 1);
	Matrix matImgAssistCoor2 = Matrix::zeros(3, 1);
	double X_1 = 0, Y_1 = 0, Z_1 = 0;
	double X_2 = 0, Y_2 = 0, Z_2 = 0;
	double N_1 = 1, N_2 = 1;
	Matrix matL = Matrix::zeros(n, 1);	// �������Q
	Matrix matX = Matrix::zeros(5, 1);	// ��Զ���Ԫ�ظ�����mu, nu, phi, omega, kappa;
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

		// ϵ����
		matA[i][0] = B_x;
		matA[i][1] = -Y_2 / Z_2 * B_x;
		matA[i][2] = -X_2 * Y_2 / Z_2 * N_2;
		matA[i][3] = -(Z_2 + pow(Y_2, 2) / Z_2) * N_2;
		matA[i][4] = X_2 * N_2;
		// ������
		matL[i][0] = N_1 * Y_1 - (N_2 * Y_2 + B_y);
	}





	FILE* fout = fopen("relativeOrientation_result.txt", "w");
	printf("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
	printf("����һ:������Զ���\n");
	fprintf(fout, "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
	fprintf(fout, "����һ:������Զ���\n");

	// ����
	bool flag = false;
	for (int i = 0; i < ITER_TIMES; i++)
	{
		matX = (matA.transpose() * matA).inverse() * (matA.transpose() * matL);

		matRelateOrientElem += matX;

		if (fabs(matX[0][0]) < 3e-5 && fabs(matX[1][0]) < 3e-5 && fabs(matX[2][0]) < 3e-5
			&& fabs(matX[3][0]) < 3e-5 && fabs(matX[4][0]) < 3e-5)
		{
			flag = true;
			std::cout << "��ɣ�����������" << i + 1 << std::endl << "��Զ���Ԫ��" << matRelateOrientElem;
			fprintf(fout, "��Զ���Ԫ��:\nmu    nu    phi    omega    kappa\n%lf    %lf    %lf    %lf    %lf\n",
				matRelateOrientElem[0][0], matRelateOrientElem[1][0], matRelateOrientElem[2][0], matRelateOrientElem[3][0], matRelateOrientElem[4][0]);

			// ����в�
			Matrix matEpsilon = matA * matX - matL;
			fprintf(fout, "�в�(mm)��\nid    V_Q\n");
			for (int i = 0; i < n; i++)
			{
				fprintf(fout, "%d\t    %f\n",(int)imgCoor[i][0], matEpsilon[i][0]);
			}
			
			// �����
			fprintf(fout, "�����(mm)��\n");
			fprintf(fout, "%f\n",sqrt((matEpsilon.transpose() * matEpsilon)[0][0] / n));

			break;
		}

		// ��ֵ
		mu = matRelateOrientElem[0][0];
		nu = matRelateOrientElem[1][0];
		phi = matRelateOrientElem[2][0];
		omega = matRelateOrientElem[3][0];
		kappa = matRelateOrientElem[4][0];

		// ������ռ丨������
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

			// ϵ����
			matA[i][0] = B_x;
			matA[i][1] = -Y_2 / Z_2 * B_x;
			matA[i][2] = -X_2 * Y_2 / Z_2 * N_2;
			matA[i][3] = -(Z_2 + pow(Y_2, 2) / Z_2) * N_2;
			matA[i][4] = X_2 * N_2;
			// ������
			matL[i][0] = N_1 * Y_1 - (N_2 * Y_2 + B_y);
		}

	}

	if (!flag)
	{
		std::cout << "����\n";
		std::cout << "���������ﵽ���ֵ";
		fprintf(fout, "���������ﵽ���ֵ!");
	}


	// ����������������Զ���
	printf("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
	printf("������:��������Զ���\n");
	fprintf(fout, "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
	fprintf(fout, "������:��������Զ���\n");
	// �ڷ�λԪ��
	f;
	x_0;;
	y_0;
	// �����ֵ
	// �ٶ����߳�
	double B = imgCoor[0][1] - imgCoor[0][3];

	// �趨��Զ������
	double phi_1 = 0, kappa_1 = 0, phi_2 = 0, omega_2 = 0, kappa_2 = 0;
	Matrix matR_1 = Matrix::zeros(3, 3);
	Matrix matR_2 = Matrix::zeros(3, 3);


	// ������ռ丨������
	matImgAssistCoor1 = Matrix::zeros(3, 1);
	matImgAssistCoor2 = Matrix::zeros(3, 1);
	X_1 = 0, Y_1 = 0, Z_1 = 0;
	X_2 = 0, Y_2 = 0, Z_2 = 0;
	matL = Matrix::zeros(n, 1);	// �������Q
	matX = Matrix::zeros(5, 1);	// ��Զ���Ԫ�ظ�����mu, nu, phi, omega, kappa;
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

		// ϵ����
		matA[i][0] = X_1 * Y_2 / Z_1;
		matA[i][1] = -X_1;
		matA[i][2] = -X_2 * Y_1 / Z_1;
		matA[i][3] = -(Z_1 + Y_1*Y_2 / Z_1);
		matA[i][4] = X_2;
		// ������q
		matL[i][0] = -f * Y_1 / Z_1 + f * Y_2 / Z_2;
	}



	// ����
	flag = false;
	for (int i = 0; i < ITER_TIMES; i++)
	{
		matX = (matA.transpose() * matA).inverse() * (matA.transpose() * matL);

		matRelateOrientElem += matX;

		if (fabs(matX[0][0]) < 3e-5 && fabs(matX[1][0]) < 3e-5 && fabs(matX[2][0]) < 3e-5
			&& fabs(matX[3][0]) < 3e-5 && fabs(matX[4][0]) < 3e-5)
		{
			flag = true;
			std::cout << "��ɣ�����������" << i + 1 << std::endl << "��Զ���Ԫ��" << matRelateOrientElem;
			fprintf(fout, "��Զ���Ԫ��:\nphi_1    kappa_1    phi_2    omega_2    kappa_2\n%lf    %lf    %lf    %lf    %lf\n",
				matRelateOrientElem[0][0], matRelateOrientElem[1][0], matRelateOrientElem[2][0], matRelateOrientElem[3][0], matRelateOrientElem[4][0]);

			// ����в�
			Matrix matEpsilon = matA * matX - matL;
			fprintf(fout, "�в�(mm)��\nid    V_Q\n");
			for (int i = 0; i < n; i++)
			{
				fprintf(fout, "%d\t    %f\n", (int)imgCoor[i][0], matEpsilon[i][0]);
			}

			// �����
			fprintf(fout, "�����(mm)��\n");
			fprintf(fout, "%f\n", sqrt((matEpsilon.transpose() * matEpsilon)[0][0] / n));

			break;
		}

		// ��ֵ
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

			// ϵ����
			matA[i][0] = X_1 * Y_2 / Z_1;
			matA[i][1] = -X_1;
			matA[i][2] = -X_2 * Y_1 / Z_1;
			matA[i][3] = -(Z_1 + Y_1 * Y_2 / Z_1);
			matA[i][4] = X_2;
			// ������q
			matL[i][0] = -f * Y_1 / Z_1 + f * Y_2 / Z_2;
		}

	}

	if (!flag)
	{
		std::cout << "����\n";
		std::cout << "���������ﵽ���ֵ";
		fprintf(fout, "���������ﵽ���ֵ!");
	}

	Matrix temp = matR.inverse() * matR_1.inverse() - matR_2.inverse();


	fclose(fout);
}