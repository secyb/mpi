// medianFilter_mpi2.cpp : �������̨Ӧ�ó������ڵ㡣
//

#include <gdal_priv.h>  
#include <iostream>
#include "mpi.h"

using namespace std;

void _medianfilter(const unsigned char* image, unsigned char* result, int width, int height);//��ֵ�˲�


int main(int argc, char *argv[])
{
	time_t start, stop, start1, stop1;//ʹ��time_t������ʱ
	start = time(NULL);
	MPI_Status status;
	const char *inPath = "/home/ubuntu/data/GF1_WFV3_E115.7_N28.9_20170729_L1A0002514125.tiff";//ԭʼӰ��·��
	//const char *inPath = "/home/ubuntu/data/1.bmp";
	//const char *outPath = "/home/ubuntu/data/2.bmp";
	const char *outPath = "/home/ubuntu/data/gdal_out.tiff";//���·��
	int nImgSizeX, nImgSizeY, bandcount;
	int rank, size;
	int interval;//ÿ�����̷ֵ���Ӱ��߶�,�ó���ֱ�Ӱ�Ӱ��߶ȷֿ�
	int i, j;
	//�����ڴ�  
	unsigned char **pImgData;//���ڵ�洢Ӱ������
	pImgData = NULL;
	unsigned char **data;//���ڵ�������
	data = NULL;
	unsigned char **result;//���ڵ�洢�������
	result = NULL;
	unsigned char **final_data;
	final_data = NULL;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);//rank�ǽ��̺�
	MPI_Comm_size(MPI_COMM_WORLD, &size);//size��ִ�г���ʱʹ�õĽ�����
	if (!rank) {//������(������0�Ž���)
		//��ͼ��  
		GDALDataset *poDataset;
		GDALAllRegister();
		poDataset = (GDALDataset *)GDALOpen(inPath, GA_ReadOnly);
		if (poDataset == NULL)
		{
			cout << "nothing" << endl;
		}

		nImgSizeX = poDataset->GetRasterXSize();     //��ȡ������Ԫ����  
		nImgSizeY = poDataset->GetRasterYSize();     //��ȡ������Ԫ����  
		bandcount = poDataset->GetRasterCount();     //��ȡ������

		pImgData = new unsigned char *[bandcount];
		for (i = 0; i < bandcount; i++)
			pImgData[i] = new unsigned char[nImgSizeX * nImgSizeY];//�����洢������ԭʼ����

		//Ӱ���и�
		if (nImgSizeY%size)
			MPI_Abort(MPI_COMM_WORLD, 1);
		else
			interval = nImgSizeY / size;

		for (i = 1; i <= bandcount; i++)
		{
			GDALRasterBand * pInRasterBand1 = poDataset->GetRasterBand(i);
			CPLErr error;
			error = pInRasterBand1->RasterIO(GF_Read, 0, 0, nImgSizeX, nImgSizeY, pImgData[i - 1], nImgSizeX, nImgSizeY, GDT_Byte, 0, 0);
			if (error == CE_Failure)
			{
				cout << "��ȡͼ������ʱ����" << endl;
				GDALDestroyDriverManager();
			}
		}
		if (poDataset != NULL)
			delete poDataset;
	}
	start1 = time(NULL);
	//����0�㲥
	MPI_Bcast(&nImgSizeX, 1, MPI_INT, 0, MPI_COMM_WORLD); 
	MPI_Bcast(&nImgSizeY, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&bandcount, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&interval, 1, MPI_INT, 0, MPI_COMM_WORLD);

	if(rank){
		pImgData = new unsigned char *[bandcount];
		for (i = 0; i < bandcount; i++)
			pImgData[i] = new unsigned char[nImgSizeX * nImgSizeY];//�����洢������ԭʼ����
	}

	//���н��̷���洢���ݼ��Ŀռ�
	data = new unsigned char *[bandcount];
	result = new unsigned char *[bandcount];
	for (i = 0; i < bandcount; i++) {
		data[i] = new unsigned char[nImgSizeX * nImgSizeY / size];
		result[i] = new unsigned char[nImgSizeX * nImgSizeY / size];
	}
	final_data = new unsigned char *[bandcount];
	for (i = 0; i < bandcount; i++)
		final_data[i] = new unsigned char[nImgSizeX * nImgSizeY];
	
	for(i=0;i<bandcount;i++)
		MPI_Scatter(&pImgData[i], nImgSizeX * nImgSizeY/size, MPI_UNSIGNED_CHAR, &data[i], nImgSizeX * nImgSizeY/size, MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);	

	MPI_Barrier(MPI_COMM_WORLD);  
	stop1 = time(NULL);
	printf("Use Time:%ld\n", (stop1 - start1)); //��������ʱ�������ڵ��������ڵ㴫���������õ�ʱ��
	for (i = 0; i<bandcount; i++)
		_medianfilter(data[i], result[i], nImgSizeX, interval);

	MPI_Barrier(MPI_COMM_WORLD);  //ͬ��һ��
	//�۽��������н������ݣ������������Լ�����˳���͸�������
	for (i = 0; i<bandcount; i++)
		MPI_Gather(result[i], nImgSizeX * nImgSizeY / size, MPI_UNSIGNED_CHAR, final_data[i], nImgSizeX * nImgSizeY / size, MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);  //ȫ�ռ��ڽ���0

	MPI_Barrier(MPI_COMM_WORLD);  

	if (!rank)
	{
		GDALDataset *poDstDS;
		const char  *pszFormat = "GTiff";
		//const char  *pszFormat = "BMP";
		GDALDriver  *poDriver;
		poDriver = GetGDALDriverManager()->GetDriverByName(pszFormat); //��ȡ����  
		if (poDriver == NULL)
			exit(1);

		//����bandcount�����ε�ͼ��  
		poDstDS = poDriver->Create(outPath, nImgSizeX, nImgSizeY, bandcount, GDT_Byte, NULL);
 
		//������pafScan�е����ݴ�����ͼ��  
		for (i = 1; i <= bandcount; i++)
		{
			poDstDS->GetRasterBand(i)->RasterIO(GF_Write, 0, 0, nImgSizeX, nImgSizeY,
				final_data[i - 1], nImgSizeX, nImgSizeY, GDT_Byte, 0, 0);
		}

		if (poDstDS != NULL)
			delete poDstDS;
	}
	MPI_Finalize();
	stop = time(NULL);
	printf("Use Time:%ld\n", (stop - start));
    return 0;
}

void _medianfilter(const unsigned char* image, unsigned char* result, int width, int height) //��ֵ�˲�
{
	//   Move window through all elements of the image
	for (int m = 1; m < height - 1; ++m)
		for (int n = 1; n < width - 1; ++n)
		{
			//   Pick up window elements
			int k = 0;
			unsigned char window[9];
			for (int j = m - 1; j < m + 2; ++j)
				for (int i = n - 1; i < n + 2; ++i)
					window[k++] = image[j * width + i];
			//   Order elements (only half of them)
			for (int j = 0; j < 5; ++j)
			{
				//   Find position of minimum element
				int min = j;
				for (int l = j + 1; l < 9; ++l)
					if (window[l] < window[min])
						min = l;
				//   Put found minimum element in its place
				const unsigned char temp = window[j];
				window[j] = window[min];
				window[min] = temp; 
			}
			//   Get result - the middle element
			result[m * width + n] = window[4];
		}
}

