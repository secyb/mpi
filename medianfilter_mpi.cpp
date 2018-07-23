// medianFilter_mpi.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"
#include <gdal_priv.h>  
#include <iostream>
#include "mpi.h"

using namespace std;

void _medianfilter(const unsigned char* image, unsigned char* result, int width, int height);//中值滤波


int main(int argc, char *argv[])
{
    time_t start, stop;
    start = time(NULL);
    MPI_Status status;
    const char *inPath = "F:\\data\\GF1_WFV3_E115.7_N28.9_20170729_L1A0002514125.tiff";
    const char *outPath = "F:\\data\\gdal_out.tiff";
    int nImgSizeX, nImgSizeY, bandcount;
    int rank, size;
    int interval;//每个进程分到的影像高度
    int i, j, k;
    //开辟内存  
    unsigned char **pImgData;//主节点存储影像数据
    pImgData = NULL;
    unsigned char **data;//各节点存放数据
    data = NULL;
    unsigned char **result;//各节点存储结果数据
    result = NULL;
    unsigned char **final_data;
    final_data = NULL;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    if (!rank) {
        //打开图像  
        GDALDataset *poDataset;
        GDALAllRegister();
        poDataset = (GDALDataset *)GDALOpen(inPath, GA_ReadOnly);
        if (poDataset == NULL)
        {
            cout << "nothing" << endl;
        }

        nImgSizeX = poDataset->GetRasterXSize();     //获取横向像元个数  
        nImgSizeY = poDataset->GetRasterYSize();     //获取纵向像元个数  
        bandcount = poDataset->GetRasterCount();     //获取波段数

        pImgData = new unsigned char *[bandcount];
        for (i = 0; i < bandcount; i++)
            pImgData[i] = new unsigned char[nImgSizeX * nImgSizeY];

        //影像切割
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
                cout << "读取图像数据时出错！" << endl;
                GDALDestroyDriverManager();
            }
        }
        if (poDataset != NULL)
            delete poDataset;
    }
    MPI_Bcast(&nImgSizeX, 1, MPI_INT, 0, MPI_COMM_WORLD);  //进程0广播
    MPI_Bcast(&nImgSizeY, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&bandcount, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&interval, 1, MPI_INT, 0, MPI_COMM_WORLD);

    //所有进程分配存储数据集的空间
    data = new unsigned char *[bandcount];
    result = new unsigned char *[bandcount];
    for (i = 0; i < bandcount; i++) {
        data[i] = new unsigned char[nImgSizeX * nImgSizeY / size];
        result[i] = new unsigned char[nImgSizeX * nImgSizeY / size];
    }
    final_data = new unsigned char *[bandcount];
    for (i = 0; i < bandcount; i++)
        final_data[i] = new unsigned char[nImgSizeX * nImgSizeY];
        
    if (!rank)//主进程分得的数据
    {
        for (i = 0; i<bandcount; i++)
        {
            for (j = 0; j<nImgSizeX * nImgSizeY / size; j++)
            {
                data[i][j] = pImgData[i][j];
                result[i][j] = data[i][j];
            }
        }
    }
    
    if (!rank) {  //进程0向其他进程分配数据集
        for (i = 0; i<bandcount; i++)
        {
            for (j = nImgSizeX*nImgSizeY/size; j<nImgSizeX * nImgSizeY; j++)
            {
                MPI_Send(&pImgData[i][j], 1, MPI_UNSIGNED_CHAR, j/(nImgSizeX*nImgSizeY/size), 99, MPI_COMM_WORLD);
                //data[i][j] = pImgData[i][j + rank*nImgSizeX * nImgSizeY / size];
                //result[i][j] = data[i][j];
            }
        }
    }
    else {  //其他进程接收进程0数据
        for (i = 0; i < bandcount; i++) {
            for (j = 0; j < nImgSizeX*nImgSizeY / size; j++) {
                MPI_Recv(&data[i][j], 1, MPI_UNSIGNED_CHAR, 0, 99, MPI_COMM_WORLD, &status);
            }
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);  //同步一下
    /*if (rank)//其他进程分得的数据
    {
        for (i = 0; i<bandcount; i++)
        {
            for (j = 0; j<nImgSizeX * nImgSizeY / size; j++)
            {
                data[i][j] = pImgData[i][j + rank*nImgSizeX * nImgSizeY / size];
                result[i][j] = data[i][j];
            }
        }
    }*/
    //MPI_Barrier(MPI_COMM_WORLD);  

    /*for (i = 0; i<bandcount; i++)
        _medianfilter(data[i], result[i], nImgSizeX, interval);*/

    //MPI_Barrier(MPI_COMM_WORLD);  //同步一下
    /*for (i = 0; i<bandcount; i++)
        MPI_Gather(result[i], nImgSizeX * nImgSizeY / size, MPI_UNSIGNED_CHAR, final_data[i], nImgSizeX * nImgSizeY / size, MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);  //全收集于进程0
*/

    for (i = 0; i<bandcount; i++)
        MPI_Gather(data[i], nImgSizeX * nImgSizeY / size, MPI_UNSIGNED_CHAR, final_data[i], nImgSizeX * nImgSizeY / size, MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);  //全收集于进程0

    MPI_Barrier(MPI_COMM_WORLD);  

    if (!rank)
    {
        GDALDataset *poDstDS;
        const char  *pszFormat = "GTiff";
        GDALDriver  *poDriver;
        poDriver = GetGDALDriverManager()->GetDriverByName(pszFormat); //获取驱动  
        if (poDriver == NULL)
            exit(1);

        //创建bandcount个波段的图像  
        poDstDS = poDriver->Create(outPath, nImgSizeX, nImgSizeY, bandcount, GDT_Byte, NULL);
 
        //将缓存pafScan中的数据存入结果图像  
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

void _medianfilter(const unsigned char* image, unsigned char* result, int width, int height)
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

