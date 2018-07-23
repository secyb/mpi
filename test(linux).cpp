// test.cpp : Defines the entry point for the console application.
//
#include <mpi.h>
#include <stdio.h>
#include <gdal.h>
#include <gdal_priv.h>
#include <ogrsf_frmts.h>
#include <time.h>
#include <stdlib.h>
#define  MASTER		0

struct Node
{
	OGRRawPoint point;
	Node *next;
};

struct Counter
{
	int num;
	Node *next;
};


//The structure of chains.
void CorrectDirection(OGRRawPoint low[], OGRRawPoint high[], int SoL, int SoH, int &StartL, int &StartH);

//Use the two points that is the most near ones as the begining of the insertion.
void FindStarts(OGRRawPoint low[], OGRRawPoint high[], int SoL, int SoH, int &StartL, int &StartH)
{
	double distance, x, y, mindis = 767297;
	for (int i = 0; i < SoL; i++)
	{
		for (int j = 0; j < SoH; j++)
		{
			x = low[i].x - high[j].x;
			y = low[i].y - high[j].y;
			distance = x * x + y * y;
			if (distance < mindis)
			{
				StartL = i;
				StartH = j;
				mindis = distance;
			}
		}
	}

	CorrectDirection(low, high, SoL, SoH, StartL, StartH);
}

Node* BeforFirstInsert(OGRRawPoint Point, OGRRawPoint TheNext)
{
	Node *target = new Node;
	target->next = NULL;
	target->point.x = Point.x + (Point.x - TheNext.x);
	target->point.y = Point.y + (Point.y - TheNext.y);
	return target;
}

bool TestSide(OGRRawPoint PointA, OGRRawPoint PointB, OGRRawPoint Point1, OGRRawPoint Point2)
{
	double xa, xb, ya, yb, x1, x2, y1, y2, S1, S2, a, b;
	xa = PointA.x;
	xb = PointB.x;
	x1 = Point1.x;
	x2 = Point2.x;
	ya = PointA.y;
	yb = PointB.y;
	y1 = Point1.y;
	y2 = Point2.y;
	a = (yb - ya) / (xb - xa);
	b = (ya*xb - yb * xa) / (xb - xa);
	S1 = y1 - a * x1 - b;
	S2 = y2 - a * x2 - b;

	if (S1*S2 < 1e-7)
	{
		return false;
	}
	else
	{
		return true;
	}

}

bool TestSideI(OGRRawPoint PointA, OGRRawPoint PointB, OGRRawPoint Point1, OGRRawPoint Point2)
{
	double xa, xb, ya, yb, x1, x2, y1, y2, S1, S2, a, b;
	xa = PointA.x;
	xb = PointB.x;
	x1 = Point1.x;
	x2 = Point2.x;
	ya = PointA.y;
	yb = PointB.y;
	y1 = Point1.y;
	y2 = Point2.y;
	a = (yb - ya) / (xb - xa);
	b = (ya*xb - yb * xa) / (xb - xa);
	S1 = y1 - a * x1 - b;
	S2 = y2 - a * x2 - b;

	if (S1*S2 > -1e-7)
	{
		return true;
	}
	else
	{
		return false;
	}

}

//Try to detect whether the two counters is going at the same direction, if not, revers the higher one.
void CorrectDirection(OGRRawPoint low[], OGRRawPoint high[], int SoL, int SoH, int &StartL, int &StartH)
{
	double flagL1, flagL2, flagH1, flagH2, a, b;
	bool test1, test2;
	a = (low[StartL].y - high[StartH].y) / (low[StartL].x - high[StartH].x);
	b = (high[StartH].y * low[StartL].x - low[StartL].y * high[StartH].x) / (low[StartL].x - high[StartH].x);

	flagL1 = low[StartL + 1].y - a * low[StartL + 1].x - b;
	flagL2 = low[StartL - 1].y - a * low[StartL - 1].x - b;
	flagH1 = high[StartH + 1].y - a * high[StartH + 1].x - b;
	flagH2 = high[StartH - 1].y - a * high[StartH - 1].x - b;
	//test1 = TestSide(low[StartL], high[StartH], low[StartL + 1], high[StartH + 1]);
	//test2 = TestSide(low[StartL], high[StartH], low[StartL - 1], high[StartH - 1]);

	if (flagL1*flagH1 > 0 || flagL2 * flagH2 > 0)
	{
		if (abs(low[0].x - low[SoL - 1].x) < 0.2 && abs(low[0].y - low[SoL - 1].y) < 0.2)
		{
			OGRRawPoint *tempL = new OGRRawPoint[SoL];
			for (int i = 0; i < SoL; i++)
			{
				tempL[i] = low[(StartL + i) % SoL];
			}
			for (int i = 0; i < SoL; i++)
			{
				low[i] = tempL[i];
			}
			StartL = 0;
			free(tempL);
		}

		if (abs(high[0].x - high[SoH - 1].x) < 0.2 && abs(high[0].y - high[SoH - 1].y) < 0.2)
		{
			OGRRawPoint *tempH = new OGRRawPoint[SoH];
			for (int i = 0; i < SoH; i++)
			{
				tempH[i] = high[(StartH + i) % SoH];
			}
			for (int i = 0; i < SoH; i++)
			{
				high[i] = tempH[i];
			}
			StartH = 0;
			free(tempH);
		}
	}
	else
	{
		OGRRawPoint *temp = new OGRRawPoint[SoH];
		for (int i = 0; i < SoH; i++)
		{
			temp[i] = high[SoH - i - 1];
		}
		for (int i = 0; i < SoH; i++)
		{
			high[i] = temp[i];
		}
		StartH = SoH - 1 - StartH;
		free(temp);

		if (abs(low[0].x - low[SoL - 1].x) < 0.2 && abs(low[0].y - low[SoL - 1].y) < 0.2)
		{
			OGRRawPoint *tempL = new OGRRawPoint[SoL];
			for (int i = 0; i < SoL; i++)
			{
				tempL[i] = low[(StartL + i) % SoL];
			}
			for (int i = 0; i < SoL; i++)
			{
				low[i] = tempL[i];
			}
			StartL = 0;
			free(tempL);
		}

		if (abs(high[0].x - high[SoH - 1].x) < 0.2 && abs(high[0].y - high[SoH - 1].y) < 0.2)
		{
			OGRRawPoint *tempH = new OGRRawPoint[SoH];
			for (int i = 0; i < SoH; i++)
			{
				tempH[i] = high[(StartH + i) % SoH];
			}
			for (int i = 0; i < SoH; i++)
			{
				high[i] = tempH[i];
			}
			StartH = 0;
			free(tempH);
		}
	}

}

Node* GetPrevious(Counter CounterAimed)
{
	Node *target = NULL;
	for (int i = 0; i < CounterAimed.num - 1; i++)
	{
		if (i == 0)
		{
			target = CounterAimed.next;
		}
		else
		{
			target = target->next;
		}
	}

	return target;
}

Node* GetLast(Counter CounterAimed)
{
	Node *target = NULL;
	for (int i = 0; i < CounterAimed.num; i++)
	{
		if (i == 0)
		{
			target = CounterAimed.next;
		}
		else
		{
			target = target->next;
		}
	}

	return target;
}

Node* FindNode(Counter CounterAll, int num)
{
	Node *target = NULL;
	if (num <= CounterAll.num)
	{
		for (int i = 0; i < num; i++)
		{
			if (i == 0)
			{
				target = CounterAll.next;
			}
			else
			{
				target = target->next;
			}
		}
	}

	return target;
}

bool TestIntersect(OGRRawPoint Start1, OGRRawPoint End1, OGRRawPoint Start2, OGRRawPoint End2)
{
	/*
	double t, u, S1_S2_x, S1_S2_y, S1_x, S1_y, S2_x, S2_y;
	S1_S2_x = Start1.x - Start2.x;
	S1_S2_y = Start1.y - Start2.y;
	S1_x = End1.x - Start1.x;
	S1_y = End1.y - Start1.y;
	S2_x = End2.x - Start2.x;
	S2_y = End2.y - Start2.y;

	t = (S1_S2_y * S2_x - S1_S2_x * S2_y) / (S1_x * S2_y - S1_y * S2_x);
	u = (S1_x * S1_S2_y - S1_y * S1_S2_x) / (S1_x * S2_y - S1_y * S2_x);

	if (t > 0&& t < 1 && u > 0 && u < 1)
	{
	return true;
	}
	else
	return false;
	*/

	bool flag1, flag2;
	flag1 = TestSideI(Start1, End1, Start2, End2);
	flag2 = TestSideI(Start2, End2, Start1, End1);
	if (flag1 || flag2)
	{
		return false;
	}
	else
	{
		return true;
	}

}

double CosLaw(OGRRawPoint PointA, OGRRawPoint PointB, OGRRawPoint PointC)
{
	double a, b, c, CosA;
	a = (PointB.x - PointC.x)*(PointB.x - PointC.x) + (PointB.y - PointC.y)*(PointB.y - PointC.y);
	b = (PointA.x - PointC.x)*(PointA.x - PointC.x) + (PointA.y - PointC.y)*(PointA.y - PointC.y);
	c = (PointB.x - PointA.x)*(PointB.x - PointA.x) + (PointB.y - PointA.y)*(PointB.y - PointA.y);
	CosA = (c + b - a) / (2 * sqrt(c) * sqrt(b));
	return CosA;
}

Counter Reverse(Counter CounterT)
{
	Node *pointT, *pointS;
	Counter Target;

	Target.num = 0;
	Target.next = NULL;

	for (int i = 0; i < CounterT.num; i++)
	{
		if (!Target.next)
		{
			pointS = GetLast(CounterT);
			Target.next = pointS;
			pointS->next = NULL;
			Target.num++;
		}
		else
		{
			pointS = FindNode(CounterT, CounterT.num - i);
			pointT = GetLast(Target);
			pointT->next = pointS;
			pointS->next = NULL;
			Target.num++;
		}
	}

	return Target;
}

Counter* insert(OGRRawPoint low[], OGRRawPoint high[], int SoL, int SoH, int number)
{
	Counter *Counter1 = new Counter[number];
	Node *compared;
	double LCos = 1, HCos = 1, Dtemp;
	bool flagL = false, flagH = false, FlagofIntersect, LChange = false, HChange = false;
	//Lnum is the next point at the lower counter. A is the current point for the lower counter. LFlag is uesed to recomande the flag point.
	int LNum = 1, HNum = 1, A = 0, B = 0, LFlag = 0, HFlag = 0;

	//Counter 1
	//initial
	for (int i = 0; i < number; i++)
	{
		Counter1[i].num = 0;
		Counter1[i].next = NULL;
	}

	Node *storage;

	//The first insert node.
	if (SoL > 0 && SoH > 0)
	{
		for (int i = 0; i < number; i++)
		{
			Counter1[i].next = new Node;
			storage = Counter1[i].next;
			storage->point.x = low[0].x + (high[0].x - low[0].x) * (i + 1) / (number + 1);
			storage->point.y = low[0].y + (high[0].y - low[0].y) * (i + 1) / (number + 1);
			Counter1[i].num++;
			storage->next = NULL;
			storage = NULL;
		}
	}


	//Second step.
	for (; LNum < SoL && HNum < SoH;)
	{

		//Calculate fot the flag of the lower counter.
		for (int temp = LNum; temp < SoL; temp++)
		{
			storage = GetPrevious(Counter1[0]);
			if (storage != NULL)
			{
				flagL = TestSide(low[A], high[B], storage->point, low[temp]);
			}
			else
			{
				storage = BeforFirstInsert(low[A], low[A + 1]);
				flagL = TestSide(low[A], high[B], storage->point, low[temp]);
				free(storage);

			}

			if (flagL == false)
			{
				if (!TestIntersect(low[A], low[A + 1], high[B], low[temp]))
				{
					if (!TestIntersect(low[A], low[temp], high[B], high[B + 1]))
					{
						Dtemp = CosLaw(low[temp], low[A], high[B]);
						if (Dtemp < LCos)
						{
							LCos = Dtemp;
							LFlag = temp;
							//LChange = true;
						}
					}
				}
			}
		}

		//Calculate fot the flag of the higher counter.
		for (int temp = HNum; temp < SoH; temp++)
		{
			storage = GetPrevious(Counter1[0]);
			if (storage != NULL)
			{
				flagH = TestSide(low[A], high[B], storage->point, high[temp]);
			}
			else
			{
				storage = BeforFirstInsert(high[B], high[B + 1]);
				flagH = TestSide(low[A], high[B], storage->point, high[temp]);
				free(storage);

			}
			if (flagH == false)
			{
				if (!TestIntersect(low[A], high[temp], high[B], high[B + 1]))
				{
					if (!TestIntersect(low[A], low[A + 1], high[B], high[temp]))
					{
						Dtemp = CosLaw(high[temp], low[A], high[B]);
						if (Dtemp < HCos)
						{
							HCos = Dtemp;
							HFlag = temp;
							//HChange = true;
						}
					}
				}
			}
		}

		//Test if ever changed.
		/*
		if (LChange == false && HChange == false)
		{
		break;
		}
		*/
		//Insert the target counter.
		if (LCos < HCos)
		{
			for (int i = A + 1; i < LFlag; i++)
			{
				for (int ii = 0; ii < SoL; ii++)
				{
					FlagofIntersect = TestIntersect(high[B], low[i], low[ii], low[ii + 1]);
					if (FlagofIntersect == true)
					{
						break;
					}
				}
				if (FlagofIntersect == false)
				{
					for (int j = 0; j < number; j++)
					{
						storage = GetLast(Counter1[j]);
						storage->next = new Node;
						storage->next->point.x = low[i].x + (high[B].x - low[i].x) * (j + 1) / (number + 1);
						storage->next->point.y = low[i].y + (high[B].y - low[i].y) * (j + 1) / (number + 1);
						Counter1[j].num++;
						storage->next->next = NULL;
					}
				}
			}
			A = LFlag;
			LNum = A + 1;
		}
		else
		{
			for (int i = B + 1; i < HFlag; i++)
			{
				for (int ii = 0; ii < SoH; ii++)
				{
					FlagofIntersect = TestIntersect(high[i], low[A], high[ii], high[ii + 1]);
					if (FlagofIntersect == true)
					{
						break;
					}
				}
				if (FlagofIntersect == false)
				{
					for (int j = 0; j < number; j++)
					{
						storage = GetLast(Counter1[j]);
						storage->next = new Node;
						storage->next->point.x = low[A].x + (high[i].x - low[A].x) * (j + 1) / (number + 1);
						storage->next->point.y = low[A].y + (high[i].y - low[A].y) * (j + 1) / (number + 1);
						Counter1[j].num++;
						storage->next->next = NULL;
					}
				}
			}
			B = HFlag;
			HNum = B + 1;
		}
		// Insert counters between A and B.

		for (int j = 0; j < number; j++)
		{
			storage = GetLast(Counter1[j]);
			storage->next = new Node;
			storage->next->point.x = low[A].x + (high[B].x - low[A].x) * (j + 1) / (number + 1);
			storage->next->point.y = low[A].y + (high[B].y - low[A].y) * (j + 1) / (number + 1);
			Counter1[j].num++;
			storage->next->next = NULL;
		}

		// Reset the initial data.

		LCos = 1;
		HCos = 1;
		HChange = false;
		LChange = false;

	}

	//Third step.
	if (LNum < SoL && HNum >= SoH)
	{
		for (int i = LNum; i < SoL; i++)
		{
			for (int j = 0; j < number; j++)
			{
				storage = new Node;
				GetLast(Counter1[j])->next = storage;
				storage->point.x = low[i].x + (high[HNum - 1].x - low[i].x) * (j + 1) / (number + 1);
				storage->point.y = low[i].y + (high[HNum - 1].y - low[i].y) * (j + 1) / (number + 1);
				Counter1[j].num++;
				storage->next = NULL;
				storage = NULL;
			}
		}
	}

	if (HNum < SoH && LNum >= SoL)
	{
		for (int i = HNum; i < SoH; i++)
		{
			for (int j = 0; j < number; j++)
			{
				storage = new Node;
				GetLast(Counter1[j])->next = storage;
				storage->point.x = low[LNum - 1].x + (high[i].x - low[LNum - 1].x) * (j + 1) / (number + 1);
				storage->point.y = low[LNum - 1].y + (high[i].y - low[LNum - 1].y) * (j + 1) / (number + 1);
				Counter1[j].num++;
				storage->next = NULL;
				storage = NULL;
			}
		}
	}

	return Counter1;
}

Counter* Merge(OGRRawPoint low[], OGRRawPoint high[], int SoL, int SoH, int number)
{
	Counter *Counter1, *Counter2, *CounterAll = new Counter[number];
	Node *Pointer1, *Pointer2;
	OGRRawPoint *lowCounter1, *lowCounter2, *highCounter1, *highCounter2;
	int StartL, StartH;
	FindStarts(low, high, SoL, SoH, StartL, StartH);


	highCounter1 = new OGRRawPoint[SoH - StartH];
	for (int i = 0; i < SoH - StartH; i++)
	{
		highCounter1[i] = high[i + StartH];
	}

	highCounter2 = new OGRRawPoint[StartH];
	for (int i = 0; i < StartH; i++)
	{
		highCounter2[i] = high[StartH - 1 - i];
	}


	lowCounter1 = new OGRRawPoint[SoL - StartL];
	for (int i = 0; i < SoL - StartL; i++)
	{
		lowCounter1[i] = low[i + StartL];
	}

	lowCounter2 = new OGRRawPoint[StartL];
	for (int i = 0; i < StartL; i++)
	{
		lowCounter2[i] = low[StartL - 1 - i];
	}

	Counter1 = insert(lowCounter1, highCounter1, SoL - StartL, SoH - StartH, number);
	Counter2 = insert(lowCounter2, highCounter2, StartL, StartH, number);

	for (int j = 0; j < number; j++)
	{
		if (Counter1[j].num == 0)
		{
			CounterAll[j] = Counter2[j];
		}
		else if (Counter2[j].num == 0)
		{
			CounterAll[j] = Counter1[j];
		}
		else
		{
			CounterAll[j] = Reverse(Counter2[j]);
			Pointer1 = GetLast(CounterAll[j]);
			Pointer1->next = Counter1[j].next;
		}
	}


	return CounterAll;
}

void Prepare(const char *Filename, OGRRawPoint *&low, OGRRawPoint *&high, int &SoL, int &SoH)
{
	GDALDataset *poVec;
	OGRLayer *poLayer;
	OGRFeature *poFeature;
	OGRGeometry *poGeo;
	OGRLineString *Line;
	OGRFieldType see;
	const char *Field;
	//OGRPoint *Low, *high;
	//OGRRawPoint *low, *high;
	//OGRSFDriver* poDriver = OGRSFDriverRegistrar::GetRegistrar()->GetDriverByName("ESRI Shapefile");


	//int PointNum;

	poVec = (GDALDataset*)GDALOpenEx("D:\\Study\\GraduatDesign\\test\\Out_Put\\Counter_out_test.shp", GDAL_OF_VECTOR, NULL, NULL, NULL);
	if (poVec == NULL)
	{
		printf("Open failed.\n");
		exit(1);
	}
	poLayer = poVec->GetLayer(0);
	poFeature = poLayer->GetFeature(0);
	Field = poFeature->GetFieldAsString("ELE");
	see = poFeature->GetFieldDefnRef(2)->GetType();
	printf("%s\n", Field);
	poGeo = poFeature->GetGeometryRef();
	Line = (OGRLineString*)poGeo;
	SoL = Line->getNumPoints();
	low = new OGRRawPoint[SoL];
	Line->getPoints(low);
	poFeature = poLayer->GetFeature(1);
	printf("%s\n", poFeature->GetFieldAsString("ELE"));
	poGeo = poFeature->GetGeometryRef();
	Line = (OGRLineString*)poGeo;
	SoH = Line->getNumPoints();
	high = new OGRRawPoint[SoH];
	Line->getPoints(high);

	GDALClose(poVec);

}

void testdata(const char *File)
{




	//read varities
	GDALDataset *poVec;
	OGRLayer *poLayer;
	OGRFeature *poFeature, *newFeature;
	OGRFieldDefn *poFD;
	OGRGeometry *poGeometry;
	OGRLineString *poLine;
	OGRPoint *poPoint = new OGRPoint;
	int count;
	double adfGeoTransform[6];

	//write varities
	GDALDriver *poDriver;
	GDALDataset *poWriDS;
	OGRLayer *poOutLayer;

	OGRSpatialReference *poSR;



	GDALAllRegister();


	//---read code---
	poVec = (GDALDataset*)GDALOpenEx(File, GDAL_OF_VECTOR, NULL, NULL, NULL);

	if (poVec == NULL)
	{
		printf("Open failed.\n");
		exit(1);
	}

	count = poVec->GetLayerCount();
	poLayer = poVec->GetLayer(0);
	poSR = poLayer->GetSpatialRef();
	poFeature = poLayer->GetFeature(55256);
	poFeature->GetFieldAsInteger(2);
	printf("%d\n", poFeature->GetFieldCount());
	poFD = poFeature->GetFieldDefnRef(0);
	poGeometry = poFeature->GetGeometryRef();
	//printf("%s\n", poFeature->GetFieldDefnRef(2));
	poVec->GetGeoTransform(adfGeoTransform);
	poLine = (OGRLineString*)poGeometry;
	poLine->StartPoint(poPoint);

	printf("%d, %d", poPoint->getX(), poPoint->getY());


	//---write code---
	poDriver = GetGDALDriverManager()->GetDriverByName("ESRI Shapefile");
	if (poDriver == NULL)
	{
		printf("Driver not available.");
		exit(1);
	}


	poWriDS = poDriver->Create("D:\\Study\\GraduatDesign\\test\\Out_Put\\Counter_out_test.shp", 0, 0, 0, GDT_Unknown, NULL);
	poOutLayer = poWriDS->CreateLayer("Conter", poSR, poGeometry->getGeometryType());
	for (int i = 0; i < poFeature->GetFieldCount(); i++)
	{
		poOutLayer->CreateField(poFeature->GetFieldDefnRef(i));
	}
	poOutLayer->CreateFeature(poFeature);
	poFeature = poLayer->GetFeature(58802);
	poOutLayer->CreateFeature(poFeature);

	GDALClose(poWriDS);
	GDALClose(poVec);

}

void WriteCounter(const char *Filename, Counter *CounterAll, int number, double EleoL, double EleoH)
{
	//read varities
	GDALDataset *poVec;
	OGRLayer *poLayer;
	OGRFeature *poFeature, *newFeature;
	OGRFieldDefn *poFD;
	OGRGeometry *poGeometry;
	OGRLineString **poLine;
	OGRPoint *poPoint = new OGRPoint;
	int count;
	double adfGeoTransform[6], elevation;

	//write varities
	GDALDriver *poDriver;
	GDALDataset *poWriDS;
	OGRLayer *poOutLayer;
	OGRErr ERR;


	OGRSpatialReference *poSR;

	Node *temp;

	GDALAllRegister();


	//---read code---
	poVec = (GDALDataset*)GDALOpenEx(Filename, GDAL_OF_VECTOR && GDAL_OF_UPDATE, NULL, NULL, NULL);

	if (poVec == NULL)
	{
		printf("Open failed.\n");
		exit(1);
	}

	count = poVec->GetLayerCount();
	poLayer = poVec->GetLayer(0);
	poSR = poLayer->GetSpatialRef();
	poFeature = poLayer->GetFeature(0);
	poGeometry = poFeature->GetGeometryRef();

	poLine = new OGRLineString *[number];

	for (int i = 0; i < number; i++)
	{
		poLine[i] = new OGRLineString;
		for (int j = 0; j < CounterAll[i].num; j++)
		{
			temp = FindNode(CounterAll[i], j + 1);
			poLine[i]->addPoint(temp->point.x, temp->point.y);
		}
	}


	//---write code---
	poDriver = GetGDALDriverManager()->GetDriverByName("ESRI Shapefile");
	if (poDriver == NULL)
	{
		printf("Driver not available.");
		exit(1);
	}


	poWriDS = poDriver->Create("..\\Out_Put\\Counter_output.shp", 0, 0, 0, GDT_Unknown, NULL);
	if (poWriDS == NULL)
	{
		printf("Write Wrong");
		exit(1);
	}
	poOutLayer = poWriDS->CreateLayer("Conter2", poSR, poGeometry->getGeometryType());


	for (int i = 0; i < poFeature->GetFieldCount(); i++)
	{
		poOutLayer->CreateField(poFeature->GetFieldDefnRef(i));
	}

	for (int i = 0; i < number; i++)
	{
		newFeature = new OGRFeature(poFeature->GetDefnRef());
		newFeature->SetGeometry(poLine[i]);
		elevation = EleoL + (i + 1.0) / (number + 1.0) * (EleoH - EleoL);
		newFeature->SetField("ELE", elevation);

		//poOutLayer->CreateFeature(newFeature);
		ERR = poLayer->CreateFeature(newFeature);
	}

	//poLayer->CreateFeature(newFeature);

	GDALClose(poVec);
	GDALClose(poWriDS);
}

void GetCounter(OGRLayer *poLayer, int NoL, int NoH, OGRRawPoint *&low, OGRRawPoint *&high, int &SoL, int &SoH, double &EleoL, double &EleoH)
{
	OGRFeature *poFeature;
	OGRGeometry *poGeo;
	OGRLineString *Line;
	const char *Field;

	poFeature = poLayer->GetFeature(NoL);
	poGeo = poFeature->GetGeometryRef();
	Line = (OGRLineString*)poGeo;
	EleoL = poFeature->GetFieldAsDouble("ELE");
	SoL = Line->getNumPoints();
	low = new OGRRawPoint[SoL];
	Line->getPoints(low);
	poFeature = poLayer->GetFeature(NoH);
	poGeo = poFeature->GetGeometryRef();
	Line = (OGRLineString*)poGeo;
	EleoH = poFeature->GetFieldAsDouble("ELE");
	SoH = Line->getNumPoints();
	high = new OGRRawPoint[SoH];
	Line->getPoints(high);

}

void GetCounter(OGRLayer *poLayer, int Order, OGRRawPoint *&Counter, int &size)
{
	OGRFeature *poFeature;
	OGRGeometry *poGeo;
	OGRLineString *Line;
	const char *Field;

	poFeature = poLayer->GetFeature(Order);
	poGeo = poFeature->GetGeometryRef();
	Line = (OGRLineString*)poGeo;
	size = Line->getNumPoints();
	Counter = new OGRRawPoint[size];
	Line->getPoints(Counter);
}

void FindLowest(OGRLayer *poLayer, int &target, double &elevation)
{
	double temp;
	elevation = 10000.0;
	;
	for (int i = 0; i < poLayer->GetFeatureCount(); i++)
	{
		temp = poLayer->GetFeature(i)->GetFieldAsDouble("ELE");
		if (temp < elevation)
		{
			elevation = temp;
			target = i;
		}
	}
}

void FindLow(OGRLayer *poLayer, int &target, double &min)
{
	target = -1;
	double limit = 100000.0, temp;
	for (int i = 0; i < poLayer->GetFeatureCount(); i++)
	{
		temp = poLayer->GetFeature(i)->GetFieldAsDouble("ELE");
		if (temp > min && temp < limit)
		{
			limit = temp;
			target = i;
		}
	}

	min = limit;

}

void FastSort(int Order[], double Elevation[], int start, int end)
{
	if (start >= end)
	{
		return;
	}
	int first = start;
	int last = end;
	double key = Elevation[first], key2 = Order[first];

	while (first < last)
	{
		while (first < last && Elevation[last] >= key)
		{
			--last;
		}

		Elevation[first] = Elevation[last];
		Order[first] = Order[last];

		while (first < last && Elevation[first] <= key)
		{
			++first;
		}

		Elevation[last] = Elevation[first];
		Order[last] = Order[first];

		Elevation[first] = key;
		Order[first] = key2;
		FastSort(Order, Elevation, start, first - 1);
		FastSort(Order, Elevation, first + 1, end);

	}

}

void SortingOrder(OGRLayer *poLayer, int *&order, double *&elevation)
{
	int num = poLayer->GetFeatureCount();
	order = new int[num];
	elevation = new double[num];

	for (int i = 0; i < num; i++)
	{
		order[i] = i;
		elevation[i] = poLayer->GetFeature(i)->GetFieldAsDouble("ELE");
	}

	FastSort(order, elevation, 0, num - 1);

}

void Chain2Arry(Counter CounterAll, OGRRawPoint *&Arry)
{
	Node *temp = NULL;
	Arry = new OGRRawPoint[CounterAll.num];
	for (int i = 1; i <= CounterAll.num; i++)
	{
		temp = FindNode(CounterAll, i);
		Arry[i - 1].x = temp->point.x;
		Arry[i - 1].y = temp->point.y;
	}
}

void Arry2Chain(OGRRawPoint *Arry, Counter &CounterAll, int num)
{
	Node *temp = NULL;
	for (int i = 0; i < num; i++)
	{

		if (i == 0)
		{
			temp = new Node;
			CounterAll.next = temp;
			CounterAll.num = 1;
			temp->point = Arry[i];
			temp->next = NULL;
		}
		else
		{
			temp = GetLast(CounterAll);
			temp->next = new Node;
			temp->next->point = Arry[i];
			temp->next->next = NULL;
			CounterAll.num++;
		}
	}
}

void CopyCounter(OGRRawPoint *&low, OGRRawPoint *high, int SizeofH)
{
	low = new OGRRawPoint[SizeofH];
	for (int i = 0; i < SizeofH; i++)
	{
		low[i] = high[i];
	}
}

void TheInsert(const char* TheFile, int number, int MyId, int MPISize)
{
	OGRRawPoint *low = NULL, *high = NULL, *middle = NULL, **CounterData = NULL;
	Counter *CounterAll = NULL;

	//OGRFieldType see;
	MPI_Datatype *MPI_Point = new MPI_Datatype;
	MPI_Status *Status = new MPI_Status;
	OGRRawPoint *Arry = NULL;
	//const char *Field;
	const char *filename = "D:\\Study\\GraduatDesign\\test\\Out_Put\\Counter_out_test.shp";
	int SoL, SoH, StartL, StartH, SizeofF, NofL = -1, Group, NofH = -1, *Order = NULL, TheStart, TheEnd, *SizeGroup, CounterSize;
	double min, temp, *Elevation = NULL;
	OGRLayer *poLayer = NULL;
	MPI_Type_contiguous(2, MPI_DOUBLE, MPI_Point);

	MPI_Type_commit(MPI_Point);

	GDALAllRegister();
	OGRRegisterAll();

	//Read the File and send the Point data out in the main process.And write the file after recieve the counter.
	if (MyId == 0)
	{
		GDALDataset *poVec;

		OGRFeature *poFeature;
		OGRGeometry *poGeo;
		OGRLineString *Line;
		int counterSize = NULL;

		poVec = (GDALDataset*)GDALOpenEx(filename, GDAL_OF_VECTOR, NULL, NULL, NULL);
		if (poVec == NULL)
		{
			printf("Open failed.\n");
			exit(1);
		}
		poLayer = poVec->GetLayer(0);
		SizeofF = poLayer->GetFeatureCount();
		Group = (SizeofF - 1) / MPISize + 1;
		MPI_Bcast(&Group, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(&SizeofF, 1, MPI_INT, 0, MPI_COMM_WORLD);
		SizeGroup = new int[Group];

		SortingOrder(poLayer, Order, Elevation);

		for (int i = 0; i < Group; i++)
		{
			poFeature = poLayer->GetFeature(Order[i]);
			poGeo = poFeature->GetGeometryRef();
			Line = (OGRLineString*)poGeo;
			SizeGroup[i] = Line->getNumPoints();
		}

		for (int i = 1; i < MPISize; i++)
		{
			TheStart = i * Group;
			TheEnd = ((i + 1) * Group) < SizeofF ? ((i + 1) * Group) : SizeofF;
			int temp = 0;
			for (int j = TheStart; j < TheEnd; j++)
			{
				GetCounter(poLayer, Order[j], low, SoL);
				MPI_Send(&SoL, 1, MPI_INT, i, temp, MPI_COMM_WORLD);
				MPI_Send(&low, SoL, *MPI_Point, i, temp, MPI_COMM_WORLD);
				temp++;
			}
		}
		TheStart = 0;
		TheEnd = Group < SizeofF ? Group : SizeofF;
	}
	else
	{
		MPI_Bcast(&Group, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(&SizeofF, 1, MPI_INT, 0, MPI_COMM_WORLD);
		SizeGroup = new int[Group];
		TheStart = MyId * Group;
		TheEnd = ((MyId + 1) * Group) < SizeofF ? ((MyId + 1) * Group) : SizeofF;
		CounterData = new OGRRawPoint*[Group];
	}

	//MPI_Barrier(MPI_COMM_WORLD);


	for (int i = TheStart; i < TheEnd; i++)
	{
		if (MyId != 0)
		{
			MPI_Recv(&SizeGroup[i], 1, MPI_INT, 0, i, MPI_COMM_WORLD, Status);
			CounterData[i] = new OGRRawPoint[SizeGroup[i]];
			MPI_Recv(&CounterData[i], SizeGroup[i], *MPI_Point, 0, i, MPI_COMM_WORLD, Status);




			if (i == 0)
			{
				low = new OGRRawPoint[SizeGroup[i]];
				MPI_Recv(&low, SizeGroup[i], *MPI_Point, 0, i, MPI_COMM_WORLD, Status);
			}
			else
			{
				high = new OGRRawPoint[SizeGroup[i]];
				MPI_Recv(&high, SizeGroup[i], *MPI_Point, 0, i, MPI_COMM_WORLD, Status);
				CounterAll = Merge(low, high, SizeGroup[i - 1], SizeGroup[i], number);
				for (int j = 0; j < number; j++)
				{
					Chain2Arry(CounterAll[j], Arry);
					MPI_Send(&CounterAll[j].num, 1, MPI_INT, 0, j, MPI_COMM_WORLD);
					MPI_Send(&Arry, CounterAll[j].num, *MPI_Point, 0, j, MPI_COMM_WORLD);
				}
				free(low);
				CopyCounter(low, high, SizeGroup[i]);
				free(high);
			}


		}
		else
		{
			if (i > 0)
			{
				GetCounter(poLayer, Order[i], Order[i + 1], low, high, SizeGroup[i - 1], SizeGroup[i], Elevation[i - 1], Elevation[i]);
				CounterAll = Merge(low, high, SizeGroup[i - 1], SizeGroup[i], number);
				//WriteCounter(TheFile, CounterAll, number, Elevation[i], Elevation[i + 1]);
				printf("No.%d succuss \n", MyId);

			}

		}


	}


	if (MyId == 0)
	{
		for (int j = 1; j < MPISize; j++)
		{
			for (int i = 0; i < number; i++)
			{
				MPI_Recv(&CounterSize, 1, MPI_INT, j, i, MPI_COMM_WORLD, Status);
				Arry = new OGRRawPoint[CounterSize];
				MPI_Recv(&Arry, CounterSize, *MPI_Point, j, i, MPI_COMM_WORLD, Status);
				Arry2Chain(Arry, CounterAll[i], CounterSize);
			}
			//WriteCounter(TheFile, CounterAll, number, Elevation[i], Elevation[i + 1]);
			printf("No.%d succuss \n", j);
			free(Arry);
		}
	}

}

void TheInsert2(const char* TheFile, int number, int MyId, int MPISize)
{
	OGRRawPoint **CounterData = NULL;
	Counter **CounterAll = NULL;
	MPI_Datatype *MPI_Point = new MPI_Datatype;
	MPI_Status *Status = new MPI_Status;
	OGRRawPoint *Arry = NULL;
	const char *filename = "./Out_Put/Counter_out_test.shp";
	int Group, SizeofF, TheStart, TheEnd, *SizeGroup = NULL;
	MPI_Type_contiguous(2, MPI_DOUBLE, MPI_Point);
	MPI_Type_commit(MPI_Point);
	GDALAllRegister();
	OGRRegisterAll();

	if (!MyId)
	{
		GDALDataset *poVec = NULL, *poWrite = NULL;
		OGRFeature *poFeature = NULL;
		OGRGeometry *poGeo = NULL;
		OGRLineString *Line = NULL;
		OGRLayer *poLayer = NULL;
		OGRRawPoint *SendData = NULL;
		Counter **Result, testResult;
		int SendSize, RecvSize, counterSize = NULL, *Order = NULL;
		double *Elevation = NULL;

		poVec = (GDALDataset*)GDALOpenEx(filename, GDAL_OF_VECTOR, NULL, NULL, NULL);
		if (poVec == NULL)
		{
			printf("Open failed.\n");
			exit(1);
		}
		poLayer = poVec->GetLayer(0);
		SizeofF = poLayer->GetFeatureCount();
		Group = (SizeofF - 1) / MPISize;
		MPI_Bcast(&Group, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(&SizeofF, 1, MPI_INT, 0, MPI_COMM_WORLD);
		printf("Bcast success of %d", MyId);
		SizeGroup = new int[Group];

		SortingOrder(poLayer, Order, Elevation);

		CounterData = new OGRRawPoint*[Group];
		TheStart = MyId * Group;
		TheEnd = ((MyId + 1) * Group) < (SizeofF - 1) ? ((MyId + 1) * Group) : (SizeofF - 1);

		for (int i = TheStart; i <= TheEnd; i++)
		{
			GetCounter(poLayer, Order[i], CounterData[i], SizeGroup[i]);
		}
		for (int i = 1; i < MPISize; i++)
		{
			TheStart = i * Group;
			TheEnd = ((i + 1) * Group) < (SizeofF - 1) ? ((i + 1) * Group) : (SizeofF - 1);
			int temp = 0;
			for (int j = TheStart; j <= TheEnd; j++)
			{
				GetCounter(poLayer, Order[j], SendData, SendSize);
				MPI_Send(&SendSize, 1, MPI_INT, i, temp, MPI_COMM_WORLD);
				MPI_Send(SendData, SendSize, *MPI_Point, i, temp, MPI_COMM_WORLD);
				temp++;
			}
		}

		TheStart = MyId * Group;
		TheEnd = ((MyId + 1) * Group) < (SizeofF - 1) ? ((MyId + 1) * Group) : (SizeofF - 1);
		CounterAll = new Counter*[TheEnd - TheStart];
		Result = new Counter*[SizeofF - 1 - (TheEnd - TheStart)];


		for (int i = TheStart; i < TheEnd; i++)
		{
			CounterAll[i] = Merge(CounterData[i], CounterData[i + 1], SizeGroup[i], SizeGroup[i + 1], number);
			printf("No.%d has merged for %d of %d.\n", MyId, i, TheEnd - TheStart);
			WriteCounter(TheFile, CounterAll[i], number, Elevation[i], Elevation[i + 1]);

		}

		int temp = 0;

		for (int i = 1; i < MPISize; i++)
		{
			TheStart = i * Group;
			TheEnd = ((i + 1) * Group) < (SizeofF - 1) ? ((i + 1) * Group) : (SizeofF - 1);
			int CounterNum = TheEnd - TheStart;
			for (int j = 0; j < CounterNum; j++)
			{
				for (int k = 0; k < number; k++)
				{
					MPI_Recv(&RecvSize, 1, MPI_INT, i, k + j * number, MPI_COMM_WORLD, Status);
					Arry = new OGRRawPoint[RecvSize];
					MPI_Recv(Arry, RecvSize, *MPI_Point, i, k + j * number, MPI_COMM_WORLD, Status);
					Arry2Chain(Arry, Result[temp][k], RecvSize);
					free(Arry);
				}
				WriteCounter(TheFile, Result[temp], number, Elevation[i], Elevation[i + 1]);
				temp++;

			}

		}


	}
	else
	{
		int RecvSize;
		MPI_Bcast(&Group, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(&SizeofF, 1, MPI_INT, 0, MPI_COMM_WORLD);
		printf("Bcast success of %d", MyId);

		TheStart = MyId * Group;
		TheEnd = ((MyId + 1) * Group) < (SizeofF - 1) ? ((MyId + 1) * Group) : (SizeofF - 1);

		RecvSize = TheEnd - TheStart + 1;

		if (RecvSize > 0)
		{
			SizeGroup = new int[RecvSize];
			CounterData = new OGRRawPoint*[RecvSize];
			CounterAll = new Counter*[RecvSize - 1];

		}

		for (int i = 0; i < RecvSize; i++)
		{
			MPI_Recv(&SizeGroup[i], 1, MPI_INT, 0, i, MPI_COMM_WORLD, Status);

			CounterData[i] = new OGRRawPoint[SizeGroup[i]];
			MPI_Recv(CounterData[i], SizeGroup[i], *MPI_Point, 0, i, MPI_COMM_WORLD, Status);

		}


		for (int i = 0; i < RecvSize - 1; i++)
		{
			CounterAll[i] = Merge(CounterData[i], CounterData[i + 1], SizeGroup[i], SizeGroup[i + 1], number);
			printf("No.%d has merged for %d of %d", MyId, i, TheEnd - TheStart);

		}

		for (int i = 0; i < RecvSize - 1; i++)
		{
			for (int j = 0; j < number; j++)
			{
				Chain2Arry(CounterAll[i][j], Arry);
				MPI_Send(&CounterAll[i][j].num, 1, MPI_INT, 0, j + i * number, MPI_COMM_WORLD);
				MPI_Send(Arry, CounterAll[i][j].num, *MPI_Point, 0, j + i * number, MPI_COMM_WORLD);
				free(Arry);
			}
		}


	}


}


void Bcast(int Id, int Size, int &test, int &test2)
{
	if (Id == 0)
	{
		test = 1;
		test2 = 2;
		MPI_Bcast(&test, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(&test2, 1, MPI_INT, 0, MPI_COMM_WORLD);
	}
	else
	{
		MPI_Bcast(&test, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(&test2, 1, MPI_INT, 0, MPI_COMM_WORLD);
	}
}


int main(int argc, char *argv[])
{
	int rank, numproces, proNum;
	int namelen, test, test2;
	double Time;
	clock_t Start_Time, End_Time;
	MPI_Datatype *MPI_Point = new MPI_Datatype;
	MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &numproces);
	
	
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);//获得进程号
	MPI_Comm_size(MPI_COMM_WORLD, &proNum);//返回通信子的进程数
	
	Start_Time = clock();


	TheInsert2("./Out_Put/Counter_out_test.shp", 2, rank, proNum);


	End_Time = clock();
	Time = End_Time - Start_Time / CLOCKS_PER_SEC;

	MPI_Finalize();
	
	return 0;
	

}
