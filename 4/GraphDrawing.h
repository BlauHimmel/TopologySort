#ifndef GRAPH_DRAWING
#define GRAPH_DRAWING

#include "stdafx.h"
#include "Graph.h"
#include <math.h>
#define RADIUS 12
#define OFFSET 15
#define PI 3.14159265354

/*
	point2：始点
	point1：末点
*/
void DrawArrow(CDC* pdc, CPoint point1, CPoint point2)
{
	double theta = PI / 12;	//箭头两侧直线与所画直线之间的夹角
	double len = 15.0;		//箭头两侧直线的长度
	double theta1;			//所画直线与水平方向之间的夹角

	theta1 = atan2(double(point2.y - point1.y), double(point2.x - point1.x));
	double dOffsetX = abs(OFFSET*cos(theta1));
	double dOffsetY = abs(OFFSET*sin(theta1));
	if (point1.x > point2.x)
	{
		point1.x -= dOffsetX;
		point2.x += dOffsetX;
	}
	else if (point1.x < point2.x)
	{
		point1.x += dOffsetX;
		point2.x -= dOffsetX;
	}

	if (point1.y > point2.y)
	{
		point1.y -= dOffsetY;
		point2.y += dOffsetY;
	}
	else if (point1.y < point2.y)
	{
		point1.y += dOffsetY;
		point2.y -= dOffsetY;
	}
	
	pdc->MoveTo(point1);	//起始点
	pdc->LineTo(point2);	//终点

	pdc->MoveTo(point1);
	pdc->LineTo((int)(point1.x + len * cos(theta + theta1)), (int)(point1.y + len * sin(theta + theta1)));
	pdc->MoveTo(point1);
	pdc->LineTo((int)(point1.x + len * cos(theta1 - theta)), (int)(point1.y + len * sin(theta1 - theta)));
}

void GraphDrawing(CRect& rect,CDC* pDC)
{
	int nHeight = rect.Height();
	int nWidth = rect.Width();
	
	double dRadiusOfGraph = (nHeight - 5 * RADIUS) / 2.0;

	CPoint CCenter;
	CCenter.SetPoint(nWidth / 2, nHeight / 2);
	int nNumOfVertices = CGraph.getNumberOfVertices();
	double dAngle = 360.0 / (double)nNumOfVertices;

	for (int i = 0; i < nNumOfVertices; i++)
	{
		CPoint CPointOfGraph;
		double dAlpha = (360.0 / nNumOfVertices)*(i + 1);
		double temp = dAlpha / 180.0 * PI;;
		double temp1 = sin(dAlpha / 180.0 * PI);
		double temp2 = cos(dAlpha / 180.0 * PI);
		CPointOfGraph.SetPoint(CCenter.x + dRadiusOfGraph*sin(dAlpha / 180.0 * PI), CCenter.y + dRadiusOfGraph*cos(dAlpha / 180.0 * PI));
		pDC->Ellipse(CPointOfGraph.x - RADIUS, CPointOfGraph.y - RADIUS, CPointOfGraph.x + RADIUS, CPointOfGraph.y + RADIUS);
		CString cstrText;
		cstrText.Format(_T("%d"), CGraph.getValue(i));
		pDC->TextOutW(CPointOfGraph.x - 6, CPointOfGraph.y - 6, cstrText);
	}

	for (int i = 0; i < nNumOfVertices; i++)
	{
		int nEdge = CGraph.getFirstNeighbor(i);
		
		if (nEdge != -1)
		{
			CPoint CPointOfStart;
			double dAlphaS = (360.0 / nNumOfVertices)*(i + 1);
			CPointOfStart.SetPoint(CCenter.x + dRadiusOfGraph*sin(dAlphaS / 180.0 * PI), CCenter.y + dRadiusOfGraph*cos(dAlphaS / 180.0 * PI));

			CPoint CPointOfEnd;
			double dAlphaE = (360.0 / nNumOfVertices)*(nEdge + 1);
			CPointOfEnd.SetPoint(CCenter.x + dRadiusOfGraph*sin(dAlphaE / 180.0 * PI), CCenter.y + dRadiusOfGraph*cos(dAlphaE / 180.0 * PI));
			DrawArrow(pDC, CPointOfEnd, CPointOfStart);
			
			int nNextEdge = CGraph.getNextNeighbor(i, nEdge);

			while (nNextEdge != -1)
			{
				dAlphaE = (360.0 / nNumOfVertices)*(nNextEdge + 1);
				CPointOfEnd.SetPoint(CCenter.x + dRadiusOfGraph*sin(dAlphaE / 180.0 * PI), CCenter.y + dRadiusOfGraph*cos(dAlphaE / 180.0 * PI));
				DrawArrow(pDC, CPointOfEnd, CPointOfStart);
				nNextEdge = CGraph.getNextNeighbor(i, nNextEdge);
			}
		}
	}
}

#endif