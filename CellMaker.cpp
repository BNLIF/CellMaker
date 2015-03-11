//
//  CellMaker - A Wire-Cell Generator and Visualization Tool for LArTPC Experiments
//
//    Author:   Michael Mooney (BNL)
//    Created:  March 7th, 2015
//    Updated:  March 10th, 2015
//

#include <iostream>
#include <iomanip>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <string>
#include <vector>

#include <TMath.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TF1.h>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TView.h>
#include <TView3D.h>
#include <TPolyLine3D.h>
#include <TPolyMarker3D.h>
#include <TRandom3.h>
#include <TStyle.h>
#include <TAxis.h>
#include <TGaxis.h>
#include <TLegend.h>
#include <TStopwatch.h>
#include <TVector3.h>
#include <TVirtualFFT.h>
#include <TSystem.h>

using namespace std;

const Double_t heightToWidthRatio = 0.50;
const Double_t wirePitch = 0.30; // in cm
const Double_t firstYWireUoffsetYval = 0.00; // in cm
const Double_t firstYWireVoffsetYval = 0.00; // in cm
const Double_t firstYWireZval = wirePitch/2.0; // in cm
const Double_t leftEdgeOffsetZval = 0.0; // in cm
const Double_t rightEdgeOffsetZval = 0.0; // in cm

const Int_t colorVec[12] = {2,3,4,5,6,7,30,36,38,40,46,48};
//const Int_t colorVec[24] = {30,36,38,40,46,48,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19}; 

const Double_t PI = 3.141592653589793;

Double_t angleU = 60.0;
Double_t angleV = 60.0;
Int_t numYWires = 10;
Int_t makePlot = 0;

struct Cell
{
  Int_t ID;
  pair<Double_t,Double_t> center;
  vector<pair<Double_t,Double_t> > vertices;
  Double_t area;
};

struct YWire
{
  Int_t ID;
  Double_t Zval;
  vector<Cell> cells;
};

vector<YWire> constructCellMap();
YWire constructYWire(Int_t wireID, Double_t wireZval, Double_t YvalOffsetU, Double_t YvalOffsetV);
Cell constructCell(Int_t cellID, Double_t YWireZval, Double_t UWireYval, Double_t VWireYval);
Bool_t formsCell(Double_t UWireYval, Double_t VWireYval);
vector<pair<Double_t,Double_t> > getCellVertices(Double_t YWireZval, Double_t UWireYval, Double_t VWireYval);
pair<Double_t,Double_t> getIntersectionUV(Double_t ZvalOffset, Double_t slope1, Double_t intercept1, Double_t slope2, Double_t intercept2);
pair<Double_t, Double_t> getIntersectionY(Double_t ZvalOffset, Double_t YWireZval, Double_t slope, Double_t intercept);
vector<pair<Double_t,Double_t> > sortVertices(vector<pair<Double_t,Double_t> > vertices);
vector<pair<Double_t,Double_t> > trimEdgeCellVertices(vector<pair<Double_t,Double_t> > vertices, Int_t edgeType, Double_t edgeVal);
Bool_t vertexOutsideBoundary(pair<Double_t,Double_t> vertex, Int_t edgeType, Double_t edgeVal);
pair<Double_t,Double_t> calcCellCenter(vector<pair<Double_t,Double_t> > vertices);
Double_t calcCellArea(vector<pair<Double_t,Double_t> > vertices);
Bool_t compareOrientedVertices(pair<Double_t,pair<Double_t,Double_t> > orientedVertex1, pair<Double_t,pair<Double_t,Double_t> > orientedVertex2);
void drawCellMap(vector<YWire> wires, Int_t numWires, Int_t numCells);

/////////////////////////////////////////////////////////////////////////////////////////////////////
// main - Main function to run program
/////////////////////////////////////////////////////////////////////////////////////////////////////
Int_t main(Int_t argc, Char_t** argv)
{
  // Setup environment
  gErrorIgnoreLevel = kError;

  // Get input parameters
  if(argc > 1)
    angleU = (Double_t) atof(argv[1]);
  if(argc > 2)
    angleV = (Double_t) atof(argv[2]);
  if(argc > 3)
    numYWires = (Int_t) atoi(argv[3]);
  if(argc > 4)
    makePlot = (Int_t) atoi(argv[4]);

  // Create and draw cell map
  vector<YWire> wires = constructCellMap();
  if(makePlot == 1)
    drawCellMap(wires,-1,-1);

  // End of program
  return 0;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
// constructCellMap - Construct map of Cells by creating array of YWires (each with associated Cells)
/////////////////////////////////////////////////////////////////////////////////////////////////////
vector<YWire> constructCellMap()
{ 
  vector<YWire> wires;

  const Double_t UdeltaY = wirePitch/tan((PI/180.0)*angleU);
  const Double_t VdeltaY = -1.0*wirePitch/tan((PI/180.0)*angleV);
  const Double_t UspacingOnWire = fabs(wirePitch/sin((PI/180.0)*angleU));
  const Double_t VspacingOnWire = fabs(wirePitch/sin((PI/180.0)*angleV));

  Double_t Zval = firstYWireZval;
  Double_t Uoffset = firstYWireUoffsetYval;
  Double_t Voffset = firstYWireVoffsetYval;

  while(Uoffset > ((UspacingOnWire-UdeltaY)/2.0)+0.0000000001)
    Uoffset -= UspacingOnWire;
  while(Voffset > ((VspacingOnWire+VdeltaY)/2.0)+0.0000000001)
    Voffset -= VspacingOnWire;

  for(Int_t i = 0; i < numYWires; i++)  
  {
    wires.push_back(constructYWire(i,Zval,Uoffset,Voffset));

    Zval += wirePitch;
    Uoffset += UdeltaY;
    Voffset += VdeltaY;

    while(Uoffset > ((UspacingOnWire-UdeltaY)/2.0)+0.0000000001)
      Uoffset -= UspacingOnWire;
    while(Voffset < ((VdeltaY-VspacingOnWire)/2.0)-0.0000000001)
      Voffset += VspacingOnWire;
  }

  return wires;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
// constructYWire - Construct YWire Cells by forming pairs of crossings from a U wire and a V wire
/////////////////////////////////////////////////////////////////////////////////////////////////////
YWire constructYWire(Int_t wireID, Double_t wireZval, Double_t YvalOffsetU, Double_t YvalOffsetV)
{ 
  YWire wire;

  wire.ID = wireID;
  wire.Zval = wireZval;
 
  const Double_t maxHeight = heightToWidthRatio*wirePitch*numYWires; 
  const Double_t UdeltaY = wirePitch/tan((PI/180.0)*angleU);
  const Double_t VdeltaY = -1.0*wirePitch/tan((PI/180.0)*angleV);
  const Double_t UspacingOnWire = fabs(wirePitch/sin((PI/180.0)*angleU));
  const Double_t VspacingOnWire = fabs(wirePitch/sin((PI/180.0)*angleV));

  Int_t numUcrosses = TMath::Ceil((maxHeight+(UdeltaY-UspacingOnWire)/2.0-YvalOffsetU)/UspacingOnWire)+1;
  Int_t numVcrosses = TMath::Ceil((maxHeight-(VdeltaY+VspacingOnWire)/2.0-YvalOffsetV)/VspacingOnWire)+1;

  Int_t cellID = 0;
  Cell cell;
  Bool_t flag1;
  Bool_t flag2;
  Int_t j;
  for(Int_t i = 0; i < numUcrosses; i++)
  {
    flag1 = false;
    flag2 = false;
    j = 0;
    while((j < numVcrosses) && (flag2 == false))
    {
      if(formsCell(YvalOffsetU+i*UspacingOnWire,YvalOffsetV+j*VspacingOnWire) == true)
      {
        flag1 = true;
        cell = constructCell(cellID,wireZval,YvalOffsetU+i*UspacingOnWire,YvalOffsetV+j*VspacingOnWire);

        if(cell.vertices.size() > 2)
	{
          wire.cells.push_back(cell);
          cellID++;
	}
      }
      else if(flag1 == true)
        flag2 = true;

      j++;
    }
  }
 
  return wire;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
// constructCell - Construct one Cell on a YWire, including vertices, area, and center point
/////////////////////////////////////////////////////////////////////////////////////////////////////
Cell constructCell(Int_t cellID, Double_t YWireZval, Double_t UWireYval, Double_t VWireYval)
{
  Cell cell;
  cell.ID = cellID;

  cell.vertices = getCellVertices(YWireZval,UWireYval,VWireYval);
  cell.center = calcCellCenter(cell.vertices);
  cell.area = calcCellArea(cell.vertices);

  return cell;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
// formsCell - Check whether or not a U/V crossing pair on a particular YWire forms a Cell
/////////////////////////////////////////////////////////////////////////////////////////////////////
Bool_t formsCell(Double_t UWireYval, Double_t VWireYval)
{
  Bool_t isCell;

  const Double_t UspacingOnWire = fabs(wirePitch/sin((PI/180.0)*angleU));
  const Double_t VspacingOnWire = fabs(wirePitch/sin((PI/180.0)*angleV));
  const Double_t tanU = tan((PI/180.0)*angleU);
  const Double_t tanV = tan((PI/180.0)*angleV);
  const Double_t maxHeight = heightToWidthRatio*wirePitch*numYWires; 

  Double_t deltaY;
  if(UWireYval > VWireYval)
    deltaY = (UWireYval-UspacingOnWire/2.0)-(VWireYval+VspacingOnWire/2.0);
  else
    deltaY = (VWireYval-VspacingOnWire/2.0)-(UWireYval+UspacingOnWire/2.0);

  if((UWireYval+UspacingOnWire/2.0 < 0.0) && (VWireYval+VspacingOnWire/2.0 < 0.0))
    isCell = false;
  else if((UWireYval-UspacingOnWire/2.0 > maxHeight) && (VWireYval-VspacingOnWire/2.0 > maxHeight))
    isCell = false;
  else if(deltaY < 0.0000000001)
    isCell = true;
  else if(UWireYval == VWireYval)
    isCell = true;
  else if(fabs((tanU*tanV*deltaY)/(tanU+tanV)) < wirePitch/2.0)
    isCell = true;
  else
    isCell = false;

  return isCell;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
// getCellVertices - Find vertices that define the boundaries of a Cell
/////////////////////////////////////////////////////////////////////////////////////////////////////
vector<pair<Double_t,Double_t> > getCellVertices(Double_t YWireZval, Double_t UWireYval, Double_t VWireYval)
{
  vector<pair<Double_t,Double_t> > cellVertices;

  const Double_t maxHeight = heightToWidthRatio*wirePitch*numYWires;
  const Double_t UspacingOnWire = fabs(wirePitch/sin((PI/180.0)*angleU));
  const Double_t VspacingOnWire = fabs(wirePitch/sin((PI/180.0)*angleV));

  const Double_t U1slope = 1.0/tan((PI/180.0)*angleU);
  const Double_t U2slope = U1slope;
  const Double_t U1intercept = UWireYval - UspacingOnWire/2.0;
  const Double_t U2intercept = UWireYval + UspacingOnWire/2.0;

  const Double_t V1slope = -1.0/tan((PI/180.0)*angleV);
  const Double_t V2slope = V1slope;
  const Double_t V1intercept = VWireYval - VspacingOnWire/2.0;
  const Double_t V2intercept = VWireYval + VspacingOnWire/2.0;

  const Double_t Y1Zval = YWireZval - wirePitch/2.0;
  const Double_t Y2Zval = YWireZval + wirePitch/2.0;

  pair<Double_t,Double_t> U1V1intersectionPoint;
  pair<Double_t,Double_t> U1V2intersectionPoint;
  pair<Double_t,Double_t> U2V1intersectionPoint;
  pair<Double_t,Double_t> U2V2intersectionPoint;

  U1V1intersectionPoint = getIntersectionUV(YWireZval,U1slope,U1intercept,V1slope,V1intercept);
  U1V2intersectionPoint = getIntersectionUV(YWireZval,U1slope,U1intercept,V2slope,V2intercept);
  U2V1intersectionPoint = getIntersectionUV(YWireZval,U2slope,U2intercept,V1slope,V1intercept);
  U2V2intersectionPoint = getIntersectionUV(YWireZval,U2slope,U2intercept,V2slope,V2intercept);

  Double_t UVminZval = U1V1intersectionPoint.first;
  Double_t UVmaxZval = U1V1intersectionPoint.first;
  Double_t UVminYval = U1V1intersectionPoint.second;
  Double_t UVmaxYval = U1V1intersectionPoint.second;

  if(U1V2intersectionPoint.first < UVminZval)
    UVminZval = U1V2intersectionPoint.first;
  if(U1V2intersectionPoint.first > UVmaxZval)
    UVmaxZval = U1V2intersectionPoint.first;
  if(U1V2intersectionPoint.second < UVminYval)
    UVminYval = U1V2intersectionPoint.second;
  if(U1V2intersectionPoint.second > UVmaxYval)
    UVmaxYval = U1V2intersectionPoint.second;

  if(U2V1intersectionPoint.first < UVminZval)
    UVminZval = U2V1intersectionPoint.first;
  if(U2V1intersectionPoint.first > UVmaxZval)
    UVmaxZval = U2V1intersectionPoint.first;
  if(U2V1intersectionPoint.second < UVminYval)
    UVminYval = U2V1intersectionPoint.second;
  if(U2V1intersectionPoint.second > UVmaxYval)
    UVmaxYval = U2V1intersectionPoint.second;

  if(U2V2intersectionPoint.first < UVminZval)
    UVminZval = U2V2intersectionPoint.first;
  if(U2V2intersectionPoint.first > UVmaxZval)
    UVmaxZval = U2V2intersectionPoint.first;
  if(U2V2intersectionPoint.second < UVminYval)
    UVminYval = U2V2intersectionPoint.second;
  if(U2V2intersectionPoint.second > UVmaxYval)
    UVmaxYval = U2V2intersectionPoint.second;

  pair<Double_t,Double_t> Y1U1intersectionPoint;
  pair<Double_t,Double_t> Y1U2intersectionPoint;
  pair<Double_t,Double_t> Y1V1intersectionPoint;
  pair<Double_t,Double_t> Y1V2intersectionPoint;
  pair<Double_t,Double_t> Y2U1intersectionPoint;
  pair<Double_t,Double_t> Y2U2intersectionPoint;
  pair<Double_t,Double_t> Y2V1intersectionPoint;
  pair<Double_t,Double_t> Y2V2intersectionPoint;

  Y1U1intersectionPoint = getIntersectionY(YWireZval,Y1Zval,U1slope,U1intercept);
  Y1U2intersectionPoint = getIntersectionY(YWireZval,Y1Zval,U2slope,U2intercept);
  Y1V1intersectionPoint = getIntersectionY(YWireZval,Y1Zval,V1slope,V1intercept);
  Y1V2intersectionPoint = getIntersectionY(YWireZval,Y1Zval,V2slope,V2intercept);
  Y2U1intersectionPoint = getIntersectionY(YWireZval,Y2Zval,U1slope,U1intercept);
  Y2U2intersectionPoint = getIntersectionY(YWireZval,Y2Zval,U2slope,U2intercept);
  Y2V1intersectionPoint = getIntersectionY(YWireZval,Y2Zval,V1slope,V1intercept);
  Y2V2intersectionPoint = getIntersectionY(YWireZval,Y2Zval,V2slope,V2intercept);

  if((Y1U1intersectionPoint.second < UVmaxYval) && (Y1U1intersectionPoint.second > UVminYval) && (Y1U1intersectionPoint.first < UVmaxZval) && (Y1U1intersectionPoint.first > UVminZval))
    cellVertices.push_back(Y1U1intersectionPoint);
  if((Y1U2intersectionPoint.second < UVmaxYval) && (Y1U2intersectionPoint.second > UVminYval) && (Y1U2intersectionPoint.first < UVmaxZval) && (Y1U2intersectionPoint.first > UVminZval))
    cellVertices.push_back(Y1U2intersectionPoint);
  if((Y1V1intersectionPoint.second < UVmaxYval) && (Y1V1intersectionPoint.second > UVminYval) && (Y1V1intersectionPoint.first < UVmaxZval) && (Y1V1intersectionPoint.first > UVminZval))
    cellVertices.push_back(Y1V1intersectionPoint);
  if((Y1V2intersectionPoint.second < UVmaxYval) && (Y1V2intersectionPoint.second > UVminYval) && (Y1V2intersectionPoint.first < UVmaxZval) && (Y1V2intersectionPoint.first > UVminZval))
    cellVertices.push_back(Y1V2intersectionPoint);

  if((Y2U1intersectionPoint.second < UVmaxYval) && (Y2U1intersectionPoint.second > UVminYval) && (Y2U1intersectionPoint.first < UVmaxZval) && (Y2U1intersectionPoint.first > UVminZval))
    cellVertices.push_back(Y2U1intersectionPoint);
  if((Y2U2intersectionPoint.second < UVmaxYval) && (Y2U2intersectionPoint.second > UVminYval) && (Y2U2intersectionPoint.first < UVmaxZval) && (Y2U2intersectionPoint.first > UVminZval))
    cellVertices.push_back(Y2U2intersectionPoint);
  if((Y2V1intersectionPoint.second < UVmaxYval) && (Y2V1intersectionPoint.second > UVminYval) && (Y2V1intersectionPoint.first < UVmaxZval) && (Y2V1intersectionPoint.first > UVminZval))
    cellVertices.push_back(Y2V1intersectionPoint);
  if((Y2V2intersectionPoint.second < UVmaxYval) && (Y2V2intersectionPoint.second > UVminYval) && (Y2V2intersectionPoint.first < UVmaxZval) && (Y2V2intersectionPoint.first > UVminZval))
    cellVertices.push_back(Y2V2intersectionPoint);
  
  if((U1V1intersectionPoint.first >= Y1Zval) && (U1V1intersectionPoint.first <= Y2Zval))
    cellVertices.push_back(U1V1intersectionPoint);
  if((U1V2intersectionPoint.first >= Y1Zval) && (U1V2intersectionPoint.first <= Y2Zval))
    cellVertices.push_back(U1V2intersectionPoint);
  if((U2V1intersectionPoint.first >= Y1Zval) && (U2V1intersectionPoint.first <= Y2Zval))
    cellVertices.push_back(U2V1intersectionPoint);
  if((U2V2intersectionPoint.first >= Y1Zval) && (U2V2intersectionPoint.first <= Y2Zval))
    cellVertices.push_back(U2V2intersectionPoint);

  vector<pair<Double_t,Double_t> > cellVerticesSorted;
  cellVerticesSorted = sortVertices(cellVertices);

  if(UVminZval < (firstYWireZval-0.5*wirePitch+leftEdgeOffsetZval))
    cellVerticesSorted = trimEdgeCellVertices(cellVerticesSorted,1,firstYWireZval-0.5*wirePitch+leftEdgeOffsetZval);
  else if(UVmaxZval > (firstYWireZval+(numYWires-0.5)*wirePitch-rightEdgeOffsetZval))
    cellVerticesSorted = trimEdgeCellVertices(cellVerticesSorted,2,(firstYWireZval+(numYWires-0.5)*wirePitch-rightEdgeOffsetZval));

  if(UVminYval < 0.0)
    cellVerticesSorted = trimEdgeCellVertices(cellVerticesSorted,3,0.0);
  else if(UVmaxYval > maxHeight)
    cellVerticesSorted = trimEdgeCellVertices(cellVerticesSorted,4,maxHeight);

  return cellVerticesSorted;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
// getIntersectionUV - Get intersection point of a U wire and a V wire
/////////////////////////////////////////////////////////////////////////////////////////////////////
pair<Double_t,Double_t> getIntersectionUV(Double_t ZvalOffset, Double_t slope1, Double_t intercept1, Double_t slope2, Double_t intercept2)
{
  Double_t Yval = (slope2*intercept1-slope1*intercept2)/(slope2-slope1);
  Double_t Zval = ZvalOffset + (Yval-intercept1)/slope1;

  pair<Double_t,Double_t> intersectionPoint;
  intersectionPoint.first = Zval;
  intersectionPoint.second = Yval;

  return intersectionPoint;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
// getIntersectionY - Get intersection point of a YWire and either a U or a V wire
/////////////////////////////////////////////////////////////////////////////////////////////////////
pair<Double_t, Double_t> getIntersectionY(Double_t ZvalOffset, Double_t YWireZval, Double_t slope, Double_t intercept)
{
  Double_t Yval = slope*(YWireZval-ZvalOffset) + intercept;
  Double_t Zval = YWireZval;

  pair<Double_t,Double_t> intersectionPoint;
  intersectionPoint.first = Zval;
  intersectionPoint.second = Yval;

  return intersectionPoint;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
// sortVertices - Sort the vertices of a particular Cell by angle, clockwise
/////////////////////////////////////////////////////////////////////////////////////////////////////
vector<pair<Double_t,Double_t> > sortVertices(vector<pair<Double_t,Double_t> > vertices)
{
  vector<pair<Double_t,Double_t> > verticesSorted;

  if(vertices.size() < 2)
    return vertices;

  pair<Double_t,Double_t> vertex;
  pair<Double_t,Double_t> center = calcCellCenter(vertices);
  pair<Double_t,pair<Double_t,Double_t> > orientedVertex;
  vector<pair<Double_t,pair<Double_t,Double_t> > > orientedVertices;
  for(Int_t i = 0; i < vertices.size(); i++)
  {
    vertex = vertices.at(i);
    orientedVertex.first = atan2(vertex.second-center.second,vertex.first-center.first);
    orientedVertex.second = vertex;
    orientedVertices.push_back(orientedVertex);
  }

  sort(orientedVertices.begin(),orientedVertices.end(),compareOrientedVertices);

  for(Int_t i = 0; i < orientedVertices.size(); i++)
    verticesSorted.push_back(orientedVertices.at(i).second);

  return verticesSorted;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
// trimEdgeCellVertices - Modify Cell at edge of map to be fully contained within map boundaries
/////////////////////////////////////////////////////////////////////////////////////////////////////
vector<pair<Double_t,Double_t> > trimEdgeCellVertices(vector<pair<Double_t,Double_t> > vertices, Int_t edgeType, Double_t edgeVal)
{
  vector<pair<Double_t,Double_t> > edgeCellVertices;

  if((edgeType < 1) && (edgeType > 4))
    return vertices;
  
  Int_t numVertOutside = 0;
  for(Int_t i = 0; i < vertices.size(); i++)
    if(vertexOutsideBoundary(vertices.at(i),edgeType,edgeVal) == true)
      numVertOutside++;

  if(numVertOutside == 0)
    return vertices;

  Int_t j;
  Double_t Zval1;
  Double_t Zval2;
  Double_t Yval1;
  Double_t Yval2;
  Double_t slope;
  Double_t intercept;
  pair<Double_t,Double_t> vertex;
  for(Int_t i = 0; i < vertices.size(); i++)
  {
    j = (i+1) % vertices.size();
    
    Zval1 = vertices.at(i).first;
    Zval2 = vertices.at(j).first;
    Yval1 = vertices.at(i).second;
    Yval2 = vertices.at(j).second;

    if(Zval1 != Zval2)
    {
      slope = (Yval2-Yval1)/(Zval2-Zval1);
      intercept = Yval1-slope*Zval1;

      if((edgeType == 1) || (edgeType == 2))
      {
        if(((edgeVal > Zval1) && (edgeVal < Zval2)) || ((edgeVal > Zval2) && (edgeVal < Zval1)))
	{
          vertex.first = edgeVal;
          vertex.second = slope*edgeVal+intercept;
          edgeCellVertices.push_back(vertex);
	}
      }
      else if((edgeType == 3) || (edgeType == 4))
      {
        if(((edgeVal > Yval1) && (edgeVal < Yval2)) || ((edgeVal > Yval2) && (edgeVal < Yval1)))
	{
          vertex.first = (edgeVal-intercept)/slope;
          vertex.second = edgeVal;
          edgeCellVertices.push_back(vertex);
	}
      }
    }
    else if((edgeType == 3) || (edgeType == 4))
    {
      if(((edgeVal > Yval1) && (edgeVal < Yval2)) || ((edgeVal > Yval2) && (edgeVal < Yval1)))
      {
        vertex.first = Zval1;
        vertex.second = edgeVal;
        edgeCellVertices.push_back(vertex);
      }
    }
  }

  for(Int_t i = 0; i < vertices.size(); i++)
    if(vertexOutsideBoundary(vertices.at(i),edgeType,edgeVal) == false)
      edgeCellVertices.push_back(vertices.at(i));

  vector<pair<Double_t,Double_t> > edgeCellVerticesSorted;
  edgeCellVerticesSorted = sortVertices(edgeCellVertices);

  return edgeCellVerticesSorted;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
// vertexOutsideBoundary - Check whether or not vertex of a Cell is contained within map boundaries
/////////////////////////////////////////////////////////////////////////////////////////////////////
Bool_t vertexOutsideBoundary(pair<Double_t,Double_t> vertex, Int_t edgeType, Double_t edgeVal)
{
  switch(edgeType)
  {
    case 1:
      return (vertex.first < edgeVal);
    case 2:
      return (vertex.first > edgeVal);
    case 3:
      return (vertex.second < edgeVal);
    case 4:
      return (vertex.second > edgeVal);
    default:
      return false;
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
// calcCellCenter - Calculate center coordinates of a Cell (units of cm)
/////////////////////////////////////////////////////////////////////////////////////////////////////
pair<Double_t,Double_t> calcCellCenter(vector<pair<Double_t,Double_t> > vertices)
{
  Double_t Zval = 0.0;
  Double_t Yval = 0.0;

  const Int_t numVertices = vertices.size();

  for(Int_t i = 0; i < numVertices; i++)
  {
    Zval += vertices.at(i).first;
    Yval += vertices.at(i).second;
  }  

  Zval /= ((Double_t) numVertices);
  Yval /= ((Double_t) numVertices);

  pair<Double_t,Double_t> centerPoint;
  centerPoint.first = Zval;
  centerPoint.second = Yval;

  return centerPoint;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
// calcCellArea - Calculate area of a cell (units of cm^2)
/////////////////////////////////////////////////////////////////////////////////////////////////////
Double_t calcCellArea(vector<pair<Double_t,Double_t> > vertices)
{ 
  Double_t cellArea = 0.0;
  const Int_t numVertices = vertices.size();

  Int_t j = numVertices-1;
  for(Int_t i = 0; i < numVertices; i++)
  {
    cellArea += (vertices.at(j).first+vertices.at(i).first)*(vertices.at(j).second-vertices.at(i).second); 
    j = i;
  }
  cellArea /= 2.0;

  return cellArea;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
// compareOrientedVertices - Function used when sorting Cell vertices by angle, clockwise
/////////////////////////////////////////////////////////////////////////////////////////////////////
Bool_t compareOrientedVertices(pair<Double_t,pair<Double_t,Double_t> > orientedVertex1, pair<Double_t,pair<Double_t,Double_t> > orientedVertex2)
{
  return (orientedVertex1.first > orientedVertex2.first);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
// drawCellMap - Draw complete Cell map including all YWires, or portion thereof
/////////////////////////////////////////////////////////////////////////////////////////////////////
void drawCellMap(vector<YWire> wires, Int_t numWires, Int_t numCells)
{
  const Double_t maxHeight = heightToWidthRatio*wirePitch*numYWires;
  const Int_t colorFactor = ((sizeof(colorVec)/sizeof(*colorVec))/2);

  numWires = min(numWires,(Int_t)wires.size());
  if(numWires < 0)
    numWires = wires.size();

  TCanvas *c1 = new TCanvas("c1","c1",1600.0,800.0*(heightToWidthRatio/0.5));

  TGraph *cellGraph;
  TMultiGraph *cellMapGraph = new TMultiGraph();

  Int_t colorNum;
  Int_t testNcell;
  Int_t testNvert;
  Double_t testX[7];
  Double_t testY[7];
  for(Int_t i = 0; i < numWires; i++)
  {
    testNcell = min(numCells,(Int_t)wires.at(i).cells.size());
    if(testNcell < 0)
      testNcell = wires.at(i).cells.size();

    for(Int_t j = 0; j < testNcell; j++)
    {
      testNvert = wires.at(i).cells.at(j).vertices.size();
      for(Int_t k = 0; k < testNvert; k++)
      {
        testX[k] = wires.at(i).cells.at(j).vertices.at(k).first;
        testY[k] = wires.at(i).cells.at(j).vertices.at(k).second;
      }
      colorNum = colorVec[colorFactor*(wires.at(i).ID % 2)+(wires.at(i).cells.at(j).ID % colorFactor)];

      cellGraph = new TGraph(testNvert,testX,testY);
      cellGraph->SetTitle("");
      cellGraph->SetLineColor(colorNum);
      cellGraph->SetFillColor(colorNum);
      cellMapGraph->Add(cellGraph);
    }
  }
  cellMapGraph->Draw("AF");
  cellMapGraph->SetTitle(Form("#theta_{U} = %.1f#circ,_{} #theta_{V} = %.1f#circ,_{} N_{wires} = %d",angleU,angleV,numYWires));
  cellMapGraph->GetXaxis()->SetTitle("Z [cm]");
  cellMapGraph->GetXaxis()->SetTitleOffset(1.0);
  cellMapGraph->GetXaxis()->SetTitleSize(0.04);
  cellMapGraph->GetYaxis()->SetTitle("Y [cm]");
  cellMapGraph->GetYaxis()->SetTitleOffset(0.8*(heightToWidthRatio/0.5));
  cellMapGraph->GetYaxis()->SetTitleSize(0.04);
  cellMapGraph->GetXaxis()->SetLimits(firstYWireZval-0.5*wirePitch+leftEdgeOffsetZval,firstYWireZval+(numYWires-0.5)*wirePitch-rightEdgeOffsetZval);
  cellMapGraph->GetHistogram()->SetMinimum(0.0);
  cellMapGraph->GetHistogram()->SetMaximum(maxHeight);

  gPad->Update();
  gPad->RedrawAxis();
  TLine borderLine;
  borderLine.DrawLine(gPad->GetUxmin(), gPad->GetUymax(), gPad->GetUxmax(), gPad->GetUymax());
  borderLine.DrawLine(gPad->GetUxmax(), gPad->GetUymin(), gPad->GetUxmax(), gPad->GetUymax());

  c1->SaveAs(Form("cellDiagram_%dUAngle_%dVAngle_%dYWires.png",(Int_t) round(angleU),(Int_t) round(angleV),numYWires));

  return;
}
