//
//  CellMaker - A Wire-Cell Generator and Visualization Tool for LArTPC Experiments
//
//    Author:   Michael Mooney (BNL)
//    Created:  March 7th, 2015
//    Updated:  March 18th, 2015
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
const Double_t wirePitchY = 0.30; // in cm
const Double_t wirePitchU = 0.30; // in cm
const Double_t wirePitchV = 0.30; // in cm
Double_t firstYwireUoffsetYval = 0.00; // in cm  // TEMP make non const so can modify after obtaining numYwires
Double_t firstYwireVoffsetYval = 0.00; // in cm  // TEMP make non const so can modify after obtaining numYwires
const Double_t firstYwireZval = wirePitchY/2.0; // in cm
const Double_t leftEdgeOffsetZval = 0.0; // in cm
const Double_t rightEdgeOffsetZval = 0.0; // in cm

const Int_t colorVec[12] = {2,3,4,5,6,7,30,36,38,40,46,48};
//const Int_t colorVec[24] = {30,36,38,40,46,48,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19}; 
const Int_t colorVec2[8] = {2,3,4,5,6,7,8,9};
//const Int_t colorVec3[12] = {22,31,32,14,28,34,35,36,24,37,39,40};
const Int_t colorVec3[12] = {11,14,35,31,24,16,22,13,34,27,15,25};

const Double_t PI = 3.141592653589793;
const Double_t epsilon = 0.0000000001;

Double_t angleU = 60.0;
Double_t angleV = 60.0;
Int_t numYwires = 10;
Int_t plotMode = 0;

enum Plane_t {kUPlane, kVPlane, kYPlane};
enum Hit_t {kNoHit, kRealHit, kFakeHit};

struct Cell
{
  Int_t ID;
  vector<pair<Double_t,Double_t> > vertices;
  pair<Double_t,Double_t> center;
  Double_t area;
  Double_t trueCharge;
  Double_t recoCharge;
  Int_t UwireID;
  Int_t VwireID;
  Int_t YwireID;
  Hit_t hitType;
};

struct Wire
{
  Int_t ID;
  Plane_t plane;
  Double_t location;
  Double_t charge;
  vector<Int_t> cellIDs;
};

struct CellMap
{
  vector<Wire> Uwires;
  vector<Wire> Vwires;
  vector<Wire> Ywires;
  vector<Cell> cells;
};

CellMap constructCellMap();
vector<Wire> constructWires(Plane_t wirePlane);
Wire constructWire(Plane_t wirePlane, Int_t wireID, Int_t wireLocation);
vector<Cell> constructCells(vector<Wire> &Uwires, vector<Wire> &Vwires, vector<Wire> &Ywires);
vector<Cell> constructCellChain(Int_t firstCellID, Double_t wireZval, Double_t YvalOffsetU, Double_t YvalOffsetV);
Cell constructCell(Int_t cellID, Double_t YwireZval, Double_t UwireYval, Double_t VwireYval);
Bool_t formsCell(Double_t UwireYval, Double_t VwireYval);
vector<pair<Double_t,Double_t> > getCellVertices(Double_t YwireZval, Double_t UwireYval, Double_t VwireYval);
pair<Double_t,Double_t> getIntersectionUV(Double_t ZvalOffset, Double_t slope1, Double_t intercept1, Double_t slope2, Double_t intercept2);
pair<Double_t, Double_t> getIntersectionY(Double_t ZvalOffset, Double_t YwireZval, Double_t slope, Double_t intercept);
vector<pair<Double_t,Double_t> > sortVertices(vector<pair<Double_t,Double_t> > vertices);
vector<pair<Double_t,Double_t> > trimEdgeCellVertices(vector<pair<Double_t,Double_t> > vertices, Int_t edgeType, Double_t edgeVal);
Bool_t vertexOutsideBoundary(pair<Double_t,Double_t> vertex, Int_t edgeType, Double_t edgeVal);
pair<Double_t,Double_t> calcCellCenter(vector<pair<Double_t,Double_t> > vertices);
Double_t calcCellArea(vector<pair<Double_t,Double_t> > vertices);
Bool_t compareOrientedVertices(pair<Double_t,pair<Double_t,Double_t> > orientedVertex1, pair<Double_t,pair<Double_t,Double_t> > orientedVertex2);
Double_t getUwireYval(Int_t IDnum, Double_t Zval);
Double_t getUwireZval(Int_t IDnum, Double_t Yval);
Int_t getUwireID(Double_t Yval, Double_t Zval);
Double_t getVwireYval(Int_t IDnum, Double_t Zval);
Double_t getVwireZval(Int_t IDnum, Double_t Yval);
Int_t getVwireID(Double_t Yval, Double_t Zval);
Double_t getYwireZval(Int_t IDnum);
Int_t getYwireID(Double_t Zval);
pair<pair<Double_t,Double_t>,pair<Double_t,Double_t> > getWireEndpoints(Int_t wireID, Plane_t wirePlane);
void addCharges(CellMap &cellMap);
void assignHitTypes(CellMap &cellMap);
void drawCellMap(CellMap const& cellMap, Int_t numWires, Int_t numCells);

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
    numYwires = (Int_t) atoi(argv[3]);
  if(argc > 4)
    plotMode = (Int_t) atoi(argv[4]);

  // Adjust offsets for MicroBooNE case
  const Double_t maxHeight = heightToWidthRatio*wirePitchY*numYwires; 
  const Double_t UspacingOnWire = fabs(wirePitchU/sin((PI/180.0)*angleU));
  const Double_t VspacingOnWire = fabs(wirePitchV/sin((PI/180.0)*angleV));
  Double_t tempUoffset = firstYwireUoffsetYval;
  while(tempUoffset > UspacingOnWire-epsilon)
    tempUoffset -= UspacingOnWire;
  firstYwireVoffsetYval = maxHeight-tempUoffset-TMath::Floor((maxHeight-tempUoffset)/UspacingOnWire)*UspacingOnWire;
  while(firstYwireVoffsetYval > VspacingOnWire-epsilon)
    firstYwireVoffsetYval -= VspacingOnWire;

  // Create and draw cell map
  CellMap globalCellMap = constructCellMap();
  addCharges(globalCellMap);
  if(plotMode > 0)
    drawCellMap(globalCellMap,-1,-1);

  // End of program
  return 0;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
// constructCellMap - Construct map of Cells formed by three Wires (one from each plane)
/////////////////////////////////////////////////////////////////////////////////////////////////////
CellMap constructCellMap()
{ 
  CellMap cellMap;

  cellMap.Uwires = constructWires(kUPlane);
  cellMap.Vwires = constructWires(kVPlane);
  cellMap.Ywires = constructWires(kYPlane);

  cellMap.cells = constructCells(cellMap.Uwires,cellMap.Vwires,cellMap.Ywires);

  return cellMap;
}
 
/////////////////////////////////////////////////////////////////////////////////////////////////////
// constructWires - Construct set of Wires for a particular plane
/////////////////////////////////////////////////////////////////////////////////////////////////////
vector<Wire> constructWires(Plane_t wirePlane)
{
  vector<Wire> wires;

  const Double_t maxZ = (numYwires-1)*wirePitchY;
  const Double_t maxY = (numYwires-1)*wirePitchY*heightToWidthRatio;
  const Double_t diagLength = sqrt(pow(maxZ,2) + pow(maxY,2));
  const Double_t diagAngle = (180.0/PI)*atan(1.0/heightToWidthRatio);

  Double_t offset = 0.0;
  Int_t maxNum;
  Double_t wirePitch;
  Double_t wireLocation;
  if(wirePlane == kUPlane)
  {
    offset = firstYwireUoffsetYval*(sin((PI/180.0)*angleU)/sin((PI/180.0)*(180.0-diagAngle-angleU)));
    maxNum = TMath::Floor(((diagLength-offset)*sin((PI/180.0)*(diagAngle+angleU)))/wirePitchU);

    wirePitch = wirePitchU;
    wireLocation = firstYwireUoffsetYval*sin((PI/180.0)*angleU);
  }
  else if(wirePlane == kVPlane)
  {
    offset = firstYwireVoffsetYval*(sin((PI/180.0)*angleV)/sin((PI/180.0)*(180.0-diagAngle-angleV)));
    maxNum = TMath::Floor(((diagLength-offset)*sin((PI/180.0)*(diagAngle+angleV)))/wirePitchV);

    wirePitch = wirePitchU;
    wireLocation = firstYwireUoffsetYval*sin((PI/180.0)*angleU);
  }
  else if(wirePlane == kYPlane)
  {
    maxNum = numYwires;
   
    wirePitch = wirePitchY;
    wireLocation = firstYwireZval;
  }
  else
    return wires;

  for(Int_t wireID = 0; wireID < maxNum; wireID++)
  {
    wires.push_back(constructWire(wirePlane,wireID,wireLocation));
    wireLocation += wirePitch;
  }

  return wires;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
// constructWire - Construct one Wire, including location and plane
/////////////////////////////////////////////////////////////////////////////////////////////////////
Wire constructWire(Plane_t wirePlane, Int_t wireID, Int_t wireLocation)
{
  Wire wire;

  wire.ID = wireID;
  wire.plane = wirePlane;
  wire.location = wireLocation;

  wire.charge = 0.0;

  return wire;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
// constructCells - Use collections of Wires to create set of all Cells
/////////////////////////////////////////////////////////////////////////////////////////////////////
vector<Cell> constructCells(vector<Wire> &Uwires, vector<Wire> &Vwires, vector<Wire> &Ywires)
{
  vector<Cell> cells;

  const Double_t maxHeight = heightToWidthRatio*wirePitchY*numYwires; 
  const Double_t UdeltaY = wirePitchY/tan((PI/180.0)*angleU);
  const Double_t VdeltaY = -1.0*wirePitchY/tan((PI/180.0)*angleV);
  const Double_t UspacingOnWire = fabs(wirePitchU/sin((PI/180.0)*angleU));
  const Double_t VspacingOnWire = fabs(wirePitchV/sin((PI/180.0)*angleV));

  Double_t Zval = firstYwireZval;
  Double_t Uoffset = maxHeight-firstYwireUoffsetYval;
  Double_t Voffset = firstYwireVoffsetYval;

  while(Uoffset < maxHeight-((UspacingOnWire-UdeltaY)/2.0)-epsilon)
    Uoffset += UspacingOnWire;
  while(Voffset > ((VspacingOnWire+VdeltaY)/2.0)+epsilon)
    Voffset -= VspacingOnWire;

  vector<Cell> cellChain;
  Int_t firstCellID = 0;
  for(Int_t i = 0; i < numYwires; i++)
  {
    cellChain = constructCellChain(firstCellID,Zval,Uoffset,Voffset); // TEMPORARY (using existing code, eventually replace with constructCell that takes three wires as input and no secondary loop)
    firstCellID += cellChain.size();

    for(Int_t j = 0; j < cellChain.size(); j++)
    {
      // NOTE:  In principal some wires that don't exist could be associated with a cell.  In the future we should merge such a cell with its adjacent cell that is formed from three wires that actually exist in the TPC.  This happens almost exclusively at the corners.  Also, the cells near the edges should change in shape due to different "closest wires"
      if(cellChain.at(j).UwireID < Uwires.size())
        Uwires.at(cellChain.at(j).UwireID).cellIDs.push_back(cellChain.at(j).ID);
      if(cellChain.at(j).VwireID < Vwires.size())
        Vwires.at(cellChain.at(j).VwireID).cellIDs.push_back(cellChain.at(j).ID);
      if(cellChain.at(j).YwireID < Ywires.size())
        Ywires.at(cellChain.at(j).YwireID).cellIDs.push_back(cellChain.at(j).ID);
  
      cells.push_back(cellChain.at(j));
    }

    Zval += wirePitchY;
    Uoffset += UdeltaY;
    Voffset += VdeltaY;

    while(Uoffset < maxHeight-((UspacingOnWire-UdeltaY)/2.0)-epsilon)
      Uoffset += UspacingOnWire;
    while(Voffset < ((VdeltaY-VspacingOnWire)/2.0)-epsilon)
      Voffset += VspacingOnWire;
  }

  return cells;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
// constructCellChain - Construct Cell chains using pairs of crossings from a U wire and a V wire
/////////////////////////////////////////////////////////////////////////////////////////////////////
vector<Cell> constructCellChain(Int_t firstCellID, Double_t wireZval, Double_t YvalOffsetU, Double_t YvalOffsetV)
{ 
  vector<Cell> cellChain;
 
  const Double_t maxHeight = heightToWidthRatio*wirePitchY*numYwires; 
  const Double_t UdeltaY = wirePitchY/tan((PI/180.0)*angleU);
  const Double_t VdeltaY = -1.0*wirePitchY/tan((PI/180.0)*angleV);
  const Double_t UspacingOnWire = fabs(wirePitchU/sin((PI/180.0)*angleU));
  const Double_t VspacingOnWire = fabs(wirePitchV/sin((PI/180.0)*angleV));

  Int_t numUcrosses = TMath::Ceil(((UdeltaY-UspacingOnWire)/2.0+YvalOffsetU)/UspacingOnWire)+1;
  Int_t numVcrosses = TMath::Ceil((maxHeight-(VdeltaY+VspacingOnWire)/2.0-YvalOffsetV)/VspacingOnWire)+1;

  Int_t cellID = firstCellID;
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
      if(formsCell(YvalOffsetU-i*UspacingOnWire,YvalOffsetV+j*VspacingOnWire) == true)
      {
        flag1 = true;
        cell = constructCell(cellID,wireZval,YvalOffsetU-i*UspacingOnWire,YvalOffsetV+j*VspacingOnWire);

        if(cell.vertices.size() > 2)
	{
          cellChain.push_back(cell);
          cellID++;
	}
      }
      else if(flag1 == true)
        flag2 = true;

      j++;
    }
  }
 
  return cellChain;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
// constructCell - Construct one Cell, including vertices, area, and center point
/////////////////////////////////////////////////////////////////////////////////////////////////////
Cell constructCell(Int_t cellID, Double_t YwireZval, Double_t UwireYval, Double_t VwireYval)
{
  Cell cell;
  cell.ID = cellID;

  cell.vertices = getCellVertices(YwireZval,UwireYval,VwireYval);
  cell.center = calcCellCenter(cell.vertices);
  cell.area = calcCellArea(cell.vertices);

  cell.trueCharge = 0.0;
  cell.recoCharge = 0.0;

  cell.UwireID = getUwireID(UwireYval,YwireZval);
  cell.VwireID = getVwireID(VwireYval,YwireZval);
  cell.YwireID = getYwireID(YwireZval);

  cell.hitType = kNoHit;

  return cell;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
// formsCell - Check whether or not a U/V crossing pair on a particular Y wire forms a Cell
/////////////////////////////////////////////////////////////////////////////////////////////////////
Bool_t formsCell(Double_t UwireYval, Double_t VwireYval)
{
  Bool_t isCell;

  const Double_t UspacingOnWire = fabs(wirePitchU/sin((PI/180.0)*angleU));
  const Double_t VspacingOnWire = fabs(wirePitchV/sin((PI/180.0)*angleV));
  const Double_t tanU = tan((PI/180.0)*angleU);
  const Double_t tanV = tan((PI/180.0)*angleV);
  const Double_t maxHeight = heightToWidthRatio*wirePitchY*numYwires; 

  Double_t deltaY;
  if(UwireYval > VwireYval)
    deltaY = (UwireYval-UspacingOnWire/2.0)-(VwireYval+VspacingOnWire/2.0);
  else
    deltaY = (VwireYval-VspacingOnWire/2.0)-(UwireYval+UspacingOnWire/2.0);

  if((UwireYval+UspacingOnWire/2.0 < -1.0*epsilon) && (VwireYval+VspacingOnWire/2.0 < -1.0*epsilon))
    isCell = false;
  else if((UwireYval-UspacingOnWire/2.0 > maxHeight) && (VwireYval-VspacingOnWire/2.0 > maxHeight+epsilon))
    isCell = false;
  else if(deltaY < epsilon)
    isCell = true;
  else if(UwireYval == VwireYval)
    isCell = true;
  else if(fabs((tanU*tanV*deltaY)/(tanU+tanV))-wirePitchY/2.0 < epsilon)
    isCell = true;
  else
    isCell = false;

  return isCell;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
// getCellVertices - Find vertices that define the boundaries of a Cell
/////////////////////////////////////////////////////////////////////////////////////////////////////
vector<pair<Double_t,Double_t> > getCellVertices(Double_t YwireZval, Double_t UwireYval, Double_t VwireYval)
{
  vector<pair<Double_t,Double_t> > cellVertices;

  const Double_t maxHeight = heightToWidthRatio*wirePitchY*numYwires;
  const Double_t UspacingOnWire = fabs(wirePitchU/sin((PI/180.0)*angleU));
  const Double_t VspacingOnWire = fabs(wirePitchV/sin((PI/180.0)*angleV));

  const Double_t U1slope = 1.0/tan((PI/180.0)*angleU);
  const Double_t U2slope = U1slope;
  const Double_t U1intercept = UwireYval - UspacingOnWire/2.0;
  const Double_t U2intercept = UwireYval + UspacingOnWire/2.0;

  const Double_t V1slope = -1.0/tan((PI/180.0)*angleV);
  const Double_t V2slope = V1slope;
  const Double_t V1intercept = VwireYval - VspacingOnWire/2.0;
  const Double_t V2intercept = VwireYval + VspacingOnWire/2.0;

  const Double_t Y1Zval = YwireZval - wirePitchY/2.0;
  const Double_t Y2Zval = YwireZval + wirePitchY/2.0;

  pair<Double_t,Double_t> U1V1intersectionPoint;
  pair<Double_t,Double_t> U1V2intersectionPoint;
  pair<Double_t,Double_t> U2V1intersectionPoint;
  pair<Double_t,Double_t> U2V2intersectionPoint;

  U1V1intersectionPoint = getIntersectionUV(YwireZval,U1slope,U1intercept,V1slope,V1intercept);
  U1V2intersectionPoint = getIntersectionUV(YwireZval,U1slope,U1intercept,V2slope,V2intercept);
  U2V1intersectionPoint = getIntersectionUV(YwireZval,U2slope,U2intercept,V1slope,V1intercept);
  U2V2intersectionPoint = getIntersectionUV(YwireZval,U2slope,U2intercept,V2slope,V2intercept);

  Double_t UVminZval = U1V1intersectionPoint.first;
  Double_t UVmaxZval = U1V1intersectionPoint.first;
  Double_t UVminYval = U1V1intersectionPoint.second;
  Double_t UVmaxYval = U1V1intersectionPoint.second;

  if(U1V2intersectionPoint.first < UVminZval+epsilon)
    UVminZval = U1V2intersectionPoint.first;
  if(U1V2intersectionPoint.first > UVmaxZval-epsilon)
    UVmaxZval = U1V2intersectionPoint.first;
  if(U1V2intersectionPoint.second < UVminYval+epsilon)
    UVminYval = U1V2intersectionPoint.second;
  if(U1V2intersectionPoint.second > UVmaxYval-epsilon)
    UVmaxYval = U1V2intersectionPoint.second;

  if(U2V1intersectionPoint.first < UVminZval+epsilon)
    UVminZval = U2V1intersectionPoint.first;
  if(U2V1intersectionPoint.first > UVmaxZval-epsilon)
    UVmaxZval = U2V1intersectionPoint.first;
  if(U2V1intersectionPoint.second < UVminYval+epsilon)
    UVminYval = U2V1intersectionPoint.second;
  if(U2V1intersectionPoint.second > UVmaxYval-epsilon)
    UVmaxYval = U2V1intersectionPoint.second;

  if(U2V2intersectionPoint.first < UVminZval+epsilon)
    UVminZval = U2V2intersectionPoint.first;
  if(U2V2intersectionPoint.first > UVmaxZval-epsilon)
    UVmaxZval = U2V2intersectionPoint.first;
  if(U2V2intersectionPoint.second < UVminYval+epsilon)
    UVminYval = U2V2intersectionPoint.second;
  if(U2V2intersectionPoint.second > UVmaxYval-epsilon)
    UVmaxYval = U2V2intersectionPoint.second;

  pair<Double_t,Double_t> Y1U1intersectionPoint;
  pair<Double_t,Double_t> Y1U2intersectionPoint;
  pair<Double_t,Double_t> Y1V1intersectionPoint;
  pair<Double_t,Double_t> Y1V2intersectionPoint;
  pair<Double_t,Double_t> Y2U1intersectionPoint;
  pair<Double_t,Double_t> Y2U2intersectionPoint;
  pair<Double_t,Double_t> Y2V1intersectionPoint;
  pair<Double_t,Double_t> Y2V2intersectionPoint;

  Y1U1intersectionPoint = getIntersectionY(YwireZval,Y1Zval,U1slope,U1intercept);
  Y1U2intersectionPoint = getIntersectionY(YwireZval,Y1Zval,U2slope,U2intercept);
  Y1V1intersectionPoint = getIntersectionY(YwireZval,Y1Zval,V1slope,V1intercept);
  Y1V2intersectionPoint = getIntersectionY(YwireZval,Y1Zval,V2slope,V2intercept);
  Y2U1intersectionPoint = getIntersectionY(YwireZval,Y2Zval,U1slope,U1intercept);
  Y2U2intersectionPoint = getIntersectionY(YwireZval,Y2Zval,U2slope,U2intercept);
  Y2V1intersectionPoint = getIntersectionY(YwireZval,Y2Zval,V1slope,V1intercept);
  Y2V2intersectionPoint = getIntersectionY(YwireZval,Y2Zval,V2slope,V2intercept);

  if((Y1U1intersectionPoint.second < UVmaxYval+epsilon) && (Y1U1intersectionPoint.second > UVminYval+epsilon) && (Y1U1intersectionPoint.first < UVmaxZval-epsilon) && (Y1U1intersectionPoint.first > UVminZval+epsilon))
    cellVertices.push_back(Y1U1intersectionPoint);
  if((Y1U2intersectionPoint.second < UVmaxYval-epsilon) && (Y1U2intersectionPoint.second > UVminYval+epsilon) && (Y1U2intersectionPoint.first < UVmaxZval-epsilon) && (Y1U2intersectionPoint.first > UVminZval+epsilon))
    cellVertices.push_back(Y1U2intersectionPoint);
  if((Y1V1intersectionPoint.second < UVmaxYval-epsilon) && (Y1V1intersectionPoint.second > UVminYval+epsilon) && (Y1V1intersectionPoint.first < UVmaxZval-epsilon) && (Y1V1intersectionPoint.first > UVminZval+epsilon))
    cellVertices.push_back(Y1V1intersectionPoint);
  if((Y1V2intersectionPoint.second < UVmaxYval-epsilon) && (Y1V2intersectionPoint.second > UVminYval+epsilon) && (Y1V2intersectionPoint.first < UVmaxZval-epsilon) && (Y1V2intersectionPoint.first > UVminZval+epsilon))
    cellVertices.push_back(Y1V2intersectionPoint);

  if((Y2U1intersectionPoint.second < UVmaxYval-epsilon) && (Y2U1intersectionPoint.second > UVminYval+epsilon) && (Y2U1intersectionPoint.first < UVmaxZval-epsilon) && (Y2U1intersectionPoint.first > UVminZval+epsilon))
    cellVertices.push_back(Y2U1intersectionPoint);
  if((Y2U2intersectionPoint.second < UVmaxYval-epsilon) && (Y2U2intersectionPoint.second > UVminYval+epsilon) && (Y2U2intersectionPoint.first < UVmaxZval-epsilon) && (Y2U2intersectionPoint.first > UVminZval+epsilon))
    cellVertices.push_back(Y2U2intersectionPoint);
  if((Y2V1intersectionPoint.second < UVmaxYval-epsilon) && (Y2V1intersectionPoint.second > UVminYval+epsilon) && (Y2V1intersectionPoint.first < UVmaxZval-epsilon) && (Y2V1intersectionPoint.first > UVminZval+epsilon))
    cellVertices.push_back(Y2V1intersectionPoint);
  if((Y2V2intersectionPoint.second < UVmaxYval-epsilon) && (Y2V2intersectionPoint.second > UVminYval+epsilon) && (Y2V2intersectionPoint.first < UVmaxZval-epsilon) && (Y2V2intersectionPoint.first > UVminZval+epsilon))
    cellVertices.push_back(Y2V2intersectionPoint);

  if((U1V1intersectionPoint.first >= Y1Zval-epsilon) && (U1V1intersectionPoint.first <= Y2Zval+epsilon))
    cellVertices.push_back(U1V1intersectionPoint);
  if((U1V2intersectionPoint.first >= Y1Zval-epsilon) && (U1V2intersectionPoint.first <= Y2Zval+epsilon))
    cellVertices.push_back(U1V2intersectionPoint);
  if((U2V1intersectionPoint.first >= Y1Zval-epsilon) && (U2V1intersectionPoint.first <= Y2Zval+epsilon))
    cellVertices.push_back(U2V1intersectionPoint);
  if((U2V2intersectionPoint.first >= Y1Zval-epsilon) && (U2V2intersectionPoint.first <= Y2Zval+epsilon))
    cellVertices.push_back(U2V2intersectionPoint);

  vector<pair<Double_t,Double_t> > cellVerticesSorted;
  cellVerticesSorted = sortVertices(cellVertices);

  if(UVminZval < (firstYwireZval-0.5*wirePitchY+leftEdgeOffsetZval)+epsilon)
    cellVerticesSorted = trimEdgeCellVertices(cellVerticesSorted,1,firstYwireZval-0.5*wirePitchY+leftEdgeOffsetZval);
  else if(UVmaxZval > (firstYwireZval+(numYwires-0.5)*wirePitchY-rightEdgeOffsetZval)-epsilon)
    cellVerticesSorted = trimEdgeCellVertices(cellVerticesSorted,2,(firstYwireZval+(numYwires-0.5)*wirePitchY-rightEdgeOffsetZval));

  if(UVminYval < epsilon)
    cellVerticesSorted = trimEdgeCellVertices(cellVerticesSorted,3,0.0);
  else if(UVmaxYval > maxHeight-epsilon)
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
// getIntersectionY - Get intersection point of a Y wire and either a U or a V wire
/////////////////////////////////////////////////////////////////////////////////////////////////////
pair<Double_t, Double_t> getIntersectionY(Double_t ZvalOffset, Double_t YwireZval, Double_t slope, Double_t intercept)
{
  Double_t Yval = slope*(YwireZval-ZvalOffset) + intercept;
  Double_t Zval = YwireZval;

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
        if(((edgeVal > Zval1+epsilon) && (edgeVal < Zval2-epsilon)) || ((edgeVal > Zval2+epsilon) && (edgeVal < Zval1-epsilon)))
	{
          vertex.first = edgeVal;
          vertex.second = slope*edgeVal+intercept;
          edgeCellVertices.push_back(vertex);
	}
      }
      else if((edgeType == 3) || (edgeType == 4))
      {
        if(((edgeVal > Yval1+epsilon) && (edgeVal < Yval2-epsilon)) || ((edgeVal > Yval2+epsilon) && (edgeVal < Yval1-epsilon)))
	{
          vertex.first = (edgeVal-intercept)/slope;
          vertex.second = edgeVal;
          edgeCellVertices.push_back(vertex);
	}
      }
    }
    else if((edgeType == 3) || (edgeType == 4))
    {
      if(((edgeVal > Yval1+epsilon) && (edgeVal < Yval2-epsilon)) || ((edgeVal > Yval2+epsilon) && (edgeVal < Yval1-epsilon)))
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
      return (vertex.first < edgeVal-epsilon);
    case 2:
      return (vertex.first > edgeVal+epsilon);
    case 3:
      return (vertex.second < edgeVal-epsilon);
    case 4:
      return (vertex.second > edgeVal+epsilon);
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
// getUwireYval - Return Y value of particular U wire at a given Z value
/////////////////////////////////////////////////////////////////////////////////////////////////////
Double_t getUwireYval(Int_t IDnum, Double_t Zval)
{
  const Double_t maxHeight = heightToWidthRatio*wirePitchY*numYwires;
  const Double_t UspacingOnWire = fabs(wirePitchU/sin((PI/180.0)*angleU));

  return Zval/tan((PI/180.0)*angleU) + maxHeight - (firstYwireZval/tan((PI/180.0)*angleU)) - firstYwireUoffsetYval - UspacingOnWire*IDnum;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
// getUwireZval - Return Z value of particular U wire at a given Y value
/////////////////////////////////////////////////////////////////////////////////////////////////////
Double_t getUwireZval(Int_t IDnum, Double_t Yval)
{
  const Double_t maxHeight = heightToWidthRatio*wirePitchY*numYwires;
  const Double_t UspacingOnWire = fabs(wirePitchU/sin((PI/180.0)*angleU));

  return tan((PI/180.0)*angleU)*(Yval - maxHeight + (firstYwireZval/tan((PI/180.0)*angleU)) + firstYwireUoffsetYval + UspacingOnWire*IDnum);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
// getUwireID - Return ID of U wire nearest to given {Y,Z} point
/////////////////////////////////////////////////////////////////////////////////////////////////////
Int_t getUwireID(Double_t Yval, Double_t Zval)
{
  const Double_t maxHeight = heightToWidthRatio*wirePitchY*numYwires;
  const Double_t UspacingOnWire = fabs(wirePitchU/sin((PI/180.0)*angleU));

  return round((Zval/tan((PI/180.0)*angleU) + maxHeight - (firstYwireZval/tan((PI/180.0)*angleU)) - firstYwireUoffsetYval - Yval)/UspacingOnWire);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
// getVwireYval - Return Y value of particular V wire at a given Z value
/////////////////////////////////////////////////////////////////////////////////////////////////////
Double_t getVwireYval(Int_t IDnum, Double_t Zval)
{
  const Double_t VspacingOnWire = fabs(wirePitchV/sin((PI/180.0)*angleV));

  return -1.0*Zval/tan((PI/180.0)*angleV) + firstYwireZval/tan((PI/180.0)*angleV) + firstYwireVoffsetYval + VspacingOnWire*IDnum;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
// getVwireZval - Return Z value of particular V wire at a given Y value
/////////////////////////////////////////////////////////////////////////////////////////////////////
Double_t getVwireZval(Int_t IDnum, Double_t Yval)
{
  const Double_t VspacingOnWire = fabs(wirePitchV/sin((PI/180.0)*angleV));

  return tan((PI/180.0)*angleV)*(firstYwireZval/tan((PI/180.0)*angleV) + firstYwireVoffsetYval + VspacingOnWire*IDnum - Yval);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
// getVwireID - Return ID of V wire nearest to given {Y,Z} point
/////////////////////////////////////////////////////////////////////////////////////////////////////
Int_t getVwireID(Double_t Yval, Double_t Zval)
{
  const Double_t VspacingOnWire = fabs(wirePitchV/sin((PI/180.0)*angleV));

  return round((Zval/tan((PI/180.0)*angleV) - firstYwireZval/tan((PI/180.0)*angleV) - firstYwireVoffsetYval + Yval)/VspacingOnWire);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
// getYwireZval - Return Z value of particular Y wire for all Y values
/////////////////////////////////////////////////////////////////////////////////////////////////////
Double_t getYwireZval(Int_t IDnum)
{
  return firstYwireZval + wirePitchY*IDnum;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
// getYwireID - Return ID of Y wire nearest to given {Y,Z} point
/////////////////////////////////////////////////////////////////////////////////////////////////////
Int_t getYwireID(Double_t Zval)
{
  return round((Zval-firstYwireZval)/wirePitchY);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
// getWireEndpoints - Return endpoints of a particular Wire on the wire plane edges
/////////////////////////////////////////////////////////////////////////////////////////////////////
pair<pair<Double_t,Double_t>,pair<Double_t,Double_t> > getWireEndpoints(Int_t wireID, Plane_t wirePlane)
{
  pair<pair<Double_t,Double_t>,pair<Double_t,Double_t> > endpoints;

  const Double_t maxHeight = heightToWidthRatio*wirePitchY*numYwires;

  if(wirePlane == kUPlane)
  {
    endpoints.first.first = max(firstYwireZval-0.5*wirePitchY+leftEdgeOffsetZval,getUwireZval(wireID,0.0));
    endpoints.first.second = max(0.0,getUwireYval(wireID,firstYwireZval-0.5*wirePitchY+leftEdgeOffsetZval));
    endpoints.second.first = min(firstYwireZval+(numYwires-0.5)*wirePitchY-rightEdgeOffsetZval,getUwireZval(wireID,maxHeight));
    endpoints.second.second = min(maxHeight,getUwireYval(wireID,firstYwireZval+(numYwires-0.5)*wirePitchY-rightEdgeOffsetZval));
  }
  else if(wirePlane == kVPlane)
  {
    endpoints.first.first = max(firstYwireZval-0.5*wirePitchY+leftEdgeOffsetZval,getVwireZval(wireID,maxHeight));
    endpoints.first.second = min(maxHeight,getVwireYval(wireID,firstYwireZval-0.5*wirePitchY+leftEdgeOffsetZval));
    endpoints.second.first = min(firstYwireZval+(numYwires-0.5)*wirePitchY-rightEdgeOffsetZval,getVwireZval(wireID,0.0));
    endpoints.second.second = max(0.0,getVwireYval(wireID,firstYwireZval+(numYwires-0.5)*wirePitchY-rightEdgeOffsetZval));
  }
  else if(wirePlane == kYPlane)
  {
    endpoints.first.first = getYwireZval(wireID);
    endpoints.first.second = 0.0;
    endpoints.second.first = endpoints.first.first;
    endpoints.second.second = maxHeight;
  }

  return endpoints;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
// addCharges - Add true charge to Cells within CellMap, consistently adding charge to relevant Wires
/////////////////////////////////////////////////////////////////////////////////////////////////////
void addCharges(CellMap &cellMap)
{
  TRandom3 *rand;
  Double_t randNum;
  Double_t tempCharge;
  for(Int_t i = 0; i < cellMap.cells.size(); i++)
  {
    // Temporarily assign random charge for visualization purposes
    rand = new TRandom3(0);
    randNum = rand->Uniform(1.0);
    if(randNum < 0.01*cellMap.cells.at(i).area/(pow(pow(wirePitchY*wirePitchU*wirePitchV,1.0/3.0),2)))
    {
      tempCharge = 100.0*randNum;
      cellMap.cells.at(i).trueCharge = tempCharge;
      // NOTE:  (SECOND COPY) In principal some wires that don't exist could be associated with a cell.  In the future we should merge such a cell with its adjacent cell that is formed from three wires that actually exist in the TPC.  This happens almost exclusively at the corners.  Also, the cells near the edges should change in shape due to different "closest wires"
      if(cellMap.cells.at(i).UwireID < cellMap.Uwires.size())
        cellMap.Uwires.at(cellMap.cells.at(i).UwireID).charge += tempCharge;
      if(cellMap.cells.at(i).VwireID < cellMap.Vwires.size())
        cellMap.Vwires.at(cellMap.cells.at(i).VwireID).charge += tempCharge;
      if(cellMap.cells.at(i).YwireID < cellMap.Ywires.size())
        cellMap.Ywires.at(cellMap.cells.at(i).YwireID).charge += tempCharge;
    }
    delete rand;
  }

  assignHitTypes(cellMap);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
// assignHitTypes - Determine whether cell contains a hit or not, and whether it is real or fake
/////////////////////////////////////////////////////////////////////////////////////////////////////
void assignHitTypes(CellMap &cellMap)
{
  Double_t cellCharge;
  Double_t UwireCharge;
  Double_t VwireCharge;
  Double_t YwireCharge;
  for(Int_t i = 0; i < cellMap.cells.size(); i++)
  {
    cellCharge = cellMap.cells.at(i).trueCharge;

    // NOTE:  (THIRD COPY) In principal some wires that don't exist could be associated with a cell.  In the future we should merge such a cell with its adjacent cell that is formed from three wires that actually exist in the TPC.  This happens almost exclusively at the corners.  Also, the cells near the edges should change in shape due to different "closest wires"
    if(cellMap.cells.at(i).UwireID < cellMap.Uwires.size())
      UwireCharge = cellMap.Uwires.at(cellMap.cells.at(i).UwireID).charge;
    else
      UwireCharge = 0.0;

    if(cellMap.cells.at(i).VwireID < cellMap.Vwires.size())
      VwireCharge = cellMap.Vwires.at(cellMap.cells.at(i).VwireID).charge;
    else
      VwireCharge = 0.0;

    if(cellMap.cells.at(i).YwireID < cellMap.Ywires.size())
      YwireCharge = cellMap.Ywires.at(cellMap.cells.at(i).YwireID).charge;
    else
      YwireCharge = 0.0;
    
    if((UwireCharge > 0.0) && (VwireCharge > 0.0) && (YwireCharge > 0.0))
    {
      if(cellCharge > 0.0)
        cellMap.cells.at(i).hitType = kRealHit;
      else
        cellMap.cells.at(i).hitType = kFakeHit;
    }
    else
    {
      cellMap.cells.at(i).hitType = kNoHit;
    }
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
// drawCellMap - Draw complete Cell map, or portion thereof
/////////////////////////////////////////////////////////////////////////////////////////////////////
void drawCellMap(CellMap const& cellMap, Int_t numWires, Int_t numCells)
{
  const Double_t shadeMinCellTrue = 0.5;
  const Double_t shadeMinCellFake = 0.3;
  const Double_t shadeMinWire = 0.1;

  Int_t maxWiresY = min(numWires,(Int_t)cellMap.Ywires.size());
  if(maxWiresY < 0)
    maxWiresY = cellMap.Ywires.size();

  const Double_t maxHeight = heightToWidthRatio*wirePitchY*numYwires;
  const Int_t colorFactor = ((sizeof(colorVec)/sizeof(*colorVec))/2);
  const Int_t colorFactor2 = ((sizeof(colorVec2)/sizeof(*colorVec2))/2);
  const Int_t colorFactor3 = ((sizeof(colorVec3)/sizeof(*colorVec3))/2);

  const Double_t maxZ = (maxWiresY-1)*wirePitchY;
  const Double_t maxY = (maxWiresY-1)*wirePitchY*heightToWidthRatio;
  const Double_t diagLength = sqrt(pow(maxZ,2) + pow(maxY,2));
  const Double_t diagAngle = (180.0/PI)*atan(1.0/heightToWidthRatio);

  const Double_t Uoffset = firstYwireUoffsetYval*(sin((PI/180.0)*angleU)/sin((PI/180.0)*(180.0-diagAngle-angleU)));
  const Double_t Voffset = firstYwireVoffsetYval*(sin((PI/180.0)*angleV)/sin((PI/180.0)*(180.0-diagAngle-angleV)));
  const Double_t maxWiresU = TMath::Floor(((diagLength-Uoffset)*sin((PI/180.0)*(diagAngle+angleU)))/wirePitchU);
  const Double_t maxWiresV = TMath::Floor(((diagLength-Voffset)*sin((PI/180.0)*(diagAngle+angleV)))/wirePitchV);

  TCanvas *c1 = new TCanvas("c1","c1",1600.0,800.0*(heightToWidthRatio/0.5));

  TGraph *cellGraph;
  TMultiGraph *cellMapGraph = new TMultiGraph();

  Int_t colorNum;
  Int_t colorNum2;
  Int_t colorNum3;
  Int_t plotNcell;
  Int_t plotNvert;
  Double_t plotX[7];
  Double_t plotY[7];

  Double_t tempCharge;
  Hit_t cellHitType;
  Double_t maxCellCharge = -999999999.0;
  Double_t minCellCharge = 999999999.0;
  Double_t maxWireCharge = -999999999.0;
  Double_t minWireCharge = 999999999.0;
  Double_t UwireCharge;
  Double_t VwireCharge;
  Double_t YwireCharge;
  Double_t maxFakeCharge = -999999999.0;
  Double_t minFakeCharge = 999999999.0;
  if(plotMode == 3)
  {
    for(Int_t i = 0; i < maxWiresY; i++)
    {
      plotNcell = min(numCells,(Int_t)cellMap.Ywires.at(i).cellIDs.size());
      if(plotNcell < 0)
        plotNcell = cellMap.Ywires.at(i).cellIDs.size();
    
      for(Int_t j = 0; j < plotNcell; j++)
      {
        tempCharge = cellMap.cells.at(cellMap.Ywires.at(i).cellIDs.at(j)).trueCharge;
	if(tempCharge > maxCellCharge)
          maxCellCharge = tempCharge;
	if(tempCharge < minCellCharge)
          minCellCharge = tempCharge;

        // NOTE:  (FOURTH COPY) In principal some wires that don't exist could be associated with a cell.  In the future we should merge such a cell with its adjacent cell that is formed from three wires that actually exist in the TPC.  This happens almost exclusively at the corners.  Also, the cells near the edges should change in shape due to different "closest wires"

        if(cellMap.cells.at(cellMap.Ywires.at(i).cellIDs.at(j)).UwireID < cellMap.Uwires.size())
          UwireCharge = cellMap.Uwires.at(cellMap.cells.at(cellMap.Ywires.at(i).cellIDs.at(j)).UwireID).charge;
        else
          UwireCharge = 0.0;
        
        if(cellMap.cells.at(cellMap.Ywires.at(i).cellIDs.at(j)).VwireID < cellMap.Vwires.size())
          VwireCharge = cellMap.Vwires.at(cellMap.cells.at(cellMap.Ywires.at(i).cellIDs.at(j)).VwireID).charge;
        else
          VwireCharge = 0.0;
        
        YwireCharge = cellMap.Ywires.at(i).charge;
        
        if(cellMap.cells.at(cellMap.Ywires.at(i).cellIDs.at(j)).hitType == kFakeHit)
        {
          if(UwireCharge+VwireCharge+YwireCharge > maxFakeCharge)
            maxFakeCharge = UwireCharge+VwireCharge+YwireCharge;
          if(UwireCharge+VwireCharge+YwireCharge < minFakeCharge)
            minFakeCharge = UwireCharge+VwireCharge+YwireCharge;
        }
      }
    }

    for(Int_t i = 0; i < maxWiresU; i++)
    {
      tempCharge = cellMap.Uwires.at(i).charge;
      if(tempCharge > maxWireCharge)
        maxWireCharge = tempCharge;
      if(tempCharge < minWireCharge)
        minWireCharge = tempCharge;
    }
    for(Int_t i = 0; i < maxWiresV; i++)
    {
      tempCharge = cellMap.Vwires.at(i).charge;
      if(tempCharge > maxWireCharge)
        maxWireCharge = tempCharge;
      if(tempCharge < minWireCharge)
        minWireCharge = tempCharge;
    }
    for(Int_t i = 0; i < maxWiresY; i++)
    {
      tempCharge = cellMap.Ywires.at(i).charge;
      if(tempCharge > maxWireCharge)
        maxWireCharge = tempCharge;
      if(tempCharge < minWireCharge)
        minWireCharge = tempCharge;
    }      
  }
  for(Int_t i = 0; i < maxWiresY; i++)
  {
    plotNcell = min(numCells,(Int_t)cellMap.Ywires.at(i).cellIDs.size());
    if(plotNcell < 0)
      plotNcell = cellMap.Ywires.at(i).cellIDs.size();

    for(Int_t j = 0; j < plotNcell; j++)
    {
      plotNvert = cellMap.cells.at(cellMap.Ywires.at(i).cellIDs.at(j)).vertices.size();

      for(Int_t k = 0; k < plotNvert; k++)
      {
        plotX[k] = cellMap.cells.at(cellMap.Ywires.at(i).cellIDs.at(j)).vertices.at(k).first;
        plotY[k] = cellMap.cells.at(cellMap.Ywires.at(i).cellIDs.at(j)).vertices.at(k).second;
      }
      colorNum = colorVec[colorFactor*(cellMap.Ywires.at(i).ID % 2)+(cellMap.cells.at(cellMap.Ywires.at(i).cellIDs.at(j)).ID % colorFactor)];
      colorNum2 = colorVec2[colorFactor2*(cellMap.Ywires.at(i).ID % 2)+(cellMap.cells.at(cellMap.Ywires.at(i).cellIDs.at(j)).ID % colorFactor2)];
      colorNum3 = colorVec3[colorFactor3*(cellMap.Ywires.at(i).ID % 2)+(cellMap.cells.at(cellMap.Ywires.at(i).cellIDs.at(j)).ID % colorFactor3)];

      cellGraph = new TGraph(plotNvert,plotX,plotY);
      cellGraph->SetTitle("");

      tempCharge = cellMap.cells.at(cellMap.Ywires.at(i).cellIDs.at(j)).trueCharge;
      cellHitType = cellMap.cells.at(cellMap.Ywires.at(i).cellIDs.at(j)).hitType;

      // NOTE:  (FIFTH COPY) In principal some wires that don't exist could be associated with a cell.  In the future we should merge such a cell with its adjacent cell that is formed from three wires that actually exist in the TPC.  This happens almost exclusively at the corners.  Also, the cells near the edges should change in shape due to different "closest wires"
      if(cellMap.cells.at(cellMap.Ywires.at(i).cellIDs.at(j)).UwireID < cellMap.Uwires.size())
        UwireCharge = cellMap.Uwires.at(cellMap.cells.at(cellMap.Ywires.at(i).cellIDs.at(j)).UwireID).charge;
      else
        UwireCharge = 0.0;
    
      if(cellMap.cells.at(cellMap.Ywires.at(i).cellIDs.at(j)).VwireID < cellMap.Vwires.size())
        VwireCharge = cellMap.Vwires.at(cellMap.cells.at(cellMap.Ywires.at(i).cellIDs.at(j)).VwireID).charge;
      else
        VwireCharge = 0.0;
    
      YwireCharge = cellMap.Ywires.at(i).charge;

      if(plotMode == 2)
      {
        if(cellHitType == kRealHit)
	{
          cellGraph->SetLineColor(colorNum2);
          cellGraph->SetFillColorAlpha(colorNum2,1.0);
	}
	else
	{
          cellGraph->SetLineColor(colorNum3);
          cellGraph->SetFillColorAlpha(colorNum3,0.2);
	}
      }
      else if(plotMode == 3)
      {
        if(cellHitType == kRealHit)
	{
          cellGraph->SetLineColor(kGreen+1);
          cellGraph->SetFillColorAlpha(kGreen+1,shadeMinCellTrue+((tempCharge-minCellCharge)/(maxCellCharge-minCellCharge))*(1.0-shadeMinCellTrue));
	}
	else if(cellHitType == kFakeHit)
	{
          cellGraph->SetLineColor(kRed-4);
          cellGraph->SetFillColorAlpha(kRed-4,shadeMinCellFake+((UwireCharge+VwireCharge+YwireCharge-minFakeCharge)/(maxFakeCharge-minFakeCharge))*(1.0-shadeMinCellFake));
	}
	else
	{
          cellGraph->SetLineColor(colorNum3);
          cellGraph->SetFillColorAlpha(colorNum3,0.2);
	}
      }
      else
      {
        cellGraph->SetLineColor(colorNum);
        cellGraph->SetFillColor(colorNum);
      }
      cellMapGraph->Add(cellGraph);
    }
  }
  cellMapGraph->Draw("AFL");

  // NOTE:  right now plotting doesn't fully support input parameters limiting number of wires/cells
  cellMapGraph->SetTitle(Form("#theta_{U} = %.1f#circ,_{} #theta_{V} = %.1f#circ,_{} N_{wires} = %d",angleU,angleV,numYwires));
  cellMapGraph->GetXaxis()->SetTitle("Z [cm]");
  cellMapGraph->GetXaxis()->SetTitleOffset(1.0);
  cellMapGraph->GetXaxis()->SetTitleSize(0.04);
  cellMapGraph->GetYaxis()->SetTitle("Y [cm]");
  cellMapGraph->GetYaxis()->SetTitleOffset(0.8*(heightToWidthRatio/0.5));
  cellMapGraph->GetYaxis()->SetTitleSize(0.04);
  cellMapGraph->GetXaxis()->SetLimits(firstYwireZval-0.5*wirePitchY+leftEdgeOffsetZval,firstYwireZval+(numYwires-0.5)*wirePitchY-rightEdgeOffsetZval);
  cellMapGraph->GetHistogram()->SetMinimum(0.0);
  cellMapGraph->GetHistogram()->SetMaximum(maxHeight);

  gPad->Update();
  gPad->RedrawAxis();

  if(plotMode == 3)
  {
    TLine wireLine;
    wireLine.SetLineWidth(2.0);

    pair<pair<Double_t,Double_t>,pair<Double_t,Double_t> > lineEnds;

    for(Int_t i = 0; i < maxWiresU; i++)
    {
      tempCharge = cellMap.Uwires.at(i).charge;
      if(tempCharge > 0.0)
      {
        lineEnds = getWireEndpoints(i,kUPlane);
        lineEnds.first.first /= numYwires*wirePitchY-leftEdgeOffsetZval-rightEdgeOffsetZval;
        lineEnds.first.second /= maxHeight;
        lineEnds.second.first /= numYwires*wirePitchY-leftEdgeOffsetZval-rightEdgeOffsetZval;
        lineEnds.second.second /= maxHeight;

        wireLine.SetLineColorAlpha(kBlue,shadeMinWire+((tempCharge-minWireCharge)/(maxWireCharge-minWireCharge))*(1.0-shadeMinWire));
        wireLine.DrawLine(gPad->GetUxmin()+(gPad->GetUxmax()-gPad->GetUxmin())*lineEnds.first.first,gPad->GetUymin()+(gPad->GetUymax()-gPad->GetUymin())*lineEnds.first.second,gPad->GetUxmin()+(gPad->GetUxmax()-gPad->GetUxmin())*lineEnds.second.first,gPad->GetUymin()+(gPad->GetUymax()-gPad->GetUymin())*lineEnds.second.second);
      }
    }
    for(Int_t i = 0; i < maxWiresV; i++)
    {
      tempCharge = cellMap.Vwires.at(i).charge;
      if(tempCharge > 0.0)
      {
        lineEnds = getWireEndpoints(i,kVPlane);
        lineEnds.first.first /= numYwires*wirePitchY-leftEdgeOffsetZval-rightEdgeOffsetZval;
        lineEnds.first.second /= maxHeight;
        lineEnds.second.first /= numYwires*wirePitchY-leftEdgeOffsetZval-rightEdgeOffsetZval;
        lineEnds.second.second /= maxHeight;
  
        wireLine.SetLineColorAlpha(kBlue,shadeMinWire+((tempCharge-minWireCharge)/(maxWireCharge-minWireCharge))*(1.0-shadeMinWire));
        wireLine.DrawLine(gPad->GetUxmin()+(gPad->GetUxmax()-gPad->GetUxmin())*lineEnds.first.first,gPad->GetUymin()+(gPad->GetUymax()-gPad->GetUymin())*lineEnds.first.second,gPad->GetUxmin()+(gPad->GetUxmax()-gPad->GetUxmin())*lineEnds.second.first,gPad->GetUymin()+(gPad->GetUymax()-gPad->GetUymin())*lineEnds.second.second);
      }
    }
    for(Int_t i = 0; i < maxWiresY; i++)
    {
      tempCharge = cellMap.Ywires.at(i).charge;
      if(tempCharge > 0.0)
      {
        lineEnds = getWireEndpoints(i,kYPlane);
        lineEnds.first.first /= numYwires*wirePitchY-leftEdgeOffsetZval-rightEdgeOffsetZval;
        lineEnds.first.second /= maxHeight;
        lineEnds.second.first /= numYwires*wirePitchY-leftEdgeOffsetZval-rightEdgeOffsetZval;
        lineEnds.second.second /= maxHeight;
        
        wireLine.SetLineColorAlpha(kBlue,shadeMinWire+((tempCharge-minWireCharge)/(maxWireCharge-minWireCharge))*(1.0-shadeMinWire));
        wireLine.DrawLine(gPad->GetUxmin()+(gPad->GetUxmax()-gPad->GetUxmin())*lineEnds.first.first,gPad->GetUymin()+(gPad->GetUymax()-gPad->GetUymin())*lineEnds.first.second,gPad->GetUxmin()+(gPad->GetUxmax()-gPad->GetUxmin())*lineEnds.second.first,gPad->GetUymin()+(gPad->GetUymax()-gPad->GetUymin())*lineEnds.second.second);
      }
    }
  }

  gPad->Update();
  gPad->RedrawAxis();

  TLine borderLine;
  borderLine.SetLineWidth(2.0);
  borderLine.DrawLine(gPad->GetUxmin(),gPad->GetUymax(),gPad->GetUxmax(),gPad->GetUymax());
  borderLine.DrawLine(gPad->GetUxmax(),gPad->GetUymin(),gPad->GetUxmax(),gPad->GetUymax());
  borderLine.DrawLine(gPad->GetUxmin(),gPad->GetUymin(),gPad->GetUxmax(),gPad->GetUymin());
  borderLine.DrawLine(gPad->GetUxmin(),gPad->GetUymin(),gPad->GetUxmin(),gPad->GetUymax());

  c1->SaveAs(Form("cellDiagram_%dUAngle_%dVAngle_%dYwires.png",(Int_t) round(angleU),(Int_t) round(angleV),numYwires));

  return;
}
