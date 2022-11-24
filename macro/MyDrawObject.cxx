#include <MyDrawObject.h>

ClassImp(MyDrawObject);

MyDrawObject::MyDrawObject():
  axisX(),
  axisY(),
  axisXerror(),
  axisYerror(),
  axisSysYerror()
{  
}

void MyDrawObject::ClearObject(){
  axisX.clear();
  axisY.clear();
  axisXerror.clear();
  axisYerror.clear();
  axisSysYerror.clear();
}