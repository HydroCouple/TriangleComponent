#include "stdafx.h"
#include "component/trianglecomponent.h"
#include "core/argument1d.h"
#include "spatial/geometryargument.h"
#include "spatial/geometryfactory.h"
#include "core/idbasedargument.h"

#include  <iostream>

using namespace std;

int main(int argc, char *argv[])
{

  for(int i = 0 ; i < argc ; i++)
  {
    cout << argv[i] << endl;
  }

  if(argc == 3)
  {
    cout << "inside" << endl;
    //ned to flesh this out later
    cout << argv[1] << endl;
    cout << argv[2] << endl;

    TriangleComponent* triComp = new TriangleComponent("test");

    Argument1DString *triangleSwitches = triComp->commandSwitchesArgument();
    (*triangleSwitches)[0] = "pqc";

    QString dataFieldName;
    QString errorMessage;

    GeometryFactory::readGeometryDataItemFromFile(QString(argv[1]).trimmed(), dataFieldName, triComp->inputVerticesArgument(),errorMessage);

    (*triComp->outputOptionsArgument())["WriteOutputVerticesToCSV"] = "Yes";
    (*triComp->outputOptionsArgument())["OutputVerticeCSVFilePath"] = QString(argv[2]);

    (*triComp->thinningOptionsArgument())["ThinningMode"] = 2;

    triComp->initialize();
    triComp->prepare();
    triComp->update();
    triComp->finish();

    delete triComp;
  }

  return 0;
}
