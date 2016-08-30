#include "stdafx.h"
#include "component/trianglecomponentinfo.h"
#include "component/trianglecomponent.h"

using namespace HydroCouple;

TriangleComponentInfo::TriangleComponentInfo(QObject *parent):
  ModelComponentInfo(parent)
{
  setId("Triangle Component 1.6");
  setCaption("Triangle ");
  setIconFilePath("./../../resources/images/superior.png");
  setDescription("A Two-Dimensional Quality Mesh Generator and Delaunay Triangulator.");
  setCategory("Computational Geometry");
  setCopyright("");
  setVendor("");
  setUrl("https://www.epa.gov/water-research/storm-water-management-model-swmm");
  setEmail("jrs@cs.berkley.gov");
  setVersion("1.6.0");

  QStringList publications;


  publications << "Shewchuk, J.R., 1996. Triangle: Engineering a 2D Quality Mesh Generator and Delaunay Triangulator."
                  "Applied Computational Geometry towards Geometric Engineering. Springer, pp. 203–222."
               << "Shewchuk, J.R., 2002. Delaunay Refinement Algorithms for Triangular Mesh Generation. Computational Geometry 22:21–74."
               << "Shewchuk, J.R., 2008. General-Dimensional Constrained Delaunay and Constrained Regular Triangulations, I: Combinatorial Properties."
                  "Discrete & Computational Geometry 39:580–637.";

  setPublications(publications);
}

TriangleComponentInfo::~TriangleComponentInfo()
{
}

IModelComponent *TriangleComponentInfo::createComponentInstance()
{
  QString id =  QUuid::createUuid().toString();
  TriangleComponent *component = new TriangleComponent(id, "Triangle Model Instance");
  component->setDescription("Triangle Model Instance");
  component->setComponentInfo(this);
  return component;
}
