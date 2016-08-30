#include "stdafx.h"
#include "component/trianglecomponent.h"
#include "core/dimension.h"
#include "core/valuedefinition.h"
#include "core/idbasedargument.h"
#include "core/argument1d.h"
#include "spatial/geometryargument.h"
#include "spatial/point.h"
#include "spatial/polyhedralsurface.h"
#include "spatial/polyhedralsurfaceexchangeitem.h"
#include "spatial/polygon.h"
#include "spatial/linestring.h"
#include "spatial/edge.h"
#include "spatial/geometryfactory.h"
#include "spatial/spatialreferencesystem.h"

#include <QString>
#include <assert.h>
#include <iostream>
#include <QDebug>
#include <math.h>
#include <QtConcurrent>
#include <set>

using namespace HydroCouple;
using namespace HydroCouple::Spatial;

TriangleComponent::TriangleComponent(const QString &id, TriangleComponent *parent):
  AbstractModelComponent(id,parent),
  m_points(nullptr),
  m_edges(nullptr),
  m_boundaries(nullptr),
  m_holes(nullptr),
  m_outputTriangulation(nullptr),
  m_identifiersArgument(nullptr),
  m_outputsArgument(nullptr),
  m_commandSwitches(nullptr),
  m_thinningArgument(nullptr),
  m_outputTIN(nullptr),
  m_maxNumberThinningIterations(1000),
  m_maxNumberOfPointsPerThinningIteration(1500),
  m_thinningCriteriaValue(0.99),
  m_thinningMode(ThinningMode::None),
  m_writePointsToFile(false),
  m_writePointsToCSV(false),
  m_writePatchsToFile(false),
  m_excludeOuterEdgesFromThinning(true),
  m_writeTriangulationToNetCDF(false),
  m_gdalVectorDriverName("ESRI Shapefile")
{
  nullifyTriangleIOObjectPtrs(&m_triangle);
  m_spatialReferenceSystem = new SpatialReferenceSystem(this);
  createArguments();
}

TriangleComponent::TriangleComponent(const QString &id, const QString &caption, TriangleComponent *parent):
  AbstractModelComponent(id,caption,parent),
  m_points(nullptr),
  m_edges(nullptr),
  m_boundaries(nullptr),
  m_holes(nullptr),
  m_outputTriangulation(nullptr),
  m_identifiersArgument(nullptr),
  m_outputsArgument(nullptr),
  m_commandSwitches(nullptr),
  m_thinningArgument(nullptr),
  m_outputTIN(nullptr),
  m_maxNumberThinningIterations(1000),
  m_maxNumberOfPointsPerThinningIteration(1500),
  m_thinningCriteriaValue(0.99),
  m_thinningMode(ThinningMode::None),
  m_writePointsToFile(false),
  m_writePointsToCSV(false),
  m_writePatchsToFile(false),
  m_excludeOuterEdgesFromThinning(true),
  m_writeTriangulationToNetCDF(false),
  m_gdalVectorDriverName("ESRI Shapefile")
{
  nullifyTriangleIOObjectPtrs(&m_triangle);
  m_spatialReferenceSystem = new SpatialReferenceSystem(this);
  createArguments();
}

TriangleComponent::~TriangleComponent()
{
  if(m_outputTIN)
    delete m_outputTIN;
}

void TriangleComponent::initialize()
{
  if(status() == HydroCouple::Created ||
     status() == HydroCouple::Initialized ||
     status() == HydroCouple::Failed)
  {
    setStatus(HydroCouple::Initializing , "Initialized Triangle Component Model");

    QString message;

    if(initializeIdentifierArgument(message) &&
       initializeCommandSwitchesArgument(message) &&
       initializeExchangeItems(message) &&
       initializeThinningArguments(message) &&
       initializeOutputsArguments(message))
    {

      setStatus(HydroCouple::Initialized , "Initialized Triangle Component Model");
      m_initialized = true;
    }
    else
    {
      setStatus(HydroCouple::Failed , message);
      m_initialized = false;
    }
  }
}

IModelComponent* TriangleComponent::clone()
{
  return nullptr;
}

QList<QString> TriangleComponent::validate()
{
  int numPoints = m_points->dimensionLength(nullptr,0);

  QList<HCGeometry*> edges = m_edges->geometries();

  for(int i = 0; i < edges.length() ; i++)
  {
    HCLineString* lineString = dynamic_cast<HCLineString*>(edges[i]);
    numPoints += lineString->pointCount();
  }

  QList<HCGeometry*> boundaries = m_boundaries->geometries();
  for(int i = 0; i < boundaries.length() ; i++)
  {
    HCPolygon* polygon = dynamic_cast<HCPolygon*>(boundaries[i]);
    numPoints += polygon->exteriorRing()->pointCount();
  }

  if(m_initialized && numPoints)
  {
    setStatus(HydroCouple::Validating,"Validating Triangle Component");


    setStatus(HydroCouple::Valid,"Triangle Component Is Valid");
  }
  else
  {
    setStatus(HydroCouple::Failed,"Triangle Component Is Valid. Inputs must have at least on point");
    return QList<QString>({"Triangle Component Is Valid. Inputs must have at least on point"});
  }

  return QList<QString>();
}

void TriangleComponent::prepare()
{
  if(m_initialized)
  {
    setStatus(HydroCouple::Preparing , "Preparing Triangle Model");

    createTriangleIOObject();

    setStatus(HydroCouple::Updated ,"Finished preparing Triangle Model");
    m_prepared = true;
  }
  else
  {
    m_prepared = false;
  }
}

void TriangleComponent::update(const QList<IOutput *> &requiredOutputs)
{
  if(status() == HydroCouple::Updated)
  {
    setStatus(HydroCouple::Updating , "Triangulation Component " + id() + " is performing triangulation");

    triangulateio outputTriangulateio, vorout;

    nullifyTriangleIOObjectPtrs(&outputTriangulateio);
    nullifyTriangleIOObjectPtrs(&vorout);

    if(!m_switches.contains("z"))
      m_switches.append('z');

    char* switches = (char*) m_switches.toStdString().c_str();
    triangulate(switches,&m_triangle,&outputTriangulateio,&vorout);

    setStatus(HydroCouple::Updating , "Triangulation Component " + id() + " has finished performing triangulation");
    setStatus(HydroCouple::Updating , "Triangulation component " + id() + " is copying triangulation outputs");

    if(m_outputTIN)
      delete m_outputTIN;

    m_outputTIN = readTINData(&outputTriangulateio);
    m_outputTIN->spatialReferenceSystemInternal()->setSrText(m_spatialReferenceSystem->srText());

    freeTriangleIOObject(&m_triangle);
    //because pointer is the same as triangle
    outputTriangulateio.holelist = NULL;
    outputTriangulateio.regionlist = NULL;
    freeTriangleIOObject(&outputTriangulateio);
    freeTriangleIOObject(&vorout);

    if(m_thinningMode != ThinningMode::None)
    {
      setStatus(HydroCouple::Updating , "Triangulation component " + id() + " is performing thinning operation");

      int pointCountBefore = m_outputTIN->vertexCount();
      int numPointsRemoved = performThinning();

      QString message;
      message   = "Thinning Operation Results: " + QString::number(numPointsRemoved) + " out of " +  QString::number(pointCountBefore)
                  +" vertices removed (" + QString::number((numPointsRemoved*100.0)/(pointCountBefore*1.0)) +" %)" ;

      setStatus(HydroCouple::Updating , "Triangulation component " + id() + " finished performing thinning operation");
      setStatus(HydroCouple::Updating , "Triangulation component " + id() + " Summary: " + message);


    }

    updateOutputExchangeItems();

    setStatus(HydroCouple::Updating , "Triangulation Component " + id() + " is done performing triangulation");
    setStatus(HydroCouple::Updating , "Triangulation Component " + id() + " Triangles: " + QString::number(m_outputTIN->patchCount()));
    setStatus(HydroCouple::Done , "Triangulation Component " + id() + " Points: " + QString::number(m_outputTIN->vertexCount()));
  }
}

void TriangleComponent::finish()
{
  if(m_prepared)
  {
    setStatus(HydroCouple::Finishing , "Triangle component simulation with component id " + id() + " is being disposed");

    writeOutput();

    m_prepared = false;
    m_initialized = false;

    setStatus(HydroCouple::Finished , "Triangle component  simulation with component id " + id() + " has been disposed");
    setStatus(HydroCouple::Created , "Triangle component  simulation with component id " + id() + " ran successfully and has been re-created");

  }
}

IdBasedArgumentQString *TriangleComponent::identifiersArgument() const
{
  return m_identifiersArgument;
}

HCGeometryArgumentDouble *TriangleComponent::inputVerticesArgument() const
{
  return m_points;
}

HCGeometryArgumentDouble *TriangleComponent::inputEdgesArgument() const
{
  return m_edges;
}

HCGeometryArgumentDouble *TriangleComponent::inputBoundariesArgument() const
{
  return m_edges;
}

HCGeometryArgumentDouble *TriangleComponent::inputHolesArgument() const
{
  return m_holes;
}

Argument1DString *TriangleComponent::commandSwitchesArgument() const
{
   return m_commandSwitches;
}

IdBasedArgumentDouble *TriangleComponent::thinningOptionsArgument() const
{
  return m_thinningArgument;
}

IdBasedArgumentQString *TriangleComponent::outputOptionsArgument() const
{
  return m_outputsArgument;
}

void TriangleComponent::createArguments()
{
  createIdentifierArgument();
  createCommandsSwitchesArgument();
  createGeometriesArguments();
  createThinningArguments();
  createOutputFileArgument();
}

void TriangleComponent::createIdentifierArgument()
{
  Dimension *identifierDimension = new Dimension("IdentifierDimension","Dimension for identifiers", this);
  QStringList identifiers;
  identifiers.append("Id");
  identifiers.append("Caption");
  identifiers.append("Description");
  Quantity* quantity = Quantity::unitLessValues("IdentifiersQuantity","", QVariant::String , this);

  m_identifiersArgument = new IdBasedArgumentQString("Identifiers", identifiers,identifierDimension,quantity,this);
  m_identifiersArgument->setCaption("Model Identifiers");

  (*m_identifiersArgument)["Id"] = id();
  (*m_identifiersArgument)["Caption"] = caption();
  (*m_identifiersArgument)["Description"] = description();


  m_identifiersArgument->addInputFileTypeFilter("Input XML File (*.xml)");
  m_identifiersArgument->setMatchIdentifiersWhenReading(true);

  addArgument(m_identifiersArgument);
}

void TriangleComponent::createCommandsSwitchesArgument()
{
  QStringList comments;

  Quantity *switchesQuantity = Quantity::unitLessValues("TriangleSwitchesQuantity","", QVariant::String , this);
  Dimension *switchesDimension = new Dimension("TriangleSwitchesDimension","", this);
  comments.append("The command switches to use for triangulation. See "
                  "https://www.cs.cmu.edu/~quake/triangle.help.html for details.");

  m_commandSwitches = new Argument1DString("TriangleSwitches",switchesDimension,1,switchesQuantity,this);
  m_commandSwitches->setCaption("Triangle Component Switches");
  m_commandSwitches->setDescription("The command switches to use for triangulation. See "
                                    "https://www.cs.cmu.edu/~quake/triangle.help.html for details.");

  m_commandSwitches->addInputFileTypeFilter("Input XML File (*.xml)");
  m_commandSwitches->setValueT(0,"q");
  m_commandSwitches->setComments(comments);

  addArgument(m_commandSwitches);

}

void TriangleComponent::createGeometriesArguments()
{
  Dimension *inputGeometryDimension = new Dimension("TnputGeometryDimension","", this);
  Quantity *inputGeometryQuantity = Quantity::unitLessValues("TnputGeometryQuantity","", QVariant::String,this);

  m_points = new HCGeometryArgumentDouble("TrigulationPoints", GeometryType::PointZ,inputGeometryDimension,inputGeometryQuantity,this);
  m_points->setCaption("Triangulation Points");
  m_points->setDescription("Input points to be used for triangulation");
  m_points->addGeometry(new HCPoint(56,67,67,m_points));
  m_points->resetDataArray();
  addArgument(m_points);

  m_edges = new HCGeometryArgumentDouble("ConstrainingEdges", GeometryType::LineStringZ,inputGeometryDimension,inputGeometryQuantity,this);
  m_edges->setCaption("Triangulation Constraining Edges");
  m_edges->setDescription("Input edges to be used for triangulation");
  m_edges->resetDataArray();
  addArgument(m_edges);

  m_boundaries = new HCGeometryArgumentDouble("BoundaryPolygonEdges", GeometryType::PolygonZ,inputGeometryDimension,inputGeometryQuantity,this);
  m_boundaries->setCaption("Triangulation Boundary Polygons");
  m_boundaries->setDescription("Input boundaries to be used for for constraining triangulation and creating holes");
  m_boundaries->resetDataArray();
  addArgument(m_boundaries);

  m_holes = new HCGeometryArgumentDouble("HoleLocations", GeometryType::PointZ,inputGeometryDimension,inputGeometryQuantity,this);
  m_holes->setCaption("Triangulation Hole Locations");
  m_holes->setDescription("Input locations to be used for holes in triangulation. These must be bounded by a boundary polygon");
  m_holes->resetDataArray();
  addArgument(m_holes);
}

void TriangleComponent::createThinningArguments()
{
  QStringList comments;
  Dimension *thinningDimension = new Dimension("ThinningOptionDimension","Dimension for thinning options", this);
  QStringList identifiers;

  comments.append("ThinningMode: Not 1 or 2 = None thinning, 1 = VolumeChangeCriteria, 2 = NormalVectorCriteria");

  identifiers.append("ThinningMode");

  comments.append("CriteriaValue: For the VolumeChangeCriteria option this value represents the minimum fraction of volume change below which surface a surface is considered flat. "
                  "For the NormalVectorCriteria option a number close to 1 signifies flat");

  identifiers.append("CriteriaValue");


  comments.append("MaxNumberThinningIterations: Maximum number of iterations per thinning cycle. "
                  "If not specified, the thinning process repeats until there is no eligible vertex candidate for removal");
  identifiers.append("MaxNumberThinningIterations");

  comments.append("MaxNumberPointsPerThinningIteration: Maximum number of eligible points to sample per iterations");
  identifiers.append("MaxNumberPointsPerThinningIteration");


  identifiers.append("ExcludeOuterEdges");

  Quantity* quantity = Quantity::unitLessValues("ThinningOptionQuantity","", QVariant::Double , this);

  m_thinningArgument = new IdBasedArgumentDouble("ThinningOptionArgument", identifiers,thinningDimension,quantity,this);
  m_thinningArgument->setCaption("Thinning Options Argument");

  (*m_thinningArgument)["ThinningMode"] = 0;
  (*m_thinningArgument)["CriteriaValue"] = m_thinningCriteriaValue;
  (*m_thinningArgument)["MaxNumberThinningIterations"] = m_maxNumberThinningIterations;
  (*m_thinningArgument)["MaxNumberPointsPerThinningIteration"] = m_maxNumberOfPointsPerThinningIteration;
  (*m_thinningArgument)["ExcludeOuterEdges"] = m_excludeOuterEdgesFromThinning ? 1:0;


  m_thinningArgument->addInputFileTypeFilter("Input XML File (*.xml)");
  m_thinningArgument->setMatchIdentifiersWhenReading(true);
  m_thinningArgument->setComments(comments);

  addArgument(m_thinningArgument);
}

void TriangleComponent::createOutputFileArgument()
{
  Dimension *outputFileDimension = new Dimension("OutputFileDimension","Dimension for outputFile", this);
  Quantity* quantity = Quantity::unitLessValues("OutputFileQuality","", QVariant::String , this);


  QStringList options;
  QStringList comments;
  options.append("WriteOutput");
  comments << "WriteOutput : Yes, True, or 1. Everything else treated as false. Valid path for NetCDF file must be specified";
  options.append("NetCDFOutputFilePath");


  options.append("GDALDriverName");
  comments << "GDALDriverName : The name of the GDAL driver to use to write vector outputs. Default is ESRI Shapefile";

  options.append("WriteOutputVertices");
  comments << "WriteOutputVertices : Yes, True, or 1. Everything else treated as false. Valid GDALDriverName must be specified. Valid OutputVerticesFilePath must also be specified";
  options.append("OutputVerticesFilePath");

  options.append("WriteOutputTriangles");
  comments << "WriteOutputTriangles : Yes, True, or 1. Everything else treated as false. Valid GDALDriverName must be specified. Valid OutputTrianglesFilePath must also be specified";
  options.append("OutputTrianglesFilePath");


  options.append("WriteOutputVerticesToCSV");
  comments << "WriteOutputVerticesToCSV : Yes, True, or 1. Everything else treated as false. Valid OutputVerticeCSVFilePath must also be specified";
  options.append("OutputVerticeCSVFilePath");


  m_outputsArgument = new IdBasedArgumentQString("OutputOptions", options,outputFileDimension,quantity,this);
  m_outputsArgument->setCaption("Output Options");

  (*m_outputsArgument)["WriteOutput"] = m_writeTriangulationToNetCDF ? "Yes" : "No";
  (*m_outputsArgument)["GDALDriverName"] = m_gdalVectorDriverName;
  (*m_outputsArgument)["WriteOutputVertices"] = m_writePointsToFile ? "Yes" : "No";
  (*m_outputsArgument)["WriteOutputTriangles"] = m_writePatchsToFile ? "Yes" : "No";
  (*m_outputsArgument)["WriteOutputVerticesToCSV"] = m_writePointsToCSV ? "Yes" : "No";


  m_outputsArgument->setComments(comments);
  addArgument(m_outputsArgument);

}

bool TriangleComponent::initializeIdentifierArgument(QString &message)
{

  QString identifier = (*m_identifiersArgument)["Id"];

  if(identifier.isNull() || identifier.isEmpty())
  {
    message = "The id provided is invalid!";
    return false;
  }
  else
  {
    setId(identifier);
  }

  setCaption((*m_identifiersArgument)["Caption"]);
  setDescription((*m_identifiersArgument)["Description"]);

  return true;
}

bool TriangleComponent::initializeCommandSwitchesArgument(QString &message)
{
  if(m_commandSwitches->dimensionLength(nullptr,0))
  {
    int index[] = {0};
    QVariant switches;
    m_commandSwitches->getValue(index,switches);
    m_switches = switches.toString();

    if(m_switches.isEmpty() || m_switches.isNull())
    {
      message = "No command switches have been specified";
      return false;
    }
  }
  else
  {
    message  = "No command switches have been specified";
    return false;
  }

  return true;
}

bool TriangleComponent::initializeExchangeItems(QString &message)
{
  clearInputExchangeItems();
  clearOutputExchangeItems();

  if(m_outputTIN)
  {
    delete m_outputTIN;
    m_outputTIN = nullptr;
  }

  m_outputTIN = new HCTIN(this);
  m_outputTIN->enable3D();

  Dimension *patchDim = new Dimension("PatchDimension","", this);
  Dimension *edgeDim = new Dimension("EdgeDimension","", this);
  Dimension *nodeDim = new Dimension("NodeDimension","", this);

  Quantity *tinQuality = Quantity::unitLessValues("OutputGeometryQuantity","", QVariant::String,this);


  m_outputTriangulation = new HCTINOutputDouble("OutputTin",PolyhedralSurfaceDataType::Patch,m_outputTIN,
                                                patchDim,edgeDim,nodeDim,tinQuality,this);

  m_outputTriangulation->setCaption("Output TIN Surface");

  addOutputExchangeItem(m_outputTriangulation);

  return true;
}

bool TriangleComponent::initializeThinningArguments(QString &message)
{
  m_maxNumberThinningIterations = 100;
  m_maxNumberOfPointsPerThinningIteration = 100;


  double  value = (*m_thinningArgument)["ThinningMode"];

  if(value  == 1.0)
  {
    m_thinningMode = ThinningMode::VolumeChangeCriteria;
    m_thinningCriteriaValue  = 0.1;
  }
  else if (value == 2.0)
  {
    m_thinningMode = ThinningMode::NormalVectorCriteria;
    m_thinningCriteriaValue = 0.95;
  }
  else
  {
    m_thinningMode = ThinningMode::None;
  }


  if(m_thinningMode != ThinningMode::None)
  {
    value = (*m_thinningArgument)["CriteriaValue"];

    if(value > 0 && value < 1.0)
    {
      m_thinningCriteriaValue = value;
    }

    value = (*m_thinningArgument)["MaxNumberThinningIterations"];

    if(value > 0 )
      m_maxNumberThinningIterations = (int)value;

    value = (*m_thinningArgument)["MaxNumberPointsPerThinningIteration"];

    if(value > 0 )
    {
      m_maxNumberOfPointsPerThinningIteration = (int)value;
    }

    value = (*m_thinningArgument)["ExcludeOuterEdges"];
    m_excludeOuterEdgesFromThinning = value == 0 ? false : true;
  }

  return true;
}

bool TriangleComponent::initializeOutputsArguments(QString &message)
{
  QString value = "";
  bool canConvertToDouble = false;
  int stride =1;

  int index = m_outputsArgument->identifiers().indexOf("WriteOutput");

  if(index > -1)
  {
    m_outputsArgument->getValues(&index, &stride , &value);

    if((m_writeTriangulationToNetCDF = getBoolFromString(value)))
    {
      index = m_outputsArgument->identifiers().indexOf("NetCDFOutputFilePath");

      if(index > -1)
      {
        m_outputsArgument->getValues(&index, &stride , &value);

        if(value.isEmpty() || value.isNull() || !(m_netCDFFile = getRelativeFilePath(value)).absoluteDir().exists())
        {
          message = "NetCDF output file directory does not exist: " + value;;
          return false;
        }
      }
    }
  }

  index = m_outputsArgument->identifiers().indexOf("WriteOutputVertices");
  if(index > -1)
  {
    m_outputsArgument->getValues(&index, &stride , &value);

    if((m_writePointsToFile = getBoolFromString(value)))
    {
      index = m_outputsArgument->identifiers().indexOf("OutputVerticesFilePath");

      if(index > -1)
      {
        m_outputsArgument->getValues(&index, &stride , &value);

        if(value.isEmpty() || value.isNull() || !(m_pointsFile = getRelativeFilePath(value)).absoluteDir().exists())
        {
          message = "Points output file directory does not exist: " + value;;
          return false;
        }
      }
    }
  }

  index = m_outputsArgument->identifiers().indexOf("WriteOutputTriangles");
  if(index > -1)
  {
    m_outputsArgument->getValues(&index, &stride , &value);

    if((m_writePatchsToFile = getBoolFromString(value)))
    {
      index = m_outputsArgument->identifiers().indexOf("OutputTrianglesFilePath");

      if(index > -1)
      {
        m_outputsArgument->getValues(&index, &stride , &value);

        if(value.isEmpty() || value.isNull() || !(m_patchesFile = getRelativeFilePath(value)).absoluteDir().exists())
        {
          message = "Polygon output file directory does not exist: " + value;;
          return false;
        }
      }
    }
  }

  index = m_outputsArgument->identifiers().indexOf("WriteOutputVerticesToCSV");
  if(index > -1)
  {
    m_outputsArgument->getValues(&index, &stride , &value);

    if((m_writePointsToCSV = getBoolFromString(value)))
    {
      index = m_outputsArgument->identifiers().indexOf("OutputVerticeCSVFilePath");

      if(index > -1)
      {
        m_outputsArgument->getValues(&index, &stride , &value);

        if(value.isEmpty() || value.isNull() || !(m_CSVFile = getRelativeFilePath(value)).absoluteDir().exists())
        {
          message = "CSV output file directory does not exist: " + value;;
          return false;
        }
      }
    }
  }

  return true;
}

void TriangleComponent::createTriangleIOObject()
{
  freeTriangleIOObject(&m_triangle);

  //check to make sure correct flags are used
  int numPoints = m_points->dimensionLength(nullptr,0);
  int numEdges = 0;

  if(m_points->geometryCount())
  {
    m_spatialReferenceSystem->setSrText(m_points->geometry(0)->spatialReferenceSystem()->srText());
  }
  else if(m_boundaries->geometryCount())
  {
    m_spatialReferenceSystem->setSrText(m_boundaries->geometry(0)->spatialReferenceSystem()->srText());

  }
  else if(m_edges->geometryCount())
  {
    m_spatialReferenceSystem->setSrText(m_edges->geometry(0)->spatialReferenceSystem()->srText());

  }
  else if(m_holes->geometryCount())
  {
    m_spatialReferenceSystem->setSrText(m_holes->geometry(0)->spatialReferenceSystem()->srText());

  }


  QList<HCGeometry*> edges = m_edges->geometries();
  for(int i = 0; i < edges.length() ; i++)
  {
    HCLineString* lineString = dynamic_cast<HCLineString*>(edges[i]);
    numPoints += lineString->pointCount();
    numEdges += lineString->pointCount() - 1;
  }

  QList<HCGeometry*> boundaries = m_boundaries->geometries();
  for(int i = 0; i < boundaries.length() ; i++)
  {
    HCPolygon* polygon = dynamic_cast<HCPolygon*>(boundaries[i]);
    numPoints += polygon->exteriorRing()->pointCount() - 1;
    numEdges += polygon->exteriorRing()->pointCount() - 1;
  }

  m_triangle.numberofpoints = numPoints;
  m_triangle.numberofpointattributes = 1;
  m_triangle.pointlist = (REAL*) malloc(numPoints * 2 * sizeof(REAL));
  m_triangle.pointmarkerlist = (int*) malloc(numPoints * sizeof(int));
  m_triangle.pointattributelist =   (REAL*) malloc(numPoints * sizeof(REAL));
  m_triangle.numberofholes = m_holes->geometries().length();
  m_triangle.holelist = (REAL*) malloc(m_holes->geometries().length() * 2 * sizeof(REAL));
  m_triangle.numberofsegments = numEdges;
  m_triangle.segmentmarkerlist = (int*) malloc(numEdges * sizeof(int));
  m_triangle.segmentlist = (int*) malloc(numEdges * 2 * sizeof(int));

  int currPoint = 0;
  int currentMarker = 0;
  int currEdge = 0;

  QList<HCGeometry*> points = m_points->geometries();
  for(int i = 0 ; i < points.length() ; i++)
  {
    int pindex = currPoint*2;

    HCPoint *p = dynamic_cast<HCPoint*>(points[i]);
    m_triangle.pointlist[pindex] = p->x();
    m_triangle.pointlist[pindex+1] = p->y();
    m_triangle.pointmarkerlist[currPoint] = currentMarker;
    m_triangle.pointattributelist[currPoint] = p->z();

    currPoint++;
  }

  for(int i = 0; i < edges.length() ; i++)
  {
    HCLineString* lineString = dynamic_cast<HCLineString*>(edges[i]);
    currentMarker++;

    for(int j = 0; j < lineString->pointCount() ; j++)
    {
      int pindex = currPoint*2;

      HCPoint *p = dynamic_cast<HCPoint*>(lineString->point(j));
      m_triangle.pointlist[pindex] = p->x();
      m_triangle.pointlist[pindex+1] = p->y();
      m_triangle.pointmarkerlist[currPoint] = currentMarker;
      m_triangle.pointattributelist[currPoint] = p->z();

      if(j < lineString->pointCount() - 1)
      {
        int edgeIndex = 2 * currEdge;
        m_triangle.segmentlist[edgeIndex] = currPoint;
        m_triangle.segmentlist[edgeIndex+1] = currPoint+1;
        m_triangle.segmentmarkerlist[currEdge] = currentMarker;
        currEdge++;
      }


      currPoint++;
    }

  }

  for(int i = 0; i < boundaries.length() ; i++)
  {
    HCPolygon* polygon = dynamic_cast<HCPolygon*>(boundaries[i]);
    HCLineString* lineString = polygon->exteriorRingInternal();
    currentMarker++;

    int initPoint = currPoint;

    for(int j = 0; j < lineString->pointCount() - 1 ; j++)
    {
      int pindex = currPoint*2;

      HCPoint *p = dynamic_cast<HCPoint*>(lineString->point(j));
      m_triangle.pointlist[pindex] = p->x();
      m_triangle.pointlist[pindex+1] = p->y();
      m_triangle.pointmarkerlist[currPoint] = currentMarker;
      m_triangle.pointattributelist[currPoint] = p->z();


      int edgeIndex = 2 * currEdge;
      m_triangle.segmentlist[edgeIndex] = currPoint;

      if(j == lineString->pointCount() - 2)
      {
        m_triangle.segmentlist[edgeIndex+1] = initPoint ;
      }
      else
      {
        m_triangle.segmentlist[edgeIndex+1] = currPoint+1 ;
      }

      m_triangle.segmentmarkerlist[currEdge] = currentMarker;

      currEdge++;
      currPoint++;

    }
  }

  QList<HCGeometry*> holes = m_holes->geometries();

  for(int i = 0; i < holes.length() ; i++)
  {
    HCPoint *point = dynamic_cast<HCPoint*>(holes[i]);

    m_triangle.holelist[i*2] = point->x();
    m_triangle.holelist[i*2 + 1] = point->y();
  }

}

void TriangleComponent::nullifyTriangleIOObjectPtrs(triangulateio *triangleObject)
{
  triangleObject->pointlist = NULL;
  triangleObject->pointattributelist = NULL;
  triangleObject->pointmarkerlist = NULL;
  triangleObject->trianglelist = NULL;
  triangleObject->triangleattributelist = NULL;
  triangleObject->trianglearealist = NULL;
  triangleObject->neighborlist = NULL;
  triangleObject->segmentlist = NULL;
  triangleObject->segmentmarkerlist = NULL;
  triangleObject->holelist = NULL;
  triangleObject->regionlist = NULL;
  triangleObject->edgelist = NULL;
  triangleObject->edgemarkerlist = NULL;

  triangleObject->numberofpoints = 0;
  triangleObject->numberofpointattributes = 0;
  triangleObject->numberoftriangles = 0;
  triangleObject->numberofcorners = 0;
  triangleObject->numberoftriangleattributes = 0;
  triangleObject->numberofsegments = 0;
  triangleObject->numberofholes = 0;
  triangleObject->numberofregions = 0;
  triangleObject->numberofedges = 0;
}

void TriangleComponent::freeTriangleIOObject(triangulateio *triangleObject)
{
  if(triangleObject->pointlist)
    free(triangleObject->pointlist);

  if(triangleObject->pointattributelist)
    free(triangleObject->pointattributelist);

  if(triangleObject->pointmarkerlist)
    free(triangleObject->pointmarkerlist);

  if(triangleObject->trianglelist)
    free(triangleObject->trianglelist);

  if(triangleObject->triangleattributelist)
    free(triangleObject->triangleattributelist);

  if(triangleObject->trianglearealist)
    free(triangleObject->trianglearealist);

  if(triangleObject->neighborlist)
    free(triangleObject->neighborlist);

  if(triangleObject->segmentlist)
    free(triangleObject->segmentlist);

  if(triangleObject->segmentmarkerlist)
    free(triangleObject->segmentmarkerlist);

  if(triangleObject->holelist)
    free(triangleObject->holelist);

  if(triangleObject->regionlist)
    free(triangleObject->regionlist);

  if(triangleObject->edgelist)
    free(triangleObject->edgelist);

  if(triangleObject->edgemarkerlist)
    free(triangleObject->edgemarkerlist);

  nullifyTriangleIOObjectPtrs(triangleObject);
}

HCTIN *TriangleComponent::readTINData(triangulateio *triangulation)
{
  HCTIN *tin = new HCTIN(nullptr);
  tin->enable3D();

  for(int i = 0 ; i < triangulation->numberofpoints ; i++)
  {
    int index = 2 *i;

    double x = triangulation->pointlist[index];
    double y = triangulation->pointlist[index + 1];
    double z = triangulation->pointattributelist[i];
    HCVertex* vertex = new HCVertex(x,y,z,tin);
  }

  QVector<HCVertex*> addedVertices = tin->vertices();

  assert(triangulation->numberofcorners == 3);

  //  QFutureSynchronizer<HCTriangle*> synchronizer;

  for(int i = 0 ; i < triangulation->numberoftriangles ; i++)
  {
    int index  = i * 3;

    HCVertex *v1 = addedVertices[triangulation->trianglelist[index]];
    HCVertex *v2 = addedVertices[triangulation->trianglelist[index+1]];
    HCVertex *v3 = addedVertices[triangulation->trianglelist[index+2]];

    //    std::set<double> fix;
    //    fix.insert(v1->z());
    //    fix.insert(v2->z());
    //    fix.insert(v3->z());

    //    std::set<double>::iterator it = fix.begin();
    //    std::advance(it,1);
    //    double mid = *it;

    //    if( fabs(v1->z() - mid)/mid > 0.1)
    //    {
    //      v1->setZ(mid);
    //    }

    //    if( fabs(v2->z() - mid)/mid > 0.1)
    //    {
    //      v2->setZ(mid);
    //    }

    //    if( fabs(v3->z() - mid)/mid > 0.1)
    //    {
    //      v3->setZ(mid);
    //    }

    HCTriangle *triangle = tin->createTriangle(v1,v2,v3);
    //    QFuture<HCTriangle*> future = QtConcurrent::run(tin,&HCTIN::createTriangle,v1,v2,v3 );
    //    synchronizer.addFuture(future);
  }

  //  synchronizer.waitForFinished();

  return tin;
}

int TriangleComponent::performThinning()
{
  int numberOfThinningOps = 0;


  if(m_outputTIN->patchCount() && m_thinningMode != ThinningMode::None)
  {
    int iterCount = 0;

    while(true)
    {
      if(iterCount >= m_maxNumberThinningIterations)
      {
        break;
      }
      else
      {

        QVector<HCVertex*> vertices = m_outputTIN->vertices();
        int vlength = vertices.length();

        //reset markers
        for(int i = 0 ; i < vlength ; i++)
        {
          vertices[i]->setMarker(0);
        }

        int eligibleCount = 0;
        int numTriangles = 0;

        for(int i = 0 ; i < vertices.length() ; i++)
        {
          HCVertex *vertex = vertices[i];

          if(eligibleCount >= m_maxNumberOfPointsPerThinningIteration)
          {
            break;
          }
          else if(isEligibleForThinning(vertex,numTriangles))
          {
            assert(numTriangles > 1);

            vertex->setMarker(numTriangles);
            setOrbitDestinationMarkers(vertex,1);
            eligibleCount++;
          }
        }

        if(eligibleCount == 0)
        {
          break;
        }


        QFutureSynchronizer<bool> synchronizer;

        for(int i = 0 ; i < vlength ; i++)
        {
          HCVertex *vertex = vertices[i];

          if(vertex->getMarker() > 1)
          {
            QFuture<bool> future = QtConcurrent::run(this,&TriangleComponent::performThinningOnVertex, vertex, m_thinningMode, m_thinningCriteriaValue);
            synchronizer.addFuture(future);
          }
        }

        synchronizer.waitForFinished();

        for(QFuture<bool> future : synchronizer.futures())
        {
          if(future.result())
          {
            numberOfThinningOps++;
          }
        }

        iterCount++;
      }
    }
  }


  return numberOfThinningOps;
}

void TriangleComponent::updateInputExchangeItems()
{

}

void TriangleComponent::updateOutputExchangeItems()
{
  m_outputTriangulation->setPolyhedralSurface(m_outputTIN);

  for(AbstractOutput *output : outputsInternal())
  {
    for(HydroCouple::IAdaptedOutput* adaptedOutput : output->adaptedOutputs())
    {
      adaptedOutput->refresh();
    }
  }
}

void TriangleComponent::writeOutput()
{
  writeOutputVerticesToCSV();

  if(m_writePointsToFile && m_pointsFile.absoluteDir().exists())
  {
    GeometryFactory::writeTINVertices(m_outputTIN,m_pointsFile.absoluteFilePath(),m_gdalVectorDriverName);
  }

  if(m_writePatchsToFile && m_patchesFile.absoluteDir().exists())
  {
    GeometryFactory::writeTINPolygons(m_outputTIN,m_patchesFile.absoluteFilePath(),m_gdalVectorDriverName);
  }

  if(m_writeTriangulationToNetCDF && m_netCDFFile.absoluteDir().exists())
  {
    QString message;
    GeometryFactory::writeTINToFile(m_outputTIN, m_netCDFFile.absoluteFilePath() ,message);
  }
}

void TriangleComponent::writeOutputVerticesToCSV()
{
  if(m_writePointsToCSV && m_CSVFile.absoluteDir().exists())
  {
    QFile file(m_CSVFile.absoluteFilePath());

    if(file.open(QIODevice::WriteOnly | QIODevice::Append))
    {
      QTextStream writer(&file);

      QVector<HCVertex*> vertices = m_outputTIN->vertices();

      for(int i = 0 ; i < m_outputTIN->vertexCount() ; i++)
      {
        HCVertex *v = vertices[i];
        writer << QString::number(v->x(), 'g' , 20) << "," << QString::number(v->y(), 'g' , 20) << "," << QString::number(v->z(), 'g' , 20) << endl;
      }


      file.flush();
      file.close();
    }
  }
}

bool TriangleComponent::isEligibleForThinning(HCVertex *vertex, int &numTriangles)
{
  numTriangles = 0;

  if(vertex->getMarker() == 0 && vertex->edge())
  {
    HCEdge *edge = vertex->edgeInternal();
    HCEdge *endEdge = edge;

    //check all edges have faces and neighbour orbit edges are internal
    do
    {
      if(edge->leftInternal() == nullptr ||
         edge->rightInternal() == nullptr ||
         //edge->destInternal()->getMarker() //||
         (m_excludeOuterEdgesFromThinning && edge->leftNextInternal()->rightInternal() == nullptr)
         )
      {
        return false;
      }

      numTriangles++;
      edge = edge->origNextInternal();

    }while(edge != endEdge) ;

    if(numTriangles > 2)
      return true;
  }

  return false;
}

bool TriangleComponent::performThinningOnVertex(HCVertex *vertex, ThinningMode thinningMode, double criteriaValue)
{
  if(thinningMode == ThinningMode::NormalVectorCriteria)
  {
    //for each triangle calculate normal
    int numTri = vertex->getMarker();

    double normals[numTri*3];
    double aveN[3] ={0};

    HCEdge *edge = vertex->edgeInternal();

    for(int i = 0 ; i < numTri; i++)
    {
      int index = i*3;

      HCVertex *p1 = edge->origInternal();
      HCVertex *p2 = edge->destInternal();
      HCVertex *p3 = edge->leftNextInternal()->destInternal();

      triangleNormalVector(p1,p2,p3, &normals[index]);

      normalizeVector(&normals[index]);

      aveN[0] += normals[index];
      aveN[1] += normals[index+1];
      aveN[2] += normals[index+2];

      edge = edge->origNextInternal();
    }

    aveN[0] = aveN[0] / (numTri*1.0);
    aveN[1] = aveN[1] / (numTri*1.0);
    aveN[2] = aveN[2] / (numTri*1.0);

    double mean = 0;

    for(int i = 0 ; i < numTri; i++)
    {
      int index = i*3;

      mean += dotProduct(aveN,&normals[index]);
    }

    mean = mean / (numTri*1.0);

    if(mean > criteriaValue)
    {
      removeAndTriangulateOrbit(vertex);
      return true;
    }
    else
    {
      return false;
    }
  }
  else
  {

  }

  return false;
}

void TriangleComponent::removeAndTriangulateOrbit(HCVertex *vertex)
{
  assert(m_outputTIN);

  int numTri = vertex->getMarker();
  HCEdge *edge = vertex->edgeInternal();

  //points to triangulate after vertex removal

  QVector<HCVertex*> points(numTri);

  //to expensive find workaround
  {
    for(int i = 0 ; i < numTri; i++)
    {
      points[i] = edge->leftNextInternal()->origInternal();
      edge = edge->origNextInternal();
    }

    m_mutex.lock();
    //delete vertex
    while (vertex->edgeInternal())
    {

      m_outputTIN->deleteFaceEdge(vertex->edgeInternal());
    }

    //delete dangling vertex
    delete vertex;

    //delete orbit polygon so that new triangulation can be inserted.
    delete points[0]->edgeInternal()->leftInternal();

    m_mutex.unlock();

  }

  triangulateio in;
  triangulateio out;
  triangulateio vrout;

  nullifyTriangleIOObjectPtrs(&in);
  nullifyTriangleIOObjectPtrs(&out);
  nullifyTriangleIOObjectPtrs(&vrout);

  in.numberofpoints = numTri;
  in.pointlist = (REAL*) malloc(numTri * 2 * sizeof(REAL));
  in.pointattributelist = (REAL*) malloc(numTri * sizeof(REAL));
  in.pointmarkerlist = (int*) malloc(numTri * sizeof(REAL));
  in.numberofsegments = numTri;
  in.segmentmarkerlist = (int*) malloc(numTri * sizeof(REAL));
  in.segmentlist = (int*) malloc(numTri * 2 * sizeof(REAL));

  for(int i = 0; i < numTri ; i++)
  {
    int index = i * 2;

    HCVertex *vertex = points[i];
    in.pointlist[index] = vertex->x();
    in.pointlist[index+1] = vertex->y();
    in.pointattributelist[i] = vertex->z();
    in.pointmarkerlist[i] = 1;

    in.segmentlist[index] = i;

    if(i == numTri - 1)
    {
      in.segmentlist[index+1]  = 0;
    }
    else
    {
      in.segmentlist[index+1] = i + 1;
    }

    in.segmentmarkerlist[i] = 1;

  }

  triangulate("pz",&in,&out,&vrout);

  freeTriangleIOObject(&in);
  out.holelist = NULL;
  out.regionlist = NULL;

  for(int i = 0 ; i < out.numberoftriangles; i++)
  {
    int index = i*3;

    HCVertex *v1 = points[out.trianglelist[index]];
    HCVertex *v2 = points[out.trianglelist[index+1]];
    HCVertex *v3 = points[out.trianglelist[index+2]];

    m_mutex.lock();
    HCTriangle *triangle = m_outputTIN->createTriangle(v1,v2,v3);
    m_mutex.unlock();
  }

  freeTriangleIOObject(&out);
  freeTriangleIOObject(&vrout);
}

void TriangleComponent::triangulateOrbit(const QVector<HCVertex*> &vertices, int length, triangulateio *out)
{

  triangulateio in;
  triangulateio vrout;

  nullifyTriangleIOObjectPtrs(&in);
  nullifyTriangleIOObjectPtrs(out);
  nullifyTriangleIOObjectPtrs(&vrout);

  in.numberofpoints = length;
  in.pointlist = (REAL*) malloc(length * 2 * sizeof(REAL));
  in.pointattributelist = (REAL*) malloc(length * sizeof(REAL));
  in.pointmarkerlist = (int*) malloc(length * sizeof(REAL));
  in.numberofsegments = length;
  in.segmentmarkerlist = (int*) malloc(length * sizeof(REAL));
  in.segmentlist = (int*) malloc(length * 2 * sizeof(REAL));

  for(int i = 0; i < length ; i++)
  {
    int index = i * 2;
    HCVertex *vertex = vertices[i];

    in.pointlist[index] = vertex->x();
    in.pointlist[index+1] = vertex->y();
    in.pointattributelist[i] = vertex->z();
    in.pointmarkerlist[i] = 1;

    in.segmentlist[index] = index;

    if(i == length -1)
    {
      in.segmentlist[index]  = 0;
    }
    else
    {
      in.segmentlist[index] = index + 1;
    }

    in.segmentmarkerlist[index] = 1;

  }

  triangulate("pz",&in,out,&vrout);

  freeTriangleIOObject(&in);
  out->holelist = NULL;
  out->regionlist = NULL;
  freeTriangleIOObject(&vrout);
}

void TriangleComponent::setOrbitDestinationMarkers(HCVertex *vertex, int marker)
{
  HCEdge *edge = vertex->edgeInternal();
  HCEdge *endEdge = edge;

  do
  {
    edge->destInternal()->setMarker(marker);
    edge = edge->origNextInternal();

  }while(edge != endEdge);
}

double TriangleComponent::volumeUnderTriangle(const HCVertex *p1, HCVertex *p2, HCVertex *p3)
{
  return 0;
}

void TriangleComponent::triangleNormalVector(const HCVertex *p1, HCVertex *p2, HCVertex *p3, double nvect[])
{
  double u[3] = {p3->x() - p2->x(),p3->y() - p2->y(),p3->z() - p2->z()};
  double v[3] = {p1->x() - p2->x(),p1->y() - p2->y(),p1->z() - p2->z()};
  crossProduct(u,v,nvect);
}

void TriangleComponent::normalizeVector(double uvect[])
{
  double x= uvect[0];
  double y= uvect[1];
  double z= uvect[2];

  double denom = sqrt(x*x + y*y + z*z);

  uvect[0] = x /denom;
  uvect[1] = y /denom;
  uvect[2] = z /denom;
}

double TriangleComponent::dotProduct(double u[], double v[])
{
  return u[0]*v[0] + u[1]*v[1] + u[2]*v[2];
}

void TriangleComponent::crossProduct(double u[], double v[], double out[])
{
  out[0] = u[1]*v[2] - u[2]*v[1];
  out[1] = u[2]*v[0] - u[0]*v[2];
  out[2] = u[0]*v[1] - u[1]*v[0];
}

bool TriangleComponent::getBoolFromString(const QString &boolString)
{
  if(!boolString.trimmed().compare("True", Qt::CaseInsensitive) ||
     !boolString.trimmed().compare("Yes", Qt::CaseInsensitive) ||
     !boolString.trimmed().compare("1", Qt::CaseInsensitive)
     )
  {
    return true;
  }

  return false;
}
