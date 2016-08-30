/*!
 *  \file   trianglecomponent.h
 *  \author Caleb Amoa Buahin <caleb.buahin@gmail.com>
 *  \version   1.0.0.0
 *  \section   Description
 *  This header files contains implementation of the IModelComponent interface definition for the Triangulation Component
 *  \section License
 *  hydrocouplespatial.h, associated files and libraries are free software;
 *  you can redistribute it and/or modify it under the terms of the
 *  Lesser GNU General Public License as published by the Free Software Foundation;
 *  either version 3 of the License, or (at your option) any later version.
 *  hydrocouplespatial.h its associated files is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.(see <http://www.gnu.org/licenses/> for details)
 *  \date 2014-2016
 *  \pre
 *  \bug
 *  \todo Enable geometry input exchange items for defining geometric inputs into triangulation.
 *  \warning
 */
#ifndef TRIANGLECOMPONENT_H
#define TRIANGLECOMPONENT_H

#include "trianglecomponent_global.h"
#include "core/abstractmodelcomponent.h"
#include "core/triangle.h"

#include <QFileInfo>

class IdBasedArgumentQString;
class IdBasedArgumentDouble;
class Argument1DString;
class HCGeometryArgumentDouble;
class HCTINOutputDouble;
class HCTIN;
class HCVertex;
class SpatialReferenceSystem;

class TRIANGLECOMPONENT_EXPORT TriangleComponent : public AbstractModelComponent
{
    friend class TriangleComponentInfo;

    Q_OBJECT

  public:

    enum ThinningMode
    {
      None,
      NormalVectorCriteria,
      VolumeChangeCriteria
    };

    TriangleComponent(const QString &id, TriangleComponent* parent = nullptr);

    TriangleComponent(const QString &id, const QString &caption, TriangleComponent* parent = nullptr);

    virtual ~TriangleComponent();

    void initialize() override;

    HydroCouple::IModelComponent* clone() override;

    QList<QString> validate() override;

    void prepare() override;

    void update(const QList<HydroCouple::IOutput*> &requiredOutputs = QList<HydroCouple::IOutput*>()) override;

    void finish() override;

    IdBasedArgumentQString *identifiersArgument() const;

    HCGeometryArgumentDouble *inputVerticesArgument() const;

    HCGeometryArgumentDouble *inputEdgesArgument() const;

    HCGeometryArgumentDouble *inputBoundariesArgument() const;

    HCGeometryArgumentDouble *inputHolesArgument() const;

    Argument1DString *commandSwitchesArgument() const;

    IdBasedArgumentDouble *thinningOptionsArgument() const;

    IdBasedArgumentQString *outputOptionsArgument() const;

  private:

    void createArguments();

    void createIdentifierArgument();

    void createCommandsSwitchesArgument();

    void createGeometriesArguments();

    void createThinningArguments();

    void createOutputFileArgument();

    bool initializeIdentifierArgument(QString& message);

    bool initializeCommandSwitchesArgument(QString &message);

    bool initializeExchangeItems(QString &message);

    bool initializeThinningArguments(QString &message);

    bool initializeOutputsArguments(QString &message);

    void createTriangleIOObject();

    void nullifyTriangleIOObjectPtrs(triangulateio *triangleObject);

    void freeTriangleIOObject(triangulateio *triangleObject);

    HCTIN* readTINData(triangulateio *triangulation);

    int performThinning();

    void updateInputExchangeItems();

    void updateOutputExchangeItems();

    void writeOutput();

  private:

    void writeOutputVerticesToCSV();

    bool isEligibleForThinning(HCVertex *vertex, int &numtriangles);

    bool performThinningOnVertex(HCVertex *vertex, ThinningMode thinningMode , double criteriaValue);

    void removeAndTriangulateOrbit(HCVertex *vertex);

    void triangulateOrbit(const QVector<HCVertex*> &vertices, int length, triangulateio *output);

    void setOrbitDestinationMarkers(HCVertex *vertex, int marker = 1);

    static double volumeUnderTriangle(const HCVertex *p1, HCVertex *p2, HCVertex *p3);

    static void triangleNormalVector(const HCVertex *p1, HCVertex *p2, HCVertex *p3 , double nvect[]);

    static void normalizeVector(double uvect[]);

    static double dotProduct(double u[] , double v[]);

    static void crossProduct(double u[], double v[], double out[]);

    static bool getBoolFromString(const QString &boolString);

  private:
    HCGeometryArgumentDouble *m_points, *m_edges, *m_boundaries, *m_holes;
    HCTINOutputDouble *m_outputTriangulation;
    IdBasedArgumentQString *m_identifiersArgument, *m_outputsArgument;
    Argument1DString *m_commandSwitches;
    IdBasedArgumentDouble *m_thinningArgument;
    HCTIN *m_outputTIN;
    SpatialReferenceSystem *m_spatialReferenceSystem;
    bool m_initialized, m_prepared;
    int m_maxNumberThinningIterations, m_maxNumberOfPointsPerThinningIteration;
    double m_thinningCriteriaValue;
    ThinningMode m_thinningMode;
    QString m_switches;
    triangulateio m_triangle;
    QMutex m_mutex;
    bool m_writePointsToFile, m_writePointsToCSV, m_writePatchsToFile, m_writeTriangulationToNetCDF , m_excludeOuterEdgesFromThinning;
    QString m_gdalVectorDriverName;
    QFileInfo m_pointsFile, m_patchesFile, m_netCDFFile, m_CSVFile;
};

#endif // TRIANGLECOMPONENT_H
