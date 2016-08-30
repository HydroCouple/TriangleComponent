#ifndef TRIANGLECOMPONENTINFO_H
#define TRIANGLECOMPONENTINFO_H

#include "trianglecomponent_global.h"
#include "core/modelcomponentinfo.h"

class TRIANGLECOMPONENT_EXPORT TriangleComponentInfo : public ModelComponentInfo
{
    Q_OBJECT
    Q_PLUGIN_METADATA(IID "TriangleComponentInfo")

  public:
    TriangleComponentInfo(QObject *parent = nullptr);

    virtual ~TriangleComponentInfo();

    HydroCouple::IModelComponent* createComponentInstance() override;
};

Q_DECLARE_METATYPE(TriangleComponentInfo*)
#endif // TRIANGLECOMPONENTINFO_H
