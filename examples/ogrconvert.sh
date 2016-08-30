#/bin/bash
ogr2ogr -f "CSV" edge.csv edge.shp -lco GEOMETRY=AS_WKT
ogr2ogr -f "CSV" hole.csv hole.shp -lco GEOMETRY=AS_WKT
ogr2ogr -f "CSV" sample_points.csv sample_points.shp -lco GEOMETRY=AS_WKT