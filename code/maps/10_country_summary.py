####################################################

import os, sys
from qgis.analysis import QgsRasterCalculator, QgsRasterCalculatorEntry,  QgsZonalStatistics

# Add polygon layer 

#specify polygon shapefile vector
polygonLayer = QgsVectorLayer('/Users/amirsiraj/GoFlex/DataNow/Admin_units/south_America.shp', 'zonepolygons', "ogr") 


vtypes = ["cumI_SA_repl_", "cumPI_SA_repl_", "cumBI_SA_repl_"]
sufs = ["I","PI","BI"]

rplst = range(1,8,1) + range(9,101,1) 
rplst = range(1,8,1)+ range(9,59,1) + range(60,92,1) + range(93,101,1)
rplst = range (501,1001,1)
rng2 = range (1,2,1) + range(3,4,1)

for vtype in rng2 :
	for repl in rplst :  
		# Add raster 1
		rasterFilePath = '/Users/amirsiraj/GoFlex/zika/outputbk/'+ vtypes[vtype-1] + str(repl) +'_stat.bil'
		zoneStat = QgsZonalStatistics (polygonLayer, rasterFilePath, sufs[vtype-1]+str(repl), 1, QgsZonalStatistics.Sum)
		zoneStat.calculateStatistics(None)
		print(vtype)
		print(repl)
		