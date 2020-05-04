#!/usr/bin/env python
# ---------------------------------------------------------------------------
# Copyright 2016, Houston Engineering Inc., All Rights Reserved
# Description: This script creates a raster whose cells measure the length of time
#              in seconds, that it takes water to flow across it and then accumulates
#              the time from the cell to the outlet of the watershed. water velocity
#              is calculated as a function of hydraulic radius, Manning's N and slope
#              Using a DEM as a source, flow direction, flow accumulation and slope (percent)
#              is calculated. Using flow accumulation, each cell is assigned one of three
#              flow regimes. Slopes are processed so that there are no cells with zero slopes.
#
#              Landcover is used to assign flow N and R values to non channelized cells. Users can
#              use various landcover sources but these must be matched to the N and R table that
#              is stored in the /data folder in this distribution. The N and R table must
#              contain one record for each of the raster valeus in the landcover grid and there
#              must be the same fields in that table for the program to work. R and N values
#              are assigned based on landcover and flow regime.
#
#              Hydraulic Radius is calucated using:
#                  (X)a + Y where:
#                       X = hydraulic radius factor 1 (User setting, 0.0032 default)
#                       a = drainage area in square miles (calculated from Flow Accumulation)
#                       Y = hydraulic radius factor 2 (User setting, 1.7255 default)
#
#              Velocity is calculated as:
#                 vel(feet/second) = (1.49 * Hydraulic Radius^0.667 * slope^0.5) / Mannings N
#
#              Travel Time across individual cells is calculated as:
#                  1 / vel * 0.3048 (if the xy units are meters, otherwise the conversion factor is not used
# Author: Paul Hedlund via Tim Loesch
# Created: April 2016
# ---------------------------------------------------------------------------
import arcpy,time,traceback,datetime, os, math, platform
from arcpy.sa import *

def main():
	try:
		#Start tracking time -------------------------------------------------------------------------
		start = time.time()

		#Toolbar settings --------------------------------------------
		elevGrid = arcpy.GetParameterAsText(0)
		fdr_surf = arcpy.GetParameterAsText(1)
		fac_surf = arcpy.GetParameterAsText(2)
		zUnits = arcpy.GetParameterAsText(3)
		minSlope = arcpy.GetParameter(4)
		MediumRetardanceAcres = arcpy.GetParameter(5)
		LowRetardanceAcres = arcpy.GetParameter(6)
		hydRad1 = arcpy.GetParameter(7)
		hydRad2 = arcpy.GetParameter(8)
		luGrid = arcpy.GetParameterAsText(9)
		Mask = arcpy.GetParameterAsText(10)
		blnChannelCreation = arcpy.GetParameter(11)
		Connectors = arcpy.GetParameterAsText(12)
		channelFC = arcpy.GetParameterAsText(13)
		nField = arcpy.GetParameterAsText(14)
		CutInterval = arcpy.GetParameter(15)
		OutputWorkspace = arcpy.GetParameterAsText(16)

		#In script settings --------------------------------------------
		if not elevGrid: #If no parameter
			#ProcessingPath = r'D:\ZMapdata\PTMA\TravelTime\Workspace.gdb'
			ProcessingPath = r'D:\ZMapdata\PTMA\Watonwan\Travel_Time\Dbogart.gdb'
			#ProcessingPath = r'D:\ZMapdata\PTMA\MC_10m_TT_inputs.gdb'
			elevGrid = ProcessingPath + os.sep + 'hyd_dem_c'
			fdr_surf = ProcessingPath + os.sep + 'fdr_total_c'
			fac_surf = ProcessingPath + os.sep + 'fac_total_c'
			zUnits = 'Meters'
			minSlope = 0.019
			MediumRetardanceAcres = 2
			LowRetardanceAcres = 40
			hydRad1 = 0.0032
			hydRad2 = 1.7255
			luGrid = ProcessingPath + os.sep + 'nlcd'
			#luGrid = ProcessingPath + os.sep + 'landusethislongthislongthislongthislongthislong'
			blnChannelCreation = False
			Mask = ProcessingPath + os.sep + 'mask'
			Connectors = ProcessingPath + os.sep + 'connectors'
			channelFC = ProcessingPath + os.sep + 'StreamsForN'
			nField = 'N_Value'
			CutInterval = 6
			OutputWorkspace = r'D:\ZMapdata\PTMA\Watonwan\Travel_Time\Dbogart2.gdb'

		#Check out Spatial Analyst Extensions --------------------------------------------------------
		blnLicSA = ExtensionManagment('spatial','Spatial Analyst')
		if blnLicSA:
			return

		# Set Geoprocessing environments
		elevDesc = arcpy.Describe(elevGrid)
		cellSize = elevDesc.meancellwidth
		if str(cellSize).endswith('.0'):
			cellSize = int(round(elevDesc.meancellwidth))
		arcpy.env.cellSize = cellSize
		#scratchDir = OutputWorkspace  #TEMP
		scratchDir = arcpy.env.scratchGDB
		arcpy.env.scratchworkspace = scratchDir
		arcpy.env.overwriteOutput = True
		arcpy.env.parallelProcessingFactor = "100%"
		arcpy.env.extent = elevDesc.extent
		arcpy.env.snapRaster = elevGrid
		#arcpy.env.mask = elevGrid

		tGrid = os.path.join(OutputWorkspace,'tt_grid')
		ds_tt = os.path.join(OutputWorkspace,'ds_tt')
		us_tt = os.path.join(OutputWorkspace,'us_tt')
		faRank = os.path.join(scratchDir,'farank')
		slope2 = os.path.join(scratchDir,'slope2')
		#luGrid = os.path.join(scratchDir,'landuse')
		#-------------------------------------------------

		#Clear scratch data
		ClearScratch(scratchDir)

		#Validate and Get Raster Cell size ----------------------------------------------------------
		arcpy.AddMessage("Validate grid cell sizes ...")
		InputRasterList = []
		InputRasterList.append(elevGrid)
		InputRasterList.append(fdr_surf)
		InputRasterList.append(fac_surf)
		InputRasterList.append(luGrid)
		blnValidate = ValidateGridSizeMatch(InputRasterList)
		if blnValidate:
			return

		arcpy.AddMessage("Validate projection ...")
		blnNoMeter = False
		facDesc = arcpy.Describe(fac_surf)
		fdrDesc = arcpy.Describe(fdr_surf)
		if facDesc.spatialReference.linearUnitName != 'Meter':
			blnNoMeter = True
		if fdrDesc.spatialReference.linearUnitName != 'Meter':
			blnNoMeter = True
		if elevDesc.spatialReference.linearUnitName != 'Meter':
			blnNoMeter = True

		if blnNoMeter:
			arcpy.AddError('All input rasters must be in a meters projection!')
			return

		#Get landuse field
		#arcpy.CopyRaster_management(landuse, luGrid)
		luField = os.path.split(luGrid)[1]

		#----------------------------------------------------------------------- Channel
		if blnChannelCreation:
			RECLSS_Poly = scratchDir + os.sep + 'RECLSS_Poly'
			Union_shp = scratchDir + os.sep + 'Union_shp'
			Wetland_Connectors_TOC_Calcs = scratchDir + os.sep + 'Wetland_Connectors_TOC_Calcs'
			Stream_2km_3m = scratchDir + os.sep + 'Stream_2km_3m'

			arcpy.AddMessage('Clip wetland connectors ...')
			valElapsed = time.time()
			arcpy.Clip_analysis(Connectors,Mask,Wetland_Connectors_TOC_Calcs,'')
			PrintTime(start,valElapsed,True)

			#--------------
			arcpy.AddMessage('Reclassify flow accumulation raster ...')
			valElapsed = time.time()

			#DA_sqmi_ras = (Raster(fac_surf) * cellSize ** 2) / 2589987.945
			DA_sqmi_ras = (Raster(fac_surf) * Power(float(cellSize),2)) / 2589987.945
			DA_sqmi_ras.save('DA_sqmi_ras')
			RECLSS_FAC = Reclassify(DA_sqmi_ras,'VALUE',"0 5 1;5 20 2;20 10000000000000 3",'DATA')
			#RECLSS_FAC = Reclassify(DA_sqmi_ras,'VALUE',RemapRange([[0,5,1],[5,20,2],[20,10000000000000,3]]),'DATA')
			RECLSS_FAC.save('RECLSS_FAC')
			PrintTime(start,valElapsed,True)
			#---------------

			arcpy.AddMessage('Selection and merge process ...')
			valElapsed = time.time()
			arcpy.RasterToPolygon_conversion(RECLSS_FAC, RECLSS_Poly,'NO_SIMPLIFY','VALUE')

			#Faster way
			arcpy.AddField_management(RECLSS_Poly,'DrainageN','DOUBLE','10','5')
			with arcpy.da.UpdateCursor(RECLSS_Poly,['GRIDCODE','DrainageN']) as cursorDrain:
				for rowDrain in cursorDrain:
					if rowDrain[0] == 1:
						rowDrain[1] = 0.045
					if rowDrain[0] == 2:
						rowDrain[1] = 0.040
					if rowDrain[0] == 3:
						rowDrain[1] = 0.035
					cursorDrain.updateRow(rowDrain)
			del cursorDrain
			PrintTime(start,valElapsed,True)

			arcpy.AddMessage('Union ...')
			valElapsed = time.time()
			arcpy.Union_analysis([Wetland_Connectors_TOC_Calcs,RECLSS_Poly], Union_shp,'ALL','','GAPS')
			#arcpy.Union_analysis([Wetland_Connectors_TOC_Calcs,Merge_shp], Union_shp,'ALL','','GAPS')
			PrintTime(start,valElapsed,True)

			arcpy.AddMessage('Second merge process ...')
			valElapsed = time.time()

			#Faster way
			arcpy.AddField_management(Union_shp,'N_Value','FLOAT','10','5')
			with arcpy.da.UpdateCursor(Union_shp,['DrainageN','N_Value','ManningsN','Connector']) as cursorDrain2:
				for rowDrain2 in cursorDrain2:
					if rowDrain2[3] == 'Lake':
						rowDrain2[1] = rowDrain2[2]
					elif rowDrain2[3] == 'Wetland':
						rowDrain2[1] = rowDrain2[2]
					else:
						rowDrain2[1] = rowDrain2[0]
					cursorDrain2.updateRow(rowDrain2)
			del cursorDrain2
			PrintTime(start,valElapsed,True)

			#-------
			arcpy.AddMessage('Stream to feature ...')
			valElapsed = time.time()

			#-----------
			#FL_km_ras = (Raster(fac_surf) * cellSize ** 2) / 1000000
			FL_km_ras = (fac_surf * Power(float(cellSize),2)) / 1000000
			Stre_null_2km_3 = SetNull(FL_km_ras,1,'VALUE < 2')
			#--------
			try:
				StreamToFeature(Stre_null_2km_3,fdr_surf, Stream_2km_3m,'NO_SIMPLIFY')
			except:
				if str(platform.architecture()[0]) == '64bit':
					arcpy.AddError('Stream to Feature failed!  You are in background processing mode.  Try removing your data from the map document to avoid a lock. Then try again.')
					return
				else:
					arcpy.AddError('Stream to Feature failed!  Try running in the foreground.')
					return

			PrintTime(start,valElapsed,True)

			arcpy.AddMessage('Intersect analysis ...')
			valElapsed = time.time()
			channelFC = OutputWorkspace + os.sep + 'StreamsForN'
			arcpy.Intersect_analysis([Union_shp,Stream_2km_3m],channelFC,'ALL','','INPUT')
			#arcpy.Intersect_analysis([Merge2_shp,Stream_2km_3m],channelFC,'ALL','','INPUT')
			nField = 'N_Value'
			PrintTime(start,valElapsed,True)

			arcpy.AddMessage('Delete temp stuff ...')
			valElapsed = time.time()
			arcpy.Delete_management(Wetland_Connectors_TOC_Calcs)
			arcpy.Delete_management(RECLSS_FAC)
			arcpy.Delete_management(Stream_2km_3m)
			arcpy.Delete_management(RECLSS_Poly)
			arcpy.Delete_management(Union_shp)
			arcpy.Delete_management(DA_sqmi_ras)
			PrintTime(start,valElapsed,True)
		else:
			desc = arcpy.Describe(Mask)
			OIDField = desc.OIDFieldName
			Stre_null_2km_3 = OutputWorkspace + os.sep + 'Strenullkm'
			arcpy.FeatureToRaster_conversion(Mask,OIDField,Stre_null_2km_3,cellSize)

		#-----------------------------------------------------------------------

		arcpy.AddMessage('Prepare hydraulic values ...')
		valElapsed = time.time()
		cellArea = elevDesc.meancellwidth * elevDesc.meancellheight
		cellAcres = cellArea * 0.0002471
		cellSqMiles = cellArea * 0.0000003861003
		xyUnits = elevDesc.SpatialReference.LinearUnitname
		xyUnits = xyUnits.upper()
		zUnits = zUnits.upper()

		if xyUnits == 'METERS' and zUnits == 'FEET':
			zFactor = 0.3048
		elif xyUnits == 'FEET' and zUnits == 'METERS':
			zFactor = 3.280839895
		else:
			zFactor = 1

		# this table stores the various values of mannings M and Channel hydraulic radius based on various landcovers and retardance
		NandMTable = {
	        11:[0.002,1,0.002,2,0.002,3],
	        12:[0.01,1,0.01,2,0.01,3],
	        21:[0.05,0.1,0.035,0.6,0.035,1],
	        22:[0.05,0.1,0.035,0.6,0.035,1],
	        23:[0.03,0.1,0.011,0.6,0.011,1],
	        24:[0.025,0.1,0.011,0.6,0.011,1],
	        31:[0.05,0.1,0.035,0.6,0.035,1],
	        41:[0.1,0.1,0.08,0.6,0.08,1],
	        42:[0.1,0.1,0.08,0.6,0.08,1],
	        43:[0.1,0.1,0.08,0.6,0.08,1],
	        52:[0.2,0.1,0.125,0.6,0.125,1],
	        71:[0.2,0.1,0.1,0.4,0.08,1],
	        81:[0.07,0.1,0.08,0.4,0.08,1],
	        82:[0.07,0.1,0.08,0.4,0.08,1],
	        90:[0.2,0.1,0.125,0.4,0.125,0.6],
	        95:[0.01,0.1,0.01,0.4,0.01,0.6]
	    }

		#Input parameter defines the minimum percent slope.
		MediumRetardanceCells = int(MediumRetardanceAcres) / cellAcres
		LowRetardanceCells = int(LowRetardanceAcres) / cellAcres
		PrintTime(start,valElapsed,True)

		arcpy.AddMessage('Processing elevation data ...')
		valElapsed = time.time()
		processElevation(scratchDir,elevGrid,fdr_surf,fac_surf,faRank,slope2,MediumRetardanceAcres,MediumRetardanceCells,LowRetardanceAcres,LowRetardanceCells,channelFC,zFactor,cellSize,minSlope,CutInterval,Stre_null_2km_3,OutputWorkspace)
		PrintTime(start,valElapsed,True)

		arcpy.AddMessage('Processing channel data ...')
		valElapsed = time.time()
		rGrid,nGrid = processChannelData(scratchDir,NandMTable,channelFC,nField,cellSize,cellSqMiles,hydRad1,hydRad2,fac_surf,luGrid,faRank,elevGrid,luField,OutputWorkspace)
		PrintTime(start,valElapsed,True)

		arcpy.AddMessage('Calculating velocity ...')
		valElapsed = time.time()
		calcVelocity(scratchDir,tGrid,ds_tt,us_tt,slope2,rGrid,nGrid,xyUnits,cellSize)
		PrintTime(start,valElapsed,True)

		arcpy.AddMessage('Calculating travel DOWNSTREAM (this may take awhile) ...')
		valElapsed = time.time()
		travTimeGridDDSV = FlowLength(fdr_surf,'DOWNSTREAM',tGrid)
		ttDDhours = Divide(travTimeGridDDSV,3600)
		ttDDhours.save(ds_tt)
		PrintTime(start,valElapsed,True)

		arcpy.AddMessage('Calculating travel UPSTREAM (this may take awhile) ...')
		valElapsed = time.time()
		travTimeGridSV = FlowLength(fdr_surf,'UPSTREAM',tGrid)
		tthours = Divide(travTimeGridSV,3600)
		tthours.save(us_tt)
		PrintTime(start,valElapsed,True)

		arcpy.AddMessage('Delete more temp stuff ...')
		valElapsed = time.time()
		arcpy.Delete_management(faRank) #TODO
		arcpy.Delete_management(slope2)
		#arcpy.Delete_management(rGrid)  #Told to keep
		arcpy.Delete_management(luGrid)
		arcpy.Delete_management(nGrid)
		arcpy.Delete_management(Stre_null_2km_3)
		arcpy.Delete_management(travTimeGridSV)
		PrintTime(start,valElapsed,True)

		arcpy.AddMessage('Complete!')
	except:
		strError = sys.exc_info()[1]
		stacktrace = sys.exc_info()[2]
		try:
			strScriptInside = str(traceback.extract_tb(stacktrace)[1][0]) #line-number inside function
			strLineInside = str(traceback.extract_tb(stacktrace)[1][1]) #line-number inside function
		except IndexError:
			strScriptInside = 'N/A'
			strLineInside = 'N/A'
		strLine = str(stacktrace.tb_lineno) #line-number in this function
		timestamp = str(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
		#Note don't want to openly display script path
		arcpy.AddError(timestamp + ': Line:' + strLine + ' - Inner Error Line:' + strLineInside + ' - Error:' + str(strError))

#//////////////// Functions//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
def TimeParser(t):
	rediv = lambda ll,b : list(divmod(ll[0],b)) + ll[1:]
	return '%d:%02d:%02d.%03d' % tuple(reduce(rediv,[[t*1000,], 1000,60,60]))

def GetCellSize(InputRaster):
	#Get grid size
	GetCellSize = arcpy.GetRasterProperties_management (InputRaster,'CELLSIZEX')
	GridSize = GetCellSize.getOutput(0)
	arcpy.AddMessage('Raster Cellsize = ' + GridSize)
	return float(GridSize)

def VeryDataFilesExist(start,BasePath,ProcessingPath):
	arcpy.AddMessage('Validate required data exists ...')
	if not arcpy.Exists(ProcessingPath + os.sep + 'catchment'):
		return ProcessingPath + os.sep + 'catchment'
	return None

def CheckSchemaLock(elevGrid):
	arcpy.AddMessage('Schema lock test ...')
	if arcpy.Exists(elevGrid):
		if not arcpy.TestSchemaLock(elevGrid):
			return elevGrid
	return None

def PrintErrors(strError,stacktrace,strTool):
	try:
		strLineInside = str(traceback.extract_tb(stacktrace)[1][1]) #line-number inside function
	except IndexError:
		strLineInside = 'N/A'
	strLine = str(stacktrace.tb_lineno) #line-number in this function
	timestamp = str(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
	#Note don't want to openly display script path
	arcpy.AddError(timestamp + ': Error ' + strTool + ' script! - ' + 'Line: ' + strLine + ' - Inner Error Line: ' + strLineInside + ' - Error:' + str(strError))

def ExtensionManagment(ExtID, ExtName):
	if arcpy.CheckExtension(ExtID) == 'Available':
		arcpy.AddMessage('Checking out the ' + ExtName + ' extension ...')
		arcpy.CheckOutExtension(ExtID)
		return False
	elif arcpy.CheckExtension(ExtID) == 'Unavailable':
		arcpy.AddError('The ' + ExtName + ' extension is currently unavailable')
		return True
	elif arcpy.CheckExtension(ExtID) == 'NotLicensed':
		arcpy.AddError('You are currently not license to use ' + ExtName + ', the tool can not be run!')
		return True
	elif arcpy.CheckExtension(ExtID) == 'Failed':
		arcpy.AddError('A system failure occurred during the request!')
		return True

def FieldExist(featureclass, fieldname):
	fieldList = arcpy.ListFields(featureclass, fieldname)

	fieldCount = len(fieldList)

	if (fieldCount == 1):
		return True
	else:
		return False

def PrintTime(start, valElapsed, blnShowElapsed):
	#Get total and elapsed time
	if not blnShowElapsed:
		arcpy.AddMessage(' - Total - ' + TimeParser(time.time() - start)[:-4])
	else:
		arcpy.AddMessage(' - Total - ' + TimeParser(time.time() - start)[:-4] + ' (Elapsed: ' + TimeParser(time.time() - valElapsed)[:-4] + ')')

def ClearScratch(LocationPath):
	#Clear Scratch
	arcpy.env.workspace = LocationPath
	try:
		featureclasses = arcpy.ListFeatureClasses()
		for fc in featureclasses:
			arcpy.Delete_management(fc)

		rasterclasses = arcpy.ListRasters()
		for fcR in rasterclasses:
			arcpy.Delete_management(fcR)

		tables = arcpy.ListTables()
		for table in tables:
			arcpy.Delete_management(table)
	except:
		arcpy.AddWarning('Old temp layer could not be deleted because it is locked - ' + str(fc))

def processElevation(scratchDir,elevDEM,fdr_surf,fac_surf,faRank,slope2,MediumRetardanceAcres,MediumRetardanceCells,LowRetardanceAcres,LowRetardanceCells,channelFC,zFactor,cellSize,minSlope,CutInterval,Stre_null_2km_3,OutputWorkspace):
	slope = os.path.join(scratchDir,'slope')
	channelSlopepre = os.path.join(scratchDir,'chanslopepre')
	#slopelandscape = os.path.join(scratchDir,'slopelandscape')
	channeldissolve = os.path.join(scratchDir,'channeldissolve')
	outexplode = os.path.join(scratchDir,'outexplode')
	outexplode2 = os.path.join(scratchDir,'outexplode2')
	outbuffer = os.path.join(scratchDir,'outbuffer')
	tempPT = os.path.join(scratchDir,'tempPT')
	tempLINE = os.path.join(scratchDir,'tempLINE')
	slopemask = os.path.join(scratchDir,'slopemask')
	elevDEMmask = os.path.join(scratchDir,'elevDEMmask')

	# Process: Reclassify the flow accumulation grid into three zones...
	arcpy.AddMessage(' - Creating flow accumulation zone grid ...')
	arcpy.AddMessage(' - Flow regime 1 = drainage area < ' + str(MediumRetardanceAcres) + ' acre(s), flow regime 2 = drainage area ' + str(MediumRetardanceAcres) + '-' +  str(LowRetardanceAcres) + ' acres, flow regime 3 = drainage area > ' + str(LowRetardanceAcres) + ' acres ...')

	reclassStr = '0 ' + str(MediumRetardanceCells) + ' 1;' + str(MediumRetardanceCells) + ' ' + str(LowRetardanceCells) + ' 2;' + str(LowRetardanceCells) + ' 1000000000000000 3'
	faRankSV = Reclassify(fac_surf,'VALUE',reclassStr,'DATA')
	faRankSV.save(faRank)

	# Process: Slope
	arcpy.AddMessage(' - Creating landscape slope from elevation DEM ...')
	slopecmd = Slope(elevDEM,'PERCENT_RISE',zFactor,'PLANAR')
	slopecmd.save(slope)

	#---New Slope workaround -----------------------------------------------------------------------------------------------------------------
	#Who needs advanced! - we don't ;p
	elevDEMmask2 = ExtractByMask(elevDEM, Stre_null_2km_3)
	elevDEMmask2.save(elevDEMmask)

	arcpy.AddMessage(' - Dissolve channel ...')
	arcpy.Dissolve_management(channelFC,channeldissolve,['GRIDCODE'],'','SINGLE_PART','DISSOLVE_LINES')  #TODO - always GRIDCODE?

	desc = arcpy.Describe(channeldissolve)
	InSR = desc.spatialReference
	arcpy.CreateFeatureclass_management(scratchDir,'tempPT','POINT','','','',InSR)
	arcpy.CreateFeatureclass_management(scratchDir,'tempLINE','POLYLINE', None, None, None, InSR)

	arcpy.AddField_management(tempPT,'MatchIDFrom','LONG')
	arcpy.AddField_management(tempPT,'MatchIDTo','LONG')
	arcpy.AddField_management(tempLINE,'MatchID','LONG')
	arcpy.AddField_management(tempPT,'PointID','LONG')

	#Create point intervals and line intervals
	arcpy.AddMessage(' - Split lines by maximum length interval ...')
	iPTcnt = 0
	iMatchID = 0
	DictPT = {}
	cursorInsertLine = arcpy.da.InsertCursor(tempLINE, ['SHAPE@','MatchID'])
	countChannel = int(arcpy.GetCount_management(channeldissolve).getOutput(0))
	with arcpy.da.SearchCursor(channeldissolve, ['SHAPE@','OID@']) as cursor:
		for row in cursor:
			blnMergeEnd = False
			blnSkipTic = False
			blnSkipMerge = False
			blnLast = False
			iCntLine = 0
			PolyLength = row[0].getLength("PLANAR")
			iCntLine = int(math.floor(PolyLength / CutInterval))
			iRemainder = PolyLength - (iCntLine * CutInterval)
			if cellSize > iRemainder:
				blnMergeEnd = True

			iPTcnt += 1
			startpt = row[0].firstPoint

			iMatchID += 1
			DictPT[iPTcnt] = [startpt.X, startpt.Y,-9999,iMatchID]

			iPT = 0
			if iCntLine > 0:
				for i in range(iCntLine + 1):
					iPTcnt += 1
					iPT += 1
					if iPT == (iCntLine + 1):
						if not blnMergeEnd:
							CutLineMeasure = row[0].segmentAlongLine(CutInterval * i,PolyLength)
						else:
							blnSkipTic = True
						blnLast = True
					elif iPT == iCntLine:
						if blnMergeEnd:
							CutLineMeasure = row[0].segmentAlongLine(CutInterval * i,((CutInterval * (i + 1)) + PolyLength))
							blnSkipMerge = True
						else:
							CutPT = row[0].positionAlongLine(CutInterval * (i + 1)).firstPoint
							CutLineMeasure = row[0].segmentAlongLine(CutInterval * i,CutInterval * (i + 1))

							DictPT[iPTcnt] = [CutPT.X, CutPT.Y,iMatchID,iMatchID + 1]
					else:
						CutPT = row[0].positionAlongLine(CutInterval * (i + 1)).firstPoint
						CutLineMeasure = row[0].segmentAlongLine(CutInterval * i,CutInterval * (i + 1))

						DictPT[iPTcnt] = [CutPT.X, CutPT.Y,iMatchID,iMatchID + 1]

					if not blnSkipTic:
						cursorInsertLine.insertRow((CutLineMeasure,iMatchID))

					if not blnSkipMerge:
						if not blnSkipTic:
							if not blnLast:
								iMatchID += 1
			else:
				cursorInsertLine.insertRow((row[0],iMatchID))

			iPTcnt += 1
			endpt = row[0].lastPoint
			DictPT[iPTcnt] = [endpt.X, endpt.Y,iMatchID,-9999]

	del cursorInsertLine

	#Create points
	cursorInsert = arcpy.da.InsertCursor(tempPT, ['SHAPE@XY','MatchIDFrom','MatchIDTo','PointID'])
	for ikey, iValue in DictPT.iteritems():
		cursorInsert.insertRow([(iValue[0], iValue[1]),iValue[2],iValue[3],ikey])
	del cursorInsert

	arcpy.AddMessage(' - Extract Point ...')
	ExtractValuesToPoints(tempPT,elevDEMmask,outexplode)
	arcpy.AddField_management(outexplode,'RASVAL1','DOUBLE')
	arcpy.CalculateField_management(outexplode,'RASVAL1','!RASTERVALU!','PYTHON_9.3')
	if FieldExist(outexplode,'RASTERVALU'):
		arcpy.DeleteField_management(outexplode,'RASTERVALU')

	arcpy.AddMessage(' - Buffer nulls and attach ...')
	arcpy.MakeFeatureLayer_management(outexplode,'explode2')
	arcpy.SelectLayerByAttribute_management('explode2','NEW_SELECTION', 'RASVAL1 IS NULL')
	counnull = int(arcpy.GetCount_management('explode2').getOutput(0))
	if counnull > 0:
		#Buffer to make sure we get a raster
		arcpy.Buffer_analysis('explode2',outbuffer,str(cellSize * 1.5) + " Meters")

		#Get minimum
		Zstat = ZonalStatisticsAsTable(outbuffer,'PointID',elevDEM,outexplode2,'DATA','MINIMUM')

		#---
		search_feats = {f[0]:f[1:] for f in arcpy.da.SearchCursor(outexplode2,['PointID','MIN'])}
		with arcpy.da.UpdateCursor('explode2',['PointID','RASVAL1']) as rowsCopyResults:
			for row in rowsCopyResults:
				if int(row[0]) in search_feats:
					row[1] = search_feats[int(row[0])][0]
					rowsCopyResults.updateRow(row)

		del rowsCopyResults
		search_feats.clear()
		del Zstat

		#---
		arcpy.SelectLayerByAttribute_management('explode2','NEW_SELECTION', 'RASVAL1 IS NULL')
		counnull3 = int(arcpy.GetCount_management('explode2').getOutput(0))
		if counnull3 > 0:
			with arcpy.da.UpdateCursor('explode2',['PointID','RASVAL1','SHAPE@']) as rowsCopyResults:
				for row in rowsCopyResults:
					with arcpy.da.SearchCursor(outexplode,['PointID','RASVAL1','SHAPE@']) as rowsSearch:
						for rowS in rowsSearch:
							if row[2].equals(rowS[2]):
								if row[1] == None:
									row[1] = rowS[1]
								rowsCopyResults.updateRow(row)

						del rowsSearch
			del rowsCopyResults
		#----

	DictPT.clear()

	arcpy.AddMessage(' - Calculate rise and slope ...')
	ElevTmpFrom = 0
	ElevTmpTo = 0
	ElevTmpPrev = 0
	sql_clause = (None,'ORDER BY PointID')
	with arcpy.da.SearchCursor(outexplode, ['MatchIDFrom','MatchIDTo','RASVAL1','PointID'],'',None,False,sql_clause) as cursorPT:
		for rowPT in cursorPT:
			blnAdd = True

			if rowPT[0] == -9999:
				ElevTmpFrom = rowPT[2]
				blnAdd = False
			else:
				ElevTmpFrom = ElevTmpPrev
				ElevTmpTo = rowPT[2]

			if blnAdd:
				if ElevTmpFrom == None:
					ElevTmpFrom = 0
				if ElevTmpTo == None:
					ElevTmpTo = 0
				rise = abs(ElevTmpFrom - ElevTmpTo)
				DictPT[rowPT[0]] = [ElevTmpFrom, ElevTmpTo,rise]

			ElevTmpPrev = rowPT[2]
	del cursorPT

	arcpy.AddField_management(tempLINE,'rise','DOUBLE')
	arcpy.AddField_management(tempLINE,'slope','DOUBLE')
	countLine = int(arcpy.GetCount_management(tempLINE).getOutput(0))
	with arcpy.da.UpdateCursor(tempLINE, ['MatchID','rise','slope','SHAPE@']) as cursor:
		for row in cursor:
			PolyLength = row[3].getLength("PLANAR")
			rise = DictPT[row[0]][2]
			row[1] = rise
			if rise == 0.0:
				row[2] = 0
			else:
				row[2] = (rise / PolyLength) * 100
			cursor.updateRow(row)
	del cursor
	DictPT.clear()

	#---------------
	# process slope along the channel, using the channel feature class as a mask...
	arcpy.AddMessage(' - Creating slope of stream channels ...')
	arcpy.PolylineToRaster_conversion(tempLINE,'slope',channelSlopepre,'MAXIMUM_LENGTH','Shape_Length',cellSize)

	themask = ExtractByMask(channelSlopepre, Stre_null_2km_3)
	themask.save(slopemask)

	#-----------------
	# now merge the channel slopes with the rest of the slopes
	arcpy.Mosaic_management(slopemask,slope,'LAST','FIRST','#','#','NONE','0','NONE')

	# make sure that all the slope values are greater than zero.
	arcpy.AddMessage(' - Setting zero slope to ' + str(minSlope) + ' in slope grid ...')
	slope2tell = Con(slope,minSlope,slope,'value <= ' + str(minSlope))
	slope2tell.save(slope2)

	#Copy outputs
	arcpy.CopyRaster_management(slopemask,OutputWorkspace + os.sep + 'channel_slope')
	arcpy.Copy_management(tempLINE,OutputWorkspace + os.sep + 'slope_line')

	arcpy.AddMessage(' - Delete temp stuff ...')
	arcpy.Delete_management(channelSlopepre) #TODO
	arcpy.Delete_management(slopemask)
	arcpy.Delete_management(channeldissolve)
	arcpy.Delete_management(outexplode)
	arcpy.Delete_management(outexplode2)
	arcpy.Delete_management(outbuffer)
	arcpy.Delete_management(tempPT)
	arcpy.Delete_management(tempLINE)
	arcpy.Delete_management(elevDEMmask)

def processChannelData(scratchDir,NandMDict,channelFC,nField,cellSize,cellSqMiles,hydRad1,hydRad2,fac_surf,luGrid,faRank,elevGrid,luField,OutputWorkspace):
	# convert the input stream polyline file to a raster for processing. Use the
	# this next step creates the stream trace (channel) hydraulic Radius Grid.....
	# this equation assumes that the drainage Area is in Square Miles...
	# from the equation HydRadius(ft) = (0.0032x + 1.7255) where X = drainage Area in Square Miles.....
	chNGrid = os.path.join(scratchDir,'channeln')
	nGrid = os.path.join(scratchDir,'nGrid')
	rGrid = os.path.join(OutputWorkspace,'rGrid')
	combGrid = os.path.join(scratchDir,'combGrid')
	chHydRadGrid2 = os.path.join(scratchDir,'chHydRadGrid2')

	arcpy.AddMessage(' - Rasterizing Channel data using field ' + nField + ' to create channelN ...')
	arcpy.PolylineToRaster_conversion(channelFC,nField,chNGrid,'MAXIMUM_LENGTH','NONE',cellSize)

	cellSqMilesStr = '%f' % cellSqMiles
	arcpy.env.mask = chNGrid

	arcpy.AddMessage(' - Creating Channel hydraulic radius grid ...')
	chHydRadGrid = (hydRad1 * (Raster(fac_surf) * float(cellSqMilesStr))) + hydRad2
	chHydRadGrid.save(chHydRadGrid2)

	arcpy.env.mask = elevGrid
	# merge the land use and flow accumulation zone grid...
	arcpy.AddMessage(' - Merging Flow Accumulation R and Land Use Grids ...')
	combGridSV = Combine([luGrid,faRank])
	combGridSV.save(combGrid)

	arcpy.AddMessage(' - Adding necessary fields, MN and HYDR, to combined grid ...')
	arcpy.AddField_management(combGrid,'MN','double')
	arcpy.AddField_management(combGrid,'hydr','double')
	arcpy.AddField_management(combGrid,'landuse','LONG')

	fields = arcpy.ListFields(combGrid)
	for field in fields:
		if not field.required:
			if field.name in luGrid:
				if field.name != 'landuse':
					arcpy.CalculateField_management(combGrid,'landuse','!' + field.name + '!','PYTHON_9.3')

	with arcpy.da.UpdateCursor(combGrid, ['landuse','farank','MN','hydr']) as cursor:
		for row in cursor:
			luVal = row[0]
			valList = NandMDict.get(luVal)
			flowVal = row[1]

			if flowVal == 1:
				n = valList[0]
				r = valList[1]
			elif flowVal == 2:
				n = valList[2]
				r = valList[3]
			elif flowVal == 3:
				n = valList[4]
				r = valList[5]
				#r = valList[3]  #DNR version error

			row[2] = n
			row[3] = r
			cursor.updateRow(row)
	del cursor

	# the next steps create the grids for the non-channel areas (based on land-use) and merge them with the grids for the Channel N and R to produce combined grids of N and R
	arcpy.AddMessage(' - Merging channel hydraulic radius with landcover hydraulic radius ...')
	rGridSV = Lookup(combGrid,'hydr')
	rGridSV.save(rGrid)
	arcpy.Mosaic_management(chHydRadGrid,rGrid)

	arcpy.AddMessage(' - Merging channel mannings N with landcover mannings N ...')
	nGridSV = Lookup(combGrid,'mn')
	nGridSV.save(nGrid)
	arcpy.Mosaic_management(chNGrid, nGrid)

	arcpy.AddMessage(' - Delete temp stuff ...')
	arcpy.Delete_management(chNGrid) #TODO
	arcpy.Delete_management(combGrid)
	arcpy.Delete_management(chHydRadGrid2)

	return rGrid,nGrid

def calcVelocity(scratchDir,tGrid,ds_tt,us_tt,slope2,rGrid,nGrid, xyUnits,cellSize):
	arcpy.AddMessage(' - Creating flow velocity grid. Values in this grid represent the time it takes for water to flow across the grid cell.')
	arcpy.AddMessage(' - Velocity for each cell is calculated based on the equation: (1.49 * Hydraulic Radius^0.667 * slope^0.5) / Mannings N')
	arcpy.AddMessage(' - Velocity in this equation is measured as feet/second.')
	velGrid = (1.49 * Power(Raster(rGrid),0.667) * Power((Raster(slope2)/ 100),0.5)) / Raster(nGrid)
	velGrid.save(scratchDir + os.sep + 'velGrid')

	arcpy.AddMessage(' - Calculating per cell travel time using the velocity grid and converting velocity to - meters/second ...')
	# the velocity grid contains values represent feet/second. To be used in the FlowLength tool the values have
	# to be converted to meters/feet per second - depending on the x,y units of the projected data.
	if xyUnits == 'METER':
		exp = 1 / (velGrid * 0.3048)
	else:
		exp = 1 / + velGrid

	exp.save(tGrid)
	arcpy.Delete_management(velGrid)

def ValidateGridSizeMatch(InputRasterList):
	#Validate grid size
	blnResults = False
	i = 0
	for iRaster in InputRasterList:
		GetCellSize = arcpy.GetRasterProperties_management (iRaster, "CELLSIZEX")
		GridSize = float(GetCellSize.getOutput(0))
		if i <> 0:
			if preGridSize <> GridSize:
				blnResults = True
				break
		preGridSize = GridSize
		i += 1

	if blnResults:
		arcpy.AddError("The grid cell size does not match among all the input rasters!")
	return blnResults

if __name__ == '__main__': main()