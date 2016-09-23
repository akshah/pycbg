from __future__ import division
from mpl_toolkits.basemap import Basemap
from matplotlib import pyplot as plt
from shapely.geometry import Point,LineString,Polygon
from shapely.ops import polygonize,cascaded_union
import shapely.ops as ops
from shapely.ops import transform
from functools import partial
import traceback
import numpy as np
import json
import pyproj
import os
import pickle
from contextlib import closing
import pymysql

class cbg():

    def __init__(self):
        return None

    def getInputConstraints(self,measurementData):
        inputPolygons=[]
        probeIDLocationDict=self.readProbeLocationInfo()
        for trace in measurementData:
            rttVals=[]
            try:
                for rttDict in trace["result"]:
                    rttVals.append(rttDict["rtt"])
                dst=self.rttToDistance(float(min(rttVals)))
                if dst<10000:#Extremly Large, Will mostly not happen
                    ptLat=float(probeIDLocationDict[trace["prb_id"]]["lat"])
                    ptLong=float(probeIDLocationDict[trace["prb_id"]]["lon"])
                    #Check lat lons
                    if ptLat<=90 and ptLat>=-90 and ptLong<=180 and ptLong>=-180:
                        inputPolygons.append(self.latlonbuffer(ptLat,ptLong,dst))
            except KeyError:
                continue
        #print('{0} input constraints'.format(len(inputPolygons)))
        return inputPolygons

    def latlonString(self,lat, lon):
        try:
            WGS84 = pyproj.Proj(init='epsg:4326')
            proj4str = '+proj=aeqd +lat_0=%s +lon_0=%s +x_0=0 +y_0=0' % (lat, lon)
            AEQD = pyproj.Proj(proj4str)
            project = partial(pyproj.transform, AEQD, WGS84)
            bufferPoint=Point(0, 0)#.buffer(0.000001)
            P = transform(project, bufferPoint)
        except:
            traceback.print_exc()
        return P


    def latlonbuffer(self,lat, lon,radius_m):
        try:
            WGS84 = pyproj.Proj(init='epsg:4326')
            proj4str = '+proj=aeqd +lat_0=%s +lon_0=%s +x_0=0 +y_0=0' % (lat, lon)
            AEQD = pyproj.Proj(proj4str)
            project = partial(pyproj.transform, AEQD, WGS84)
            P = transform(project, Point(0, 0).buffer(radius_m))
        except:
            traceback.print_exc()
        return P

    def checkIfCrossesBoundaries(self,inPoly):
        ##
        #Map boundaries
        rtTopCorner=self.latlonString(88,178)
        rtBottomCorner=self.latlonString(-88,178)
        ltTopCorner=self.latlonString(88,-178)
        ltBottomCorner=self.latlonString(-88,-178)
        ltMidPt=self.latlonString(0,-178)
        rtMidPt=self.latlonString(0,178)
        topMidPt=self.latlonString(178,0)
        bottomMidPt=self.latlonString(-178,0)
        rtLine=LineString([rtBottomCorner,rtTopCorner])
        topLine=LineString([ltTopCorner,rtTopCorner])
        bottomLine=LineString([rtBottomCorner,ltBottomCorner])
        ltLine=LineString([ltBottomCorner,ltTopCorner])
        horMidLine=LineString([ltMidPt,rtMidPt])
        verMidLine=LineString([topMidPt,bottomMidPt])
        boundaryLinesList=[{'rtLine':rtLine},{'topLine':topLine},{'bottomLine':bottomLine},{'ltLine':ltLine},{'horMidLine':horMidLine},{'verMidLine':verMidLine}]
        ##
        linesIntersectsWith=[]
        intersectionFlag=False
        for boundaryLineDict in boundaryLinesList:
            name=list(boundaryLineDict.keys())[0]
            line=list(boundaryLineDict.values())[0]
            if inPoly.intersects(line):
                linesIntersectsWith.append(name)
                intersectionFlag=True
        return intersectionFlag,linesIntersectsWith

    def checkQuadrants(self,inPoly):
        quadrantsSet=set()
        xs,ys = inPoly.exterior.xy
        for iter in range(0,len(xs)):
            xval=xs[iter]
            yval=ys[iter]
            if xval >=0 and yval>=0:
                quadrantsSet.add(1)
            if xval <0 and yval>0:
                quadrantsSet.add(2)
            if xval <=0 and yval<=0:
                quadrantsSet.add(3)
            if xval >0 and yval<0:
                quadrantsSet.add(4)
        return quadrantsSet

    def plotPolygon(self,A,m):
        #Show A
        lonsA=[]
        latsA=[]
        ax,ay=A.exterior.xy
        for iter in range(0,len(ax)):
            lat=ay[iter]
            lon=ax[iter]
            lonsA.append(lon)
            latsA.append(lat)
        xs,ys = m(lonsA, latsA)
        m.plot(xs, ys,'m-', lw=1.3)

    def plotBoundaries(self,bList,m):
        for boundaryLineDict in bList:
            name=list(boundaryLineDict.keys())[0]
            ln=list(boundaryLineDict.values())[0]
            lonsrtLine=[]
            latsrtLine=[]
            #Show Line
            rlx=ln.xy[0]
            rly=ln.xy[1]
            for iter in range(0,len(rlx)):
                lat=rly[iter]
                lon=rlx[iter]
                lonsrtLine.append(lon)
                latsrtLine.append(lat)
            xs,ys = m(lonsrtLine, latsrtLine)
            #print(xs,ys)
            m.plot(xs, ys,'-', lw=0.5)

    def rttToDistance(self,rtt):
        #Returns KM
        c=3*(10**8)
        #return (4/9)*c*rtt*(10**-3)/1000/2
        return (2/3)*c*rtt*(10**-3)/1000/2

    def showPolygons(self,polyList):
        #Set up figure
        fig=plt.figure(figsize=(15,10))
        ax=fig.add_axes([0.1,0.1,0.8,0.8])
        m=Basemap()
        #m.bluemarble()

        m.drawcoastlines()
        m.fillcontinents()

        #plotBoundaries(boundaryLinesList,m)
        for ply in polyList:
            plotPolygon(ply,m)

        plt.show()

    def haversine(self,lon1, lat1, lon2, lat2):
        """
        Calculate the great circle distance between two points
        on the earth (specified in decimal degrees)
        """
        try:
            # convert decimal degrees to radians
            lon1, lat1, lon2, lat2 = map(np.radians, [lon1, lat1, lon2, lat2])
            # haversine formula
            dlon = lon2 - lon1
            dlat = lat2 - lat1
            a = np.sin(dlat/2)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon/2)**2
            c = 2 * np.arcsin(np.sqrt(a))
            km = 6367 * c
            return km
        except:
            print(lon1, lat1, lon2, lat2)
            traceback.print_exc()

    def geoArea(self,geom):
        geom_area = ops.transform(
        partial(
        pyproj.transform,
        pyproj.Proj(init='EPSG:4326'),
        pyproj.Proj(
            proj='aea',
            lat1=geom.bounds[1],
            lat2=geom.bounds[3])),
    geom)
        return geom_area.area/10**6 #Return in Sq. KM

    def solConstraints(self,inputPolygons):
        unionPolys=[]
        for inPoly in inputPolygons:
            if not inPoly.is_valid:
                exterior = inPoly.exterior
                segments = cascaded_union([exterior,exterior])
                polyParts = list(polygonize(segments))
                for pp in polyParts:
                    unionPolys.append(pp)
            else:
                unionPolys.append(inPoly)

        polyList=[]
        #print('{0} SoL constraints'.format(len(unionPolys)))
        for inPoly in unionPolys:
            crossFlag,linesIntersectsWith=self.checkIfCrossesBoundaries(inPoly)
            quadrants=self.checkQuadrants(inPoly)
            if crossFlag:
                xs,ys = inPoly.exterior.xy

                #1st quadrant
                if 1 in quadrants:
                    lsArray=[]
                    for it in range(0,len(xs)):
                        xval=xs[it]
                        yval=ys[it]
                        if 'rtLine' in linesIntersectsWith or 'ltLine' in linesIntersectsWith:
                            if not xval>=0:
                                xval=180
                        lsArray.append((xval,yval))
                    polyList.append(Polygon(lsArray))

                #2nd quadrant
                if 2 in quadrants:
                    lsArray=[]
                    for it in range(0,len(xs)):
                        xval=xs[it]
                        yval=ys[it]
                        if 'rtLine' in linesIntersectsWith or 'ltLine' in linesIntersectsWith:
                            if not xval<0:
                                xval=-180
                        lsArray.append((xval,yval))
                    polyList.append(Polygon(lsArray))

                #3rd quadrant
                if 3 in quadrants:
                    lsArray=[]
                    for it in range(0,len(xs)):
                        xval=xs[it]
                        yval=ys[it]
                        if 'rtLine' in linesIntersectsWith or 'ltLine' in linesIntersectsWith:
                            if not xval<=0:
                                xval=-180
                        lsArray.append((xval,yval))
                    polyList.append(Polygon(lsArray))

                #4th quadrant
                if 4 in quadrants:
                    if 'rtLine' in linesIntersectsWith:
                        print('Intersects with right line')
                    lsArray=[]
                    for it in range(0,len(xs)):
                        xval=xs[it]
                        yval=ys[it]
                        if 'rtLine' in linesIntersectsWith or 'ltLine' in linesIntersectsWith:
                            if not xval>0:
                                xval=180
                        lsArray.append((xval,yval))
                    polyList.append(Polygon(lsArray))
            else:
                polyList.append(inPoly)

        return polyList

    def loadAllCities(self):
        retDict={}
        if os.path.exists('data/cityInfo.pickle'):
            retDict=pickle.load(open('data/cityInfo.pickle','rb'))
        else:
            #Prepare DB info
            db = pymysql.connect(host="proton.netsec.colostate.edu",
                             user="netsecstudent",
                             passwd="n3ts3cL@bs",#Read only account
                             db="caida_geonames")
            with closing( db.cursor() ) as cur:
                try:
                    cur.execute('select city,region,country,latitude,longitude from Location')
                    row=cur.fetchone()

                    while row is not None:
                        city=row[0]
                        country=row[1]
                        region=row[2]
                        lat=row[3]
                        lon=row[4]
                        retDict[city+'|'+region+'|'+country]={"lat":lat,"lon":lon}
                        row=cur.fetchone()
                except Exception:
                   raise Exception('Select Query Failed')
            pickle.dump(retDict,open('data/cityInfo.pickle','wb'))
        return retDict

    def readProbeLocationInfo(self):
        probeIDLocationDict={}
        probesJson=json.load(open('data/ripeProbes.json'))
        for pEntry in probesJson['objects']:
            id=pEntry['id']
            lat=pEntry['latitude']
            long=pEntry['longitude']
            if pEntry['status'] == 1 or pEntry['status'] == 2:
                if isinstance(lat, float) and isinstance(long, float):
                    if id not in probeIDLocationDict.keys():
                        probeIDLocationDict[id]={}
                    locDict={'lat':lat,'lon':long}
                    probeIDLocationDict[id]=locDict
        return probeIDLocationDict

    def getMaxIntersectionRegions(self,polyList):
        intersecDict=[[]]*len(polyList)
        for iter in range(0,len(polyList)):
            p1=polyList[iter]
            intersecDict[iter]=[]
            intersecDict[iter].append(p1)
            for iterIn in range(0,len(polyList)):
                p2=polyList[iterIn]
                insertFlag=False
                for pInSet in intersecDict[iter]:
                    if p2.intersects(pInSet):
                        insertFlag=True
                    else:
                        insertFlag=False
                        break
                if insertFlag:
                    if p2 not in intersecDict[iter]:
                        intersecDict[iter].append(p2)

        maxIntersectionVal=0
        for entry in intersecDict:
            if maxIntersectionVal < len(entry):
                maxIntersectionVal=len(entry)

        polysWithMaxIntersection=[]
        for entry in intersecDict:
            if len(entry)==maxIntersectionVal:
                polysWithMaxIntersection.append(entry)

        intersectionRegions=[]
        for entry in polysWithMaxIntersection:
            interSecRegion=entry[0]
            for itr in range(1,len(entry)):
                interSecRegion=interSecRegion.intersection(entry[itr])
            intersectionRegions.append(interSecRegion)

        '''
        if len(intersectionRegions)>1:
            while True:
                intersectionRegions=self.getMaxIntersectionRegions(intersectionRegions)
                if len(intersectionRegions)<=2:
                    break
        '''
        return intersectionRegions

    def getCities(self,inputPolygons,kmThreshold=50):
        geoCityDict={}
        cityLocs=self.loadAllCities()
        polyList=self.solConstraints(inputPolygons)

        #result = [poly1.intersection(poly2) for poly1,poly2 in  itertools.combinations(polyList, 2) if poly1.intersects(poly2)]
        if len(polyList)>1:
            intersectionRegions=self.getMaxIntersectionRegions(polyList)

            for pp in intersectionRegions:
                if pp:
                    centroidPoint=pp.centroid
                    centroidLon=centroidPoint.x
                    centroidLat=centroidPoint.y
                    for city,locDict in cityLocs.items():
                        cityLat=float(locDict["lat"])
                        cityLon=float(locDict["lon"])
                        #print(cityLat,cityLon)
                        dst=self.haversine(cityLon,cityLat,centroidLon,centroidLat)
                        if dst<=kmThreshold:
                            city,country,region=city.split('|')
                            if country not in geoCityDict.keys():
                                geoCityDict[country]={}
                            if region not in geoCityDict[country].keys():
                                geoCityDict[country][region]={}
                            if city not in geoCityDict[country][region].keys():
                                geoCityDict[country][region][city]={}
                            geoCityDict[country][region][city]=locDict
        return geoCityDict


    def getEstimatedLocation(self,inputPolygons):
        polyList=self.solConstraints(inputPolygons)
        if len(polyList)==0:
            return None,None,None
        centroidsLatList=[]
        centroidsLonList=[]
        areaList=[]
        if len(polyList)>1:
            intersectionRegions=self.getMaxIntersectionRegions(polyList)
            for pp in intersectionRegions[:1]:#For now work with top intersection
                if pp:
                    centroidPoint=pp.centroid
                    centroidLon=centroidPoint.x
                    centroidLat=centroidPoint.y
                    centroidsLatList.append(centroidLat)
                    centroidsLonList.append(centroidLon)
                    areaList.append(self.geoArea(pp))
        if len(centroidsLatList)>0:
            return np.average(np.array(centroidsLatList)),np.average(np.array(centroidsLonList)),np.average(np.array(areaList))
        else:
            return None,None,None