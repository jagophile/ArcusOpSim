import numpy
import datetime
import calendar
import pylab
import scipy.interpolate as spint
import astropy.coordinates as coord
import astropy.units as u
from targetClass import *

# Some constants

NOACTIONSET = 0

def calc_distance(x, y, startx, starty):
  # calculate the distance in degrees between from startx, starty to
  # x,y. x,y, can be arrays or individual points.

  dx = numpy.abs(x-startx)
  cosd = numpy.cos(numpy.deg2rad(90-y))*\
         numpy.cos(numpy.deg2rad(90-starty)) +\
         numpy.sin(numpy.deg2rad(90-y))*\
         numpy.sin(numpy.deg2rad(90-starty))*\
         numpy.cos(numpy.deg2rad(x-startx))
  d = numpy.arccos(cosd)

  r= numpy.rad2deg(d)
  for i in range(len(r)):
    if not(numpy.isfinite(r[i])):
      r[i] = 1.0
#  print(max(r), min(r))
#  zzz=input('hmm1')
  return r


def toTimestamp(d):
  return calendar.timegm(d.timetuple())

def read_targetlist(filename):
  dtype1 = numpy.dtype({'names':['target','ra','dec','obstime'],\
                        'formats':['|S20',float,float,float]})
  dat = numpy.genfromtxt('targetlist_v207.csv',dtype=dtype1, autostrip=True,\
                         filling_values=numpy.nan, invalid_raise=False)



  dtype2 = numpy.dtype({'names':['target','ra','dec','obstime','observed','observable'],\
                        'formats':['|S20',float,float,float, bool, bool]})

  obslist = numpy.zeros(len(dat), dtype=dtype2)
  for name in dtype1.names:
    obslist[name]=dat[name]
  obslist['observed'][:]=False
  obslist['observable'][:]=False

class ObsPlan():
  """
  The main class. Make one of these!
  """

  def __init__(self):
    self.targetList = []
    self.restrictionList = []

  def set_timestep(self, timestep, tstart, tend):
    """
    Set the timestep of the observation plan

    timestep: int
      Timestep in seconds
    tstart : datetime
      Time to start plan
    tend : datetime
      Time to end plan
    """
    ts = int(timestep * 1000)
    self.timestep = numpy.timedelta64(ts, 'ms')
    ts = int(timestep * 1)
    self.timestep = numpy.timedelta64(ts, 's')
    # check inputs
    print(type(tstart))
    if type(tstart) is numpy.datetime64:
      pass
    else:
      raise Exception("Error: tstart is not a datetime")

    if type(tend) is numpy.datetime64:
      pass
    else:
      raise Exception("Error: tend is not a datetime")

    self.tstart = tstart
    self.tend = tend

    # set up an array of times
    self.timeList = numpy.arange(tstart, tend, self.timestep)

    self.nextTimeList = 0

    # set up data arrays to store the satellite's actions
    self.observatoryActions = numpy.zeros(len(self.timeList), dtype = int)
    self.pointing = numpy.zeros(len(self.timeList), dtype = int)

    # arrays to store data
    self.dataRates = numpy.zeros(len(self.timeList), dtype = float)
    self.dataBuffer = numpy.zeros(len(self.timeList), dtype = float)

    # XXX PLACEHOLDER: Add instruments
    self.SLEWRATE=10/60 # deg/min


  def add_restriction(self, restrictionType, restrictionName, restrictionData):
    """
    Add a restriction to the ObsPlan's restrictionList

    PARAMETERS
    ----------
    RestrictionType : class
      The class of restriction you are adding. Should be a subclass of Restriction
    restrictionName : string
      The name of the restriction
    restrictionData : dict
      Dictionary of infomation to be passed to restrictionType function

    RETURNS
    -------
    None
    """

    self.restrictionList.append(restrictionType(restrictionName, restrictionData, observationPlan = self))

  def add_target(self, Target):
    """
    Add a target to the ObsPlan's targetList

    PARAMETERS
    ----------

    Target: Target
      A Target class object

    RETURNS
    -------
    None
    """
    
    Target.ObsPlan = self
    self.targetList.append(Target)


  def calc_target_distances(self):
    t = []
    for target in self.targetList:
      t.append(target.Coord)
    self.TargetCoords = coord.SkyCoord(t)
    for target in self.targetList:
      target.distances = numpy.array(numpy.ceil(target.Coord.separation(self.TargetCoords).value/self.SLEWRATE), dtype=int)


  def recalc_target_restrictions(self, startTimeIndex=0, firstCalc=False):
    for target in self.targetList:
      target.calc_visibility(self.restrictionList, startTimeIndex=startTimeIndex, firstCalc=firstCalc)


  def initialize_run(self):

    # calculate the distance between all targets
    self.calc_target_distances()
    
    # calculate when each target is visible
    self.recalc_target_restrictions(startTimeIndex=0, firstCalc=True) 


  def iterate_timestep(self, timeIndex=0):
    # so what happens in here???

    

    for target in self.targetList:
      # update the visibility of all targets
      target.update_soft_visibility(self, timeIndex)

      





class Restriction():

  def __init__(self, restrictionName, restrictionData, observationPlan = None):
    """
    An observatory/instrument restriction

    restrictionName : string
      The name of the restriction
    restrictionData : dict
      A dictionary containing any of the restriction information required
    setupRoutine: func
      A function for handling any initial setup
    implementationRoutine: func
      A function for calculating the restriction. Should return True or False for all
      time periods.
    observationPlan : ObsPlan
      The parent observation plan object
    """

    self.restrictionName = restrictionName
    self.restrictionData = restrictionData

    self.observationPlan=observationPlan


    # call the setup routine
    self.setupRoutine(self)


  def applyRestriction(target):
    """
    Apply restrictions to the target

    target : Target
      The target we are applying restrictions to.
    """
    visibility = self.implementationRoutine(self, target)

    return(visibility)


class ObservatoryMode():

  def __init__(self, obsModeName, modeData, observationPlan = None):
    """
    An observatory/instrument mode

    restrictionName : string
      The name of the restriction
    restrictionData : dict
      A dictionary containing any of the restriction information required
    setupRoutine: func
      A function for handling any initial setup
    implementationRoutine: func
      A function for calculating the restriction. Should return True or False for all
      time periods.
    observationPlan : ObsPlan
      The parent observation plan object
    """

    self.obsModeName = obsModeName
    self.modeData = modeData

    self.observationPlan=observationPlan


    # call the setup routine
    self.setupRoutine(self)


  def applyMode(target):
    """
    Apply restrictions to the target

    target : Target
      The target we are applying restrictions to.
    """
    visibility = self.implementationRoutine(self, target)

    return(visibility)

class KeepOutRestriction(Restriction):

  def __init__(self, restrictionName, restrictionData, observationPlan = None):
    """
    An observatory/instrument restriction

    restrictionName : string
      The name of the restriction
    restrictionData : dict
      A dictionary containing any of the restriction information required
    setupRoutine: func
      A function for handling any initial setup
    implementationRoutine: func
      A function for calculating the restriction. Should return True or False for all
      time periods.
    observationPlan : ObsPlan
      The parent observation plan object
    """

    self.restrictionName = restrictionName
    self.restrictionData = restrictionData

    self.observationPlan=observationPlan

 #   fig = pylab.figure()
 #   fig.show()
 #   ax1 = fig.add_subplot(211)
 #   ax2 = fig.add_subplot(212, sharex=ax1)

    # call the setup routine
    if observationPlan is not None:
      timeList_out = numpy.array(observationPlan.timeList, dtype=float)
      timeList_in = numpy.array(restrictionData['Time'], dtype=float)
      print(timeList_out[:5])
      print(timeList_in[:5])

      tmpRA = self.restrictionData['RA']*1.0



      # corrections to deal with periodicity (RA going past 360 deg)
      # keep adding/subtracting factors of 360 until it's smooth


      i = numpy.where(tmpRA[1:] < tmpRA[:-1]-270)[0]
      for ii in i:
        tmpRA[ii+1:]+=360

      i = numpy.where(tmpRA[1:] > tmpRA[:-1]+270)[0]
      for ii in i:
        tmpRA[ii+1:]-=360

      RA_out = numpy.interp(timeList_out, timeList_in, tmpRA)

      # remove the 360 deg factor
      RA_out = RA_out % 360
      Dec_out = numpy.interp(timeList_out, timeList_in, self.restrictionData['Dec'])
      self.RA = RA_out
      self.Dec = Dec_out
      self.Coord = coord.FK5(self.RA*u.deg, self.Dec*u.deg)

    self.Time = self.restrictionData['Time']
    self.keepOutAngle = self.restrictionData['Anglelimits']

    self.needsRecalculated = True


  def applyRestriction(self, target):
    """
    Apply restrictions to the target

    target : Target
      The target we are applying restrictions to.
    """
    

    tmp = calc_distance(self.RA, self.Dec, target.RA, target.Dec)
    tt = coord.FK5(target.RA*u.deg, target.Dec*u.deg)
    tmp2 = tt.separation(self.Coord)



#    fig = pylab.figure()
#    fig.show()
#    ax1 = fig.add_subplot(211)
#    ax2 = fig.add_subplot(212, sharex=ax1)

#    ax1.plot(self.observationPlan.timeList, tmp)
#    ax1.plot(self.observationPlan.timeList, tmp2)


    self.tmp = tmp

    vis = (tmp > self.keepOutAngle[0]) &\
                      (tmp < self.keepOutAngle[1])
#
#
#    vis2 = (tmp2 > self.keepOutAngle[0]*u.deg) &\
#                      (tmp2 < self.keepOutAngle[1]*u.deg)
#
#    ax2.plot(self.observationPlan.timeList, vis)
#    ax2.plot(self.observationPlan.timeList, vis2)

#    pylab.draw()
#    zzz=input()

    self.needsRecalculated = False

    return(vis)




