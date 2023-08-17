import numpy
import datetime
import calendar
import pylab
import scipy.interpolate as spint
import astropy.coordinates as coord
import astropy.units as u


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
  print(max(r), min(r))
  zzz=input('hmm1')
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
    self.restrictions = []

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
    self.dataRates = numpy.zeros(len(self.timeList), dtype = numpy.float)
    self.dataBuffer = numpy.zeros(len(self.timeList), dtype = numpy.float)

    self.SLEWRATE=10/60 # deg/min


    # XXX PLACEHOLDER: Add instruments

  def add_restriction(self, RestrictionType, restrictionName, restrictionData, observationPlan = None):

    self.restrictions.append(RestrictionType(restrictionName, restrictionData, observationPlan = self))

  def add_target(self, Target):
    Target.ObsPlan = self
    self.targetList.append(Target)


  def calc_target_distances(self):
    t = []
    for target in self.targetList:
      t.append(target.coord)
    self.TargetCoords = coord.SkyCoord(t)
    for target in self.targetList:
      target.distances = numpy.array(numpy.ceil(target.coord.separation(self.TargetCoords).values/self.slewrate), dtype=int)


  def recalc_target_restrictions(self, startTimeIndex=0):
    for target in self.targetList:
      r = calc_visibility(self.restrictions, startTimeIndex=startTimeIndex)
      r[:startTimeIndex]=False
      target.visibility=r
      target.totalVisibility = sum(target.visiblity)


class Target():
  """
  Individual targets
  """

  def __init__(self, RA, Dec, reqObsTime, name="None", scienceGoal=-1, ObsPlan=False):
    """
    RA: float
      Right Acension (J2000, decimal)
    Dec: float
      Declination (J2000, decimal)
    reqObsTime : float
      Total required observation time (s)
    name : string
      target name
    scienceGoal : int or [int, int,...]
      Science goal target(s) target applies to
    ObsPlan : ObsPlan
      The parent observation plan this target is within
    """

    self.RA = RA
    self.Dec = Dec
    self.Coord = coord.SkyCoord(RA*u.deg, Dec*u.deg)
    self.reqObsTime = reqObsTime
    self.name = name
    try:
      tmp = iter(scieneGoal)

    except:
      scienceGoal = [scienceGoal]
    self.scienceGoal = scienceGoal
    self.ObsPlan = ObsPlan


  def set_time_intervals(timeList):
    # all the times go here
    self.timeList = timelist
    self.visibility = numpy.zeros(len(timelist), dtype=bool)

  def calc_visibility(restrictions, startTimeIndex=0):
    """
    Calculate the times on which the target is visible

    restrictions : observatory_restrictions
      The various restrictions which are enforced

    startTimeIndex : int
      Only calculate visibility after this time, all before is
      set to False (used to ignore visibility in the past)
    """

    isgood = numpy.ones(len(self.timeList), dtype=bool)


    # apply general restrictions

    for r in restrictions:
      tmp = restrictions.applyRestrictions(self)
      isgood[tmp==False] = False

    # set the visibility to match the results
    self.visibility= isgood

    # set the visibility to False for all the earlier observations
    self.visibility[:startTimeIndex] = False


    # now calculate target specific restrictions
    self.visibility[self.ObsPlan.observatoryActions != NOACTIONSET] = False





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



  def applyRestriction(self, target):
    """
    Apply restrictions to the target

    target : Target
      The target we are applying restrictions to.
    """

    tmp = calc_distance(self.RA, self.Dec, target.RA, target.Dec)
    tt = coord.FK5(target.RA*u.deg, target.Dec*u.deg)
    tmp2 = tt.separation(self.Coord)

    fig = pylab.figure()
    fig.show()
    ax1 = fig.add_subplot(211)
    ax2 = fig.add_subplot(212, sharex=ax1)

    ax1.plot(self.observationPlan.timeList, tmp)
    ax1.plot(self.observationPlan.timeList, tmp2)


    self.tmp = tmp

    vis = (tmp > self.keepOutAngle[0]) &\
                      (tmp < self.keepOutAngle[1])


    vis2 = (tmp2 > self.keepOutAngle[0]*u.deg) &\
                      (tmp2 < self.keepOutAngle[1]*u.deg)

    ax2.plot(self.observationPlan.timeList, vis)
    ax2.plot(self.observationPlan.timeList, vis2)

    pylab.draw()
    zzz=input()
    return(vis)




