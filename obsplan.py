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
MODEWAIT = -1
MODESLEW = -2
MODESETTLE = -3

  

  
  
    
    

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
    self.targetClasses = []
    self.restrictionList = []
    self.STARTTIMECUTOFF = 5*3600 # 5 hours
    self.modes = []

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
    ts = int(timestep * 1)
    self.timestep = numpy.timedelta64(ts, 's')

    self.STARTTIMECUTOFFSTEPS = int(numpy.ceil(numpy.timedelta64(self.STARTTIMECUTOFF, 's')/ self.timestep))
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
    self.RA = numpy.zeros(len(self.timeList), dtype = float)
    self.Dec = numpy.zeros(len(self.timeList), dtype = float)

    # arrays to store data
    self.dataRates = numpy.zeros(len(self.timeList), dtype = float)
    self.dataBuffer = numpy.zeros(len(self.timeList), dtype = float)

    # XXX PLACEHOLDER: Add instruments
    self.SLEWRATE=5/60 # deg/min
    self.SETTLETIME=3600 # s
    self.SETTLETIMESTEPS = int(numpy.ceil(numpy.timedelta64(self.SETTLETIME, 's')/ self.timestep))


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

  def add_compulsory_action(self, actionType, restrictionName, restrictionData):
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


  def calc_target_restrictions(self, startTimeIndex=0, firstCalc=False):
    for target in self.targetList:
      target.calc_visibility(self.restrictionList, startTimeIndex=startTimeIndex, firstCalc=firstCalc)
      target.calc_start_visibility(startTimeIndex)

  def initialize_run(self):

    for target in self.targetList:
      target.set_time_intervals(self.timeList)

    # calculate the distance between all targets
    self.calc_target_distances()
    
    # calculate when each target is visible
    self.calc_target_restrictions(startTimeIndex=0, firstCalc=True) 

    # add compulsory actions
    # PLACEHOLDER FOR e.g. SDDL, ECLIPSE
#    self.calc_required_actions(startTimeIndex=0, firstCalc=True) 


  def iterate_timestep(self, timeIndex=0):
    # so what happens in here???

    total_starttime=numpy.zeros(len(self.targetList), dtype=int)
    total_available_obstime=numpy.zeros(len(self.targetList), dtype=int)
    total_remaining_obstime_required=numpy.zeros(len(self.targetList), dtype=int)
    next_starttime=numpy.zeros(len(self.targetList), dtype=int)
    next_obslength=numpy.zeros(len(self.targetList), dtype=int)
    
    

    for itarget, target in enumerate(self.targetList):
      # update the visibility of all targets

      # ignore completed targets
      if target.targetComplete:
        continue

      # slew onto target
      target.update_start_visibility(self, timeIndex) # this will take ages!
      total_starttime[itarget] = sum(target.startVisibility)
      total_available_obstime[itarget] = sum(target.visibility)
      total_remaining_obstime[itarget] = target.reqObsTime-target.completeObsTime
      ist = numpy.where(target.startVisibility[itarget:]==True)[0][0]+timeIndex
      next_starttime[itarget] = ist
      
      #
      for i in range(ist, len(self.timeList)):
        if target.startVisibility[i]==False:
          break
      next_obslength[itarget] = i-ist

    ## IN HERE CALL RANKING ALGORITHM

    rank = self.rank_targets(total_starttime, total_available_obstime,
                        total_remaining_obstime, next_starttime, next_obslength, len(self.timeList)-timeIndex)
                        
    if len(rank)==0:
      self.observatoryActions[timeIndex:timeIndex+self.STARTTIMECUTOFFSTEPS]=MODEWAIT
      timeOut = timeIndex+self.STARTTIMECUTOFFSTEPS
    else:
      # pick something that sucks eggs.
      if len(rank)==1:
        chosenTarget = rank[0]['index']
      else:
        rng = numpy.random.default_rng()
        s=sum(rank['score'])
        while True:
          exp = rng.exponential(s)
          if exp < 2*s:
            break
        s =0
        i_s = 0
        while exp > 2*s:
          s+=2*rank['score'][i_s]*2
          i_s+=1
        chosenTarget = rank[i_s-1]['index']

      # implement the observation
      self.make_observation(chosenTarget, timeIndex)
        
        
  def make_observation(self, chosenTarget, timeIndex):
    # find the matching target

    target = self.targetList(chosenTarget)

    # slew to target
    # find the last assigned pointing

    if timeIndex > 0:
      pointing = self.pointing[timeIndex-1]
      if pointing == 0:
        tmptimeIndex = timeIndex-1
        while tmptimeIndex > 0:
          tmptimeIndex-1
          pointing = self.pointing[tmptimeIndex]
          if pointing != 0: break

        # fill in the pointing up until now
        pointing[tmptimeIndex+1:timeindex-1] = pointing[tmptimeIndex]

      # calculate slew time
      slewtime = int(numpy.ceil(pointing[timeIndex-1].separation(target.Coord)))

      # Put in the slew
      self.observatoryActions[timeIndex:timeIndex+slewtime] = MODESLEW
      timeIndex += slewtime
      
      # Put in the settle
    self.observatoryActions[timeIndex:timeIndex+self.SETTLETIMESTEPS] = MODESETTLE
    timeIndex += self.SETTLETIMESTEPS
      
      # Make the observation

    while( (target.completeObsTime < target.reqObsTime) & #target not complete
           (target.visibility[timeIndex]==True) & # target is visible
           (self.observatoryActions[timeIndex] == 0) & # observatory not doing something else
           (timeIndex < len(self.observatoryActions))): # haven't overrun end of mission
      self.observatoryActions[timeIndex] = target.idNumber
      self.pointing[timeIndex] = target.Coord
      timeIndex+=1

    # update the visibility chart for every target
    if target.completeObsTime >= target.reqObsTime : target.targetComplete=True
    
    for target in self.targetList:
      target.calc_start_visibility(timeIndex)
      
#    for target in self.targetList:
#      t.append(target.Coord)
#    self.TargetCoords = coord.SkyCoord(t)
#    for target in self.targetList:
#      target.distances = numpy.array(numpy.ceil(target.Coord.separation(self.TargetCoords).value/self.SLEWRATE), dtype=int)

          
        
  def rank_targets(self, total_starttime, total_available_obstime,
                   total_remaining_obstime, next_starttime, next_obslength, remaining_time):
  
    """
    This is the ranking algorithm. Whee. All parameters are ntarget arrays
  
    PARAMETERS
    ----------
  
    total_starttime : int
      number of timesteps at which observations could start
    total_available_obstime : int
      number of timesteps at which target is observable
    total_remaining_obstime : float
      amount of time remaining to observe (s)
    next_starttime : int
      when the target is next observable
    next_obslength : int
      number of timesteps next observation could take place for
    remaining_time : int
      total number of timesteps left
  
    RETURNS
    -------
    array(int, float) Ranked list of potential targets and theis scores
    """
  
    ### check for emergency targets
    iurgent = numpy.where(total_remaining_obstime/total_available_obstime < 2.0)[0]
  
    score = total_remaining_obstime/total_starttime
    score[total_remaining_obstime==0]= 0.0
  
    ###
    scoredtype = numpy.dtype({'names':['index','score'],
                              'formats':[int, float]})
  
    if len(iurgent) > 0:
      score[~iurgent]=0
  
      rank = numpy.zeros(len(iurgent), dtype=scoredtype)
  
      for iu in iurgent:
        rank[iu]['index'] = iu
        rank[iu]['score'] = score[iu]
  
      rank.sort(score)
      rank=rank[::-1]
      return(rank)
  
    # check for constrained observations
    iconst = numpy.where(total_starttime < 0.05*remaining_time)[0]
  
    if len(iconst) > 0:
      score[~iconst] = 0
      rank = numpy.zeros(len(iconst), dtype=scoredtype)
      for ic in iconst:
        rank[ic]['index'] = ic
        rank[ic]['score'] = score[ic]
      rank.sort(score)
      rank=rank[::-1]
      return(rank)
  
    # score the remaining observations
  
  
    ir = numpy.where((score > 0) & (next_starttime < self.STARTTIMECUTOFFSTEPS))[0]
  
    if len(ir) == 0:
      rank = numpy.zeros(0, dtype=scoredtype)
      return(rank)
  
    score[~ir] = 0
    rank = numpy.zeros(len(ir), dtype=scoredtype)
    for iir in ir:
      rank[iir]['index'] = iir
      rank[iir]['score'] = score[iir]
    rank.sort(score)
    rank=rank[::-1]
    return(rank)
    
    
      

      





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


  def apply_restriction(target):
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


  def apply_mode(target, tstart, tend):
    """
    Apply restrictions to the target

    target : Target
      The target we are applying restrictions to.
    """
    visibility = self.implementationRoutine(self, target)

    return(visibility)



class ObservatoryModeWait(ObservatoryMode):

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
    self.modeIdentifier = MODEWAIT


  def apply_mode(istart, ilen):
    """
    Apply restrictions to the target

    target : Target
      The target we are applying restrictions to.
    """
    self.observationPlan.observatoryActions[istart:istart+ilen] = self.modeIdentifier

    
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

class SlewRestriction(Restriction):

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

    self.slewRate = restrictionData['slewRate']
    self.settle = restrictionData['settleTime']
    self.needsRecalculated = True


  def applyRestriction(self, target):
    """
    Apply restrictions to the target

    target : Target
      The target we are applying restrictions to.
    """
    timestep = self.ObsPlan.nextTimeList
#    if timestep > 0:
      # look backwards
      
      
    tt = coord.FK5(target.RA*u.deg, target.Dec*u.deg)
#    tmp2 = target.Coord.separation(self.observationPlan.)

#
#    tmp = calc_distance(self.RA, self.Dec, target.RA, target.Dec)
#    tt = coord.FK5(target.RA*u.deg, target.Dec*u.deg)
#    tmp2 = tt.separation(self.Coord)
    



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

class MinObsRestriction(Restriction):

  def __init__(self, restrictionName, restrictionData, observationPlan = None):
    """
    Minimum required observing time

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

    self.minObsTime = restrictionData['minObsTime']
    self.minObsTimeSteps = int(numpy.ceil(numpy.timedelta64(self.minObsTime,'s')/self.observationPlan.timestep))
    
  def applyRestriction(self, target):
    """
    Apply restrictions to the target

    target : Target
      The target we are applying restrictions to.
    """

    isgood = numpy.ones(len(target.visibility), dtype=bool)
    
    isgood[target.visibility!=True] = False

    indexes = numpy.where(((isgood[1:]+1) - (isgood[:-1]+1)) < 0)[0]
    
    for ii in indexes:
      isgood[ max([0, ii+1-self.minObsTimeSteps]) : ii+1]=False

    return(isgood)



