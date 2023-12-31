import numpy
import datetime
import calendar
import pylab
import scipy.interpolate as spint
import astropy.coordinates as coord
import astropy.units as u

NOACTIONSET = 0

class Target():
  """
  Individual targets
  """

  def __init__(self, RA, Dec, reqObsTime, name="None", scienceGoal=-1, ObsPlan=False, idnumber=0):
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
    self.copmleteObsTime = 0
    self.targetComplete=False
    
    self.name = name
    try:
      tmp = iter(scieneGoal)

    except:
      scienceGoal = [scienceGoal]
    self.scienceGoal = scienceGoal
    self.ObsPlan = ObsPlan
    self.idNumber=idnumber
    self.restrictions=[]


  def set_time_intervals(self, timeList):
    # all the times go here
    self.timeList = timeList
    # this is for the hardwired, permanent, what can be seen
    self.visibility = numpy.ones(len(timeList), dtype=bool)

    # this is the visibilty, plus periods ruled out for other activities
    self.startVisibility = numpy.ones(len(timeList), dtype=bool)

  def calc_visibility(self, restrictions, startTimeIndex=0, firstCalc=False):
    """
    Calculate the times on which the target is visible

    restrictions : observatory_restrictions
      The various restrictions which are enforced

    startTimeIndex : int
      Only calculate visibility after this time, all before is
      set to False (used to ignore visibility in the past)
    """

    

    # apply general restrictions

    for r in restrictions:
      if ((r.needsRecalculated != False) | (firstCalc==True)):
        tmp = r.applyRestriction(self)
        self.visibility[tmp==False] = False
    
    # set the visibility to match the results
#    self.visibility= isgood

    # set the visibility to False for all the earlier observations
    self.visibility[:startTimeIndex] = False

    # now calculate target specific restrictions
#    self.visibility[self.ObsPlan.observatoryActions != NOACTIONSET] = False

  def calc_start_visibility(self, timeIndex):
    """
    This updates the softVisibility based on the current observation plan

    The idea is this reflects potential starting points of observations instead
    of hardwired times
    """
    self.startVisibility[:timeIndex]=False
    
    
    for r in self.restrictions:
      tmp = r.applyRestriction(self)
      self.startVisibility[tmp==False] = False

#    self.startVisibility = numpy.zeros(len(timelist), dtype=bool)

  def add_restriction(self, restriction):
    self.restrictions.append(restriction)

class TargetType():
  """
  Class for groups of targets
  """

  def __init__(self, name="None", scienceGoal=-1, ObsPlan=False):
    """
    name : string
      target name
    scienceGoal : int or [int, int,...]
      Science goal target(s) target applies to
    ObsPlan : ObsPlan
      The parent observation plan this target is within
    """

    self.name=name
    self.scienceGoal = scienceGoal
    self.ObsPlan = ObsPlan
    self.restrictions = []
    seld.targets = []

  def add_target(self, target):
    self.targets.append(target)
  
  def add_restriction(self, restriction):

    self.restrictions.append(restrictions)

  def apply_restrictions():
    for r in self.restrictions:
      for t in self.targets:
        r.apply_restriction(t)


    
  
