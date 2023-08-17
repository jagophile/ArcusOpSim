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


  def set_time_intervals(self, timeList):
    # all the times go here
    self.timeList = timelist
    self.visibility = numpy.zeros(len(timelist), dtype=bool)
    self.softVisibility = numpy.zeros(len(timelist), dtype=bool)

  def calc_visibility(self, restrictions, startTimeIndex=0, firstCalc=False):
    """
    Calculate the times on which the target is visible

    restrictions : observatory_restrictions
      The various restrictions which are enforced

    startTimeIndex : int
      Only calculate visibility after this time, all before is
      set to False (used to ignore visibility in the past)
    """

    isgood = numpy.ones(len(self.ObsPlan.timeList), dtype=bool)


    # apply general restrictions

    for r in restrictions:
      if ((r.needsRecalculated != False) | (firstCalc==True)):
        tmp = r.applyRestriction(self)
        isgood[tmp==False] = False
    
    # set the visibility to match the results
    self.visibility= isgood

    # set the visibility to False for all the earlier observations
    self.visibility[:startTimeIndex] = False

    # now calculate target specific restrictions
    self.visibility[self.ObsPlan.observatoryActions != NOACTIONSET] = False

  def calc_soft_visibility(self, timeIndex):
    """
    This updates the softVisibility based on the current observation plan

    The idea is this reflects potential starting points of observations instead
    of hardwired times
    """
