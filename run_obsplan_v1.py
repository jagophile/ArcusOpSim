import obsplan
import datetime
import time
import numpy
import pylab
def read_bright_objects(filename):


  def parsetime(v):
    return(numpy.datetime64(datetime.datetime.strptime(v.decode(), '%d %b %Y %H:%M:%S.%f')))

  a = numpy.dtype( {'names': ['Time','Sun_RA','Sun_Dec','SunProbeEarthAngle', 'Radius',\
                                     'Earth_RA','Earth_Dec','SunProbeMooonAngle', \
                                     'Moon_RA','Moon_Dec'],\
                    'formats': ['datetime64[s]', float, float, float, float,\
                                     float, float, float,\
                                     float, float]})

  d = numpy.genfromtxt(filename,\
                       delimiter = [24, 32, 25, 29, 20,34,28,27, 33, 35],\
                       dtype = a,
                       skip_header=6,
                       converters= {0:parsetime})

  return(d)



## INPUT SETTINGS
TIMESTEP = 1000.0 # in s

# declare the observation plan
o = obsplan.ObsPlan()

# read in the bright object locations
brightObjects = read_bright_objects('inputs/RA_Dec Bright Objects.txt')
print(brightObjects[:5])


# declare the timing
o.set_timestep(TIMESTEP, brightObjects['Time'][0], brightObjects['Time'][-1])


# add the restrictions
print(brightObjects.dtype.names)

o.add_restriction(obsplan.KeepOutRestriction, 'SunKeepOut', {'Time':brightObjects['Time'],
                                   'RA':brightObjects['Sun_RA'],
                                   'Dec':brightObjects['Sun_Dec'],\
                                   'Anglelimits':[45.0,150.0]},)

o.add_restriction(obsplan.KeepOutRestriction, 'MoonKeepOut', {'Time':brightObjects['Time'],
                                   'RA':brightObjects['Moon_RA'],
                                   'Dec':brightObjects['Moon_Dec'],\
                                   'Anglelimits':[45.0,360.0]},)

o.add_restriction(obsplan.KeepOutRestriction, 'EarthKeepOut', {'Time':brightObjects['Time'],
                                   'RA':brightObjects['Earth_RA'],
                                   'Dec':brightObjects['Earth_Dec'],\
                                   'Anglelimits':[45.0,360.0]},)




# MAKE FAKE TARGETLIST

NTARG = 50
TOTOBS=10e6
MINOBS=5000
MAXOBS=1e6
NOACTIONSET = 0
totobs = 0

rng = numpy.random.default_rng(seed=1234)
rntimes = rng.lognormal(size=NTARG)
rntimes = rntimes *(TOTOBS/sum(rntimes))
rntimes[rntimes<MINOBS] = MINOBS
rntimes[rntimes>MAXOBS] = MAXOBS

#targetList = o.TargetType(self, name="Generic", scienceGoal=1, ObsPlan=o)

shortObsRestriction = obsplan.MinObsRestriction("ShortObs", {'minObsTime':5000}, observationPlan = o)

for i in range(len(rntimes)):
  RA = numpy.random.random()*360
  Dec = (numpy.random.random()-0.5) * 180
  targtime = rntimes[i]
  name = "Target_%03i"%(i)
  idnumber = i+1
  t = obsplan.Target(RA, Dec, targtime, name=name, idnumber = idnumber)
  t.add_restriction(shortObsRestriction)
  o.add_target(t)

###
print("Added all the targets")

o.initialize_run()

# NOW DOOOO IT
o.iterate_timestep()
fig = pylab.figure()
fig.show()
ax1 = fig.add_subplot(111)

for t in o.targetList:
  it = sum(t.visibility)/len(t.visibility)
  print(it)

for io in range(1,10):
  ax1.plot(o.timeList, (io*2)+o.targetList[io].visibility, drawstyle='steps')
  ax1.plot(o.timeList, (io*2)+o.targetList[io].startVisibility, '--', linewidth=3,drawstyle='steps')
  print(o.targetList[io].RA, o.targetList[io].Dec)


print(o.observatoryActions)

#ax2 = fig.add_subplot(412, sharex=ax1)
#ax3 = fig.add_subplot(413, sharex=ax1)
#ax4 = fig.add_subplot(414, sharex=ax1)

#ax1.plot(brightObjects['Time'], brightObjects['Sun_RA'])
#ax1.plot(brightObjects['Time'], brightObjects['Moon_RA'])
#ax2.plot(brightObjects['Time'], brightObjects['Sun_Dec'])
#ax2.plot(brightObjects['Time'], brightObjects['Moon_Dec'])
#ax3.plot(o.timeList, vis)
#ax3.plot(o.timeList, vis2)
#ax3.plot(o.timeList, vis&vis2)
#ax4.plot(o.timeList, o.restrictions[0].tmp)
#ax4.plot(o.timeList, o.restrictions[1].tmp)

pylab.draw()
zzz=input()



#def __init__(self, restrictionName, restrictionData, setupRoutine, implementationRoutine, observationPlan = None):
