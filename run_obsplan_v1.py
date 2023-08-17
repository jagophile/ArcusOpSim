import obsplan
import datetime
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

brightObjects = read_bright_objects('inputs/RA_Dec Bright Objects.txt')
print(brightObjects[:5])

## INPUT SETTINGS
TIMESTEP = 100.0

# declare the observation plan
o = obsplan.ObsPlan()

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

NTARG = 100
TOTOBS=42e6
MINOBS=5000
MAXOBS=1e6

totobs = 0

rng = numpy.random.default_rng(seed=1234)
rntimes = rng.lognormal(size=NTARG)
rntimes = rntimes *(TOTOBS/sum(rntimes))
rntimes[rntimes<MINOBS] = MINOBS
rntimes[rntimes>MAXOBS] = MAXOBS


for i in range(len(rntimes)):
  RA = numpy.random.random()*360
  Dec = (numpy.random.random()-0.5) * 180
  time = rntimes[i]
  name = "Target_%03i"%(i)
  t = obsplan.Target(RA, Dec, time, name=name)

  o.add_target(t)

###
print("Added all the targets")


vis = o.restrictions[0].applyRestriction(target)
vis2 = o.restrictions[1].applyRestriction(target)


fig = pylab.figure()
fig.show()
ax1 = fig.add_subplot(411)
ax2 = fig.add_subplot(412, sharex=ax1)
ax3 = fig.add_subplot(413, sharex=ax1)
ax4 = fig.add_subplot(414, sharex=ax1)

ax1.plot(brightObjects['Time'], brightObjects['Sun_RA'])
ax1.plot(brightObjects['Time'], brightObjects['Moon_RA'])
ax2.plot(brightObjects['Time'], brightObjects['Sun_Dec'])
ax2.plot(brightObjects['Time'], brightObjects['Moon_Dec'])
ax3.plot(o.timeList, vis)
ax3.plot(o.timeList, vis2)
ax3.plot(o.timeList, vis&vis2)
ax4.plot(o.timeList, o.restrictions[0].tmp)
ax4.plot(o.timeList, o.restrictions[1].tmp)

pylab.draw()
zzz=input()



#def __init__(self, restrictionName, restrictionData, setupRoutine, implementationRoutine, observationPlan = None):
