# A simple python resource profiling designed to
#     Use only core python packages (since this was written for code running on large government clusters where there is root access and installing packages is annoying)
#     Work on Windows
#
# The two way to use this class.
#
# First way (look at net memery allocation within function):
#   just change:
#     def function_toatch(dsf,fd,fd):
#   to:
#     @memcounter
#     def function_toatch(dsf,fd,fd):
#
# Second way.. sample usage peoridcally. Just import this file then
#   import memcounter
#   memcounter.start()

import threading
import resource


"""
--------------------- Basic IO Handling Abstractions ------------------------
"""
# Resource usage using the builtin Python resource package
class MemUsageResource:
  def __init__(self):
    self.resultsStale = True
    self.usageData = None
    self.mode = resource.RUSAGE_SELF
    #self.mode = resource.RUSAGE_CHILDREN
    #self.mode = resource.RUSAGE_BOTH
  def queryUsage(self):
    if self.resultsStale:
      self.usageData = resource.getrusage(self.mode)
      self.resultsStale = False
    return
  def getResidentMemUsage(self):
    self.queryUsage()
    # 2 is the same as the famous Unixy 'ru_maxrss' returned by the command 'getrusage' (see man page of getrusage for more infos)
    return self.usageData[2]
  def getSharedMemUsage(self):
    # index 3 is 'ru_ixrss'
    return self.usageData[3]
  def getUnsharedMemUsage(self):
    # index 3 is 'ru_idrss'
    return self.usageData[4]
  def getStackMemUsage(self):
    # index 3 is 'ru_isrss'
    return self.usageData[5]

# The following is an implmentation sent to me by David Stozzi from LLNL that he uses (given to him by Joseph M. Koning also at LLNL)
#   It uses a direct 'cat' of the linux kernal's /proc info. I put this in because I am courrious about how well Python built 'resource' pakckage
#   will match... Also it gives acces to more/different info.
import os
class MemUsageLinuxProcQuery:
  def __init__(self):
    self._proc_status = '/proc/%d/status' % os.getpid()
    self._scale = {'kB': 1024.0, 'mB': 1024.0*1024.0, 'KB': 1024.0, 'MB': 1024.0*1024.0}

  def _VmB(self, VmKey):
    # get pseudo file  /proc/<pid>/status
    try:
      t = open(self._proc_status)
      v = t.read()
      t.close()
    except:
      return 0.0  # non-Linux?
    
    # get VmKey line e.g. 'VmRSS:  9999  kB\n ...'
    i = v.index(VmKey)
    v = v[i:].split(None, 3)  # whitespace
    if len(v) < 3:
      return 0.0  # invalid format?
    # convert Vm value to bytes
    return float(v[1]) * self._scale[v[2]]

  def memory(self,since=0.0):
     '''Return memory usage in bytes.
     '''
     return _VmB('VmSize:') - since

  def resident(self,since=0.0):
     '''Return resident memory usage in bytes.
     '''
     return _VmB('VmRSS:') - since


  def stacksize(self,since=0.0):
     '''Return stack size in bytes.
     '''
     return _VmB('VmStk:') - since

  def hwm(self,since=0.0):
     '''Return high water mark size in bytes.
     '''
     return _VmB('VmHWM:') - since

  def peak(self,since=0.0):
     '''Return peak memory size in bytes.
     '''
     return _VmB('VmPeak:') - since

  def getResidentMemUsage(self):
    return resident()
  def getStackMemUsage(self):
    return stacksize()


"""
--------------------- Higher level niceities ------------------------------------------------
"""
import math

# global state..
memcounterTimer = None
pollingInterval = 100.0

baseIO = MemUsageResource()
baseIOLinux = MemUsageLinuxProcQuery()
units = ['bytes', 'kb', 'mb','gb','you messed up big time', 'you messed up big really time or this is the magical future. enjoy your exobyte computer and hourly shuttles to the moon where they have legalized all vice. speak fondly of me. i sit at UCLA at 3:30am. Do they still have HR departments in your time. God I hate HR. And Jews. I hate them too.']

def scaleMemorySize( mem ):
  powerIdx = 0
  try:
    if(mem > 0):
      powerIdx = int( math.log10(float(mem))) / 3
  except:
    print "log died on the incoming value of"
    print mem
    exit(-1)

  d = "%.3f %s" % ( float(mem)/math.pow(1000.0, powerIdx ), units[powerIdx] )
  return d


"""
--------------------- Integration in OSHUN Event System ------------------------------------------------
"""


previousResidentMem  = 0
# add this event type to the event queue.
def register(profileIncludeFilters=None,profileExcludeFilters=None):
  from ..OshunEventQueue import Event, add_event
  start_event = Event( "main_loop_resource_profiler", reportResourceUsage)
  add_event("main_loop_tick_start", start_event)

def reportResourceUsage( local_variables, simulation_state ):
  global previousResidentMem
  baseIO.resultsStale = True
  memUsed = baseIO.getResidentMemUsage()
  try:
    memUsedLinux = baseIO.getResidentMemUsage(),
  except:
    memUsedLinux = 0

  delta = memUsed - previousResidentMem
  print "Resident Memory Usage: %s (%s from Linux Kernal /proc) (delta %s from previous)" % (scaleMemorySize(memUsed),scaleMemorySize(memUsedLinux[0]), scaleMemorySize(delta) )
  previousResidentMem = memUsed

def timerCallback():
  baseIO.resultsStale = True
  reportResourceUsage(None,None)
  threading.Timer( pollingInterval, timerCallback)

def start():
  timerCallback()
