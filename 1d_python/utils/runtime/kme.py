import glob
import re
from subprocess import call

jls = glob.glob('*.joblog')
pattern = r'Job (.*), queue'

f = open( jls[-1] )
lines = f.readlines()

job_id = None
for line in lines:
	
	stuff = re.findall( pattern , line)
	if len(stuff) > 0:
		job_id = stuff[0]
	pass

if job_id != None:
	print "killing job with id %s" % job_id
	call(["qdel", job_id] )
	print "done"
	
else:
	print "Job id not found... doing nothing..."
