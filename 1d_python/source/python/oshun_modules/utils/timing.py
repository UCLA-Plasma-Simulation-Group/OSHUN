
# timing libraries
import time
import time as time_utils
import inspect


def timeit_basic(f):

    def timed(*args, **kw):
        ts = time.time()
        result = f(*args, **kw)
        te = time.time()

        print 'func:%r args:[%r, %r] took: %2.4f sec' % \
          (f.__name__, args, kw, te-ts)
        return result

    return timed

class timings_info:
	
	def __init__(self, fragment_name):
		self.name = fragment_name
		self._num_times_run = 0
		self._total_time = 0.0
		self._min_time = 1e12
	def report_timing(self, time_taken):
		self._num_times_run += 1
		self._total_time += time_taken
		if time_taken < self._min_time:
			self._min_time = time_taken
	def min_time(self):
		return self._min_time
	def average_time(self):
		return self._total_time / float(self._num_times_run)
	def total_time(self):
		return self._total_time
	
	
class timings:
	def __init__(self):
		self.fragments = {}
		self.order_list = []
	def add_timing(self, func, time_taken):
		class_of_func =  self.get_class_that_defined_method(func)
		if class_of_func == None:
			name = func.__name__
		else:
			name = "%s->%s" % (func.__name__, class_of_func.__name__)
		
		if name in self.fragments:
			timing_info = self.fragments[name]
		else:
			timing_info = timings_info(name)
			self.order_list.append( name)
			self.fragments[name] = timing_info
		timing_info.report_timing(time_taken)
	def add_timing_by_name( self, name, time_taken):
		if name in self.fragments:
			timing_info = self.fragments[name]
		else:
			timing_info = timings_info(name)
			self.order_list.append( name)
			self.fragments[name] = timing_info
		timing_info.report_timing(time_taken)
		
	
	def get_class_that_defined_method1(self, func):
		# see if it is a method call of a class..
		if func.__self__:
			classes = [func.__self__.__class__]
		else:
			#in this case, it is an unbound function...
			classes = [method.im_class]
		
		while classes:
			c = classes.pop()
			if func.__name__ in c.__dict__:
				return c
			else:
				classes = list(c.__bases__) + classes
		return None	
	def get_class_that_defined_method(self, meth):
		for cls in inspect.getmro(meth.im_class):
			if meth.__name__ in cls.__dict__: return cls
		return None
	
	def print_results(self):
		for name in self.order_list:
			entry = self.fragments[name]
			print "------------------\n\t\t%s\n\t\t\t: min: %e ave %e total %e " % (name, entry.min_time(), entry.average_time(), entry.total_time())
		print ""
	

TIMINGS = timings()

# function decorator version of the timing...
def timeit(f, id=None):
	
	def time_clousure(*args, **kw):
		_id = id
		if _id == None:
			stack_frames = inspect.getouterframes(inspect.currentframe())
			try:
				calling_line_index = stack_frames[1][5]
				lineNumber = stack_frames[1][2]
				calling_line_source_code = stack_frames[1][4][calling_line_index]
				_id = "line %d: %s" % ( lineNumber, calling_line_source_code.strip() )
			finally:
				del stack_frames
		
		start_time = time.time()
		result = f(*args, **kw)
		end_time   = time.time()
		
		#TIMINGS.add_timing(f, end_time-start_time)
		TIMINGS.add_timing_by_name(_id, end_time-start_time)
		
		return result
		
	return time_clousure
