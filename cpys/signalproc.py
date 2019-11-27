#
## MODULES 
#
from scipy import ndimage as nd # for gaussian filter
import sys
from my_lib import *
#
## FUNCTIONS
#
# RETURN VALUE OF NORMALIZED GAUSSIAN DISTRIBUTION (MEAN: mu, VARIANCE: sigma) IN POINT x
def gaussian(x, mu, sigma):
  return 1 / (sigma * sqrt(2 * pi)) * exp(-(x - mu)**2 / (2 * sigma**2))

# FILTER data WITH GAUSSIAN FILTER OF VARIANCE sigma
def filter_gaussian(data, sigma):
  return nd.filters.gaussian_filter1d(data, sigma) * sqrt(2 * pi) # last term necessary provided Python uses unitary FT

# COMPUTE GAUSSIAN VARIANCE sigma (IN SAMPLES UNITS, DIMENSIONLESS) FROM GAUSSIAN PARAMETER alpha (USED IN CPS)
def alpha2sigma(alpha, sampling_freq):
  sigma = float(sampling_freq) / ( float(alpha) * sqrt(2) ) 
  # multiplied by sampling freq so that it's dimensionless, for we operate usually in 'samples domain')
  return sigma

# FILTER RF WITH GAUSSIAN FILTER
def filter_RF(data, alpha, sampling_freq):
  sigma = alpha2sigma(alpha, sampling_freq)
  return (filter_gaussian(data, sigma), alpha) # ALPHA IN FILE NAMES (MEANINGFUL AND SHORTER)

# TURN LIST li OF SUBLISTS INTO A PLAIN LIST
def flatten_list(li):
  return list(itertools.chain.from_iterable(li)) # make flat (single) list out of list of sublists

# DIVIDE PLAIN LIST li INTO LIST OF SUBLISTS, EACH CONTAINING n ELEMENTS  
def get_chunks(li, n):
  for i in range(0, len(li), n):
    yield li[i : i + n]

# CREATE A SAC FILE: WRITE HEADER AND DATA (CHUNK -> LINE CONSISTED OF 5 COLUMNS)
def write_SAC(fname, header, chunks):
  f = open(fname, 'w')
  for line in header:
    f.write(line) # write SAC header into f2
  for chunk in chunks:
    line = ''
    for number in chunk:
      number = "{:15.7E}".format(number)
      line = line + str(number)
    f.write(line + '\n') 
  f.close()

# CALCULATE A MEAN OF THE TIME-WINDOW BEFORE THE 1st PEAK AND MOVE THE RF VERTICALLY BY THIS MEAN
def move_RF_vertically(data, samples_before_1st_peak):
  mean =  np.mean(data[ : samples_before_1st_peak])
  return ([d - mean for d in data], '') # empty pattern value (because float in the filename would be too messy)

# PROCESS (process IS AN EXTERNAL FUNCTION!!!) RF FILE DECIDING - SAC OR PLAIN FILE
def process_RF(file_name, pattern, process, *args): 
  content = []
  # READ FILE
  f = open(file_name, 'r')
  for line in f:
    content.append(line) # DONT SPLIT (WE WANT PRESERVE HEADER FORMAT UNTOUCHED)
  f.close()
  ## FOR SAC FILES
  if len(content[0].split(None)) > 1: 
    header = []
    data = []
    header = content[ : sac_header_length]
    data = content[sac_header_length : ]
    del content # DELETE THIS VARIABLE, IT WONT BE USFUL ANY MORE 
    # TRANSFORM DATA INTO A SINGLE LIST AND FILTER IT
    data = [line.split(None) for line in data]
    data = flatten_list(data)
    data = [float(i) for i in data]
    processed, pattern_value = process(data, *args) # process is a variable-function
    chunks = get_chunks(processed, nr_of_sac_columns) # 5 is a number of columns in sac data part
    # WRITE TO FILE
    stem = core(file_name, ' ', '_a.sac')
    fname = stem + pattern + pattern_value + '_a.sac'
    write_SAC(fname, header, chunks)
  ## FOR PLAIN 1-COLUMN FILES
  else: 
    data = [float(line) for line in content]
    processed, pattern_value = process(data, *args)
    # WRITE TO FILE
    stem = core(file_name, ' ', '.txt')
    f = open(stem + '_' + pattern + pattern_value + '.txt', 'w')
    for line in processed:
      f.write(str(line) + '\n')
    f.close()