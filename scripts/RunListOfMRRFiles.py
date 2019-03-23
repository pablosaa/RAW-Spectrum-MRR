#! /usr/bin/python

import glob, os, sys

# Default values:
CMD = '../MRR4ADMI'
UZIP = 'bunzip2 -q'
ext = 'raw'
pathin = '.'
pathout = ''

if len(sys.argv)==4:
	pathin = str(sys.argv[1])
	ext = str(sys.argv[2])
	pathout = str(sys.argv[3])
elif len(sys.argv)==3:
	ext = str(sys.argv[2])
	pathin = str(sys.argv[1])
elif len(sys.argv)==2:
	ext = '7z'
	pathin = str(sys.argv[1])
else:
	print 'USAGE:	'
	print ' > ./RunListOfMRRFiles /Disk/path_to_MRRfiles raw /Path_to_output/data'
	print '	> ./RunListOfMRRFiles /Disk/path_to_MRRfiles raw.bz2'
	print ' > ./RunListOfMRRFiles /Disk/path_to_MRRfiles'

for infile in glob.glob(pathin+'/*.'+ext):
	os.system(UZIP+' '+infile)
	end = infile.find(ext)-1
	mrrfile = infile[0:end]
	#print CMD+' '+mrrfile+' '+pathout
	os.system(CMD+' '+mrrfile+' '+pathout)
