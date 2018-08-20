#!/usr/bin/python
# get-config.py -- get a python config variable.

import os
import distutils
import distutils.sysconfig
#from setuptools import setup
import sys
import datetime

if len(sys.argv)!=2:
	print("ERROR -- one argument required.")
else:
	# f = open('py-log.txt', 'a')
	# now=str(datetime.datetime.now())
	# f.write("Hello "+sys.argv[1]+" "+now+"\n")
	# f.close()
	print(distutils.sysconfig.get_config_var(sys.argv[1]))
