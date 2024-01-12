#!/usr/bin/env python  
#
# Copyright 2023 Shreyas Sunil Gaikwad

# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions
# are met:

# 1. Redistributions of source code must retain the above copyright
# notice, this list of conditions and the following disclaimer.

# 2. Redistributions in binary form must reproduce the above copyright
# notice, this list of conditions and the following disclaimer in the
# documentation and/or other materials provided with the distribution.

# 3. Neither the name of the copyright holder nor the names of its
# contributors may be used to endorse or promote products derived from
# this software without specific prior written permission.

# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
# FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
# COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
# INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
# BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
# LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
# ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.

# /**********************************************************************/
# /*! \file  
#
#  \brief Deal with metadata in some data files during the preprocessing
#         step for SICOPOLIS-AD v2.
# */
# /**********************************************************************/



import numpy as np
import sys
import argparse
import os
import re


### Deal with metadata for TSURFACE==4 case
def TSURFACE_4(sico_specs_file, sico_variables_file):

	'''
	sico_specs_file: Location of the header file
	sico_variables_file: Location of the sico_variables_m.F90 file
	'''

	try:
		with open(sico_specs_file) as f:
			file_lines = f.readlines()
			for line in file_lines:
				if ('#define TSURFACE' in line):
					
					l = []
					for t in line.split():
						try:
							l.append(int(t))
						except ValueError:
							pass
					TSURFACE = l[0]
					break
			for line in file_lines:
				if ('#define GRIP_TEMP_FILE' in line):
					
					l = []
					for t in line.split():
						try:
							l.append(t)
						except ValueError:
							pass
					Tsurface_data_file = '../sico_in/general/' + l[2][1:-1]
					break
	except FileNotFoundError :
		print(f'{sico_specs_file} does not exist')
		sys.exit(1)

	except Exception as err :
		print('Some error in checking sico specs for TSURFACE value')
		print(err)
		sys.exit(1)

	if TSURFACE==4:
		try:
			with open(Tsurface_data_file) as f:
	
				file_lines = f.readlines()
	
				if ('#' in file_lines[0]):
					l = []
					for t in file_lines[0].split():
						try:
							l.append(int(t))
						except ValueError:
							pass
					
					[grip_time_min, grip_time_stp, grip_time_max] = [entry for entry in l]
				else:
					"Metadata missing in T surface data file!"
					sys.exit(1)
	
		except FileNotFoundError :
			print(f'{Tsurface_data_file} does not exist')
			sys.exit(1)
	
		except Exception as err :
			print('Some error in using metadata from TSURFACE == 4 data file.')
			print(err)
			sys.exit(1)
	
		try:
			with open(sico_variables_file) as f:
	
				file_lines = f.readlines()
				new_file_lines = []
				for line in file_lines:
					if ('!@ python_automated_metadata GRIP_TIME_MIN @' in line):
						line = f'   integer(i4b) :: grip_time_min = {grip_time_min} !@ python_automated_metadata GRIP_TIME_MIN @ \n'
					if ('!@ python_automated_metadata GRIP_TIME_STP @' in line):
						line = f'   integer(i4b) :: grip_time_stp = {grip_time_stp} !@ python_automated_metadata GRIP_TIME_STP @ \n'
					if ('!@ python_automated_metadata GRIP_TIME_MAX @' in line):
						line = f'   integer(i4b) :: grip_time_max = {grip_time_max} !@ python_automated_metadata GRIP_TIME_MAX @ \n'
					if ('!@ python_automated_metadata NDATA_GRIP @' in line):
						ndata_grip = int((grip_time_max-grip_time_min)/grip_time_stp)
						line = f'   integer(i4b), parameter :: ndata_grip = {ndata_grip} !@ python_automated_metadata NDATA_GRIP @ \n'
	
					new_file_lines.append(line)
	
			with open(sico_variables_file, "w") as f:
				f.write(''.join(new_file_lines))
	
		except FileNotFoundError :
			print(f'{sico_variables_file} does not exist')
			sys.exit(1)
	
		except Exception as err :
			print('Some error in copying TSURFACE metadata into sico_variables.')
			print(err)
			sys.exit(1)

		return None

	else:
		return None

### Deal with metadata for SEA_LEVEL==3 case
def SEA_LEVEL_3(sico_specs_file, sico_variables_file):

	'''
	sico_specs_file: Location of the header file
	sico_variables_file: Location of the sico_variables_m.F90 file
	'''

	try:
		with open(sico_specs_file) as f:
			file_lines = f.readlines()
			for line in file_lines:
				if ('#define SEA_LEVEL' in line):
					
					l = []
					for t in line.split():
						try:
							l.append(int(t))
						except ValueError:
							pass
					SEA_LEVEL = l[0]
					break
			for line in file_lines:
				if ('#define SEA_LEVEL_FILE' in line):
					
					l = []
					for t in line.split():
						try:
							l.append(t)
						except ValueError:
							pass
					sea_level_data_file = '../sico_in/general/' + l[2][1:-1]
					break
	except FileNotFoundError :
		print(f'{sico_specs_file} does not exist')
		sys.exit(1)

	except Exception as err :
		print('Some error in checking sico specs for SEA_LEVEL value')
		print(err)
		sys.exit(1)

	if SEA_LEVEL==3:
		try:
			with open(sea_level_data_file) as f:
	
				file_lines = f.readlines()
	
				if ('#' in file_lines[0]):
					l = []
					for t in file_lines[0].split():
						try:
							l.append(int(t))
						except ValueError:
							pass
					
					[specmap_time_min, specmap_time_stp, specmap_time_max] = [entry for entry in l]
				else:
					"Metadata missing in sea level data file!"
					sys.exit(1)
	
		except FileNotFoundError :
			print(f'{sea_level_data_file} does not exist')
			sys.exit(1)
	
		except Exception as err :
			print('Some error in using metadata from SEA_LEVEL == 3 data file.')
			print(err)
			sys.exit(1)
	
		try:
			with open(sico_variables_file) as f:
	
				file_lines = f.readlines()
				new_file_lines = []
				for line in file_lines:
					if ('!@ python_automated_metadata SPECMAP_TIME_MIN @' in line):
						line = f'   integer(i4b) :: specmap_time_min = {specmap_time_min} !@ python_automated_metadata SPECMAP_TIME_MIN @ \n'
					if ('!@ python_automated_metadata SPECMAP_TIME_STP @' in line):
						line = f'   integer(i4b) :: specmap_time_stp = {specmap_time_stp} !@ python_automated_metadata SPECMAP_TIME_STP @ \n'
					if ('!@ python_automated_metadata SPECMAP_TIME_MAX @' in line):
						line = f'   integer(i4b) :: specmap_time_max = {specmap_time_max} !@ python_automated_metadata SPECMAP_TIME_MAX @ \n'
					if ('!@ python_automated_metadata NDATA_SPECMAP @' in line):
						ndata_specmap = int((specmap_time_max-specmap_time_min)/specmap_time_stp)
						line = f'   integer(i4b), parameter :: ndata_specmap = {ndata_specmap} !@ python_automated_metadata NDATA_SPECMAP @ \n'
	
					new_file_lines.append(line)
	
			with open(sico_variables_file, "w") as f:
				f.write(''.join(new_file_lines))
	
		except FileNotFoundError :
			print(f'{sico_variables_file} does not exist')
			sys.exit(1)
	
		except Exception as err :
			print('Some error in copying sea level metadata into sico_variables.')
			print(err)
			sys.exit(1)

		return None

	else:
		return None

if __name__ == '__main__':

	## Example Usage from src/ dir -
	## python subroutines/tapenade/deal_with_metadata.py -SS ../headers/sico_specs_${HEADER}.h -SV subroutines/general/sico_variables_m.F90

	parser = argparse.ArgumentParser()
	parser.add_argument("-SV", "--sico_variables", help="name of sico variables file", type=str, required=True)
	parser.add_argument("-SS", "--sico_specs", help="name of sico specs file", type=str, required=True)
	args = parser.parse_args()

	SEA_LEVEL_3(args.sico_specs, args.sico_variables)
	TSURFACE_4(args.sico_specs, args.sico_variables)
