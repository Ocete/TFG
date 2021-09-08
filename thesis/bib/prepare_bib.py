#!/usr/bin/env python3

filename = 'bib/library.bib'

# Read in the file
with open(filename, 'r') as file:
	filedata = file.read()

# Replace the target string
filedata = filedata.replace('{\_}', '\_')
filedata = filedata.replace('{\#}', '\#')
filedata = filedata.replace('{\%}', '\%')
filedata = filedata.replace('{~}', '~')
filedata = filedata.replace('{\^}', '^')
filedata = filedata.replace('{\"}', '"')

# Add not-found citation in Mendeley manually
if 'Manin1980' not in filedata:
  filedata += """
  @article{Manin1980,
    title={Computable and uncomputable},
    author={Manin, Yuri},
    journal={Sovetskoye Radio, Moscow},
    volume={128},
    year={1980}
  }
  """

# Write the file out again
with open(filename, 'w') as file:
	file.write(filedata)