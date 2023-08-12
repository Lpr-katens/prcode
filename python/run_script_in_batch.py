import subprocess as sp

sp.run('python -W ignore find_companion.py',shell=True,check=True,)
sp.run('python -W ignore find_pair.py',shell=True,check=True)#
# sp.run('python -W ignore pair_img.py',shell=True,check=True)