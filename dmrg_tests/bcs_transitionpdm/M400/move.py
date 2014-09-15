import os
import shutil

files = os.listdir(".")

for f in files:
  if f.startswith("Rotation") and "state0" in f:
    shutil.copy(f, "../"+f.replace("state0", "state1"))
  if f.startswith("wave") and "0.tmp" in f:
    shutil.copy(f, "../"+f.replace("0.tmp", "1.tmp"))
shutil.copy("statefile.0.tmp", "../statefile.1.tmp")
#print files
