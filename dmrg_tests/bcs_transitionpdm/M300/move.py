import os
import shutil

files = os.listdir(".")

for f in files:
  if f.startswith("Rotation") and "state0" in f:
    shutil.copy(f, "../"+f.replace("state0", "state0"))
  if f.startswith("wave") and "0.tmp" in f:
    shutil.copy(f, "../"+f.replace("0.tmp", "0.tmp"))
shutil.copy("statefile.0.tmp", "../statefile.0.tmp")
#print files
