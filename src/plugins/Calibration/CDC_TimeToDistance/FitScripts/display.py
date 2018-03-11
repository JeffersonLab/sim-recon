import os
import subprocess

images = subprocess.check_output(["ls", "Monitoring/"]).splitlines()

for image in images:
    print("Displaying " + image + ". Close the window to move to the next image")
    subprocess.call(["display", "Monitoring/" + image])
    raw_input("Hit Enter")
