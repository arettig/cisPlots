import sys
import io
from os import listdir
import re



title = "CH3CH3_x"
ranges = [(1.0, 5.0, 0.01)]
reverse = True




inFile = io.open(sys.argv[1], "r")
template = inFile.read()
inFile.close()
outFile = io.open("%s.in" % title, "wb")


def stop(dist, minDist, maxDist, reverse):
    if(reverse):
        return dist < minDist
    else:
        return dist > maxDist

for rangeTuple in ranges:
    minDist = rangeTuple[0]
    maxDist = rangeTuple[1]
    distStep = rangeTuple[2]

    if(reverse):
        dist = maxDist
    else:
        dist = minDist

    while(not stop(dist, minDist, maxDist, reverse)):
        toWrite = template.replace("dist1", str(dist / 2.0))
        toWrite = toWrite.replace("dist2", str(dist / 2.0 + 0.3893))

        # don't read in guess if first job
        if((not reverse and dist == minDist) or (reverse and dist == maxDist)):
            toWrite = toWrite.replace("scf_guess read", "")

        # separate jobs unless last geometry
        if((not reverse and dist != maxDist) or (reverse and dist - minDist > 1e-8)):
            toWrite += "\n\n\n@@@@@@\n\n"

        outFile.write(toWrite)

        if(reverse):
            dist -= distStep
        else:
            dist += distStep

outFile.close()
