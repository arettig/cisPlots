import sys
import pathlib2
import io
from os import listdir
import re



title = "H2"
ranges = [(0.5, 5.0, 0.01)]
reverse = True

spinFlipParams = ["", "spin_flip 1", "spin_flip_xcis true"]
spinParams = ["1", "3"]


inFile = io.open(sys.argv[1], "r")
template = inFile.read()
inFile.close()



def stop(dist, minDist, maxDist, reverse):
    if(reverse):
        return dist < minDist
    else:
        return dist > maxDist

for sfParam in spinFlipParams:
    for spin in spinParams:
        # create outfile
        outTitle = title
        if(sfParam == "spin_flip 1"):
            outTitle += "_sf"
        elif(sfParam == "spin_flip_xcis true"):
            outTitle += "_x"

        if(spin == "3"):
            outTitle += "_t"
        outFile = io.open("%s.in" % outTitle, "wb")
        
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
                toWrite = toWrite.replace("spinParam", spin).replace("sfParam", sfParam)
                if(sfParam == "spin_flip_xcis true"):
                    toWrite = toWrite.replace("unrestricted true", "unrestricted false")


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
