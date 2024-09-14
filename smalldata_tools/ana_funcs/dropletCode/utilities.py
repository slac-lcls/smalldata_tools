import resource
import time


def printMsg(eventNr, run, rank=0, size=1):
    printFreq = 10
    # if eventNr > 10000:
    #  printFreq = 10000
    if eventNr > 1000:
        printFreq = 1000
    elif eventNr > 120:
        printFreq = 100

    if eventNr % printFreq == 0:
        if rank == 0:
            usage = resource.getrusage(resource.RUSAGE_SELF)
            print(
                "*** In Event: run",
                run,
                ",event# in single job =",
                eventNr,
                ", total about ",
                eventNr * size,
                " memory used: ",
                usage[2] * resource.getpagesize() / 1000000.0,
                " at ",
                time.strftime("%X"),
            )
