#! /usr/bin/env python3
# process the submit command from Cluster in one of the following ways
#   1. just execute the command directly by calling popen
#   2. execute a pre-submit python script to support resume
#
import     sys
import     os.path
import     os
import     subprocess

def popen(cmd, sout=subprocess.PIPE, serr=subprocess.PIPE, pshell=True):
    #
    try:
        process = subprocess.Popen(cmd, stdout=sout, stderr=serr, shell=pshell)
    except Exception as e:
        emsg = "Error: Popen exception: " + str(e) + "\ncommand: " + cmd
        print(emsg)
        raise Exception(emsg)

    status = process.wait()
    if status != 0:
        print("Popen: Error executing cmd \n\t" + cmd)
        print("Popen: Error status " + str(status))
        if serr == subprocess.PIPE:
            pipe = process.stderr
            eMsg = pipe.readline()
            # compatibility p2/p3: byte seq or string converts to string
            eMsg = bytes(eMsg).decode()
            print("Popen: Error msg: " + eMsg)
        sys.exit(status)
    # if pipes, read the results of popen
    if sout == subprocess.PIPE:
        pipe = process.stdout
        popen_results = pipe.readline()
        # compatibility p2/p3: byte seq or string converts to string
        popen_results = bytes(popen_results).decode()
        popen_results = popen_results.strip(' \t\n\r')
    # else results have been sent to stdout
    else:
        popen_results = ""
    return popen_results

def popen_stdout(cmd, logfile=None, jobname=None, shell=True):
    sout = sys.stdout
    serr = sys.stderr
    # redirect stdout/stderr to file
    if logfile != None:
        flog = open ( logfile, 'w' )
        sout = flog
        serr = flog
    try:
        popen(cmd, sout=sout, serr=serr, pshell=shell)
    except Exception as e:
        # redirect stdout back to print message to stdout
        print("Popen: Exception executing cmd \n\t" + str(cmd))
        print("Popen: Exception is: " + str(e))
        sys.exit(2)
    if logfile != None:
        flog.close()
