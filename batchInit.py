import json
import os
import sys
import time
try:
    import boto3
    import cecontext
except ImportError:
    print ("AWS batch not supported.")


def createEnv(batchC_a, queue_a, pricing_a, instancetypes_a, cename_a, tag_a, profile_a, verbose=False):
    # In order to place a new tag in the ce the following steps are done:
    #    1. If the queue exist (and has two ce's - ag_suport and cename),
    #          a. queue is disabled
    #          b. the cename is removed from queue
    #    2. cename is deleted
    #    3. a new cename is created (with a new tag)
    #    4. if the queue exists,
    #          a. cename is added to the queue
    #          b. queue is enabled
    #       else,
    #          a. the queue is created with two ce's: ag_support and cename
    #
    supportCE = "ag_support"
    # remove the ce from queue (if queue exists)
    if verbose:
        print("createEnv: deleting queue " + queue_a)
    cenames = removeCE(batchC_a, queue_a, cename_a, supportCE, verbose)
    # delete the ce if it exists
    if verbose:
        print("createEnv: delete ce " + cename_a)
    deleteCE(batchC_a, cename_a, verbose)
    # create the ce
    if verbose:
        print("createEnv: create ce " + cename_a)
    cenameArn = createCE(batchC_a, pricing_a, instancetypes_a, cename_a, tag_a, profile_a, verbose)
    cenames[cename_a] = cenameArn
    # update (or create) the queue
    if verbose:
        print("createEnv: create queue " + queue_a)
    updateQueue(batchC_a, queue_a, cenames, cename_a, supportCE, verbose)

def removeCE(batchC, queue, cename, supportCE, verbose):
    # get the cenames (supportCE and cename)
    ceinfo = batchC.describe_compute_environments(computeEnvironments = [cename, supportCE])
    ces = ceinfo['computeEnvironments']   # list of dicts
    cenames = {}
    for d in ces:
        cenames[d['computeEnvironmentName']]=d['computeEnvironmentArn']
    if len(cenames) == 0 or supportCE not in cenames.keys():
        print("removeCE: support ce " + supportCE + " not found.")
        sys.exit(2)
    # see if the queue exists
    queueD = batchC.describe_job_queues(jobQueues = [queue])
    if len(queueD['jobQueues']) > 0:
        if verbose:
            print("removeCE: update_job_queue to Disabled ...")
        # check for cename (it must exist)
        if cename not in cenames:
            print("removeCE: required ce " + cename + " not found.")
            sys.exit(2)
        # disable the queue
        batchC.update_job_queue(jobQueue = queue, state='DISABLED')
        # wait for state to change to disabled
        maxTime = 60
        timeW = 0
        sTime = 2
        while True:
            response = batchC.describe_job_queues(jobQueues = [queue])
            if response['jobQueues'][0]['status'] == 'VALID':
                break
            time.sleep(sTime)
            timeW += sTime
            if timeW > maxTime:
                print("Error removeCE: queue did not become DISABLED before " + str(maxTime) + " seconds")
                sys.exit(2)
        # queue should have two ce - supportCE (always there) and
        # cename (which is removed so tag can be updated)
        if verbose:
            print("removeCE: remove ce from queue ...")
        theQueue = queueD['jobQueues'][0]
        ceo = theQueue['computeEnvironmentOrder']
        if len(ceo) != 2:
            print("Error removeCE: queue " + queue + " has incorrect number of ce's")
            sys.exit(2)
        # insure cename is there
        cel = [d['computeEnvironment'] for d in ceo]
        if verbose:
            print("removeCE: ce list: " + str(cel))
        if cenames[cename] not in cel:
            print("Error removeCE: queue " + queue + " does not contain ce " + cename)
            sys.exit(2)
        # remove the cenname
        ceo[:] = [d for d in ceo if d.get('computeEnvironment') != cenames[cename]]
        # update the queue
        batchC.update_job_queue(jobQueue = queue, computeEnvironmentOrder = ceo)
        # wait for status to be valid
        maxTime = 60
        timeW = 0
        sTime = 2
        while True:
            response = batchC.describe_job_queues(jobQueues = [queue])
            if response['jobQueues'][0]['status'] == 'VALID':
                break
            time.sleep(sTime)
            timeW += sTime
            if timeW > maxTime:
                print("Error removeCE: queue was not updated before " + str(maxTime) + " seconds")
                sys.exit(2)
        if verbose:
            print("removeCE: ce removed from queue")
    else:
        if verbose:
            print("removeCE: queue does not exist")
    return cenames

def deleteCE(batchC, cename, verbose):
    # see if the ce exists
    response = batchC.describe_compute_environments(computeEnvironments=[cename])
    if len(response['computeEnvironments']) > 0:
        # disable
        if verbose:
            print("deleteCE: update_compute_environment to Disabled ...")
        batchC.update_compute_environment(computeEnvironment = cename, state = 'DISABLED')
        # wait for state to change to disabled
        maxTime = 60
        timeW = 0
        sTime = 2
        while True:
            response = batchC.describe_compute_environments(computeEnvironments=[cename])
            if response['computeEnvironments'][0]['status'] == 'VALID':
                break
            time.sleep(sTime)
            timeW += sTime
            if timeW > maxTime:
                print("Error deleteCE: ce did not become DISABLED before " + str(maxTime) + " seconds")
                sys.exit(2)
        # delete the ce
        batchC.delete_compute_environment(computeEnvironment = cename)
        # wait for status to be deleted
        maxTime = 60
        timeW = 0
        sTime = 2
        while True:
            response = batchC.describe_compute_environments(computeEnvironments = [cename])
            if len(response['computeEnvironments']) == 0:
                break
            if response['computeEnvironments'][0]['status'] == 'DELETED':
                break
            time.sleep(sTime)
            timeW += sTime
            if timeW > maxTime:
                print("Error deleteCE: ce was not deleted before " + str(maxTime) + " seconds")
                sys.exit(2)
        if verbose:
            print("deleteCE: ce deleted")
    else:
        if verbose:
            print("deleteCE: ce does not exist")

def updateQueue(batchC, queue, cenames, cename, supportCE, verbose):
    # see if the queue exists
    queueD = batchC.describe_job_queues(jobQueues = [queue])
    if len(queueD['jobQueues']) > 0:
        # update the queue
        if verbose:
            print("updateQueue: update_job_queue ...")
        theQueue = queueD['jobQueues'][0]
        ceo = theQueue['computeEnvironmentOrder']
        # append the ce to the ceo
        thece = {'order': 1, 'computeEnvironment': cenames[cename]}
        ceo.append(thece)   # note: just a reference which is ok in this case
        # make the call
        batchC.update_job_queue(jobQueue = queue, computeEnvironmentOrder = ceo,
                                state = 'ENABLED')
    else:
        # create the queue
        if verbose:
            print("updateQueue: create_job_queue ...")
        batchC.create_job_queue(jobQueueName = queue,
                                priority = 10,
                                state='ENABLED',
                                computeEnvironmentOrder=[
                                    {'computeEnvironment': cenames[cename],
                                     'order': 1},
                                    {'computeEnvironment': cenames[supportCE],
                                     'order': 2}  ]
                                )
    # wait for for queue to be created or updated and is in a valid state
    maxTime = 60
    timeW = 0
    sTime = 2
    while True:
        response = batchC.describe_job_queues(jobQueues = [queue])
        if (len(response['jobQueues']) > 0 and response['jobQueues'][0]['status']) == 'VALID':
            break
        time.sleep(sTime)
        timeW += sTime
        if timeW > maxTime:
            print("Error updateQueue: queue was not created or updated before " + str(maxTime) + " seconds")
            sys.exit(2)
    if verbose:
        print("updateQueue: queue created or updated")

def createCE(batchC, pricing_a, instancetypes_a, cename_a, tag_a, profile_a, verbose):
    # get the ce attributes from json file
    cectx = cecontext.cecontext()
    # get all the ce resources based on accnt ctxt
    ce_resources = cectx.allceresources(profile_a)
    if ce_resources == None:
        print("Invalid aws profile " + profile_a + " [" + str(cectx.accntnames()) + "]")
        sys.exit(2)
    ce_resources['instanceTypes'] = instancetypes_a.replace(' ','').split(',')
    ce_resources['type'] = pricing_a
    # service role
    ce_servicerole = cectx.accntservice(profile_a)
    if ce_servicerole == None:
        pError("Cannot find the serviceRole for " + profile_a)
        sys.exit(2)
    # create tag dictionary
    tagdict = {}
    tagdict['name'] = 'noname'
    tagdict['mode'] = 'autogen'
    tagdict['analysis'] = tag_a
    ce_resources['tags'] = tagdict
    # create the ce
    if verbose:
        print("createCE: create_compute_environment ...")
    ceinfo = batchC.create_compute_environment(
                computeEnvironmentName = cename_a,
                type = cectx.ctype(),
                state = cectx.cstate(),
                computeResources = ce_resources,
                serviceRole = ce_servicerole)
    ceArn = ceinfo['computeEnvironmentArn']
    # wait for the ce to be created
    maxTime = 90
    timeW = 0
    sTime = 2
    while True:
        response = batchC.describe_compute_environments(computeEnvironments = [cename_a])
        if (len(response['computeEnvironments']) > 0 and
            response['computeEnvironments'][0]['status'] == 'VALID'):
            break
        time.sleep(sTime)
        timeW += sTime
        if timeW > maxTime:
            print("Error createCE: ce was not created before " + str(maxTime) + " seconds")
            sys.exit(2)
    if verbose:
        print("createCE: ce created")
    return ceArn
