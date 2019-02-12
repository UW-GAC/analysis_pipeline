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
    # delete the queue if it exists
    if verbose:
        print("createEnv: deleting queue " + queue_a)
    deleteQueue(batchC_a, queue_a, verbose)
    # delete the ce if it exists
    if verbose:
        print("createEnv: delete ce " + cename_a)
    deleteCE(batchC_a, cename_a, verbose)
    # create the ce
    if verbose:
        print("createEnv: create ce " + cename_a)
    createCE(batchC_a, pricing_a, instancetypes_a, cename_a, tag_a, profile_a, verbose)
    # create the queue
    if verbose:
        print("createEnv: create queue " + queue_a)
    createQueue(batchC_a, queue_a, cename_a, verbose)

def deleteQueue(batchC, queue, verbose):
    # see if the queue exists
    response = batchC.describe_job_queues(jobQueues = [queue])
    if len(response['jobQueues']) > 0:
        # disable the queue
        if verbose:
            print("deleteQueue: update_job_queue to Disabled ...")
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
                print("Error deleteQueue: queue did not become DISABLED before " + str(maxTime) + " seconds")
                sys.exit(2)
        # delete the queue
        if verbose:
            print("deleteQueue: delete_job_queue ...")
        batchC.delete_job_queue(jobQueue = queue)
        # wait for status to be deleted
        maxTime = 90
        timeW = 0
        sTime = 5
        while True:
            response = batchC.describe_job_queues(jobQueues = [queue])
            if len(response['jobQueues']) == 0:
                break
            if response['jobQueues'][0]['status'] == 'DELETED':
                break
            time.sleep(sTime)
            timeW += sTime
            if timeW > maxTime:
                print("Error deleteQueue: queue was not deleted before " + str(maxTime) + " seconds")
                sys.exit(2)
        if verbose:
            print("deleteQueue: delete_job_que completed")
    else:
        if verbose:
            print("deleteQueue: queue does not exist")


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

def createQueue(batchC, queue, cename, verbose):
    # create the queue
    if verbose:
        print("createQueue: create_job_queue ...")
    batchC.create_job_queue(jobQueueName = queue,
                            priority = 10,
                            state='ENABLED',
                            computeEnvironmentOrder=[
                                {'computeEnvironment': cename,
                                 'order': 1} ]
                            )
    # wait for for queue to be created
    maxTime = 60
    timeW = 0
    sTime = 2
    while True:
        response = batchC.describe_job_queues(jobQueues = [queue])
        if len(response['jobQueues']) > 0:
            break
        if response['jobQueues'][0]['status'] == 'VALID':
            break
        time.sleep(sTime)
        timeW += sTime
        if timeW > maxTime:
            print("Error createQueue: queue was not created before " + str(maxTime) + " seconds")
            sys.exit(2)
    if verbose:
        print("createQueue: queue created")

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
    batchC.create_compute_environment(
                computeEnvironmentName = cename_a,
                type = cectx.ctype(),
                state = cectx.cstate(),
                computeResources = ce_resources,
                serviceRole = ce_servicerole)
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
