import json
import os
import sys
import time
try:
    import boto3
    import cecontext
except ImportError:
    print ("AWS batch not supported.")
stime = 2
sloops = 30
def createEnv(batchC_a, queue_a, pricing_a, instancetypes_a, cename_a, tag_a, profile_a, verbose=False):
    # delete the queue if it exists
    if verbose:
        print("createEnv: deleting queue " + queue_a)
    deleteQueue(batchC_a, queue_a)
    # delete the ce if it exists
    if verbose:
        print("createEnv: delete ce " + cename_a)
    deleteCE(batchC_a, cename_a)
    # create the ce
    if verbose:
        print("createEnv: create ce " + cename_a)
    createCE(batchC_a, pricing_a, instancetypes_a, cename_a, tag_a, profile_a)
    # create the queue
    if verbose:
        print("createEnv: create queue " + queue_a)
    createQueue(batchC_a, queue_a, cename_a)

def deleteQueue(batchC, queue):
    # see if the queue exists
    response = batchC.describe_job_queues(jobQueues = [queue])
    if len(response['jobQueues']) > 0:
        # disable the queue
        batchC.update_job_queue(jobQueue = queue, state='DISABLED')
        # delete the queue
        # since update_job_queue is async and there are no "waiters" for batch
        # we'll need to loop and sleep
        msg = ""
        for ctr in range(sloops):
            status = True
            time.sleep(stime)
            try:
                batchC.delete_job_queue(jobQueue = queue)
            except Exception as e:
                status = False
                msg = str(e)
            if status:
                break
        if not status:
            print("deleteQueue: Error " + msg)
            sys.exit(2)



def deleteCE(batchC, cename):
    # see if the ce exists
    response = batchC.describe_compute_environments(computeEnvironments=[cename])
    if len(response['computeEnvironments']) > 0:
        # disable
        msg = ""
        for ctr in range(sloops):
            status = True
            time.sleep(stime)
            try:
                batchC.update_compute_environment(computeEnvironment = cename,
                                                  state = 'DISABLED')
            except Exception as e:
                status = False
                msg = str(e)
            if status:
                break
        if not status:
            print("deleteCE: Error " + msg)
            sys.exit(2)

        # delete the ce
        for ctr in range(sloops):
            status = True
            time.sleep(stime)
            try:
                batchC.delete_compute_environment(computeEnvironment = cename)
            except Exception as e:
                status = False
                msg = str(e)
            if status:
                break
        if not status:
            print("deleteCE: Error " + msg)
            sys.exit(2)

def createQueue(batchC, queue, cename):
    msg = ""
    for ctr in range(sloops):
        status = True
        time.sleep(stime)
        try:
            batchC.create_job_queue(jobQueueName = queue,
                                      priority = 10,
                                      state='ENABLED',
                                      computeEnvironmentOrder=[
                                        {'computeEnvironment': cename,
                                         'order': 1}
                                        ]
                                      )
        except Exception as e:
            status = False
            msg = str(e)
        if status:
            break
    if not status:
        print("createQueue: Error " + msg)
        sys.exit(2)
    # wait for queue to be ready
    for ctr in range(sloops):
        status = True
        time.sleep(stime)
        try:
            response = batchC.describe_job_queues(jobQueues = [queue])
        except Exception as e:
            status = False
            msg = str(e)
        if status:
            break
    if not status:
        print("createQueue - waiting for queue to complete: Error  " + msg)
        sys.exit(2)


def createCE(batchC, pricing_a, instancetypes_a, cename_a, tag_a, profile_a):
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
    stime = 5
    status = True
    msg = ""
    for ctr in range(10):
        time.sleep(stime)
        try:
            res = batchC.create_compute_environment(
                        computeEnvironmentName = cename_a,
                        type = cectx.ctype(),
                        state = cectx.cstate(),
                        computeResources = ce_resources,
                        serviceRole = ce_servicerole)
        except Exception as e:
            status = False
            msg = str(e)
        if status:
            break
    if not status:
        print("createCE: Error " + msg)
        sys.exit(2)
