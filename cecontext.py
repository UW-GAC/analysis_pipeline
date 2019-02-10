import json
import os
import sys

class cecontext(object):
    def __init__(self,ctx_file=None, ctx_version="1.0", verbose = False):
        self.verbose = verbose
        self.ctx_file = os.path.join(os.path.dirname(os.path.abspath(sys.argv[0])),
                                     "cecontext.json")
        if self.verbose:
            print(">>>cecontext: ctx file is: " + self.ctx_file)
        if ctx_file != None:
            self.ctx_file = ctx_file
        # open the json ctx file
        with open(self.ctx_file) as cfh:
            ctxinfo = json.load(cfh)
        # check version
        key = "version"
        if key in ctxinfo:
            if ctxinfo[key] != ctx_version:
                print("Error: version of : " + self.ctx_file + " should be " + ctx_version +
                      " not " + ctxinfo[key])
                sys.exit(2)
        else:
            print("Error: " + key + " key not found in " + self.ctx_file)
            sys.exit(2)
        # type and state
        key = "type"
        if key in ctxinfo:
            self.type = ctxinfo[key]
        else:
            print("Error: " + key + " key not found in " + self.ctx_file)
            sys.exit(2)
        key = "state"
        if key in ctxinfo:
            self.state = ctxinfo[key]
        else:
            print("Error: " + key + " key not found in " + self.ctx_file)
            sys.exit(2)

        self.accnt_key = "accnt_ctx"
        # get the account resources
        if self.accnt_key in ctxinfo:
            self.accnt_ctx = ctxinfo[self.accnt_key]
        else:
            print("Error: " + key + " key not found in " + self.ctx_file)
            sys.exit(2)
        # get the common resources
        self.resource_key = "resources"
        if self.resource_key in ctxinfo:
            self.resources = ctxinfo[self.resource_key]
        else:
            print("Error: " + key + " key not found in " + self.ctx_file)
            sys.exit(2)
        #
        self.resource_names = self.resources.keys()
        self.accntctx_names = self.accnt_ctx.keys()
        if self.verbose:
            for name,resources in self.accnt_ctx.iteritems():
                print( "\t>>>account ctx: name: " + name + "  resources: " + str(resources))
            for name,value in self.resources.iteritems():
                print( "\t>>>resources: name: " + name + "  value: " + str(value))
    def cstate(self):
        return self.state
    def ctype(self):
        return self.type
    def accntnames(self):
        return self.accntctx_names
    def resourcenames(self):
        return self.resource_names
    def accntctx(self, accntname_a):
        ctx = None
        if accntname_a in self.accntctx_names:
            ctx = self.accnt_ctx[accntname_a]
        return ctx
    def commonresource(self, resname_a):
        rsrc = None
        if resname_a in self.resource_names:
            rsrc = self.resources[resname_a]
        return rsrc
    def allaccntresources(self, accntname_a):
        arsrc = None
        if accntname_a in self.accntctx_names:
            # get resources in account
            if "resources" in self.accnt_ctx[accntname_a].keys():
                arsrc = self.accnt_ctx[accntname_a]["resources"]
        return arsrc
    def accntresource(self, accntname_a,resname_a):
        rsrc = None
        arsrc = self.allaccntresources(accntname_a)
        if arsrc != None:
            if resname_a in arsrc.keys():
                rsrc = arsrc[resname_a]
        return rsrc
    # get all resources for an account
    def allceresources(self, accntname_a):
        arsc = None
        if accntname_a in self.accntctx_names:
            # combine the accnt resources with other resources
            arsc = dict(self.allaccntresources(accntname_a).items() + self.resources.items())
        return arsc
    def accntservice(self, accntname_a):
        svc = None
        if accntname_a in self.accntctx_names:
            if "serviceRole" in self.accnt_ctx[accntname_a].keys():
                svc = self.accnt_ctx[accntname_a]["serviceRole"]
        return svc
    def accntprofile(self, accntname_a):
        prof = None
        if accntname_a in self.accntctx_names:
            if "profile" in self.accnt_ctx[accntname_a].keys():
                prof = self.accnt_ctx[accntname_a]["profile"]
        return prof
