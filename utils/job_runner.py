########################################################################
# psychip - ChIP-seq pipeline for psychENCODE project
# Copyright (C) 2016  Weng's Lab @ UMASS medical school - Worcester (MA)
#
# Developers:
# Michael Purcaro : michael.purcaro@umassmed.edu 
# Arjan van der Velde : arjan.vandervelde@umassmed.edu
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
########################################################################

from __future__ import print_function

import os, sys, json, traceback
try:
    from joblib import Parallel, delayed
except:
    from files_and_paths import Dirs
    sys.path.append(Dirs.ToolsFnp('python2.7'))
    from joblib import Parallel, delayed
import multiprocessing

from utils import Utils

def runJob(job, idx, total, debug, func):
    try:
        rets = job(idx, total, debug, func)
        outs = ""
        for row in rets:
            for r in row:
                outs += r
        print(idx + 1, '/', total, "OK:", outs)
    except Exception, e:
        print(idx + 1, '/', total, "ERROR: {e}".format(e=e))
        print(traceback.print_exc())
        print(str(job)[:1024])

class PythonJob:
    def __init__(self, arr):
        self.arr = arr

    def writeJobFile(self, fnp):
        with open(fnp, 'w') as f:
            f.write(json.dumps(self.arr))

    def __call__(self, idx, total, debug, func):
        rets = []
        for args in self.arr:
            ret = func(idx, args)
            if ret:
                rets.append(ret)
                #print(idx, '/', total, "".join(ret).strip())
        return rets

    def __repr__(self):
        return "job: " + "\n\t" + str(self.arr)

    @staticmethod
    def ClusterBsubCmd(jobFolder, scriptFnp):
        return scriptFnp + " --job $LSB_JOBINDEX"

class BashJob:
    def __init__(self, arr):
        self.arr = arr

    def writeJobFile(self, fnp):
        with open(fnp, 'w') as f:
            f.write("#!/bin/bash" + "\n\n")
            for cmd in self.arr:
                f.write(' '.join(map(str, cmd)) + "\n")
        st = os.stat(fnp)
        os.chmod(fnp, st.st_mode | 0111) # chmod +x

    def __call__(self, idx, total, debug, func = None):
        rets = []
        for cmds in self.arr:
            ret = Utils.runCmds(cmds, debug)
            if ret:
                rets.append(ret)
                #print(idx, '/', total, "".join(ret).strip())
        return rets

    def __repr__(self):
        s = "bash job:\n\t"
        for c in self.arr:
            s += "\n\t" + str(c)
            try:
                s += "\n\t" + " ".join(c)
            except:
                pass
        return s

    @staticmethod
    def ClusterBsubCmd(jobFolder, scriptFnp):
        return os.path.join(jobFolder, "$LSB_JOBINDEX")

class JobRunner:
    def __init__(self, fraction = 0.5, cpus = None, scriptFnp = None, jobType = BashJob):
        if cpus:
            self.num_cores = cpus
        else:
            self.num_cores = int(float(multiprocessing.cpu_count()) * fraction)
        self.scriptFnp = scriptFnp
        self.jobType = jobType

        self.jobs = [] # array of array of cmds

    def append(self, job):
        self.jobs.append(self.jobType(job))

    def run(self, func = None):
        total = len(self.jobs)
        ret = Parallel(n_jobs = self.num_cores)(delayed(runJob)(job, idx, total, False, func)
                                                for idx, job in enumerate(self.jobs))
        self.jobs = []
        return ret

    def dump(self):
        for job in self.jobs:
            for cmd in job:
                print(' '.join(cmd) + "\n")

    def runOne(self, func = None):
        return runJob(self.jobs[0], 0, len(self.jobs), True, func)

    def runIdx(self, idx, func = None):
        return runJob(self.jobs[idx], idx, len(self.jobs), True, func)

    def dumpJobFiles(self, path):
        Utils.mkdir_p(path)
        print("writing", len(self.jobs), "job files to", path)
        for idx, job in enumerate(self.jobs):
            if (idx+1) % 1000 == 0:
                print("\twrote", idx + 1, "of", len(self.jobs))
            fnp = os.path.join(path, str(idx + 1))
            job.writeJobFile(fnp)
        print("dumped", len(self.jobs), "jobs to", path)

#    def cluster(self, folder, jobOptionsUsr = {}):
#        jobsFolder = os.path.join(folder, "jobs")
#        self.dumpJobFiles(jobsFolder)
#
#        Utils.mkdir_p(os.path.join(folder, "out"))
#        Utils.mkdir_p(os.path.join(folder, "err"))
#
#        jo = {"mem" : 4000, "time" : "1:00", "cores" : 1, "queue" : "short"}
#        for key in ["mem", "time", "cores", "queue"]:
#            if key in jobOptionsUsr:
#                jo[key] = jobOptionsUsr[key]
#
#        bins = []
#        maxJobs = 5000
#        for i in xrange(1, len(self.jobs), maxJobs):
#            start = i
#            end = i + maxJobs - 1 # inclusive
#            if end > len(self.jobs):
#                end = len(self.jobs)
#            bins.append([start, end])
#
#        outFnps = []
#        for b in bins:
#            outFnp = os.path.join(folder, "run_cluster_{}_{}.sh".format(b[0], b[1]))
#            outFnps.append(outFnp)
#            with open(outFnp, 'w') as f:
#                f.write(ClusterUtils.batchCmd(b[0], b[1],
#                                              os.path.join(folder, "err"),
#                                              os.path.join(folder, "out"),
#                                              self.jobType.ClusterBsubCmd(jobsFolder, self.scriptFnp),
#                                              jo))
#            print("wrote", outFnp)
#
#        if ClusterUtils.onCluster():
#            if GetYesNoToQuestion.immediate("run?"):
#                for outFnp in outFnps:
#                    cmds = ["bsub", "<", outFnp]
#                    print(Utils.runCmds(cmds)[0])
#        else:
#            ClusterUtils.rsyncToCluster(folder)
#            ClusterUtils.mkdirOutErrFolders(folder)
#            print("to run jobs, run this command on ghpcc:")
#            for outFnp in outFnps:
#                print("bsub <", outFnp)

class SequentialRunner:
    
    def __init__(self, sequence_length, cpus):
        self._runners = []
        for i in range(0, sequence_length):
            self._runners.append(JobRunner(cpus=cpus))

    def append_sublist(self, joblist, start_idx):
        assert len(joblist) + start_idx <= len(self._runners)
        for i in range(start_idx, start_idx + len(joblist)):
            self._runners[i].append([joblist[i - start_idx]])

    def append(self, joblist):
        assert len(joblist) == len(self._runners)
        for i in range(0, len(joblist)):
            self._runners[i].append([joblist[i]])

    def run(self):
        for i in range(0, len(self._runners)):
            self._runners[i].run()

    def clear(self):
        for i in range(0, len(self._runners)):
            self._runners[i] = JobRunner(cpus=cpus)
