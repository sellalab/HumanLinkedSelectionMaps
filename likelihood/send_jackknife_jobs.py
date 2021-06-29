__author__ = 'davidmurphy'


import re
import os
import sys
import glob
import time
import shutil
import subprocess
from classes.runstruct import ChromStruct, root_dir


run_path = '/ifs/data/c2b2/gs_lab/dam2214/run'
final_dir = root_dir + '/result/final_files'


reidx = re.compile('iprm_(\d\d)')


def jobcount():
    """count the current cluster jobs"""
    return int(subprocess.check_output('qstat|wc', shell=1).split()[0])


def errmsg(msg):
    sys.stderr.write(msg+'\n')
    sys.stdout.flush()

class JobSender:
    """a class that sends jackknife jobs and moniters job counts"""
    def __init__(self, f_init, istart, iend):
        """initialize the JobSender and define the range of jobs to send"""
        self.f_init = os.path.abspath(f_init)
        self.istart = int(istart)
        self.icur = int(istart)
        self.iend = int(iend)
        self.cst = ChromStruct('chr1', init=self.f_init)
        self.lbl = self.cst.tkn

    def _sendjob(self, idx):
        """send job for current index"""
        command = self.cmmd_string1.format(idx=idx) + ' > /dev/null'
        subprocess.call(command, shell=True)

    def _sendbatch(self):
        """send batch of jobs for the current jackknife index"""
        for idx in xrange(15):
            self._sendjob(idx)

    def _movefiles(self):
        """move jackknife step 1 files to a jkidx labeled folder"""
        # make new folder for jkidx
        fpath = self.jkfolder
        if not os.path.isdir(fpath):
            os.mkdir(fpath)
        # move all jkidx files to the new folder using glob patterns
        globpattern = '{}/*.{}.*'.format(final_dir, self.jkstring)
        jkfiles = glob.glob(globpattern)
        msg = '{} files moved for {}.'.format(len(jkfiles), self.jkstring)
        errmsg(msg)
        for f in jkfiles:
            shutil.move(f, fpath)

    def _sendmissing(self):
        """send missing step 1 jobs"""
        # check that all expected indices are in the set
        for i in range(15):
            if i not in self.step1indices:
                jobwaiting = True
                while jobwaiting:
                    njb = jobcount()
                    # if there are too many jobs, pause for 5 minutes
                    if njb > 900:
                        msg = '--There are {} jobs. Pausing 5min.'.format(njb)
                        errmsg(msg)
                        time.sleep(300)
                    else:
                        self._sendjob(i)
                        msg = '--Sent jkidx{}:idx{} job.'.format(self.icur, i)
                        errmsg(msg)
                        jobwaiting = False

    def process_step1_jobs(self):
        """process a set of jobs from jackknife range istart to iend"""
        while self.icur < self.iend:
            # check how many jobs are currently running
            njb = jobcount()
            # if there are too many jobs, pause for 5 minutes
            if njb > 900:
                msg = 'There are {} jobs. Pausing for 5min.'.format(njb)
                errmsg(msg)
                time.sleep(300)
            # otherwise, send the batch of jobs for the current jackknife
            else:
                self._sendbatch()
                msg = 'Job batch for jkidx {} sent.'.format(self.icur)
                errmsg(msg)
                self.icur += 1

    def process_step2_jobs(self):
        """process a set of step 2 jobs from jackknife range istart:iend"""
        # create JackKnifeJob object
        while self.icur < self.iend:
            # check how many jobs are currently running
            njb = jobcount()
            # if there are too many jobs, pause for 5 minutes
            if njb > 900:
                msg = 'There are {} jobs. Pausing for 5min.'.format(njb)
                errmsg(msg)
                time.sleep(300)
            # otherwise, send inf step 2 jobs if folder is complete
            else:
                # folder is complete and ready for step 2
                if self.sendstep2ok:
                    subprocess.call(self.cmmd_string2, shell=True)
                    msg = 'Inf step 2 for jkidx {} sent.'.format(self.icur)
                    errmsg(msg)
                # folder is empty, probably a genome gap
                elif len(self.step1files) == 0:
                    msg = 'jkidx {} folder empty. Skipping.'.format(self.icur)
                    errmsg(msg)
                # some jobs are missing - send missing
                else:
                    msg = 'Sending missing jobs for jkidx {}.'.format(self.icur)
                    errmsg(msg)
                    self._sendmissing()
                self.icur += 1

    def move_all_files(self):
        """move all files in the istart:iend range to their folders"""
        for jkidx in range(self.istart, self.iend+1):
            self.icur = jkidx
            self._movefiles()

    @property
    def cmmd_string1(self):
        """get command string to send inf step 1 job using qsub"""
        # program path
        prog = run_path + '/sh/jkinf1.sh'
        # log file name
        lfmt = '{}/{}.idx{{idx}}.jkidx{}.inf1.log'
        flog = lfmt.format(run_path, self.lbl, self.icur)
        # qsub command string
        sfmt = 'qsub -l mem=24G,time=10:: -cwd -j y -o {} {} {} {{idx}} {}'
        cmmd = sfmt.format(flog, prog, self.f_init, self.icur)

        return cmmd

    @property
    def cmmd_string2(self):
        """get command string to send inf step 2 job using qsub"""
        # program path
        prog = run_path + '/sh/jkinf2.sh'
        # log file name
        flog = run_path + '/{}_{}.jkinf2.log'.format(self.lbl, self.jkstring)
        # qsub command string
        sfmt = 'qsub -l mem=24G,time=10:: -cwd -j y -o {} {} {} > /dev/null'
        cmmd = sfmt.format(flog, prog, self.jkfolder.split('/')[-1])

        return cmmd

    @property
    def jkstring(self):
        """return jacknife index string based on current jkidx"""
        return 'jkidx_{:04}'.format(self.icur)

    @property
    def jkfolder(self):
        """get the folder path for a given jkidx"""
        return '{}/{}_{}'.format(final_dir, self.lbl, self.jkstring)

    @property
    def allfolderfiles(self):
        """list of all files in the current jkidx folder"""
        return os.listdir(self.jkfolder)

    @property
    def step1files(self):
        """return a list of the files in the completed step 1 folder"""
        fl = [f for f in self.allfolderfiles if f.endswith('.final.txt')]
        return fl

    @property
    def step1indices(self):
        """return the indices of all the files in step1files"""
        return [int(reidx.search(f).group(1)) for f in self.step1files]

    @property
    def sendstep2ok(self):
        """True the folder is ready for step 2"""
        # check that all 15 step 1 files are there
        check_1 = all(i in self.step1indices for i in range(15))
        # check final composite step 2 file is NOT there
        check_2 = not any('.composite.txt' in f for f in self.allfolderfiles)
        return check_1 and check_2


def main():
    if len(sys.argv) != 4:
        print 'usage: send_jackknife_jobs <f_init> <istart> <iend>'
        exit(1)

    # get command line args
    f_init, istart, iend = sys.argv[1:]

    # initialize JobSender and process batch of jobs that are given
    jbs = JobSender(f_init, istart, iend)
    # jbs.process_step1_jobs()
    jbs.move_all_files()
    jbs.icur = int(istart)
    jbs.process_step2_jobs()


if __name__ == '__main__':
    main()
