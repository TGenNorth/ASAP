'''
Created on Jul 8, 2015

@author: Darrin Lemmer
'''

import logging

job_manager = "SLURM"
job_manager_args = ""

def _submit_job(job_submitter, command, job_parms, waitfor_id=None, hold=False, notify=False):
    import subprocess
    import re
    import os

    # TODO(jtravis): remove unused output variable
    output = jobid = None
    logging.info("command = %s" % command)
    args = job_parms["args"] or job_manager_args or ""
    if job_submitter == "PBS":
        waitfor = ""
        if waitfor_id and waitfor_id[0]:
            dependency_string = waitfor_id[1] if len(waitfor_id) > 1 else 'afterok'
            waitfor = "-W depend=%s:%s" % (dependency_string, waitfor_id[0])
        queue = ""
        if job_parms["queue"]:
            queue = "-q %s" % job_parms["queue"]
        if hold:
            args += " -h"
        if notify:
            args += " -m e"
        submit_command = "qsub -V -d \'%s\' -w \'%s\' -l ncpus=%s,mem=%sgb,walltime=%s:00:00 -m a -N \'%s\' %s %s %s" % (
            job_parms["work_dir"], job_parms["work_dir"], job_parms['num_cpus'], job_parms['mem_requested'],
            job_parms['walltime'], job_parms['name'], waitfor, queue, args)
        logging.debug("submit_command = %s", submit_command)
        output = subprocess.getoutput("echo \"%s\" | %s - " % (command, submit_command))
        logging.debug("output = %s" % output)
        job_match = re.search('^(\d+)\..*$', output)
        if job_match:
            jobid = job_match.group(1)
        else:
            logging.warning("Job not submitted!!")
            print("WARNING: Job not submitted: %s" % output)
    elif job_submitter == "SLURM":
        waitfor = ""
        if waitfor_id:
            dependency_string = waitfor_id[1] if len(waitfor_id) > 1 else 'afterok'
            waitfor = "-d %s:%s" % (dependency_string, waitfor_id[0])
        queue = ""
        if job_parms["queue"]:
            queue = "-p %s" % job_parms["queue"]
        if hold:
            args += " -H"
        if notify:
            args += " --mail-type=END"
        #submit_command = "sbatch -D \'%s\' -c%s --mem=%s000 --time=%s:00:00 --mail-type=FAIL -J \'%s\' %s %s %s" % (
        #    job_parms["work_dir"], job_parms['num_cpus'], job_parms['mem_requested'], job_parms['walltime'],
        #    job_parms['name'], waitfor, queue, args)
        submit_command = "sbatch -D \'%s\' -c%s --mem=%s000 --mail-type=FAIL -J \'%s\' %s %s %s" % (
            job_parms["work_dir"], job_parms['num_cpus'], job_parms['mem_requested'],
            job_parms['name'], waitfor, queue, args)
        logging.debug("submit_command = %s" % submit_command)
        output = subprocess.getoutput("%s --wrap=\"%s\"" % (submit_command, command))
        logging.debug("output = %s" % output)
        job_match = re.search('^Submitted batch job (\d+)$', output)
        if job_match:
            jobid = job_match.group(1)
        else:
            logging.warning("Job not submitted!!")
            print("WARNING: Job not submitted: %s" % output)
    elif job_submitter == "SGE":
        waitfor = ""
        if waitfor_id:
            waitfor = "-hold_jid %s" % (re.sub(":", ",", waitfor_id[0]))
        queue = ""
        if job_parms["queue"]:
            queue = "-q %s" % job_parms["queue"]
        if hold:
            args += " -h"
        if notify:
            args += " -m e"
        mem_needed = float(job_parms['mem_requested']) * 1024 * 1024
        # Apparently the number of processors a job uses is controlled by the queue it is running on in SGE, so there is no way to request a specific number of CPUs??
        submit_command = "qsub -V -cwd \'%s\' -wd \'%s\' -l h_data=%sgb,h_rt=%s:00:00 -m a -N \'%s\' %s %s %s" % (
            job_parms["work_dir"], job_parms["work_dir"], mem_needed, job_parms['walltime'], job_parms['name'], waitfor,
            queue, args)
        logging.debug("submit_command = %s", submit_command)
        output = subprocess.getoutput("echo \"%s\" | %s - " % (command, submit_command))
        logging.debug("output = %s" % output)
        job_match = re.search('^(\d+)\..*$', output)
        if job_match:
            jobid = job_match.group(1)
        else:
            logging.warning("Job not submitted!!")
            print("WARNING: Job not submitted: %s" % output)
    else:
        command = re.sub('\n', '; ', command)
        work_dir = job_parms['work_dir']
        dependency_check = ""
        if waitfor_id:
            if re.search(":", waitfor_id[0]):
                pid_filename = os.path.join(work_dir, "%s_dependent_pids" % job_parms['name'])
                pid_file = open(pid_filename, 'w')
                pids = waitfor_id[0].split(":")
                pid_file.write("\n".join(pids))
                pid_file.close()
                dependency_check = "while [ -s %s ]; do sleep 600; for pid in `cat %s`; do kill -0 \"$pid\" 2>/dev/null || sed -i \"/^$pid$/d\" %s; done; done; rm %s; " % (
                    pid_filename, pid_filename, pid_filename, pid_filename)
            else:
                pid = waitfor_id[0]
                dependency_check = "while kill -0 %s; do sleep 300; done; " % pid
        # cpu_check = "grep 'cpu ' /proc/stat | awk '{usage=($2+$4)/($2+$4+$5)} END {print 1-usage}'"
        total_mem = subprocess.getoutput("free -g | grep Mem: | awk '{ print $2 }'")
        mem_requested = job_parms['mem_requested']
        if mem_requested > total_mem:
            mem_requested = total_mem
        print("total_mem = %s" % total_mem)
        mem_needed = float(mem_requested) * 750
        mem_check = "while [ `free -m | grep cache: | awk '{ print $4 }'` -lt %d ]; do sleep 300; done; " % mem_needed
        submit_command = "%s%s%s" % (dependency_check, mem_check, command)
        logging.debug("submit_command = %s" % submit_command)
        output_log = os.open(os.path.join(work_dir, "%s.out" % job_parms['name']), os.O_WRONLY | os.O_CREAT)
        proc = subprocess.Popen(submit_command, stderr=subprocess.STDOUT, stdout=output_log, shell=True, cwd=work_dir)
        jobid = str(proc.pid)
    logging.info("jobid = %s" % jobid)
    return jobid

def _run_bwa(sample, reads, reference, outdir='', dependency=None, sampath='samtools', bwapath='bwa', ncpus=4, args=''):
    import os
    read1 = reads[0]
    read2 = reads[1] if len(reads) > 1 else ""
    bam_string = "\'@RG\\tID:%s\\tSM:%s\'" % (sample, sample)
    job_params = {'queue':'', 'mem_requested':4, 'num_cpus':ncpus, 'walltime':8, 'args':''}
    job_params['name'] = "asap_bwa_%s" % sample
    aligner_name = "bwamem"
    aligner_command = "%s mem -R %s %s -t %s %s %s %s" % (bwapath, bam_string, args, ncpus, reference, read1, read2)
    bam_nickname = "%s-%s" % (sample, aligner_name)
    samview_command = "%s view -S -b -h -" % sampath
    #samsort_command = "%s sort - %s" % (sampath, bam_nickname)
    #samtools 1.3 version
    samsort_command = "%s sort -T %s -o %s.bam -" % (sampath, bam_nickname, bam_nickname)
    samindex_command = "%s index %s.bam" % (sampath, bam_nickname)
    command = "%s | %s | %s \n %s" % (aligner_command, samview_command, samsort_command, samindex_command)
    work_dir = os.path.join(outdir, aligner_name)
    job_params['work_dir'] = work_dir
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    final_file = os.path.join(work_dir, "%s.bam" % bam_nickname)
    job_id = _submit_job(job_manager, command, job_params, (dependency,))
    return (final_file, job_id)

def _run_bowtie2(sample, reads, reference, outdir='', dependency=None, sampath='samtools', bt2path='bowtie2', ncpus=4, args=''):
    import os
    read1 = reads[0]
    read2 = reads[1] if len(reads) > 1 else ""
    read_string = "-1 %s -2 %s" % (read1, read2) if read2 else "-U %s" % read1
    ref_string = os.path.splitext(reference)[0]
    bam_string = "--rg-id \'%s\' --rg \'SM:%s\'" % (sample, sample)
    job_params = {'queue':'', 'mem_requested':4, 'num_cpus':ncpus, 'walltime':8, 'args':''}
    job_params['name'] = "asap_bt2_%s" % sample
    aligner_name = "bowtie2"
    aligner_command = "%s %s -p %s %s -x %s %s" % (bt2path, args, ncpus, bam_string, ref_string, read_string)
    bam_nickname = "%s-%s" % (sample, aligner_name)
    samview_command = "%s view -S -b -h -" % sampath
    #samsort_command = "%s sort - %s" % (sampath, bam_nickname)
    #samtools 1.3 version
    samsort_command = "%s sort -T %s -o %s.bam -" % (sampath, bam_nickname, bam_nickname)
    samindex_command = "%s index %s.bam" % (sampath, bam_nickname)
    command = "%s | %s | %s \n %s" % (aligner_command, samview_command, samsort_command, samindex_command)
    work_dir = os.path.join(outdir, aligner_name)
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    final_file = os.path.join(work_dir, "%s.bam" % bam_nickname)
    job_params['work_dir'] = work_dir
    job_id = _submit_job(job_manager, command, job_params, (dependency,))
    return (final_file, job_id)

def _run_novoalign(sample, reads, reference, outdir='', dependency=None, sampath='samtools', novopath='novoalign', ncpus=4, args=''):
    import os
    read1 = reads[0]
    read2 = reads[1] if len(reads) > 1 else ""
    #paired_string = "-i PE 500,100" if read2 else ""
    paired_string = ""
    bam_string = "\'@RG\\tID:%s\\tSM:%s\'" % (sample, sample)
    job_params = {'queue':'', 'mem_requested':4, 'num_cpus':ncpus, 'walltime':24, 'args':''}
    job_params['name'] = "asap_novo_%s" % sample
    aligner_name = "novo"
    aligner_command = "%s -f %s %s %s -c %s -o SAM %s -d %s.idx %s" % (novopath, read1, read2, paired_string, ncpus, bam_string, reference, args)
    bam_nickname = "%s-%s" % (sample, aligner_name)
    samview_command = "%s view -S -b -h -" % sampath
    #samsort_command = "%s sort - %s" % (sampath, bam_nickname)
    #samtools 1.3 version
    samsort_command = "%s sort -T %s -o %s.bam -" % (sampath, bam_nickname, bam_nickname)
    samindex_command = "%s index %s.bam" % (sampath, bam_nickname)
    command = "%s | %s | %s \n %s" % (aligner_command, samview_command, samsort_command, samindex_command)
    work_dir = os.path.join(outdir, aligner_name)
    job_params['work_dir'] = work_dir
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    final_file = os.path.join(work_dir, "%s.bam" % bam_nickname)
    job_id = _submit_job(job_manager, command, job_params, (dependency,))
    return (final_file, job_id)

def findReads(path):
    import os
    import re
    import subprocess
    from collections import namedtuple
    read_list = []
    Read = namedtuple('Read', ['sample', 'reads'])
    for file in os.listdir(path):
        is_read = re.search('(.*)(\.f(?:ast)?q(\.gz)?)$', file, re.IGNORECASE)
        if is_read:
            sample_name = is_read.group(1)
            #print("Found read: %s" % sample_name)
            full_file = os.path.join(path, file)
            if os.path.getsize(full_file) == 0 or (is_read.group(3) and subprocess.getoutput("gzip -l %s | awk 'NR > 1{print $2}'" % full_file) == '0'):
                logging.warning("Read file %s has no data, skipping..." % file)
                read_list.append(Read(sample_name, None))
                continue
            is_merged = re.search('^(.*?)(?:[_\.](?:assembled|merged))+$', sample_name, re.IGNORECASE)
            if is_merged:
                sample_name = is_merged.group(1)
                read = Read(sample_name, [os.path.join(path, file)])
                read_list.append(read)
                logging.info(read)
            else:
                is_paired = re.search('^(?:((.*?)(?:_L\d\d\d)?(?:(?:[_\.](?:R(?:ead)?)?)))([12])([_\.])?)(?!.*[_\.](?:R(?:ead)?)?[12][_\.])(.*)$', sample_name, re.IGNORECASE)
                if is_paired:
                    if is_paired.group(3) == '1':  # If paired, only process read 1, so we don't double count the pair, see TODO below
                        sample_name = is_paired.group(2)
                        read1 = file
                        read2 = "%s2%s%s%s" % (is_paired.group(1), is_paired.group(4) or '', is_paired.group(5), is_read.group(2))
                        #print("\t%s\t%s\t%s" % (sample_name, read1, read2))
                        if os.path.exists(os.path.join(path, read2)):
                            read = Read(sample_name, [os.path.join(path, read1), os.path.join(path, read2)])
                            read_list.append(read)
                            logging.info(read)
                        else:
                            # TODO: If only R2 exists, it won't be included
                            logging.warning("Cannot find %s, the matching read to %s. Including as unpaired..." % (read2, read1))
                            read = Read(sample_name, [os.path.join(path, read1)])
                            read_list.append(read)
                            logging.info(read)
                else: #Read is unpaired
                    sample_name = is_read.group(1)
                    read = Read(sample_name, [os.path.join(path, file)])
                    read_list.append(read)
                    logging.info(read)
    return read_list

def findBams(path):
    import os
    import re
    from collections import namedtuple
    bam_list = []
    Bam = namedtuple('Bam', ['sample', 'bam', 'job'])
    for file in os.listdir(path):
        is_bam = re.search('(.*)-(.*)(\.bam)$', file, re.IGNORECASE)
        if is_bam:
            sample_name = is_bam.group(1)
            bam = Bam(sample_name, os.path.join(path, file), None)
            bam_list.append(bam)
            logging.info(bam)
    return bam_list

def expandPath(path):
    import re
    import os
    user_match = re.match('^(~)(.*)$', path)
    if user_match:
        path = os.path.expanduser(path)
    return os.path.abspath(path)

def indexFasta(fasta, aligner="bwa"):
    import os
    import re
    job_params = {'queue':'', 'mem_requested':2, 'num_cpus':1, 'walltime':4, 'args':''}
    job_params['name'] = "asap_index_%s" % os.path.basename(fasta)
    job_params['work_dir'] = os.path.dirname(fasta)
    if re.search('novo', aligner, re.IGNORECASE):
        command = "novoindex %s.idx %s" % (fasta, fasta)
    elif re.search('b(ow)?t(ie)?2', aligner, re.IGNORECASE):
        command = "bowtie2-build %s %s" % (fasta, os.path.splitext(fasta)[0])
    else:
        command = "bwa index %s" % (fasta)
    return _submit_job(job_manager, command, job_params)

def trimAdapters(sample, reads, outdir, quality=None, adapters="../illumina_adapters_all.fasta", minlen=80, dependency=None):
    from collections import namedtuple
    import os
    read1 = reads[0]
    read2 = reads[1] if len(reads) > 1 else None
    TrimmedRead = namedtuple('TrimmedRead', ['sample', 'jobid', 'reads'])
    trim_dir = os.path.join(outdir, 'trimmed')
    if not os.path.exists(trim_dir):
        os.makedirs(trim_dir)
    job_params = {'queue':'', 'mem_requested':6, 'num_cpus':4, 'walltime':8, 'args':''}
    job_params['name'] = "asap_trim_%s" % sample
    job_params['work_dir'] = trim_dir
    qual_string = quality if quality else ''
    if read2:
        out_reads1 = [sample + "_R1_trimmed.fastq.gz", sample + "_R1_unpaired.fastq.gz"]
        out_reads2 = [sample + "_R2_trimmed.fastq.gz", sample + "_R2_unpaired.fastq.gz"]
        out_reads = [os.path.join(trim_dir, out_reads1[0]), os.path.join(trim_dir, out_reads2[0])]
        command = "java -Xmx%sg org.usadellab.trimmomatic.Trimmomatic PE -threads %d %s %s %s %s %s %s ILLUMINACLIP:%s:4:30:10:1:true %s MINLEN:%d" % (job_params['mem_requested'], job_params['num_cpus'], read1, read2, out_reads1[0], out_reads1[1], out_reads2[0], out_reads2[1], adapters, qual_string, minlen)
    else:
        out_reads = [os.path.join(trim_dir, sample + "_trimmed.fastq.gz")]
        command = "java -Xmx%sg org.usadellab.trimmomatic.Trimmomatic SE -threads %d %s %s ILLUMINACLIP:%s:2:30:10 %s MINLEN:%d" % (job_params['mem_requested'], job_params['num_cpus'], read1, out_reads[0], adapters, qual_string, minlen)
    jobid = _submit_job(job_manager, command, job_params, (dependency,))
    return TrimmedRead(sample, jobid, out_reads)

def alignReadsToReference(sample, reads, reference, outdir, jobid=None, aligner="bowtie2", args=None):
    import re
    if re.search('novo', aligner, re.IGNORECASE):
        return _run_novoalign(sample, reads, reference, outdir, jobid, args=args)
    elif re.search('bwa', aligner, re.IGNORECASE):
        return _run_bwa(sample, reads, reference, outdir, jobid, args=args)
    else: #re.search('b(ow)?t(ie)?2', aligner, re.IGNORECASE)
        return _run_bowtie2(sample, reads, reference, outdir, jobid, args=args)

def processBam(sample_name, json_file, bam_file, xml_dir, dependency, depth, breadth, proportion, percid, mutdepth, smor=False):
    import os
    job_params = {'queue':'', 'mem_requested':4, 'num_cpus':2, 'walltime':24, 'args':''}
    job_params['name'] = "asap_bamprocesser_%s" % sample_name
    job_params['work_dir'] = xml_dir
    out_file = os.path.join(xml_dir, sample_name+".xml")
    smor_option = " -s" if smor else ""
    command = "bamProcessor -j %s -b %s -o %s -d %d --breadth %f -p %f -m %d -i %f%s" % (json_file, bam_file, out_file, depth, breadth, proportion, mutdepth, percid, smor_option)
    jobid = _submit_job(job_manager, command, job_params, (dependency,)) if dependency else _submit_job(job_manager, command, job_params)
    return (out_file, jobid)

def combineOutputFiles(run_name, xml_dir, out_dir, dependencies):
    import os
    job_params = {'queue':'', 'mem_requested':4, 'num_cpus':2, 'walltime':4, 'args':''}
    job_params['name'] = "asap_final_%s" % run_name
    job_params['work_dir'] = out_dir
    out_file = os.path.join(out_dir, run_name+"_analysis.xml")
    dependency_string = ":".join(dependencies)
    command = "outputCombiner -n %s -x %s -o %s" % (run_name, xml_dir, out_file)
    jobid = _submit_job(job_manager, command, job_params, (dependency_string, 'afterany'), notify=True)
    return (out_file, jobid)

if __name__ == '__main__':
    pass
