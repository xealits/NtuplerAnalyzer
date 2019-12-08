
#foreach d (`cat merge-sets/jobs_hadd.dtags.large_ab`)
#foreach d (`cat merge-sets/jobs_hadd.dtags.large_aa`)
#foreach d (`cat merge-sets/tt-sys/updowns.dtags`)
# merge-sets/jobs_hadd.dtags.single-muon
# merge-sets/jobs_hadd.dtags.single-electron

#echo foo & \
#


submit: job_file=cur_jobs
submit:
	for c in `cat ${job_file}`; \
	do \
	crab submit -c $$c ; \
	done

resubmit: grep=.
resubmit:
	for js in `ls crab_projects/crab_Ntupler_${nt}* -d | grep ${grep}`; \
	do \
	echo $$js ; \
	crab resubmit -d $$js ; \
	done

status: grep=.
status:
	for js in `ls crab_projects/crab_Ntupler_${nt}* -d | grep ${grep}`; \
	do \
	echo $$js ; \
	crab status -d $$js ; \
	done


