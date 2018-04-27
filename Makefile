
sub_proc:
	ls /exper-sw/cmst3/cmssw/users/olek/CMSSW_8_0_26_patch1/src/UserCode/NtuplerAnalyzer/queue_dir/${nt}/${proc}/${node}/com
	ssh fermi0${node} 'screen -d -m tcsh -c "source /exper-sw/cmst3/cmssw/users/olek/CMSSW_8_0_26_patch1/src/UserCode/NtuplerAnalyzer/queue_dir/${nt}/${proc}/${node}/com "'

#foreach d (`cat merge-sets/jobs_hadd.dtags.large_ab`)
#foreach d (`cat merge-sets/jobs_hadd.dtags.large_aa`)
#foreach d (`cat merge-sets/tt-sys/updowns.dtags`)
# merge-sets/jobs_hadd.dtags.single-muon
# merge-sets/jobs_hadd.dtags.single-electron

#echo foo & \
#

#make merge_proc dtags=merge-sets/jobs_hadd.dtags.large_ab nt=v25 proc=p3
merge_proc:
	for d in `cat ${dtags}`; \
	do \
	hadd merge-sets/${nt}/${proc}/$$d.root /lstore/cms/olek/outdirs/${nt}/${proc}/$$d/*.root & \
	done

merge_tt:
	for fl in `ls tt_files_*`; \
	do \
	hadd merge-sets/${nt}/${proc}/$$fl.root `sed 's,^,/lstore/cms/olek/outdirs/${nt}/${proc}/${dtag}/,' $$fl` & \
	done


resubmit:
	for js in `ls crab_projects/crab_Ntupler_${nt}* -d`; \
	do \
	echo $$js ; \
	crab resubmit -d $$js ; \
	done


status:
	for js in `ls crab_projects/crab_Ntupler_${nt}* -d`; \
	do \
	echo $$js ; \
	crab status -d $$js ; \
	done


