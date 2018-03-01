
sub_proc:
	ls /exper-sw/cmst3/cmssw/users/olek/CMSSW_8_0_26_patch1/src/UserCode/NtuplerAnalyzer/queue_dir/${nt}/${proc}/${node}/com
	ssh fermi0${node} 'screen -d -m tcsh -c "source /exper-sw/cmst3/cmssw/users/olek/CMSSW_8_0_26_patch1/src/UserCode/NtuplerAnalyzer/queue_dir/${nt}/${proc}/${node}/com "'

