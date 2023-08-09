#!usr/bin/env nextflow
nextflow.enable.dsl=2


params.proteins='/project/arcc-students/protein_folding/EffectorWorkflow/Proteins/Asplac1/Asplac1_GeneCatalog_proteins_20140612_short.aa.fasta'	// Fasta input (defaults to shortened Aspergillus niger proteome)
params.output='/project/arcc-students/ecrane3/protfold/SigCazEff_out'					// Output directory
params.db_path='/project/arcc-students/FungalPathogens/Installs/db'					// Path to CAZyme database for run_dbcan
params.EffP_script='/project/arcc-students/FungalPathogens/Installs/EffectorP-3.0-main/EffectorP.py'	// Path to the EffectorP program


// Create a copy of original input fasta w/o special characters or other substrings that may interfere with sed editing of program outputs
process ReformatFasta {
	input:
	val params_proteins

	output:
	env new_fasta

	shell:
	'''
	x=!{params_proteins}
	x2=$(echo "${x##*/}")
	x3=$(echo "${x2%%.*}")
	new_fasta="$(pwd)/${x3}_numID.fasta"
	bash /project/arcc-students/protein_folding/EffectorWorkflow/bash_scripts/ReformatFasta.sh !{params_proteins} ${new_fasta}
	'''
}

// Run SignalP & output a string containing the path to its output file
process SignalP {
	input:
	val fasta_file
	val params_out

	output:
	env output
	env sig_data

	shell:
	'''
	x=!{fasta_file}
	x2=$(echo "${x##*/}")
	x3=$(echo "${x2%%_numID.fasta}")
	output="!{params_out}/${x3}"
	mkdir -p "${output}"

	signalp6 -ff !{fasta_file} -org eukarya -od ${output} -fmt none
	sig_data="${output}/prediction_results.txt"
	'''
}

// Create copy of input fasta w/o the non-secreted proteins identified by SignalP
process RemoveNonSecreted {
	input:
	val sig_data
	val fasta_file

	output:
	path 'nonsec_numIDs.txt'
	path 'lowconfsec_numIDs.txt'
	env filename
	env hc_filename

	shell:
	'''
	## Copy input fasta to new file 
	x=!{fasta_file}
	x2=$(echo "${x##*/}")
	x3=$(echo "${x2%%.*}")
	filename="$(pwd)/${x3}_secreted.fasta"
	hc_filename="$(pwd)/${x3}_hiconf_secreted.fasta"
	cp ${x} ${filename} 

	cp !{sig_data} nonsec_numIDs.txt 
	cp !{sig_data} lowconfsec_numIDs.txt				## Copy output from SignalP into new files
	sed -i '1,2d;/OTHER/!d;s/[[:blank:]].*//' nonsec_numIDs.txt      ## and pare these files down to just the non-secreted protein sequence IDs
	sed -i '1,2d;/OTHER/d;s/..........CS.*//' lowconfsec_numIDs.txt ## and the low-confidence predictions of secreted protein sequence IDs

	max=$(sed -n "$=" lowconfsec_numIDs.txt)
	for ((j=1;j<=${max};j++))
	do
		str=$(head -n $j lowconfsec_numIDs.txt | tail -1)
		conf=${str: -6}
		if [ ${conf} -lt 001000 ]
		then
			sed -i "$(echo $j)d" lowconfsec_numIDs.txt
			j=$(($j - 1))
			max=$((${max} - 1))
		 fi
	done
	sed -i 's/[[:blank:]].*//' lowconfsec_numIDs.txt	

	## Remove each non-secreted protein sequence from the copy of the input fasta file
	max=$(sed -n "$=" nonsec_numIDs.txt)
	for ((j=1;j<=${max};j++))
	do
	        sed -i "/^>$(head -n $j nonsec_numIDs.txt | tail -1)$/,/^>/{/^>/!d};/>$(head -n $j nonsec_numIDs.txt | tail -1)$/d" ${filename}
	done

	## Create a copy of the shortened protein list above & remove the low-confidence predictions by SignalP from it
	cp ${filename} ${hc_filename}
	max=$(sed -n "$=" lowconfsec_numIDs.txt)
	for ((j=1;j<=${max};j++))
	do
		sed -i "/^>$(head -n $j lowconfsec_numIDs.txt | tail -1)$/,/^>/{/^>/!d};/>$(head -n $j lowconfsec_numIDs.txt | tail -1)$/d" ${hc_filename}
	done
	'''
}

// Run run_dbcan (cazyme identification software) & output a string containing the path to its output file
process run_dbcan {
	input:
	path fasta_file
	val output
	val db_path

	output:
	env Caz_data

	shell:
	'''
	run_dbcan !{fasta_file} protein --out_dir !{output} --db_dir !{db_path}
	Caz_data="!{output}/overview.txt"
	'''
}

// Create copy of input fasta w/o the cazymes identified by run_dbcan
process RemoveCazymes {
	input:
	val Caz_data
	val fasta_file
	val hc_fasta_file

	output:
	path 'cazyme_numIDs.txt'
	env filename
	env hc_filename

	shell:
	''' 
	cp !{Caz_data} cazyme_numIDs.txt		## Copy output of run_dbcan into two new files
	sed -i '1d;s/	.*//' cazyme_numIDs.txt		## Pare the new files down to just CAZyme sequence IDs

	## Create copies of input fasta files
	x=!{fasta_file}
	x2=$(echo "${x##*/}")
	x3=$(echo "${x2%%.*}")
	filename="$(pwd)/${x3}_nocazymes.fasta"
	cp ${x} ${filename}

	y=!{hc_fasta_file}
        y2=$(echo "${y##*/}")
        y3=$(echo "${y2%%.*}")
	hc_filename="$(pwd)/${y3}_nocazymes.fasta"
	cp ${y} ${hc_filename}

	## Remove each CAZyme sequence from the copy
	a=$(sed -n "$=" cazyme_numIDs.txt)
	for ((j=1;j<=${a};j++))
	do
	        sed -i "/^>$(head -n $j cazyme_numIDs.txt | tail -1)$/,/^>/{/^>/!d};/>$(head -n $j cazyme_numIDs.txt | tail -1)$/d" ${filename} ${hc_filename}
	done
	'''
}

// Run EffectorP
process EffectorP {
	input:
	path ep_script
	val fasta_file

	output:
	stdout emit: table

	shell:
	'''
	python !{ep_script} -i !{fasta_file}
	'''
}

// Create a copy of the input fasta w/o the effector proteins identified by EffectorP
process RemoveNonEffectors {
	input:
	val TableData
	val fasta_file
	val hc_fasta_file

	output:
	path 'EffectorP_output_num.txt'
	path 'noneff_numIDs.txt'
	path 'lowconfeff_numIDs.txt'
	env filename
	env hc_filename

	shell:
	'''
	echo "!{TableData}" >> EffectorP_output_num.txt		## Put output data from EffectorP into a file
	cp EffectorP_output_num.txt noneff_numIDs.txt		## Copy the output into two new files
	cp EffectorP_output_num.txt lowconfeff_numIDs.txt
	sed -i -e :a -e '$d;N;2,12ba' -e 'P;D' noneff_numIDs.txt lowconfeff_numIDs.txt	## and pare the new files down to just the non-effectors
	sed -i '1,11d;/Non-effector/!d;s/- .*//;s/[[:blank:]]*$//' noneff_numIDs.txt	## & low-confidence effectors, respectively
	sed -i '1,11d;/Non-effector/d;s/Cyt.*//;s/Apo.*//;s/)//g;s/(//g;s/-/ /g;s/[[:blank:]]*$//' lowconfeff_numIDs.txt

	## Create copy of the input fasta file for regular effector list
	x=!{fasta_file}
	x2=$(echo "${x##*/}")
	x3=$(echo "${x2%%.*}")
	filename="$(pwd)/${x3}_effectors.fasta"
	cp ${x} ${filename}

	## Remove each non-effector from the first copy
	a=$(sed -n "$=" noneff_numIDs.txt)
	for ((j=1;j<=${a};j++))
	do
	        sed -i "/^>$(head -n $j noneff_numIDs.txt | tail -1)$/,/^>/{/^>/!d};/>$(head -n $j noneff_numIDs.txt | tail -1)$/d" ${filename}
	done

	## Create copy of the input fasta file for high-confidence effector list
	y=!{hc_fasta_file}
	y2=$(echo "${y##*/}")
	y3=$(echo "${y2%%.*}")
	hc_filename="$(pwd)/${y3}_effectors.fasta"
	cp ${filename} ${hc_filename}

	## Remove each low-confidence effector from the high-confidence effector sequence ID list
	max=$(sed -n "$=" lowconfeff_numIDs.txt)
	for ((j=1;j<=${max};j++))
	do
		str=$(awk "NR==$j" lowconfeff_numIDs.txt | awk -F "." '{print $NF}')
		for ((k=${#str};k<=2;k++))
		do
			str="${str}0"
		done
		alt_str=$(awk "NR==$j" lowconfeff_numIDs.txt | awk -F "." '{print $2}')
		alt_str=${alt_str%% *}
		for ((k=${#alt_str};k<=2;k++))
		do
			alt_str="${alt_str}0"
		done
		
		if [ $str -gt 800 ] || [ $alt_str -gt 800 ]
		then
			sed -i "$(echo $j)d" lowconfeff_numIDs.txt
			j=$(($j - 1))
			max=$(($max - 1))
		fi
	done
	sed -i "s/[[:blank:]].*//" lowconfeff_numIDs.txt

	## Remove each non-effector from the first copy
	a=$(sed -n "$=" noneff_numIDs.txt)
	for ((j=1;j<=${a};j++))
	do
		sed -i "/^>$(head -n $j lowconfeff_numIDs.txt | tail -1)$/,/^>/{/^>/!d};/>$(head -n $j lowconfeff_numIDs.txt | tail -1)$/d" ${hc_filename}
	done
	'''
}

process RestoreSequenceIDs {
	input:
	val params_proteins
	val numID_file

	output:
	stdout

	shell:
	'''
	bash /project/arcc-students/protein_folding/EffectorWorkflow/bash_scripts/RestoreSeqID.sh !{params_proteins} !{numID_file}
	'''
}


// Move a file from its current directory to a new path
process RelocateFiles {
	input:
	path oldpath
	path newpath

	output:
	stdout

	shell:
	'''
	cp !{oldpath} !{newpath}
	'''
}

// Move a file from its current directory to a new path
process RelocateMiscFiles {
	input:
	path oldpath
	path newpath

	output:
	stdout

	shell:
	'''
	mkdir -p !{newpath}/numID_files
	cp !{oldpath} !{newpath}/numID_files
	'''
}


workflow {
	proteins_prepped=ReformatFasta(params.proteins)

	(new_output,SigP_data)=SignalP(proteins_prepped,params.output)
	(nonsecIDs,lowconfsecIDs,proteins_2,proteins_2hc)=RemoveNonSecreted(SigP_data,proteins_prepped)
	
	Caz_data=run_dbcan(proteins_2,new_output,params.db_path)
	(cazymeIDs,proteins_3,proteins_3hc)=RemoveCazymes(Caz_data,proteins_2,proteins_2hc)
	
	EffectorP(params.EffP_script,proteins_3)
	(EffP_out,noneffIDs,lowconfeffIDs,proteins_final,proteins_final_hc)=RemoveNonEffectors(EffectorP.out.table,proteins_3,proteins_3hc)

	numID_files=proteins_2.concat(nonsecIDs,lowconfsecIDs,proteins_2hc,cazymeIDs,proteins_3,proteins_3hc,noneffIDs,lowconfeffIDs,EffP_out,proteins_final,proteins_final_hc)
	bread=RestoreSequenceIDs(params.proteins,numID_files)
	RelocateFiles(bread,new_output)
	RelocateMiscFiles(numID_files,new_output)

}
