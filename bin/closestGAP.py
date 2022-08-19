#!/bin/python

import re, argparse, operator
from collections import defaultdict

def parseOptions():
	#use argparse to setup/get options from cmd line
	parser = argparse.ArgumentParser(description="In order to fragment scaffolds after HiC data assessment, find clostest gap of a scaff using fragments and gaps coordinates as inputs")
	parser.add_argument("-i", "--input", help="TAB table of scaff fragments coordinates on source assembly: /#scaff/#scaff_fragment/scaff_length/start/stop/", required=True)
	parser.add_argument("-d", "--directory", help="Directory path to find nCount.pl output files with gaps description", required=True)
	parser.add_argument("-s", "--suffix", help="Suffix of nCount.pl output files in the -d directory (scaffold_name+suffix)", required=True)
	return(parser.parse_args())

parameters=parseOptions()
# table of fragments coordinates on source assembly
with open(parameters.input, 'r') as f1:
	
	scaff_length=defaultdict(int)
	#scaff_fragment=defaultdict(list)

	scaff_start_fragment=defaultdict(lambda:defaultdict(str))
	scaff_starts=defaultdict(list)

	scaff_fragment_start=defaultdict(lambda:defaultdict(int))
	scaff_fragment_end=defaultdict(lambda:defaultdict(int))

	for li in f1.readlines():
		li=li.rstrip().split()
		scaff=li[0]
		fragment=li[1]
		length=int(li[2])
		fragment_start=int(li[3])
		fragment_end=int(li[4])
		
		#save HiC reviewed scaffolds data in dictionaries
		scaff_length[scaff]=length
		#scaff_fragment[scaff].append(fragment)

		scaff_start_fragment[scaff][fragment_start]=fragment
		scaff_starts[scaff].append(fragment_start)

		scaff_fragment_start[scaff][fragment]=fragment_start
		scaff_fragment_end[scaff][fragment]=fragment_end

#print output
output=open(re.sub(".txt", "_closestGAP.out", parameters.input), 'w')

#loop on scaffolds
for s in scaff_starts:
	print(" ")
	print("=================================================================")
	print(s+"   scaff_length: "+str(scaff_length[s]))

	#loop on fragments sorted by their start coordinate
	for start in sorted(scaff_starts[s]):
		#frag name given by its parent scaffold and its start coordinate
		frag=scaff_start_fragment[s][start]
		print("-----------------------------------------------------------------")
		print("OLD    "+frag+"    "+str(scaff_fragment_start[s][frag])+"    "+str(scaff_fragment_end[s][frag]))

		# table of gaps coordinates on source assembly
		with open(parameters.directory+"/"+s+parameters.suffix, 'r') as f2:

			scaff_newFRAG_start=defaultdict(lambda:defaultdict(int))
			scaff_newFRAG_end=defaultdict(lambda:defaultdict(int))

			#initialisation: diff(fragment_end-gap_start) set to a max
			diff=1000000000
			for li in f2.readlines():
				li=li.rstrip().split()
				gap=li[0]
				startGAP=int(li[2])
				endGAP=int(li[3])
				#print(gap+"    "+str(startGAP)+"    "+str(endGAP))
				
				# CASE 1 fragment=beginning of the scaff
				if scaff_fragment_start[s][frag]==1:
					# searching closest gap (min diff) to determine best end coord of the fragment
					current_diff=abs(scaff_fragment_end[s][frag]-startGAP)

					if current_diff<diff:
						diff=current_diff
						#new end is the nucl before start of closest GAP
						scaff_newFRAG_end[s][frag]=startGAP-1
						#start stay 1
						scaff_newFRAG_start[s][frag]=scaff_fragment_start[s][frag]
						#following fragment new start is nucl after end of closest GAP
						middle_fragment_start=endGAP+1
						last_fragment_start=endGAP+1
						targetGAP=li

				# CASE 2 fragment in middle of the scaff
				if scaff_fragment_start[s][frag]!=1 and scaff_fragment_end[s][frag]!=scaff_length[s] and len(scaff_starts[s])>2:
					# searching closest gap (min diff) to determine best end coord of the fragment
					current_diff=abs(scaff_fragment_end[s][frag]-startGAP)

					if current_diff<diff:
						diff=current_diff
						scaff_newFRAG_end[s][frag]=startGAP-1
						#new start is nucl after end of closest GAP found for former fragment
						scaff_newFRAG_start[s][frag]=middle_fragment_start
						#following fragment new start is nucl after end of closest GAP
						last_fragment_start=endGAP+1
						targetGAP=li

				# CASE 3 fragment=end of the scaff
				if scaff_fragment_end[s][frag]==scaff_length[s]:
					#new start is nucl after end of closest GAP found for former fragment
					scaff_newFRAG_start[s][frag]=last_fragment_start
					#end stay scaff length
					scaff_newFRAG_end[s][frag]=scaff_fragment_end[s][frag]

		print("NEW    "+frag+"    "+str(scaff_newFRAG_start[s][frag])+"    "+str(scaff_newFRAG_end[s][frag]))
		print(targetGAP)

		output.write("\t".join([s, frag, str(scaff_newFRAG_start[s][frag]), str(scaff_newFRAG_end[s][frag])])+"\n")