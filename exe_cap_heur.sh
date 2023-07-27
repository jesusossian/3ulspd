#!/bin/bash
maxt=120
bal=0  #0 unbal 1 bal
P=1
for fixsize in 2 3 4
do
    for freeint in 2 3 4
    do
	for model in rf
	do
	    for cap in 1.75 # 1.50 1.75 2.00
	    do
		for capw in 0.0 # 1.50 1.75 2.00
		do
		    for capr in 0.0 # 1.50 1.75 2.00
		    do
			for N in 50 100 200
			do
			    for T in 15 30
			    do
				for W in 5 10 15 20
				do
				    for V in SD_SF SD_DF DD_SF DD_DF
				    do
					echo "instance;NR;NT;NW;cap;capw;capr;bal;fixsizerf;freeintervalrf;maxtimerf;bestsol;time" >> saida.txt
					for id in $(seq 5)
					do
					    julia 3lspd.jl --inst instances/N${N}T${T}/N${N}T${T}P${P}W${W}${V}${id}.dat --form ${model} --balanced ${bal} --maxtimerf ${maxt} --fixsizerf ${fixsize} --freeintervalrf ${freeint} --capacity ${cap} --capacityw ${capw} --capacityr ${capr} >> report/bal${bal}_${model}_fixsize${fixsize}_freeint${freeint}_${model}_cap${cap}_capw${capw}_capr${capr}_N${N}T${T}P${P}W${W}${V}${id}.txt
					done
					mv saida.txt result/bal${bal}_fixsize${fixsize}_freeint${freeint}_${model}_cap${cap}_capw${capw}_capr${capr}_N${N}_T${T}_W${W}_${V}.csv
				    done
				done
			    done
			done
		    done
		done
	    done
	done
    done
done
