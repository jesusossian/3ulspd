#!/bin/bash
for model in randdpheurCC #esn std mc 3level partialmc echelon esn estp esls
do
	for bal in 0 1
	do
		for cap in 1.5 1.75 2.0
		do
			for N in 200 #50 100 200
			do
				for T in 15 30
				do
					for P in 1
					do
						for W in 5 10 #15 20
						do
							for V in DD_SF #SD_SF SD_DF DD_SF DD_DF
							do
								for id in 6
								do
									for dpalpha in 0.20 0.25
									do
										for dpredperiods in 0.7 0.85
										do
											for dpcostreductionret in 0.7 0.85
 											do
												for dpcostreductionware in 0.7 0.85
												do
													julia 3lspd.jl --inst instances/N${N}T${T}/N${N}T${T}P${P}W${W}${V}${id}.dat --form ${model} --capacity ${cap} --balanced ${bal} --dpalpha ${dpalpha} --dpredperiods ${dpredperiods} --dpcostreductionred ${dpcostreductionret} --dpcostreductionware ${dpcostreductionware}
												done
											done
										done
									done								
								done
								#mv saida.txt result/${version}_${style}_${model}_N${N}_T${T}_result.txt		
							done
						done
					done
				done
			done
		done
	done
done
