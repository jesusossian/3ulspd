module BranchAndCut

using JuMP
using Gurobi
using CPLEX
using Data
using Parameters
using MathProgBase

import DPHeuristics

mutable struct stdFormVars
	xp
	xr
	xw
	yp
	yr
	yw
	sp
	sr
	sw
end


export standardFormulation, stdFormVars, schelonStockFormulation

function bcstandardFormulation(inst::InstanceData, params::ParameterData)

    println("Running BranchAndCut.bcstandardFormulation")

	if params.solver == "Gurobi"
		if params.disablesolver == 1 #Disable gurobi cuts and presolve
			if params.maxnodes < 999.0
				model = Model(solver = GurobiSolver(TimeLimit=params.maxtime,MIPGap=params.tolgap,CliqueCuts=0, CoverCuts=0, FlowCoverCuts=0,FlowPathCuts=0,MIRCuts=0,NetworkCuts=0,GomoryPasses=0, PreCrush=1,NodeLimit=params.maxnodes))
			else
				model = Model(solver = GurobiSolver(TimeLimit=params.maxtime,MIPGap=params.tolgap,CliqueCuts=0, CoverCuts=0, FlowCoverCuts=0,FlowPathCuts=0,MIRCuts=0,NetworkCuts=0,GomoryPasses=0, PreCrush=1))
			end

			#model = Model(solver = GurobiSolver(TimeLimit=params.maxtime,MIPGap=params.tolgap,CliqueCuts=0, CoverCuts=0, FlowCoverCuts=0,FlowPathCuts=0,MIRCuts=0,NetworkCuts=0,GomoryPasses=0, Presolve=0,NodeLimit=1,CutPasses=5,OutputFlag=0))
			#model = Model(solver = GurobiSolver(TimeLimit=params.maxtime,MIPGap=params.tolgap,CliqueCuts=0, CoverCuts=0, FlowCoverCuts=0,FlowPathCuts=0,MIRCuts=0,NetworkCuts=0,GomoryPasses=1,ImpliedCuts=0, Presolve=0))
			#model = Model(solver = GurobiSolver(TimeLimit=params.maxtime,MIPGap=params.tolgap,CliqueCuts=0, CoverCuts=0, FlowCoverCuts=0,FlowPathCuts=0,MIRCuts=0,NetworkCuts=0,GomoryPasses=0,ImpliedCuts=0, Presolve=0))
			#model = Model(solver = GurobiSolver(MIPGap=params.tolgap,Cuts=0, Presolve=0))
		else
			if params.maxnodes < 999.0
        		model = Model(solver = GurobiSolver(TimeLimit=params.maxtime,MIPGap=params.tolgap,PreCrush=1,NodeLimit=params.maxnodes))
			else
				model = Model(solver = GurobiSolver(TimeLimit=params.maxtime,MIPGap=params.tolgap,PreCrush=1))
			end
		end
	elseif params.solver == "Cplex"
        model = Model(solver = CplexSolver(CPX_PARAM_TILIM=params.maxtime,CPX_PARAM_EPGAP=params.tolgap,CPX_PARAM_PREIND=0))
	else
		println("No solver selected")
		return 0
	end

	### Defining variables ###
    @variable(model,0 <= xp[i=1:inst.NP,t=1:inst.NT] <= Inf) #sum(inst.DP[i,t] for t=1:inst.NT))
    @variable(model,0 <= xr[i=1:inst.NR,t=1:inst.NT] <= Inf) #sum(inst.DR[i,t] for t=1:inst.NT))
    @variable(model,0 <= xw[i=1:inst.NW,t=1:inst.NT] <= Inf) #sum(inst.D[i,t] for t=1:inst.NT))
    @variable(model, yp[i=1:inst.NP,t=1:inst.NT], Bin)
    @variable(model, yr[i=1:inst.NR,t=1:inst.NT], Bin)
    @variable(model, yw[i=1:inst.NW,t=1:inst.NT], Bin)
    @variable(model,0 <= sp[i=1:inst.NP,t=0:inst.NT] <= Inf) #sum(inst.DP[i,t] for t=1:inst.NT))
    @variable(model,0 <= sr[i=1:inst.NR,t=0:inst.NT] <= Inf) #sum(inst.DR[i,t] for t=1:inst.NT))
    @variable(model,0 <= sw[i=1:inst.NW,t=0:inst.NT] <= Inf)


	allretcumdem = zeros(Int,inst.NT,inst.NT)
	for t in 1:inst.NT,t2 in t:inst.NT
		allretcumdem[t,t2] = sum(inst.cumdem[:,t,t2])
	end


	warecumdem = zeros(Int,inst.NW,inst.NT,inst.NT)
	for w in 1:inst.NW
		for t in 1:inst.NT,t2 in t:inst.NT
			for k in 1:length(inst.DeltaW[w])
				r = inst.DeltaW[w][k]
				warecumdem[w,t,t2] += inst.cumdem[r,t,t2]
			end
		end
	end

println("Defined variables")
	### Objective function ###
	@objective(model, Min,
        sum(inst.SCP[i,t]*yp[i,t] + inst.HCP[i]*sp[i,t] for i=1:inst.NP, t=1:inst.NT)
        + sum(inst.SCR[i,t]*yr[i,t] + inst.HCR[i]*sr[i,t] for i=1:inst.NR, t=1:inst.NT)
        + sum(inst.SCW[i,t]*yw[i,t] + inst.HCW[i]*sw[i,t] for i=1:inst.NW, t=1:inst.NT)
	)

println("Finished objective function")

	### Setup constraints ###
	@constraint(model,
				setupP[t=1:inst.NT], xp[1,t] <= sum(inst.DP[1,k] for k in t:inst.NT)*yp[1,t]
				)
	@constraint(model,
				setupW[i=1:inst.NW, t=1:inst.NT], xw[i,t] <= sum(inst.DW[i,k] for k in t:inst.NT)*yw[i,t]
				)
	@constraint(model,
				setupR[i=1:inst.NR, t=1:inst.NT], xr[i,t] <= sum(inst.D[i,k] for k in t:inst.NT)*yr[i,t]
				)

	println("Finished setup constraints")

	### no initial inventory
	@constraint(model,
				zeroinvP[p=1:inst.NP],
				sp[p,0] == 0
				)
	@constraint(model,
				zeroinvW[w=1:inst.NW],
				sw[w,0] == 0
				)
	@constraint(model,
				zeroinvR[r=1:inst.NR],
				sr[r,0] == 0
				)

println("Finished zero inventory constraints")

	### Balance constraints ###
	@constraint(model,
				balanceP[i=1:inst.NP, t=1:inst.NT],
				sp[i,t-1] + xp[i,t] == sum(xw[w,t] for w in 1:inst.NW) + sp[i,t]
				)
	@constraint(model,
				balanceW[w=1:inst.NW,t=1:inst.NT],
				sw[w,t-1] + xw[w,t] == sum(xr[inst.DeltaW[w][k],t]   for k in 1:length(inst.DeltaW[w])) + sw[w,t]
				)
	@constraint(model,
				balanceR[i=1:inst.NR,t=1:inst.NT],
				sr[i,t-1] + xr[i,t] == inst.D[i,t] + sr[i,t]
				)


	if params.capacity != 0.0
		@constraint(model,
					capP[t=1:inst.NT], xp[1,t] <= min(inst.C[t],sum(inst.DP[1,k] for k in t:inst.NT))*yp[1,t]
					)
	end

	### Capacity constraints ###
#	@constraint(model, capacity[i=1:inst.NP, t=1:inst.NT], xp[i,t] <= min(inst.DP[i,t],inst.C[t])*yp[i,t])

	#writeLP(model,"modelo.lp",genericnames=false)


	if params.form == "randdpheurbcstd"
		println("HERE")
		SETR,SETW,SETP = DPHeuristics.RandomizedDPHeuristicBottomUp(inst,params)

		ypsol = zeros(inst.NP,inst.NT)
		ywsol = zeros(inst.NW,inst.NT)
		yrsol = zeros(inst.NR,inst.NT)

		for t in 1:length(SETP)
			ypsol[1,SETP[t]] = 1
		end
		for w in 1:inst.NW, t in 1:length(SETW[w])
			ywsol[w,SETW[w][t]] = 1
		end
		for r in 1:inst.NR, t in 1:length(SETR[r])
			yrsol[r,SETR[r][t]] = 1
		end

		for t = 1:inst.NT
			JuMP.setvalue(yp[1,t],ypsol[1,t])
		end
		for w=1:inst.NW, t=1:inst.NT
			JuMP.setvalue(yw[w,t], ywsol[w,t])
		end
		for r=1:inst.NR, t=1:inst.NT
			JuMP.setvalue(yr[r,t],yrsol[r,t])
		end
	end


	S = zeros(Bool,inst.NR,inst.NT)
	function slowsingleretailerseparator(cb)
		node = MathProgBase.cbgetexplorednodes(cb)

		addedcuts = 0
		#println("Nodes = ",node)

		if node < 1
			xp_vals = getvalue(xp)
			yp_vals = getvalue(yp)
			xw_vals = getvalue(xw)
			yw_vals = getvalue(yw)
			xr_vals = getvalue(xr)
			yr_vals = getvalue(yr)

			EPSTOL = params.tolsep

			for r in 1:inst.NR

				for l in 1:inst.NT
					alpha = 0.0
					for k in 1:l
						if yp_vals[1,k]*inst.cumdem[r,k,l] < xp_vals[1,k]
							S[r,k] = true
							alpha += yp_vals[1,k]*inst.cumdem[r,k,l]
						else
							S[r,k] = false
							alpha +=  xp_vals[1,k]
						end
					end
					#println("viol = ",inst.cumdem[r,1,l]-alpha)
					if alpha < inst.cumdem[r,1,l] - EPSTOL
#						println("r = $(r), l = $(l), violation = ",inst.cumdem[r,1,l]-alpha)

						expr = 0
						for k in 1:l
							if S[r,k]
								expr += yp[1,k]*inst.cumdem[r,k,l]
							else
								expr += xp[1,k]
							end
						end
						@usercut(cb, expr >= inst.cumdem[r,1,l])
						addedcuts += 1

					end

				end
			end


			for r in 1:inst.NR

				for l in 1:inst.NT
					alpha = 0.0
					for k in 1:l
						if yw_vals[inst.DeltamR[r],k]*inst.cumdem[r,k,l] < xw_vals[inst.DeltamR[r],k]
							S[r,k] = true
							alpha += yw_vals[inst.DeltamR[r],k]*inst.cumdem[r,k,l]
						else
							S[r,k] = false
							alpha +=  xw_vals[inst.DeltamR[r],k]
						end
					end
					#println("R viol = ",inst.cumdem[r,1,l]-alpha)
					if alpha < inst.cumdem[r,1,l] - EPSTOL
#						println("W r = $(r), l = $(l), violation = ",inst.cumdem[r,1,l]-alpha)

						expr = 0
						for k in 1:l
							if S[r,k]
								expr += yw[inst.DeltamR[r],k]*inst.cumdem[r,k,l]
							else
								expr += xw[inst.DeltamR[r],k]
							end
						end
						@usercut(cb, expr >= inst.cumdem[r,1,l])
						addedcuts += 1
					end

				end
			end


			for r in 1:inst.NR

				for l in 1:inst.NT
					alpha = 0.0
					for k in 1:l
						if yr_vals[r,k]*inst.cumdem[r,k,l] < xr_vals[r,k]
							S[r,k] = true
							alpha += yr_vals[r,k]*inst.cumdem[r,k,l]
						else
							S[r,k] = false
							alpha +=  xr_vals[r,k]
						end
					end
					#println("R viol = ",inst.cumdem[r,1,l]-alpha)
					if alpha < inst.cumdem[r,1,l] - EPSTOL
#						println("R r = $(r), l = $(l), violation = ",inst.cumdem[r,1,l]-alpha)

						expr = 0
						for k in 1:l
							if S[r,k]
								expr += yr[r,k]*inst.cumdem[r,k,l]
							else
								expr += xr[r,k]
							end
						end
						@usercut(cb, expr >= inst.cumdem[r,1,l])
						addedcuts += 1

					end

				end
			end

			println("Added cuts: $(addedcuts)")

		end

	end  # End of callback function



	S = zeros(Bool,inst.NR,inst.NT)
	SP = zeros(Bool,inst.NT)
	SW = zeros(Bool,inst.NW,inst.NT)





	function slowallretailerseparatornew(cb)
		node = MathProgBase.cbgetexplorednodes(cb)

		addedcuts = 0
		#println("Nodes = ",node)

		if node < 1 && cuttingplanerounds <= params.cpmaxrounds


			cuttingplanerounds += 1

			#println("cuttingplanerounds = $(cuttingplanerounds), bestbound = $(MathProgBase.cbgetbestbound(cb))")

			xp_vals = getvalue(xp)
			yp_vals = getvalue(yp)
			xw_vals = getvalue(xw)
			yw_vals = getvalue(yw)
			xr_vals = getvalue(xr)
			yr_vals = getvalue(yr)

			EPSTOL = params.tolsep

			if true

				#P
				for l in 1:inst.NT
					alpha = 0.0
					for k in 1:l
						if yp_vals[1,k]*allretcumdem[k,l] < xp_vals[1,k]
							SP[k] = true
							alpha += yp_vals[1,k]*allretcumdem[k,l]
						else
							SP[k] = false
							alpha +=  xp_vals[1,k]
						end
					end
					#println("viol = ",inst.cumdem[r,1,l]-alpha)
					if alpha < allretcumdem[1,l] - EPSTOL
					#		println("HERE P -> l = $(l), violation = ",allretcumdem[1,l]-alpha)

						expr = 0
						for k in 1:l
							if SP[k]
								expr += yp[1,k]*allretcumdem[k,l]
							else
								expr += xp[1,k]
							end
						end
						@usercut(cb, expr >= allretcumdem[1,l])
						addedcuts += 1

					end

				end



				#W
				for w in 1:inst.NW

					for l in 1:inst.NT
						alpha = 0.0
						for k in 1:l
							if yw_vals[w,k]*warecumdem[w,k,l] < xw_vals[w,k]
								SW[w,k] = true
								alpha += yw_vals[w,k]*warecumdem[w,k,l]
							else
								SW[w,k] = false
								alpha +=  xw_vals[w,k]
							end
						end
						#println("R viol = ",inst.cumdem[r,1,l]-alpha)
						if alpha < warecumdem[w,1,l] - EPSTOL
		#						println("W r = $(r), l = $(l), violation = ",inst.cumdem[r,1,l]-alpha)

							expr = 0
							for k in 1:l
								if SW[w,k]
									expr += yw[w,k]*warecumdem[w,k,l]
								else
									expr += xw[w,k]
								end
							end
							@usercut(cb, expr >= warecumdem[w,1,l])
							addedcuts += 1
						end

					end
				end




				# R
				for r in 1:inst.NR

					for l in 1:inst.NT
						alpha = 0.0
						for k in 1:l
							if yr_vals[r,k]*inst.cumdem[r,k,l] < xr_vals[r,k]
								S[r,k] = true
								alpha += yr_vals[r,k]*inst.cumdem[r,k,l]
							else
								S[r,k] = false
								alpha +=  xr_vals[r,k]
							end
						end
						#println("R viol = ",inst.cumdem[r,1,l]-alpha)
						if alpha < inst.cumdem[r,1,l] - EPSTOL
		#						println("R r = $(r), l = $(l), violation = ",inst.cumdem[r,1,l]-alpha)

							expr = 0
							for k in 1:l
								if S[r,k]
									expr += yr[r,k]*inst.cumdem[r,k,l]
								else
									expr += xr[r,k]
								end
							end
							@usercut(cb, expr >= inst.cumdem[r,1,l])
							addedcuts += 1

						end

					end
				end
			end


			if mod(cuttingplanerounds,5)==0

				# P+W
				for l in 2:inst.NT
					for lp in 1:l-1
						alpha = 0.0
						for k in 1:lp
							if yp_vals[1,k]*allretcumdem[k,l] < xp_vals[1,k]
								SP[k] = true
								alpha += yp_vals[1,k]*allretcumdem[k,l]
							else
								SP[k] = false
								alpha +=  xp_vals[1,k]
							end
						end
						#println("viol = ",inst.cumdem[r,1,l]-alpha)

						for k in lp+1:l
							for w in 1:inst.NW
								if yw_vals[w,k]*warecumdem[w,k,l] < xw_vals[w,k]
									SW[w,k] = true
									alpha += yw_vals[w,k]*warecumdem[w,k,l]
								else
									SW[w,k] = false
									alpha +=  xw_vals[w,k]
								end
							end
						end




						if alpha < allretcumdem[1,l] - EPSTOL
		#						println("r = $(r), l = $(l), violation = ",inst.cumdem[r,1,l]-alpha)

							expr = 0
							for k in 1:lp
								if SP[k]
									expr += yp[1,k]*allretcumdem[k,l]
								else
									expr += xp[1,k]
								end
							end

							for k in lp+1:l
								for w in 1:inst.NW
									if SW[w,k]
										expr += yw[w,k]*warecumdem[w,k,l]
									else
										expr += xw[w,k]
									end
								end
							end


							@usercut(cb, expr >= allretcumdem[1,l])
							addedcuts += 1

						end
					end
				end



				# P+R #apparently, these do not help at all
				for l in 2:inst.NT
					for lp in 1:l-1
						alpha = 0.0
						for k in 1:lp
							if yp_vals[1,k]*allretcumdem[k,l] < xp_vals[1,k]
								SP[k] = true
								alpha += yp_vals[1,k]*allretcumdem[k,l]
							else
								SP[k] = false
								alpha +=  xp_vals[1,k]
							end
						end
						#println("viol = ",inst.cumdem[r,1,l]-alpha)

						for k in lp+1:l
							for r in 1:inst.NR
								if yr_vals[r,k]*inst.cumdem[r,k,l] < xr_vals[r,k]
									S[r,k] = true
									alpha += yr_vals[r,k]*inst.cumdem[r,k,l]
								else
									S[r,k] = false
									alpha +=  xr_vals[r,k]
								end
							end
						end


						if alpha < allretcumdem[1,l] - EPSTOL
							#println("l = $(l), violation = ",allretcumdem[1,l]-alpha)

							expr = 0
							for k in 1:lp
								if SP[k]
									expr += yp[1,k]*allretcumdem[k,l]
								else
									expr += xp[1,k]
								end
							end

							for k in lp+1:l
								for r in 1:inst.NR
									if S[r,k]
										expr += yr[r,k]*inst.cumdem[r,k,l]
									else
										expr += xr[r,k]
									end
								end
							end


							@usercut(cb, expr >= allretcumdem[1,l])
							addedcuts += 1

						end
					end
				end








				#println("********************* W+R ***********************")

				for w in 1:inst.NW

					for l in 2:inst.NT

						for lw in 1:l-1
							alpha = 0.0
							for k in 1:lw
				#				println("r= $(r), lw= $(lw), l= $(l)  -> k = $(k)")
								if yw_vals[w,k]*warecumdem[w,k,l] < xw_vals[w,k]
									SW[w,k] = true
									alpha += yw_vals[w,k]*warecumdem[w,k,l]
								else
									SW[w,k] = false
									alpha +=  xw_vals[w,k]
								end
							end
							#println("R viol = ",inst.cumdem[r,1,l]-alpha)
							for ret in 1:length(inst.DeltaW[w])
								r = inst.DeltaW[w][ret]
								for k in lw+1:l
									if yr_vals[r,k]*inst.cumdem[r,k,l] < xr_vals[r,k]
										S[r,k] = true
										alpha += yr_vals[r,k]*inst.cumdem[r,k,l]
									else
										S[r,k] = false
										alpha +=  xr_vals[r,k]
									end
								end
							end
							if alpha < warecumdem[w,1,l] - EPSTOL
				#				println("WR r = $(r), l = $(l), violation = ",inst.cumdem[r,1,l]-alpha)

								expr = 0
								for k in 1:lw
									if SW[w,k]
										expr += yw[w,k]*warecumdem[w,k,l]
									else
										expr += xw[w,k]
									end
								end

								#println("R r = $(r), l = $(l), violation = ",inst.cumdem[r,1,l]-alpha)
								for ret in 1:length(inst.DeltaW[w])
									r = inst.DeltaW[w][ret]
									for k in lw+1:l
										if S[r,k]
											expr += yr[r,k]*inst.cumdem[r,k,l]
										else
											expr += xr[r,k]
										end
									end
								end
								@usercut(cb, expr >= warecumdem[w,1,l])
								addedcuts += 1

							end
						end
					end
				end


			end


			if mod(cuttingplanerounds,10)==0

				#println("********************* P+W+R ***********************")

				for l in 3:inst.NT
					for lp in 1:l-2
						for lw in lp+1:l-1
							alpha = 0.0

							for k in 1:lp
								if yp_vals[1,k]*allretcumdem[k,l] < xp_vals[1,k]
									SP[k] = true
									alpha += yp_vals[1,k]*allretcumdem[k,l]
								else
									SP[k] = false
									alpha +=  xp_vals[1,k]
								end
							end

							for w in 1:inst.NW
								for k in lp+1:lw
									if yw_vals[w,k]*warecumdem[w,k,l] < xw_vals[w,k]
										SW[w,k] = true
										alpha += yw_vals[w,k]*warecumdem[w,k,l]
									else
										SW[w,k] = false
										alpha +=  xw_vals[w,k]
									end
								end
							end
							for r in 1:inst.NR
								for k in lw+1:l
									if yr_vals[r,k]*inst.cumdem[r,k,l] < xr_vals[r,k]
										S[r,k] = true
										alpha += yr_vals[r,k]*inst.cumdem[r,k,l]
									else
										S[r,k] = false
										alpha +=  xr_vals[r,k]
									end
								end
							end
							#println("R viol = ",inst.cumdem[r,1,l]-alpha)


							if alpha < allretcumdem[1,l] - EPSTOL
								#println("PWR r = $(r), l = $(l), violation = ",inst.cumdem[r,1,l]-alpha)

								expr = 0
								for k in 1:lp
									if SP[k]
										expr += yp[1,k]*allretcumdem[k,l]
									else
										expr += xp[1,k]
									end
								end

								#println("R r = $(r), l = $(l), violation = ",inst.cumdem[r,1,l]-alpha)

								for w in 1:inst.NW
									for k in lp+1:l
										if SW[w,k]
											expr += yw[w,k]*warecumdem[w,k,l]
										else
											expr += xw[w,k]
										end
									end
								end

								for r in 1:inst.NR
									for k in lw+1:l
										if S[r,k]
											expr += yr[r,k]*inst.cumdem[r,k,l]
										else
											expr += xr[r,k]
										end
									end
								end

								@usercut(cb, expr >= allretcumdem[1,l])
								addedcuts += 1

							end
						end
					end
				end

			end #if true three level inequalities



			#println("Separated cuts: $(addedcuts)")

		end

	end  # End of callback function slowallretailerseparatornew








	function slowallretailerseparator(cb)
		node = MathProgBase.cbgetexplorednodes(cb)

		addedcuts = 0
		#println("Nodes = ",node)

		if node < 1
			xp_vals = getvalue(xp)
			yp_vals = getvalue(yp)
			xw_vals = getvalue(xw)
			yw_vals = getvalue(yw)
			xr_vals = getvalue(xr)
			yr_vals = getvalue(yr)

			EPSTOL = params.tolsep

			#P
			for l in 1:inst.NT
				alpha = 0.0
				for k in 1:l
					if yp_vals[1,k]*allretcumdem[k,l] < xp_vals[1,k]
						SP[k] = true
						alpha += yp_vals[1,k]*allretcumdem[k,l]
					else
						SP[k] = false
						alpha +=  xp_vals[1,k]
					end
				end
				#println("viol = ",inst.cumdem[r,1,l]-alpha)
				if alpha < allretcumdem[1,l] - EPSTOL
				#		println("HERE P -> l = $(l), violation = ",allretcumdem[1,l]-alpha)

					expr = 0
					for k in 1:l
						if SP[k]
							expr += yp[1,k]*allretcumdem[k,l]
						else
							expr += xp[1,k]
						end
					end
					@usercut(cb, expr >= allretcumdem[1,l])
					addedcuts += 1

				end

			end



			#W
			for w in 1:inst.NW

				for l in 1:inst.NT
					alpha = 0.0
					for k in 1:l
						if yw_vals[w,k]*warecumdem[w,k,l] < xw_vals[w,k]
							SW[w,k] = true
							alpha += yw_vals[w,k]*warecumdem[w,k,l]
						else
							SW[w,k] = false
							alpha +=  xw_vals[w,k]
						end
					end
					#println("R viol = ",inst.cumdem[r,1,l]-alpha)
					if alpha < warecumdem[w,1,l] - EPSTOL
	#						println("W r = $(r), l = $(l), violation = ",inst.cumdem[r,1,l]-alpha)

						expr = 0
						for k in 1:l
							if SW[w,k]
								expr += yw[w,k]*warecumdem[w,k,l]
							else
								expr += xw[w,k]
							end
						end
						@usercut(cb, expr >= warecumdem[w,1,l])
						addedcuts += 1
					end

				end
			end




			# R
			for r in 1:inst.NR

				for l in 1:inst.NT
					alpha = 0.0
					for k in 1:l
						if yr_vals[r,k]*inst.cumdem[r,k,l] < xr_vals[r,k]
							S[r,k] = true
							alpha += yr_vals[r,k]*inst.cumdem[r,k,l]
						else
							S[r,k] = false
							alpha +=  xr_vals[r,k]
						end
					end
					#println("R viol = ",inst.cumdem[r,1,l]-alpha)
					if alpha < inst.cumdem[r,1,l] - EPSTOL
	#						println("R r = $(r), l = $(l), violation = ",inst.cumdem[r,1,l]-alpha)

						expr = 0
						for k in 1:l
							if S[r,k]
								expr += yr[r,k]*inst.cumdem[r,k,l]
							else
								expr += xr[r,k]
							end
						end
						@usercut(cb, expr >= inst.cumdem[r,1,l])
						addedcuts += 1

					end

				end
			end


			# P+W
			for l in 2:inst.NT
				for lp in 1:l-1
					alpha = 0.0
					for k in 1:lp
						if yp_vals[1,k]*allretcumdem[k,l] < xp_vals[1,k]
							SP[k] = true
							alpha += yp_vals[1,k]*allretcumdem[k,l]
						else
							SP[k] = false
							alpha +=  xp_vals[1,k]
						end
					end
					#println("viol = ",inst.cumdem[r,1,l]-alpha)

					for k in lp+1:l
						for w in 1:inst.NW
							if yw_vals[w,k]*warecumdem[w,k,l] < xw_vals[w,k]
								SW[w,k] = true
								alpha += yw_vals[w,k]*warecumdem[w,k,l]
							else
								SW[w,k] = false
								alpha +=  xw_vals[w,k]
							end
						end
					end




					if alpha < allretcumdem[1,l] - EPSTOL
	#						println("r = $(r), l = $(l), violation = ",inst.cumdem[r,1,l]-alpha)

						expr = 0
						for k in 1:lp
							if SP[k]
								expr += yp[1,k]*allretcumdem[k,l]
							else
								expr += xp[1,k]
							end
						end

						for k in lp+1:l
							for w in 1:inst.NW
								if SW[w,k]
									expr += yw[w,k]*warecumdem[w,k,l]
								else
									expr += xw[w,k]
								end
							end
						end


						@usercut(cb, expr >= allretcumdem[1,l])
						addedcuts += 1

					end
				end
			end



			# P+R #apparently, these do not help at all
			for l in 2:inst.NT
				for lp in 1:l-1
					alpha = 0.0
					for k in 1:lp
						if yp_vals[1,k]*allretcumdem[k,l] < xp_vals[1,k]
							SP[k] = true
							alpha += yp_vals[1,k]*allretcumdem[k,l]
						else
							SP[k] = false
							alpha +=  xp_vals[1,k]
						end
					end
					#println("viol = ",inst.cumdem[r,1,l]-alpha)

					for k in lp+1:l
						for r in 1:inst.NR
							if yr_vals[r,k]*inst.cumdem[r,k,l] < xr_vals[r,k]
								S[r,k] = true
								alpha += yr_vals[r,k]*inst.cumdem[r,k,l]
							else
								S[r,k] = false
								alpha +=  xr_vals[r,k]
							end
						end
					end


					if alpha < allretcumdem[1,l] - EPSTOL
						#println("l = $(l), violation = ",allretcumdem[1,l]-alpha)

						expr = 0
						for k in 1:lp
							if SP[k]
								expr += yp[1,k]*allretcumdem[k,l]
							else
								expr += xp[1,k]
							end
						end

						for k in lp+1:l
							for r in 1:inst.NR
								if S[r,k]
									expr += yr[r,k]*inst.cumdem[r,k,l]
								else
									expr += xr[r,k]
								end
							end
						end


						@usercut(cb, expr >= allretcumdem[1,l])
						addedcuts += 1

					end
				end
			end








			#println("********************* W+R ***********************")

			for w in 1:inst.NW

				for l in 2:inst.NT

					for lw in 1:l-1
						alpha = 0.0
						for k in 1:lw
			#				println("r= $(r), lw= $(lw), l= $(l)  -> k = $(k)")
							if yw_vals[w,k]*warecumdem[w,k,l] < xw_vals[w,k]
								SW[w,k] = true
								alpha += yw_vals[w,k]*warecumdem[w,k,l]
							else
								SW[w,k] = false
								alpha +=  xw_vals[w,k]
							end
						end
						#println("R viol = ",inst.cumdem[r,1,l]-alpha)
						for ret in 1:length(inst.DeltaW[w])
							r = inst.DeltaW[w][ret]
							for k in lw+1:l
								if yr_vals[r,k]*inst.cumdem[r,k,l] < xr_vals[r,k]
									S[r,k] = true
									alpha += yr_vals[r,k]*inst.cumdem[r,k,l]
								else
									S[r,k] = false
									alpha +=  xr_vals[r,k]
								end
							end
						end
						if alpha < warecumdem[w,1,l] - EPSTOL
			#				println("WR r = $(r), l = $(l), violation = ",inst.cumdem[r,1,l]-alpha)

							expr = 0
							for k in 1:lw
								if SW[w,k]
									expr += yw[w,k]*warecumdem[w,k,l]
								else
									expr += xw[w,k]
								end
							end

							#println("R r = $(r), l = $(l), violation = ",inst.cumdem[r,1,l]-alpha)
							for ret in 1:length(inst.DeltaW[w])
								r = inst.DeltaW[w][ret]
								for k in lw+1:l
									if S[r,k]
										expr += yr[r,k]*inst.cumdem[r,k,l]
									else
										expr += xr[r,k]
									end
								end
							end
							@usercut(cb, expr >= warecumdem[w,1,l])
							addedcuts += 1

						end
					end
				end
			end





			#println("Separated cuts: $(addedcuts)")

		end

	end  # End of callback function



	function fastallretailerseparator(cb)
		node = MathProgBase.cbgetexplorednodes(cb)
		lk = zeros(Int,inst.NT)
		addedcuts = 0
		#println("Nodes = ",node)

		alpha = zeros(Float64,inst.NT)

		if node < 1
			xp_vals = getvalue(xp)
			yp_vals = getvalue(yp)
			xw_vals = getvalue(xw)
			yw_vals = getvalue(yw)
			xr_vals = getvalue(xr)
			yr_vals = getvalue(yr)

			EPSTOL = params.tolsep

			changel = Vector{Vector{Int}}(undef,inst.NT)
			for t in 1:inst.NT
				changel[t]=[]
			end
		#	println("\n\n\n !!!!!!!!!!!!!!!!!!!!!!!! Production site")
			for k in 1:inst.NT
				#println("\n\n k is ",k)
				#println("yk is ",yp_vals[1,k])
				#println("xk will be ",xp_vals[1,k])
				#println("will search in ",yp_vals[1,k]*allretcumdem[k,k:inst.NT])
				if yp_vals[1,k] > params.tolgap
					lk[k] = binarysearchfindlsinglelevel(k,xp_vals[1,k],yp_vals[1,k],allretcumdem,k,inst.NT)
					if lk[k] != -1
						push!(changel[lk[k]],k)
					end
				else
					lk[k] = -1
				end
				#println("lk[$(k)] is ",lk[k])
				#println("\n")

			end

			#println("changel = ",changel)

			sumY = 0.0
			sumX = 0.0
			bestviol = 0.0
			bestviolindex = 0
			if lk[1] == 1
				alpha[1] = xp_vals[1,1]
				sumX += xp_vals[1,1]
			else
				alpha[1] = yp_vals[1,1]*allretcumdem[1,1]
				sumY += yp_vals[1,1]
			end

			#println("L = ",1)
			#println("alpha[l] = $(alpha[1])")
			#println("violation = ",alpha[1] - allretcumdem[1,1])

			if alpha[1] - allretcumdem[1,1] < bestviol
				bestviol = alpha[1] - allretcumdem[1,1]
				bestviolindex = 1
			end

			for l in 2:inst.NT
			#	println("L = ",l)
				alpha[l] = alpha[l-1]

				sumY += yp_vals[1,l]
				if length(changel[l]) > 0
			#		println("length changel = $(length(changel[l]))")
					for perind in 1:length(changel[l])
						per = changel[l][perind]
			#			println("removing $(per) from Y_l")
						sumX += xp_vals[1,per]
						sumY -= yp_vals[1,per]
						alpha[l] += xp_vals[1,per]
						alpha[l] -= yp_vals[1,per]*allretcumdem[per,l-1]
					end
				end
				alpha[l] += sumY*allretcumdem[l,l]

				if alpha[l] - allretcumdem[1,l] < bestviol
					bestviol = alpha[l] - allretcumdem[1,l]
					bestviolindex = l
				end

			#	println("alpha[l] = $(alpha[l])")
			#	println("violation = ",alpha[l] - allretcumdem[1,l])


			end


			if bestviol < -EPSTOL
				l = bestviolindex
				expr = 0
				for k in 1:l
					if yp_vals[1,k]*allretcumdem[k,l] < xp_vals[1,k]
						expr += yp[1,k]*allretcumdem[k,l]
					else
						expr += xp[1,k]
					end
				end
				#@usercut(cb, expr >= allretcumdem[1,l])
				addedcuts += 1
			end

			#println("addedcuts = $(addedcuts)")
			#println("bestviol = $(bestviol)   ---->   period = $(bestviolindex)")
			for w in 1:inst.NW
#				println("\n\n\n ################ Warehouse : $(w)")
				for t in 1:inst.NT
					changel[t]=[]
				end
				for k in 1:inst.NT
					#println("\n k is ",k)
					#println("yk is ",yw_vals[w,k])
					#println("xk will be ",xw_vals[w,k])
					#println("will search in ",yw_vals[w,k]*warecumdem[w,k,k:inst.NT])
					if yw_vals[w,k] > params.tolgap
						lk[k] = binarysearchfindlsinglelevel(k,xw_vals[w,k],yw_vals[w,k],warecumdem[w,:,:],k,inst.NT)
						if lk[k] != -1
							push!(changel[lk[k]],k)
						end
					else
						lk[k] = -1
					end
					#println("lk[$(k)] is ",lk[k])
					#println("\n")
				end
				#println("changel = ",changel)


				sumY = 0.0
				sumX = 0.0
				bestviol = 0.0
				bestviolindex = 0
				if lk[1] == 1
					alpha[1] = xw_vals[w,1]
					sumX += xw_vals[w,1]
				else
					alpha[1] = yw_vals[w,1]*warecumdem[w,1,1]
					sumY += yw_vals[w,1]
				end

				#println("L = ",1)
				#println("alpha[l] = $(alpha[1])")
				#println("violation = ",alpha[1] - allretcumdem[1,1])

				if alpha[1] - warecumdem[w,1,1] < bestviol
					bestviol = alpha[1] - warecumdem[w,1,1]
					bestviolindex = 1
				end

				for l in 2:inst.NT
				#	println("L = ",l)
					alpha[l] = alpha[l-1]

					sumY += yw_vals[w,l]
					if length(changel[l]) > 0
				#		println("length changel = $(length(changel[l]))")
						for perind in 1:length(changel[l])
							per = changel[l][perind]
				#			println("removing $(per) from Y_l")
							sumX += xw_vals[w,per]
							sumY -= yw_vals[w,per]
							alpha[l] += xw_vals[w,per]
							alpha[l] -= yw_vals[w,per]*warecumdem[w,per,l-1]
						end
					end
					alpha[l] += sumY*warecumdem[w,l,l]

					if alpha[l] - warecumdem[w,1,l] < bestviol
						bestviol = alpha[l] - warecumdem[w,1,l]
						bestviolindex = l
					end

				#	println("alpha[l] = $(alpha[l])")
				#	println("violation = ",alpha[l] - allretcumdem[1,l])


				end

				if bestviol < -EPSTOL
					l = bestviolindex
					expr = 0
					for k in 1:l
						if yw_vals[w,k]*warecumdem[w,k,l] < xw_vals[w,k]
							expr += yw[w,k]*warecumdem[w,k,l]
						else
							expr += xw[w,k]
						end
					end
					#@usercut(cb, expr >= warecumdem[w,1,l])
					addedcuts += 1
				end


			end


			for r in 1:inst.NR
				#println("\n\n\n ***************** Retailer : $(r)")
				for t in 1:inst.NT
					changel[t]=[]
				end
				for k in 1:inst.NT
					#println("\n k is ",k)
					#println("yr is ",yr_vals[r,k])
					#println("xk will be ",xr_vals[r,k])
					#println("will search in ",yr_vals[r,k]*inst.cumdem[r,k,k:inst.NT])
					if yr_vals[r,k] > params.tolgap
						lk[k] = binarysearchfindlsinglelevel(k,xr_vals[r,k],yr_vals[r,k],inst.cumdem[r,:,:],k,inst.NT)
						if lk[k] != -1
							push!(changel[lk[k]],k)
						end
					else
						lk[k] = -1
					end
					#println("lk[$(k)] is ",lk[k])
					#println("\n")
				end
				#println("changel = ",changel)


				sumY = 0.0
				sumX = 0.0
				bestviol = 0.0
				bestviolindex = 0
				if lk[1] == 1
					alpha[1] = xr_vals[r,1]
					sumX += xr_vals[r,1]
				else
					alpha[1] = yr_vals[r,1]*inst.cumdem[r,1,1]
					sumY += yr_vals[r,1]
				end

				#println("L = ",1)
				#println("alpha[l] = $(alpha[1])")
				#println("violation = ",alpha[1] - allretcumdem[1,1])

				if alpha[1] - inst.cumdem[r,1,1] < bestviol
					bestviol = alpha[1] - inst.cumdem[r,1,1]
					bestviolindex = 1
				end

				for l in 2:inst.NT
				#	println("L = ",l)
					alpha[l] = alpha[l-1]

					sumY += yr_vals[r,l]
					if length(changel[l]) > 0
				#		println("length changel = $(length(changel[l]))")
						for perind in 1:length(changel[l])
							per = changel[l][perind]
				#			println("removing $(per) from Y_l")
							sumX += xr_vals[r,per]
							sumY -= yr_vals[r,per]
							alpha[l] += xr_vals[r,per]
							alpha[l] -= yr_vals[r,per]*inst.cumdem[r,per,l-1]
						end
					end
					alpha[l] += sumY*inst.cumdem[r,l,l]

					if alpha[l] - inst.cumdem[r,1,l] < bestviol
						bestviol = alpha[l] - inst.cumdem[r,1,l]
						bestviolindex = l
					end

				#	println("alpha[l] = $(alpha[l])")
				#	println("violation = ",alpha[l] - allretcumdem[1,l])


				end

				#println("bestviol = $(bestviol)")
				if bestviol < -EPSTOL
					l = bestviolindex
					expr = 0
					for k in 1:l
						if yr_vals[r,k]*inst.cumdem[r,k,l] < xr_vals[r,k]
							expr += yr[r,k]*inst.cumdem[r,k,l]
						else
							expr += xr[r,k]
						end
					end
					@usercut(cb, expr >= inst.cumdem[r,1,l])
					addedcuts += 1
				end






			end

			#println("Array lk = ",lk)

		end

	end  # End of callback function






	cuttingplanerounds = 0
	#addcutcallback(model,slowsingleretailerseparator)
	#addcutcallback(model,slowallretailerseparator)

	addcutcallback(model,slowallretailerseparatornew)

	#addcutcallback(model,fastallretailerseparator)

	t1 = time_ns()
	status = solve(model)
	t2 = time_ns()
	elapsedtime = (t2-t1)/1.0e9

	bestsol = getobjectivevalue(model)
	bestbound = getobjbound(model)
	gap = 100*(bestsol-bestbound)/bestsol
	numnodes = getnodecount(model)
	time = getsolvetime(model)
	#gap = getobjgap(model)
	println("bestsol = ", bestsol)
	println("bestbound = ", bestbound)
    println("gap = ", gap)
    println("time = ", time)
    println("nodes = ", numnodes)

	open("saida.txt","a") do f
		write(f,";$(params.form);$bestbound;$bestsol;$gap;$time;$numnodes;$(params.disablesolver) \n")
	end

#	if params.printsol == 1
#		printStandardFormulationSolution(inst,xp,xr,xw,yp,yr,yw,sp,sr,sw)
#	end









end #function bcstandardFormulation()



function bcthreelevelFormulation(inst::InstanceData, params::ParameterData)
	println("Running NewFormulations.threelevelFormulation")

	if params.solver == "Gurobi"
		if params.disablesolver == 1 #Disable gurobi cuts and presolve
			if params.maxnodes < 999.0
				model = Model(solver = GurobiSolver(TimeLimit=params.maxtime,MIPGap=params.tolgap,CliqueCuts=0, CoverCuts=0, FlowCoverCuts=0,FlowPathCuts=0,MIRCuts=0,NetworkCuts=0,GomoryPasses=0, PreCrush=1,NodeLimit=params.maxnodes))
			else
				model = Model(solver = GurobiSolver(TimeLimit=params.maxtime,MIPGap=params.tolgap,CliqueCuts=0, CoverCuts=0, FlowCoverCuts=0,FlowPathCuts=0,MIRCuts=0,NetworkCuts=0,GomoryPasses=0, PreCrush=1))
			end

			#model = Model(solver = GurobiSolver(TimeLimit=params.maxtime,MIPGap=params.tolgap,CliqueCuts=0, CoverCuts=0, FlowCoverCuts=0,FlowPathCuts=0,MIRCuts=0,NetworkCuts=0,GomoryPasses=0, Presolve=0,NodeLimit=1,CutPasses=5,OutputFlag=0))
			#model = Model(solver = GurobiSolver(TimeLimit=params.maxtime,MIPGap=params.tolgap,CliqueCuts=0, CoverCuts=0, FlowCoverCuts=0,FlowPathCuts=0,MIRCuts=0,NetworkCuts=0,GomoryPasses=1,ImpliedCuts=0, Presolve=0))
			#model = Model(solver = GurobiSolver(TimeLimit=params.maxtime,MIPGap=params.tolgap,CliqueCuts=0, CoverCuts=0, FlowCoverCuts=0,FlowPathCuts=0,MIRCuts=0,NetworkCuts=0,GomoryPasses=0,ImpliedCuts=0, Presolve=0))
			#model = Model(solver = GurobiSolver(MIPGap=params.tolgap,Cuts=0, Presolve=0))
		else
			if params.maxnodes < 999.0
        		model = Model(solver = GurobiSolver(TimeLimit=params.maxtime,MIPGap=params.tolgap,PreCrush=1,NodeLimit=params.maxnodes))
			else
				model = Model(solver = GurobiSolver(TimeLimit=params.maxtime,MIPGap=params.tolgap,PreCrush=1))
			end
		end
	elseif params.solver == "Cplex"
		model = Model(solver = CplexSolver(CPX_PARAM_TILIM=params.maxtime,CPX_PARAM_EPGAP=params.tolgap,CPX_PARAM_PREIND=0))
	else
		println("No solver selected")
		return 0
	end

	### Defining variables ###
	@variable(model,0 <= x0[i=1:inst.NR,t=1:inst.NT] <= Inf) #sum(inst.DP[i,t] for t=1:inst.NT))
	@variable(model,0 <= x2[i=1:inst.NR,t=1:inst.NT] <= Inf) #sum(inst.DR[i,t] for t=1:inst.NT))
	@variable(model,0 <= x1[i=1:inst.NR,t=1:inst.NT] <= Inf) #sum(inst.D[i,t] for t=1:inst.NT))
	@variable(model, yp[i=1:inst.NP,t=1:inst.NT], Bin)
	@variable(model, yr[i=1:inst.NR,t=1:inst.NT], Bin)
	@variable(model, yw[i=1:inst.NW,t=1:inst.NT], Bin)
	@variable(model,0 <= s0[i=1:inst.NR,t=0:inst.NT] <= Inf) #sum(inst.DP[i,t] for t=1:inst.NT))
	@variable(model,0 <= s2[i=1:inst.NR,t=0:inst.NT] <= Inf) #sum(inst.DR[i,t] for t=1:inst.NT))
	@variable(model,0 <= s1[i=1:inst.NR,t=0:inst.NT] <= Inf)

cumdem = zeros(Int,inst.NR,inst.NT,inst.NT)

for r in 1:inst.NT, t in 1:inst.NT
	cumdem[r,t,t] = inst.D[r,t]
	for k in t+1:inst.NT
		cumdem[r,t,k] = cumdem[r,t,k-1] + inst.D[r,k]
	end
end

println("Defined variables")
	### Objective function ###
	@objective(model, Min,
		sum(inst.SCP[i,t]*yp[i,t]  for i=1:inst.NP, t=1:inst.NT)
		+ sum(inst.SCR[i,t]*yr[i,t] for i=1:inst.NR, t=1:inst.NT)
		+ sum(inst.SCW[i,t]*yw[i,t] for i=1:inst.NW, t=1:inst.NT)
		+ sum( inst.HCP[1]*s0[r,t] for r=1:inst.NR, t=1:inst.NT)
		+ sum( inst.HCR[r]*s2[r,t] for r=1:inst.NR, t=1:inst.NT)
		+ sum( inst.HCW[inst.DeltamR[r]]*s1[r,t] for r=1:inst.NR, t=1:inst.NT)
	)

println("Finished objective function")

	### Setup constraints ###
	@constraint(model,
				setupP[r=1:inst.NR,t=1:inst.NT], x0[r,t] <= sum(inst.D[r,k] for k in t:inst.NT)*yp[1,t]
				)
	@constraint(model,
				setupW[r=1:inst.NR, t=1:inst.NT], x1[r,t] <= sum(inst.D[r,k] for k in t:inst.NT)*yw[inst.DeltamR[r],t]
				)
	@constraint(model,
				setupR[r=1:inst.NR, t=1:inst.NT], x2[r,t] <= sum(inst.D[r,k] for k in t:inst.NT)*yr[r,t]
				)

	if params.capacity != 0
		@constraint(model,
					capP[t=1:inst.NT], sum(x0[r,t] for r in 1:inst.NR) <= inst.C[t]*yp[1,t]
					)
	end

	println("Finished setup constraints")

	### no initial inventory
	@constraint(model,
				zeroinvP[r=1:inst.NR],
				s0[r,0] == 0
				)
	@constraint(model,
				zeroinvW[r=1:inst.NR],
				s1[r,0] == 0
				)
	@constraint(model,
				zeroinvR[r=1:inst.NR],
				s2[r,0] == 0
				)

println("Finished zero inventory constraints")

	### Balance constraints ###
	@constraint(model,
				balanceP[r=1:inst.NR, t=1:inst.NT],
				s0[r,t-1] + x0[r,t] == x1[r,t] + s0[r,t]
				)
	@constraint(model,
				balanceW[r=1:inst.NR,t=1:inst.NT],
				s1[r,t-1] + x1[r,t] == x2[r,t] + s1[r,t]
				)
	@constraint(model,
				balanceR[r=1:inst.NR,t=1:inst.NT],
				s2[r,t-1] + x2[r,t] == inst.D[r,t] + s2[r,t]
				)

	if params.form == "randdpheurbc3level"
		#println("HERE")
		SETR,SETW,SETP = DPHeuristics.RandomizedDPHeuristicBottomUp(inst,params)

		ypsol = zeros(inst.NP,inst.NT)
		ywsol = zeros(inst.NW,inst.NT)
		yrsol = zeros(inst.NR,inst.NT)

		for t in 1:length(SETP)
			ypsol[1,SETP[t]] = 1
		end
		for w in 1:inst.NW, t in 1:length(SETW[w])
			ywsol[w,SETW[w][t]] = 1
		end
		for r in 1:inst.NR, t in 1:length(SETR[r])
			yrsol[r,SETR[r][t]] = 1
		end

		for t = 1:inst.NT
			JuMP.setvalue(yp[1,t],ypsol[1,t])
		end
		for w=1:inst.NW, t=1:inst.NT
			JuMP.setvalue(yw[w,t], ywsol[w,t])
		end
		for r=1:inst.NR, t=1:inst.NT
			JuMP.setvalue(yr[r,t],yrsol[r,t])
		end
	end

	S = zeros(Bool,inst.NR,inst.NT)
	function slowsingleretailerseparator3level(cb)
		node = MathProgBase.cbgetexplorednodes(cb)

		addedcuts = 0
		#println("Nodes = ",node)

		if node < 1
			x0_vals = getvalue(x0)
			yp_vals = getvalue(yp)
			x1_vals = getvalue(x1)
			yw_vals = getvalue(yw)
			x2_vals = getvalue(x2)
			yr_vals = getvalue(yr)

			EPSTOL = params.tolsep

			#println("********************* P ***********************")
			#println("length lw_vals = ",length(yw_vals))
			for r in 1:inst.NR

				for l in 1:inst.NT
					alpha = 0.0
					for k in 1:l
						if yp_vals[1,k]*inst.cumdem[r,k,l] < x0_vals[r,k]
							S[r,k] = true
							alpha += yp_vals[1,k]*inst.cumdem[r,k,l]
						else
							S[r,k] = false
							alpha +=  x0_vals[r,k]
						end
					end
					#println("viol = ",inst.cumdem[r,1,l]-alpha)
					if alpha < inst.cumdem[r,1,l] - EPSTOL
						#println("r = $(r), l = $(l), violation = ",inst.cumdem[r,1,l]-alpha)

						expr = 0
						for k in 1:l
							if S[r,k]
								expr += yp[1,k]*inst.cumdem[r,k,l]
							else
								expr += x0[r,k]
							end
						end
						@usercut(cb, expr >= inst.cumdem[r,1,l])
						addedcuts += 1

					end

				end
			end

			#println("********************* W ***********************")

			for r in 1:inst.NR

				for l in 1:inst.NT
					alpha = 0.0
					for k in 1:l
						if yw_vals[inst.DeltamR[r],k]*inst.cumdem[r,k,l] < x1_vals[r,k]
							S[r,k] = true
							alpha += yw_vals[inst.DeltamR[r],k]*inst.cumdem[r,k,l]
						else
							S[r,k] = false
							alpha +=  x1_vals[r,k]
						end
					end
					#println("R viol = ",inst.cumdem[r,1,l]-alpha)
					if alpha < inst.cumdem[r,1,l] - EPSTOL
						#println("W r = $(r), l = $(l), violation = ",inst.cumdem[r,1,l]-alpha)

						expr = 0
						for k in 1:l
							if S[r,k]
								expr += yw[inst.DeltamR[r],k]*inst.cumdem[r,k,l]
							else
								expr += x1[r,k]
							end
						end
						@usercut(cb, expr >= inst.cumdem[r,1,l])
						addedcuts += 1
					end

				end
			end

			#println("********************* R ***********************")

			for r in 1:inst.NR

				for l in 1:inst.NT
					alpha = 0.0
					for k in 1:l
						if yr_vals[r,k]*inst.cumdem[r,k,l] < x2_vals[r,k]
							S[r,k] = true
							alpha += yr_vals[r,k]*inst.cumdem[r,k,l]
						else
							S[r,k] = false
							alpha +=  x2_vals[r,k]
						end
					end
					#println("R viol = ",inst.cumdem[r,1,l]-alpha)
					if alpha < inst.cumdem[r,1,l] - EPSTOL
						#println("R r = $(r), l = $(l), violation = ",inst.cumdem[r,1,l]-alpha)

						expr = 0
						for k in 1:l
							if S[r,k]
								expr += yr[r,k]*inst.cumdem[r,k,l]
							else
								expr += x2[r,k]
							end
						end
						@usercut(cb, expr >= inst.cumdem[r,1,l])
						addedcuts += 1

					end

				end
			end


			#println("********************* W+R ***********************")

			for r in 1:inst.NR

				for l in 2:inst.NT

					for lw in 1:l-1
						alpha = 0.0
						for k in 1:lw
			#				println("r= $(r), lw= $(lw), l= $(l)  -> k = $(k)")
							if yw_vals[inst.DeltamR[r],k]*inst.cumdem[r,k,l] < x1_vals[r,k]
								S[r,k] = true
								alpha += yw_vals[inst.DeltamR[r],k]*inst.cumdem[r,k,l]
							else
								S[r,k] = false
								alpha +=  x1_vals[r,k]
							end
						end
						#println("R viol = ",inst.cumdem[r,1,l]-alpha)

						for k in lw+1:l
							if yr_vals[r,k]*inst.cumdem[r,k,l] < x2_vals[r,k]
								S[r,k] = true
								alpha += yr_vals[r,k]*inst.cumdem[r,k,l]
							else
								S[r,k] = false
								alpha +=  x2_vals[r,k]
							end
						end

						if alpha < inst.cumdem[r,1,l] - EPSTOL
							#println("WR r = $(r), l = $(l), violation = ",inst.cumdem[r,1,l]-alpha)

							expr = 0
							for k in 1:lw
								if S[r,k]
									expr += yw[inst.DeltamR[r],k]*inst.cumdem[r,k,l]
								else
									expr += x1[r,k]
								end
							end

							#println("R r = $(r), l = $(l), violation = ",inst.cumdem[r,1,l]-alpha)

							for k in lw+1:l
								if S[r,k]
									expr += yr[r,k]*inst.cumdem[r,k,l]
								else
									expr += x2[r,k]
								end
							end
							@usercut(cb, expr >= inst.cumdem[r,1,l])
							addedcuts += 1

						end
					end
				end
			end


			#println("********************* P+W ***********************")

			for r in 1:inst.NR

				for l in 2:inst.NT

					for lp in 1:l-1
						alpha = 0.0

						for k in 1:lp
							if yp_vals[1,k]*inst.cumdem[r,k,l] < x0_vals[r,k]
								S[r,k] = true
								alpha += yp_vals[1,k]*inst.cumdem[r,k,l]
							else
								S[r,k] = false
								alpha +=  x0_vals[r,k]
							end
						end

						for k in lp+1:l
							if yw_vals[inst.DeltamR[r],k]*inst.cumdem[r,k,l] < x1_vals[r,k]
								S[r,k] = true
								alpha += yw_vals[inst.DeltamR[r],k]*inst.cumdem[r,k,l]
							else
								S[r,k] = false
								alpha +=  x1_vals[r,k]
							end
						end
						#println("R viol = ",inst.cumdem[r,1,l]-alpha)


						if alpha < inst.cumdem[r,1,l] - EPSTOL
							#println("PW r = $(r), l = $(l), violation = ",inst.cumdem[r,1,l]-alpha)

							expr = 0
							for k in 1:lp
								if S[r,k]
									expr += yp[1,k]*inst.cumdem[r,k,l]
								else
									expr += x0[r,k]
								end
							end

							#println("R r = $(r), l = $(l), violation = ",inst.cumdem[r,1,l]-alpha)

							for k in lp+1:l
								if S[r,k]
									expr += yw[inst.DeltamR[r],k]*inst.cumdem[r,k,l]
								else
									expr += x1[r,k]
								end
							end
							@usercut(cb, expr >= inst.cumdem[r,1,l])
							addedcuts += 1

						end
					end
				end
			end


			#println("********************* P+R ***********************")

			for r in 1:inst.NR

				for l in 2:inst.NT

					for lp in 1:l-1
						alpha = 0.0

						for k in 1:lp
							if yp_vals[1,k]*inst.cumdem[r,k,l] < x0_vals[r,k]
								S[r,k] = true
								alpha += yp_vals[1,k]*inst.cumdem[r,k,l]
							else
								S[r,k] = false
								alpha +=  x0_vals[r,k]
							end
						end

						for k in lp+1:l
							if yr_vals[r,k]*inst.cumdem[r,k,l] < x2_vals[r,k]
								S[r,k] = true
								alpha += yr_vals[r,k]*inst.cumdem[r,k,l]
							else
								S[r,k] = false
								alpha +=  x2_vals[r,k]
							end
						end
						#println("R viol = ",inst.cumdem[r,1,l]-alpha)


						if alpha < inst.cumdem[r,1,l] - EPSTOL
							#println("PR r = $(r), l = $(l), violation = ",inst.cumdem[r,1,l]-alpha)

							expr = 0
							for k in 1:lp
								if S[r,k]
									expr += yp[1,k]*inst.cumdem[r,k,l]
								else
									expr += x0[r,k]
								end
							end

							#println("R r = $(r), l = $(l), violation = ",inst.cumdem[r,1,l]-alpha)

							for k in lp+1:l
								if S[r,k]
									expr += yr[r,k]*inst.cumdem[r,k,l]
								else
									expr += x2[r,k]
								end
							end
							@usercut(cb, expr >= inst.cumdem[r,1,l])
							addedcuts += 1

						end
					end
				end
			end


			#println("********************* P+W+R ***********************")

			for r in 1:inst.NR
				for l in 3:inst.NT
					for lp in 1:l-2
						for lw in lp+1:l-1
							alpha = 0.0

							for k in 1:lp
								if yp_vals[1,k]*inst.cumdem[r,k,l] < x0_vals[r,k]
									S[r,k] = true
									alpha += yp_vals[1,k]*inst.cumdem[r,k,l]
								else
									S[r,k] = false
									alpha +=  x0_vals[r,k]
								end
							end

							for k in lp+1:lw
								if yw_vals[inst.DeltamR[r],k]*inst.cumdem[r,k,l] < x1_vals[r,k]
									S[r,k] = true
									alpha += yw_vals[inst.DeltamR[r],k]*inst.cumdem[r,k,l]
								else
									S[r,k] = false
									alpha +=  x1_vals[r,k]
								end
							end

							for k in lw+1:l
								if yr_vals[r,k]*inst.cumdem[r,k,l] < x2_vals[r,k]
									S[r,k] = true
									alpha += yr_vals[r,k]*inst.cumdem[r,k,l]
								else
									S[r,k] = false
									alpha +=  x2_vals[r,k]
								end
							end
							#println("R viol = ",inst.cumdem[r,1,l]-alpha)


							if alpha < inst.cumdem[r,1,l] - EPSTOL
								#println("PWR r = $(r), l = $(l), violation = ",inst.cumdem[r,1,l]-alpha)

								expr = 0
								for k in 1:lp
									if S[r,k]
										expr += yp[1,k]*inst.cumdem[r,k,l]
									else
										expr += x0[r,k]
									end
								end

								#println("R r = $(r), l = $(l), violation = ",inst.cumdem[r,1,l]-alpha)

								for k in lp+1:l
									if S[r,k]
										expr += yw[inst.DeltamR[r],k]*inst.cumdem[r,k,l]
									else
										expr += x1[r,k]
									end
								end

								for k in lw+1:l
									if S[r,k]
										expr += yr[r,k]*inst.cumdem[r,k,l]
									else
										expr += x2[r,k]
									end
								end

								@usercut(cb, expr >= inst.cumdem[r,1,l])
								addedcuts += 1

							end
						end
					end
				end
			end

			#println("Separated cuts: $(addedcuts)")
		end




	end  # End of callback function


	S = zeros(Bool,inst.NR,inst.NT)
	function slowsingleretailerseparator3levelnew(cb)
		node = MathProgBase.cbgetexplorednodes(cb)

		addedcuts = 0
		#println("Nodes = ",node)

		if node < 1 && cuttingplanerounds <= params.cpmaxrounds

			cuttingplanerounds += 1

			#println("cuttingplanerounds = $(cuttingplanerounds), bestbound = $(MathProgBase.cbgetbestbound(cb))")

			x0_vals = getvalue(x0)
			yp_vals = getvalue(yp)
			x1_vals = getvalue(x1)
			yw_vals = getvalue(yw)
			x2_vals = getvalue(x2)
			yr_vals = getvalue(yr)

			EPSTOL = params.tolsep

			if true
				#println("********************* P ***********************")
				#println("length lw_vals = ",length(yw_vals))
				for r in 1:inst.NR

					for l in 1:inst.NT
						alpha = 0.0
						for k in 1:l
							if yp_vals[1,k]*inst.cumdem[r,k,l] < x0_vals[r,k]
								S[r,k] = true
								alpha += yp_vals[1,k]*inst.cumdem[r,k,l]
							else
								S[r,k] = false
								alpha +=  x0_vals[r,k]
							end
						end
						#println("viol = ",inst.cumdem[r,1,l]-alpha)
						if alpha < inst.cumdem[r,1,l] - EPSTOL
							#println("r = $(r), l = $(l), violation = ",inst.cumdem[r,1,l]-alpha)

							expr = 0
							for k in 1:l
								if S[r,k]
									expr += yp[1,k]*inst.cumdem[r,k,l]
								else
									expr += x0[r,k]
								end
							end
							@usercut(cb, expr >= inst.cumdem[r,1,l])
							addedcuts += 1

						end

					end
				end

				#println("********************* W ***********************")

				for r in 1:inst.NR

					for l in 1:inst.NT
						alpha = 0.0
						for k in 1:l
							if yw_vals[inst.DeltamR[r],k]*inst.cumdem[r,k,l] < x1_vals[r,k]
								S[r,k] = true
								alpha += yw_vals[inst.DeltamR[r],k]*inst.cumdem[r,k,l]
							else
								S[r,k] = false
								alpha +=  x1_vals[r,k]
							end
						end
						#println("R viol = ",inst.cumdem[r,1,l]-alpha)
						if alpha < inst.cumdem[r,1,l] - EPSTOL
							#println("W r = $(r), l = $(l), violation = ",inst.cumdem[r,1,l]-alpha)

							expr = 0
							for k in 1:l
								if S[r,k]
									expr += yw[inst.DeltamR[r],k]*inst.cumdem[r,k,l]
								else
									expr += x1[r,k]
								end
							end
							@usercut(cb, expr >= inst.cumdem[r,1,l])
							addedcuts += 1
						end

					end
				end

				#println("********************* R ***********************")

				for r in 1:inst.NR

					for l in 1:inst.NT
						alpha = 0.0
						for k in 1:l
							if yr_vals[r,k]*inst.cumdem[r,k,l] < x2_vals[r,k]
								S[r,k] = true
								alpha += yr_vals[r,k]*inst.cumdem[r,k,l]
							else
								S[r,k] = false
								alpha +=  x2_vals[r,k]
							end
						end
						#println("R viol = ",inst.cumdem[r,1,l]-alpha)
						if alpha < inst.cumdem[r,1,l] - EPSTOL
							#println("R r = $(r), l = $(l), violation = ",inst.cumdem[r,1,l]-alpha)

							expr = 0
							for k in 1:l
								if S[r,k]
									expr += yr[r,k]*inst.cumdem[r,k,l]
								else
									expr += x2[r,k]
								end
							end
							@usercut(cb, expr >= inst.cumdem[r,1,l])
							addedcuts += 1

						end

					end
				end
			end


			#return
			if mod(cuttingplanerounds,5)== 0
				#println("********************* W+R ***********************")

				for r in 1:inst.NR

					for l in 2:inst.NT

						for lw in 1:l-1
							alpha = 0.0
							for k in 1:lw
				#				println("r= $(r), lw= $(lw), l= $(l)  -> k = $(k)")
								if yw_vals[inst.DeltamR[r],k]*inst.cumdem[r,k,l] < x1_vals[r,k]
									S[r,k] = true
									alpha += yw_vals[inst.DeltamR[r],k]*inst.cumdem[r,k,l]
								else
									S[r,k] = false
									alpha +=  x1_vals[r,k]
								end
							end
							#println("R viol = ",inst.cumdem[r,1,l]-alpha)

							for k in lw+1:l
								if yr_vals[r,k]*inst.cumdem[r,k,l] < x2_vals[r,k]
									S[r,k] = true
									alpha += yr_vals[r,k]*inst.cumdem[r,k,l]
								else
									S[r,k] = false
									alpha +=  x2_vals[r,k]
								end
							end

							if alpha < inst.cumdem[r,1,l] - EPSTOL
								#println("WR r = $(r), l = $(l), violation = ",inst.cumdem[r,1,l]-alpha)

								expr = 0
								for k in 1:lw
									if S[r,k]
										expr += yw[inst.DeltamR[r],k]*inst.cumdem[r,k,l]
									else
										expr += x1[r,k]
									end
								end

								#println("R r = $(r), l = $(l), violation = ",inst.cumdem[r,1,l]-alpha)

								for k in lw+1:l
									if S[r,k]
										expr += yr[r,k]*inst.cumdem[r,k,l]
									else
										expr += x2[r,k]
									end
								end
								@usercut(cb, expr >= inst.cumdem[r,1,l])
								addedcuts += 1

							end
						end
					end
				end


				#println("********************* P+W ***********************")

				for r in 1:inst.NR

					for l in 2:inst.NT

						for lp in 1:l-1
							alpha = 0.0

							for k in 1:lp
								if yp_vals[1,k]*inst.cumdem[r,k,l] < x0_vals[r,k]
									S[r,k] = true
									alpha += yp_vals[1,k]*inst.cumdem[r,k,l]
								else
									S[r,k] = false
									alpha +=  x0_vals[r,k]
								end
							end

							for k in lp+1:l
								if yw_vals[inst.DeltamR[r],k]*inst.cumdem[r,k,l] < x1_vals[r,k]
									S[r,k] = true
									alpha += yw_vals[inst.DeltamR[r],k]*inst.cumdem[r,k,l]
								else
									S[r,k] = false
									alpha +=  x1_vals[r,k]
								end
							end
							#println("R viol = ",inst.cumdem[r,1,l]-alpha)


							if alpha < inst.cumdem[r,1,l] - EPSTOL
								#println("PW r = $(r), l = $(l), violation = ",inst.cumdem[r,1,l]-alpha)

								expr = 0
								for k in 1:lp
									if S[r,k]
										expr += yp[1,k]*inst.cumdem[r,k,l]
									else
										expr += x0[r,k]
									end
								end

								#println("R r = $(r), l = $(l), violation = ",inst.cumdem[r,1,l]-alpha)

								for k in lp+1:l
									if S[r,k]
										expr += yw[inst.DeltamR[r],k]*inst.cumdem[r,k,l]
									else
										expr += x1[r,k]
									end
								end
								@usercut(cb, expr >= inst.cumdem[r,1,l])
								addedcuts += 1

							end
						end
					end
				end


				#println("********************* P+R ***********************")

				for r in 1:inst.NR

					for l in 2:inst.NT

						for lp in 1:l-1
							alpha = 0.0

							for k in 1:lp
								if yp_vals[1,k]*inst.cumdem[r,k,l] < x0_vals[r,k]
									S[r,k] = true
									alpha += yp_vals[1,k]*inst.cumdem[r,k,l]
								else
									S[r,k] = false
									alpha +=  x0_vals[r,k]
								end
							end

							for k in lp+1:l
								if yr_vals[r,k]*inst.cumdem[r,k,l] < x2_vals[r,k]
									S[r,k] = true
									alpha += yr_vals[r,k]*inst.cumdem[r,k,l]
								else
									S[r,k] = false
									alpha +=  x2_vals[r,k]
								end
							end
							#println("R viol = ",inst.cumdem[r,1,l]-alpha)


							if alpha < inst.cumdem[r,1,l] - EPSTOL
								#println("PR r = $(r), l = $(l), violation = ",inst.cumdem[r,1,l]-alpha)

								expr = 0
								for k in 1:lp
									if S[r,k]
										expr += yp[1,k]*inst.cumdem[r,k,l]
									else
										expr += x0[r,k]
									end
								end

								#println("R r = $(r), l = $(l), violation = ",inst.cumdem[r,1,l]-alpha)

								for k in lp+1:l
									if S[r,k]
										expr += yr[r,k]*inst.cumdem[r,k,l]
									else
										expr += x2[r,k]
									end
								end
								@usercut(cb, expr >= inst.cumdem[r,1,l])
								addedcuts += 1

							end
						end
					end
				end
			end



			if mod(cuttingplanerounds,10)== 0
				#println("********************* P+W+R ***********************")

				for r in 1:inst.NR
					for l in 3:inst.NT
						for lp in 1:l-2
							for lw in lp+1:l-1
								alpha = 0.0

								for k in 1:lp
									if yp_vals[1,k]*inst.cumdem[r,k,l] < x0_vals[r,k]
										S[r,k] = true
										alpha += yp_vals[1,k]*inst.cumdem[r,k,l]
									else
										S[r,k] = false
										alpha +=  x0_vals[r,k]
									end
								end

								for k in lp+1:lw
									if yw_vals[inst.DeltamR[r],k]*inst.cumdem[r,k,l] < x1_vals[r,k]
										S[r,k] = true
										alpha += yw_vals[inst.DeltamR[r],k]*inst.cumdem[r,k,l]
									else
										S[r,k] = false
										alpha +=  x1_vals[r,k]
									end
								end

								for k in lw+1:l
									if yr_vals[r,k]*inst.cumdem[r,k,l] < x2_vals[r,k]
										S[r,k] = true
										alpha += yr_vals[r,k]*inst.cumdem[r,k,l]
									else
										S[r,k] = false
										alpha +=  x2_vals[r,k]
									end
								end
								#println("R viol = ",inst.cumdem[r,1,l]-alpha)


								if alpha < inst.cumdem[r,1,l] - EPSTOL
									#println("PWR r = $(r), l = $(l), violation = ",inst.cumdem[r,1,l]-alpha)

									expr = 0
									for k in 1:lp
										if S[r,k]
											expr += yp[1,k]*inst.cumdem[r,k,l]
										else
											expr += x0[r,k]
										end
									end

									#println("R r = $(r), l = $(l), violation = ",inst.cumdem[r,1,l]-alpha)

									for k in lp+1:l
										if S[r,k]
											expr += yw[inst.DeltamR[r],k]*inst.cumdem[r,k,l]
										else
											expr += x1[r,k]
										end
									end

									for k in lw+1:l
										if S[r,k]
											expr += yr[r,k]*inst.cumdem[r,k,l]
										else
											expr += x2[r,k]
										end
									end

									@usercut(cb, expr >= inst.cumdem[r,1,l])
									addedcuts += 1

								end
							end
						end
					end
				end
			end
			#println("Separated cuts: $(addedcuts)")
		end




	end  # End of callback function slowsingleretailerseparator3levelnew





	cuttingplanerounds = 0
	addcutcallback(model,slowsingleretailerseparator3levelnew)
	#addcutcallback(model,slowsingleretailerseparator3level)


	t1 = time_ns()
	status = solve(model)
	t2 = time_ns()
	elapsedtime = (t2-t1)/1.0e9

	bestsol = getobjectivevalue(model)
	bestbound = getobjbound(model)
	gap = 100*(bestsol-bestbound)/bestsol
	numnodes = getnodecount(model)
	time = getsolvetime(model)
	#gap = getobjgap(model)
	println("bestsol = ", bestsol)
	println("bestbound = ", bestbound)
    println("gap = ", gap)
    println("time = ", time)
    println("nodes = ", numnodes)

	open("saida.txt","a") do f
		write(f,";$(params.form);$bestbound;$bestsol;$gap;$time;$numnodes;$(params.disablesolver) \n")
	end

#	if params.printsol == 1
#		printthreelevelFormulationSolution(inst,xp,xr,xw,yp,yr,yw,sp,sr,sw)
#	end

end #function bcthreelevelFormulation()



function binarysearchfindlsinglelevel(k,xk,yk,cumdem,st,en)
	#println(st," ",en)
	if st>en
		return -1
	elseif st<en

		med = trunc(Int,ceil((st+en)/2))
		if yk*cumdem[k,med-1] <  xk && xk <= yk*cumdem[k,med]
			return med
		elseif yk*cumdem[k,med] < xk
			binarysearchfindlsinglelevel(k,xk,yk,cumdem,med+1,en)
		else
			binarysearchfindlsinglelevel(k,xk,yk,cumdem,st,med-1)
		end

	else
		return st
	end

end #function binarysearchfindlsinglelevel









function cuttingplanestandardFormulation(inst::InstanceData, params::ParameterData)

    println("Running BranchAndCut.cuttingplanestandardFormulation")

	if params.solver == "Gurobi"
		if params.disablesolver == 1 #Disable gurobi cuts and presolve
			model = Model(solver = GurobiSolver(TimeLimit=params.maxtime,MIPGap=params.tolgap,CliqueCuts=0, CoverCuts=0, FlowCoverCuts=0,FlowPathCuts=0,MIRCuts=0,NetworkCuts=0,GomoryPasses=0, PreCrush=1))
			#model = Model(solver = GurobiSolver(TimeLimit=params.maxtime,MIPGap=params.tolgap,CliqueCuts=0, CoverCuts=0, FlowCoverCuts=0,FlowPathCuts=0,MIRCuts=0,NetworkCuts=0,GomoryPasses=0, Presolve=0,NodeLimit=1,CutPasses=5,OutputFlag=0))
			#model = Model(solver = GurobiSolver(TimeLimit=params.maxtime,MIPGap=params.tolgap,CliqueCuts=0, CoverCuts=0, FlowCoverCuts=0,FlowPathCuts=0,MIRCuts=0,NetworkCuts=0,GomoryPasses=1,ImpliedCuts=0, Presolve=0))
			#model = Model(solver = GurobiSolver(TimeLimit=params.maxtime,MIPGap=params.tolgap,CliqueCuts=0, CoverCuts=0, FlowCoverCuts=0,FlowPathCuts=0,MIRCuts=0,NetworkCuts=0,GomoryPasses=0,ImpliedCuts=0, Presolve=0))
			#model = Model(solver = GurobiSolver(MIPGap=params.tolgap,Cuts=0, Presolve=0))
		else
        	model = Model(solver = GurobiSolver(TimeLimit=params.maxtime,MIPGap=params.tolgap,PreCrush=1))
		end
	elseif params.solver == "Cplex"
        model = Model(solver = CplexSolver(CPX_PARAM_TILIM=params.maxtime,CPX_PARAM_EPGAP=params.tolgap,CPX_PARAM_PREIND=0))
	else
		println("No solver selected")
		return 0
	end

	### Defining variables ###
    @variable(model,0 <= xp[i=1:inst.NP,t=1:inst.NT] <= Inf) #sum(inst.DP[i,t] for t=1:inst.NT))
    @variable(model,0 <= xr[i=1:inst.NR,t=1:inst.NT] <= Inf) #sum(inst.DR[i,t] for t=1:inst.NT))
    @variable(model,0 <= xw[i=1:inst.NW,t=1:inst.NT] <= Inf) #sum(inst.D[i,t] for t=1:inst.NT))
    @variable(model, yp[i=1:inst.NP,t=1:inst.NT], Bin)
    @variable(model, yr[i=1:inst.NR,t=1:inst.NT], Bin)
    @variable(model, yw[i=1:inst.NW,t=1:inst.NT], Bin)
    @variable(model,0 <= sp[i=1:inst.NP,t=0:inst.NT] <= Inf) #sum(inst.DP[i,t] for t=1:inst.NT))
    @variable(model,0 <= sr[i=1:inst.NR,t=0:inst.NT] <= Inf) #sum(inst.DR[i,t] for t=1:inst.NT))
    @variable(model,0 <= sw[i=1:inst.NW,t=0:inst.NT] <= Inf)


	for t in 1:inst.NT
		setcategory(yp[1,t],:Cont)
		for w in 1:inst.NW
			setcategory(yw[w,t],:Cont)
		end

		for r in 1:inst.NR
			setcategory(yr[r,t],:Cont)
		end

	end


	allretcumdem = zeros(Int,inst.NT,inst.NT)
	for t in 1:inst.NT,t2 in t:inst.NT
		allretcumdem[t,t2] = sum(inst.cumdem[:,t,t2])
	end


	warecumdem = zeros(Int,inst.NW,inst.NT,inst.NT)
	for w in 1:inst.NW
		for t in 1:inst.NT,t2 in t:inst.NT
			for k in 1:length(inst.DeltaW[w])
				r = inst.DeltaW[w][k]
				warecumdem[w,t,t2] += inst.cumdem[r,t,t2]
			end
		end
	end

println("Defined variables")
	### Objective function ###
	@objective(model, Min,
        sum(inst.SCP[i,t]*yp[i,t] + inst.HCP[i]*sp[i,t] for i=1:inst.NP, t=1:inst.NT)
        + sum(inst.SCR[i,t]*yr[i,t] + inst.HCR[i]*sr[i,t] for i=1:inst.NR, t=1:inst.NT)
        + sum(inst.SCW[i,t]*yw[i,t] + inst.HCW[i]*sw[i,t] for i=1:inst.NW, t=1:inst.NT)
	)

println("Finished objective function")

	### Setup constraints ###
	@constraint(model,
				setupP[t=1:inst.NT], xp[1,t] <= sum(inst.DP[1,k] for k in t:inst.NT)*yp[1,t]
				)
	@constraint(model,
				setupW[i=1:inst.NW, t=1:inst.NT], xw[i,t] <= sum(inst.DW[i,k] for k in t:inst.NT)*yw[i,t]
				)
	@constraint(model,
				setupR[i=1:inst.NR, t=1:inst.NT], xr[i,t] <= sum(inst.D[i,k] for k in t:inst.NT)*yr[i,t]
				)

	println("Finished setup constraints")

	### no initial inventory
	@constraint(model,
				zeroinvP[p=1:inst.NP],
				sp[p,0] == 0
				)
	@constraint(model,
				zeroinvW[w=1:inst.NW],
				sw[w,0] == 0
				)
	@constraint(model,
				zeroinvR[r=1:inst.NR],
				sr[r,0] == 0
				)

println("Finished zero inventory constraints")

	### Balance constraints ###
	@constraint(model,
				balanceP[i=1:inst.NP, t=1:inst.NT],
				sp[i,t-1] + xp[i,t] == sum(xw[w,t] for w in 1:inst.NW) + sp[i,t]
				)
	@constraint(model,
				balanceW[w=1:inst.NW,t=1:inst.NT],
				sw[w,t-1] + xw[w,t] == sum(xr[inst.DeltaW[w][k],t]   for k in 1:length(inst.DeltaW[w])) + sw[w,t]
				)
	@constraint(model,
				balanceR[i=1:inst.NR,t=1:inst.NT],
				sr[i,t-1] + xr[i,t] == inst.D[i,t] + sr[i,t]
				)


	if params.capacity != 0.0
		@constraint(model,
					capP[t=1:inst.NT], xp[1,t] <= min(inst.C[t],sum(inst.DP[1,k] for k in t:inst.NT))*yp[1,t]
					)
	end

	### Capacity constraints ###
#	@constraint(model, capacity[i=1:inst.NP, t=1:inst.NT], xp[i,t] <= min(inst.DP[i,t],inst.C[t])*yp[i,t])

	#writeLP(model,"modelo.lp",genericnames=false)



	S = zeros(Bool,inst.NR,inst.NT)
	function slowsingleretailerseparator(model)


		addedcuts = 0
		#println("Nodes = ",node)


		xp_vals = getvalue(xp)
		yp_vals = getvalue(yp)
		xw_vals = getvalue(xw)
		yw_vals = getvalue(yw)
		xr_vals = getvalue(xr)
		yr_vals = getvalue(yr)

		EPSTOL = params.tolsep

		for r in 1:inst.NR

			for l in 1:inst.NT
				alpha = 0.0
				for k in 1:l
					if yp_vals[1,k]*inst.cumdem[r,k,l] < xp_vals[1,k]
						S[r,k] = true
						alpha += yp_vals[1,k]*inst.cumdem[r,k,l]
					else
						S[r,k] = false
						alpha +=  xp_vals[1,k]
					end
				end
				#println("viol = ",inst.cumdem[r,1,l]-alpha)
				if alpha < inst.cumdem[r,1,l] - EPSTOL
#						println("r = $(r), l = $(l), violation = ",inst.cumdem[r,1,l]-alpha)

					expr = 0
					for k in 1:l
						if S[r,k]
							expr += yp[1,k]*inst.cumdem[r,k,l]
						else
							expr += xp[1,k]
						end
					end
					@constraint(model, expr >= inst.cumdem[r,1,l])
					addedcuts += 1

				end

			end
		end


		for r in 1:inst.NR

			for l in 1:inst.NT
				alpha = 0.0
				for k in 1:l
					if yw_vals[inst.DeltamR[r],k]*inst.cumdem[r,k,l] < xw_vals[inst.DeltamR[r],k]
						S[r,k] = true
						alpha += yw_vals[inst.DeltamR[r],k]*inst.cumdem[r,k,l]
					else
						S[r,k] = false
						alpha +=  xw_vals[inst.DeltamR[r],k]
					end
				end
				#println("R viol = ",inst.cumdem[r,1,l]-alpha)
				if alpha < inst.cumdem[r,1,l] - EPSTOL
#						println("W r = $(r), l = $(l), violation = ",inst.cumdem[r,1,l]-alpha)

					expr = 0
					for k in 1:l
						if S[r,k]
							expr += yw[inst.DeltamR[r],k]*inst.cumdem[r,k,l]
						else
							expr += xw[inst.DeltamR[r],k]
						end
					end
					@usercut(cb, expr >= inst.cumdem[r,1,l])
					addedcuts += 1
				end

			end
		end


		for r in 1:inst.NR

			for l in 1:inst.NT
				alpha = 0.0
				for k in 1:l
					if yr_vals[r,k]*inst.cumdem[r,k,l] < xr_vals[r,k]
						S[r,k] = true
						alpha += yr_vals[r,k]*inst.cumdem[r,k,l]
					else
						S[r,k] = false
						alpha +=  xr_vals[r,k]
					end
				end
				#println("R viol = ",inst.cumdem[r,1,l]-alpha)
				if alpha < inst.cumdem[r,1,l] - EPSTOL
#						println("R r = $(r), l = $(l), violation = ",inst.cumdem[r,1,l]-alpha)

					expr = 0
					for k in 1:l
						if S[r,k]
							expr += yr[r,k]*inst.cumdem[r,k,l]
						else
							expr += xr[r,k]
						end
					end
					@constraint(model, expr >= inst.cumdem[r,1,l])
					addedcuts += 1

				end

			end
		end

		#println("Added cuts: $(addedcuts)")


	end  # End of callback function



	S = zeros(Bool,inst.NR,inst.NT)
	SP = zeros(Bool,inst.NT)
	SW = zeros(Bool,inst.NW,inst.NT)

	function slowallretailerseparator(model)

		addedcuts = 0
		#println("Nodes = ",node)


		xp_vals = getvalue(xp)
		yp_vals = getvalue(yp)
		xw_vals = getvalue(xw)
		yw_vals = getvalue(yw)
		xr_vals = getvalue(xr)
		yr_vals = getvalue(yr)

		EPSTOL = params.tolsep

		#P  #apparently, these do not help at all
		for l in 1:inst.NT
			alpha = 0.0
			for k in 1:l
				if yp_vals[1,k]*allretcumdem[k,l] < xp_vals[1,k]
					SP[k] = true
					alpha += yp_vals[1,k]*allretcumdem[k,l]
				else
					SP[k] = false
					alpha +=  xp_vals[1,k]
				end
			end
			#println("viol = ",inst.cumdem[r,1,l]-alpha)
			if alpha < allretcumdem[1,l] - EPSTOL
			#		println("HERE P -> l = $(l), violation = ",allretcumdem[1,l]-alpha)

				expr = 0
				for k in 1:l
					if SP[k]
						expr += yp[1,k]*allretcumdem[k,l]
					else
						expr += xp[1,k]
					end
				end
				@constraint(model, expr >= allretcumdem[1,l])
				addedcuts += 1

			end

		end



		#W
		for w in 1:inst.NW

			for l in 1:inst.NT
				alpha = 0.0
				for k in 1:l
					if yw_vals[w,k]*warecumdem[w,k,l] < xw_vals[w,k]
						SW[w,k] = true
						alpha += yw_vals[w,k]*warecumdem[w,k,l]
					else
						SW[w,k] = false
						alpha +=  xw_vals[w,k]
					end
				end
				#println("R viol = ",inst.cumdem[r,1,l]-alpha)
				if alpha < warecumdem[w,1,l] - EPSTOL
#						println("W r = $(r), l = $(l), violation = ",inst.cumdem[r,1,l]-alpha)

					expr = 0
					for k in 1:l
						if SW[w,k]
							expr += yw[w,k]*warecumdem[w,k,l]
						else
							expr += xw[w,k]
						end
					end
					@constraint(model, expr >= warecumdem[w,1,l])
					addedcuts += 1
				end

			end
		end




		# R
		for r in 1:inst.NR

			for l in 1:inst.NT
				alpha = 0.0
				for k in 1:l
					if yr_vals[r,k]*inst.cumdem[r,k,l] < xr_vals[r,k]
						S[r,k] = true
						alpha += yr_vals[r,k]*inst.cumdem[r,k,l]
					else
						S[r,k] = false
						alpha +=  xr_vals[r,k]
					end
				end
				#println("R viol = ",inst.cumdem[r,1,l]-alpha)
				if alpha < inst.cumdem[r,1,l] - EPSTOL
#						println("R r = $(r), l = $(l), violation = ",inst.cumdem[r,1,l]-alpha)

					expr = 0
					for k in 1:l
						if S[r,k]
							expr += yr[r,k]*inst.cumdem[r,k,l]
						else
							expr += xr[r,k]
						end
					end
					@constraint(model, expr >= inst.cumdem[r,1,l])
					addedcuts += 1

				end

			end
		end


		return addedcuts

		# P+W
		for l in 2:inst.NT
			for lp in 1:l-1
				alpha = 0.0
				for k in 1:lp
					if yp_vals[1,k]*allretcumdem[k,l] < xp_vals[1,k]
						SP[k] = true
						alpha += yp_vals[1,k]*allretcumdem[k,l]
					else
						SP[k] = false
						alpha +=  xp_vals[1,k]
					end
				end
				#println("viol = ",inst.cumdem[r,1,l]-alpha)

				for k in lp+1:l
					for w in 1:inst.NW
						if yw_vals[w,k]*warecumdem[w,k,l] < xw_vals[w,k]
							SW[w,k] = true
							alpha += yw_vals[w,k]*warecumdem[w,k,l]
						else
							SW[w,k] = false
							alpha +=  xw_vals[w,k]
						end
					end
				end




				if alpha < allretcumdem[1,l] - EPSTOL
#						println("r = $(r), l = $(l), violation = ",inst.cumdem[r,1,l]-alpha)

					expr = 0
					for k in 1:lp
						if SP[k]
							expr += yp[1,k]*allretcumdem[k,l]
						else
							expr += xp[1,k]
						end
					end

					for k in lp+1:l
						for w in 1:inst.NW
							if SW[w,k]
								expr += yw[w,k]*warecumdem[w,k,l]
							else
								expr += xw[w,k]
							end
						end
					end


					#@constraint(model, expr >= allretcumdem[1,l])
					addedcuts += 1

				end
			end
		end



		# P+R #apparently, these do not help at all
		for l in 2:inst.NT
			for lp in 1:l-1
				alpha = 0.0
				for k in 1:lp
					if yp_vals[1,k]*allretcumdem[k,l] < xp_vals[1,k]
						SP[k] = true
						alpha += yp_vals[1,k]*allretcumdem[k,l]
					else
						SP[k] = false
						alpha +=  xp_vals[1,k]
					end
				end
				#println("viol = ",inst.cumdem[r,1,l]-alpha)

				for k in lp+1:l
					for r in 1:inst.NR
						if yr_vals[r,k]*inst.cumdem[r,k,l] < xr_vals[r,k]
							S[r,k] = true
							alpha += yr_vals[r,k]*inst.cumdem[r,k,l]
						else
							S[r,k] = false
							alpha +=  xr_vals[r,k]
						end
					end
				end


				if alpha < allretcumdem[1,l] - EPSTOL
					#println("l = $(l), violation = ",allretcumdem[1,l]-alpha)

					expr = 0
					for k in 1:lp
						if SP[k]
							expr += yp[1,k]*allretcumdem[k,l]
						else
							expr += xp[1,k]
						end
					end

					for k in lp+1:l
						for r in 1:inst.NR
							if S[r,k]
								expr += yr[r,k]*inst.cumdem[r,k,l]
							else
								expr += xr[r,k]
							end
						end
					end


					@constraint(model, expr >= allretcumdem[1,l])
					addedcuts += 1

				end
			end
		end








		#println("********************* W+R ***********************")

		for w in 1:inst.NW

			for l in 2:inst.NT

				for lw in 1:l-1
					alpha = 0.0
					for k in 1:lw
		#				println("r= $(r), lw= $(lw), l= $(l)  -> k = $(k)")
						if yw_vals[w,k]*warecumdem[w,k,l] < xw_vals[w,k]
							SW[w,k] = true
							alpha += yw_vals[w,k]*warecumdem[w,k,l]
						else
							SW[w,k] = false
							alpha +=  xw_vals[w,k]
						end
					end
					#println("R viol = ",inst.cumdem[r,1,l]-alpha)
					for ret in 1:length(inst.DeltaW[w])
						r = inst.DeltaW[w][ret]
						for k in lw+1:l
							if yr_vals[r,k]*inst.cumdem[r,k,l] < xr_vals[r,k]
								S[r,k] = true
								alpha += yr_vals[r,k]*inst.cumdem[r,k,l]
							else
								S[r,k] = false
								alpha +=  xr_vals[r,k]
							end
						end
					end
					if alpha < warecumdem[w,1,l] - EPSTOL
		#				println("WR r = $(r), l = $(l), violation = ",inst.cumdem[r,1,l]-alpha)

						expr = 0
						for k in 1:lw
							if SW[w,k]
								expr += yw[w,k]*warecumdem[w,k,l]
							else
								expr += xw[w,k]
							end
						end

						#println("R r = $(r), l = $(l), violation = ",inst.cumdem[r,1,l]-alpha)
						for ret in 1:length(inst.DeltaW[w])
							r = inst.DeltaW[w][ret]
							for k in lw+1:l
								if S[r,k]
									expr += yr[r,k]*inst.cumdem[r,k,l]
								else
									expr += xr[r,k]
								end
							end
						end
						#@constraint(model, expr >= warecumdem[w,1,l])
						addedcuts += 1

					end
				end
			end
		end



		return addedcuts

		#println("Added cuts: $(addedcuts)")


	end  # End of callback function



	function fastallretailerseparator(model)
		lk = zeros(Int,inst.NT)
		addedcuts = 0
		#println("Nodes = ",node)

		alpha = zeros(Float64,inst.NT)


		xp_vals = getvalue(xp)
		yp_vals = getvalue(yp)
		xw_vals = getvalue(xw)
		yw_vals = getvalue(yw)
		xr_vals = getvalue(xr)
		yr_vals = getvalue(yr)

		EPSTOL = params.tolsep

		changel = Vector{Vector{Int}}(undef,inst.NT)
		for t in 1:inst.NT
			changel[t]=[]
		end
	#	println("\n\n\n !!!!!!!!!!!!!!!!!!!!!!!! Production site")
		for k in 1:inst.NT
			#println("\n\n k is ",k)
			#println("yk is ",yp_vals[1,k])
			#println("xk will be ",xp_vals[1,k])
			#println("will search in ",yp_vals[1,k]*allretcumdem[k,k:inst.NT])
			if yp_vals[1,k] > params.tolgap
				lk[k] = binarysearchfindlsinglelevel(k,xp_vals[1,k],yp_vals[1,k],allretcumdem,k,inst.NT)
				if lk[k] != -1
					push!(changel[lk[k]],k)
				end
			else
				lk[k] = -1
			end
			#println("lk[$(k)] is ",lk[k])
			#println("\n")

		end

		#println("changel = ",changel)

		sumY = 0.0
		sumX = 0.0
		bestviol = 0.0
		bestviolindex = 0
		if lk[1] == 1
			alpha[1] = xp_vals[1,1]
			sumX += xp_vals[1,1]
		else
			alpha[1] = yp_vals[1,1]*allretcumdem[1,1]
			sumY += yp_vals[1,1]
		end

		#println("L = ",1)
		#println("alpha[l] = $(alpha[1])")
		#println("violation = ",alpha[1] - allretcumdem[1,1])

		if alpha[1] - allretcumdem[1,1] < bestviol
			bestviol = alpha[1] - allretcumdem[1,1]
			bestviolindex = 1
		end

		for l in 2:inst.NT
		#	println("L = ",l)
			alpha[l] = alpha[l-1]

			sumY += yp_vals[1,l]
			if length(changel[l]) > 0
		#		println("length changel = $(length(changel[l]))")
				for perind in 1:length(changel[l])
					per = changel[l][perind]
		#			println("removing $(per) from Y_l")
					sumX += xp_vals[1,per]
					sumY -= yp_vals[1,per]
					alpha[l] += xp_vals[1,per]
					alpha[l] -= yp_vals[1,per]*allretcumdem[per,l-1]
				end
			end
			alpha[l] += sumY*allretcumdem[l,l]

			if alpha[l] - allretcumdem[1,l] < bestviol
				bestviol = alpha[l] - allretcumdem[1,l]
				bestviolindex = l
			end

		#	println("alpha[l] = $(alpha[l])")
		#	println("violation = ",alpha[l] - allretcumdem[1,l])


		end


		if bestviol < -EPSTOL
			l = bestviolindex
			expr = 0
			for k in 1:l
				if yp_vals[1,k]*allretcumdem[k,l] < xp_vals[1,k]
					expr += yp[1,k]*allretcumdem[k,l]
				else
					expr += xp[1,k]
				end
			end
			@constraint(model, expr >= allretcumdem[1,l])
			addedcuts += 1
		end

		#println("addedcuts = $(addedcuts)")
		#println("bestviol = $(bestviol)   ---->   period = $(bestviolindex)")
		for w in 1:inst.NW
#				println("\n\n\n ################ Warehouse : $(w)")
			for t in 1:inst.NT
				changel[t]=[]
			end
			for k in 1:inst.NT
				#println("\n k is ",k)
				#println("yk is ",yw_vals[w,k])
				#println("xk will be ",xw_vals[w,k])
				#println("will search in ",yw_vals[w,k]*warecumdem[w,k,k:inst.NT])
				if yw_vals[w,k] > params.tolgap
					lk[k] = binarysearchfindlsinglelevel(k,xw_vals[w,k],yw_vals[w,k],warecumdem[w,:,:],k,inst.NT)
					if lk[k] != -1
						push!(changel[lk[k]],k)
					end
				else
					lk[k] = -1
				end
				#println("lk[$(k)] is ",lk[k])
				#println("\n")
			end
			#println("changel = ",changel)


			sumY = 0.0
			sumX = 0.0
			bestviol = 0.0
			bestviolindex = 0
			if lk[1] == 1
				alpha[1] = xw_vals[w,1]
				sumX += xw_vals[w,1]
			else
				alpha[1] = yw_vals[w,1]*warecumdem[w,1,1]
				sumY += yw_vals[w,1]
			end

			#println("L = ",1)
			#println("alpha[l] = $(alpha[1])")
			#println("violation = ",alpha[1] - allretcumdem[1,1])

			if alpha[1] - warecumdem[w,1,1] < bestviol
				bestviol = alpha[1] - warecumdem[w,1,1]
				bestviolindex = 1
			end

			for l in 2:inst.NT
			#	println("L = ",l)
				alpha[l] = alpha[l-1]

				sumY += yw_vals[w,l]
				if length(changel[l]) > 0
			#		println("length changel = $(length(changel[l]))")
					for perind in 1:length(changel[l])
						per = changel[l][perind]
			#			println("removing $(per) from Y_l")
						sumX += xw_vals[w,per]
						sumY -= yw_vals[w,per]
						alpha[l] += xw_vals[w,per]
						alpha[l] -= yw_vals[w,per]*warecumdem[w,per,l-1]
					end
				end
				alpha[l] += sumY*warecumdem[w,l,l]

				if alpha[l] - warecumdem[w,1,l] < bestviol
					bestviol = alpha[l] - warecumdem[w,1,l]
					bestviolindex = l
				end

			#	println("alpha[l] = $(alpha[l])")
			#	println("violation = ",alpha[l] - allretcumdem[1,l])


			end

			if bestviol < -EPSTOL
				l = bestviolindex
				expr = 0
				for k in 1:l
					if yw_vals[w,k]*warecumdem[w,k,l] < xw_vals[w,k]
						expr += yw[w,k]*warecumdem[w,k,l]
					else
						expr += xw[w,k]
					end
				end
				@constraint(model, expr >= warecumdem[w,1,l])
				addedcuts += 1
			end


		end


		for r in 1:inst.NR
			#println("\n\n\n ***************** Retailer : $(r)")
			for t in 1:inst.NT
				changel[t]=[]
			end
			for k in 1:inst.NT
				#println("\n k is ",k)
				#println("yr is ",yr_vals[r,k])
				#println("xk will be ",xr_vals[r,k])
				#println("will search in ",yr_vals[r,k]*inst.cumdem[r,k,k:inst.NT])
				if yr_vals[r,k] > params.tolgap
					lk[k] = binarysearchfindlsinglelevel(k,xr_vals[r,k],yr_vals[r,k],inst.cumdem[r,:,:],k,inst.NT)
					if lk[k] != -1
						push!(changel[lk[k]],k)
					end
				else
					lk[k] = -1
				end
				#println("lk[$(k)] is ",lk[k])
				#println("\n")
			end
			#println("changel = ",changel)


			sumY = 0.0
			sumX = 0.0
			bestviol = 0.0
			bestviolindex = 0
			if lk[1] == 1
				alpha[1] = xr_vals[r,1]
				sumX += xr_vals[r,1]
			else
				alpha[1] = yr_vals[r,1]*inst.cumdem[r,1,1]
				sumY += yr_vals[r,1]
			end

			#println("L = ",1)
			#println("alpha[l] = $(alpha[1])")
			#println("violation = ",alpha[1] - allretcumdem[1,1])

			if alpha[1] - inst.cumdem[r,1,1] < bestviol
				bestviol = alpha[1] - inst.cumdem[r,1,1]
				bestviolindex = 1
			end

			for l in 2:inst.NT
			#	println("L = ",l)
				alpha[l] = alpha[l-1]

				sumY += yr_vals[r,l]
				if length(changel[l]) > 0
			#		println("length changel = $(length(changel[l]))")
					for perind in 1:length(changel[l])
						per = changel[l][perind]
			#			println("removing $(per) from Y_l")
						sumX += xr_vals[r,per]
						sumY -= yr_vals[r,per]
						alpha[l] += xr_vals[r,per]
						alpha[l] -= yr_vals[r,per]*inst.cumdem[r,per,l-1]
					end
				end
				alpha[l] += sumY*inst.cumdem[r,l,l]

				if alpha[l] - inst.cumdem[r,1,l] < bestviol
					bestviol = alpha[l] - inst.cumdem[r,1,l]
					bestviolindex = l
				end

			#	println("alpha[l] = $(alpha[l])")
			#	println("violation = ",alpha[l] - allretcumdem[1,l])


			end

			#println("bestviol = $(bestviol)")
			if bestviol < -EPSTOL
				l = bestviolindex
				expr = 0
				for k in 1:l
					if yr_vals[r,k]*inst.cumdem[r,k,l] < xr_vals[r,k]
						expr += yr[r,k]*inst.cumdem[r,k,l]
					else
						expr += xr[r,k]
					end
				end
				@constraint(model, expr >= inst.cumdem[r,1,l])
				addedcuts += 1
			end






		end


		return addedcuts
		#println("Array lk = ",lk)

	end  # End of callback function







	#addcutcallback(model,slowsingleretailerseparator)
	#addcutcallback(model,slowallretailerseparator)

	#addcutcallback(model,fastallretailerseparator)

	t1 = time_ns()

	keepdoing = true
	rounds = 0
	while keepdoing
		rounds += 1
		status = solve(model)
		#addedcuts = slowallretailerseparator(model)
		addedcuts = fastallretailerseparator(model)
		#println("addedcuts = ",addedcuts)
		#println("newobjective = ",getobjectivevalue(model))
		if addedcuts > 0
			keepdoing = true
		else
			keepdoing = false
		end

	end


	t2 = time_ns()
	elapsedtime = (t2-t1)/1.0e9
	newobjective = getobjectivevalue(model)
	println("newobjective = ",newobjective)
	println("time = ", elapsedtime)
	println("rounds = ",rounds)
	timeperround = elapsedtime/rounds
	println("time/round = ",timeperround)

	open("saida.txt","a") do f
		write(f,";$newobjective;$elapsedtime;$rounds;$timeperround;$(params.disablesolver) \n")
	end



#	if params.printsol == 1
#		printStandardFormulationSolution(inst,xp,xr,xw,yp,yr,yw,sp,sr,sw)
#	end









end #function cuttingplanestandardFormulation()













end #module
