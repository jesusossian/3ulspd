for i in instances/Preliminary/*
do
	echo Solving $i with multi	
	julia 3lspd.jl --form cuttingplane --inst $i  --balanced 1 --solver Gurobi --disablesolver 1
done
