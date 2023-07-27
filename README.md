# 3ulspd

Three-level lot sizing and replenishment problem

for run one instances with default configuration
> julia 3ulspd.jl

for run one instance with formulation solver
> julia 3ulspd.jl --inst XXXX --form XXXX --solver XXXX

There a .sh to run set of instances
./executar.sh

for run one instance with formulation, solver and capacity
> julia 3ulspd.jl --inst XXXX --form XXXX --solver XXXX --capacity XXXX

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

Instances group

## Uncapacitated instances

P = { 1 }

N = {50, 100, 200}

T = {15, 30}

W = {5, 10, 15, 20}

V = {SD_SF, SD_DF, DD_SF, DD_DF}
SD_SF = static demand and static fixed cost
SD_DF = static demand and dynamic fixed cost
DD_SF = dynamic demand and static fixed cost
DD_DF = dynamic demand and dynamic fixed cost

instances = {1,2, \ldots, 5}

style = {bal, unbal}

Number instances uncapacitated = 3 * 2 * 4 * 4 * 5 * 2 = 960 instances 
