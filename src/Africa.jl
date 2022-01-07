module Africa

include("Tools.jl")
export emptypopulate!, reverseBurnin!, db_to_array, startingArray

include("RunSim.jl")
export runsim!, runburnin!, keepYear!

end # module
