using EcoSISTEM
using Unitful
using Unitful.DefaultSymbols
using Distributions
using EcoSISTEM.Units
using Phylo
using DataFrames
using Diversity

function create_eco(paramDict::Dict, abenv::A; bound::B = Torus(), size = (mean = 1.0m^2, std = 0.1m^2)) where {A <: EcoSISTEM.AbstractAbiotic, B <: EcoSISTEM.BoundaryCondition}
    # Set up initial parameters for ecosystem
    birth = haskey(paramDict, "birth") ? paramDict["birth"] : 0.6/year
    death = haskey(paramDict, "death") ? paramDict["death"] : 0.6/year
    l = haskey(paramDict, "l") ? paramDict["l"] : 1.0
    s = haskey(paramDict, "s") ? paramDict["s"] : 0.2
    boost = haskey(paramDict, "boost") ? paramDict["boost"] : 1.0

    numSpecies = paramDict["numSpecies"]
    numInvasive = paramDict["numInvasive"]
    individuals = paramDict["numIndiv"]
    req = paramDict["reqs"]
    opts = paramDict["opts"]
    vars = paramDict["vars"]
    kernel = paramDict["kernel"]

    # Set up how much energy each species consumes
    #energy_vec = SolarRequirement(fill(req, numSpecies))
    names = map(x -> "$x", 1:(numSpecies + numInvasive))
    tree = rand(Ultrametric{BinaryTree{DataFrame, DataFrame}}(names))
    #rand(Ultrametric{RootedTree}(names))
    units = unit(size.mean)
    trts = ContinuousEvolve(uconvert(NoUnits, size.mean / units),
                            uconvert(NoUnits, size.std / units), tree)
    # rand(BrownianTrait(tree, "Size",
    #                    size.mean, σ=size.std, f=abs))
    # trts = getnodedata.(tree, names, "Size")
    # energy_vec1 = SolarRequirement(trts .* req[1])
    # energy_vec2 = WaterRequirement(trts .* req[2])

    energy_vec1 = SolarRequirement(abs.(trts.mean) .* (req[1] * units))
    energy_vec2 = WaterRequirement(abs.(trts.mean) .* (req[2] * units))

    energy_vec = ReqCollection2(energy_vec1, energy_vec2)
    # Collect model parameters together (in this order!!)
    if length(birth) > 1
        param = PopGrowth{typeof(unit(birth[1]))}(birth, death, l, s , boost)
    else
        param = EqualPop(birth, death, l, s, boost)
    end

    # Create ecosystem

    movement = BirthOnlyMovement(kernel, bound)

    traits = GaussTrait(opts, vars)
    native = fill(true, numSpecies + numInvasive)
    native[(numSpecies+1):end] .= false
    abun = [rand(Multinomial(individuals, numSpecies)); fill(0, numInvasive)]
    sppl = SpeciesList(numSpecies + numInvasive, traits, abun, energy_vec,
        movement, param, native)
    rel = Gauss{typeof(first(paramDict["opts"]))}()
    eco = Ecosystem(traitpopulate!, sppl, abenv, rel)
end
