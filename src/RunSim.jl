using EcoSISTEM
using Unitful
using EcoSISTEM.Units
using JLD2
using Printf

function keepYear!(eco::Ecosystem)
  newhabchange1 = HabitatUpdate(eraChange, eco.abenv.habitat.h1.change.rate)
  newhabchange2 = HabitatUpdate(eraChange, eco.abenv.habitat.h2.change.rate)
  eco.abenv.habitat.h1.change = newhabchange1
  eco.abenv.habitat.h2.change = newhabchange2
end

function forwardTime!(eco::Ecosystem)
  newhabchange1 = HabitatUpdate(eraChange, eco.abenv.habitat.h1.change.rate)
  newhabchange2 = HabitatUpdate(eraChange, eco.abenv.habitat.h2.change.rate)
  eco.abenv.habitat.h1.change = newhabchange1
  eco.abenv.habitat.h2.change = newhabchange2
  eco.abenv.habitat.h1.time = 1
  eco.abenv.habitat.h2.time = 1
end

function eraSteady(eco::Ecosystem, hab::ContinuousTimeHab, timestep::Unitful.Time)
    monthstep = uconvert(month, timestep)
    hab.time += round(Int64, monthstep/month)
    if hab.time > 12
        hab.time = 1
        @warn "More timesteps than available, have repeated"
    end
end

function runburnin!(eco::Ecosystem, duration::Unitful.Time, timestep::Unitful.Time)
  times = length(0s:timestep:duration)
  keepYear!(eco)
  for i in 1:times
    update!(eco, timestep)
  end
  forwardTime!(eco)
end

function simulate_record_diversity!(storage::AbstractArray, eco::Ecosystem, times::Unitful.Time, interval::Unitful.Time, timestep::Unitful.Time, divfuns::Array{Function, 1}, qs::Vector{Float64}, cacheInterval::Unitful.Time, cacheFolder::String, scenario_name::String)
  mod(interval,timestep) == 0.0year || error("Interval must be a multiple of timestep")
  record_seq = 0s:interval:times
  time_seq = 0s:timestep:times
  counting = 0
  for i in 1:length(time_seq)
      update!(eco, timestep);
      # Record diversity profiles for each measure
      if time_seq[i] in record_seq
          counting = counting + 1
        for j in eachindex(divfuns)
          storage[j, :, counting, :] .= reshape(divfuns[j](eco, qs)[:diversity], size(storage, 2), size(storage, 4))
        end
      end
      # Save cache of abundances
      if mod(time_seq[i], cacheInterval) == 0year
          JLD.save(joinpath(cacheFolder, scenario_name * (@sprintf "%02d.jld" uconvert(NoUnits,time_seq[i]/cacheInterval))), "abun", eco.abundances.matrix)
      end
  end
  storage
end

function runsim!(div::Array{Float64, 4}, eco::Ecosystem, simDict::Dict, folder::String = "./") where AB <: EcoSISTEM.AbstractAbiotic
    cacheFolder = joinpath(folder, "cache")
    isdir(cacheFolder) || mkdir(cacheFolder)
    divfuns = simDict["divfuns"]
    q = simDict["q"]
    simulate!(eco, simDict["burnin"], simDict["timestep"])
    simulate_record_diversity!(div, eco, simDict["times"], simDict["interval"], simDict["timestep"], divfuns, q, simDict["cacheInterval"], cacheFolder, simDict["fileName"])
    JLD.save(joinpath(folder, simDict["fileName"] * ".jld"), "div", div)
end
