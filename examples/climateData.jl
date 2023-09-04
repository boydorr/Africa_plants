using Africa
using Diversity
using JLD2
using Unitful
using Unitful.DefaultSymbols
using EcoSISTEM.Units
using EcoSISTEM.ClimatePref
using EcoSISTEM
using Distributions
using AxisArrays
using Random

## Load CERA climates
dir1 = "data/CERA"
tempax1 = readCERA(dir1, "cera_20c_temp2m", "t2m")
precax1 = readCERA(dir1, "cera_20c_totalprec", "tp")
soilwaterax1 = readCERA(dir1, "cera_20c_soilwater1", "swvl1")
solarradax1 = readCERA(dir1, "cera_20c_surfacenetsolar", "ssr")

times = [collect(1979year:1month:(1980year - 1.0month)),
    collect(1980year:1month:(1990year - 1.0month)),
    collect(1990year:1month:(2000year - 1.0month)),
    collect(2000year:1month:(2010year - 1.0month)),
    collect(2010year:1month:(2019year- 1.0month))]

# Load ERA climates
dir2 = "data/ERA"
tempax2 = readERA(dir2, "era_int_temp2m_", "t2m", times)
precax2 = readERA(dir2, "era_int_totalprec_", "tp", times)
soilwaterax2 = readERA(dir2, "era_int_soilwater1", "swvl1", times)
solarradax2 = readERA(dir2, "era_int_netsolar", "ssr", times)  

# Concatenate and crop to Africa
temp = cat(tempax1.array[:, :, 1901year .. (1979year-1month)], tempax2.array, dims = 3)
prec = cat(precax1.array[:, :, 1901year .. (1979year-1month)], precax2.array, dims = 3)
prec[prec .< 0m] *= 0
soilwater = cat(soilwaterax1.array[:, :, 1901year .. (1979year-1month)], soilwaterax2.array, dims = 3)
soilwater[soilwater .< 0m^3] *= 0
solarrad = cat(solarradax1.array[:, :, 1901year .. (1979year-1month)], solarradax2.array, dims = 3)
solarrad[solarrad .< 0J/m^2] *= 0

africa_temp = ERA(temp[-25°.. 50°, -35°.. 40°, :])
africa_prec = ERA(prec[-25°.. 50°, -35°.. 40°, :])
africa_prec.array = AxisArray(uconvert.(mm, africa_prec.array), africa_prec.array.axes)
JLD2.@save "Africa_temp.jld2" africa_temp
JLD2.@save "Africa_prec.jld2" africa_prec

# Convert water and solar to the right area/volume for one grid square
africa_water = uconvert.(m^3, soilwater[-25°.. 50°, -35°.. 40°, :] .* (80km * 80km * 7cm) ./ m^3)
africa_solar = uconvert.(kJ, solarrad[-25°.. 50°, -35°.. 40°, :] .* (80km * 80km))
africa_solar = SolarTimeBudget(africa_solar, 1)
africa_water =VolWaterTimeBudget(africa_water, 1)
JLD2.@save "Africa_water.jld2" africa_water
JLD2.@save "Africa_solar.jld2" africa_solar