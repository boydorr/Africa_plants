using JLD
using Diversity
using Plots
using Unitful
using Unitful.DefaultSymbols
using EcoSISTEM.Units
using AxisArrays
using EcoSISTEM.ClimatePref
using Statistics
plotlyjs()

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
dir1 = "data/ERA"
tempax2 = readERA(dir1, "era_int_temp2m_", "t2m", times)
precax2 = readERA(dir1, "era_int_totalprec_", "tp", times)
soilwaterax2 = readERA(dir1, "era_int_soilwater1", "swvl1", times)
solarradax2 = readERA(dir1, "era_int_netsolar", "ssr", times)

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
africa_water = uconvert.(m^3, soilwater[-25°.. 50°, -35°.. 40°, :] .* (80km * 80km * 7cm) ./ m^3)
africa_solar = uconvert.(kJ, solarrad[-25°.. 50°, -35°.. 40°, :] .* (80km * 80km))
sa_temps = ERA(temp[-85°.. -30.25°, -60.25°.. 20°, :])
sa_prec = ERA(prec[-85°.. -30.25°, -60.25°.. 20°, :])
sa_prec.array = AxisArray(uconvert.(mm, sa_prec.array), sa_prec.array.axes)
sa_water = uconvert.(m^3, soilwater[-85°.. -30.25°, -60.25°.. 20°, :] .* (80km * 80km * 7cm) ./ m^3)
sa_solar = uconvert.(kJ, solarrad[-85°.. -30.25°, -60.25°.. 20°, :] .* (80km * 80km))


africa = readfile("data/Africa.tif", -25°, 50°, -35°, 40°)[:, end:-1:1]
sa = readfile("data/SouthAmerica.tif", -85°, -30.25°, -60.25°, 20°)[:, end:-1:1]
active_a =  Array{Bool, 2}(.!isnan.(africa))
active_sa =  Array{Bool, 2}(.!isnan.(sa))

a_temp = ustrip.(uconvert.(°C, mapslices(mean, Array(africa_temp.array), dims = 3)[:, :, 1]))
a_temp[.!(active_a)] .= NaN
sa_temp = ustrip.(uconvert.(°C, mapslices(mean, Array(sa_temps.array), dims = 3)[:, :, 1]))
sa_temp[.!(active_sa)] .= NaN

a_times = map(0:117) do i
    ustrip.(uconvert.(°C, mapslices(mean, Array(africa_temp.array[:,:,(i*12 + 1):min(1415,(i* 12 +12))]), dims = 3)[:, :, 1]))
end
a_time_temp = mean.(a_times)
sa_times = map(0:117) do i
    ustrip.(uconvert.(°C, mapslices(mean, Array(sa_temps.array[:,:,(i*12 + 1):min(1415,(i* 12 +12))]), dims = 3)[:, :, 1]))
end
a_time_temp = mean.(a_times)
sa_time_temp = mean.(sa_times)

heatmap(-25:0.75:49.25, -35:0.75:39.25, transpose(a_temp), background_color = :lightblue, background_color_outside=:white, grid = false, colorbar_title = "Temperature (°C)", guidefontsize =12,tickfontsize=12, titlefontsize=16, margin = 1.0*Plots.mm, layout = (@layout [a{0.5w} b{0.5w}]),size = (1000, 300), aspect_ratio = 1, clim = (5, 35),legendfontsize= 12)
heatmap!(-84.25:0.75:-31, -59.75:0.75:19.25, transpose(sa_temp), subplot = 2, background_color = :lightblue, background_color_outside=:white, grid = false, colorbar_title = "Temperature (°C)", clim = (5, 35), guidefontsize =12,tickfontsize=12, titlefontsize=16, margin = 1.0*Plots.mm, aspect_ratio = 1, legendfontsize= 12, right_margin = 0.0*Plots.mm)
Plots.pdf("plots/Climates1.pdf")

plot(1901:2018, a_time_temp, background_color = :white, background_color_outside=:white, grid = false, guidefontsize =12,tickfontsize=12, titlefontsize=16, margin = 1.0*Plots.mm, legend= false, ylim = (6, 28), size = (1000, 200), layout = (@layout [a{0.45w} b{0.45w}]), color = 2, ribbon = std.(a_times))
plot!(1901:2018, sa_time_temp, subplot=2, background_color = :white, background_color_outside=:white, grid = false, colorbar_title = "Temperature (°C)", guidefontsize =12,tickfontsize=12, titlefontsize=16, margin = 1.0*Plots.mm, legend = false, ylim = (6, 28), ribbon = std.(sa_times), color = 2)
Plots.pdf("plots/ClimatesMeans1.pdf")

a_temp = ustrip.(uconvert.(°C, mapslices(minimum, Array(africa_temp.array), dims = 3)[:, :, 1]))
a_temp[.!(active_a)] .= NaN
sa_temp = ustrip.(uconvert.(°C, mapslices(minimum, Array(sa_temps.array), dims = 3)[:, :, 1]))
sa_temp[.!(active_sa)] .= NaN

a_times = map(0:117) do i
    ustrip.(uconvert.(°C, mapslices(minimum, Array(africa_temp.array[:,:,(i*12 + 1):min(1415,(i* 12 +12))]), dims = 3)[:, :, 1]))
end
a_time_temp = mean.(a_times)
sa_times = map(0:117) do i
    ustrip.(uconvert.(°C, mapslices(minimum, Array(sa_temps.array[:,:,(i*12 + 1):min(1415,(i* 12 +12))]), dims = 3)[:, :, 1]))
end
a_time_temp = mean.(a_times)
sa_time_temp = mean.(sa_times)

heatmap(-25:0.75:49.25, -35:0.75:39.25, transpose(a_temp), background_color = :lightblue, background_color_outside=:white, grid = false, colorbar_title = "Minimum Temperature (°C)", guidefontsize =12,tickfontsize=12, titlefontsize=16, margin = 1.0*Plots.mm, layout = (@layout [a{0.5w} b{0.5w}]),size = (1000, 300), aspect_ratio = 1, clim = (5, 35),legendfontsize= 12)
heatmap!(-84.25:0.75:-31, -59.75:0.75:19.25, transpose(sa_temp), subplot = 2, background_color = :lightblue, background_color_outside=:white, grid = false, colorbar_title = "Minimum Temperature (°C)", clim = (5, 35), guidefontsize =12,tickfontsize=12, titlefontsize=16, margin = 1.0*Plots.mm, aspect_ratio = 1, legendfontsize= 12, right_margin = 0.0*Plots.mm)
Plots.pdf("plots/Climates1.5.pdf")

plot(1901:2018, a_time_temp, background_color = :white, background_color_outside=:white, grid = false, guidefontsize =12,tickfontsize=12, titlefontsize=16, margin = 1.0*Plots.mm, legend= false, ylim = (0, 28), size = (1000, 200), layout = (@layout [a{0.45w} b{0.45w}]), color = 2, ribbon = std.(a_times))
plot!(1901:2018, sa_time_temp, subplot=2, background_color = :white, background_color_outside=:white, grid = false, colorbar_title = "Temperature (°C)", guidefontsize =12,tickfontsize=12, titlefontsize=16, margin = 1.0*Plots.mm, legend = false, ylim = (0, 28), ribbon = std.(sa_times), color = 2)
Plots.pdf("plots/ClimatesMeans1.5.pdf")


a_rains = map(0:117) do i
    mapslices(sum, Array(africa_prec.array[:,:,(i*12 + 1):min(1415,(i* 12 +12))]), dims = 3)[:, :, 1] ./ mm
end
a_rain = mapslices(mean, reshape(hcat(a_rains...), 100, 100, 118),dims = 3)[:, :, 1]
a_rain[.!(active_a)] .= NaN

sa_rains = map(0:117) do i
    mapslices(sum, Array(sa_prec.array[:,:,(i*12 + 1):min(1415,(i* 12 +12))]), dims = 3)[:, :, 1] ./ mm
end
sa_rain = mapslices(mean, reshape(hcat(sa_rains...), 73, 107, 118),dims = 3)[:, :, 1]
sa_rain[.!(active_sa)] .= NaN

a_time_rain = mean.(a_rains)
sa_time_rain = mean.(sa_rains)

heatmap(-25:0.75:49.25, -35:0.75:39.25, transpose(a_rain), background_color = :lightblue, background_color_outside=:white, grid = false, clim = (0, 140), color = :viridis_r, colorbar_entry = false, guidefontsize =12,tickfontsize=12, titlefontsize=16, margin = 1.0*Plots.mm, layout = (@layout [a{0.5w} b{0.5w}]),size = (1000, 300), aspect_ratio = 1, legendfontsize= 12)
heatmap!(-84.25:0.75:-31, -59.75:0.75:19.25, transpose(sa_rain), background_color = :lightblue, background_color_outside=:white, grid = false, colorbar_title = "Rainfall (mm)", subplot = 2, clim = (0, 140), color = :viridis_r,guidefontsize =12,tickfontsize=12, titlefontsize=16, margin = 1.0*Plots.mm, aspect_ratio = 1, legendfontsize= 12, right_margin = 0.0*Plots.mm)

Plots.pdf("plots/Climates2.pdf")

plot(1901:2018, a_time_rain, background_color = :white, background_color_outside=:white, grid = false, guidefontsize =12,tickfontsize=12, titlefontsize=16, margin = 1.0*Plots.mm, legend= false, ylim = (0, 80), size = (1000, 200), layout = (@layout [a{0.45w} b{0.45w}]), ribbon = std.(a_rains), color = 3)
plot!(1901:2018, sa_time_rain, subplot=2, background_color = :white, background_color_outside=:white, grid = false, colorbar_title = "Temperature (°C)", guidefontsize =12,tickfontsize=12, titlefontsize=16, margin = 1.0*Plots.mm, legend = false, ylim = (0, 80), ribbon = std.(sa_rains), color = 3)
Plots.pdf("plots/ClimatesMeans2.pdf")

a_solar = uconvert.(kW, mapslices(mean, africa_solar, dims = 3)[:, :, 1] ./ month) ./ kW
a_solar[.!(active_a)] .= NaN
sa_solar = uconvert.(kW, mapslices(mean, sa_solar, dims = 3)[:, :, 1] ./ month) ./ kW
sa_solar[.!(active_sa)] .= NaN

heatmap(-25:0.75:49.25, -35:0.75:39.25, transpose(a_solar), background_color = :lightblue, background_color_outside=:white, grid = false, layout = (@layout [a{0.5w} b{0.5w}]), clim = (2e7, 6e7), color = :solar, colorbar_entry = false, size = (1000, 300), guidefontsize =12,tickfontsize=12, titlefontsize=16, margin = 1.0*Plots.mm, aspect_ratio=1)
heatmap!(-84.25:0.75:-31, -59.75:0.75:19.25, transpose(sa_solar), background_color = :lightblue, background_color_outside=:white, grid = false, colorbar_title = "Solar radiation (kW)", subplot = 2, clim = (2e7, 6e7), color = :solar, guidefontsize =12,tickfontsize=12, titlefontsize=16, margin = 1.0*Plots.mm, aspect_ratio=1, right_margin = 0.0 *Plots.mm)
Plots.pdf("plots/Resources1.pdf")

a_times = map(0:117) do i
    uconvert.(kW, mapslices(mean, Array(africa_solar[:,:,(i*12 + 1):min(1415,(i* 12 +12))]), dims = 3)[:, :, 1]./ month) ./ kW
end
sa_times = map(0:117) do i
    uconvert.(kW, mapslices(mean, Array(sa_solar[:,:,(i*12 + 1):min(1415,(i* 12 +12))]), dims = 3)[:, :, 1] ./ month) ./ kW
end
a_time_sol = mean.(a_times)
sa_time_sol = mean.(sa_times)
plot(1901:2018, a_time_sol./1e6, background_color = :white, background_color_outside=:white, grid = false, guidefontsize =12,tickfontsize=12, titlefontsize=16, margin = 1.0*Plots.mm, legend= false, ylim = (20, 50), size = (1000, 200), layout = (@layout [a{0.45w} b{0.45w}]), ribbon = std.(a_times)./1e6, color = 5)
plot!(1901:2018, sa_time_sol./1e6, subplot=2, background_color = :white, background_color_outside=:white, grid = false, guidefontsize =12,tickfontsize=12, titlefontsize=16, margin = 1.0*Plots.mm, legend = false, ylim = (20, 50), ribbon = std.(sa_times)./1e6, color = 5)
Plots.pdf("plots/ResourceMeans1.pdf")

a_waters = map(0:117) do i
    mapslices(sum, Array(africa_water[:,:,(i*12 + 1):min(1415,(i* 12 +12))]), dims = 3)[:, :, 1] ./ m^3
end
a_water = mapslices(mean, reshape(hcat(a_waters...), 100, 100, 118),dims = 3)[:, :, 1]
a_water[.!(active_a)] .= NaN

sa_waters = map(0:117) do i
    mapslices(sum, Array(sa_water[:,:,(i*12 + 1):min(1415,(i* 12 +12))]), dims = 3)[:, :, 1] ./ m^3
end
sa_water = mapslices(mean, reshape(hcat(sa_waters...), 73, 107, 118),dims = 3)[:, :, 1]
sa_water[.!(active_sa)] .= NaN

heatmap(-25:0.75:49.25, -35:0.75:39.25, transpose(a_water), background_color = :lightblue, background_color_outside=:white, grid = false, clim = (0, 2.5e9), color = :ice_r, colorbar_entry = false, guidefontsize =12,tickfontsize=12, titlefontsize=16, margin = 1.0*Plots.mm, layout = (@layout [a{0.5w} b{0.5w}]),size = (1000, 300), aspect_ratio = 1, legendfontsize= 12)
heatmap!(-84.25:0.75:-31, -59.75:0.75:19.25, transpose(sa_water), background_color = :lightblue, background_color_outside=:white, grid = false, colorbar_title = "Soil water volume (m³)",  color = :ice_r, subplot = 2, clim = (0, 2.5e9), guidefontsize =12,tickfontsize=12, titlefontsize=16, margin = 1.0*Plots.mm, right_margin = 0.0*Plots.mm, aspect_ratio = 1, legendfontsize= 12)
Plots.pdf("plots/Resources2.pdf")

plot(1901:2018, mean.(a_waters)./1e8, background_color = :white, background_color_outside=:white, grid = false, guidefontsize =12,tickfontsize=12, titlefontsize=16, margin = 1.0*Plots.mm, legend= false, ylim = (0, 18), size = (1000, 200), layout = (@layout [a{0.45w} b{0.45w}]), ribbon = std.(a_waters)./1e8, color = 13)
plot!(1901:2018, mean.(sa_waters)./1e8, subplot=2, background_color = :white, background_color_outside=:white, grid = false, colorbar_title = "Temperature (°C)", guidefontsize =12,tickfontsize=12, titlefontsize=16, margin = 1.0*Plots.mm, legend = false, ylim = (0, 18), ribbon = std.(sa_waters)./1e8, color = 13)
Plots.pdf("plots/ResourceMeans2.pdf")


plot(ustrip.(uconvert.(°C, mapslices(mean, Array(sa_temps.array), dims = (1,2))[1, 1, :]))[1:12])
plot!(ustrip.(uconvert.(°C, mapslices(mean, Array(sa_temps.array), dims = (1,2))[1, 1, :]))[end-11:end])

plot(ustrip.(uconvert.(°C, mapslices(mean, Array(africa_temp.array), dims = (1,2))[1, 1, :]))[1:12], clim = (18, 24))
plot!(ustrip.(uconvert.(°C, mapslices(mean, Array(africa_temp.array), dims = (1,2))[1, 1, :]))[end-11:end], clim = (18, 24))
