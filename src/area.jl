# SPDX-License-Identifier: LGPL-3.0-or-later
#
# The continent as a study area: where it is, which projection it is simulated in, and the
# finest grid the memory available can hold.

"""
    AFRICA_CRS

The projected coordinate system every Africa run is simulated in: EPSG:10592, WGS 84 / GLANCE
Africa, an equal-area projection for the continent and EcoSISTEM's own suggestion for this
extent. A simulation needs a projected grid because dispersal assumes one uniform cell size.
"""
const AFRICA_CRS = EPSG(10592)

"""
    africa_within(; outline = true, coverage = AllTerritories())

Return the `within` specification positioning a study area over Africa: Natural Earth's
continent-level "Africa", downloaded the first time it is built, activating every cell the outline
touches at all ([`AnyOverlap`](@ref)), so that a coastal cell holding any land is simulated and the
available-ground supply scales it to that land; or with `outline = false` its bounding box alone at
no download.

`coverage` says how much of the name to take. The default keeps every territory, Madagascar
and the islands included; `LargestLandmass()` keeps the mainland only.
"""
function africa_within(; outline::Bool = true, coverage = AllTerritories())
    outline || return EcoSISTEM.boundingbox("Africa", level = "CONTINENT",
                                 coverage = coverage)
    return ShapeMaskSpec(NaturalEarthSpec("Africa", level = "CONTINENT",
                                          coverage = coverage),
                         AnyOverlap())
end

"""
    perspecies_bytes(cellsize; regime, supply, within = africa_within(), crs = AFRICA_CRS)

Return the bytes one species' abundances would occupy on the Africa grid at `cellsize`,
over the whole grid, active cells or not. Asked of EcoSISTEM's own grid analysis
(`investigate_study_area`, which resolves the grid without building it), so it cannot disagree
with what a run allocates. `regime` and `supply` are the layer specifications that will shape
the grid.
"""
function perspecies_bytes(cellsize::Unitful.Length; regime, supply,
                          within = africa_within(), crs = AFRICA_CRS)
    report = investigate_study_area(regime = regime, supply = supply,
                                    within = within, crs = crs,
                                    cellsize = cellsize)
    return EcoSISTEM.getspeciesstorage(report)
end

"""
    choose_cellsize(numspecies; regime, supply, within = africa_within(), crs = AFRICA_CRS,
                    candidates = CANDIDATE_CELLSIZES, budget = MemoryGuidance.memory_budget())

Return the finest of `candidates` that `numspecies` species fit in, as a named tuple of
`cellsize`, `perspecies` bytes, the `required` and `usable` bytes and whether it `fits`.

An explicit `AFRICA_CELLSIZE` is honoured, and checked through the same path so that one that
does not fit is reported rather than discovered by an out-of-memory kill; a small-scale run
(`AFRICA_SCALE=small` or `medium`) takes the coarsest candidate outright. Under MPI `budget` is the sum
over the nodes the job landed on, and its first measurement is collective: reach this on every
rank.
"""
function choose_cellsize(numspecies::Integer; regime, supply,
                         within = africa_within(), crs = AFRICA_CRS,
                         candidates = CANDIDATE_CELLSIZES,
                         budget = MemoryGuidance.memory_budget())
    want = requested_cellsize()
    pool = (smallscale() || mediumscale()) ? (maximum(candidates),) :
           isnothing(want) ? candidates : (want,)
    return MemoryGuidance.choose_cellsize(numspecies, pool,
                                          budget = budget) do cellsize
        return perspecies_bytes(cellsize, regime = regime, supply = supply,
                                within = within, crs = crs)
    end
end

"""
    africa_area(cellsize; regime, supply, within = africa_within(), crs = AFRICA_CRS,
                verbosity = :normal)

Return the [`StudyArea`](@ref) for Africa at `cellsize`: the same inputs
[`perspecies_bytes`](@ref) priced, so the grid built is the one that was costed.
"""
function africa_area(cellsize::Unitful.Length; regime, supply,
                     within = africa_within(), crs = AFRICA_CRS,
                     verbosity::Symbol = :normal)
    return StudyArea(regime = regime, supply = supply, within = within,
                     crs = crs, cellsize = cellsize, verbosity = verbosity)
end
