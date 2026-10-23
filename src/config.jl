# SPDX-License-Identifier: LGPL-3.0-or-later
#
# Run configuration, read from the environment so that a batch script needs no edits:
#
# | variable            | default  | meaning                                             |
# |---------------------|----------|-----------------------------------------------------|
# | `AFRICA_SCALE`      | `large`  | `small` for a smoke test: coarsest grid, few species; `medium` for a trial: coarsest grid, a sample of the species, the full span |
# | `AFRICA_CELLSIZE`   | chosen   | `50km` to override the memory-based choice           |
#
# Read at call time, not at load time, so a test can set them after `using AfricaPlants`.

"""
    smallscale()

Return `true` when `AFRICA_SCALE` is `small`, the smoke-test setting: the coarsest candidate
grid and a small species pool, so that a run finishes in seconds.
"""
smallscale() = get(ENV, "AFRICA_SCALE", "large") == "small"

"""
    mediumscale()

Return `true` when `AFRICA_SCALE` is `medium`, the setting for trying a change to the model: the
coarsest candidate grid and a sample of the species drawn across every level of recording, run over
the full span of the climate, so that a whole experiment finishes in minutes and still holds the
rarely recorded species a change may bear on.
"""
mediumscale() = get(ENV, "AFRICA_SCALE", "large") == "medium"

"""
    SCALED_SPECIES

How many species a run keeps at each reduced scale.
"""
const SCALED_SPECIES = (small = 200, medium = 3000)

"""
    CANDIDATE_CELLSIZES

The grid resolutions a run may choose between, finest first. What each costs is asked of
EcoSISTEM at run time (see [`perspecies_bytes`](@ref)), never written down here.
"""
const CANDIDATE_CELLSIZES = (5.0km, 10.0km, 20.0km, 50.0km, 100.0km)

"""
    requested_cellsize()

Return the cell size named by `AFRICA_CELLSIZE`, one of [`CANDIDATE_CELLSIZES`](@ref) written
as `"50km"` or `"50.0 km"`, or `nothing` when the variable is unset and the size is to be chosen
from the memory available. A value outside the candidates is an error.
"""
function requested_cellsize()
    want = get(ENV, "AFRICA_CELLSIZE", "")
    isempty(want) && return nothing
    for cs in CANDIDATE_CELLSIZES
        (string(cs) == want || "$(Int(ustrip(cs)))km" == want) && return cs
    end
    return error("AFRICA_CELLSIZE=$want is not one of " *
                 join(CANDIDATE_CELLSIZES, ", "))
end
