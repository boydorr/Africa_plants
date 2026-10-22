# SPDX-License-Identifier: LGPL-3.0-or-later
#
# The whole experiment in one process, on as many threads as it was started with: the backwards
# burn-in and the years of held climate that settle it, then the forward run recording diversity.
# `examples/africa_mpi.jl` is the same experiment across MPI ranks.
#
#     julia -t 12 --project examples/africa_threads.jl
#     AFRICA_SCALE=small julia -t 4 --project examples/africa_threads.jl
#
# The environment variables are `examples/africa_burnin.jl`'s and `examples/africa_forward.jl`'s,
# which this runs in turn; `AFRICA_DISTRIBUTED` is set here and must not be.
#
# Threads are what this file is for: the step loop is threaded, so the run takes about as long as
# one thread's share of it. Started with one thread it says so and runs anyway.

haskey(ENV, "AFRICA_DISTRIBUTED") &&
    error("`examples/africa_threads.jl` sets AFRICA_DISTRIBUTED itself; unset it, or run " *
          "`examples/africa_burnin.jl` and `examples/africa_forward.jl` directly.")
ENV["AFRICA_DISTRIBUTED"] = "serial"

Threads.nthreads() > 1 ||
    @warn "Running on one thread. Start Julia with `-t 12` (this machine's performance cores) for a " *
          "many-fold quicker run."
println("Africa on ", Threads.nthreads(), " thread",
        Threads.nthreads() == 1 ? "" : "s")
flush(stdout)

# Each stage is a script of its own, so each goes in its own module: they name the same constants,
# and a bare `include` of both would redefine them.
module Burnin
include("africa_burnin.jl")
end

module Forward
include("africa_forward.jl")
end
