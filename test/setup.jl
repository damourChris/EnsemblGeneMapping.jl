# Script to setup the conda environment for the test suite
# This script is called by the test/runtests.jl script

using CondaPkg
using UUIDs
using Preferences
using Libdl

# Setup the channels
CondaPkg.add_channel("conda-forge")
CondaPkg.add_channel("bioconda")

# Install the required packages
CondaPkg.add(["r-base",
              "r-dplyr",
              "bioconductor-biomaRt"])

target_rhome = joinpath(CondaPkg.envdir(), "lib", "R")
if Sys.iswindows()
    target_libr = joinpath(target_rhome, "bin", Sys.WORD_SIZE == 64 ? "x64" : "i386",
                           "R.dll")
else
    target_libr = joinpath(target_rhome, "lib", "libR.$(Libdl.dlext)")
end

const RCALL_UUID = UUID("6f49c342-dc21-5d91-9882-a32aef131414")

# Check if LocalPreferences.toml exists
if !isfile("LocalPreferences.toml")
    # Create a LocalPreferences.toml file
    open("LocalPreferences.toml", "w") do io
        return write(io, """
               [RCall]
               Rhome = "$target_rhome"
               libR = "$target_libr"
               """)
    end
else
    set_preferences!(RCALL_UUID, "Rhome" => target_rhome, "libR" => target_libr;
                     force=true)
end