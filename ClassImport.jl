module ClassImport

function get_root()
    return dirname(Base.active_project())
end

const FUNCTION_PATH = joinpath(get_root(), "src", "functions")
const HDF5_PATH = joinpath(get_root(), "data", "hdf5")
const JSON_PATH = joinpath(get_root(), "data", "json")

include(joinpath(FUNCTION_PATH, "GreedyRegression.jl"))
include(joinpath(FUNCTION_PATH, "Integrate.jl"))
include(joinpath(FUNCTION_PATH, "ReportModel.jl"))
include(joinpath(FUNCTION_PATH, "PickSubdomains.jl"))
include(joinpath(FUNCTION_PATH, "EnvelopePol.jl"))
include(joinpath(FUNCTION_PATH, "OpEval.jl"))

end