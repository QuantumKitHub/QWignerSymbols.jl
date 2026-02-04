import TensorKitSectors

testsuite_path = joinpath(
    dirname(dirname(pathof(TensorKitSectors))), # TensorKitSectors root
    "test", "testsuite.jl"
)
include(testsuite_path)

sectorlist = (
    SU2qIrrep{1.11}, SU2qIrrep{1.3}, ntuple(x -> SU2qIrrep{RootOfUnity(x + 1)}, 8)...,
)

for sector in sectorlist
    SectorTestSuite.test_sector(sector)
end
