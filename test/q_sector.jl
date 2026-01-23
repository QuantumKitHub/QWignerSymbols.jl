import TensorKitSectors
testsuite_path = joinpath(
    dirname(dirname(pathof(TensorKitSectors))), # TensorKitSectors root
    "test", "testsuite.jl"
)
include(testsuite_path)

sectorlist = (SU2qIrrep{1.11}, (SU2kIrrep{level} for level in 2:2)...)

for sector in sectorlist
    SectorTestSuite.test_sector(sector)
end
