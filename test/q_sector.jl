import TensorKitSectors

testsuite_path = joinpath(
    dirname(dirname(pathof(TensorKitSectors))), # TensorKitSectors root
    "test", "testsuite.jl"
)
include(testsuite_path)

sectorlist = (
    SU2qIrrep{1.11}, SU2qIrrep{1.3}, SU2kIrrep{2}, SU2kIrrep{3}, SU2kIrrep{4},
    SU2kIrrep{5}, SU2kIrrep{6}, SU2kIrrep{7}, SU2kIrrep{8}, SU2kIrrep{9}, SU2kIrrep{10},
    # SU2qIrrep{convert(ComplexF64, RootOfUnity{3}())}, SU2qIrrep{1.42 + 0.4im}
)

for sector in sectorlist
    SectorTestSuite.test_sector(sector)
end
