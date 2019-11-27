push!(LOAD_PATH,"../src/")

using Documenter
using Iris

makedocs(
    sitename = "Iris.jl",
    doctest = false,
    modules = [Iris, Iris.Common],
    pages = [
            "index.md",
            "Building a simulation" => "Building.md",
            "Library" => [
                "Iris.md",
                "Common" => [
                    "Common/Common.md",
                    "Common/Points.md",
                    "Common/ElectricFields.md",
                    "Common/Shapes.md",
                    "Common/BoundaryLayers.md",
                    "Common/BoundaryConditions.md",
                    "Common/Boundaries.md",
                    "Common/DielectricFunctions.md",
                    "Common/PumpFunctions.md",
                    "Common/Dispersions.md",
                    "Dispersions" => [
                        "Common/Dispersions/Kerrs.md",
                        "Common/Dispersions/TwoLevelSystems.md",
                        ],
                    "Common/Lattices.md",
                    "Common/Domains.md",
                    "Common/Curlcurls.md",
                    "Common/SelfEnergies.md",
                    "Common/Simulations.md",
                    "Common/LU_Factorizations.md",
                    "Common/Plotting.md",
                    ],
                "Spectral.md",
                "Floquet.md",
                "Scattering.md",
                "SMatrices.md",
                "Lasing.md",
                "SaturableCPA.md",
                "TimeDomain.md"
                ]
            ]
    )
