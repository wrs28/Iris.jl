push!(LOAD_PATH,"../src/")

using Documenter
using Iris

makedocs(
    sitename = "Iris.jl",
    doctest = false,
    modules = [Iris, Iris.Common, Iris.Spectral],
    pages = [
            "index.md",
            "Manual" => [
                "Building a simulation" => "Common.md",
                "Spectral Analysis" => "Spectral.md",
                "Band dispersion (under construction)" => "Floquet.md",
                "Scattering" => "Scattering.md",
                "S-Matrix" => "SMatrices.md",
                "Lasing" => "Lasing.md",
                "Saturating CPA" => "SaturableCPA.md",
                "FDTD" => "TimeDomain.md",
                ],
            "Examples.md",
            "Library" => [
                "Library/Iris.md",
                "Common" => [
                    "Library/Common/Common.md",
                    "Library/Common/Points.md",
                    "Library/Common/ElectricFields.md",
                    "Library/Common/Shapes.md",
                    "Library/Common/BoundaryLayers.md",
                    "Library/Common/BoundaryConditions.md",
                    "Library/Common/Boundaries.md",
                    "Library/Common/DielectricFunctions.md",
                    "Library/Common/PumpFunctions.md",
                    "Library/Common/Dispersions.md",
                    "Library/Dispersions" => [
                        "Library/Common/Dispersions/Kerrs.md",
                        "Library/Common/Dispersions/TwoLevelSystems.md",
                        ],
                    "Library/Common/Lattices.md",
                    "Library/Common/Domains.md",
                    "Library/Common/Curlcurls.md",
                    "Library/Common/SelfEnergies.md",
                    "Library/Common/Simulations.md",
                    "Library/Common/LU_Factorizations.md",
                    "Library/Common/Plotting.md",
                    ],
                "Library/Spectral.md",
                "Library/Floquet.md",
                "Library/Scattering.md",
                "Library/SMatrices.md",
                "Library/Lasing.md",
                "Library/SaturableCPA.md",
                "Library/TimeDomain.md"
                ]
            ]
    )
