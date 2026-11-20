# Recommendation ITU-R P.528-5 - U.S. Reference Implementation #

Persistent Identifier: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5081071.svg)](https://doi.org/10.5281/zenodo.5081071)

This code repository contains the U.S. Reference Software Implementation of Recommendation ITU-R P.528. This Recommendation contains a method for predicting basic transmission loss in the frequency range 100-30 000 MHz for aeronautical services.  The software implements Annex 2 of Rec P.528, the Step-by-Step method to computing propagation loss for air-to-ground paths.

Companion software ([p528-gui](https://github.com/NTIA/p528-gui)) provides a Graphical User Interface (GUI) that can be used with this software implementation. 

## About This Fork ##

This repository is a fork of [NTIA/p528](https://github.com/NTIA/p528). It adds a much faster way to compute whole curves, a WebAssembly build, and SIMD experiments. It is not an official ITU or NTIA distribution; results are informative.

**Try it in the browser: [miro77.github.io](https://miro77.github.io/)** – the P.528 Explorer. Sliders set frequency, terminal heights, time percentage, radiated power, distance range and the number of points, and the field strength (or basic transmission loss) curve follows while a slider moves. The model runs locally as WebAssembly; a 401-point curve takes about 15 ms.

### Branches

The work lives on branches. `master` is the upstream state at the time of forking plus this README.

| Branch | Content |
|--------|---------|
| [`feature/profile-performance`](https://github.com/miro77/p528-wasm/tree/feature/profile-performance) | `P528_Profile()`: many distances in one call, plus the `P528Bench` benchmark and unit tests. Based on upstream's CMake branch (`cmake_basic`). The web page runs this code, compiled to WebAssembly through a small wrapper that is not in the repository yet. |
| [`cmake_simd`](https://github.com/miro77/p528-wasm/tree/cmake_simd) | Upstream's CMake build merged with SIMD versions of the P.676 line-by-line calculation (`std::experimental::simd`). SIMD is active with GCC/libstdc++; other toolchains such as MSVC use the scalar code. Keeps the P.528 test data in the repository. |
| [`develop-miro`](https://github.com/miro77/p528-wasm/tree/develop-miro) | The earlier Vue web GUI with a multi-threaded WebAssembly build. Superseded by the P.528 Explorer. |

### Computing a Whole Curve: `P528_Profile()`

`P528()` spends almost all of its time on work that does not depend on the path distance, mainly the atmospheric ray traces of Recommendation ITU-R P.676. `P528_Profile()` evaluates many distances that share the terminal heights, frequency, polarization and time percentage, and does that work once:

```cpp
std::vector<double> d__km = { /* distances, in any order */ };
std::vector<Result> results(d__km.size());
std::vector<int> return_codes(d__km.size());

int rtn = P528_Profile(d__km.data(), (int)d__km.size(), h_1__meter, h_2__meter,
                       f__mhz, T_pol, time, results.data(), return_codes.data());
```

Every point gets the result and return code `P528()` would give for it.

 * Terminal geometry, the smooth earth diffraction line, the line-of-sight setup, the K-value at the radio horizon and the diffraction/troposcatter crossover are computed once per call. These values are identical to those of `P528()`.
 * Line-of-sight rays share the atmospheric layers between the two terminals.
 * For rays whose layers change with the distance, the specific gaseous attenuation and the refractivity are interpolated versus height from a table built on demand for the frequency (piecewise Chebyshev with error control), instead of being summed line by line for every layer of every ray.

Measured with the benchmark of that branch, five scenarios with 1 000 distances each, on one desktop machine:

| Build | `P528()` per point | `P528_Profile()` | Speed-up | Largest deviation |
|-------|--------------------|------------------|----------|-------------------|
| MSVC, x64 | 73.7 s | 0.25 s | about 290x | 1.4e-9 dB |

For large profiles the remaining cost is roughly 10 to 40 microseconds per point, depending on compiler and scenario. `P528()` itself returns bit-identical results to upstream on that branch.

The benchmark compares both functions on every run and checks them against stored reference results, see the README of the branch:

```
P528Bench --points 1000
P528Bench --points 100000 --profile-only
P528Bench --points 200 --check-reference bench/reference_200.txt
```

### Does SIMD Help?

Measured with GCC 13 on a machine with AVX-512, on top of `P528_Profile()`:

 * Single calls of `P528()` become about 1.75x faster with the SIMD line-by-line code of `cmake_simd`.
 * Short profiles (300 points) gain about 1.4x, because the fixed setup per profile gets shorter.
 * Large profiles (20 000 points) gain about 1 %: they interpolate from the table and hardly run the line-by-line code.

### Building the Branches

Both CMake branches build with the presets in their root directory:

```
git submodule update --init extern/googletest
cmake --preset release64 -DBUILD_DOCS=OFF
cmake --build --preset release64
ctest --preset release64
```

`-DBUILD_DOCS=OFF` skips the Doxygen site. The upstream test data submodule (`extern/test-data`) is not publicly reachable; `cmake_simd` therefore carries its own copy of the data, while on `feature/profile-performance` the two tests that need it fail.

## Inputs ##

| Variable          | Type   | Units | Limits       | Description  |
|-------------------|--------|-------|--------------|--------------|
| `d__km`           | double | km    | 0 <= `d__km` | Great circle path distance between terminals |
| `h_1__meter`      | double | meter | 1.5 <= `h_1__meter` <= 20 000 | Height of the low terminal |
| `h_2__meter`      | double | meter | 1.5 <= `h_2__meter` <= 20 000 | Height of the high terminal |
| `f__mhz`          | double | MHz   | 100 <= `f__mhz` <= 30 000 | Frequency |
| `T_pol`           | int    |       |              | Polarization <ul><li>0 = Horizontal</li><li>1 = Vertical</li></ul> |
| `time`            | double |       | 1 <= `time` <= 99 | Time percentage |
 
## Outputs ##

Outputs to P.528 are contained within a defined `Results` structure.

| Variable   | Type   | Units | Description |
|------------|--------|-------|-------------|
| `d__km`    | double | km    | Great circle path distance.  Could be slightly different than specified in input variable if within LOS region |
| `A__db`    | double | dB    | Basic transmission loss |
| `A_fs__db` | double | dB    | Free space basic transmission loss |
| `A_a__db`  | double | dB    | Median atmospheric absorption loss |
| `theta_h1__rad` | double | rad | Elevation angle of the ray at the low terminal |
| `propagation_mode` | int |  | Mode of propagation <ul><li>1 = Line of Sight</li><li>2 = Diffraction</li><li>3 = Troposcatter</li></ul> |
| `warnings` | int    |       | Warning flags |

## Error Codes and Warning Flags ##

P.528 supports a defined list of error codes and warning flags.  A complete list can be found [here](ERRORS_AND_WARNINGS.md).

## Example Values ##

The below table includes a select set of example inputs and outputs for testing purposes. The [P.528 CSV data files](https://www.itu.int/rec/R-REC-P.528/en) contain an extensive set of validation example values.

| `d__km` | `h_1__meter` | `h_2__meter` | `f__mhz` | `T_pol` | `time` | `A__db` |
| --------|--------------|--------------|----------|---------|--------|---------|
|      15 |           10 |        1 000 |      500 |       0 |     50 |   110.0 |
|     100 |          100 |       15 000 |    3 600 |       0 |     90 |   151.6 |
|   1 500 |           15 |       10 000 |    5 700 |       0 |     10 |   293.4 |
|      30 |            8 |       20 000 |   22 000 |       1 |     50 |   151.1 |

## Notes on Code Style ##

 * In general, variables follow the naming convention in which a single underscore denotes a subscript (pseudo-LaTeX format), where a double underscore is followed by the units, i.e. h_1__meter.
 * Variables are named to match their corresponding mathematical variables in the underlying Recommendation text.
 * Wherever possible, equation numbers are provided.  It is assumed that a user reviewing this source code would have a copy of the Recommendation's text available as a primary reference.

## Configure and Build ##

### C++ Software

The software is designed to be built into a DLL (or corresponding library for non-Windows systems).  The source code can be built for any OS that supports the standard C++ libraries.  A Visual Studio 2019 project file is provided for Windows users to support the build process and configuration.

### C#/.NET Wrapper Software

The .NET support of P.528 consists of a simple pass-through wrapper around the native DLL.  It is compiled to target .NET Framework 4.8.  Distribution and updates are provided through the published [NuGet package](https://github.com/NTIA/p528/packages).

## References ##

 * [Recommendation ITU-R P.528](https://www.itu.int/rec/R-REC-P.528/en)
 * [Rec ITU-R P.528 GUI](https://github.com/NTIA/p528-gui)
 * [Upstream repository NTIA/p528](https://github.com/NTIA/p528)
 * [P.528 Explorer](https://miro77.github.io/), the web page of this fork

## Contact ##

The fork-specific parts described under "About This Fork" are not maintained by NTIA.

For questions about the reference implementation itself, contact the upstream maintainers: Billy Kozma, (303) 497-6082, wkozma@ntia.gov
