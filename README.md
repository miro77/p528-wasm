# Recommendation ITU-R P.528-5 - U.S. Reference Implementation #

Persistent Identifier: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5081071.svg)](https://doi.org/10.5281/zenodo.5081071)

This code repository contains the U.S. Reference Software Implementation of Recommendation ITU-R P.528. This Recommendation contains a method for predicting basic transmission loss in the frequency range 100-30 000 MHz for aeronautical services.  The software implements Annex 2 of Rec P.528, the Step-by-Step method to computing propagation loss for air-to-ground paths.

Companion software ([p528-gui](https://github.com/NTIA/p528-gui)) provides a Graphical User Interface (GUI) that can be used with this software implementation. 

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

## Profiles: Many Distances at Once ##

`P528()` spends almost all of its time on quantities that do not depend on the path distance, above all on the atmospheric ray traces of Recommendation ITU-R P.676. When the loss is needed for many distances that share the terminal heights, frequency, polarization and time percentage (the points of a profile), use `P528_Profile()` instead of calling `P528()` per point:

```cpp
std::vector<double> d__km = { /* distances along the profile, in any order */ };
std::vector<Result> results(d__km.size());
std::vector<int> return_codes(d__km.size());

int rtn = P528_Profile(d__km.data(), (int)d__km.size(), h_1__meter, h_2__meter,
                       f__mhz, T_pol, time, results.data(), return_codes.data());
```

Every point gets the result and return code `P528()` would give for it. `return_codes` may be null; the function's own return value is `SUCCESS`, the first error code, or `SUCCESS_WITH_WARNINGS`.

The function is faster for these reasons:

 * The terminal geometries, the smooth earth diffraction line, the line-of-sight setup, the K-value at the radio horizon and the diffraction/troposcatter crossover are computed once per profile. The geometry is identical to that of `P528()`, bit for bit: it decides where line-of-sight ends and which steps the searches of the line-of-sight region take.
 * The specific gaseous attenuation is interpolated versus height from a table that is built on demand for the frequency, instead of being summed line by line for every layer of every ray.
 * Beyond the horizon the ray to the common volume is interpolated instead of traced. Its loss is not a smooth function of the height of the common volume, because the number of layers changes with it. But for a given layer scale all these rays share their layers, and the scale varies by a few percent only, so the sums up to each layer are tabulated versus the scale.
 * Upward-looking line-of-sight rays share the layers between the terminals. Their loss and length are interpolated versus the elevation angle.
 * Downward-looking line-of-sight rays cannot be interpolated: their grazing height is the result of a search with a 1 m stopping criterion, and the loss jumps between its possible values. The search and the rays of each grazing height are remembered instead. The same holds for the search of the reflection angle of a distance.

Every table is verified when it is built, and is only built once enough distances use it. In the benchmark the results stay within 1e-8 dB of `P528()`, with the same propagation mode, warnings and return code for every point.

`P528_Profile()` runs on the calling thread and may be called from several threads at once, each call keeps its own state. For large profiles `P528_ProfileParallel()` takes one more argument, the largest number of threads to use (0: one per processor), and evaluates contiguous ranges of the distances at the same time. It pays off from about 4000 distances; smaller profiles and WebAssembly builds are evaluated on the calling thread.

```cpp
int rtn = P528_ProfileParallel(d__km.data(), (int)d__km.size(), h_1__meter, h_2__meter,
                               f__mhz, T_pol, time, results.data(), return_codes.data(), 0);
```

### Measuring the Calculation Time

The `P528Bench` executable (CMake option `BUILD_BENCHMARK`, on by default) times a profile for several scenarios with both functions and compares the results:

```
P528Bench --points 1000            # P528() per point versus P528_Profile()
P528Bench --points 100000 --profile-only   # large profiles, P528_Profile() only
P528Bench --points 100000 --profile-only --threads 0   # P528_ProfileParallel()
P528Bench --points 2 --profile-only --repeat 5   # the fixed cost of a profile
P528Bench --points 50 --stages     # additionally time the building blocks
P528Bench --points 200 --check-reference bench/reference_200.txt
```

It exits with a non-zero code if the two functions differ by more than 1e-6 dB, or disagree on a propagation mode, warning or return code. `bench/reference_200.txt` holds the results of `P528()` from before `P528_Profile()` was added; checking against it is part of the test suite and guards changes that affect both functions. Use `--save-reference` to store new reference results, and `--csv` to log timings.

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

## Contact ##

For questions, contact Billy Kozma, (303) 497-6082, wkozma@ntia.gov
