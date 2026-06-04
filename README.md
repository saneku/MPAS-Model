MPAS-v8.4.0
====

The Model for Prediction Across Scales (MPAS) is a collaborative project for
developing atmosphere, ocean, and other earth-system simulation components for
use in climate, regional climate, and weather studies. The primary development
partners are the climate modeling group at Los Alamos National Laboratory
(COSIM) and the National Center for Atmospheric Research. Both primary
partners are responsible for the MPAS framework, operators, and tools common to
the applications; LANL has primary responsibility for the ocean model, and NCAR
has primary responsibility for the atmospheric model.

The MPAS framework facilitates the rapid development and prototyping of models
by providing infrastructure typically required by model developers, including
high-level data types, communication routines, and I/O routines. By using MPAS,
developers can leverage pre-existing code and focus more on development of
their model.

BUILDING
========

This README is provided as a brief introduction to the MPAS framework. It does
not provide details about each specific model, nor does it provide building
instructions.

For information about building and running each core, please refer to each
core's user's guide, which can be found at the following web sites:

[MPAS-Atmosphere](http://mpas-dev.github.io/atmosphere/atmosphere_download.html)

[MPAS-Albany Land Ice](http://mpas-dev.github.io/land_ice/download.html)

[MPAS-Ocean](http://mpas-dev.github.io/ocean/releases.html)

[MPAS-Seaice](http://mpas-dev.github.io/sea_ice/releases.html)


Code Layout
----------

Within the MPAS repository, code is laid out as follows. Sub-directories are
only described below the src directory.

	MPAS-Model
	├── src
	│   ├── driver -- Main driver for MPAS in stand-alone mode (Shared)
	│   ├── external -- External software for MPAS (Shared)
	│   ├── framework -- MPAS Framework (Includes DDT Descriptions, and shared routines. Shared)
	│   ├── operators -- MPAS Opeartors (Includes Operators for MPAS meshes. Shared)
	│   ├── tools -- Empty directory for include files that Registry generates (Shared)
	│   │   ├── registry -- Code for building Registry.xml parser (Shared)
	│   │   └── input_gen -- Code for generating streams and namelist files (Shared)
	│   └── core_* -- Individual model cores.
	│       └── inc -- Empty directory for include files that Registry generates
	├── testing_and_setup -- Tools for setting up configurations and test cases (Shared)
	└── default_inputs -- Copies of default stream and namelists files (Shared)

Model cores are typically developed independently. For information about
building and running a particular core, please refer to that core's user's
guide.


MPAS-Chem / MUSICA Development
------------------------------

This branch includes experimental MPAS-Chem development for coupling the MPAS
atmosphere core with MUSICA chemistry components.

Current chemistry capabilities include:

* MICM chemistry configured with `config_micm_file` in the `&musica` namelist.
* TUV-x photolysis configured with `config_tuvx_config_file`.
* Dynamic emission species selected with `config_chem_emission_species`.
* Extra passive or non-MICM tracers selected with
  `config_chem_additional_species`.
* Hourly anthropogenic emission input through the `chem_emissions` stream,
  using files named like `emi.$Y-$M-$D_$h.$m.$s.nc`.

Emission files should contain fields named `emi_<species>` for the species in
`config_chem_emission_species`. The emission fields are treated as surface or
layer fluxes in `kg m-2 s-1`; during chemistry stepping they are converted to
mixing-ratio increments with the local air-column mass and added to the
matching MPAS chemistry tracers.

Regional chemistry boundary conditions use the standard MPAS limited-area LBC
infrastructure. Dynamic chemistry species are added to the `lbc_scalars`
var-array as `lbc_<species>` fields, and regional runs require compatible
`lbc.*.nc` files containing those fields. Global runs do not use chemistry LBC
files. IC/LBC chemistry fields must be prepared by Merra2BC or another
preprocessor before running MPAS-Chem.

Useful diagnostics include:

* `cos_sza`, the cosine of the solar zenith angle used by chemistry
  photolysis.
* `photolysis_rates`, a dynamic var-array containing TUV-x photolysis rates
  discovered from the selected TUV-x configuration.

TUV-x also returns heating rates through its solver API. These are currently
kept only as internal solver workspace; they are not written as diagnostics and
are not coupled to MPAS thermodynamic or radiation tendencies.

Relevant chemistry modules include:

* `src/core_atmosphere/chemistry/mpas_chem_photolysis_driver.F`
* `src/core_atmosphere/chemistry/mpas_chem_photolysis_tuvx.F`
* `src/core_atmosphere/chemistry/mpas_chem_solar_geometry.F`
* `src/core_atmosphere/chemistry/musica/mpas_musica.F`
