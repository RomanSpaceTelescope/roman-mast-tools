roman-mast-tools
Copyright 2026 United States Government as represented by the
Administrator of the National Aeronautics and Space Administration.
All Rights Reserved.

Developed at NASA Goddard Space Flight Center (GSFC).

This project is a lightweight collection of Python utilities that wrap
publicly available Roman Space Telescope data services (MAST archive,
MAST Engineering Database, and public S3 buckets under s3://stpubdata/).
It does not include or redistribute any restricted or proprietary data.

--------------------------------------------------------------------------
Third-Party Software
--------------------------------------------------------------------------

This project depends on, but does not redistribute, the following
open-source packages. Each is used under its own license; see the
respective project for details.

  numpy, scipy              (BSD-3-Clause)
  astropy, photutils        (BSD-3-Clause)
  matplotlib                (Matplotlib License, BSD-style)
  pandas, pyarrow           (BSD-3-Clause / Apache-2.0)
  roman_datamodels,
      romancal, rad, gwcs,
      asdf                  (BSD-3-Clause; Space Telescope Science
                             Institute)
  astroquery                (BSD-3-Clause)
  s3fs, fsspec, requests    (BSD-3-Clause / Apache-2.0)
  pyyaml                    (MIT)
  keyring, python-dotenv    (MIT / BSD-3-Clause)

--------------------------------------------------------------------------
Acknowledgment of AI-Assisted Development
--------------------------------------------------------------------------

Portions of this software, including code scaffolding, documentation,
and README content, were developed with the assistance of Anthropic's
Claude large language model. All AI-generated content was reviewed,
edited, and tested by the human author(s) prior to release. The human