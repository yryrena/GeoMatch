# GeoMatch.jl

Minimal, working proof-of-concept for spatial matching in Julia.

Given a **treated** table and a **control** table with point geometries (WGS84 lon-lat), GeoMatch builds matched pairs or weighted sets based on great-circle distance, with optional covariate balancing diagnostics.

**MVP scope**: haversine distance, K-nearest neighbors or kernel matching, basic balance table, and a tiny HTML report.

------

## Features

- Haversine distance on lon-lat coordinates  
- K-nearest neighbor matching with or without replacement  
- Kernel matching with Gaussian or triangular kernels  
- Optional radius caliper and same-region constraint  
- Simple balance diagnostics (standardized mean differences)  
- Minimal HTML report generator

