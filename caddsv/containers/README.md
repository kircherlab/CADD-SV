# CADD-SV rule containers

The same Dockerfile builds the four dependency images used by the Snakemake
workflow. It uses the conda-forge Miniforge image and `conda env create` with
Conda's libmamba solver to install each environment in an isolated image
environment. Use the `caddsv` package directory as the Docker build context:

```bash
docker build \
  --file caddsv/containers/Dockerfile \
  --build-arg ENV_FILE=workflow/envs/SV.yml \
  --build-arg CADD_SV_VERSION=2.0 \
  --tag cadd-sv-sv:2.0 \
  caddsv
```

Valid environment files are `preprocessing.yml`, `SV.yml`, `NT.yml`, and
`training.yml`. Release tags publish the corresponding images to GHCR. The
workflow defaults can be overridden through the `containers` mapping in a
CADD-SV config file, including with paths to prebuilt local SIF images.
