"""Default OCI images used by the packaged CADD-SV workflow."""

DEFAULT_CONTAINER_IMAGES = {
    "preprocessing": "docker://docker.io/ocatona/cadd-sv-envs:preprocessing-2.0",
    "sv": "docker://docker.io/ocatona/cadd-sv-envs:sv-2.0",
    "nt": "docker://docker.io/ocatona/cadd-sv-envs:nt-2.0",
    "training": "docker://docker.io/ocatona/cadd-sv-envs:training-2.0",
}

COORDINATE_BASED_CONTAINER_ENVIRONMENTS = (
    "preprocessing",
    "sv",
    "training",
)
