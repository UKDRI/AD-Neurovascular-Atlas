#!/bin/bash

# Variables
VERSION="1.0.0"
GITHUB_REPO="UKDRI/AD-Neurovascular-Atlas"
ZENODO_DOI="10.5281/zenodo.17737332"
BUILD_DATE=$(date -u +'%Y-%m-%dT%H:%M:%SZ')
GIT_COMMIT=$(git rev-parse --short HEAD)

# Build image
docker build \
  --build-arg VERSION="${VERSION}" \
  --build-arg BUILD_DATE="${BUILD_DATE}" \
  --build-arg GIT_COMMIT="${GIT_COMMIT}" \
  --build-arg GITHUB_REPO="${GITHUB_REPO}" \
  --build-arg ZENODO_DOI="${ZENODO_DOI}" \
  -t ad-bbb-analysis:${VERSION} \
  -t ad-bbb-analysis:latest \
  .

echo "Built image: ad-bbb-analysis:${VERSION}"
