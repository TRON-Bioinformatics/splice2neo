#!/bin/bash
#
# Maintainer: Johannes Hausmann <johannes.hausmann@posteo.de>
# Build OCI images for splice2neo

cp ../DESCRIPTION DESCRIPTION

export TAG=$(grep 'Version:' DESCRIPTION | awk '{print $2}')

podman build -t tronbioinformatics/splice2neo:"${TAG}" .
#podman build -t tronbioinformatics/splice2neo:latest .

rm DESCRIPTION
