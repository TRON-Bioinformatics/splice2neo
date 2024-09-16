# Container support

We provide a Singularity and a Docker image that contains all software pre-installed for running the splice2neo analysis.
However, singularity is our preferred method for running splice2neo as it allows unprivileged users to run the application in a container.


## Singularity

To use the the singularity image, you have to first build it from the definition file. Please make sure you have `fakeroot` permissions or ask your sysadmin to build the image.

```
singularity build --fakeroot splice2neo.sif splice2neo.def
```


To run the analysis, bind your input directory into the container. The container ships the HG19 version of BSgenome.

### Script usage

```
singularity exec -e -B `pwd` -B /path/to/your/inputs splice2neo.sif Rscript splice2neo.R
```

### Interactive usage

```
singularity shell -B `pwd` -B /path/to/your/inputs splice2neo.sif
```

Within the container you can now open a R-Session and execute the commands.


## Docker / podman

We use podman to build and manage our Docker images. Podman is API compatible with Docker and provides unprivileged execution of images.


```
podman build -t tronbioinformatics/splice2neo:"${TAG}" .
```
