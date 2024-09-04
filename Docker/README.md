# Container support

We provide a Singularity and soon a Docker image that contains all software pre-installed for running a splice2neo analysis.
However, singularity is our preferred method for running splice2neo and allows unprivileged users to run the application in a container.

To use the the singularity image, you have to first build. Please make sure you have `fakeroot` permissions or ask your admin to build the image.

```
singularity build --fakeroot splice2neo.sif splice2neo.def
```


To run an analysis, bind your input directory in the container and run the analysis code either as script or interactively.

### Script usage

```
singularity exec -e -B `pwd` -B /path/to/your/inputs splice2neo.sif Rscript splice2neo.R
```

### Interactive usage

```
singularity shell -B `pwd` -B /path/to/your/inputs splice2neo.sif
```

Within the container you can now open a R-Session and execute the commands.
