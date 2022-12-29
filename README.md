# CFL_Blogpost
This repo contains the accompanying code to our blogpost about "Causal Feature Learning" by Chalupka et al.

TODO: add the link to the blogpost when online


You can run the code by cloning this repo locally and launching Julia within the folder. Afterwards, you can run the following commands to install the necessary packages:

```julia
(YOUR_JULIA_VERSION) pkg> activate .
  Activating project at `YOUR_REPO_PATH/CFL_Blogpost`

(CFL_Blogpost) pkg> instantiate
```

`pkg>` is a special REPL mode you can enter from the julia REPL by pressing `]`.

------

The blogpost and the code in this repo (although written from scratch) are mainly based upon the following publications:

Chalupka, K., Eberhardt, F., Perona, P., 2017. Causal feature learning: an overview. Behaviormetrika 44, 137–164. https://doi.org/10.1007/s41237-016-0008-2
Chalupka, K., Bischoff, T., Perona, P., Eberhardt, F., 2016. Unsupervised Discovery of El Nino Using Causal Feature Learning on Microlevel Climate Data.
Chalupka, K., Perona, P., Eberhardt, F., 2015a. Visual Causal Feature Learning.
Chalupka, K., Perona, P., Eberhardt, F., 2015b. Multi-Level Cause-Effect Systems.
