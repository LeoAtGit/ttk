# Installation

The software works with the following packages:

* Paraview v5.9.0, commit `0aaf903ba3e8d012a608a466f798a041b4aebe17`
* TTK from this repository, branch `streamtree`
* websocketpp 0.8.2

To be able to use the Streamtree code, you have to build websocketpp, TTK and Paraview from scratch.
Build TTK and Paraview with websocket support enabled.
To install and build Paraview with TTK, follow the guide at https://topology-tool-kit.github.io/installation.html .

# Overview

The backend of the application is Paraview/TTK.
There, you have to do the data manipulation, calculate the merge tree, and so on.
The Paraview statefiles are provided (TODO: where?) for the use cases in the paper.
The datasets are also provided there (TODO: where?).
The last step of the data pipeline in Paraview will always be a Websocket, that sends the data to the frontend.

The frontend can be found in `TTK/standalone/MergeTreeIntegration`.
To use the application, open `KDERenderer.html` in your browser.
The main code for the tree visualization can be found in `TreeRenderer.js`.
The main code for the map can be found in `KDERenderer.js`.

# Different Use Cases

The different use cases require different code, like different labels on the y-axis.
In its current state, the code is set to the Chicago use case.
For the asteroid use case, search in the `TreeRenderer.js` file for `asteroid` and read the comments, to see what you have to change.
