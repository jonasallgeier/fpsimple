# FPsimple

## Summary
A simplified semi-analytical Matlab model of groundwater flow in floodplain aquifers.

## Documentation
This project is consists of two **classes**, which allow the effortless creation of simplified floodplain models, that can be solved semi-analytically. The two classes are:

  * `fpSimple()` &rarr; abstract description of simplified floodplain aquifer
  * `fpAna()` &rarr; semi-analytical model of simplified floodplain aquifer; inherits from `fpSimple()`

The decision to separate the two classes makes it easier to understand the semi-analytical solution, as all helper methods and property definitions are part of `fpSimple`. It also allows the potential development of other classes inheriting from fpSimple, without inheriting the specifics of the semi-analytical solution (e.g., this can be used to build a numerical model class `fpNum()` for any numerical modeling software).

The classes are **documented** by means of Matlab's `help` and `doc` utilities. To explore the functionality of `fpAna()` for example, use the following commands:
```matlab
  fpAna.help        % or "help fpAna" 
  fpAna.doc         % or "doc fpAna"
  methods(fpAna)
  properties(fpAna)
```

Tab completion is another way to access a list of available properties and methods (or sometimes even allowed values of a property). For example, type "`fpAna.`&#8633;" or "`fpAna("shape","`&#8633;".

You can also operate the ` help` and `doc` functions in individual **properties**/**methods**:
```matlab
  % investigate a property
  help fpAna.L
  doc fpAna.L
  % investigate a method
  help fpAna.plot
  doc fpAna.plot
```

## Usage

In quick summary, you can create model entities simply by calling the class **constructor**:
```matlab
    myModel = fpAna(); % create a model with default properties
```

The constructor can also be called with arguments that adjust the model geometry/parameters:
```matlab
    myModel = fpAna("L",500,"Tx",1e-5,"shape","composite"); % create a custom model
```

You can then **modify**, for example the model geometry:

```matlab
    myModel.L = 1200; % change domain length of an existing model
```
...or similarly the other parameters:
```matlab
    myModel.Tx(1) = 1e-5; % change transmissivity in x-direction of existing model
```

To investigate quantitative model outcomes, access the respective properties:
```matlab
    exchangeFlux = myModel.Qex; % get hyporheic exchange flux of model
```

Finally, results can be visualized by using the `plot` method:
```matlab
    plot(myModel); % myModel.plot() also works
```

For more details, please consult the Matlab `help` and `doc` commands.

## Peculiarities
To make the interaction with the models as smooth and user-friendly as possible, the classes include some peculiarities.
* **automatic solution**:
  Any change to a property of a model that has an effect on the flow field automatically triggers the automated re-evaluation of the semi-analytical solution, if the `autoSolve` property is set to `true`.  If it is set to `false`, an update of the solution can be triggered by the method `solve`.
* **automatic plot updating**:
  If you use the `plot()` command to plot a model, the resulting plot will automatically be updated if you change properties of the model, if the `autoUpdate` property is set to `true` and the plot is still *active* (i.e., its axis equals `gca` and it was not deleted).
* **figure storage**:
  The `plot()` command produces several plots at once that are superimposed (e.g., an outline of the geometry and a flow net). To avoid confusion of these components, it does *not* provide any output handle. Instead, all handles are stored with a *description identifier* in the `stor` property. It has the fields `go` (*graphics objects*) and `id` (the description). If you want to modify the plot's components, use this storage.

## Figures Directory
The `figures` directory contains Matlab functions/scripts to (re-)create all auto-generated figures of the main manuscript and supporting information. The file names correspond with figure number in the following way:
* `figArea.m` &rarr; Figure 6
* `figExample.m` &rarr; Figure 8
* `figFlowNetExamples.m` &rarr; Figure 3
* `figGeometry.m` &rarr; Figure 4
* `figInfluxNorth.m` &rarr; Figure 5
* `figSuppQnorth.m` &rarr; Figure S1 (supporting information)
* `figTravelTimes.m` &rarr; Figure 7 (+ Figure S2 in supporting information)

All other files in this directory are helper functions or data files.
One of the scripts (`figExample.m`) requires `pdftk` (https://en.wikipedia.org/wiki/PDFtk) to be installed.


## Troubleshooting

This version of the code is developed & tested with **Matlab R2021a** on **Ubuntu 20.04**. Consider updating, if you are operating an older version of Matlab and errors/warnings appear. If you have further problems or questions, feel free to contact the author of this code.