# FPsimple - Python

## Summary
This is the interactive Python Bokeh variant of the simplified floodplain model.

## Usage
Run `fpsimple.py` to create an HTML file called `fpsimple.html`. This is a standalone interactive model that can be directly executed in the browser. If you want to move this file to another place without losing its functionality, you need to be careful to copy the following files and directories too, as they are linked:
* `fpsimple.js` (javascript helper functions)
* `marchingsquares.min.js` (javascript function to generate contour lines)
* `themodel.js` (javascript callback for the model)
* `static` (directory with static website files)

The `bokeh` command is part of the Python Bokeh library, which can be found at https://bokeh.org/.

## Troubleshooting
This version of the code is developed & tested with **Mozilla Firefox 89.0.1**, **Bokeh 2.3.2** and **Python 3.8.5** on **Ubuntu 20.04**. In case you encounter errors or warnings or have further problems or questions, feel free to contact the author of this code.