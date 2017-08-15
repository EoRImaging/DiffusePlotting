## Synopsis

The purpose of these plotting modules is to be able to analyze the linearly polarized components
of a set of Stokes data files, I, Q, and U.

## Functions

The run_data module is a wrapper that allows you to create multiple kinds of plots of the same
data.  It expects to be given fits files.

```
run_data.run_data(filename_Q, filename_U, filename_I, plot_variable=None,
             save_show='show', graph_selection='all', file_extension='.png', projection='cyl', transparency=1
             histogram_file_basename=None, map_file_basename=None, drapery_file_basename=None, directory=None)
```

  Filename_Q, filename_U, filename_I:  3 files of a typical 4pol data set.  They are configured to
  use .fits files.

  plot_variable: the variable to be plotted for the plothealpix_map module.  The other two plotting
  modules plot specific variables.  Default is None.

  save_show: Options are 'save' and 'show'. 'Save' option saves the plots you've chosen, as whatever type
  of file has been specified.

  graph_selection: Options are 'map', 'stokes_histogram', 'drapery', and 'all'.

  file_extension: Type of file to be saved. Not necessary if save_show option 'show' is selected.
  Default is .png.

  projection: The map projection for the mapping module plothealpix_map. The default for run_data is 'cyl'.
  'ortho' is a good projection if you want to see curved RA and Dec lines, while 'cyl' matches up with the drapery plot.
  Other projections include 'ortho', 'cyl', 'merc', 'moll' and more, detailed in the matplotlib package Basemap.

  transparency: this dictates the transparency of the drapery plot.  This is useful for overlaying the map with the drapery of the same data.
  histogram_file_basename, map_file_basename, drapery_file_basename: Defaults are None.
  If not given an input, run_data will extract the ObsID from filename_Q and append it with variable information.
  (The ObsID extraction finds an integer in the filename.  If the file does not contain an integer number,
  you must set a basename for all plots you're creating.)

  directory: If you want to put the plots in a different directory than the current working directory.
  If a directory with the name inputted doesn't exist, one will be created.


plothealpix_map.mapping(ra, dec, data_vals, var_name, obsID=None, map_file_name=None,
            projection='ortho', save_show='show', full_image=True)

  This module creates a map projected onto a sphere of some variable.

  ra, dec: The coordinates of the data points.

  data_vals: The 3rd dimension. This info will be plotted with color.

  var_name: The name of the 3rd dimension.  This will be used to title the plot.
  If 'sources' is selected, the generated map will be a scatterplot, instead of a contourf plot.

  obsID: The Obs ID of the area of data

  map_file_name: The name of the file to be written out.

  projection: The map projection for the mapping module plothealpix_map. The default for run_data is 'cyl'.
  'ortho' is a good projection if you want to see curved RA and Dec lines, while 'cyl' matches up with the drapery plot.
  Other projections include 'ortho', 'cyl', 'merc', 'moll' and more, detailed in the matplotlib package Basemap.

  save_show: Options are 'save' and 'show'. 'Save' option saves the plots you've chosen, as whatever type
  of file has been specified.

  full_image: If True, all the data will be displayed.  If False, the maximum rectangular shape that fits inside the circle of data is selected. This area is selected with the module fits_data_extraction.cutout_square(ra, dec).

stokes_tester.plotstokesQU(obsID, x_stokes, y_stokes, ra, dec,
                 new_histogram_filename=None, directory=None,
                 save_show='show', projection='ortho')

  This module takes data lists 'x_stokes' and 'y_stokes', which are obtained using the module stokes_math.math(signal_I, signal_Q, signal_U), and creates a 2D histogram of the data.

  x_stokes and y_stokes are the projections onto a single axis, of the relationship between the Q and U axes.  These axes are at 45-degree angles to each other, making this transform necessary.

  The output plot shows the density of points at some angle between 0 and pi, and some magnitude (distance from zero).  Ideally, the data should appear as a roughly semicircular shape.

  obsID: obsID: The Obs ID of the area of data

  ra, dec: The coordinates of the data points.

  new_histogram_filename: The name of the written-out file

  directory: location of the file.  Default is current working directory.

  save_show: Options are 'save' and 'show'. 'Save' option saves the plot as whatever type
  of file has been specified.

  projection: The map projection for the mapping module plothealpix_map. The default for run_data is 'cyl'.
  'ortho' is a good projection if you want to see curved RA and Dec lines, while 'cyl' matches up with the drapery plot.
  Other projections include 'ortho', 'cyl', 'merc', 'moll' and more, detailed in the matplotlib package Basemap.

lic_maps.LIC(obsID, x_stokes, y_stokes, ra, dec, dpi, size, length=31, full_image=True,
        disp_drapery='save', name_of_plot='flow-image.png', transparency=1)

  This module creates a 'drapery' plot of stokes linear polarizations.  The direction of polarization at every point is represented visually as a texture or draping of the angles.
