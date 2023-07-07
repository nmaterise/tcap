# -*- encoding: utf-8 -*-
"""
Post processing tools for plotting, viewing states, etc.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import cm
from scipy.optimize import curve_fit
import scipy.sparse as scsp
import datetime
import pickle as pk
import subprocess
import h5py as hdf
from mpl_toolkits.mplot3d import Axes3D
import re as regex
# import vtk


def get_linestyle_cycler():
    """
    Returns a linestyle cycler for plotting
    """

    # Different types of dashing styles
    linestyle_cycle = [
     (0, (1, 10)),
     (0, (1, 1)),
     (0, (1, 1)),
     (0, (5, 10)),
     (0, (5, 5)),
     (0, (5, 1)),
     (0, (3, 10, 1, 10)),
     (0, (3, 5, 1, 5)),
     (0, (3, 1, 1, 1)),
     (0, (3, 5, 1, 5, 1, 5)),
     (0, (3, 10, 1, 10, 1, 10)),
     (0, (3, 1, 1, 1, 1, 1))]

    return linestyle_cycle

 
def get_alpha_color_cycler(alpha=0.5):
    """
    Returns color_cycler default with transparency fraction set to alpha
    """

    # Get the color cycler as a hex
    color_cycle_hex = plt.rcParams['axes.prop_cycle'].by_key()['color']
    hex2rgb = lambda hx: [int(hx[0:2],16)/256., \
                          int(hx[2:4],16)/256., \
                          int(hx[4:6],16)/256.]
    color_cycle_rgb = [hex2rgb(cc[1:]) for cc in color_cycle_hex]

    return [(*cc, alpha) for cc in color_cycle_rgb]


def get_marker_cycler():
    """
    Returns a marker style cycler for plotting
    """

    # Different marker icons
    marker_cycle = ['D', 'x', 'o', 's', 'v', '^', '*', '>', '<', 'p']

    return marker_cycle

def set_axes_num_format(ax, fmt, axis='x'):
    """
    Sets the number format for the x and y axes in a 1D plot
    """
    
    # Check for the axis x or y and act accordingly
    if axis == 'x':
        ax.xaxis.set_major_formatter(mpl.ticker.StrMethodFormatter(fmt))
    elif axis == 'y':
        ax.yaxis.set_major_formatter(mpl.ticker.StrMethodFormatter(fmt))
    elif axis == 'xy':
        ax.xaxis.set_major_formatter(mpl.ticker.StrMethodFormatter(fmt))
        ax.yaxis.set_major_formatter(mpl.ticker.StrMethodFormatter(fmt))
    else:
        raise ValueError('axis (%s) not a valid input.' % axis)


def init_subplots(fsize=20, tight_layout=True):
    """
    Initializes the figure and axes objects to default values and font sizes
    """

    # Get the figure and axes objects
    fig, ax = plt.subplots(1, 1, tight_layout=tight_layout)

    # Set the axes fonts
    set_axes_fonts(ax, fsize)

    return fig, ax


def set_leg_outside(ax, fsize):
    """
    Sets the legend location outside
    """
    
    # Shrink current axis by 20%
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    
    # Put a legend to the right of the current axis
    hdls, legs = ax.get_legend_handles_labels()
    leg = ax.legend(hdls, legs, fontsize=fsize,\
                    loc='center left', bbox_to_anchor=(1, 0.5))

    return leg


def set_leg_hdls_lbs(ax, fsize, loc='best'):
    """
    Set the legend handles and labels
    """

    hdls, legs = ax.get_legend_handles_labels()
    ax.legend(hdls, legs, loc=loc, fontsize=fsize)


def set_axes_fonts(ax, fsize):
    """
    Set axes font sizes because it should be abstracted away
    """
    
    for tick in ax.get_xticklabels():
        tick.set_fontsize(fsize)
    for tick in ax.get_yticklabels():
        tick.set_fontsize(fsize)


def set_axes_fonts_3d(ax, fsize):
    """
    Set axes font sizes because it should be abstracted away
    """
    
    for tick in ax.get_xticklabels():
        tick.set_fontsize(fsize)
    for tick in ax.get_yticklabels():
        tick.set_fontsize(fsize)
    for tick in ax.get_zticklabels():
        tick.set_fontsize(fsize)


def set_xaxis_rot(ax, angle=45):
    """
    Rotate the x-axis labels
    """
        
    for tick in ax.get_xticklabels():
        tick.set_rotation(angle)


def write_fig_to_file(fig, fname, leg=None, ax=None,\
                      is_leg_outside=True, \
                      format='pdf', fsize=20, \
                      leg_loc='best', transparent=True):
    """
    Writes a figure object to file, sets legends accordingly
    
    Parameters:
    ----------

    fig:                    matplotlib figure object
    fname:                  path to output file 
    leg:                    matplotlib legends object
    ax:                     matplotlib axes object
    is_leg_outside:         True if legends set outside of figure, else false
    format:                 output file format
    fsize:                  legend fontsize
    leg_loc:                legend location
    transparent:            set background to transparent

    """

    # Check for no legends
    if leg is None:
        fig.savefig(fname, format=format, transparent=transparent)
    
    # Otherwise save with legends
    else:
        ## Check for setting legend outside
        if is_leg_outside:
            fig.savefig(fname, format=format, \
                  bbox_extra_artists=(leg, ), bbox_inches='tight', \
                  transparent=transparent)
        else:
            ### Check for the ax object to set the legends
            if ax is None:
                pass
            else:
                set_leg_hdls_lbs(ax, fsize, loc=leg_loc)
            
            fig.savefig(fname, format=format, transparent=transparent)


def setup_3d_surf_plot(x, y, z, labels=[], xyzstrs={'x':'','y':'','z':''},
                       fsize=20, zscale='linear', show_leg=False):
    """
    Setup the 3d surface plot for a set of data ordered by expectation for
    plotting with plt.plot_surface(0
    
    Parameters:
    ----------

    x, y:       1d arrays for the abscissa
    z:          3d array for the ordinates
    labels:     1d list of legend labels for the various 2d slices of z
    xyzstrs:    dictionary of x,y,z axes labels
    fsize:      axes labels and numbers fontsize
    zscale:     scale on the z-axis data, applied as log10 or none
    show_leg:   display legend if True

    Returns:
    -------

    fig, ax:    figure and axes objects

    """

    # Get the figure and axes objects
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    
    # Set the axes fonts sizes
    set_axes_fonts_3d(ax, fsize)
    ax.set_xlabel(xyzstrs['x'], fontsize=fsize, labelpad=15)
    ax.set_ylabel(xyzstrs['y'], fontsize=fsize, labelpad=15)
    ax.set_zlabel(xyzstrs['z'], fontsize=fsize, labelpad=15)
    
    # Get the meshgrid points
    X, Y = np.meshgrid(x, y)

    # Plot the results
    alpha_cycle = get_alpha_color_cycler(alpha=0.5)
    if zscale == 'log':
        surfs = [ax.plot_surface(X, Y, np.log10(zz), \
                                 lw=0, label=r'%s' % labels[idx], \
                                 color=alpha_cycle[idx%len(alpha_cycle)]) \
            for idx, zz in enumerate(z)]

    elif zscale == 'linear':
        surfs = [ax.plot_surface(X, Y, zz, lw=0, label=r'%s' % ll) \
            for zz, ll in zip(z, labels)]

    else:
        raise ValueError('zscale (%s) no recognized.' % zscale)

    # Trying this hack
    for surf in surfs:
        surf._facecolors2d=surf._facecolors3d
        surf._edgecolors2d=surf._edgecolors3d

    # Set the legend outside
    if show_leg:
        leg = set_leg_outside(ax, fsize)

        return fig, ax, leg
    else:
        return fig, ax, None


def complex_str_itoj(fname):
    """
    Replaces all instances of i with j in complex data generated by COMSOL or
    other software
    """

    # Open the file for reading
    with open(fname, 'r') as iid:
        sdata = iid.read()
    iid.close()

    # Find all instances of i, if any proceed, otherwise close the file
    finds = regex.findall('i', sdata) 

    # Replace the i's with j's
    if len(finds) > 0:
        print('Replacing i with j in (%s) ...' % fname)

        cleaned_data = regex.sub('i', 'j', sdata) 

        # Reopen the file and overwrite the contents with the replaced string 
        with open(fname, 'w') as oid:
            oid.write(cleaned_data)
        oid.close()

def stderr_fill(ax, xval, yval, yerr, fill_color, alpha=0.5):
    """
    Shaded region for standard deviation on a linear plot
    """

    # Shaded region command
    ax.fill_between(xval, yval + yerr, yval - yerr, \
            facecolor=fill_color, alpha=alpha)


def moving_avg(x, stride=2):
    """
    Computes the moving average of x, replacing x_i with
    (x_i-stride + x_i+stride)/2 by padding the ends with stride-zeros
    """

    # Read in the data as a numpy array
    x = np.asarray(x)
    
    # Generate the zero padding and hstack with the input
    zpad = np.zeros(stride)
    xpadded = np.hstack((zpad, x, zpad))
    xout = np.zeros(x.size)

    # Compute the moving averages
    for xidx in range(x.size):
        xout[xidx] = (xpadded[xidx-stride] + xpadded[xidx+stride]) / 2

    return xout


def plot_2d_cmap(x, y, z, fname,
                 xstr='', ystr='',
                 tstr='', cbar_str='',
                 cmap=cm.viridis, nlevels=100,
                 use_imshow=False, show_cbar=False,
                 zlims=None, ndisp=3):
    """
    Plot 2D colormap data such that 

         -----------------------
         |                     |
         |                     |
    y    |          z          |
         |                     |
         |                     |
         ----------------------- 
    
                    x
    
    Parameters:
    ----------

    x, y:       independent variables
    z:          resulting data, z = z(x, y) 
    fname:      output figure filename, relative path with file extension
    xstr:       x-axis label
    ystr:       y-axis label
    tstr:       title label
    cbar_str:   lable for color bar
    cmap:       colormap to use on the 2D figure output

    """

    # Setup the color map, normalizations, etc
    if zlims is not None:
        norm = mpl.colors.Normalize(zlims[0], zlims[1])
    else:
        norm = mpl.colors.Normalize(z.min(), z.max())

    # Fontsize
    fsize = 36; tsize = 26;
    fsize2 = 32

    # Setup the figure and axes
    fig, ax = plt.subplots(1, 1, figsize=(10, 8))

    # Read in the data as numpy arrays
    x = np.asarray(x)
    y = np.asarray(y)
    z = np.asarray(z)

    # Use imshow to display results
    if use_imshow:
        plt1 = ax.imshow(z, cmap=cmap, norm=norm) #, cmap=cmap, norm=norm, \
                         #extent=(x.min(), x.max(), y.min(), y.max()))
        if show_cbar:
            # Set the x,y ticks to the appropriate values
            ax.set_xticks(np.arange(len(x), dtype=int))
            ax.set_yticks(np.arange(len(y), dtype=int))
            ax.set_xticklabels([int(xx) for xx in x])
            ax.set_yticklabels([int(xx) for xx in y])

            # Iterate over all of the values and plot on their pixels
            for i in range(z.shape[0]):
                for j in range(z.shape[1]):
                    text = ax.text(j, i, f'{z[i, j]:.{ndisp}g}', ha='center',
                            # va='center', color='tab:gray', fontsize=fsize2)
                            va='center', color='k', fontsize=fsize2)
            cbar = fig.colorbar(plt1, ax=ax, format=f'%.{ndisp}g')
            cbar.ax.set_title(cbar_str, fontsize=fsize, y=1.025)
            cbar.ax.tick_params(labelsize=fsize)
        else:
            # Iterate over all of the values and plot on their pixels
            ax.set_xticks(np.arange(len(x), dtype=int))
            ax.set_yticks(np.arange(len(y), dtype=int))
            ax.set_xticklabels([int(xx) for xx in x])
            ax.set_yticklabels([int(xx) for xx in y])
            for i in range(z.shape[0]):
                for j in range(z.shape[1]):
                    text = ax.text(j, i, f'{z[i, j]:.{ndisp}g}', ha="center",
                            va="center", color="tab:gray", fontsize=fsize2)
    else:
        plt1 = ax.contourf(x, y, z, nlevels, cmap=cmap, norm=norm)
        # Set the color bar, offset the title slightly from top
        cbar = fig.colorbar(plt1, ax=ax)
        cbar.ax.set_title(cbar_str, fontsize=fsize, y=1.025)
        cbar.ax.tick_params(labelsize=fsize)

    ax.set_xlabel(xstr, fontsize=fsize)
    ax.set_ylabel(ystr, fontsize=fsize)
    ax.set_title(tstr, fontsize=tsize)

    # Set the axis tick labels to a reasonable size
    set_axes_fonts(ax, fsize)

    # Write the results to file
    fig.tight_layout()
    fig.savefig(fname, format='pdf', transparent=True) 
        

def gnuplot_term_trace(x, y, title=None):
    """
    Uses gnuplot to plot data in the terminal
    """

    # Open a process and pipe commands to gnuplot
    gnuplot = subprocess.Popen(["/usr/bin/gnuplot"], stdin=subprocess.PIPE)

    # Set the terminal to dumb
    gnuplot.stdin.write(b"set term dumb 79 25\n")
    
    # Plot the data, starting with the lines command
    tstr = str.encode(title) if title != None else b"Line1"
    gnuplot.stdin.write(b"plot '-' using 1:2 title '%b' with dots\n" % tstr)
    
    # Iterate over the input data
    for i, j in zip(x, y):
        gnuplot.stdin.write(b"%f %f\n" % (i, j))

    # Write the execute command, then flush 
    gnuplot.stdin.write(b"e\n")
    gnuplot.stdin.flush()
    gnuplot.stdin.write(b"quit")


def gnuplot_dumb_traces(xlist, ylist, 
                        tlist=None, linespoints='points',
                        xlim=[], ylim=[], xscale='', yscale=''):
    """
    Plot multiple traces, each input is assumed a list

    Example:
    -------

    gnuplot_dumb_traces([x1, x2], [y1, y2], tlist=['y1(x1)', 'y2(x2)'])

    """

    # Open a process and pipe commands to gnuplot
    gnuplot = subprocess.Popen(["/usr/bin/gnuplot"], stdin=subprocess.PIPE)

    # Set the terminal to dumb
    gnuplot.stdin.write(b"set term dumb 79 25\n")

    # Set the x and y ranges if not empty
    if xlim != []:
        xrange_str = 'set xrange [%g:%g]\n' % (xlim[0], xlim[1])
        gnuplot.stdin.write(str.encode(xrange_str))

    # Set the x and y ranges if not empty
    if ylim != []:
        yrange_str = 'set yrange [%g:%g]\n' % (ylim[0], ylim[1])
        gnuplot.stdin.write(str.encode(yrange_str))

    # Set the scale of x and y axes
    if xscale == 'log':
        xscale_str = 'set logscale x\n'
        gnuplot.stdin.write(str.encode(xscale_str))
    if yscale == 'log':
        yscale_str = 'set logscale y\n'
        gnuplot.stdin.write(str.encode(yscale_str))

    # Get the number of entries
    Nitems = len(xlist)

    # Set the linespoints option
    lp = str.encode(linespoints)

    # Set the title strings
    tstrs = [str.encode(tt) for tt in tlist] if tlist != None \
            else [b'Line%b' % (str.encode('%d' % i)) for i in range(Nitems)]

    # Set the titles in one shot
    plot_str = b"plot %b" % (b",".join([b"'-' u 1:2 t '%b' w %b" \
                % (tt, lp) for tt in tstrs]))
    plot_str += b"\n"
    gnuplot.stdin.write(plot_str) 

    # Iterate over the inputs
    for x, y in zip(xlist, ylist):
    
        # Iterate over the input data
        for i, j in zip(x, y):
            gnuplot.stdin.write(b"%f %f\n" % (i, j))

        # Write the execute command, then flush 
        gnuplot.stdin.write(b"e\n")

    # Flush gnuplot to view
    gnuplot.stdin.flush()
    gnuplot.stdin.write(b"quit")


def gnuplot_dumb_boxes(x, y, title=None):
    """
    Plot a single bar chart with x, y data
    """

    # Open a process and pipe commands to gnuplot
    gnuplot = subprocess.Popen(["/usr/bin/gnuplot"], stdin=subprocess.PIPE)

    # Set the terminal to dumb
    gnuplot.stdin.write(b"set term dumb 79 25\n")
    
    # Plot the data, starting with the lines command
    gnuplot.stdin.write(b"plot '-' using 1:2 t '' w boxes\n")
    
    # Iterate over the input data
    for i, j in zip(x, y):
        gnuplot.stdin.write(b"%f %f\n" % (i, j))

    # Write the execute command, then flush 
    gnuplot.stdin.write(b"e\n")
    if title != None:
        print('\n%s\n' % title)
    gnuplot.stdin.flush()
    gnuplot.stdin.write(b"quit")


def write_to_hdf(data_list, key_list, fname):
    """
    Writes a list of arrays to data sets with names in key_list to an hdf5 at
    the path given by fname
    
    Parameters:
    ----------

    data_list:      list of arrays to write to data sets
    key_list:       list of string keys for each data set in data_list in the
                    same order as data_list 
    fname:          path to the hdf5 file output

    """

    # Open the new hdf5 file and get the file handle
    fid = hdf.File(fname, 'w')

    # Iterate over the data_list and key_list
    for d, k in zip(data_list, key_list):
        fid.create_dataset(k, data=d)

    # Close the file handle
    fid.close()


def twinx_plot(xlist, ylist, xstr, ystrs, fname):
    """
    Plots data (x1, y1) and (x2, y2) on shared x-axis using differnt y-axes
    
    Parameters:
    ----------

    xlist:      [x1_array, x2_array]
    ylist:      [y1_array, y2_array]
    xstr:       x-axis label
    ystrs:      ['yleft-label', 'yright-label']
    fname:      output file name

    """
    # Default font sizes
    fsize = 20
    
    # Create a separate figure for each key
    fig, ax1 = plt.subplots(1, 1, tight_layout=True)
    
    # Set the axes fontsizes
    set_axes_fonts(ax1, fsize)
    
    # Set the colors for each line
    color1 = 'tab:blue'
    color2 = 'tab:red'
    
    # Plot the result of the key
    ax1.plot(xlist[0], ylist[0], color=color1)
    ax1.set_xlabel(xstr, fontsize=fsize)
    ax1.set_ylabel(ystrs[0], fontsize=fsize, color=color1)
    
    # Rotate the x-axis labels
    set_xaxis_rot(ax1)
    
    # Set the axis color
    ax1.tick_params(axis='y', labelcolor=color1)
    
    # Get the twin-axis
    ax2 = ax1.twinx()
    set_axes_fonts(ax2, fsize)
    ax2.plot(xlist[1], ylist[1], color=color2)
    ax2.set_xlabel(xstr, fontsize=fsize)
    ax2.set_ylabel(ystrs[1], fontsize=fsize, color=color2)
    ax2.tick_params(axis='y', labelcolor=color2)


    # Write the figure to file
    # Assuming the format is the same as the file extension
    fig.savefig(fname, format=fname.split('.')[-1], transparent=True)

# def import_vti(fname):
#     """
#     Imports a vtk file (.vti) and renders it for interactive manipulation
#     """
# 
#     # Get the vtk colors
#     colors = vtk.vtkNamedColors()
# 
#     # Read the source file.
#     reader = vtk.vtkXMLImageDataReader()
#     reader.SetFileName(fname)
# 
#     # Create the mapper that creates graphics elements
#     mapper = vtk.vtkDataSetMapper()
#     mapper.SetInputConnection(reader.GetOutputPort())
# 
#     # Create the Actor
#     actor = vtk.vtkActor()
#     actor.SetMapper(mapper)
#     # show the edges of the image grid
#     actor.GetProperty().SetRepresentationToWireframe()
#     actor.GetProperty().SetColor(colors.GetColor3d("DarkSalmon"))
# 
#     # Create the Renderer
#     renderer = vtk.vtkRenderer()
#     renderer.AddActor(actor)
#     renderer.ResetCamera()
#     renderer.SetBackground(colors.GetColor3d("Silver"))
# 
#     # Create the RendererWindow
#     renderer_window = vtk.vtkRenderWindow()
#     renderer_window.AddRenderer(renderer)
# 
#     # Create the RendererWindowInteractor and display the vti file
#     interactor = vtk.vtkRenderWindowInteractor()
#     interactor.SetRenderWindow(renderer_window)
#     interactor.Initialize()
#     interactor.Start()
# 
# 
# def import_vtu(fname):
#     """
#     Imports a vtk file (.vtu) unstructured grid file
#     """
# 
#     # Get the colors defined in vtk
#     colors = vtk.vtkNamedColors()
# 
#     # Start the reader
#     reader = vtk.vtkUnstructuredGridReader()
#     reader.SetFileName(fname)
#     reader.Update()
# 
#     # Extract the edges
#     extractEdges = vtk.vtkExtractEdges()
#     extractEdges.SetInputConnection(reader.GetOutputPort())
# 
#     # Get the legend values
#     legendValues = vtk.vtkVariantArray()
#     it = reader.GetOutput().NewCellIterator()
#     it.InitTraversal()
# 
#     # Start traversing the cells
#     while not it.IsDoneWithTraversal():
#         cell = vtk.vtkGenericCell()
#         it.GetCell(cell)
#         cellName = vtk.vtkCellTypes.GetClassNameFromTypeId(cell.GetCellType())
#         print(cellName, "NumberOfPoints:", 
#                 cell.GetNumberOfPoints(), "CellDimension:",
#         cell.GetCellDimension())
#         legendValues.InsertNextValue(cellName)
#         it.GoToNextCell()
# 
#     # Tube the edges
#     tubes = vtk.vtkTubeFilter()
#     tubes.SetInputConnection(extractEdges.GetOutputPort())
#     tubes.SetRadius(.05)
#     tubes.SetNumberOfSides(21)
# 
#     edgeMapper = vtk.vtkPolyDataMapper()
#     edgeMapper.SetInputConnection(tubes.GetOutputPort())
#     edgeMapper.SetScalarRange(0, 26)
# 
#     edgeActor = vtk.vtkActor()
#     edgeActor.SetMapper(edgeMapper)
#     edgeActor.GetProperty().SetSpecular(.6)
#     edgeActor.GetProperty().SetSpecularPower(30)
# 
#     # Glyph the points
#     sphere = vtk.vtkSphereSource()
#     sphere.SetPhiResolution(21)
#     sphere.SetThetaResolution(21)
#     sphere.SetRadius(.08)
# 
#     pointMapper = vtk.vtkGlyph3DMapper()
#     pointMapper.SetInputConnection(reader.GetOutputPort())
#     pointMapper.SetSourceConnection(sphere.GetOutputPort())
#     pointMapper.ScalingOff()
#     pointMapper.ScalarVisibilityOff()
# 
#     pointActor = vtk.vtkActor()
#     pointActor.SetMapper(pointMapper)
#     pointActor.GetProperty().SetDiffuseColor(colors.GetColor3d("Banana"))
#     pointActor.GetProperty().SetSpecular(.6)
#     pointActor.GetProperty().SetSpecularColor(1.0, 1.0, 1.0)
#     pointActor.GetProperty().SetSpecularPower(100)
# 
#     # Label the points
#     labelMapper = vtk.vtkLabeledDataMapper()
#     labelMapper.SetInputConnection(reader.GetOutputPort())
#     labelActor = vtk.vtkActor2D()
#     labelActor.SetMapper(labelMapper)
# 
#     # The geometry
#     geometryShrink = vtk.vtkShrinkFilter()
#     geometryShrink.SetInputConnection(reader.GetOutputPort())
#     geometryShrink.SetShrinkFactor(.8)
# 
#     # NOTE: We must copy the originalLut because the CategoricalLegend
#     # needs an indexed lookup table, but the geometryMapper uses a
#     # non-index lookup table
#     # categoricalLut = vtk.vtkLookupTable()
#     # originalLut = reader.GetOutput().GetCellData().GetScalars().GetLookupTable()
# 
#     # categoricalLut.DeepCopy(originalLut)
#     # categoricalLut.IndexedLookupOn()
# 
#     geometryMapper = vtk.vtkDataSetMapper()
#     geometryMapper.SetInputConnection(geometryShrink.GetOutputPort())
#     geometryMapper.SetScalarModeToUseCellData()
#     geometryMapper.SetScalarRange(0, 11)
# 
#     geometryActor = vtk.vtkActor()
#     geometryActor.SetMapper(geometryMapper)
#     geometryActor.GetProperty().SetLineWidth(3)
#     geometryActor.GetProperty().EdgeVisibilityOn()
#     geometryActor.GetProperty().SetEdgeColor(0, 0, 0)
# 
#     # Legend
#     # for v in range(0, legendValues.GetNumberOfTuples()):
#     #     categoricalLut.SetAnnotation(legendValues.GetValue(v),
#     #             legendValues.GetValue(v).ToString())
#     # legend = vtk.vtkCategoryLegend()
#     # legend.SetScalarsToColors(categoricalLut)
#     # legend.SetValues(legendValues)
#     # legend.SetTitle("Cell Type")
#     # legend.GetBrush().SetColor(colors.GetColor4ub("Silver"))
# 
#     # placeLegend = vtk.vtkContextTransform()
#     # placeLegend.AddItem(legend)
#     # placeLegend.Translate(640 - 20, 480 - 12 * 16)
# 
#     contextView = vtk.vtkContextView()
#     # contextView.GetScene().AddItem(placeLegend)
# 
#     renderer = contextView.GetRenderer()
# 
#     renderWindow = contextView.GetRenderWindow()
# 
#     renderWindowInteractor = vtk.vtkRenderWindowInteractor()
#     renderWindowInteractor.SetRenderWindow(renderWindow)
# 
#     renderer.AddActor(geometryActor)
#     renderer.AddActor(labelActor)
#     renderer.AddActor(edgeActor)
#     renderer.AddActor(pointActor)
#     renderer.SetBackground(colors.GetColor3d("SlateGray"))
# 
#     aCamera = vtk.vtkCamera()
#     aCamera.Azimuth(-40.0)
#     aCamera.Elevation(50.0)
# 
#     renderer.SetActiveCamera(aCamera)
#     renderer.ResetCamera()
# 
#     renderWindow.SetSize(640, 480)
#     renderWindow.Render()
# 
#     renderWindowInteractor.Start()

